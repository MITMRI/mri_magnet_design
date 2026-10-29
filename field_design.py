# -*- coding: utf-8 -*-
"""
Field-driven design of an aligned-ring Halbach magnet, in the spirit of CoilGen (target field -> continuous
source distribution -> discrete layout -> refinement) instead of a genetic search from random designs.

1. Target field: Bx = B0 on check points covering the target cylinder (or sphere, --sphere-d). Bx is harmonic, so its extremes over the
   region are on its surface and checking the surface is enough.
2. Sensitivity matrix: field per magnet of a Halbach ring at every candidate (layer radius, z) on a fine grid,
   the magnet equivalent of CoilGen's stream-function basis functions.
3. Continuous source (LP): magnet density along z for each layer that minimizes the peak field deviation, with
   density limited by the ring pitch and the total magnet budget. This is the "theoretical" magnet distribution.
4. Discrete placement (MILP): whole ring slices on the grid (all layers of a slice share z, slices >= cube + 3 mm
   apart) with an integer magnet count per ring, minimizing the same peak deviation. This is the initial design.
5. Gradient refinement: slice z, per-ring radius and rotation offset. Gradients are analytic for point dipoles; the
   difference to the exact cuboid field (magpylib) is added as a correction refreshed every few steps, so the
   objective is exact and the gradient is the dipole one (cube corrections are ~(a/d)^4). First Adam on a smooth
   peak-to-peak (log-sum-exp), then SLSQP on the exact minimax form, then screened +-1 / move-one-magnet changes.
6. --ga instead refines the same initial design with a genetic algorithm (pymoo) over the same variables and
   constraints, for comparison.

All runs stop at --time-limit seconds (default 600).

usage:
    python field_design.py out.xlsx [--layers 3] [--max-magnets 650] [--min-b0 45] [--time-limit 600] [--ga]
"""

import argparse
import time
import numpy as np
import pandas as pd
import magpylib as magpy
from scipy.optimize import linprog, milp, LinearConstraint, Bounds, minimize
from scipy.spatial.transform import Rotation as R
import magsimulator
from refine_design import COLS

T0 = time.time()


def log(msg):
    print(f'[{time.time()-T0:6.1f} s] {msg}', flush=True)


def ppm(B):
    return 1e6*(B.max()-B.min())/B.mean()


def ring_field(r, z, n, phase, sensors, a, br):
    """Exact Bx (mT) on the sensors from a Halbach ring of n cubes at radius r; a mirror ring at -z unless z == 0."""
    th = phase + 2*np.pi*np.arange(n)/n
    zs = [0.0] if abs(z) < 1e-9 else [z, -z]
    pos = np.vstack([np.stack([r*np.cos(th), r*np.sin(th), np.full(n, zz)], axis=1) for zz in zs])
    rotz = np.concatenate([2*th for _ in zs])
    n_m, n_s = len(pos), len(sensors)
    B = magpy.getB('Cuboid', observers=np.tile(sensors, (n_m, 1)), position=np.repeat(pos, n_s, axis=0),
                   orientation=R.from_euler('z', np.repeat(rotz, n_s)[:, None]), magnetization=(br, 0, 0),
                   dimension=(a,)*3)
    return B[:, 0].reshape(n_m, n_s).sum(axis=0)


class Geometry:
    def __init__(self, args):
        self.a, self.br = args.cube, args.br
        self.pitch = args.cube + args.slice_wall                       # slice spacing along z
        self.dr = np.sqrt(2)*args.cube + args.layer_wall              # layer spacing (in-plane cube diagonal + resin)
        self.rmin = args.bore/2 + args.bore_wall + args.cube/np.sqrt(2)
        self.rmax = args.max_od/2 - args.od_wall - args.cube/np.sqrt(2)
        self.zmax = (args.max_length - args.cube)/2
        self.radii = [self.rmin + k*self.dr for k in range(args.layers)]
        assert self.radii[-1] <= self.rmax + 1e-9, f'{args.layers} layers do not fit within OD {args.max_od} mm'
        self.max_magnets, self.min_b0 = args.max_magnets, args.min_b0
        self.fitc = np.sqrt(3)*args.cube/2    # magnets fit around a ring if 2*pi/n >= 2*asin(fitc/r) (as in magsimulator)

    def nmax(self, r):
        return int(np.floor(2*np.pi/(2*np.arcsin(self.fitc/r)) + 1e-9))


# ---------------------------------------------------------------------------------------------------------------
# steps 2-4: sensitivity matrix, continuous LP, discrete MILP

def sensitivity(geo, zgrid, sensors):
    """A[l, j, :] = Bx per magnet of a ring in layer l at z = zgrid[j] (with its mirror ring unless z = 0)."""
    A = np.zeros((len(geo.radii), len(zgrid), len(sensors)))
    for l, r in enumerate(geo.radii):
        nref = geo.nmax(r)
        for j, z in enumerate(zgrid):
            A[l, j] = ring_field(r, z, nref, 0.0, sensors, geo.a, geo.br)/nref
    return A


def continuous_source(geo, sensors, b0, dz=1.0):
    """LP: magnet density per layer along z (magnets per dz cell on each side) minimizing max |Bx - B0|."""
    zc = np.arange(dz/2, geo.zmax, dz)
    A = sensitivity(geo, zc, sensors)                         # mirror pairs; variable = magnets at +z in the cell
    L, nz, ns = A.shape
    nv = L*nz + 1
    M = A.reshape(L*nz, ns).T                                  # (ns, L*nz)
    ub = np.concatenate([np.full(nz, geo.nmax(r)*dz/geo.pitch) for r in geo.radii])
    c = np.zeros(nv); c[-1] = 1
    A_ub = np.vstack([np.hstack([M, -np.ones((ns, 1))]), np.hstack([-M, -np.ones((ns, 1))]),
                      np.concatenate([2*np.ones(L*nz), [0]])[None, :]])
    b_ub = np.concatenate([np.full(ns, b0), np.full(ns, -b0), [geo.max_magnets]])
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=list(zip(np.zeros(L*nz), ub)) + [(0, None)], method='highs')
    rho = res.x[:-1].reshape(L, nz)
    return zc, rho/dz, 2e6*res.x[-1]/b0, M @ res.x[:-1]       # density in magnets/mm (one side), ppm, field


def discrete_placement(geo, sensors, b0, dz, nmin, time_limit):
    """MILP: ring slices on a z grid with integer magnets per ring, minimizing max |Bx - B0|."""
    zg = np.concatenate([[0.0], np.arange(geo.pitch/2, geo.zmax + 1e-9, dz)])
    A = sensitivity(geo, zg, sensors)
    L, nz, ns = A.shape
    mult = np.where(zg == 0, 1, 2)
    nm = np.array([geo.nmax(r) for r in geo.radii])
    # variables: n[l, j] (int), y[l, j] (ring present), s[j] (slice present), t
    i_n, i_y, i_s, i_t = 0, L*nz, 2*L*nz, 2*L*nz + nz
    nv = i_t + 1
    rows, lo, hi = [], [], []

    def add(row, l_, h_):
        rows.append(row); lo.append(l_); hi.append(h_)

    M = A.reshape(L*nz, ns).T
    for k in range(ns):                                         # |B_k - B0| <= t
        row = np.zeros(nv); row[i_n:i_y] = M[k]
        row[i_t] = -1; add(row, -np.inf, b0)
        row = row.copy(); row[i_t] = 1; add(row, b0, np.inf)
    for l in range(L):
        for j in range(nz):
            q = l*nz + j
            row = np.zeros(nv); row[i_n+q] = 1; row[i_y+q] = -nm[l]; add(row, -np.inf, 0)   # n <= nmax*y
            row = np.zeros(nv); row[i_n+q] = 1; row[i_y+q] = -nmin; add(row, 0, np.inf)     # n >= nmin*y
            row = np.zeros(nv); row[i_y+q] = 1; row[i_s+j] = -1; add(row, -np.inf, 0)       # y <= s
    for j in range(nz):                                         # slices at least one pitch apart
        win = np.where((zg >= zg[j]) & (zg < zg[j] + geo.pitch - 1e-9))[0]
        if len(win) > 1:
            row = np.zeros(nv); row[i_s+win] = 1; add(row, -np.inf, 1)
    row = np.zeros(nv); row[i_n:i_y] = np.tile(mult, L); add(row, -np.inf, geo.max_magnets)
    c = np.zeros(nv); c[i_t] = 1
    integrality = np.concatenate([np.ones(2*L*nz + nz), [0]])
    ub = np.concatenate([np.repeat(nm, nz), np.ones(L*nz + nz), [np.inf]])
    res = milp(c, constraints=LinearConstraint(np.array(rows), lo, hi), integrality=integrality,
               bounds=Bounds(np.zeros(nv), ub), options={'time_limit': time_limit, 'disp': False, 'mip_rel_gap': 1e-3})
    x = res.x
    n = np.round(x[i_n:i_y]).astype(int).reshape(L, nz)
    rings = [dict(layer=l, r=geo.radii[l], z=zg[j], n=int(n[l, j]), phase=0.0)
             for j in range(nz) for l in range(L) if n[l, j] > 0]
    return rings, 2e6*x[i_t]/b0, res


# ---------------------------------------------------------------------------------------------------------------
# step 5: ring design with analytic dipole gradients

class Design:
    """Rings grouped into slices; slice 0 is at z = 0 if present. Variables: slice z (> 0), ring radius, ring phase."""

    PH = 10.0   # phase scaling so that all variables move by similar amounts (mm)

    def __init__(self, rings, geo):
        self.geo = geo
        zs = sorted({round(abs(rg['z']), 6) for rg in rings})
        self.has_center = zs[0] == 0
        self.zpos = [z for z in zs if z > 0]
        self.rings = []
        for rg in sorted(rings, key=lambda q: (abs(q['z']), q['r'])):
            z = round(abs(rg['z']), 6)
            k = -1 if z == 0 else self.zpos.index(z)
            self.rings.append(dict(k=k, r=float(rg['r']), n=int(rg['n']), phase=float(rg['phase'])))

    def copy(self):
        d = object.__new__(Design)
        d.geo, d.has_center, d.zpos = self.geo, self.has_center, list(self.zpos)
        d.rings = [dict(rg) for rg in self.rings]
        return d

    def z(self, rg):
        return 0.0 if rg['k'] < 0 else self.zpos[rg['k']]

    def count(self):
        return sum(rg['n']*(1 if rg['k'] < 0 else 2) for rg in self.rings)

    # variable vector
    def pack(self):
        return np.array(self.zpos + [rg['r'] for rg in self.rings] + [rg['phase']*self.PH for rg in self.rings])

    def unpack(self, x):
        nz, nr = len(self.zpos), len(self.rings)
        self.zpos = list(x[:nz])
        for i, rg in enumerate(self.rings):
            rg['r'], rg['phase'] = x[nz+i], x[nz+nr+i]/self.PH

    def magnets(self):
        """Positions, rotation angles (rad), and per magnet: ring index and z sign (for gradients)."""
        pos, rot, ring, sgn = [], [], [], []
        for i, rg in enumerate(self.rings):
            th = rg['phase'] + 2*np.pi*np.arange(rg['n'])/rg['n']
            z = self.z(rg)
            for s in ([1.0] if rg['k'] < 0 else [1.0, -1.0]):
                pos.append(np.stack([rg['r']*np.cos(th), rg['r']*np.sin(th), np.full(len(th), s*z)], axis=1))
                rot.append(2*th); ring.append(np.full(len(th), i)); sgn.append(np.full(len(th), s))
        return np.vstack(pos), np.concatenate(rot), np.concatenate(ring), np.concatenate(sgn)

    def field_exact(self, sensors):
        return sum(ring_field(rg['r'], self.z(rg), rg['n'], rg['phase'], sensors, self.geo.a, self.geo.br)
                   for rg in self.rings)

    def field_dipole(self, sensors, grad=False):
        """Point-dipole Bx and, with grad, its Jacobian wrt the variable vector (ns x nv)."""
        pos, rot, ring, sgn = self.magnets()
        k = self.geo.br*self.geo.a**3/(4*np.pi)
        d = sensors[None, :, :] - pos[:, None, :]                       # (nm, ns, 3)
        r2 = (d**2).sum(-1); ir = 1/np.sqrt(r2); ir3 = ir**3; ir5 = ir3/r2; ir7 = ir5/r2
        mx, my = np.cos(rot)[:, None], np.sin(rot)[:, None]
        dx, dy, dz = d[..., 0], d[..., 1], d[..., 2]
        md = mx*dx + my*dy
        Bx = k*(3*md*dx*ir5 - mx*ir3)
        if not grad:
            return Bx.sum(0)
        # dBx/dd_j, then dBx/dp = -dBx/dd
        g = [k*(3*m*dx*ir5 + 3*md*(j == 0)*ir5 - 15*md*dx*dj*ir7 + 3*mx*dj*ir5)
             for j, (m, dj) in enumerate([(mx, dx), (my, dy), (0.0, dz)])]
        gpx, gpy, gpz = -g[0], -g[1], -g[2]
        dBdmx, dBdmy = k*(3*dx*dx*ir5 - ir3), k*(3*dy*dx*ir5)
        th = (rot/2)[:, None]
        r = np.array([self.rings[i]['r'] for i in ring])[:, None]
        dBdr = gpx*np.cos(th) + gpy*np.sin(th)
        dBdth = r*(-gpx*np.sin(th) + gpy*np.cos(th)) + 2*(-dBdmx*np.sin(2*th) + dBdmy*np.cos(2*th))
        dBdz = gpz*sgn[:, None]
        nz, nr = len(self.zpos), len(self.rings)
        J = np.zeros((len(sensors), nz + 2*nr))
        for i, rg in enumerate(self.rings):
            sel = ring == i
            J[:, nz+i] = dBdr[sel].sum(0)
            J[:, nz+nr+i] = dBdth[sel].sum(0)/self.PH
            if rg['k'] >= 0:
                J[:, rg['k']] += dBdz[sel].sum(0)
        return Bx.sum(0), J

    # constraints g(x) >= 0 and their Jacobian
    def constraints(self, grad=False):
        geo, a = self.geo, self.geo.a
        nz, nr = len(self.zpos), len(self.rings)
        nv = nz + 2*nr
        g, J = [], []

        def add(val, coef):
            row = np.zeros(nv)
            for j, c in coef:
                row[j] += c
            g.append(val); J.append(row)

        z = self.zpos
        if nz:
            add(z[0] - (geo.pitch if self.has_center else geo.pitch/2), [(0, 1)])
            for i in range(nz-1):
                add(z[i+1] - z[i] - geo.pitch, [(i+1, 1), (i, -1)])
            add(geo.zmax - z[-1], [(nz-1, -1)])
        for i, rg in enumerate(self.rings):
            r = rg['r']
            add(r - geo.rmin, [(nz+i, 1)])
            add(geo.rmax - r, [(nz+i, -1)])
            c = geo.fitc/r
            add(2*np.pi/rg['n'] - 2*np.arcsin(c), [(nz+i, 2*geo.fitc/(r*r*np.sqrt(1-c*c)))])
            nxt = [j for j, q in enumerate(self.rings) if q['k'] == rg['k'] and q['r'] > r]
            if nxt:
                j = min(nxt, key=lambda j: self.rings[j]['r'])
                add(self.rings[j]['r'] - r - geo.dr, [(nz+j, 1), (nz+i, -1)])
        return (np.array(g), np.array(J)) if grad else np.array(g)

    def ledger(self):
        pos, rot, _, _ = self.magnets()
        rows = [[p[0], p[1], p[2], 0, 0, np.rad2deg(t), 0, 0, 0, 0, 0, self.geo.a, ''] for p, t in zip(pos, rot)]
        return pd.DataFrame(rows, columns=COLS)

    def summary(self):
        out = []
        for k in ([-1] if self.has_center else []) + list(range(len(self.zpos))):
            rs = [rg for rg in self.rings if rg['k'] == k]
            z = 0.0 if k < 0 else self.zpos[k]
            out.append(f'z={z:6.2f}: ' + ', '.join(f"{rg['n']}@r{rg['r']:.1f}" for rg in rs))
        return '\n'.join(out)


def score(d, sensors):
    """Exact ppm, or inf if a constraint or B0 is violated."""
    B = d.field_exact(sensors)
    if B.mean() < d.geo.min_b0 - 1e-6 or d.constraints().min() < -1e-6:
        return np.inf
    return ppm(B)


def refine_adam(d, sensors, steps, deadline, lr=0.05, refresh=25):
    """Gradient descent (Adam) on a smooth peak-to-peak with quadratic penalties for constraints and B0."""
    x = d.pack(); m = np.zeros_like(x); v = np.zeros_like(x)
    best_x, best = x.copy(), score(d, sensors)
    corr = None
    for it in range(steps):
        if time.time() > deadline:
            break
        d.unpack(x)
        if it % refresh == 0:
            corr = d.field_exact(sensors) - d.field_dipole(sensors)
        B, J = d.field_dipole(sensors, grad=True)
        B = B + corr
        mean = B.mean(); dev = 1e6*(B - mean)/mean
        dmean = J.mean(0)
        ddev = 1e6*(J*mean - B[:, None]*dmean[None, :])/mean**2
        beta = 0.05 + 0.45*it/max(1, steps-1)                   # anneal the log-sum-exp towards a hard max
        wp = np.exp(beta*(dev - dev.max())); wp /= wp.sum()
        wm = np.exp(-beta*(dev - dev.min())); wm /= wm.sum()
        grad = wp @ ddev - wm @ ddev
        g, Jg = d.constraints(grad=True)
        viol = np.minimum(g, 0)
        grad += 2*1e4*(viol @ Jg)
        db = mean - d.geo.min_b0
        if db < 0.02:                                            # keep a little above the minimum B0
            grad += 2*1e6*(db - 0.02)*dmean
        m = 0.9*m + 0.1*grad; v = 0.999*v + 0.001*grad**2
        x = x - lr*(m/(1 - 0.9**(it+1)))/(np.sqrt(v/(1 - 0.999**(it+1))) + 1e-8)
        if it % refresh == refresh - 1:
            d.unpack(x)
            c = score(d, sensors)
            if c < best:
                best, best_x = c, x.copy()
    d.unpack(best_x)
    return best


def refine_minimax(d, sensors, deadline, outer=6, maxiter=100):
    """SLSQP on the exact minimax problem: min t s.t. |dev_i| <= t, B0 >= min, geometry; analytic dipole Jacobian."""
    best = score(d, sensors)
    for _ in range(outer):
        if time.time() > deadline:
            break
        x0 = d.pack()
        corr = d.field_exact(sensors) - d.field_dipole(sensors)
        cache = {}

        def fld(x):
            key = x.tobytes()
            if key not in cache:
                d.unpack(x[:-1])
                B, J = d.field_dipole(sensors, grad=True)
                cache.clear(); cache[key] = (B + corr, J)
            return cache[key]

        def cons(x):
            B, _ = fld(x)
            mean = B.mean(); dev = 1e6*(B - mean)/mean; t = x[-1]
            d.unpack(x[:-1])
            return np.concatenate([t - dev, t + dev, [mean - d.geo.min_b0], d.constraints()])

        def cjac(x):
            B, J = fld(x)
            mean = B.mean(); dmean = J.mean(0)
            ddev = 1e6*(J*mean - B[:, None]*dmean[None, :])/mean**2
            d.unpack(x[:-1])
            _, Jg = d.constraints(grad=True)
            ns, nv = J.shape
            top = np.vstack([-ddev, ddev, dmean[None, :], Jg])
            tcol = np.concatenate([np.ones(2*ns), np.zeros(1 + len(Jg))])[:, None]
            return np.hstack([top, tcol])

        B0 = d.field_exact(sensors)
        xs = np.concatenate([x0, [0.5*ppm(B0)]])
        res = minimize(lambda x: x[-1], xs, jac=lambda x: np.eye(len(x))[-1], method='SLSQP',
                       constraints=[{'type': 'ineq', 'fun': cons, 'jac': cjac}],
                       options={'maxiter': maxiter, 'ftol': 1e-6})
        d.unpack(res.x[:-1])
        c = score(d, sensors)
        if c < best - 0.5:
            best = c
        else:
            d.unpack(x0)
            break
    return best


def count_moves(d, sensors, deadline, top=6):
    """Screened +-1 magnet changes per ring and single-magnet moves between rings, each re-tuned by SLSQP."""
    best = score(d, sensors)
    while time.time() < deadline:
        cands = []
        idx = range(len(d.rings))
        moves = [((i, s),) for i in idx for s in (-1, 1)]
        moves += [((i, -1), (j, 1)) for i in idx for j in idx if i != j and (d.rings[i]['k'] < 0) == (d.rings[j]['k'] < 0)]
        for mv in moves:
            e = d.copy()
            for i, s in mv:
                e.rings[i]['n'] += s
            if min(rg['n'] for rg in e.rings) < 4 or e.count() > d.geo.max_magnets:
                continue
            cands.append((score(e, sensors), mv, e))
        cands.sort(key=lambda c: c[0])
        found = None
        for _, mv, e in cands[:top]:
            if time.time() > deadline:
                break
            c = refine_minimax(e, sensors, deadline, outer=2, maxiter=40)
            if c < best - 1:
                best, found = c, (mv, e)
        if found is None:
            break
        mv, e = found
        d.zpos, d.rings = e.zpos, e.rings
        log('  ' + ', '.join(f"ring {i}: {'+' if s > 0 else '-'}1" for i, s in mv) + f' -> {d.count()} magnets, {best:.0f} ppm')
    return best


# ---------------------------------------------------------------------------------------------------------------
# step 6: GA refinement of the same initial design, for comparison

def refine_ga(d, sensors, deadline, pop=60, seed=1, pool=None):
    from pymoo.core.problem import Problem
    from pymoo.algorithms.soo.nonconvex.ga import GA
    from pymoo.optimize import minimize as pminimize
    from pymoo.core.population import Population
    from pymoo.termination import get_termination

    # genome: the continuous variables plus each ring's magnet count (rounded), so the GA can do what the
    # gradient stages and count moves do together
    n0 = np.array([rg['n'] for rg in d.rings], float)
    x0 = np.concatenate([d.pack(), n0]); nz, nr = len(d.zpos), len(d.rings)
    lo = np.concatenate([np.maximum(np.array(d.zpos) - 6, d.geo.pitch/2), np.full(nr, d.geo.rmin), np.full(nr, -np.pi)*Design.PH,
                         np.maximum(n0 - 4, 4)])
    hi = np.concatenate([np.minimum(np.array(d.zpos) + 6, d.geo.zmax), np.full(nr, d.geo.rmax), np.full(nr, np.pi)*Design.PH,
                         n0 + 4])
    lo, hi = np.minimum(lo, x0), np.maximum(hi, x0)
    base = d.copy()

    class P(Problem):
        def __init__(self):
            super().__init__(n_var=len(x0), n_obj=1, n_ieq_constr=1, xl=lo, xu=hi)

        def _evaluate(self, X, out, *a, **kw):
            res = (pool.map if pool else map)(_ga_eval, [(base, x, sensors) for x in X])
            res = np.array(list(res))
            out['F'], out['G'] = res[:, :1], res[:, 1:]

    rng = np.random.default_rng(seed)
    sig = np.concatenate([np.full(len(x0) - nr, 0.3), np.full(nr, 0.7)])
    X0 = np.vstack([x0] + [np.clip(x0 + rng.normal(0, 1, len(x0))*sig, lo, hi) for _ in range(pop - 1)])
    algo = GA(pop_size=pop, sampling=Population.new('X', X0), eliminate_duplicates=True)
    secs = max(1, int(deadline - time.time()))
    res = pminimize(P(), algo, get_termination('time', f'{secs//3600:02d}:{secs//60%60:02d}:{secs%60:02d}'),
                    seed=seed, verbose=False)
    x = np.atleast_2d(res.X if res.X is not None else x0)[0]
    _ga_apply(d, x)
    return score(d, sensors), res.algorithm.n_gen, res.algorithm.evaluator.n_eval


def _ga_apply(d, x):
    nr = len(d.rings)
    d.unpack(x[:-nr])
    for rg, n in zip(d.rings, np.round(x[-nr:]).astype(int)):
        rg['n'] = int(n)


def _ga_eval(args):
    base, x, sensors = args
    e = base.copy(); _ga_apply(e, x)
    B = e.field_exact(sensors)
    viol = max(0, e.geo.min_b0 - B.mean()) + np.maximum(0, -e.constraints()).sum() + max(0, e.count() - e.geo.max_magnets)
    return ppm(B), viol


# ---------------------------------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--layers', type=int, default=3)
    ap.add_argument('--max-magnets', type=int, default=650)
    ap.add_argument('--min-b0', type=float, default=45)
    ap.add_argument('--max-length', type=float, default=150)
    ap.add_argument('--bore', type=float, default=100)
    ap.add_argument('--max-od', type=float, default=200)
    ap.add_argument('--cube', type=float, default=6.35)
    ap.add_argument('--br', type=float, default=1320)
    ap.add_argument('--cyl-d', type=float, default=40); ap.add_argument('--cyl-h', type=float, default=40)
    ap.add_argument('--sphere-d', type=float, help='target a sphere of this diameter instead of the cylinder')
    ap.add_argument('--bore-wall', type=float, default=3); ap.add_argument('--od-wall', type=float, default=3)
    ap.add_argument('--slice-wall', type=float, default=3); ap.add_argument('--layer-wall', type=float, default=2)
    ap.add_argument('--sensors', type=int, default=400, help='check points on the target cylinder (z >= 0 half used)')
    ap.add_argument('--milp-dz', type=float, default=1.0, help='z grid of candidate slices (mm)')
    ap.add_argument('--milp-time', type=float, default=120)
    ap.add_argument('--nmin', type=int, default=8, help='fewest magnets in a ring that is used')
    ap.add_argument('--b0-margin', type=float, default=0.1, help='mT above --min-b0 targeted by LP/MILP')
    ap.add_argument('--adam-steps', type=int, default=400)
    ap.add_argument('--time-limit', type=float, default=600)
    ap.add_argument('--ga', action='store_true', help='refine the initial design with a GA instead of gradients')
    ap.add_argument('--pool', type=int, default=0, help='worker processes for the GA')
    ap.add_argument('--start', help='skip steps 1-4 and start from this ledger')
    args = ap.parse_args()

    deadline = T0 + args.time_limit
    geo = Geometry(args)
    if args.sphere_d:
        sensors = magsimulator.define_sensor_points_on_sphere(args.sensors, args.sphere_d/2, [0, 0, 0]).position
    else:
        sensors = magsimulator.define_sensor_points_on_cylinder(args.sensors, args.cyl_d/2, args.cyl_h, [0, 0, 0]).position
    sensors = sensors[sensors[:, 2] >= 0]
    b0 = args.min_b0 + args.b0_margin
    log(f'{len(sensors)} check points, layer radii ' + ', '.join(f'{r:.2f}' for r in geo.radii) +
        f' mm (max {[geo.nmax(r) for r in geo.radii]} magnets/ring), slice pitch {geo.pitch} mm, z <= {geo.zmax:.2f} mm')
    stem = args.out.rsplit('.', 1)[0]
    times = {}

    if args.start:
        from refine_design import rings_from_ledger
        rings = rings_from_ledger(pd.read_excel(args.start))
        d = Design([dict(r=rg['r'], z=rg['z'], n=rg['n'], phase=rg['phase']) for rg in rings if rg['z'] >= -1e-6], geo)
        log(f'start {args.start}: {d.count()} magnets, {score(d, sensors):.0f} ppm')
    else:
        t = time.time()
        zc, rho, lp_ppm, _ = continuous_source(geo, sensors, b0)
        np.savez(stem + '_source.npz', z=zc, density=rho, radii=geo.radii, ppm=lp_ppm)
        times['lp'] = time.time() - t
        log(f'continuous source (LP): {lp_ppm:.1f} ppm at {b0} mT with {2*rho.sum():.0f} magnets '
            f'(per layer {np.round(2*rho.sum(1)).astype(int).tolist()})')
        t = time.time()
        rings, milp_ppm, res = discrete_placement(geo, sensors, b0, args.milp_dz, args.nmin, args.milp_time)
        times['milp'] = time.time() - t
        d = Design(rings, geo)
        B = d.field_exact(sensors)
        log(f'discrete placement (MILP, {res.message.strip()}): {milp_ppm:.0f} ppm linearized, exact {ppm(B):.0f} ppm, '
            f'B0 {B.mean():.2f} mT, {d.count()} magnets, {len(d.zpos) + d.has_center} slices\n' + d.summary())
        d.ledger().to_excel(stem + '_initial.xlsx', index=False)

    if args.ga:
        from multiprocessing import Pool
        pool = Pool(args.pool) if args.pool > 1 else None
        t = time.time()
        c, ngen, neval = refine_ga(d, sensors, deadline, pool=pool)
        times['ga'] = time.time() - t
        log(f'GA: {c:.0f} ppm after {ngen} generations, {neval} evaluations')
    else:
        t = time.time()
        c = refine_adam(d, sensors, args.adam_steps, deadline)
        times['adam'] = time.time() - t
        log(f'gradient descent (Adam): {c:.0f} ppm')
        t = time.time()
        c = refine_minimax(d, sensors, deadline)
        times['slsqp'] = time.time() - t
        log(f'minimax (SLSQP): {c:.0f} ppm')
        t = time.time()
        c = count_moves(d, sensors, deadline)
        c = refine_minimax(d, sensors, deadline)
        times['moves'] = time.time() - t
        log(f'magnet count moves: {c:.0f} ppm')

    d.ledger().to_excel(args.out, index=False)
    B = d.field_exact(sensors)
    log(f'done: {d.count()} magnets, {ppm(B):.0f} ppm on the optimization points, B0 {B.mean():.2f} mT -> {args.out}\n'
        + d.summary() + '\nstage times: ' + ', '.join(f'{k} {v:.0f} s' for k, v in times.items()))


if __name__ == '__main__':
    main()
