# -*- coding: utf-8 -*-
"""
Passive shim design for a slide-in shim sleeve: small cubes in trays inside the bore, placed with the same pipeline as
field_design.py (target field -> continuous source density -> discrete magnets -> gradient refinement).

Hardware (see SHIMMING.md and make_shim_cad.py): a printed sleeve slides into the bore. Its outside has --trays axial
grooves, closed by the bore wall of the holder; a printed tray slides into each groove from the +z end. Each tray is
printed for one shim solution, so a cube can sit at any z and, lying flat, can turn to any angle about the radial axis
(magnetization in the tangential-z plane), or stand with its magnetization radially in or out. The magnet array never
comes apart.

Steps, per field map:
  1. LP: continuous shim magnetization density on a fine grid of tray positions (6 directions per point), minimizing
     the peak-to-peak field over the check points (minimax, as in field_design.continuous_source) with B0 >= --min-b0
     and at most --max-cubes cubes.
  2. Placement: walk each tray and put a cube wherever the accumulated |density| passes k + offset; it points along the
     local density (flat at the matching angle, or standing radial if the density is mostly radial); beyond
     --max-cubes the cubes with the least density around them are dropped. Several offsets are tried (--starts).
  3. Gradient refinement: SLSQP minimax over every cube's z and angle, with tray spacing, z range and B0 constraints.
     The cubes are 25+ mm from the check points, so the field model is a point dipole (analytic, cheap); the final
     result is checked with the exact cuboid field (magpylib).
  4. The shim cubes' own errors (Br, magnetization angle) are propagated linearly and reported as the typical and 90th
     percentile shimmed result.

Field maps:
  --map FILE      measured map: .csv with columns x, y, z (mm) and Bx (mT; B0 along x), or the MRI4ALL field mapper
                  format (two .mat files, --map pos.mat --map-field field.mat: 'pos' (n,3) mm, 'field_matrix' (n,>=1) mT)
  --simulate N    N simulated as-built builds of the ledger (Br, angle and position errors as in field_design --robust),
                  checked on the target cylinder; for studying what shimming buys before anything is built

usage:
  python shim_design.py results/field_design/FD_robust.xlsx results/shim/FD_robust_sim --simulate 3
  python shim_design.py results/field_design/FD_robust.xlsx results/shim/measured --map map.csv
Outputs: <out>_b<i>.csv (one row per shim cube: tray, z, orientation, position, magnetization), <out>_b<i>.npz,
<out>.log.
"""

import argparse
import os
import time
import numpy as np
import pandas as pd
import magpylib as magpy
from scipy.optimize import linprog, minimize
from scipy.spatial.transform import Rotation as R
import magsimulator


class Sleeve:
    """tray geometry: cube centre radius, tray axes and the spacing a cube needs along its tray"""
    def __init__(self, args):
        self.a, self.n, self.br, self.zmax = args.cube, args.trays, args.br, args.zmax
        # tray outer face just inside the sleeve OD (sleeve OD = bore - 2 fit); pocket open outwards, cube on the floor
        self.r_face = args.bore/2 - args.fit - args.tray_gap
        self.r = self.r_face - args.pocket_tol - self.a/2
        self.wall = args.tray_wall
        self.mom = self.br*1e-3/(4e-7*np.pi)*(self.a*1e-3)**3          # dipole moment, A m^2
        th = 2*np.pi*(np.arange(self.n) + 0.5)/self.n
        self.th = th
        self.er = np.stack([np.cos(th), np.sin(th), 0*th], 1)
        self.et = np.stack([-np.sin(th), np.cos(th), 0*th], 1)

    def spacing(self, kind):
        # flat cubes may turn to any angle: allow the full diagonal along the tray
        return np.where(kind == 0, self.a*np.sqrt(2), self.a) + self.wall

    def mdir(self, tr, kind, psi):
        """unit magnetization: kind 0 = flat (cos psi along z + sin psi along the tangent), +1/-1 = radially out/in"""
        flat = np.cos(psi)[:, None]*np.array([0, 0, 1.0]) + np.sin(psi)[:, None]*self.et[tr]
        return np.where((kind == 0)[:, None], flat, kind[:, None]*self.er[tr])

    def pos(self, tr, z):
        return self.r*self.er[tr] + z[:, None]*np.array([0, 0, 1.0])


def dipole_bx(pos, mvec, P):
    """Bx (mT) at points P (m,3) mm from point dipoles at pos (n,3) mm with moments mvec (n,3) A m^2 -> (n, m)"""
    r = P[None] - pos[:, None]
    d = np.linalg.norm(r, axis=2)
    rh = r/d[..., None]
    mr = (mvec[:, None]*rh).sum(2)
    return 1e5*(3*rh[..., 0]*mr - mvec[:, None, 0])/d**3


def exact_bx(sl, tr, z, kind, psi, P):
    """exact cuboid field of the shim cubes (faces aligned with the tray, flat cubes turned by psi about the radial axis)"""
    if len(z) == 0:
        return np.zeros(len(P))
    cs = []
    for t, zz, k, p in zip(tr, z, kind, psi):
        base = R.from_matrix(np.stack([sl.er[t], sl.et[t], [0, 0, 1.0]], 1))   # local x radial, y tangential, z axial
        if k == 0:
            rot, mloc = base*R.from_euler('x', -p), (0, 0, sl.br)
        else:
            rot, mloc = base, (k*sl.br, 0, 0)
        cs.append(magpy.magnet.Cuboid(magnetization=mloc, dimension=(sl.a,)*3, position=sl.pos(np.array([t]), np.array([zz]))[0],
                                      orientation=rot))
    return magpy.getB(cs, P, sumup=True)[:, 0]


def ppm(b):
    return 1e6*np.ptp(b)/b.mean()


def lp_density(sl, b0, P, nmax, b0min):
    """continuous shim density: points every a/2 along each tray, 6 directions, minimax ppm. Returns the LP ppm,
    the z grid and the moment density per tray and point (in cubes)"""
    zs = np.arange(-sl.zmax, sl.zmax + 1e-9, sl.a/2)
    tr = np.repeat(np.arange(sl.n), len(zs))
    pos = sl.pos(tr, np.tile(zs, sl.n))
    ez = np.tile([0, 0, 1.0], (len(tr), 1))
    D = np.stack([sl.er[tr], -sl.er[tr], sl.et[tr], -sl.et[tr], ez, -ez], 1)
    A = np.stack([dipole_bx(pos, sl.mom*D[:, o], P) for o in range(6)], 1).reshape(-1, len(P))
    nv, ns = A.shape
    sc = b0.mean()/1e6
    Ac = (A - A.mean(1, keepdims=True)).T
    dev0 = b0 - b0.mean()
    one = np.zeros((len(pos), nv + 1))     # a point every a/2 holds at most half a cube
    for s in range(len(pos)):
        one[s, 6*s:6*s+6] = 1
    Aub = np.vstack([np.hstack([Ac, -sc*np.ones((ns, 1))]), np.hstack([-Ac, -sc*np.ones((ns, 1))]), one,
                     np.r_[-A.mean(1), 0][None], np.r_[np.ones(nv), 0][None]])
    bub = np.r_[-dev0, dev0, 0.5*np.ones(len(pos)), b0.mean() - b0min, nmax]
    res = linprog(np.r_[np.zeros(nv), 1], A_ub=Aub, b_ub=bub, bounds=(0, None), method='highs')
    x = res.x[:-1].reshape(len(pos), 6)
    return res.x[-1], zs, (x[..., None]*D).sum(1).reshape(sl.n, len(zs), 3)


def place(sl, zs, V, offset, nmax):
    """cubes from the density: per tray, one wherever the accumulated |moment| passes k + offset; at most nmax, keeping
    the ones with the most density around them"""
    tr, zc, kind, psi, wt = [], [], [], [], []
    for t in range(sl.n):
        v = V[t]
        cum = np.cumsum(np.linalg.norm(v, axis=1))
        for k in range(int(np.floor(cum[-1] + 1 - offset))):
            j = min(np.searchsorted(cum, k + offset), len(zs) - 1)
            m = v[max(j - 2, 0):j + 3].sum(0)
            if np.linalg.norm(m) < 1e-9:
                continue
            cr, ct, cz = m @ sl.er[t], m @ sl.et[t], m[2]
            tr.append(t); zc.append(zs[j]); wt.append(np.linalg.norm(m))
            if abs(cr) > np.hypot(ct, cz):
                kind.append(np.sign(cr)); psi.append(0.0)
            else:
                kind.append(0.0); psi.append(np.arctan2(ct, cz))
    tr, zc, kind, psi = np.array(tr, int), np.array(zc, float), np.array(kind, float), np.array(psi, float)
    top = np.sort(np.argsort(wt)[::-1][:nmax])
    tr, zc, kind, psi = tr[top], zc[top], kind[top], psi[top]
    # spread cubes that landed too close together along the tray (the refinement keeps them apart from here on)
    for t in np.unique(tr):
        idx = np.flatnonzero(tr == t)
        idx = idx[np.argsort(zc[idx])]
        for i0, i1 in zip(idx[:-1], idx[1:]):
            need = (sl.spacing(kind[i0:i0+1]) + sl.spacing(kind[i1:i1+1]))[0]/2
            zc[i1] = max(zc[i1], zc[i0] + need)
    keep = np.abs(zc) <= sl.zmax
    return tr[keep], zc[keep], kind[keep], psi[keep]


def refine(sl, b0, P, tr, z, kind, psi, b0min, iters):
    """SLSQP minimax over z (every cube) and psi (flat cubes): min t s.t. |dev| <= t, B0 >= b0min, tray spacing, |z| <= zmax"""
    o = np.lexsort((z, tr))
    tr, z, kind, psi = tr[o], z[o], kind[o], psi[o]
    n = len(z)
    if n == 0:
        return tr, z, kind, psi
    flat = np.flatnonzero(kind == 0)
    sp = sl.spacing(kind)
    pairs = [(i, i + 1) for i in range(n - 1) if tr[i] == tr[i + 1]]

    def unpack(x):
        ps = psi.copy()
        ps[flat] = x[n:n + len(flat)]
        return x[:n], ps

    def fields(x):
        zz, ps = unpack(x)
        return dipole_bx(sl.pos(tr, zz), sl.mom*sl.mdir(tr, kind, ps), P)

    def cons(x):
        Bt = b0 + fields(x).sum(0)
        m = Bt.mean()
        dev = 1e6*(Bt - m)/m
        zz, _ = unpack(x)
        g = [zz[j] - zz[i] - (sp[i] + sp[j])/2 for i, j in pairs]
        return np.r_[x[-1] - dev, x[-1] + dev, m - b0min, g, sl.zmax - zz, zz + sl.zmax]

    def cjac(x):
        Bi = fields(x)
        Bt = b0 + Bi.sum(0)
        m = Bt.mean()
        zz, ps = unpack(x)
        h = 1e-4
        J = []
        for k in range(n + len(flat)):          # each variable moves one cube
            i = k if k < n else flat[k - n]
            z1, p1 = zz[i:i+1].copy(), ps[i:i+1].copy()
            if k < n:
                z1 += h
            else:
                p1 += h
            J.append((dipole_bx(sl.pos(tr[i:i+1], z1), sl.mom*sl.mdir(tr[i:i+1], kind[i:i+1], p1), P)[0] - Bi[i])/h)
        J = np.array(J)
        dm = J.mean(1)
        dd = 1e6*(J*m - Bt[None]*dm[:, None])/m**2
        ns, nv = len(P), len(x)
        top = np.zeros((2*ns + 1, nv))
        top[:ns, :-1], top[ns:2*ns, :-1] = -dd.T, dd.T
        top[:2*ns, -1] = 1
        top[-1, :-1] = dm
        g = np.zeros((len(pairs) + 2*n, nv))
        for r_, (i, j) in enumerate(pairs):
            g[r_, i], g[r_, j] = -1, 1
        g[len(pairs) + np.arange(n), np.arange(n)] = -1
        g[len(pairs) + n + np.arange(n), np.arange(n)] = 1
        return np.vstack([top, g])

    x0 = np.r_[z, psi[flat], 0.0]
    x0[-1] = 0.5*ppm(b0 + fields(x0).sum(0))*1.05
    res = minimize(lambda x: x[-1], x0, jac=lambda x: np.eye(len(x))[-1], method='SLSQP',
                   constraints=[{'type': 'ineq', 'fun': cons, 'jac': cjac}], options={'maxiter': iters, 'ftol': 1e-6})
    zz, ps = unpack(res.x)
    return tr, zz, kind, np.mod(ps + np.pi, 2*np.pi) - np.pi


def shim_error_spread(sl, tr, z, kind, psi, B_after, P, sigma_br, sigma_ang, n=400, seed=0):
    """typical / p90 ppm after shimming when the shim cubes themselves have Br and angle errors (linear, point dipoles)"""
    if len(z) == 0:
        return ppm(B_after), ppm(B_after)
    pos = sl.pos(tr, z)
    m = sl.mdir(tr, kind, psi)
    # two directions perpendicular to each cube's magnetization
    ref = np.where(np.abs(m[:, 2:3]) < 0.9, np.array([[0, 0, 1.0]]), np.array([[1.0, 0, 0]]))
    p1 = np.cross(m, ref); p1 /= np.linalg.norm(p1, axis=1, keepdims=True)
    p2 = np.cross(m, p1)
    G = np.concatenate([dipole_bx(pos, sl.mom*m, P)*sigma_br, dipole_bx(pos, sl.mom*p1, P)*np.deg2rad(sigma_ang),
                        dipole_bx(pos, sl.mom*p2, P)*np.deg2rad(sigma_ang)])
    s = B_after[None] + np.random.default_rng(seed).standard_normal((n, len(G))) @ G
    e = 1e6*np.ptp(s, 1)/s.mean(1)
    return np.median(e), np.percentile(e, 90)


def design(sl, Bopt, Popt, Bchk, Pchk, args, log):
    """one shim solution for one field map: Bopt on Popt drives the optimization, Bchk on Pchk is the check"""
    t0 = time.time()
    lp, zs, V = lp_density(sl, Bopt, Popt, args.max_cubes, args.min_b0)
    best = None
    for off in np.linspace(0.5, 0.25, args.starts) if args.starts > 1 else [0.5]:
        sol = place(sl, zs, V, off, args.max_cubes)
        sol = refine(sl, Bopt, Popt, *sol, args.min_b0, args.iters)
        after = Bchk + exact_bx(sl, *sol, Pchk)
        ok = after.mean() >= args.min_b0 - 0.05 and len(sol[1]) <= args.max_cubes
        log(f'    start {off:.2f}: {len(sol[1])} cubes, {ppm(after):.0f} ppm, B0 {after.mean():.2f} mT' +
            ('' if ok else ' (rejected: B0 too low or too many cubes)'))
        if ok and (best is None or ppm(after) < best[0]):
            best = (ppm(after), sol, after)
    return lp, best, time.time() - t0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('ledger', help='magnet ledger of the main array (for --simulate, and recorded in the log)')
    ap.add_argument('out', help='output prefix')
    ap.add_argument('--map', help='measured field map: .csv (x, y, z, Bx) or MRI4ALL pos .mat (with --map-field)')
    ap.add_argument('--map-field', help='MRI4ALL field .mat (field_matrix, Bx in the first column) for --map pos.mat')
    ap.add_argument('--simulate', type=int, default=0, help='number of simulated as-built builds to shim')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--sigma-br', type=float, default=0.8, help='%%, main array (simulated builds) and shim cubes')
    ap.add_argument('--sigma-ang', type=float, default=0.7, help='deg, main array (simulated builds) and shim cubes')
    ap.add_argument('--pos-tol', type=float, default=0.25, help='+- mm, main array magnet position (simulated builds)')
    ap.add_argument('--cyl-d', type=float, default=40); ap.add_argument('--cyl-h', type=float, default=40)
    ap.add_argument('--sensors', type=int, default=2000, help='check points on the target cylinder (simulated builds)')
    ap.add_argument('--cube', type=float, default=3.175, help='shim cube size, mm (1/8")')
    ap.add_argument('--br', type=float, default=1320, help='shim cube remanence, mT')
    ap.add_argument('--trays', type=int, default=36)
    ap.add_argument('--zmax', type=float, default=70, help='shim cubes within +- this z, mm')
    ap.add_argument('--bore', type=float, default=100)
    ap.add_argument('--fit', type=float, default=0.2, help='radial clearance between the sleeve and the bore, mm')
    ap.add_argument('--tray-gap', type=float, default=0.1, help='tray outer face below the sleeve OD, mm')
    ap.add_argument('--pocket-tol', type=float, default=0.2, help='clearance added to the cube in the tray pockets, mm')
    ap.add_argument('--tray-wall', type=float, default=0.8, help='resin between pockets along a tray, mm')
    ap.add_argument('--max-cubes', type=int, default=60)
    ap.add_argument('--min-b0', type=float, default=44.8, help='mT, after shimming')
    ap.add_argument('--starts', type=int, default=3, help='placement offsets tried before refinement')
    ap.add_argument('--iters', type=int, default=200, help='SLSQP iterations per start')
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    logf = open(args.out + '.log', 'w')
    def log(s):
        print(s, flush=True); logf.write(s + '\n'); logf.flush()

    sl = Sleeve(args)
    log(f'shim sleeve: {sl.n} trays, cube centres at r {sl.r:.2f} mm, {sl.a} mm cubes (Br {sl.br:.0f} mT, {sl.mom:.4f} A m2), '
        f'|z| <= {sl.zmax}, max {args.max_cubes} cubes, B0 >= {args.min_b0} mT')

    maps = []       # (name, Bopt, Popt, Bchk, Pchk)
    if args.map:
        if args.map_field:
            from scipy.io import loadmat
            P = np.asarray(loadmat(args.map)['pos'], float)
            B = np.asarray(loadmat(args.map_field)['field_matrix'], float)[:, 0]
        else:
            M = pd.read_csv(args.map)
            P, B = M[['x', 'y', 'z']].values.astype(float), M['Bx'].values.astype(float)
        maps.append(('map', B, P, B, P))
        log(f'field map {args.map}: {len(P)} points, B0 {B.mean():.3f} mT, {ppm(B):.0f} ppm')
    if args.simulate:
        L = pd.read_excel(args.ledger)
        S = magsimulator.define_sensor_points_on_cylinder(args.sensors, args.cyl_d/2, args.cyl_h, [0, 0, 0]).position
        cubes = magsimulator.build_magnet_collection(L, [args.br, 0, 0]).children
        B = magpy.getB(cubes, S, sumup=False)[..., 0]
        sens = [B*args.sigma_br/100]
        for mloc in ((0, args.br, 0), (0, 0, args.br)):          # magnetization tipped towards local y, z
            for c in cubes:
                c.magnetization = mloc
            sens.append(magpy.getB(cubes, S, sumup=False)[..., 0]*np.deg2rad(args.sigma_ang))
        for c in cubes:
            c.magnetization = (args.br, 0, 0)
        for e in np.eye(3):                                      # moving a magnet by +d = moving the sensors by -d
            sens.append((magpy.getB(cubes, S - 0.01*e, sumup=False)[..., 0] - B)/0.01*args.pos_tol/np.sqrt(3))
        G = np.concatenate(sens)
        builds = B.sum(0)[None] + np.random.default_rng(args.seed).standard_normal((args.simulate, len(G))) @ G
        for i, Bb in enumerate(builds):
            maps.append((f'b{i}', Bb[::2], S[::2], Bb, S))
        log(f'{args.simulate} simulated builds of {args.ledger} (Br {args.sigma_br}%, {args.sigma_ang} deg, '
            f'+-{args.pos_tol} mm), {len(S)} check points on the {args.cyl_d:g} x {args.cyl_h:g} mm cylinder')

    for name, Bopt, Popt, Bchk, Pchk in maps:
        log(f'{name}: unshimmed {ppm(Bchk):.0f} ppm, B0 {Bchk.mean():.2f} mT')
        lp, best, dt = design(sl, Bopt, Popt, Bchk, Pchk, args, log)
        if best is None:
            log(f'{name}: no start kept B0 >= {args.min_b0} mT')
            continue
        e, (tr, z, kind, psi), after = best
        typ, p90 = shim_error_spread(sl, tr, z, kind, psi, after, Pchk, args.sigma_br/100, args.sigma_ang)
        log(f'{name}: {ppm(Bchk):.0f} -> {e:.0f} ppm with {len(z)} cubes ({np.sum(kind == 0)} flat, {np.sum(kind > 0)} out, '
            f'{np.sum(kind < 0)} in) in {len(np.unique(tr))} trays; LP bound {lp:.0f} ppm; B0 {after.mean():.2f} mT; '
            f'with shim-cube errors {typ:.0f} typical / {p90:.0f} p90  ({dt:.0f} s)')
        pos, m = sl.pos(tr, z), sl.mdir(tr, kind, psi)
        pd.DataFrame(dict(tray=tr, tray_angle_deg=np.degrees(sl.th[tr]), z=z,
                          orientation=np.where(kind == 0, 'flat', np.where(kind > 0, 'radial_out', 'radial_in')),
                          psi_deg=np.where(kind == 0, np.degrees(psi), 0.0), x=pos[:, 0], y=pos[:, 1],
                          mx=m[:, 0], my=m[:, 1], mz=m[:, 2])).round(4).to_csv(f'{args.out}_{name}.csv', index=False)
        np.savez(f'{args.out}_{name}.npz', trs=tr, z=z, kind=kind, psi=psi, rtr=sl.r, a=sl.a, ntr=sl.n,
                 before=ppm(Bchk), after=e, ledger=args.ledger)
    log(f'-> {args.out}_*.csv')


if __name__ == '__main__':
    main()
