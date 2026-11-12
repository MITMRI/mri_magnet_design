"""
Passive shim insert for a 12 cm design: a stack of printed rings that slides into the bore, one ring per main plate
(boundaries halfway between plates, the end rings out to the magnet ends), each holding a few shim cubes (default
1/8" N42, K&J B222) lying flat in pockets open at the ring's +z (top) face and closed by the next ring -
magnetization in the ring plane at any angle, at any angle around the ring. The stack bolts to the magnet's end
cover. Same pipeline as the main array (field_design_12cm.py):
continuous optimum -> discrete magnets -> refinement.

1. As-built builds: every block of the design gets its own remanence error (uniform within +-br_tol),
   magnetization tilt (uniform on the cap within ang_tol, random azimuth) and position error (+-pos_tol per axis); exact cuboid fields of all blocks, B0 along
   x. Homogeneity = peak-to-peak Bx over points on the DSV surface / mean. (With a measured map this step is replaced
   by the map.)
2. Continuous optimum (Gurobi barrier): magnetization density on every ring plane, `grid` angular positions x 8
   in-plane directions, in cubes, with the packing limit per position. For each target ppm, the least total magnet
   that reaches it. Minimizing the total (not just the ppm) keeps the density sparse, so it discretizes well; a ppm
   minimum with a large budget fills the rings with cancelling density that does not.
3. Discretization: walking around each ring, one cube wherever the accumulated |density| passes k + offset, pointing
   along the local density; several offsets are tried.
4. Refinement: SLSQP minimax over every cube's angle around its ring and its magnetization angle, exact cuboid fields
   (finite-difference Jacobian), neighbours on a ring kept apart so any rotation fits; best offset on a denser check
   point set.
5. Shim-cube errors: the cubes' own Br, angle and position tolerances, exact fields over random draws: the as-built result of one installation.
6. Loads on the best layout (pole-face surface-charge integral, as forces_12cm.py): force and torque on every cube in
   the field of the as-built blocks plus the other cubes; net force and torque per ring and on the whole stack (what
   the bolts carry); axial force while the stack slides in from +z; reverse field on each cube.

Why small cubes: a cube at r ~ 96 mm is ~35 mm from the 120 mm DSV surface. One 1/4" cube puts ~1.4 mT (~14,000 ppm)
on the nearest point, far coarser than a typical as-built error (~0.2-0.4 mT), and its own 3% error is worth
~400 ppm; a 1/8" cube has 1/8 of both. 1/4" cubes in full rings with only the angles designed reached ~1,000 ppm as
built with 216 cubes; 1/8" cubes placed by this pipeline reach ~300-360 ppm with 20-50.

    python utilities/shim_12cm.py designs/12cm_bore200/placements.csv --targets 100 50
    python utilities/shim_12cm.py designs/12cm_bore200/placements.csv --map map.csv --targets 100 50
(--map: x, y, z in mm and Bx in mT, e.g. a few hundred points on the 120 mm sphere surface; the extremes of a
harmonic field lie on the surface, so that is enough.)
Builds run in parallel (--parallel, ~1.5-2 GB each at peak). Writes <design>/shim_<cube>mm/ (shim_<cube>mm_map/ with
--map): shim.log, summary.csv (per build and target), rings.csv (z span, height and cube plane of every ring, cubes per
ring), b<i>_layout.csv (ring, z, phi, magnetization angle psi of every cube), b<i>_build.npz (simulated block errors),
b<i>_cubes.csv (loads per cube), b<i>_insertion.csv (axial force vs insertion). Pictures: plot_shim_12cm.py.
"""
import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor

import gurobipy as gp
import magpylib as magpy
import numpy as np
import pandas as pd
from gurobipy import GRB
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation as R

MU0 = 4e-7 * np.pi
BLOCK = np.array([12.7, 12.7, 25.4])                       # main blocks, mm; local x = magnetization
T0 = time.time()


def ppm(B):
    return 1e6 * np.ptp(B) / B.mean()


def sphere_points(n, r):
    i = np.arange(n) + 0.5
    phi, th = np.arccos(1 - 2 * i / n), np.pi * (1 + 5 ** 0.5) * i
    return r * np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], axis=1)


def cuboid_B(pos, quat, mag, dim, P, rows=1_500_000):
    """Total B (mT) at points P (m,3) of cuboids (pos (n,3) mm, quat (n,4), local magnetization mag (n,3) mT,
    dimension dim (3,)), chunked over sources to bound memory."""
    B = np.zeros((len(P), 3))
    step = max(1, rows // len(P))
    for s in range(0, len(pos), step):
        k = min(step, len(pos) - s)
        B += magpy.getB('Cuboid', observers=np.tile(P, (k, 1)), position=np.repeat(pos[s:s + k], len(P), 0),
                        orientation=R.from_quat(np.repeat(quat[s:s + k], len(P), 0)),
                        magnetization=np.repeat(mag[s:s + k], len(P), 0),
                        dimension=dim).reshape(k, len(P), 3).sum(0)
    return B


def each_bx(pos, psi, br, a, P, rows=1_500_000):
    """Bx (mT) at P of each cube separately: cubes at pos (n,3) lying flat, magnetized in-plane at angle psi -> (n, m)"""
    n, m = len(pos), len(P)
    out = np.empty((n, m))
    step = max(1, rows // m)
    for s in range(0, n, step):
        k = min(step, n - s)
        out[s:s + k] = magpy.getB('Cuboid', observers=np.tile(P, (k, 1)), position=np.repeat(pos[s:s + k], m, 0),
                                  orientation=R.from_euler('z', np.repeat(psi[s:s + k], m)[:, None]),
                                  magnetization=(br, 0, 0), dimension=(a, a, a))[:, 0].reshape(k, m)
    return out


def random_magnetization(n, br, br_tol, ang_tol, rng):
    """Local magnetization vectors (n,3), mT: nominal along local x, Br scaled by 1+u, tilted within ang_tol."""
    u = rng.uniform(-br_tol, br_tol, n)
    ca = rng.uniform(np.cos(np.deg2rad(ang_tol)), 1, n)
    sa, b = np.sqrt(1 - ca ** 2), rng.uniform(0, 2 * np.pi, n)
    return br * (1 + u)[:, None] * np.c_[ca, sa * np.cos(b), sa * np.sin(b)]


def as_built(csv, br, br_tol, ang_tol, pos_tol, P_list, seed, builds):
    """Exact Bx (mT) of `builds` random as-built versions of the design on each point set in P_list (Br, angle and
    position errors; position uniform within +-pos_tol mm on each axis)."""
    df = pd.read_csv(csv)
    pos = df[['x', 'y', 'z']].to_numpy()
    quat = R.from_euler('z', df.rot_z_deg.to_numpy()[:, None], degrees=True).as_quat()
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(builds):
        mag = random_magnetization(len(df), br, br_tol / 100, ang_tol, rng)
        p = pos + rng.uniform(-pos_tol, pos_tol, pos.shape)
        out.append(([cuboid_B(p, quat, mag, BLOCK, P)[:, 0] for P in P_list], mag, p))
    return df, pos, quat, out


class Stack:
    """Ring geometry: cube centres on radius r at the ring planes z; neighbours at least `gap` apart (centre to
    centre, so any rotation of either cube fits)."""

    def __init__(self, cube, shim_br, ring_z, bore, fit=0.2, ring_wall=1.0, pocket_wall=1.0):
        self.a, self.br = cube, shim_br
        self.z = np.asarray(ring_z, float)
        self.od = bore - 2 * fit
        half_diag = cube / np.sqrt(2)
        self.r = self.od / 2 - ring_wall - half_diag
        self.clear = 2 * (self.r - half_diag - ring_wall)
        self.gap = 2 * half_diag + pocket_wall
        self.dphi = 2 * np.arcsin(self.gap / (2 * self.r))
        self.per_ring = int(2 * np.pi / self.dphi)

    def pos(self, ring, phi):
        return np.c_[self.r * np.cos(phi), self.r * np.sin(phi), self.z[ring]]


# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
# placement: continuous optimum -> cubes -> refinement

class Density:
    """Continuous optimum: magnetization density on every ring plane, `grid` positions x 8 in-plane directions, in
    units of cubes, at most the packing limit per grid point and max_cubes in total. The model is built once; solve()
    returns, for a target ppm, the least total magnet that reaches it (or the best reachable)."""

    def __init__(self, st, B0, P, grid, max_cubes, threads):
        self.phi = 2 * np.pi * np.arange(grid) / grid
        self.nr, self.nd, self.grid = len(st.z), 8, grid
        self.dirs = 2 * np.pi * np.arange(self.nd) / self.nd                       # relative to radial
        ring = np.repeat(np.arange(self.nr), grid * self.nd)
        ph = np.tile(np.repeat(self.phi, self.nd), self.nr)
        A = each_bx(st.pos(ring, ph), ph + np.tile(self.dirs, self.nr * grid), st.br, st.a, P, rows=400_000)
        Ac = (A - A.mean(1, keepdims=True)).T
        del A
        d0 = B0 - B0.mean()
        self.scale = B0.mean() / 2e6
        m = gp.Model()
        m.Params.OutputFlag, m.Params.Threads, m.Params.Method, m.Params.Crossover = 0, threads, 2, 0
        self.x = m.addMVar(Ac.shape[1], lb=0)
        self.t = m.addVar(lb=0)
        m.addConstr(Ac @ self.x - self.t <= -d0)
        m.addConstr(-Ac @ self.x - self.t <= d0)
        m.addConstr(self.x.reshape(self.nr * grid, self.nd).sum(1) <= (2 * np.pi * st.r / grid) / st.gap)
        m.addConstr(self.x.sum() <= max_cubes)
        m.setObjective(self.t, GRB.MINIMIZE)
        m.optimize()
        self.m, self.best = m, self.t.X / self.scale
        m.setObjective(self.x.sum(), GRB.MINIMIZE)

    def solve(self, target):
        self.t.UB = max(target, self.best * 1.0001) * self.scale
        self.m.optimize()
        X = self.x.X.reshape(self.nr, self.grid, self.nd)
        return self.t.X / self.scale, X.sum(), self.phi, X @ np.c_[np.cos(self.dirs), np.sin(self.dirs)]


def place(st, phi, V, offset, max_cubes):
    """Cubes from the density: walking around each ring, one cube wherever the accumulated |density| passes
    k + offset, pointing along the local density (summed over +-3 grid steps); close neighbours are pushed apart."""
    grid = len(phi)
    rings, ph, ps, wt = [], [], [], []
    for k in range(len(st.z)):
        v = V[k]
        cum = np.cumsum(np.linalg.norm(v, axis=1))
        sel = []
        for c in range(int(np.floor(cum[-1] + 1 - offset))):
            j = min(np.searchsorted(cum, c + offset), grid - 1)
            w = v[np.arange(j - 3, j + 4) % grid].sum(0)
            if np.linalg.norm(w) > 1e-12:
                sel.append((phi[j], phi[j] + np.arctan2(w[1], w[0]), np.linalg.norm(w)))
        sel = sorted(sorted(sel, key=lambda s_: -s_[2])[:st.per_ring])
        if not sel:
            continue
        p = np.array([s_[0] for s_ in sel])
        for _ in range(500):
            short = np.maximum(st.dphi - np.diff(np.r_[p, p[0] + 2 * np.pi]), 0)
            if len(p) == 1 or short.max() < 1e-6:
                break
            p = p - short / 2 + np.roll(short, 1) / 2
        for (_, s_psi, s_w), pp in zip(sel, p):
            rings.append(k); ph.append(np.mod(pp, 2 * np.pi)); ps.append(s_psi); wt.append(s_w)
    rings, ph, ps, wt = (np.array(v_) for v_ in (rings, ph, ps, wt))
    keep = np.sort(np.argsort(-wt)[:max_cubes])
    return rings[keep].astype(int), ph[keep], ps[keep]


def refine_cubes(st, B0, P, ring, phi, psi, iters):
    """SLSQP minimax ppm over every cube's angle around its ring (phi) and magnetization angle (psi), exact cuboid
    fields with a finite-difference Jacobian; neighbours on a ring stay at least st.dphi apart."""
    n = len(phi)
    if n == 0:
        return ring, phi, psi
    o = np.lexsort((phi, ring))
    ring, phi, psi = ring[o], phi[o], psi[o]
    nxt = []
    for k in np.unique(ring):
        idx = np.flatnonzero(ring == k)
        if len(idx) > 1:
            nxt += [(i, j, 0.0) for i, j in zip(idx[:-1], idx[1:])] + [(idx[-1], idx[0], 2 * np.pi)]
    h = 1e-5

    def fields(ph, ps):
        return each_bx(st.pos(ring, ph), ps, st.br, st.a, P)

    def cons(x):
        ph, ps = x[:n], x[n:2 * n]
        Bt = B0 + fields(ph, ps).sum(0)
        mu = Bt.mean()
        dev = 1e6 * (Bt - mu) / mu
        return np.r_[x[-1] - dev, x[-1] + dev, [ph[j] + w - ph[i] - st.dphi for i, j, w in nxt]]

    def cjac(x):
        ph, ps = x[:n], x[n:2 * n]
        Bi = fields(ph, ps)
        Bt = B0 + Bi.sum(0)
        mu = Bt.mean()
        J = np.vstack([(fields(ph + h, ps) - Bi) / h, (fields(ph, ps + h) - Bi) / h])
        dd = 1e6 * (J * mu - Bt[None] * J.mean(1)[:, None]) / mu ** 2
        npt = len(P)
        top = np.zeros((2 * npt, 2 * n + 1))
        top[:npt, :-1], top[npt:, :-1], top[:, -1] = -dd.T, dd.T, 1
        g = np.zeros((len(nxt), 2 * n + 1))
        for r_, (i, j, _) in enumerate(nxt):
            g[r_, i], g[r_, j] = -1, 1
        return np.vstack([top, g])

    x0 = np.r_[phi, psi, 0.5 * ppm(B0 + fields(phi, psi).sum(0)) * 1.05]
    res = minimize(lambda x: x[-1], x0, jac=lambda x: np.eye(len(x))[-1], method='SLSQP',
                   constraints=[{'type': 'ineq', 'fun': cons, 'jac': cjac}], options={'maxiter': iters, 'ftol': 1e-7})
    return ring, np.mod(res.x[:n], 2 * np.pi), np.mod(res.x[n:2 * n], 2 * np.pi)


def design_density(st, dens, Bopt, Popt, Bchk, Pchk, target, max_cubes, starts, iters, log=print):
    """Continuous optimum for the target -> cubes (several walking offsets) -> refinement; best on the check points.
    Returns (continuous ppm, continuous cubes, (ring, phi, psi, Bx after))."""
    lp, used, phi_g, V = dens.solve(target)
    log(f'    target {target:g} ppm: continuous optimum {lp:.0f} ppm with {used:.1f} cubes of magnet')
    best = None
    for off in np.linspace(0.5, 0.0, starts, endpoint=False):
        ring, ph, ps = place(st, phi_g, V, off, max_cubes)
        ring, ph, ps = refine_cubes(st, Bopt, Popt, ring, ph, ps, iters)
        after = Bchk + (each_bx(st.pos(ring, ph), ps, st.br, st.a, Pchk).sum(0) if len(ph) else 0)
        if best is None or ppm(after) < ppm(best[3]):
            best = (ring, ph, ps, after)
    log(f'      best of {starts} offsets: {len(best[1])} cubes -> {ppm(best[3]):.0f} ppm')
    return lp, used, best


def run_build(job):
    """One as-built build: continuous model once, then every target; returns rows and the best layout."""
    a, b, Bopt, Bchk, st, Popt, Pchk = job
    rng = np.random.default_rng(a.seed + 1000 + b)
    lines = []

    def log(s):
        lines.append(s)
        print(f'[{time.time() - T0:6.0f} s] build {b}: {s.strip()}', flush=True)

    dens = Density(st, Bopt, Popt, a.grid, a.max_cubes, a.threads)
    log(f'continuous model ready; best reachable {dens.best:.0f} ppm')
    rows, best = [], None
    for target in a.targets:
        lp, used, (ring, phi, psi, after) = design_density(st, dens, Bopt, Popt, Bchk, Pchk, target, a.max_cubes,
                                                          a.starts, a.iters, log)
        med, p90 = cube_error_spread(st, ring, phi, psi, Bchk, Pchk, a.shim_br_tol, a.shim_ang_tol, a.shim_pos_tol,
                                     a.error_draws, rng) if len(phi) else (ppm(Bchk), ppm(Bchk))
        rows.append(dict(build=b, target=target, continuous_ppm=lp, continuous_cubes=used, cubes=len(phi),
                         unshimmed=ppm(Bchk), perfect_cubes=ppm(after), real_median=med, real_p90=p90,
                         mass_g=len(phi) * a.mass, cost=len(phi) * a.price))
        log(f'target {target:g}: {len(phi)} cubes -> {ppm(after):.0f} ppm perfect, {med:.0f} median / {p90:.0f} p90 '
            f'real (one installation)')
        if best is None or med < best[0]:
            best = (med, ring, phi, psi)
    return rows, best, lines


def cube_error_spread(st, ring, phi, psi, Bbase, P, br_tol, ang_tol, pos_tol, draws, rng):
    """median / p90 ppm when every shim cube has its own Br, angle and position error (position uniform within
    +-pos_tol mm on each axis), exact fields per draw"""
    spos, quat = st.pos(ring, phi), R.from_euler('z', psi[:, None]).as_quat()
    n, m = len(spos), len(P)
    ori = R.from_quat(np.repeat(quat, m, 0))
    e = []
    for _ in range(draws):
        p = spos + rng.uniform(-pos_tol, pos_tol, spos.shape)
        mag = random_magnetization(n, st.br, br_tol / 100, ang_tol, rng)
        B = magpy.getB('Cuboid', observers=np.tile(P, (n, 1)), position=np.repeat(p, m, 0), orientation=ori,
                       magnetization=np.repeat(mag, m, 0), dimension=(st.a,) * 3)[:, 0].reshape(n, m)
        e.append(ppm(Bbase + B.sum(0)))
    return np.median(e), np.percentile(e, 90)


# ----------------------------------------------------------------------------------------------------------------
# loads

def face_quadrature(a, nq):
    """Gauss points on the 6 faces of a cube of side a (local), outward normals, weights (mm^2)."""
    g, w = np.polynomial.legendre.leggauss(nq)
    U, V = np.meshgrid(g * a / 2, g * a / 2, indexing='ij')
    W = np.outer(w, w).ravel() * (a / 2) ** 2
    pts, nrm, wts = [], [], []
    for ax in range(3):
        o = [i for i in range(3) if i != ax]
        for s in (1, -1):
            p = np.zeros((U.size, 3)); p[:, ax] = s * (a / 2 + 1e-3); p[:, o[0]] = U.ravel(); p[:, o[1]] = V.ravel()
            n = np.zeros(3); n[ax] = s
            pts.append(p); nrm.append(np.tile(n, (U.size, 1))); wts.append(W)
    return np.vstack(pts), np.vstack(nrm), np.concatenate(wts)


def cube_loads(pos, quat, mag, a, field_fn, nq=4):
    """Force (N, global) and torque about the cube centre (N m) on cubes (pos, quat, local mag mT) from the external
    field field_fn(points, cube index) -> B (mT); only faces carrying charge are integrated."""
    lp, ln, lw = face_quadrature(a, nq)
    rot = R.from_quat(quat)
    sig = mag @ ln.T
    idx, q = np.nonzero(np.abs(sig) > 1e-9)
    P = rot[idx].apply(lp[q]) + pos[idx]
    B = field_fn(P, idx)
    dF = ((sig[idx, q] * 1e-3 / MU0) * lw[q] * 1e-6)[:, None] * B * 1e-3
    F = np.zeros((len(pos), 3)); np.add.at(F, idx, dF)
    Tq = np.zeros((len(pos), 3)); np.add.at(Tq, idx, np.cross((P - pos[idx]) * 1e-3, dF))
    return F, Tq


_W = {}


def _init_worker(state):
    _W.update(state)


def _insertion_row(dz):
    """axial force (N) on the stack and on each ring with every cube shifted dz along +z, blocks only"""
    w = _W
    Fd, _ = cube_loads(w['spos'] + [0, 0, dz], w['squat'], w['smag'], w['a'],
                       lambda P, idx: cuboid_B(w['bpos'], w['bquat'], w['bmag'], BLOCK, P, rows=w['rows']))
    return np.r_[dz, Fd[:, 2].sum(), [Fd[w['ring'] == k, 2].sum() for k in w['rings']]]


def loads(st, ring, phi, psi, bpos, bquat, bmag, zmax_blocks, jobs, out, b, log):
    spos = st.pos(ring, phi)
    squat = R.from_euler('z', psi[:, None]).as_quat()
    smag = np.tile([st.br, 0, 0], (len(spos), 1))

    def main_field(P, idx=None):
        return cuboid_B(bpos, bquat, bmag, BLOCK, P)

    def total_field(P, idx):
        Bs = cuboid_B(spos, squat, smag, (st.a,) * 3, P)
        self_ = magpy.getB('Cuboid', observers=P, position=spos[idx], orientation=R.from_quat(squat[idx]),
                           magnetization=smag[idx], dimension=(st.a,) * 3)
        return main_field(P) + Bs - self_

    F, Tq = cube_loads(spos, squat, smag, st.a, total_field)
    Fm, Tm = cube_loads(spos, squat, smag, st.a, lambda P, idx: main_field(P))
    Bc = main_field(spos) + cuboid_B(spos, squat, smag, (st.a,) * 3, spos) - \
        magpy.getB('Cuboid', observers=spos, position=spos, orientation=R.from_quat(squat), magnetization=smag,
                   dimension=(st.a,) * 3)
    mhat = np.c_[np.cos(psi), np.sin(psi), 0 * psi]
    rev = np.maximum(-(Bc * mhat).sum(1), 0)
    er = np.c_[np.cos(phi), np.sin(phi), 0 * phi]
    et = np.c_[-np.sin(phi), np.cos(phi), 0 * phi]
    cubes = pd.DataFrame(dict(ring=ring, z=st.z[ring], phi_deg=np.degrees(phi), psi_deg=np.degrees(psi),
                              x=spos[:, 0], y=spos[:, 1], Fr=(F * er).sum(1), Ft=(F * et).sum(1), Fz=F[:, 2],
                              F=np.linalg.norm(F, axis=1), torque_Nmm=1e3 * np.linalg.norm(Tq, axis=1),
                              Tz_Nmm=1e3 * Tq[:, 2], B_local_mT=np.linalg.norm(Bc, axis=1), reverse_mT=rev))
    cubes.round(4).to_csv(os.path.join(out, f'b{b}_cubes.csv'), index=False)
    i, j = cubes.F.idxmax(), cubes.torque_Nmm.idxmax()
    log(f'build {b} loads (as-built blocks + other cubes):')
    log(f'    per cube: max force {cubes.F[i]:.3f} N (ring z {cubes.z[i]:+.0f}); max |Fr| {cubes.Fr.abs().max():.3f}, '
        f'|Ft| {cubes.Ft.abs().max():.3f}, |Fz| {cubes.Fz.abs().max():.3f} N; max torque {cubes.torque_Nmm[j]:.1f} '
        f'N mm (field {cubes.B_local_mT[j]:.0f} mT); field at the cubes {cubes.B_local_mT.min():.0f}-'
        f'{cubes.B_local_mT.max():.0f} mT')
    log(f'    worst reverse field on a cube {rev.max():.0f} mT = {rev.max() * 1e-3 / MU0 / 1e3:.0f} kA/m '
        f'(N42 intrinsic coercivity >= 955 kA/m)')
    per = pd.DataFrame(dict(ring=ring, Fx=Fm[:, 0], Fy=Fm[:, 1], Fz=Fm[:, 2],
                            Tz=np.cross(spos * 1e-3, Fm)[:, 2] + Tm[:, 2])).groupby('ring').agg(
        n=('Fx', 'size'), Fx=('Fx', 'sum'), Fy=('Fy', 'sum'), Fz=('Fz', 'sum'), Tz=('Tz', 'sum'))
    for k, rw in per.iterrows():
        log(f'    ring z {st.z[k]:+5.0f} ({int(rw.n):2d} cubes): net F ({rw.Fx:+.2f}, {rw.Fy:+.2f}, {rw.Fz:+.2f}) N, '
            f'torque about the bore axis {rw.Tz:+.3f} N m')
    net = Fm.sum(0)
    Tnet = (np.cross(spos * 1e-3, Fm) + Tm).sum(0)
    log(f'    whole stack (what the bolts carry): F ({net[0]:+.2f}, {net[1]:+.2f}, {net[2]:+.2f}) N, '
        f'torque ({Tnet[0]:+.3f}, {Tnet[1]:+.3f}, {Tnet[2]:+.3f}) N m')

    span = st.z.max() + st.a / 2
    dzs = np.r_[np.arange(0, 20, 2), np.arange(20, 100, 5), np.arange(100, zmax_blocks + 12.7 + 2 * span + 40, 10)]
    state = dict(spos=spos, squat=squat, smag=smag, a=st.a, bpos=bpos, bquat=bquat, bmag=bmag, ring=ring,
                 rings=np.unique(ring), rows=600_000)
    with ProcessPoolExecutor(jobs, initializer=_init_worker, initargs=(state,)) as ex:
        rows = list(ex.map(_insertion_row, dzs))
    ins = pd.DataFrame(rows, columns=['dz', 'stack_Fz'] + [f'ring{k}_Fz' for k in np.unique(ring)])
    ins.round(4).to_csv(os.path.join(out, f'b{b}_insertion.csv'), index=False)
    k = ins.stack_Fz.abs().idxmax()
    log(f'    sliding the stack in from +z: peak axial force {ins.stack_Fz[k]:+.2f} N at {ins.dz[k]:.0f} mm to go '
        f'(+ = pushed back out); seated {ins.stack_Fz[0]:+.2f} N')


# ----------------------------------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('csv')
    ap.add_argument('--builds', type=int, default=3)
    ap.add_argument('--map', help='measured field map csv: x, y, z (mm), Bx (mT), e.g. points on the DSV surface; '
                                  'replaces the simulated builds')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--br', type=float, default=1300, help='main block remanence, mT')
    ap.add_argument('--br-tol', type=float, default=3.0, help='+-%% remanence error of the main blocks')
    ap.add_argument('--ang-tol', type=float, default=1.5, help='+-deg magnetization angle error of the main blocks')
    ap.add_argument('--pos-tol', type=float, default=0.1, help='+-mm position error of the main blocks (each axis)')
    ap.add_argument('--dsv', type=float, default=120)
    ap.add_argument('--opt-points', type=int, default=500)
    ap.add_argument('--check-points', type=int, default=2000)
    ap.add_argument('--cube', type=float, default=3.175, help='shim cube side, mm (K&J B222: 1/8")')
    ap.add_argument('--shim-br', type=float, default=1300, help='shim cube remanence, mT (N42)')
    ap.add_argument('--shim-br-tol', type=float, default=3.0, help='+-%%, shim cubes')
    ap.add_argument('--shim-ang-tol', type=float, default=1.5, help='+-deg, shim cubes')
    ap.add_argument('--shim-pos-tol', type=float, default=0.1, help='+-mm position error, shim cubes (each axis)')
    ap.add_argument('--price', type=float, default=0.22, help='$ per shim cube (B222)')
    ap.add_argument('--mass', type=float, default=0.24, help='g per shim cube (B222)')
    ap.add_argument('--top-gap', type=float, default=0.3, help='cube top below the ring top face, mm')
    ap.add_argument('--bore', type=float, default=200, help='clear bore of the magnet, mm (ring OD = bore - 2 fit)')
    ap.add_argument('--targets', type=float, nargs='+', default=[100, 50],
                    help='ppm targets for the continuous optimum (fewest cubes that reach each)')
    ap.add_argument('--max-cubes', type=int, default=200)
    ap.add_argument('--grid', type=int, default=240, help='continuous positions per ring')
    ap.add_argument('--starts', type=int, default=8, help='discretization offsets tried per target')
    ap.add_argument('--iters', type=int, default=300, help='SLSQP iterations')
    ap.add_argument('--error-draws', type=int, default=300)
    ap.add_argument('--no-loads', action='store_true')
    ap.add_argument('--parallel', type=int, default=3, help='builds in parallel (~1.5-2 GB each)')
    ap.add_argument('--threads', type=int, default=4, help='Gurobi threads per build')
    ap.add_argument('--jobs', type=int, default=max(1, min(8, (os.cpu_count() or 2) - 1)),
                    help='worker processes for the insertion sweep (~350 MB each)')
    a = ap.parse_args()

    # one shim ring per main plate: boundaries halfway between plates, ends at the magnet ends; cubes at each ring's
    # +z (top) face, the pocket open to the top and closed by the next ring
    pz = np.sort(pd.read_csv(a.csv).z.round(1).unique())
    bounds = np.r_[pz[0] - BLOCK[2] / 2, (pz[:-1] + pz[1:]) / 2, pz[-1] + BLOCK[2] / 2]
    st = Stack(a.cube, a.shim_br, bounds[1:] - a.cube / 2 - a.top_gap, a.bore)
    out = os.path.join(os.path.dirname(a.csv), f'shim_{a.cube:g}mm' + ('_map' if a.map else ''))
    os.makedirs(out, exist_ok=True)
    logf = open(os.path.join(out, 'shim.log'), 'w')

    def log(s):
        s = f'[{time.time() - T0:6.0f} s] {s}'
        print(s, flush=True); logf.write(s + '\n'); logf.flush()

    if a.map:
        # measured map: x, y, z (mm), Bx (mT; B0 along x), e.g. points on the DSV surface. Bx is harmonic, so a map of
        # the sphere surface holds the extremes of the whole sphere. Optimize on up to opt_points of it, check on all.
        M = pd.read_csv(a.map)
        Pchk, Bmap = M[['x', 'y', 'z']].to_numpy(float), np.abs(M['Bx'].to_numpy(float))
        sub = np.linspace(0, len(Pchk) - 1, min(a.opt_points, len(Pchk))).round().astype(int)
        Popt = Pchk[sub]
        df = pd.read_csv(a.csv)
        bpos = df[['x', 'y', 'z']].to_numpy()
        bquat = R.from_euler('z', df.rot_z_deg.to_numpy()[:, None], degrees=True).as_quat()
        builds = [((Bmap[sub], Bmap), np.tile([a.br, 0, 0], (len(df), 1)), bpos)]   # loads use nominal blocks
        rad = np.linalg.norm(Pchk - Pchk.mean(0), axis=1)
        log(f'{a.map}: {len(Pchk)} measured points (radius {rad.min():.1f}-{rad.max():.1f} mm about '
            f'{np.round(Pchk.mean(0), 1).tolist()}), optimizing on {len(Popt)}; B0 {Bmap.mean():.3f} mT')
    else:
        Popt, Pchk = sphere_points(a.opt_points, a.dsv / 2), sphere_points(a.check_points, a.dsv / 2)
        df, bpos, bquat, builds = as_built(a.csv, a.br, a.br_tol, a.ang_tol, a.pos_tol, [Popt, Pchk], a.seed, a.builds)
        log(f'{a.csv}: {len(df)} blocks; block errors Br +-{a.br_tol}%, angle +-{a.ang_tol} deg, '
            f'position +-{a.pos_tol} mm (uniform)')
    log(f'shim stack: {len(st.z)} ring planes at z = {", ".join(f"{z:.0f}" for z in st.z)} mm; ring OD {st.od:.1f}, '
        f'cube centres at r {st.r:.2f}, clear bore {st.clear:.1f} mm; {a.cube} mm cubes lying flat '
        f'(Br {a.shim_br:.0f} mT, errors +-{a.shim_br_tol}% / +-{a.shim_ang_tol} deg / +-{a.shim_pos_tol} mm), at most {st.per_ring} per ring, '
        f'{a.max_cubes} in total')
    for b, ((Bopt, Bchk), bmag, bp) in enumerate(builds):
        log(f'build {b}: unshimmed {ppm(Bchk):.0f} ppm, B0 {Bchk.mean():.2f} mT')
        if not a.map:
            np.savez(os.path.join(out, f'b{b}_build.npz'), bmag=bmag, bpos=bp)   # the build's block magnetizations, positions
    jobs = [(a, b, Bopt, Bchk, st, Popt, Pchk) for b, ((Bopt, Bchk), _, _) in enumerate(builds)]
    with ProcessPoolExecutor(min(a.parallel, len(jobs))) as ex:
        results = list(ex.map(run_build, jobs))
    rows = [r for res_, _, _ in results for r in res_]
    for b, (_, best, lines) in enumerate(results):
        for line in lines:
            logf.write(f'build {b}: {line}\n')
        pd.DataFrame(dict(ring=best[1], z=st.z[best[1]], phi_deg=np.degrees(best[2]), psi_deg=np.degrees(best[3]))
                     ).round(3).to_csv(os.path.join(out, f'b{b}_layout.csv'), index=False)
        if not a.no_loads:
            loads(st, best[1], best[2], best[3], builds[b][2], bquat, builds[b][1], df.z.max(), a.jobs, out, b, log)
    rt = pd.DataFrame(dict(ring=np.arange(len(st.z)), z_bottom=bounds[:-1], z_top=bounds[1:],
                           height=np.diff(bounds), cube_plane_z=st.z))
    for b, (_, best, _) in enumerate(results):
        rt[f'cubes_b{b}'] = np.bincount(best[1], minlength=len(st.z))
    rt.round(2).to_csv(os.path.join(out, 'rings.csv'), index=False)
    log('rings (one per main plate, cubes at the top face):\n' + rt.round(1).to_string(index=False))
    res = pd.DataFrame(rows)
    res.round(1).to_csv(os.path.join(out, 'summary.csv'), index=False)
    log('summary (median over builds):\n' + res.groupby('target').median(numeric_only=True).drop(columns='build')
        .round(0).to_string())
    log(f'-> {out}/')


if __name__ == '__main__':
    main()
