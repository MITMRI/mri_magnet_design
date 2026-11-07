"""
Maximum-B0 Halbach ring design for a 120 mm imaging sphere with N42SH 25.4 x 12.7 x 12.7 mm blocks
(KJ Magnetics BX088SH, magnetized through a 12.7 mm face).

Each block lies with its 25.4 mm length along z; its magnetization lies in the ring plane at angle 2*theta
(k = 1 Halbach). Rings sit in slices (shared z planes, mirrored at -z); a slice holds one or more rings of
different radius. The design is a list of slices, each a list of rings (radius, magnet count).

Pipeline:
  1. Sensitivity: exact Bx per magnet of a full ring, for every candidate (radius, z) on a grid.
  2. MILP: pick rings and integer magnet counts that maximize B0 subject to |Bx - B0| <= eps*B0 on the
     sphere surface (Bx is harmonic, so its extremes lie on the surface), rings in a slice not overlapping radially,
     slices not overlapping axially.
  3. Gradient refinement (SLSQP): move slice z and ring radii continuously, same objective and constraints.
     Gradient from a point-dipole model (3 dipoles per block); the exact-minus-dipole field is added as a
     correction, so constraint values are exact.
  4. Greedy +1 magnet per ring, re-refined, kept if B0 rises.
  5. Check on 20,000 points on the sphere surface.

    python field_design_12cm.py results/12cm/bore180.xlsx --bore 180
"""

import argparse
import json
import time

import magpylib as magpy
import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, milp, minimize
from scipy.sparse import coo_matrix
from scipy.spatial.transform import Rotation

T0 = time.time()


def log(msg):
    print(f'[{time.time() - T0:6.1f} s] {msg}', flush=True)


def ppm(B):
    return 1e6 * (B.max() - B.min()) / B.mean()


def sphere_points(n, r):
    i = np.arange(n) + 0.5
    phi, th = np.arccos(1 - 2 * i / n), np.pi * (1 + 5 ** 0.5) * i
    return r * np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], axis=1)


class Geometry:
    def __init__(self, a):
        self.a, self.h, self.br = a.a, a.h, a.br
        half_diag = a.a / np.sqrt(2)                                  # block footprint radius in the ring plane
        self.pitch = a.h + a.slice_wall                               # min z spacing of slices
        self.dr = 2 * half_diag + a.layer_wall                        # min radial spacing of rings in a slice
        self.rmin = a.bore / 2 + a.bore_wall + half_diag
        self.rmax = a.max_od / 2 - a.od_wall - half_diag
        self.zmax = (a.max_length - a.h) / 2
        self.fit = half_diag + a.ring_wall / 2                        # n blocks fit at r if pi/n >= asin(fit/r)
        self.rs, self.ripple = a.sphere_d / 2, a.ripple

    def nmax(self, r):
        return int(np.pi / np.arcsin(self.fit / r) + 1e-9)

    def nmin(self, r):
        """Fewest blocks in a ring at r: an n-block ring's field ripple on the sphere scales as (rs/r)^(n-1)."""
        return int(np.ceil(1 + np.log(self.ripple) / np.log(self.rs / r)))


def ring_positions(r, z, n):
    th = 2 * np.pi * np.arange(n) / n
    zs = [z, -z] if z > 1e-9 else [0.0]
    pos = np.vstack([np.stack([r * np.cos(th), r * np.sin(th), np.full(n, zz)], 1) for zz in zs])
    return pos, np.tile(2 * th, len(zs))


def ring_field(geo, r, z, n, sensors):
    """Exact Bx (mT) of a ring (and its mirror at -z) on the sensors."""
    pos, rot = ring_positions(r, z, n)
    nm, ns = len(pos), len(sensors)
    B = magpy.getB('Cuboid', observers=np.tile(sensors, (nm, 1)), position=np.repeat(pos, ns, 0),
                   orientation=Rotation.from_euler('z', np.repeat(rot, ns)[:, None]), magnetization=(geo.br, 0, 0),
                   dimension=(geo.a, geo.a, geo.h))
    return B[:, 0].reshape(nm, ns).sum(0)


# --------------------------------------------------------------------------------------------------------------
# MILP: x = [n (L*Z ints), y (L*Z ring used), s (Z slice used), b]

def ring_table(geo, r, n, zs, sensors):
    """Bx per magnet (mT) on the sensors of a ring of n blocks at radius r, for every z in zs: (len(zs), S)."""
    pos, rot = zip(*[ring_positions(r, z, n) for z in zs])
    sizes = [len(p) for p in pos]
    pos, rot = np.vstack(pos), np.concatenate(rot)
    nm, ns = len(pos), len(sensors)
    B = magpy.getB('Cuboid', observers=np.tile(sensors, (nm, 1)), position=np.repeat(pos, ns, 0),
                   orientation=Rotation.from_euler('z', np.repeat(rot, ns)[:, None]), magnetization=(geo.br, 0, 0),
                   dimension=(geo.a, geo.a, geo.h))[:, 0].reshape(nm, ns)
    return np.array([blk.sum(0) for blk in np.split(B, np.cumsum(sizes)[:-1])]) / n


def milp_design(geo, sensors, eps, dr_grid, dz, time_limit):
    radii = np.arange(geo.rmax, geo.rmin - 1e-9, -dr_grid)[::-1]        # grid anchored at the OD limit
    radii = np.array([r for r in radii if geo.nmin(r) <= geo.nmax(r)])
    zs = np.concatenate([[0.0], np.arange(geo.pitch / 2, geo.zmax + 1e-9, dz)])
    L, Z, S = len(radii), len(zs), len(sensors)
    nmax = np.array([geo.nmax(r) for r in radii])
    nmin = np.array([geo.nmin(r) for r in radii])
    A = np.array([ring_table(geo, r, n, zs, sensors) for r, n in zip(radii, nmax)])     # (L, Z, S)
    log(f'MILP: {L} radii ({radii[0]:.0f}-{radii[-1]:.0f} mm) x {Z} z positions, {S} check points')

    N, Y, Sl, b = 0, L * Z, 2 * L * Z, 2 * L * Z + Z
    nv = b + 1
    q = lambda l, j: l * Z + j
    ri, ci, vals, lo, hi = [], [], [], [], []                           # sparse constraint matrix

    def add(entries, l_, h_):
        for i, c in entries:
            ri.append(len(lo)); ci.append(i); vals.append(c)
        lo.append(l_); hi.append(h_)

    M = A.reshape(L * Z, S).T
    for k in range(S):                                                  # (1-eps) b <= B_k <= (1+eps) b
        add(list(enumerate(M[k])) + [(b, -(1 + eps))], -np.inf, 0)
        add(list(enumerate(M[k])) + [(b, -(1 - eps))], 0, np.inf)
    for l in range(L):
        for j in range(Z):
            add([(N + q(l, j), 1), (Y + q(l, j), -nmax[l])], -np.inf, 0)  # n <= nmax * y
            add([(N + q(l, j), 1), (Y + q(l, j), -nmin[l])], 0, np.inf)   # n >= nmin * y (no ripple)
            add([(Y + q(l, j), 1), (Sl + j, -1)], -np.inf, 0)             # ring only in a used slice
        win = [m for m in range(l, L) if radii[m] < radii[l] + geo.dr]
        for j in range(Z):                                              # rings in a slice do not overlap
            if len(win) > 1:
                add([(Y + q(m, j), 1) for m in win], -np.inf, 1)
    for j in range(Z):                                                  # slices do not overlap
        win = np.where((zs >= zs[j]) & (zs < zs[j] + geo.pitch))[0]
        if len(win) > 1:
            add([(Sl + w, 1) for w in win], -np.inf, 1)

    c = np.zeros(nv); c[b] = -1
    res = milp(c, constraints=LinearConstraint(coo_matrix((vals, (ri, ci)), shape=(len(lo), nv)), lo, hi),
               integrality=np.r_[np.ones(nv - 1), 0],
               bounds=Bounds(0, np.r_[np.repeat(nmax, Z), np.ones(L * Z + Z), np.inf]),
               options={'time_limit': time_limit, 'mip_rel_gap': 1e-3})
    if res.x is None:
        raise SystemExit(f'MILP found no solution: {res.message}')
    n = np.round(res.x[:L * Z]).astype(int).reshape(L, Z)
    rings = [dict(r=radii[l], z=zs[j], n=n[l, j]) for j in range(Z) for l in range(L) if n[l, j] > 0]
    if not rings:
        raise SystemExit(f'MILP returned no rings: {res.message}')
    if res.mip_dual_bound is not None:
        log(f'MILP upper bound {-res.mip_dual_bound:.1f} mT')
    return rings, res.x[b], res.message


# --------------------------------------------------------------------------------------------------------------
# design: slices at z (a z = 0 slice if present, the others mirrored at -z), each holding rings (radius, count)

class Design:
    def __init__(self, geo, rings):
        self.geo = geo
        self.z = sorted({float(rg['z']) for rg in rings})
        self.center = self.z[0] == 0
        self.rings = [dict(k=self.z.index(float(rg['z'])), r=float(rg['r']), n=int(rg['n']))
                      for rg in sorted(rings, key=lambda g: (g['z'], g['r']))]

    def copy(self):
        d = Design.__new__(Design)
        d.geo, d.z, d.center, d.rings = self.geo, list(self.z), self.center, [dict(rg) for rg in self.rings]
        return d

    def count(self):
        return sum(rg['n'] * (2 if self.z[rg['k']] > 1e-9 else 1) for rg in self.rings)

    # variables: z of each slice (a z = 0 slice stays there), then radius of each ring
    def pack(self):
        return np.array(self.z + [rg['r'] for rg in self.rings])

    def unpack(self, x):
        self.z = list(x[:len(self.z)])
        for rg, r in zip(self.rings, x[len(self.z):]):
            rg['r'] = r

    def field(self, sensors):
        return sum(ring_field(self.geo, rg['r'], self.z[rg['k']], rg['n'], sensors) for rg in self.rings)

    def field_dipole(self, sensors):
        """Bx and its Jacobian wrt (slice z, ring r), each block as 3 point dipoles along its length."""
        nz = len(self.z)
        k = self.geo.br * self.geo.a ** 2 * self.geo.h / 3 / (4 * np.pi)
        B, J = np.zeros(len(sensors)), np.zeros((len(sensors), nz + len(self.rings)))
        for i, rg in enumerate(self.rings):
            pos, rot = ring_positions(rg['r'], self.z[rg['k']], rg['n'])
            sign = np.sign(pos[:, 2])[:, None]
            th = (rot / 2)[:, None]
            mx, my = np.cos(rot)[:, None], np.sin(rot)[:, None]
            for off in (-self.geo.h / 3, 0, self.geo.h / 3):
                d = sensors[None] - (pos + [0, 0, off])[:, None]
                r2 = (d ** 2).sum(-1); ir5 = r2 ** -2.5; ir7 = ir5 / r2
                dx, dy, dz = d[..., 0], d[..., 1], d[..., 2]
                md = mx * dx + my * dy
                B += k * (3 * md * dx * ir5 - mx * ir5 * r2).sum(0)
                # dBx/d(source position) = -dBx/dd
                gx = -k * (6 * mx * dx * ir5 + 3 * md * ir5 - 15 * md * dx * dx * ir7)
                gy = -k * (3 * my * dx * ir5 - 15 * md * dx * dy * ir7 + 3 * mx * dy * ir5)
                gz = -k * (-15 * md * dx * dz * ir7 + 3 * mx * dz * ir5)
                J[:, nz + i] += (gx * np.cos(th) + gy * np.sin(th)).sum(0)
                J[:, rg['k']] += (gz * sign).sum(0)
        return B, J

    def constraints(self):
        """Geometry constraints g >= 0 and their Jacobian."""
        geo, nz = self.geo, len(self.z)
        nv = nz + len(self.rings)
        g, J = [], []

        def add(val, coef):
            row = np.zeros(nv)
            for j, c in coef:
                row[j] += c
            g.append(val); J.append(row)

        if not self.center:                                             # clear of its own mirror slice
            add(self.z[0] - geo.pitch / 2, [(0, 1)])
        for j in range(nz - 1):                                         # slices do not overlap
            add(self.z[j + 1] - self.z[j] - geo.pitch, [(j + 1, 1), (j, -1)])
        add(geo.zmax - self.z[-1], [(nz - 1, -1)])
        for i, rg in enumerate(self.rings):
            r, c = rg['r'], geo.fit / rg['r']
            add(r - geo.rmin, [(nz + i, 1)])
            add(geo.rmax - r, [(nz + i, -1)])
            add(np.pi / rg['n'] - np.arcsin(min(c, 1)), [(nz + i, c / r / np.sqrt(max(1 - c * c, 1e-12)))])
            rlo = geo.rs / geo.ripple ** (1 / (rg['n'] - 1))              # ripple: n >= nmin(r) <=> r >= rlo
            add(r - rlo, [(nz + i, 1)])
            if i + 1 < len(self.rings) and self.rings[i + 1]['k'] == rg['k']:   # rings in a slice do not overlap
                add(self.rings[i + 1]['r'] - r - geo.dr, [(nz + i + 1, 1), (nz + i, -1)])
        return np.array(g), np.array(J)

    def ok(self, sensors, cap):
        B = self.field(sensors)
        return self.constraints()[0].min() >= -1e-6 and ppm(B) <= cap, B.mean()

    def table(self):
        return [dict(z=round(z, 2), rings=[dict(n=rg['n'], r=round(rg['r'], 2))
                                           for rg in self.rings if rg['k'] == k]) for k, z in enumerate(self.z)]

    def summary(self):
        return '\n'.join(f"  z = +-{t['z']:6.2f}: " + ', '.join(f"{q['n']} @ r {q['r']:.1f}" for q in t['rings'])
                         for t in self.table())

    def ledger(self):
        rows = []
        for rg in self.rings:
            pos, rot = ring_positions(rg['r'], self.z[rg['k']], rg['n'])
            rows += [dict(x=p[0], y=p[1], z=p[2], rot_z_deg=np.rad2deg(t) % 360) for p, t in zip(pos, rot)]
        df = pd.DataFrame(rows)
        df['size_xyz_mm'] = f'{self.geo.a} x {self.geo.a} x {self.geo.h}'
        return df


def refine(d, sensors, eps, iters=10, step=2.0):
    """SLSQP: maximize mean Bx s.t. |B_k - mean| <= eps*mean and geometry, each call within +-step mm of its start
    (a trust region, so the dipole-model correction stays valid). Returns the best feasible B0, or None."""
    cap = 2e6 * eps + 1e-6
    feas, b = d.ok(sensors, cap)
    best = (b, d.pack()) if feas else None
    for _ in range(iters):
        if best is not None:
            d.unpack(best[1])                                           # always continue from the best feasible
        x0 = d.pack()
        lo, hi = x0 - step, x0 + step
        if d.center:
            lo[0] = hi[0] = 0.0
        corr = d.field(sensors) - d.field_dipole(sensors)[0]
        cache = {}

        def fld(x):
            if x.tobytes() not in cache:
                d.unpack(x)
                B, J = d.field_dipole(sensors)
                cache.clear(); cache[x.tobytes()] = (B + corr, J)
            return cache[x.tobytes()]

        def obj(x):
            B, J = fld(x)
            return -B.mean(), -J.mean(0)

        def cons(x):
            B, _ = fld(x); m = B.mean(); d.unpack(x)
            return np.r_[1e6 * ((1 + eps) * m - B) / m, 1e6 * (B - (1 - eps) * m) / m, d.constraints()[0]]

        def cjac(x):
            B, J = fld(x); m, dm = B.mean(), J.mean(0); d.unpack(x)
            t = 1e6 * (B[:, None] * dm[None] - J * m) / m ** 2               # d/dx of 1e6*(m - B)/m
            return np.vstack([t, -t, d.constraints()[1]])

        res = minimize(obj, x0, jac=True, method='SLSQP', bounds=list(zip(lo, hi)),
                       constraints=[dict(type='ineq', fun=cons, jac=cjac)], options=dict(maxiter=100, ftol=1e-9))
        d.unpack(res.x)
        feas, b = d.ok(sensors, cap)
        log(f'    SLSQP: B0 {b:.2f} mT, {ppm(d.field(sensors)):.1f} ppm, geometry {d.constraints()[0].min():+.3f}, {res.message}')
        if feas and (best is None or b > best[0] + 0.01):
            best = (b, d.pack())
        elif best is not None:
            break
    if best is None:
        return None
    d.unpack(best[1])
    return best[0]


def add_magnets(d, sensors, eps, deadline):
    """Greedy: add one magnet to a ring (smallest ppm hit first), re-refine, keep if B0 rises."""
    b = refine(d, sensors, eps)
    if b is None:
        return None
    improved = True
    while improved and time.time() < deadline:
        improved = False
        cands = []
        for i in range(len(d.rings)):
            e = d.copy(); e.rings[i]['n'] += 1
            if e.constraints()[0].min() >= -1e-6:
                cands.append((ppm(e.field(sensors)), e))
        for _, e in sorted(cands, key=lambda c: c[0])[:6]:
            be = refine(e, sensors, eps, iters=3)
            if be is not None and be > b + 0.05:
                d.z, d.rings, b, improved = e.z, e.rings, be, True
                log(f'  +1 magnet: {d.count()} magnets, B0 {b:.2f} mT')
                break
    return b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--sphere-d', type=float, default=120)
    ap.add_argument('--ppm', type=float, default=60, help='peak-to-peak limit on the sphere')
    ap.add_argument('--margin', type=float, default=8, help='ppm held back during optimization')
    ap.add_argument('--bore', type=float, default=140, help='clear bore diameter, mm')
    ap.add_argument('--max-od', type=float, default=325, help='holder OD incl. plastic, mm')
    ap.add_argument('--max-length', type=float, default=360, help='end-to-end magnet length, mm')
    ap.add_argument('--a', type=float, default=12.7, help='block cross-section in the ring plane, mm')
    ap.add_argument('--h', type=float, default=25.4, help='block length along z, mm')
    ap.add_argument('--br', type=float, default=1300, help='remanence, mT')
    ap.add_argument('--bore-wall', type=float, default=3)
    ap.add_argument('--od-wall', type=float, default=3)
    ap.add_argument('--slice-wall', type=float, default=3)
    ap.add_argument('--layer-wall', type=float, default=2)
    ap.add_argument('--ring-wall', type=float, default=1.5)
    ap.add_argument('--milp-slack', type=float, default=1.0, help='MILP uses this fraction of the ppm budget')
    ap.add_argument('--dr-grid', type=float, default=10)
    ap.add_argument('--dz', type=float, default=4)
    ap.add_argument('--sensors', type=int, default=1200, help='points on the full sphere; one quarter is used')
    ap.add_argument('--ripple', type=float, default=2e-6, help='max (R/r)^(n-1) per ring')
    ap.add_argument('--milp-time', type=float, default=180)
    ap.add_argument('--time-limit', type=float, default=900)
    a = ap.parse_args()

    geo = Geometry(a)
    sensors = sphere_points(a.sensors, a.sphere_d / 2)
    sensors = sensors[(sensors[:, 1] >= 0) & (sensors[:, 2] >= 0)]     # Bx is mirror symmetric in y and z
    eps = (a.ppm - a.margin) * 1e-6 / 2
    log(f'r {geo.rmin:.1f}-{geo.rmax:.1f} mm, axial pitch {geo.pitch} mm, z <= {geo.zmax:.1f} mm')

    rings, b_lin, msg = milp_design(geo, sensors, a.milp_slack * eps, a.dr_grid, a.dz, a.milp_time)
    d = Design(geo, rings)
    B = d.field(sensors)
    log(f'MILP ({msg}): {b_lin:.1f} mT linear, exact {B.mean():.1f} mT / {ppm(B):.0f} ppm, {d.count()} magnets\n'
        + d.summary())
    start = d.copy()
    b = add_magnets(d, sensors, eps, T0 + a.time_limit)
    if b is None:                                   # refinement never got under the cap: keep the MILP design
        log('refinement found no feasible design, keeping the MILP design')
        d = start

    dense = sphere_points(20000, a.sphere_d / 2)
    Bd = d.field(dense)
    info = dict(B0_mT=round(float(Bd.mean()), 2), B_center_mT=round(float(d.field(np.zeros((1, 3)))[0]), 2),
                ppm_120mm_sphere=round(float(ppm(Bd)), 1), magnets=d.count(), bore=a.bore, br_mT=a.br,
                milp_linear_mT=round(float(b_lin), 1), slices=d.table())
    stem = a.out.rsplit('.', 1)[0]
    d.ledger().to_excel(a.out, index=False)
    with open(stem + '.json', 'w') as f:
        json.dump(info, f, indent=1)
    log(f'done: B0 {Bd.mean():.2f} mT, {ppm(Bd):.1f} ppm (20k points), {d.count()} magnets\n' + d.summary())


if __name__ == '__main__':
    main()
