"""
Axial forces while blocks are pushed into one plate's pockets, for pockets with a floor on the outward face
(away from the midplane) and an open face toward the midplane, held only by friction.

Seated blocks in one plate share a z centre and put no axial force on each other. While block b is being slid in
through the open face it sits d mm out of its seat, and then
  * every seated block feels an axial force from b; the part pointing toward the open face is a pull-out
    load the friction fit must resist,
  * b itself feels the field of all its seated neighbours (worst case: b is the last block of the plate);
    "resist" > 0 = b is pushed back out (you must press it in), < 0 = it is sucked into its pocket.
Same pole-face surface-charge integral as forces_12cm.py; other plates are assumed absent (plates are filled
before stacking).

    python utilities/insertion_12cm.py designs/12cm_bore200/placements.csv
"""
import argparse
from concurrent.futures import ProcessPoolExecutor

import magpylib as magpy
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation as R

MU0 = 4e-7 * np.pi
DIM = np.array([12.7, 12.7, 25.4])
BR = 1300.0
NQ = 5
D = np.array([0.5, 1, 2, 3, 5, 8, 12.7, 18, 25.4, 35])   # mm out of the seat, toward the open face


def face_points():
    gy, wy = np.polynomial.legendre.leggauss(NQ)
    gz, wz = np.polynomial.legendre.leggauss(2 * NQ)
    Y, Z = np.meshgrid(gy * DIM[1] / 2, gz * DIM[2] / 2, indexing='ij')
    W = np.outer(wy * DIM[1] / 2, wz * DIM[2] / 2).ravel()
    loc = np.vstack([np.c_[np.full(Y.size, s * (DIM[0] / 2 + 1e-3)), Y.ravel(), Z.ravel()] for s in (1, -1)])
    return loc, np.r_[W, -W]


LOC, SW = face_points()
K = BR * 1e-3 / MU0 * 1e-9                                 # (A/m) * mm^2 -> m^2 * mT -> T


def pts_of(p, q):
    return R.from_quat(np.repeat(q, len(LOC), 0)).apply(np.tile(LOC, (len(p), 1))) + np.repeat(p, len(LOC), 0)


def fz_on(obs_p, obs_q, src_p, src_q):
    """Axial force (N) on each observer block from all source blocks (summed)."""
    pts = pts_of(obs_p, obs_q)
    m = len(pts)
    Bz = np.zeros(m)
    for s in range(len(src_p)):
        Bz += magpy.getB('Cuboid', observers=pts, position=src_p[s], orientation=R.from_quat(src_q[s]),
                         magnetization=(BR, 0, 0), dimension=DIM)[:, 2]
    return K * (np.tile(SW, len(obs_p)) * Bz).reshape(len(obs_p), -1).sum(1)


def plate_job(args):
    z, p, q = args
    u = -np.sign(z)                                        # open face points toward the midplane
    n = len(p)
    pull = np.zeros((n, len(D)))                           # worst pull-out on any seated neighbour while b goes in
    who = np.zeros((n, len(D)), int)
    resist = np.zeros((n, len(D)))                         # force on b itself, + = pushed back out
    for b in range(n):
        others = np.r_[0:b, b + 1:n]
        for k, d in enumerate(D):
            pb = p[b] + [0, 0, u * d]
            f = u * fz_on(p[others], q[others], pb[None], q[b][None])
            j = int(np.argmax(f))
            pull[b, k], who[b, k] = f[j], others[j]
            resist[b, k] = u * fz_on(pb[None], q[b][None], p[others], q[others])[0]
    return z, pull, resist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('csv')
    ap.add_argument('--jobs', type=int, default=6)
    a = ap.parse_args()
    df = pd.read_csv(a.csv)
    df['zk'] = df.z.round(1)
    jobs = []
    for z, g in df.groupby('zk'):
        jobs.append((z, g[['x', 'y', 'z']].to_numpy(),
                     R.from_euler('z', g.rot_z_deg.to_numpy()[:, None], degrees=True).as_quat()))
    with ProcessPoolExecutor(a.jobs) as ex:
        res = list(ex.map(plate_job, jobs))

    print(f'insertion d (mm out of seat): {", ".join(f"{d:g}" for d in D)}')
    print('  plate z  blocks | worst pull-out on a seated neighbour N (at d) | inserted block: max push-back N (at d)'
          ' / max suck-in N (at d)')
    worst = (0, None)
    for z, pull, resist in sorted(res):
        i = np.unravel_index(np.argmax(pull), pull.shape)
        r = np.unravel_index(np.argmax(resist), resist.shape)
        s = np.unravel_index(np.argmin(resist), resist.shape)
        print(f'  {z:+7.1f} {len(pull):6d} | {pull[i]:10.1f} ({D[i[1]]:g})                          | '
              f'{resist[r]:8.1f} ({D[r[1]]:g}) / {-resist[s]:6.1f} ({D[s[1]]:g})')
        if pull[i] > worst[0]:
            worst = (pull[i], z)
    print(f'worst pull-out on a seated block while a neighbour is inserted: {worst[0]:.1f} N '
          f'({worst[0] / 9.81:.2f} kgf), plate z={worst[1]:+.1f}')
    print('per-d maximum pull-out over all plates (N): ' +
          ', '.join(f'{d:g}: {max(r[1][:, k].max() for r in res):.1f}' for k, d in enumerate(D)))


if __name__ == '__main__':
    main()
