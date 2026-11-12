"""
Magnetic forces and torques on every block of a design, i.e. the loads the printed ring holders must carry.

Each N42SH block is a uniformly magnetized cuboid (magnetization along its local x, rotated by rot_z about z).
The force on it from all the other blocks is the surface-charge integral over its two pole faces,
    F = sum_faces  sigma * B_ext dA,   sigma = +-Br/mu0,
and the torque about its centre is sum r x dF. Each face uses Gauss-Legendre quadrature; B_ext is the total
field of the design minus the block's own field at the same points (its self force is zero).

    python utilities/forces_12cm.py designs/12cm_bore200/placements.csv [--br 1300]
"""
import argparse
from concurrent.futures import ProcessPoolExecutor

import magpylib as magpy
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation as R

MU0 = 4e-7 * np.pi
DIM = np.array([12.7, 12.7, 25.4])                  # mm, local x = magnetization
_G = {}


def _init(flat, pos, quat, br):
    _G.update(flat=flat, pos=pos, quat=quat, br=br)


def _partial(chunk):
    """Field (mT) of source blocks s..s+k at every observer point, summed over those sources."""
    s, k = chunk
    flat, m = _G['flat'], len(_G['flat'])
    B = magpy.getB('Cuboid', observers=np.tile(flat, (k, 1)), position=np.repeat(_G['pos'][s:s + k], m, 0),
                   orientation=R.from_quat(np.repeat(_G['quat'][s:s + k], m, 0)),
                   magnetization=(_G['br'], 0, 0), dimension=DIM)
    return B.reshape(k, m, 3).sum(0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('csv')
    ap.add_argument('--br', type=float, default=1300, help='remanence, mT')
    ap.add_argument('--nq', type=int, default=6, help='quadrature points per face edge (y), 2x along z')
    ap.add_argument('--chunk', type=int, default=2, help='source blocks per getB call (~200 MB per block per worker)')
    ap.add_argument('--jobs', type=int, default=7, help='worker processes (total RAM ~ jobs x chunk x 200 MB)')
    a = ap.parse_args()

    df = pd.read_csv(a.csv)
    pos = df[['x', 'y', 'z']].to_numpy()
    rot = R.from_euler('z', df.rot_z_deg.to_numpy()[:, None], degrees=True)
    n = len(df)

    # quadrature on a +x face in local coords (mm), weights in mm^2
    gy, wy = np.polynomial.legendre.leggauss(a.nq)
    gz, wz = np.polynomial.legendre.leggauss(2 * a.nq)
    Y, Z = np.meshgrid(gy * DIM[1] / 2, gz * DIM[2] / 2, indexing='ij')
    W = np.outer(wy * DIM[1] / 2, wz * DIM[2] / 2).ravel()
    eps = 1e-3                                       # mm, just outside the face
    face = []
    for s in (+1, -1):
        face.append(np.c_[np.full(Y.size, s * (DIM[0] / 2 + eps)), Y.ravel(), Z.ravel()])
    loc = np.vstack(face)                            # (q, 3) local points
    sig = np.r_[np.ones(Y.size), -np.ones(Y.size)]   # charge sign per point
    w = np.r_[W, W]
    q = len(loc)

    # observer points for every block, in global coords
    pts = np.stack([rot[i].apply(loc) + pos[i] for i in range(n)])          # (n, q, 3)
    flat = pts.reshape(-1, 3)

    # total field, a few sources at a time so memory stays ~100s of MB (a Collection getB tiles all n x n*q at once)
    quat = rot.as_quat()
    chunks = [(s, min(a.chunk, n - s)) for s in range(0, n, a.chunk)]
    with ProcessPoolExecutor(a.jobs, initializer=_init, initargs=(flat, pos, quat, a.br)) as ex:
        Btot = sum(ex.map(_partial, chunks, chunksize=4))                  # mT
    Bself = magpy.getB('Cuboid', observers=flat, position=np.repeat(pos, q, 0),
                       orientation=R.from_quat(np.repeat(rot.as_quat(), q, 0)),
                       magnetization=(a.br, 0, 0), dimension=DIM)
    Bext = (Btot - Bself).reshape(n, q, 3) * 1e-3                          # T

    sigma = a.br * 1e-3 / MU0                                               # A/m
    dF = sigma * (sig * w * 1e-6)[None, :, None] * Bext                     # N per point
    F = dF.sum(1)                                                           # (n, 3) global
    r = (pts - pos[:, None, :]) * 1e-3                                      # m
    T = np.cross(r, dF).sum(1)                                              # N m

    phi = np.arctan2(pos[:, 1], pos[:, 0])
    er = np.c_[np.cos(phi), np.sin(phi), np.zeros(n)]
    et = np.c_[-np.sin(phi), np.cos(phi), np.zeros(n)]
    out = df.copy()
    out['r'] = np.hypot(pos[:, 0], pos[:, 1]).round(2)
    out['Fr'], out['Ft'], out['Fz'] = (F * er).sum(1), (F * et).sum(1), F[:, 2]
    out['F'] = np.linalg.norm(F, axis=1)
    out['Tz'] = T[:, 2]
    out['T'] = np.linalg.norm(T, axis=1)
    out.to_csv(a.csv.rsplit('/', 1)[0] + '/forces.csv', index=False)

    kg = 1 / 9.81
    print(f'{n} blocks, Br {a.br:.0f} mT, {q} quadrature points per block')
    print(f'net force on whole magnet (should be ~0): {F.sum(0).round(2)} N')
    print('\nper block (N, N m):')
    for c in ['F', 'Fr', 'Ft', 'Fz', 'T', 'Tz']:
        i = out[c].abs().idxmax()
        print(f'  max |{c:2}| = {abs(out[c][i]):7.2f}  ({abs(out[c][i]) * (kg if c[0] == "F" else 1):6.2f}'
              f'{" kgf" if c[0] == "F" else " N m"})  at r={out.r[i]:.1f} z={out.z[i]:+.1f} rot={out.rot_z_deg[i]:.1f}')

    print('\nper ring (holder ring = one radius at one z):')
    print('     z      r    n   sum Fz N   max|F| N  max|Fr| N  max|Ft| N  max|Fz| N  max|T| Nm  mean Fr N')
    for (z, rr), g in out.groupby([out.z.round(1), out.r]):
        print(f'{z:+7.1f} {rr:6.1f} {len(g):4d} {g.Fz.sum():+10.1f} {g.F.max():10.2f} {g.Fr.abs().max():10.2f}'
              f' {g.Ft.abs().max():10.2f} {g.Fz.abs().max():10.2f} {g["T"].max():10.3f} {g.Fr.mean():+10.2f}')

    print('\nper slice (all rings at one z) and cumulative axial load toward the midplane:')
    zs = np.sort(out.z.round(1).unique())
    sl = out.groupby(out.z.round(1)).Fz.sum()
    for z in zs:
        cum = sl[zs[zs >= z]].sum() if z > 0 else sl[zs[zs <= z]].sum()
        print(f'  z {z:+7.1f}: slice Fz {sl[z]:+9.1f} N ({sl[z] * kg:+7.1f} kgf), '
              f'outboard stack incl. this slice {cum:+9.1f} N ({cum * kg:+7.1f} kgf)')
    half = out[out.z > 0].Fz.sum()
    print(f'\nforce between the two halves (sum Fz for z>0): {half:+.1f} N ({half * kg:+.1f} kgf);'
          f' negative = halves attract')


if __name__ == '__main__':
    main()
