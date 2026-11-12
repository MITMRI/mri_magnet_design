"""
Monte Carlo as-built homogeneity of a 12 cm design (placements.csv from field_design_12cm.py).

Each block gets independent random errors, all uniform within +-tol:
  - position: x, y, z each uniform in +-pos_tol (mm)
  - magnetization direction: tilted off its nominal axis by an angle uniform-on-the-cap within tol_deg, random azimuth
  - magnetization strength: Br * (1 + u), u uniform in +-tol_br
ppm = peak-to-peak of Bx over points on the sphere surface / mean Bx (Bx is harmonic, extremes on the surface).

Field model: exact cuboid fields (magpylib), linearized in the errors. Per block we store Bx on the sensors for unit
magnetization along each local axis (3 basis fields) and dBx/d(position) (central differences), so a trial is a
matrix product. The linearization is checked against exact magpylib on a few trials (--check).

    python utilities/tolerance_12cm.py designs/12cm_bore200_84mT/placements.csv
"""

import argparse
import time

import magpylib as magpy
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

BR, A, H = 1300.0, 12.7, 25.4


def sphere_points(n, r):
    i = np.arange(n) + 0.5
    phi, th = np.arccos(1 - 2 * i / n), np.pi * (1 + 5 ** 0.5) * i
    return r * np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], axis=1)


def bx(pos, rot, mag, sensors):
    """Bx (mT) on sensors of each block: (N, S). mag: (N, 3) local magnetization (mT)."""
    nm, ns = len(pos), len(sensors)
    B = magpy.getB('Cuboid', observers=np.tile(sensors, (nm, 1)), position=np.repeat(pos, ns, 0),
                   orientation=Rotation.from_euler('z', np.repeat(rot, ns)[:, None]), magnetization=np.repeat(mag, ns, 0),
                   dimension=(A, A, H))
    return B[:, 0].reshape(nm, ns)


def chunked(f, n, *args, step=100):
    return np.vstack([f(*[a[i:i + step] for a in args[:3]], *args[3:]) for i in range(0, n, step)])


def tilt_dirs(rng, n, tol_deg):
    """Unit vectors tilted from +x by an angle uniform on the spherical cap of half-angle tol_deg."""
    c = 1 - rng.random(n) * (1 - np.cos(np.deg2rad(tol_deg)))
    s, az = np.sqrt(1 - c * c), 2 * np.pi * rng.random(n)
    return np.stack([c, s * np.cos(az), s * np.sin(az)], 1)


def ppm(B):
    return 1e6 * (B.max(-1) - B.min(-1)) / B.mean(-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('placements')
    ap.add_argument('--sensors', type=int, default=4000)
    ap.add_argument('--trials', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--pos-tol', type=float, default=0.05, help='placement tolerance, +-mm per axis')
    ap.add_argument('--check', type=int, default=5, help='exact magpylib trials to validate the linear model')
    a = ap.parse_args()

    df = pd.read_csv(a.placements)
    pos, rot = df[['x', 'y', 'z']].to_numpy(), np.deg2rad(df['rot_z_deg'].to_numpy())
    n = len(pos)
    S = sphere_points(a.sensors, 60.0)
    t0 = time.time()
    eye = np.eye(3)
    basis = [chunked(bx, n, pos, rot, np.tile(BR * eye[j], (n, 1)), S) for j in range(3)]  # local Mx, My, Mz
    h = 0.01
    grad = []
    for j in range(3):
        dp = h * eye[j]
        grad.append((chunked(bx, n, pos + dp, rot, np.tile([BR, 0, 0], (n, 1)), S)
                     - chunked(bx, n, pos - dp, rot, np.tile([BR, 0, 0], (n, 1)), S)) / (2 * h))
    basis, grad = np.array(basis), np.array(grad)          # (3, N, S) each
    B_ideal = basis[0].sum(0)
    print(f'{n} blocks, {a.sensors} sensor points, tables in {time.time() - t0:.0f} s')
    print(f'ideal: B0 {B_ideal.mean():.2f} mT, {ppm(B_ideal):.1f} ppm')

    def trials(rng, k, pos_tol, ang_tol, br_tol):
        m = np.ones((k, n, 3)) * [1, 0, 0]
        if ang_tol > 0:
            m = tilt_dirs(rng, k * n, ang_tol).reshape(k, n, 3)
        m = m * (1 + br_tol * (2 * rng.random((k, n, 1)) - 1))
        dp = pos_tol * (2 * rng.random((k, n, 3)) - 1)
        B = np.einsum('knj,jns->ks', m, basis, optimize=True) + np.einsum('knj,jns->ks', dp, grad, optimize=True)
        return B, m, dp

    p = a.pos_tol
    cases = [(f'position +-{p} mm only', p, 0, 0),
             ('+ angle +-1 deg', p, 1, 0),
             ('+ angle +-3 deg', p, 3, 0),
             ('+ Br +-0.5%', p, 0, 0.005),
             ('+ Br +-1.5%', p, 0, 0.015),
             (f'all: pos {p}, 1 deg, 0.5%', p, 1, 0.005),
             (f'all: pos {p}, 3 deg, 1.5%', p, 3, 0.015)]

    # validate the linear model against exact magpylib on the loosest case
    rng = np.random.default_rng(a.seed + 99)
    B_lin, m, dp = trials(rng, a.check, a.pos_tol, 3, 0.015)
    for i in range(a.check):
        B_ex = chunked(bx, n, pos + dp[i], rot, BR * m[i], S).sum(0)
        print(f'  check {i}: linear {ppm(B_lin[i]):7.1f} ppm, exact {ppm(B_ex):7.1f} ppm')

    print(f'\n{a.trials} random builds per case; ppm over the 120 mm sphere surface (ideal {ppm(B_ideal):.1f})')
    print(f'{"case":32s} {"median":>8s} {"mean":>8s} {"p90":>8s} {"p99":>8s} {"max":>8s}   B0 sd (mT)')
    for name, pt, at, bt in cases:
        rng = np.random.default_rng(a.seed)
        P, B0 = [], []
        for _ in range(0, a.trials, 100):
            B = trials(rng, 100, pt, at, bt)[0]
            P.append(ppm(B)); B0.append(B.mean(-1))
        P, B0 = np.concatenate(P), np.concatenate(B0)
        q = np.percentile(P, [50, 90, 99])
        print(f'{name:32s} {q[0]:8.1f} {P.mean():8.1f} {q[1]:8.1f} {q[2]:8.1f} {P.max():8.1f}   {B0.std():.3f}')


if __name__ == '__main__':
    main()
