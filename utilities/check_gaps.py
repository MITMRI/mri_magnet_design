"""Exact minimum gap between every pair of blocks in a design ledger (rotated 12.7 x 12.7 squares in the ring
plane, 25.4 mm along z). Prints the smallest gap within a ring, between rings of one slice, and between slices.

    python utilities/check_gaps.py results/12cm/final/<design>.xlsx
"""
import sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

A, H = 12.7, 25.4


def square(x, y, deg):
    t = np.deg2rad(deg)
    c = np.array([[1, 1], [-1, 1], [-1, -1], [1, -1]]) * A / 2
    rot = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    return c @ rot.T + [x, y]


def seg_dist(p, a, b):
    ab = b - a
    s = np.clip(np.dot(p - a, ab) / np.dot(ab, ab), 0, 1)
    return np.linalg.norm(p - a - s * ab)


def poly_dist(P, Q):
    """Distance between two convex polygons; 0 if they overlap (separating-axis test)."""
    for poly in (P, Q):
        for i in range(4):
            e = poly[(i + 1) % 4] - poly[i]; n = np.array([-e[1], e[0]])
            if (P @ n).max() < (Q @ n).min() or (Q @ n).max() < (P @ n).min():
                break
        else:
            continue
        break
    else:
        return 0.0
    return min(min(seg_dist(p, Q[i], Q[(i + 1) % 4]) for p in P for i in range(4)),
               min(seg_dist(q, P[i], P[(i + 1) % 4]) for q in Q for i in range(4)))


df = pd.read_excel(sys.argv[1])
xyz = df[['x', 'y', 'z']].values; rot = df['rot_z_deg'].values
r = np.hypot(xyz[:, 0], xyz[:, 1])
sq = [square(x, y, d) for (x, y, _), d in zip(xyz, rot)]
best = {}
for i, j in cKDTree(xyz).query_pairs(40.0):
    dz = max(0.0, abs(xyz[i, 2] - xyz[j, 2]) - H)
    gap = np.hypot(dz, poly_dist(sq[i], sq[j]))
    kind = ('same ring' if abs(r[i] - r[j]) < 1 else 'rings in one slice') if abs(xyz[i, 2] - xyz[j, 2]) < 1 else 'between slices'
    best[kind] = min(best.get(kind, np.inf), gap)
print(f'{len(df)} blocks; outer edge r {max(np.abs(s).max() for s in []) if False else (np.hypot(*np.vstack(sq).T)).max():.2f} mm, '
      f'inner edge r {min(np.hypot(*s.T).min() for s in sq):.2f} mm, z span {2 * (np.abs(xyz[:, 2]).max() + H / 2):.1f} mm')
for k, v in best.items():
    print(f'  min gap {k}: {v:.2f} mm')
print('OVERLAP' if min(best.values()) <= 0 else 'no overlaps')
