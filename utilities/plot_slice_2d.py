"""
2D cross-section of one slice (plate) of a design: the 12.7 x 12.7 mm block footprints with magnetization arrows,
the holder annulus (bore ID to OD), the 120 mm DSV and each ring's radius and block count.

    python utilities/plot_slice_2d.py designs/12cm_bore200/placements.csv [--z 166.2] [--id 201.1 --od 325]

Without --z it draws the slice with the most blocks. Writes slice_z<z>.png next to the csv.
"""
import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, Polygon, Wedge

matplotlib.use('Agg')
A = 12.7

ap = argparse.ArgumentParser()
ap.add_argument('csv')
ap.add_argument('--z', type=float)
ap.add_argument('--id', type=float, default=201.1, help='holder inner diameter, mm')
ap.add_argument('--od', type=float, default=325.0, help='holder outer diameter, mm')
ap.add_argument('--bore', type=float, default=200.0, help='clear bore, mm')
ap.add_argument('--dsv', type=float, default=120.0)
a = ap.parse_args()

df = pd.read_csv(a.csv)
df['zk'] = df.z.round(1)
z = df.zk.value_counts().idxmax() if a.z is None else a.z
s = df[df.zk == round(z, 1)].copy()
s['r'] = np.hypot(s.x, s.y).round(1)

INK, MUTED, HOLDER, BLOCK, NORTH = '#1f2328', '#6e7781', '#e8e4dc', '#9aa4ae', '#c4314b'
fig, ax = plt.subplots(figsize=(9, 9), dpi=150)
ax.add_patch(Wedge((0, 0), a.od / 2, 0, 360, width=(a.od - a.id) / 2, facecolor=HOLDER, edgecolor=MUTED, lw=1))
ax.add_patch(Circle((0, 0), a.bore / 2, fill=False, ls=(0, (6, 4)), ec=INK, lw=1))
ax.add_patch(Circle((0, 0), a.dsv / 2, fill=False, ls=(0, (2, 3)), ec=MUTED, lw=1))

half = A / 2
for row in s.itertuples():
    t = np.deg2rad(row.rot_z_deg)
    R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    c = np.array([row.x, row.y])
    sq = (R @ np.array([[-half, -half], [half, -half], [half, half], [-half, half]]).T).T + c
    north = (R @ np.array([[0, -half], [half, -half], [half, half], [0, half]]).T).T + c    # +x half = north side
    ax.add_patch(Polygon(sq, closed=True, facecolor=BLOCK, edgecolor='white', lw=0.6))
    ax.add_patch(Polygon(north, closed=True, facecolor=NORTH, edgecolor='none', alpha=0.85))
    ax.add_patch(Polygon(sq, closed=True, fill=False, edgecolor=INK, lw=0.4))
    d = R @ np.array([4.2, 0])
    ax.annotate('', xy=c + d, xytext=c - d, arrowprops=dict(arrowstyle='-|>', lw=0.7, color='white',
                                                             mutation_scale=6))

# ring labels to the right, outermost on top, leaders to each ring at 35 deg (no crossings)
rings = list(s.groupby('r'))[::-1]
ang = np.deg2rad(35)
for i, (r, g) in enumerate(rings):
    ax.annotate(f'r = {r:g} mm, {len(g)} blocks', xy=(r * np.cos(ang), r * np.sin(ang)),
                xytext=(a.od / 2 + 12, 150 - 24 * i), fontsize=9, color=INK, va='center',
                arrowprops=dict(arrowstyle='-', lw=0.6, color=MUTED, shrinkB=0))
ax.text(0, a.bore / 2 - 18, f'{a.bore:g} mm clear bore', ha='center', fontsize=9, color=INK)
ax.text(0, a.dsv / 2 - 12, f'{a.dsv:g} mm DSV', ha='center', fontsize=8, color=MUTED)
ax.text(0, -a.od / 2 - 12, f'holder ID {a.id:g} mm  /  OD {a.od:g} mm', ha='center', fontsize=9, color=INK)
ax.annotate('', xy=(18, 0), xytext=(-18, 0), arrowprops=dict(arrowstyle='-|>', lw=1.4, color=INK))
ax.text(0, 6, 'B0 (+x)', fontsize=10, color=INK, ha='center')

ax.set_aspect('equal')
lim = a.od / 2 + 30
ax.set_xlim(-lim, lim + 95)
ax.set_ylim(-lim, lim)
ax.axis('off')
ax.set_title(f'Slice z = {z:+g} mm: {len(s)} blocks (12.7 x 12.7 mm footprint, red = north half)',
             fontsize=11, color=INK)
out = a.csv.rsplit('/', 1)[0] + f'/slice_z{z:+g}.png'
fig.savefig(out, bbox_inches='tight', facecolor='white')
print(out)
