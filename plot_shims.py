# -*- coding: utf-8 -*-
"""
Picture of a shim_design.py solution in the slide-in sleeve: 3D view of the main array (front half cut away) with the
trays and shim cubes and one tray pulled out, an end view of the sleeve, and the busiest tray seen from the bore axis.
Sleeve dimensions are the make_shim_cad.py defaults.

usage:
    python plot_shims.py ledger.xlsx shim_solution.npz out.png [holder OD, mm]
"""
import sys, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Circle, Polygon, Patch, Rectangle

L = pd.read_excel(sys.argv[1]); sol = np.load(sys.argv[2]); out = sys.argv[3]
od = float(sys.argv[4]) if len(sys.argv) > 4 else 185.9
bore, a_m = 100.0, float(L['Magnet_length'].iloc[0])
trs, z, kind, psi = sol['trs'], sol['z'], sol['kind'], sol['psi']
rtr, a, ntr = float(sol['rtr']), float(sol['a']), int(sol['ntr'])
th = 2*np.pi*(np.arange(ntr) + 0.5)/ntr
tw = a*np.sqrt(2) + 0.2 + 2*0.8                    # tray width: a flat cube can turn to any angle (make_shim_cad.py)
r_clear = rtr - a/2 - 0.8 - 0.15 - 1.2             # sleeve inner radius: tray floor, slide clearance, sleeve wall
tray_len = 2*70 + a*np.sqrt(2) + 4
COLS = {'flat': '#009E73', 'out': '#D55E00', 'in': '#0072B2'}   # Okabe-Ito
def cat(k): return 'flat' if k == 0 else ('out' if k > 0 else 'in')

def frame(t):
    return np.array([np.cos(t), np.sin(t), 0]), np.array([-np.sin(t), np.cos(t), 0]), np.array([0, 0, 1.0])

def box(c, ex, ey, ez, sx, sy, sz):
    v = np.array([c + i*sx*ex + j*sy*ey + k*sz*ez for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)])
    return [v[q] for q in [[0, 1, 3, 2], [4, 5, 7, 6], [0, 1, 5, 4], [2, 3, 7, 6], [0, 2, 6, 4], [1, 3, 7, 5]]]

def cube_axes(t, k, p):
    er, et, ez = frame(t)
    if k != 0: return er, et, ez
    return er, np.cos(p)*et - np.sin(p)*ez, np.sin(p)*et + np.cos(p)*ez   # turned about the radial axis by psi

counts = np.bincount(trs, minlength=ntr); pull = int(np.argmax(counts)); PULL = 120.0; AZ = -55.0
vdir = np.array([np.cos(np.deg2rad(AZ)), np.sin(np.deg2rad(AZ))])

fig = plt.figure(figsize=(19, 15))
gs = fig.add_gridspec(3, 3, height_ratios=(1, 1, 0.5))
ax = fig.add_subplot(gs[0:2, 0:2], projection='3d')
polys = []
for _, row in L.iterrows():
    if row[['X-pos', 'Y-pos']].values.astype(float) @ vdir > 0: continue     # front half cut away
    polys += box(row[['X-pos', 'Y-pos', 'Z-pos']].values.astype(float), *frame(np.deg2rad(row['Z-rot'])), a_m/2, a_m/2, a_m/2)
ax.add_collection3d(Poly3DCollection(polys, facecolors='#c4c4c4', edgecolors='#8a8a8a', linewidths=0.1))
for t in range(ntr):
    er, et, ez = frame(th[t]); dz = PULL if t == pull else 0.0
    ax.add_collection3d(Poly3DCollection(box(rtr*er + dz*ez, er, et, ez, a/2 + 0.6, tw/2, tray_len/2),
                        facecolors=(0.55, 0.5, 0.75, 0.06 if dz == 0 else 0.25), edgecolors=(0.35, 0.3, 0.55, 0.25), linewidths=0.3))
sp, sc = [], []
for t, zz, k, p in zip(trs, z, kind, psi):
    er, _, ez = frame(th[t]); dz = PULL if t == pull else 0.0
    sp += box(rtr*er + (zz + dz)*ez, *cube_axes(th[t], k, p), a/2, a/2, a/2); sc += [COLS[cat(k)]]*6
ax.add_collection3d(Poly3DCollection(sp, facecolors=sc, edgecolors='k', linewidths=0.3))
phi = np.linspace(0, 2*np.pi, 200)
for R_ in (od/2, bore/2, r_clear):
    for zz in (-75, 75):
        ax.plot(R_*np.cos(phi), R_*np.sin(phi), zz, color='#555', lw=0.8, alpha=0.6)
er = frame(th[pull])[0]
ax.text(*(rtr*er*1.6 + [0, 0, PULL + tray_len/2 + 8]), 'tray slides out\nthrough the end cap', fontsize=10, ha='center')
lim = od/2 + 5
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, PULL + lim/2)
ax.set_box_aspect((2*lim, 2*lim, PULL + 1.5*lim)); ax.view_init(elev=22, azim=AZ)
ax.set_xlabel('x (mm), B0'); ax.set_ylabel('y (mm)'); ax.set_zlabel('z (mm), bore axis')
ax.text2D(0.02, 0.02, 'front half of the main array cut away', transform=ax.transAxes, fontsize=10, color='#555')
ax.set_title(f'Main array: {sys.argv[1].split("/")[-1]} ({len(L)} magnets)\n{ntr} shim trays in a sleeve inside the bore (r = {rtr:.1f} mm); {len(z)} shim cubes ({a:g} mm), '
             f'placed by LP density + gradient descent\nas-built {float(sol["before"]):.0f} ppm -> {float(sol["after"]):.0f} ppm shimmed (simulated build)', fontsize=12)

ax2 = fig.add_subplot(gs[0, 2])
ax2.add_patch(Circle((0, 0), od/2, fc='#f1eef8', ec='#555')); ax2.add_patch(Circle((0, 0), bore/2, fc='white', ec='#555'))
for _, row in L.iterrows():
    er, et, _ = frame(np.deg2rad(row['Z-rot'])); c = row[['X-pos', 'Y-pos']].values.astype(float)
    ax2.add_patch(Polygon([c + s1*a_m/2*er[:2] + s2*a_m/2*et[:2] for s1, s2 in ((-1, -1), (1, -1), (1, 1), (-1, 1))], fc='#c8c8c8', ec='#999', lw=0.3))
for t in range(ntr):
    er, et, _ = frame(th[t]); c = rtr*er[:2]; hw = a/2 + 0.6
    ax2.add_patch(Polygon([c + s1*hw*er[:2] + s2*tw/2*et[:2] for s1, s2 in ((-1, -1), (1, -1), (1, 1), (-1, 1))],
                          fc='#cfc6e8' if counts[t] else 'white', ec='#4b3f7a', lw=0.6))
ax2.add_patch(Circle((0, 0), r_clear, fc='none', ec='#4b3f7a', lw=0.8, ls='--'))
ax2.annotate('+x (B0)', (od/2, 0), (od/2 - 30, -12), fontsize=9, arrowprops=dict(arrowstyle='->'))
ax2.set_xlim(-od/2 - 5, od/2 + 5); ax2.set_ylim(-od/2 - 5, od/2 + 5); ax2.set_aspect('equal')
ax2.set_title(f'End view: sleeve slots (filled = tray has cubes, {np.sum(counts > 0)}/{ntr});\n'
              f'clear bore {2*(r_clear):.1f} mm (dashed)', fontsize=11)
ax2.set_xlabel('x (mm)'); ax2.set_ylabel('y (mm)')

# the busiest tray, seen from the bore axis looking out (z along, tangential up)
ax3 = fig.add_subplot(gs[2, :])
z0 = -tray_len/2
ax3.add_patch(Rectangle((z0, -tw/2), tray_len, tw, fc='#e4def3', ec='#4b3f7a', lw=1.2))
ax3.add_patch(Rectangle((z0 + tray_len, -tw/4), 10, tw/2, fc='#cfc6e8', ec='#4b3f7a', lw=1.2))
ax3.text(z0 + tray_len + 12, 0, 'handle', va='center', fontsize=10)
sel = np.flatnonzero(trs == pull)
for i in sel:
    k, p, zz = kind[i], psi[i], z[i]
    ang = -p if k == 0 else 0.0            # square turned by psi in the (z, tangential) plane
    ca, sa = np.cos(ang), np.sin(ang)
    corners = [(zz + ca*u - sa*v, -(sa*u + ca*v)) for u, v in ((-a/2, -a/2), (a/2, -a/2), (a/2, a/2), (-a/2, a/2))]
    ax3.add_patch(Polygon(corners, fc=COLS[cat(k)], ec='k', lw=0.8))
    if k == 0:
        dz_, dt_ = np.cos(p), np.sin(p)
        ax3.annotate('', (zz + dz_*a*0.45, dt_*a*0.45), (zz - dz_*a*0.45, -dt_*a*0.45),
                     arrowprops=dict(arrowstyle='-|>', lw=1.8, color='k', mutation_scale=14))
    else:
        ax3.plot(zz, 0, marker='o', ms=9, mfc='white', mec='k', mew=1.3)
        ax3.plot(zz, 0, marker='.' if k < 0 else 'x', ms=6, color='k', mew=1.6)
ax3.set_xlim(z0 - 3, z0 + tray_len + 30); ax3.set_ylim(-tw, tw); ax3.set_aspect('equal')
ax3.set_xlabel('z along the bore (mm)'); ax3.set_ylabel('tangential (mm)')
ax3.set_title(f'The pulled-out tray, seen from the bore axis: {len(sel)} cubes at free z positions; arrows = magnetization '
              f'(flat cubes, any angle); dot = radially inward (toward you), cross = outward', fontsize=11)
fig.legend(handles=[Patch(fc=COLS['flat'], ec='k', label='flat: magnetization in the tangential-z plane, any angle'),
                    Patch(fc=COLS['out'], ec='k', label='standing: radially outward'),
                    Patch(fc=COLS['in'], ec='k', label='standing: radially inward'),
                    Patch(fc='#c4c4c4', label='main array cube')], loc='lower center', ncol=4, fontsize=10)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(out, dpi=120)
print('->', out)
