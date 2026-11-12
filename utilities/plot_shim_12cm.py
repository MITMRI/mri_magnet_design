"""
Pictures of a shim solution from shim_12cm.py: the field over the DSV before and after shimming, the central planes,
and the shim cubes ring by ring.

    python utilities/plot_shim_12cm.py designs/12cm_bore200/placements.csv designs/12cm_bore200/shim_3.175mm_plates
        [--build N]   (default: the build with the worst unshimmed field)

"After" uses one draw of real shim cubes (their own Br and angle errors, as in shim_12cm.py), not ideal cubes.
Writes field_sphere_b<N>.png, field_planes_b<N>.png and shim_layout_b<N>.png into the shim folder.
"""
import argparse
import os

import magpylib as magpy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, Polygon
from scipy.spatial.transform import Rotation as R

import shim_12cm as S

matplotlib.use('Agg')
INK, MUTED, NORTH, BODY = '#1f2328', '#6e7781', '#c4314b', '#57606a'
CMAP = 'RdBu_r'                                            # diverging, neutral at 0 ppm


def ppm_map(B, mean):
    return 1e6 * (B - mean) / mean


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('csv')
    ap.add_argument('shim_dir')
    ap.add_argument('--build', type=int)
    ap.add_argument('--dsv', type=float, default=120)
    ap.add_argument('--cube', type=float, default=3.175)
    ap.add_argument('--shim-br', type=float, default=1300)
    ap.add_argument('--shim-br-tol', type=float, default=3.0)
    ap.add_argument('--shim-ang-tol', type=float, default=1.5)
    ap.add_argument('--bore', type=float, default=200)
    ap.add_argument('--seed', type=int, default=7)
    a = ap.parse_args()

    summ = pd.read_csv(os.path.join(a.shim_dir, 'summary.csv'))
    b = a.build if a.build is not None else int(summ.groupby('build').unshimmed.first().idxmax())
    df = pd.read_csv(a.csv)
    bpos = df[['x', 'y', 'z']].to_numpy()
    bquat = R.from_euler('z', df.rot_z_deg.to_numpy()[:, None], degrees=True).as_quat()
    bmag = np.load(os.path.join(a.shim_dir, f'b{b}_build.npz'))['bmag']
    lay = pd.read_csv(os.path.join(a.shim_dir, f'b{b}_layout.csv'))
    rings = pd.read_csv(os.path.join(a.shim_dir, 'rings.csv')) if os.path.exists(
        os.path.join(a.shim_dir, 'rings.csv')) else None
    r = S.Stack(a.cube, a.shim_br, [0], a.bore).r
    phi, psi = np.deg2rad(lay.phi_deg.to_numpy()), np.deg2rad(lay.psi_deg.to_numpy())
    spos = np.c_[r * np.cos(phi), r * np.sin(phi), lay.z.to_numpy()]
    squat = R.from_euler('z', psi[:, None]).as_quat()
    smag = S.random_magnetization(len(lay), a.shim_br, a.shim_br_tol / 100, a.shim_ang_tol,
                                  np.random.default_rng(a.seed))

    def fields(P):
        before = S.cuboid_B(bpos, bquat, bmag, S.BLOCK, P)[:, 0]
        return before, before + S.cuboid_B(spos, squat, smag, (a.cube,) * 3, P)[:, 0]

    # ---- field over the DSV surface (equirectangular map) ----
    rs = a.dsv / 2
    lon, lat = np.linspace(-180, 180, 181), np.linspace(-90, 90, 91)
    LO, LA = np.meshgrid(np.deg2rad(lon), np.deg2rad(lat))
    P = rs * np.c_[(np.cos(LA) * np.cos(LO)).ravel(), (np.cos(LA) * np.sin(LO)).ravel(), np.sin(LA).ravel()]
    B0, B1 = fields(P)
    m0, m1 = ppm_map(B0, B0.mean()).reshape(LA.shape), ppm_map(B1, B1.mean()).reshape(LA.shape)
    lim0 = np.abs(m0).max()
    lim1 = np.abs(m1).max()
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6), dpi=140, constrained_layout=True)
    panels = [(m0, lim0, f'Before: {np.ptp(m0):.0f} ppm peak-to-peak'),
              (m1, lim0, f'After, same scale: {np.ptp(m1):.0f} ppm'),
              (m1, lim1, f'After, own scale: {np.ptp(m1):.0f} ppm')]
    for ax, (m, lim, title) in zip(axs, panels):
        im = ax.pcolormesh(lon, lat, m, cmap=CMAP, vmin=-lim, vmax=lim, shading='auto')
        ax.contour(lon, lat, m, levels=np.linspace(-lim, lim, 11), colors=INK, linewidths=0.3, alpha=0.4)
        ax.set_title(title, fontsize=11, color=INK)
        ax.set_xlabel('longitude about z (deg; 0 = +x = B0)', fontsize=9, color=MUTED)
        ax.set_ylabel('latitude (deg; +90 = +z)', fontsize=9, color=MUTED)
        ax.set_aspect('equal')
        fig.colorbar(im, ax=ax, shrink=0.8, label='ppm from the mean')
    fig.suptitle(f'Bx over the {a.dsv:g} mm DSV surface, build {b} (blocks +-3% Br, +-1.5 deg), {len(lay)} shim cubes '
                 f'(real cubes, one installation)', fontsize=12, color=INK)
    f1 = os.path.join(a.shim_dir, f'field_sphere_b{b}.png')
    fig.savefig(f1, facecolor='white'); plt.close(fig)

    # ---- central planes ----
    g = np.linspace(-rs, rs, 121)
    U, V = np.meshgrid(g, g)
    inside = U ** 2 + V ** 2 <= rs ** 2
    planes = {'xz plane (y = 0)': np.c_[U.ravel(), 0 * U.ravel(), V.ravel()],
              'xy plane (z = 0)': np.c_[U.ravel(), V.ravel(), 0 * U.ravel()]}
    fig, axs = plt.subplots(2, 2, figsize=(10.5, 9.5), dpi=140, constrained_layout=True)
    for row, (name, Pp) in enumerate(planes.items()):
        Bb, Ba = fields(Pp[inside.ravel()])
        for col, (Bv, mean, label) in enumerate([(Bb, B0.mean(), 'before'), (Ba, B1.mean(), 'after')]):
            M = np.full(U.shape, np.nan)
            M[inside] = ppm_map(Bv, mean)
            lim = lim0 if col == 0 else max(np.nanmax(np.abs(M)), 1)
            ax = axs[row, col]
            im = ax.pcolormesh(g, g, M, cmap=CMAP, vmin=-lim, vmax=lim, shading='auto')
            ax.contour(g, g, M, levels=np.linspace(-lim, lim, 11), colors=INK, linewidths=0.3, alpha=0.4)
            ax.add_patch(Circle((0, 0), rs, fill=False, ec=MUTED, lw=1))
            ax.set_aspect('equal')
            xl, yl = name[0], name[1]
            ax.set_xlabel(f'{xl} (mm)', fontsize=9, color=MUTED); ax.set_ylabel(f'{yl} (mm)', fontsize=9, color=MUTED)
            ax.set_title(f'{name}, {label}: {np.nanmax(M) - np.nanmin(M):.0f} ppm in this plane', fontsize=10,
                         color=INK)
            fig.colorbar(im, ax=ax, shrink=0.85, label='ppm from the DSV mean')
    fig.suptitle(f'Bx in the central planes, build {b} (each panel on its own symmetric scale)', fontsize=12,
                 color=INK)
    f2 = os.path.join(a.shim_dir, f'field_planes_b{b}.png')
    fig.savefig(f2, facecolor='white'); plt.close(fig)

    # ---- shim layout ----
    zs = np.sort(lay.z.unique()) if rings is None else rings.cube_plane_z.to_numpy()
    n = len(zs)
    cols = 6
    fig = plt.figure(figsize=(3.0 * cols, 3.2 * int(np.ceil(n / cols)) + 4.2), dpi=140)
    gs = fig.add_gridspec(int(np.ceil(n / cols)) + 1, cols, height_ratios=[1] * int(np.ceil(n / cols)) + [1.35])
    half = a.cube / 2 * 4                                   # cubes drawn 4x size so they are visible
    for k, z in enumerate(zs):
        ax = fig.add_subplot(gs[k // cols, k % cols])
        sel = np.isclose(lay.z, z, atol=0.05)
        ax.add_patch(Circle((0, 0), r, fill=False, ec=MUTED, lw=0.8, ls=(0, (3, 3))))
        ax.add_patch(Circle((0, 0), a.bore / 2, fill=False, ec=MUTED, lw=0.6))
        for (pp, ss) in zip(phi[sel], psi[sel]):
            c = r * np.array([np.cos(pp), np.sin(pp)])
            Rm = np.array([[np.cos(ss), -np.sin(ss)], [np.sin(ss), np.cos(ss)]])
            sq = (Rm @ np.array([[-half, -half], [half, -half], [half, half], [-half, half]]).T).T + c
            nh = (Rm @ np.array([[0, -half], [half, -half], [half, half], [0, half]]).T).T + c
            ax.add_patch(Polygon(sq, facecolor='#d0d7de', edgecolor=BODY, lw=0.6))
            ax.add_patch(Polygon(nh, facecolor=NORTH, edgecolor='none'))
            ax.annotate('', xy=c + Rm @ [half * 2.2, 0], xytext=c, arrowprops=dict(arrowstyle='-|>', lw=0.8,
                                                                                    color=INK, mutation_scale=7))
        span = '' if rings is None else f'ring z {rings.z_bottom[k]:+.0f}..{rings.z_top[k]:+.0f}, '
        ax.set_title(f'{span}cubes at z {z:+.1f}\n{int(sel.sum())} cubes', fontsize=8.5, color=INK)
        ax.set_xlim(-112, 112); ax.set_ylim(-112, 112); ax.set_aspect('equal'); ax.axis('off')
        ax.annotate('', xy=(30, 0), xytext=(-30, 0), arrowprops=dict(arrowstyle='-|>', lw=0.8, color=MUTED))
        ax.text(0, 6, 'B0', ha='center', fontsize=7, color=MUTED)
    ax = fig.add_subplot(gs[-1, :])
    ax.scatter(np.degrees(phi) % 360, lay.z, s=40, marker='s', facecolor='#d0d7de', edgecolor=BODY, zorder=2)
    ax.quiver(np.degrees(phi) % 360, lay.z, np.cos(psi - phi), np.sin(psi - phi), angles='uv', scale=70, width=0.0022,
              color=NORTH, zorder=3)
    for z in zs:
        ax.axhline(z, color=MUTED, lw=0.4, ls=(0, (3, 3)))
    if rings is not None:
        for zb in np.r_[rings.z_bottom, rings.z_top.iloc[-1]]:
            ax.axhline(zb, color=INK, lw=0.6)
    ax.set_xlim(0, 360); ax.set_xticks(range(0, 361, 45))
    ax.set_xlabel('angle around the bore (deg, 0 = +x)', fontsize=9, color=MUTED)
    ax.set_ylabel('z (mm)', fontsize=9, color=MUTED)
    ax.set_title('Unrolled: every cube (arrow = magnetization relative to the radial direction: right = outward, '
                 'up = +tangential); solid lines = ring boundaries, dashed = cube planes', fontsize=9, color=INK)
    fig.suptitle(f'Shim cubes, build {b}: {len(lay)} x {a.cube} mm cubes at r {r:.1f} mm (squares drawn 4x size, '
                 f'red = north half)', fontsize=12, color=INK)
    fig.tight_layout()
    f3 = os.path.join(a.shim_dir, f'shim_layout_b{b}.png')
    fig.savefig(f3, facecolor='white'); plt.close(fig)
    print(f'build {b}: {np.ptp(m0):.0f} -> {np.ptp(m1):.0f} ppm on the map grid')
    print(f1, f2, f3, sep='\n')


if __name__ == '__main__':
    main()
