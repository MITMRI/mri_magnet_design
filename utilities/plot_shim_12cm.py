"""
Pictures of a shim solution from shim_12cm.py: the field over the DSV before and after shimming, the central planes,
and the shim cubes ring by ring.

    python utilities/plot_shim_12cm.py designs/12cm_bore200/placements.csv designs/12cm_bore200/shim_3.175mm
        [--build N]   (default: the build with the worst unshimmed field)

"After" uses one draw of real shim cubes (their own Br, angle and position errors, as in shim_12cm.py), not ideal cubes.
Writes field_cylinder_b<N>.png (field on a cylinder in the imaging region: side unrolled at true size, end caps),
field_planes_b<N>.png and shim_layout_b<N>.png into the shim folder.
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
    ap.add_argument('--shim-pos-tol', type=float, default=0.1)
    ap.add_argument('--cyl-d', type=float, help='cylinder diameter, mm (default: DSV / sqrt 2, fits in the DSV)')
    ap.add_argument('--cyl-len', type=float, help='cylinder length, mm (default: DSV / sqrt 2)')
    ap.add_argument('--bore', type=float, default=200)
    ap.add_argument('--seed', type=int, default=7)
    a = ap.parse_args()

    summ = pd.read_csv(os.path.join(a.shim_dir, 'summary.csv'))
    b = a.build if a.build is not None else int(summ.groupby('build').unshimmed.first().idxmax())
    df = pd.read_csv(a.csv)
    bpos = df[['x', 'y', 'z']].to_numpy()
    bquat = R.from_euler('z', df.rot_z_deg.to_numpy()[:, None], degrees=True).as_quat()
    npz = np.load(os.path.join(a.shim_dir, f'b{b}_build.npz'))
    bmag = npz['bmag']
    if 'bpos' in npz:
        bpos = npz['bpos']                                  # as-built block positions
    lay = pd.read_csv(os.path.join(a.shim_dir, f'b{b}_layout.csv'))
    rings = pd.read_csv(os.path.join(a.shim_dir, 'rings.csv')) if os.path.exists(
        os.path.join(a.shim_dir, 'rings.csv')) else None
    r = S.Stack(a.cube, a.shim_br, [0], a.bore).r
    phi, psi = np.deg2rad(lay.phi_deg.to_numpy()), np.deg2rad(lay.psi_deg.to_numpy())
    spos = np.c_[r * np.cos(phi), r * np.sin(phi), lay.z.to_numpy()]
    squat = R.from_euler('z', psi[:, None]).as_quat()
    rng = np.random.default_rng(a.seed)
    smag = S.random_magnetization(len(lay), a.shim_br, a.shim_br_tol / 100, a.shim_ang_tol, rng)
    spos = spos + rng.uniform(-a.shim_pos_tol, a.shim_pos_tol, spos.shape)

    def fields(P):
        before = S.cuboid_B(bpos, bquat, bmag, S.BLOCK, P)[:, 0]
        return before, before + S.cuboid_B(spos, squat, smag, (a.cube,) * 3, P)[:, 0]

    # ---- field on a cylinder in the imaging region: side unrolled (no distortion) + both end caps ----
    rs = a.dsv / 2
    rc = a.cyl_d / 2 if a.cyl_d else rs / np.sqrt(2)        # default: largest D = L cylinder inside the DSV
    hc = a.cyl_len / 2 if a.cyl_len else rs / np.sqrt(2)
    s_arc = np.linspace(-np.pi * rc, np.pi * rc, 241)       # arc length around the side, 0 = +x
    zc = np.linspace(-hc, hc, 81)
    SA, ZC = np.meshgrid(s_arc, zc)
    Pside = np.c_[(rc * np.cos(SA / rc)).ravel(), (rc * np.sin(SA / rc)).ravel(), ZC.ravel()]
    gc = np.linspace(-rc, rc, 81)
    CX, CY = np.meshgrid(gc, gc)
    disc = (CX ** 2 + CY ** 2 <= rc ** 2).ravel()
    Pcap = [np.c_[CX.ravel()[disc], CY.ravel()[disc], np.full(disc.sum(), z)] for z in (-hc, hc)]
    n1, n2 = len(Pside), disc.sum()
    Bb, Ba = fields(np.r_[Pside, Pcap[0], Pcap[1]])
    B0, B1 = Bb, Ba                                          # means used by the plane plots too

    def split(B):
        m = ppm_map(B, B.mean())
        caps = []
        for k in range(2):
            c = np.full(CX.size, np.nan)
            c[disc] = m[n1 + k * n2:n1 + (k + 1) * n2]
            caps.append(c.reshape(CX.shape))
        return m[:n1].reshape(SA.shape), caps, np.ptp(m)

    side0, caps0, p0 = split(Bb)
    side1, caps1, p1 = split(Ba)
    lim0 = max(np.abs(side0).max(), *(np.nanmax(np.abs(c)) for c in caps0))
    lim1 = max(np.abs(side1).max(), *(np.nanmax(np.abs(c)) for c in caps1))
    fig = plt.figure(figsize=(18, 7.2), dpi=140, constrained_layout=True)
    gsf = fig.add_gridspec(2, 3, height_ratios=[hc * 2, rc * 2 * 0.75])
    panels = [(side0, caps0, lim0, f'Before: {p0:.0f} ppm peak-to-peak'),
              (side1, caps1, lim0, f'After, same scale: {p1:.0f} ppm'),
              (side1, caps1, lim1, f'After, own scale: {p1:.0f} ppm')]
    off = rc * 1.3
    for col, (side, caps, lim, title) in enumerate(panels):
        ax = fig.add_subplot(gsf[0, col])
        im = ax.pcolormesh(s_arc, zc, side, cmap=CMAP, vmin=-lim, vmax=lim, shading='auto')
        ax.contour(s_arc, zc, side, levels=np.linspace(-lim, lim, 11), colors=INK, linewidths=0.3, alpha=0.4)
        ax.set_aspect('equal')
        ax.set_title(f'{title}\nside, unrolled (true size)', fontsize=10.5, color=INK)
        ticks = np.array([-180, -90, 0, 90, 180])
        ax.set_xticks(np.deg2rad(ticks) * rc)
        ax.set_xticklabels(['-x', '-y', '+x (B0)', '+y', '-x'], fontsize=8)
        ax.set_xlabel(f'around the cylinder ({2 * np.pi * rc:.0f} mm of arc)', fontsize=9, color=MUTED)
        ax.set_ylabel('z (mm)', fontsize=9, color=MUTED)
        ax2 = fig.add_subplot(gsf[1, col])
        for k, (c, lab) in enumerate(zip(caps, (f'-z end cap (z = {-hc:.0f})', f'+z end cap (z = {hc:+.0f})'))):
            xs = gc + (2 * k - 1) * off
            ax2.pcolormesh(xs, gc, c, cmap=CMAP, vmin=-lim, vmax=lim, shading='auto')
            ax2.contour(xs, gc, c, levels=np.linspace(-lim, lim, 11), colors=INK, linewidths=0.3, alpha=0.4)
            ax2.add_patch(Circle(((2 * k - 1) * off, 0), rc, fill=False, ec=MUTED, lw=0.8))
            ax2.text((2 * k - 1) * off, rc + 4, lab, ha='center', fontsize=8.5, color=INK)
            ax2.annotate('', xy=((2 * k - 1) * off + 14, -rc - 5), xytext=((2 * k - 1) * off - 2, -rc - 5),
                         arrowprops=dict(arrowstyle='-|>', lw=0.8, color=MUTED))
            ax2.text((2 * k - 1) * off + 17, -rc - 5, '+x', fontsize=7.5, color=MUTED, va='center')
        ax2.set_xlim(-off - rc - 3, off + rc + 3); ax2.set_ylim(-rc - 10, rc + 9)
        ax2.set_aspect('equal'); ax2.axis('off')
        ax2.text(0, -rc - 5, 'both seen from +z', ha='center', va='center', fontsize=8, color=MUTED)
        fig.colorbar(im, ax=[ax, ax2], shrink=0.7, label='ppm from the cylinder mean')
    fig.suptitle(f'Bx on a {2 * rc:.0f} mm dia x {2 * hc:.0f} mm long cylinder inside the {a.dsv:g} mm DSV, build {b} '
                 f'(blocks +-3% Br, +-1.5 deg, +-0.1 mm), {len(lay)} shim cubes (real cubes, one installation)',
                 fontsize=12, color=INK)
    f1 = os.path.join(a.shim_dir, f'field_cylinder_b{b}.png')
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
            lim = max(np.nanmax(np.abs(M)), 1)
            ax = axs[row, col]
            im = ax.pcolormesh(g, g, M, cmap=CMAP, vmin=-lim, vmax=lim, shading='auto')
            ax.contour(g, g, M, levels=np.linspace(-lim, lim, 11), colors=INK, linewidths=0.3, alpha=0.4)
            ax.add_patch(Circle((0, 0), rs, fill=False, ec=MUTED, lw=1))
            ax.set_aspect('equal')
            xl, yl = name[0], name[1]
            ax.set_xlabel(f'{xl} (mm)', fontsize=9, color=MUTED); ax.set_ylabel(f'{yl} (mm)', fontsize=9, color=MUTED)
            ax.set_title(f'{name}, {label}: {np.nanmax(M) - np.nanmin(M):.0f} ppm in this plane', fontsize=10,
                         color=INK)
            fig.colorbar(im, ax=ax, shrink=0.85, label='ppm from the cylinder mean')
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
    half = a.cube / 2                                       # true size

    def draw_cubes(ax, sel, arrow):
        for (pp, ss) in zip(phi[sel], psi[sel]):
            c = r * np.array([np.cos(pp), np.sin(pp)])
            Rm = np.array([[np.cos(ss), -np.sin(ss)], [np.sin(ss), np.cos(ss)]])
            sq = (Rm @ np.array([[-half, -half], [half, -half], [half, half], [-half, half]]).T).T + c
            nh = (Rm @ np.array([[0, -half], [half, -half], [half, half], [0, half]]).T).T + c
            ax.add_patch(Polygon(sq, facecolor='#d0d7de', edgecolor=BODY, lw=0.6))
            ax.add_patch(Polygon(nh, facecolor=NORTH, edgecolor='none'))
            ax.annotate('', xy=c + Rm @ [arrow, 0], xytext=c, arrowprops=dict(arrowstyle='-|>', lw=0.8, color=INK,
                                                                               mutation_scale=7))

    for k, z in enumerate(zs):
        ax = fig.add_subplot(gs[k // cols, k % cols])
        sel = np.isclose(lay.z, z, atol=0.05)
        ax.add_patch(Circle((0, 0), r, fill=False, ec=MUTED, lw=0.8, ls=(0, (3, 3))))
        ax.add_patch(Circle((0, 0), a.bore / 2, fill=False, ec=MUTED, lw=0.6))
        draw_cubes(ax, sel, half * 5)
        span = '' if rings is None else f'ring z {rings.z_bottom[k]:+.0f}..{rings.z_top[k]:+.0f}, '
        ax.set_title(f'{span}cubes at z {z:+.1f}\n{int(sel.sum())} cubes', fontsize=8.5, color=INK)
        ax.set_xlim(-112, 112); ax.set_ylim(-112, 112); ax.set_aspect('equal'); ax.axis('off')
        ax.annotate('', xy=(30, 0), xytext=(-30, 0), arrowprops=dict(arrowstyle='-|>', lw=0.8, color=MUTED))
        ax.text(0, 6, 'B0', ha='center', fontsize=7, color=MUTED)
    # zoom on the most crowded spot: the cube whose neighbour on its ring is closest
    best = (np.inf, None, None)
    for z in zs:
        idx = np.flatnonzero(np.isclose(lay.z, z, atol=0.05))
        if len(idx) < 2:
            continue
        o = idx[np.argsort(phi[idx])]
        gaps = np.diff(np.r_[phi[o], phi[o[0]] + 2 * np.pi])
        j = int(np.argmin(gaps))
        if gaps[j] < best[0]:
            best = (gaps[j], z, phi[o[j]] + gaps[j] / 2)
    ax = fig.add_subplot(gs[-1, cols - 2:])
    if best[1] is not None:
        sel = np.isclose(lay.z, best[1], atol=0.05)
        ax.add_patch(Circle((0, 0), a.bore / 2, fill=False, ec=MUTED, lw=0.8))
        ax.add_patch(Circle((0, 0), r, fill=False, ec=MUTED, lw=0.6, ls=(0, (3, 3))))
        draw_cubes(ax, sel, half * 2.2)
        cx, cy = r * np.cos(best[2]), r * np.sin(best[2])
        ax.set_xlim(cx - 18, cx + 18); ax.set_ylim(cy - 18, cy + 18); ax.set_aspect('equal')
        ax.set_title(f'True size, most crowded spot (cubes at z {best[1]:+.1f}): closest centres '
                     f'{2 * r * np.sin(best[0] / 2):.2f} mm,\n>= {2 * r * np.sin(best[0] / 2) - a.cube * np.sqrt(2):.1f}'
                     f' mm of resin between cubes at any rotation; solid arc = ring OD', fontsize=8.5, color=INK)
        ax.tick_params(labelsize=7)
        ax.set_xlabel('x (mm)', fontsize=8, color=MUTED); ax.set_ylabel('y (mm)', fontsize=8, color=MUTED)
    ax = fig.add_subplot(gs[-1, :cols - 2])
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
    fig.suptitle(f'Shim cubes, build {b}: {len(lay)} x {a.cube} mm cubes at r {r:.1f} mm (true size; arrows = magnetization, '
                 f'red = north half)', fontsize=12, color=INK)
    fig.tight_layout()
    f3 = os.path.join(a.shim_dir, f'shim_layout_b{b}.png')
    fig.savefig(f3, facecolor='white'); plt.close(fig)
    print(f'build {b}: {p0:.0f} -> {p1:.0f} ppm on the cylinder')
    print(f1, f2, f3, sep='\n')


if __name__ == '__main__':
    main()
