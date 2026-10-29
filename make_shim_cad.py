# -*- coding: utf-8 -*-
"""
CAD for the slide-in shim sleeve and the shim trays of one shim solution (from shim_design.py).

  sleeve: a tube that slides into the bore of the holder (OD = bore - 2 fit). Its outside has one axial groove per tray,
          closed by the bore wall of the holder, so the grooves are the tray channels. A flange at the +z end sits on
          the outer face of the top end cap: it sets the sleeve's z position and, with a notch at +x, its rotation.
          The flange has clearance holes for the inner (M2) rods and nuts of the holder. The grooves run through the
          flange, so trays slide in and out from the +z end with the magnet fully assembled.
  trays:  one bar per groove that holds cubes, printed for this solution. Pockets open on the outer face (towards the
          holder's bore wall, which keeps the cubes in), each with a push-out hole through the floor. Flat cubes lie in
          a pocket turned to their angle psi, with a V notch at the N face; standing cubes have a round notch if their N
          face points out (towards the bore wall), a square notch if it points in (towards the bore axis). The head at
          the +z end rests on the flange and sets the tray's z position; the tray number is engraved on it.

Geometry options must match the shim_design.py run (cube, trays, bore, fit, tray gap, pocket tolerance); the script
checks the cube centre radius it builds against the one in the solution.

usage:
  python make_shim_cad.py results/shim/FD_robust_sim_b0.csv --holder cad/fd_robust --out cad/fd_robust_shim
Outputs in --out: sleeve.step/.stl, tray_<k>.step/.stl (trays with cubes), shim_assembly.step (sleeve, trays and
shim cubes in place, same coordinates as the holder's assembly.step), tray_loading.png, fit_section.png, trays.csv.
"""

import argparse
import os
import numpy as np
import pandas as pd
import cadquery as cq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('solution', help='shim solution csv from shim_design.py')
    ap.add_argument('--holder', required=True, help='make_ring_cad.py output folder of the holder (slabs.csv, holes.csv)')
    ap.add_argument('--out', default='cad/shim')
    ap.add_argument('--cube', type=float, default=3.175); ap.add_argument('--trays', type=int, default=36)
    ap.add_argument('--bore', type=float, default=100); ap.add_argument('--fit', type=float, default=0.2)
    ap.add_argument('--tray-gap', type=float, default=0.1); ap.add_argument('--pocket-tol', type=float, default=0.2)
    ap.add_argument('--tray-floor', type=float, default=0.8, help='tray floor under the pockets, mm')
    ap.add_argument('--tray-side', type=float, default=0.8, help='tray side wall next to a pocket at its widest, mm')
    ap.add_argument('--slide', type=float, default=0.15, help='clearance between a tray and its groove, each side, mm')
    ap.add_argument('--sleeve-wall', type=float, default=1.2, help='sleeve wall under the grooves, mm')
    ap.add_argument('--endcap', type=float, default=5.0, help='thickness of the end caps (your design), mm')
    ap.add_argument('--flange', type=float, default=3.0, help='flange thickness, mm')
    ap.add_argument('--flange-od', type=float, default=124.0, help='flange outer diameter, mm')
    ap.add_argument('--nut-clear', type=float, default=7.0, help='flange hole diameter over the inner rods and nuts, mm')
    ap.add_argument('--head', type=float, default=6.0, help='tray head height above the flange, mm')
    ap.add_argument('--no-step', action='store_true')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    sol = pd.read_csv(args.solution)
    slabs = pd.read_csv(os.path.join(args.holder, 'slabs.csv'))
    holes = pd.read_csv(os.path.join(args.holder, 'holes.csv'))
    a, n = args.cube, args.trays
    p = a + args.pocket_tol                                   # pocket size
    r_od = args.bore/2 - args.fit                             # sleeve outer radius
    r_face = r_od - args.tray_gap                             # tray outer face
    r_floor = r_face - p - args.tray_floor                    # tray inner face
    r_cube = r_floor + args.tray_floor + a/2
    r_groove = r_floor - args.slide                           # groove bottom
    r_in = r_groove - args.sleeve_wall                        # sleeve inner radius = clear bore / 2
    w = a*np.sqrt(2) + args.pocket_tol + 2*args.tray_side     # tray width: a flat cube may turn to any angle
    gw = w + 2*args.slide
    r_sol = np.hypot(sol['x'], sol['y']).mean()
    assert abs(r_sol - r_cube) < 0.02, f'solution has cubes at r {r_sol:.3f}, this geometry at r {r_cube:.3f}: options differ'
    rib = 2*np.pi*r_groove/n - gw
    z_bot = slabs['bottom'].min() - args.endcap               # outer face of the bottom end cap
    z_top = slabs['top'].max() + args.endcap                  # outer face of the top end cap (flange sits here)
    z_fl = z_top + args.flange                                # flange top = tray stop
    zt0 = sol['z'].min() if len(sol) else 0
    tray_z0 = min(-70 - a - 2, zt0 - a - 2)
    th = 2*np.pi*(np.arange(n) + 0.5)/n
    print(f'sleeve: OD {2*r_od:.2f}, clear bore {2*r_in:.2f} mm, {n} grooves {gw:.2f} wide x {r_od - r_groove:.2f} deep '
          f'(rib {rib:.2f} mm at the groove bottom), z {z_bot:+.2f} to {z_fl:+.2f} mm incl. a {args.flange} mm flange on '
          f'the top end cap ({args.endcap} mm end caps assumed)')
    print(f'trays: {w:.2f} wide x {r_face - r_floor:.2f} thick, z {tray_z0:+.2f} to {z_fl:+.2f} + {args.head} mm head; '
          f'cube centres at r {r_cube:.3f} mm (solution {r_sol:.3f})')
    if rib < 0.8:
        print(f'  WARNING: ribs between grooves only {rib:.2f} mm')

    # inner rods/nuts against the flange; outer rods are outside the flange
    inner = holes[holes['kind'] == 'inner']
    flange_holes = [(x, y) for x, y, r in zip(inner['x'], inner['y'], inner['r']) if r - args.nut_clear/2 < args.flange_od/2]
    outer_r = holes[holes['kind'] == 'outer']['r'].min() if (holes['kind'] == 'outer').any() else np.inf
    if outer_r - 3.5 < args.flange_od/2:
        print(f'  WARNING: flange OD {args.flange_od} reaches the outer rods (r {outer_r:.1f})')

    # ---- sleeve
    L_sl = z_top - z_bot
    sleeve = cq.Workplane('XY').workplane(offset=z_bot).circle(r_od).circle(r_in).extrude(L_sl).val()
    flange = cq.Workplane('XY').workplane(offset=z_top).circle(args.flange_od/2).circle(r_in).extrude(args.flange).val()
    sleeve = sleeve.fuse(flange)
    cut = []
    for t in th:
        g = cq.Solid.makeBox(r_od + 1 - r_groove, gw, z_fl - z_bot + 2, pnt=cq.Vector(r_groove, -gw/2, z_bot - 1))
        cut.append(g.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), np.degrees(t)))
    for x, y in flange_holes:
        cut.append(cq.Solid.makeCylinder(args.nut_clear/2, args.flange + 2, cq.Vector(x, y, z_top - 1)))
    cut.append(cq.Workplane('XY').workplane(offset=z_top - 1)
               .polyline([(args.flange_od/2 - 1.5, 0), (args.flange_od/2 + 1, -2.5), (args.flange_od/2 + 1, 2.5)]).close()
               .extrude(args.flange + 2).val())                                 # +x (B0) notch
    sleeve = sleeve.cut(*cut)
    cq.exporters.export(cq.Workplane().add(sleeve), os.path.join(args.out, 'sleeve.stl'), tolerance=0.01, angularTolerance=0.1)
    if not args.no_step:
        cq.exporters.export(cq.Workplane().add(sleeve), os.path.join(args.out, 'sleeve.step'))

    # ---- trays, built along +x (radial) at angle 0, then turned to their groove
    def tray_solid(k, rows):
        body = cq.Solid.makeBox(r_face - r_floor, w, z_fl - tray_z0, pnt=cq.Vector(r_floor, -w/2, tray_z0))
        head = cq.Solid.makeBox(r_face - (r_groove - 3.0), w, args.head, pnt=cq.Vector(r_groove - 3.0, -w/2, z_fl))
        body = body.fuse(head)
        cutters = []
        for _, c in rows.iterrows():
            zc = c['z']
            pocket = cq.Solid.makeBox(p, p, p + 1, pnt=cq.Vector(r_face - p, -p/2, -p/2))   # local: x radial, y tang, z axial
            # push-out hole through the floor
            pocket = pocket.fuse(cq.Solid.makeCylinder(0.6, args.tray_floor + 2, cq.Vector(r_floor - 1, 0, 0), cq.Vector(1, 0, 0)))
            if c['orientation'] == 'flat':
                # V notch at the N face (+z of the unturned pocket), then the pocket turned by psi about the radial axis
                notch = (cq.Workplane('YZ').workplane(offset=r_face - 1.2)
                         .polyline([(-0.8, p/2 - 0.01), (0.8, p/2 - 0.01), (0, p/2 + 0.8)]).close().extrude(2).val())
                pocket = pocket.fuse(notch).rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), -c['psi_deg'])
            elif c['orientation'] == 'radial_out':
                pocket = pocket.fuse(cq.Solid.makeCylinder(0.7, 2.2, cq.Vector(r_face - 1.2, 0, p/2), cq.Vector(1, 0, 0)))
            else:
                pocket = pocket.fuse(cq.Solid.makeBox(2.2, 1.4, 1.4, pnt=cq.Vector(r_face - 1.2, -0.7, p/2 - 0.7)))
            cutters.append(pocket.translate(cq.Vector(0, 0, zc)))
        body = body.cut(*cutters)
        try:
            txt = (cq.Workplane('XY').workplane(offset=z_fl + args.head - 0.5).center((r_face + r_groove - 3.0)/2, 0)
                   .text(str(k), 3.0, 1.0, kind='bold').val())
            body = body.cut(txt)
        except Exception:
            pass
        return body.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), np.degrees(th[k]))

    trays, rows_out = {}, []
    for k in sorted(sol['tray'].unique()):
        rows = sol[sol['tray'] == k].sort_values('z')
        s = tray_solid(k, rows)
        trays[k] = s
        cq.exporters.export(cq.Workplane().add(s.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), -np.degrees(th[k]))),
                            os.path.join(args.out, f'tray_{k:02d}.stl'), tolerance=0.01, angularTolerance=0.1)
        if not args.no_step:
            cq.exporters.export(cq.Workplane().add(s), os.path.join(args.out, f'tray_{k:02d}.step'))
        for _, c in rows.iterrows():
            rows_out.append(dict(tray=k, groove_angle_deg=np.degrees(th[k]), z=c['z'], from_flange_mm=z_fl - c['z'],
                                 orientation=c['orientation'], psi_deg=c['psi_deg']))
    pd.DataFrame(rows_out).round(3).to_csv(os.path.join(args.out, 'trays.csv'), index=False)
    print(f'{len(trays)} trays with {len(sol)} cubes; empty grooves: {n - len(trays)}')

    # ---- interference check: trays vs sleeve (the sleeve OD clears the bore by --fit by construction)
    worst = 0.0
    for k, s in trays.items():
        v = s.intersect(sleeve).Volume()
        worst = max(worst, v)
    print(f'tray/sleeve overlap: {worst:.4f} mm3 (should be 0)')
    print(f'flange: OD {args.flange_od}, {len(flange_holes)} holes of {args.nut_clear} mm over the inner rods (r {inner["r"].mean():.1f})')

    # ---- assembly: sleeve, trays and shim cubes in the holder's coordinates
    if not args.no_step:
        assy = cq.Assembly(name='shim_sleeve')
        assy.add(sleeve, name='sleeve', color=cq.Color(0.55, 0.5, 0.75, 1))
        ta = cq.Assembly(name='trays')
        for k, s in trays.items():
            ta.add(s, name=f'tray_{k:02d}', color=cq.Color(0.85, 0.8, 0.95, 1))
        assy.add(ta)
        cubes = cq.Assembly(name='shim_cubes')
        for j, c in sol.iterrows():
            t = th[int(c['tray'])]
            cube = cq.Solid.makeBox(a, a, a, pnt=cq.Vector(-a/2, -a/2, -a/2))
            if c['orientation'] == 'flat':
                cube = cube.rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), -c['psi_deg'])
            cube = cube.translate(cq.Vector(r_cube, 0, c['z'])).rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), np.degrees(t))
            cubes.add(cube, name=f'shim_{j}', color=cq.Color(0.8, 0.1, 0.1, 1))
        assy.add(cubes)
        assy.export(os.path.join(args.out, 'shim_assembly.step'))

    # ---- loading guide: every tray seen from outside (the open pocket side), z increasing to the right
    ks = sorted(trays)
    fig, axs = plt.subplots(len(ks), 1, figsize=(16, 0.75*len(ks) + 1.2), squeeze=False)
    for ax, k in zip(axs[:, 0], ks):
        ax.add_patch(plt.Rectangle((tray_z0, -w/2), z_fl - tray_z0, w, fc='#e4def3', ec='#4b3f7a'))
        ax.add_patch(plt.Rectangle((z_fl, -w/2), args.head, w, fc='#cfc6e8', ec='#4b3f7a'))
        for _, c in sol[sol['tray'] == k].iterrows():
            ang = np.deg2rad(c['psi_deg']) if c['orientation'] == 'flat' else 0.0
            ca, sa = np.cos(ang), np.sin(ang)
            P = [(c['z'] + ca*u - sa*v, sa*u + ca*v) for u, v in ((-a/2, -a/2), (a/2, -a/2), (a/2, a/2), (-a/2, a/2))]
            col = {'flat': '#009E73', 'radial_out': '#D55E00', 'radial_in': '#0072B2'}[c['orientation']]
            ax.add_patch(plt.Polygon(P, fc=col, ec='k', lw=0.6))
            if c['orientation'] == 'flat':
                ax.annotate('', (c['z'] + ca*a*0.45, sa*a*0.45), (c['z'] - ca*a*0.45, -sa*a*0.45),
                            arrowprops=dict(arrowstyle='-|>', lw=1.3, color='k', mutation_scale=10))
            else:
                ax.text(c['z'], 0, 'N↑' if c['orientation'] == 'radial_out' else 'N↓', ha='center', va='center', fontsize=6)
        ax.set_xlim(tray_z0 - 2, z_fl + args.head + 2); ax.set_ylim(-w, w); ax.set_aspect('equal')
        ax.set_yticks([]); ax.set_ylabel(f'tray {k}\n{np.degrees(th[k]):.0f}°', rotation=0, ha='right', va='center', fontsize=8)
        ax.tick_params(labelsize=7)
    axs[-1, 0].set_xlabel('z (mm, magnet centre = 0; the head rests on the flange)')
    fig.suptitle(f'Shim tray loading guide ({os.path.basename(args.solution)}), seen from outside the tray (open pockets): '
                 'green = flat, arrow = N; orange = standing, N out (N↑, round notch); blue = standing, N in '
                 '(N↓, square notch)', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'tray_loading.png'), dpi=130)
    plt.close(fig)

    # ---- fit drawing: end view near the bore and a side section through one tray (outlines from the same dimensions)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 8.5), gridspec_kw=dict(width_ratios=(1, 1.25)))
    phi = np.linspace(0, 2*np.pi, 400)
    for R_, st in ((args.bore/2, '-'), (r_od, '-'), (r_in, '-'), (args.flange_od/2, ':')):
        ax1.plot(R_*np.cos(phi), R_*np.sin(phi), color='#4b3f7a' if R_ != args.bore/2 else 'k', ls=st, lw=1)
    for t in th:
        er, et = np.array([np.cos(t), np.sin(t)]), np.array([-np.sin(t), np.cos(t)])
        G = [r_groove*er - gw/2*et, r_od*er - gw/2*et, r_od*er + gw/2*et, r_groove*er + gw/2*et]
        ax1.add_patch(plt.Polygon(G, fc='white', ec='#4b3f7a', lw=0.6))
        k = int(np.argmin(np.abs(th - t)))
        if k in trays:
            T = [r_floor*er - w/2*et, r_face*er - w/2*et, r_face*er + w/2*et, r_floor*er + w/2*et]
            ax1.add_patch(plt.Polygon(T, fc='#cfc6e8', ec='#4b3f7a', lw=0.6))
    for _, h in holes.iterrows():
        ax1.add_patch(plt.Circle((h['x'], h['y']), h['hole_d']/2, fc='k' if h['kind'] == 'inner' else '0.5'))
    for x, y in flange_holes:
        ax1.add_patch(plt.Circle((x, y), args.nut_clear/2, fc='none', ec='#4b3f7a', ls='--', lw=0.8))
    ax1.set_xlim(-62, 62); ax1.set_ylim(-62, 62); ax1.set_aspect('equal')
    ax1.set_title(f'End view near the bore: holder bore (black, d {args.bore:g}), sleeve OD {2*r_od:.1f} / clear bore '
                  f'{2*r_in:.1f} (purple),\n{n} grooves, trays with cubes filled; M2 rods black at r {inner["r"].mean():.0f}, '
                  f'flange holes over them dashed, flange OD dotted', fontsize=9)
    # side section through the busiest tray: holder slabs (grey), end caps (hatched), sleeve, tray, cubes
    kb = int(sol['tray'].value_counts().idxmax())
    for _, s in slabs.iterrows():
        ax2.add_patch(plt.Rectangle((s['bottom'], args.bore/2), s['top'] - s['bottom'], 10, fc='0.85', ec='0.5', lw=0.5))
    sp = os.path.join(args.holder, 'spacers.csv')
    if os.path.exists(sp):
        for _, s in pd.read_csv(sp).iterrows():
            ax2.add_patch(plt.Rectangle((s['bottom'], args.bore/2), s['top'] - s['bottom'], 10, fc='#bcd0ea', ec='0.5', lw=0.5))
    for z0 in (slabs['bottom'].min() - args.endcap, slabs['top'].max()):
        ax2.add_patch(plt.Rectangle((z0, args.bore/2 - 0.0), args.endcap, 10, fc='white', ec='0.3', hatch='///', lw=0.5))
    ax2.add_patch(plt.Rectangle((z_bot, r_in), z_top - z_bot, r_groove - r_in, fc='#b7acd9', ec='#4b3f7a', lw=0.6))
    ax2.add_patch(plt.Rectangle((z_top, r_in), args.flange, args.flange_od/2 - r_in, fc='#b7acd9', ec='#4b3f7a', lw=0.6))
    ax2.add_patch(plt.Rectangle((tray_z0, r_floor), z_fl - tray_z0, r_face - r_floor, fc='#e4def3', ec='#4b3f7a', lw=0.8))
    ax2.add_patch(plt.Rectangle((z_fl, r_groove - 3.0), args.head, r_face - r_groove + 3.0, fc='#cfc6e8', ec='#4b3f7a', lw=0.8))
    for _, c in sol[sol['tray'] == kb].iterrows():
        ax2.add_patch(plt.Rectangle((c['z'] - a/2, r_cube - a/2), a, a, fc='#c0392b', ec='k', lw=0.5))
    ax2.set_xlim(z_bot - 4, z_fl + args.head + 4); ax2.set_ylim(r_in - 4, args.bore/2 + 12)
    ax2.set_xlabel('z (mm)'); ax2.set_ylabel('radius (mm)')
    ax2.set_title(f'Side section through tray {kb} (radius exaggerated): holder slabs grey, FDM spacers blue, end caps '
                  f'hatched ({args.endcap} mm assumed),\nsleeve dark purple with the flange on the top end cap, tray light '
                  f'purple with its head on the flange, shim cubes red', fontsize=9)
    ax2.set_aspect(2.0)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'fit_section.png'), dpi=130)
    plt.close(fig)
    print(f'-> {args.out}')


if __name__ == '__main__':
    main()
