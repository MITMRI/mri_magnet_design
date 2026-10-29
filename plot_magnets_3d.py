# -*- coding: utf-8 -*-
"""
3D views of the magnets in one or more ledgers (the figures/best_magnets3D_views.png layout): top down, side x-z,
side y-z and low angle, drawn with magpylib like magsimulator.plot_magnets3D. With several ledgers each design is
one row, on the same axes, so they can be compared side by side.

usage:
    python plot_magnets_3d.py out.png ledger.xlsx[:label] [ledger2.xlsx[:label] ...]
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import magpylib as magpy
import magsimulator

VIEWS = [('Top down (looking along the bore, z axis)', 90, -90),
         ('Side view (x-z plane)', 0, -90),
         ('Side view (y-z plane)', 0, 0),
         ('Low angle', 10, -60)]


def main():
    out, specs = sys.argv[1], sys.argv[2:]
    designs = []
    for s in specs:
        path, _, label = s.partition(':')
        L = pd.read_excel(path)
        designs.append((L, label or path.split('/')[-1]))
    lim = max(np.abs(L[['X-pos', 'Y-pos', 'Z-pos']].values).max() + L['Magnet_length'].iloc[0] for L, _ in designs)

    fig = plt.figure(figsize=(5*len(VIEWS), 5.2*len(designs)))
    for i, (L, label) in enumerate(designs):
        col = magsimulator.build_magnet_collection(L, [1320, 0, 0])
        for j, (title, elev, azim) in enumerate(VIEWS):
            ax = fig.add_subplot(len(designs), len(VIEWS), i*len(VIEWS) + j + 1, projection='3d')
            magpy.show(col, canvas=ax, backend='matplotlib')
            ax.view_init(elev=elev, azim=azim)
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
            ax.set_box_aspect((1, 1, 1))
            ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)'); ax.set_zlabel('z (mm)')
            leg = ax.get_legend()
            if leg:
                leg.remove()
            ax.set_title(title, fontsize=11)
        fig.text(0.01, 1 - (i + 0.5)/len(designs), f'{label}\n{len(L)} magnets', rotation=90,
                 va='center', ha='left', fontsize=14, weight='bold')
    fig.tight_layout(rect=(0.03, 0, 1, 1))
    fig.savefig(out, dpi=110)
    print('->', out)


if __name__ == '__main__':
    main()
