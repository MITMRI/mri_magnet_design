"""
3D views of a magnet design, drawn the way MRI4ALL's magsimulator.plot_magnets3D does it (a magpylib Collection of
Cuboids shown with magpy.show), but with the real 12.7 x 12.7 x 25.4 mm blocks instead of cubes, saved as PNGs.

    python plot_design_3d.py designs/12cm_bore200/placements.csv
"""
import sys

import magpylib as magpy
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial.transform import Rotation as R

matplotlib.use('Agg')

path = sys.argv[1]
out = path.rsplit('/', 1)[0] if '/' in path else '.'
df = pd.read_csv(path)


def collection(rows):
    magnets = magpy.Collection(style_label='magnets')
    for row in rows.itertuples():
        magnets.add(magpy.magnet.Cuboid(magnetization=(1300, 0, 0), dimension=(12.7, 12.7, 25.4),
                                        position=(row.x, row.y, row.z),
                                        orientation=R.from_euler('z', row.rot_z_deg, degrees=True)))
    return magnets


# (rows shown, elevation, azimuth); bore axis = z
views = {'iso': (df, 25, -50), 'end': (df, 90, -90), 'side': (df, 0, -90), 'cutaway': (df[df.y > 0], 20, -90)}
for name, (rows, elev, azim) in views.items():
    fig = plt.figure(figsize=(9, 9), dpi=130)
    ax = fig.add_subplot(projection='3d')
    magpy.show(collection(rows), canvas=ax, backend='matplotlib', style_magnetization_show=True)
    ax.set_box_aspect((1, 1, (df.z.max() - df.z.min() + 25.4) / 325))     # true proportions
    ax.view_init(elev=elev, azim=azim)
    ax.set_title(f'{len(df)} blocks, {name} view (bore along z, B0 along x)')
    fig.savefig(f'{out}/magnets_{name}.png', bbox_inches='tight')
    plt.close(fig)
    print(f'{out}/magnets_{name}.png')
