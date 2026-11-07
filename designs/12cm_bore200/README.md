# 12cm_bore200: 200 mm bore, 50 ppm, 58.1 mT

| | |
|---|---|
| B0 | 58.11 mT (Br = 1.30 T; scales linearly with Br) |
| Homogeneity | 46.7 ppm peak-to-peak of Bx over the 120 mm sphere (20,000 surface points, ideal magnets) |
| Magnets | 610 x KJ BX088SH (N42SH, 25.4 x 12.7 x 12.7 mm), ~18.7 kg |
| Bore ID (plastic) | 200.0 mm (magnet inner faces at r = 102.0 mm, 2 mm wall) |
| OD (plastic) | 324.9 mm (magnet outer edge r = 160.45 mm, 2 mm wall) |
| Magnet length | 343.2 mm end to end |
| Min gaps | 1.69 mm within a ring, 3.00 mm between slices, no overlaps ([gap_check.txt](gap_check.txt)) |

Orientation: each block's 25.4 mm side along the bore (z), magnetized through a 12.7 mm face, in the ring plane at
angle 2*theta for a block at azimuth theta; B0 along x.

## Layout (slices at +-z, rings: n blocks @ radius of block centers, mm)

| Slice z | Rings |
|---|---|
| +-14.20 | 43 @ 147.55 |
| +-43.62 | 48 @ 148.77 |
| +-72.02 | 48 @ 148.77 |
| +-100.42 | 42 @ 141.76 |
| +-128.82 | 42 @ 151.48 |
| +-158.88 | 35 @ 110.98, 47 @ 145.68 |

Every block: [placements.csv](placements.csv) (x, y, z of the block center in mm, rot_z_deg = magnetization angle
from +x). Specs as JSON: [design.json](design.json). Full optimizer log: [run.log](run.log).

## Settings that produced it

```bash
python field_design_12cm.py results/12cm/ppm50/dr10_dz3.xlsx --bore 200 --bore-wall 2 --od-wall 2 --ppm 50 \
    --dr-grid 10 --dz 3 --milp-time 220 --time-limit 330
```

Defaults otherwise: margin 8 ppm, 1200 sphere points (quarter used), ripple 2e-6, slice wall 3 mm, layer wall 2 mm,
ring wall 1.5 mm, max OD 325 mm, max length 360 mm, Br 1300 mT.

## Known limitation

This is far below what the geometry allows. The continuous (LP) upper bound for this envelope at 50 ppm is ~168 mT;
the integer ring-selection MILP on a coarse grid only reaches ~50-58 mT (and ~81-95 mT at 60 ppm). The gap is the
optimizer (coarse grid, whole-ring choices, single-threaded HiGHS hitting time limits), not physics. Better designs
should come from a stronger formulation/solver (e.g. Gurobi, finer grids, continuous design then rounding).
