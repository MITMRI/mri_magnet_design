# 12cm_bore200: 200 mm bore, ~50 ppm, 96.0 mT

| | |
|---|---|
| B0 | **96.01 mT** (Br = 1.30 T; scales linearly with Br) |
| Homogeneity | 50.07 ppm peak-to-peak of Bx over the 120 mm sphere (50,000 surface points; 50.0 on 20,000), ideal magnets |
| Magnets | **970** x KJ BX088SH (N42SH, 25.4 x 12.7 x 12.7 mm) |
| Magnet mass | 29.8 kg (30.73 g each) |
| Bore ID (plastic) | **201.1 mm** (closest block corner at r = 102.54 mm, 2 mm wall) |
| OD (plastic) | **325.0 mm** (outermost block corner at r = 160.50 mm, 2 mm wall) |
| Magnet length | **357.8 mm** end to end (limit 360) |
| Min gaps | 1.87 mm within a ring, 2.17 mm between rings of a slice, 3.00 mm between slices; no overlaps ([gap_check.txt](gap_check.txt)) |

Orientation: each block's 25.4 mm side along the bore (z), magnetized through a 12.7 mm face, in the ring plane at
angle 2*theta for a block at azimuth theta; B0 along x.

## Layout (6 slice pairs at +-z; rings: n blocks @ radius of block centers, mm)

| Slice z | Inner | Middle | Outer |
|---|---|---|---|
| +-14.2 | | 23 @ 131.52 | 48 @ 151.52 |
| +-44.2 | | 31 @ 131.52 | 48 @ 151.52 |
| +-76.2 | | 21 @ 121.52 | 47 @ 151.52 |
| +-106.2 | | 18 @ 131.52 | 46 @ 151.52 |
| +-136.2 | | 39 @ 121.52 | 48 @ 151.52 |
| +-166.2 | 35 @ 111.52 | 38 @ 131.52 | 43 @ 151.52 |

Every block: [placements.csv](placements.csv) (x, y, z of the block center in mm; rot_z_deg = magnetization angle from
+x). Specs as JSON: [design.json](design.json). Optimizer output: [run.log](run.log).

## Settings that produced it

```bash
python field_design_12cm.py results/12cm/m1_dr5_dz2.xlsx --bore 200 --bore-wall 2 --od-wall 2 --ppm 50 --margin 1 \
    --dr-grid 5 --dz 2 --milp-time 150 --time-limit 360
```

Gurobi 13.0.3 (academic license), MILP proven optimal for this grid (gap 0.00%) in 57 s on 12 cores. Defaults
otherwise: 1200 sphere points (quarter used), ripple 2e-6, slice wall 3 mm, layer wall 2 mm, ring wall 1.5 mm, max OD
325 mm, max length 360 mm, Br 1300 mT.

## Notes

- B0 vs ppm is a staircase, not a slope: on this grid the optimum is ~60 mT up to ~42 ppm and jumps to ~88-96 mT from
  ~46 ppm, because above that point the inner ring layer can be mostly filled. Small changes to the ppm target near
  the step change B0 a lot.
- Ideal magnets only. As built, magnet errors put real builds at thousands of ppm before shimming.
