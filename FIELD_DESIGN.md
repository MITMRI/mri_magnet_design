# Field-driven design: target field -> LP/MILP placement -> gradient refinement

`field_design.py` designs the aligned-ring Halbach magnet the way CoilGen designs coils: start from the field you
want, solve for the source distribution that makes it, turn that into buildable parts, then refine. It replaces the
GA + Powell + minimax pipeline in [DESIGN_NOTES.md](DESIGN_NOTES.md) as the design method going forward: in about 3 minutes it
reaches 45 ppm (ideal magnets) where that pipeline reached 430 ppm after hours.

```bash
source .venv/bin/activate
python field_design.py results/field_design/FD_gd.xlsx      # ~3 min: LP 2 s, MILP 2 min, gradient stages ~1 min
```

## Result

`results/field_design/FD_gd.xlsx` (ledger format, same as the other designs)

| | Field-driven (`FD_gd.xlsx`) | Previous BEST |
|---|---|---|
| Magnets | 650 x 1/4" N42 cubes | 649 |
| B0 | 45.00 mT | 45.00 mT |
| Homogeneity, ideal magnets, 40 x 40 mm cylinder | **45 ppm** (44 ppm in the independent MATLAB model, 51,537 points) | 430 ppm |
| 40 mm sphere | 35 ppm | 302 ppm |
| As built, p90 (1% Br, 0.05 mm, 1 degree, 400 builds, no shimming) | 4,644 ppm | 4,423 ppm |
| Length of magnets / holder OD | 148.0 mm / 174 mm | 143.9 mm / 184 mm |
| Ring slices | 6 pairs (no center ring) | 4 pairs + center |
| Required iHc at 40 C, 1.2x safety | 953 kA/m: N42 (955) passes, with almost no margin | 941 kA/m |
| Run time | ~3 min | hours (GA) + refinement |

Ring layout (each row is a ring slice, present at +z and -z; `n @ r` = magnets at center radius r in mm):

| Slice z (mm) | Inner layer | Middle layer | Outer layer |
|---|---|---|---|
| +-7.9 | 26 @ 57.5 | 26 @ 68.5 | - |
| +-18.9 | 13 @ 59.0 | - | - |
| +-29.3 | 32 @ 57.5 | 32 @ 68.5 | - |
| +-52.1 | 24 @ 58.2 | 16 @ 69.2 | - |
| +-61.5 | 32 @ 57.5 | 39 @ 68.5 | - |
| +-70.8 | 30 @ 57.5 | 37 @ 68.5 | 18 @ 79.5 |

**As built it is not better.** Below ~1,000 ppm on paper, magnet errors dominate (mostly each cube's magnetization
angle), and every design lands at ~4,400-4,700 ppm without shimming. The method's value is getting to the paper
optimum in minutes; reaching better as-built numbers needs measured magnets or shimming, or an as-built objective in
the gradient stage (see Next steps).

Figures: `figures/compare_magnets3D_views.png` (3D magnet views, previous BEST on top) and
`figures/compare/fig1-5_*.png` (MATLAB field figures side by side). The ppm maps (`fig3`) use separate color
scales: about -190 to +210 ppm for the previous BEST, -20 to +13 ppm for this design.

## How it works

1. **Target field.** Bx = B0 (min B0 + 0.1 mT) on check points covering the target cylinder (`--sensors 400`; only
   the z >= 0 half is used since designs are mirror symmetric). Bx is harmonic in the bore, so its extremes over
   the region are on the surface and checking the surface is enough.
2. **Sensitivity matrix.** For every candidate ring (each layer radius x z on a 1 mm grid), the exact cuboid field
   one magnet in that ring adds at every check point (a full ring computed with magpylib, divided by its magnet
   count). These are the magnet equivalent of CoilGen's stream-function basis functions. Layer radii are fixed by
   the bore and spacing rules: 57.5, 68.5, 79.5 mm for 3 layers.
3. **Continuous source (LP).** Magnet density along z for each layer that minimizes the peak field deviation, with
   density capped by the ring pitch (one full ring per 9.35 mm) and the total magnet budget. It is convex, so this
   is the global optimum of the relaxed problem: 59 ppm with 650 magnets. Saved as `<out>_source.npz`.
4. **Discrete placement (MILP).** Whole ring slices on a z grid (`--milp-dz 1`), all layers of a slice at the same
   z, slices >= 9.35 mm apart, an integer magnet count per ring (0, or `--nmin 8` up to what fits), total <=
   budget, minimizing the same peak deviation. HiGHS, capped at `--milp-time 120` s. This alone gave 64 ppm (exact
   field) and is saved as `<out>_initial.xlsx`.
5. **Gradient refinement.** Variables: slice z positions, each ring's radius and rotation offset. Gradients are
   analytic for point dipoles. The difference between dipole and exact cuboid field (~20 ppm) is added as a
   correction refreshed every 25 steps, so the objective is exact while the gradient comes from the dipole formula
   (checked against finite differences of the exact field). Stages:
   - Adam on a smooth peak-to-peak (log-sum-exp, sharpened over the run) with penalties for constraints (~6 s).
   - SLSQP on the exact minimax form: min t subject to |deviation_i| <= t, B0 >= min, geometry (~2 s).
   - Screened magnet-count changes: +-1 per ring and moving one magnet between rings, each re-tuned by SLSQP;
     kept only if the exact field improves (~40 s).

Constraints are the same as `refine_minimax.py`: <= 650 magnets, B0 >= 45 mT, magnets <= 150 mm long, 3 mm wall to
the 100 mm bore, holder OD <= 200 mm, slices >= cube + 3 mm apart, layers >= sqrt(2)*cube + 2 mm apart, magnets fit
around each ring (3D diagonal rule from `magsimulator`). Each is a command-line option; see `python field_design.py -h`.
`--time-limit` (default 600 s) stops any stage that is still running.

## Why the magnets end up at the ends

On the axis, a ring near the middle (|z| < ~25 mm) makes a field that peaks at the center and falls towards the
ends of the region; a ring near the end (|z| > ~30 mm) makes a field that is weak at the center but rises towards the
ends of the region (inner layer, change from z = 0 to z = 20 mm: -20% for a ring at 10 mm, +41% at 70 mm). A flat field needs the two to
cancel. End magnets are ~4x weaker at the center but bend the field strongly, so many are needed: the LP puts
~60% of the magnets beyond z = 45 mm and leaves 35-50 mm empty. The GA never made that jump from its evenly spread
starting layouts.

Radially the opposite holds: an inner-layer magnet adds ~2.5x more field at the center than an outer-layer one, so
the budget left after flattening goes to the inner layers (18 outer-layer magnets per side vs 66 in the previous
BEST), which is how the design still reaches 45 mT.

## GA comparison

Same constraints and evaluation, all 3 layers, 40 x 40 mm cylinder:

| Method | ppm, ideal (20k points) | As built p90 | Wall time |
|---|---|---|---|
| `optimize_robust_aligned.py` GA from scratch (paper objective, 11 processes) | 3,669 (at 47.6 mT) | 6,886 | 10 min |
| LP + MILP initial placement only | 64 | 4,716 | 2 min |
| Initial placement + GA refinement (`--ga --pool 11`, same variables plus counts) | 189 | 4,619 | 19 min* |
| **Initial placement + gradient refinement** | **45** | 4,644 | **3 min** |

\* The GA was given 10 min but pymoo's time termination overran to 19 min. It also ended worse than its 64 ppm
start: the reloaded start broke a constraint by a rounding error (radii are stored to 0.01 mm), so the GA first had
to find valid designs. The gradient stages handle the same start in 95 s (`FD_gdonly.log`).
Logs and ledgers: `results/field_design/`.

```bash
python field_design.py out.xlsx --start results/field_design/FD_gd_initial.xlsx --ga --pool 11   # GA refinement
python field_design.py out.xlsx --start <any ledger.xlsx>                                       # gradient refinement only
TIME_LIMIT=00:10:00 POOL=11 ALIGN=1 LAYER_WALL=2 python optimize_robust_aligned.py 40 100000 650 45 b0 150 100 40 3 paper
```

## Checking a new design

The optimizer only sees ~200 check points. Small rings (13 or 16 magnets) make azimuthal ripple that sparse points
can miss, so always check on a dense set, and then the rest of the toolchain:

```bash
LAYOUT=$PWD/results/field_design/FD_gd.xlsx OUTDIR=$PWD/verify/FD_gd matlab -batch "cd verify; simulate_halbach_field"
python plot_magnets_3d.py figures/compare_magnets3D_views.png results/BEST_bore100_cyl40x40_649mag_427ppm.xlsx:"Previous BEST" results/field_design/FD_gd.xlsx:"Field-driven"
python check_demag.py results/field_design/FD_gd.xlsx --Br 1320 --max-temp 40
python make_ring_cad.py results/field_design/FD_gd.xlsx --out cad/fd_gd --magnets     # not yet run for this design
```

The MATLAB script reports homogeneity on a 51,537-point grid with an independent cuboid model; the Python numbers
above use 10-20k points on the cylinder surface (`results/field_design/comparison.txt`).

## Notes and limits

- The MILP treats each ring's field as proportional to its magnet count. That is exact for large rings, but small
  rings have azimuthal ripple and depend on their rotation, so the exact field after the MILP can differ from its
  linear estimate (61 vs 74 ppm here; 466 vs 87 ppm in a test run with 8-magnet rings). The gradient stages fix
  this; `--nmin` raises the smallest ring allowed.
- The MILP usually hits its time limit before proving optimality; a longer `--milp-time` can give a better start.
- Phase (ring rotation) only matters for small rings; for n >= ~15 the field is nearly independent of it.
- Holder CAD (`make_ring_cad.py`) has not been generated or checked for this layout yet. The rings are thinner
  radially (mostly 2 layers), which should make printing easier.

## Next steps

1. Optimize the as-built metric in the gradient stage: `magsimulator.tolerance_homogeneity` already has the
   per-magnet error sensitivities, so the expected spread from magnet errors can be added to the objective.
2. Assign measured magnets (Br, angle) to positions with the same gradient machinery.
3. Generate and check the ring CAD for `FD_gd.xlsx`.
