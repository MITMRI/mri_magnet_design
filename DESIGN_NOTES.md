# Design notes: compact Halbach magnet (100 mm bore, 40 x 40 mm imaging region)

This fork adds tools to the MRI4ALL magnet scripts for a smaller magnet built from 1/4" N42 cubes in resin-printed rings,
and records the design study that produced the recommended layout. All numbers are simulated with magpylib 4.

## Recommended design

`results/BEST_bore100_cyl40x40_649mag_427ppm.xlsx` (ledger format; readable by `magcadexporter.py` and `check_demag.py`)

| | |
|---|---|
| Magnets | 649 x 6.35 mm (1/4") N42 cubes, Br 1.32 T |
| B0 | 45.0 mT (1.92 MHz proton) |
| Magnet bore (holder inner diameter) | 100 mm, with a 3 mm wall to the nearest magnet corner |
| Holder outer diameter | 184 mm (largest ring); rings print upright on a Formlabs Form 4 |
| Length | 144 mm of magnets, 9 ring slices, mirror symmetric in z |
| Layers | 3, all layers of a ring share the same z (one stack of whole rings) |
| Target region | 40 mm diameter x 40 mm long cylinder |
| Homogeneity, ideal magnets | **427 ppm** on the cylinder (302 ppm on a 40 mm sphere) |
| Homogeneity, as built | 3,222 ppm typical, 4,413 ppm worst case (90th percentile), no shimming |
| Demagnetization (N42, iHc 955 kA/m) | corner damage from ~66 C; 1.2x safety margin up to ~42 C |

"As built" assumes random errors on every magnet: 1% Br spread, 0.05 mm placement and 1 degree magnetization
angle (sigma), evaluated at 400 random builds with `magsimulator.tolerance_homogeneity`.

Ring layout (each row is one printed ring, present at +z and -z):

| Ring z (mm) | Inner layer | Middle layer | Outer layer |
|---|---|---|---|
| 0 | 23 @ r 57.5 | 22 @ 68.5 | 6 @ 84.6 |
| +-18.9 | 28 @ 57.7 | 39 @ 68.7 | 23 @ 79.7 |
| +-44.7 | 31 @ 58.2 | 39 @ 69.2 | 22 @ 80.2 |
| +-59.4 | 32 @ 57.5 | 39 @ 68.5 | 15 @ 79.5 |
| +-68.8 | 31 @ 57.5 | - | - |

r is the radius of the magnet centers; each ring places its pockets at its own radii. Magnets follow the Halbach
rule: a magnet at azimuth theta is rotated 2*theta about z, magnetization along its local x.

Figures: `figures/best_fieldmap.png` (field slices), `figures/best_magnets3D_views.png` (top, two sides, 3D), and the
same for the previous design (`final_*`). In the field slice plots from `extract_3Dfields`, the panel titled
"xz plane" actually shows the yz plane and the one titled "yz plane" shows the xz plane (see Known issues).

## Requirements that set the design

- Imaging region 40 x 40 mm cylinder inside a 40 mm coil bore; 100 mm clear magnet bore leaves 30 mm per side for coils.
- B0 >= 45 mT, at most 650 magnets, total length <= 150 mm.
- Whole printed rings only (no split segments), which must fit the Form 4 (200 x 125 x 210 mm) standing upright.
- Ring slices >= 9.35 mm apart (6.35 mm cube + 3 mm resin), >= 2 mm resin between layers.

## What was learned

**Paper vs as built.** Below about 2,000 ppm on paper, magnet errors dominate. For this magnet size they add roughly
3,000-4,000 ppm regardless of design, mostly from the magnetization angle of each cube (a property of the magnets, not the
printer) and placement. Going from 2,057 to 427 ppm on paper improved the as-built worst case only from 4,823 to 4,413 ppm.
Reaching < 1,000 ppm as built requires measuring the magnets and assigning positions/rotations, or shimming.

**Model check against MRI4ALL.** Their saved design (`optimization_after_neonate_...maxmag990.xlsx`) gives 3,843 ppm on
paper over its 140 mm sphere and 738 ppm over a 100 mm sphere (their quoted ~740). With the same error model it
predicts ~6,600 ppm as built over 140 mm, close to their reported ~6,000 ppm measurement. Scaled to a 100 mm bore it scores
1,106 ppm on the 40 x 40 cylinder; the recommended design reaches 427 ppm with fewer magnets.

**Geometry.** Homogeneity depends on region size relative to the magnet radius and on magnet length relative to radius.
A 50 mm region in a 100 mm bore could not get below ~8,000 ppm as built; a magnet only 70-100 mm long around a 100 mm
bore could not flatten a 40 mm region. Length of ~2.5x the magnet radius (~140-150 mm here) was needed.

**Layer spacing.** Cubes only rotate about z, so the clearance between layers is set by the in-plane diagonal,
sqrt(2)*a = 8.98 mm, not the 3D diagonal sqrt(3)*a = 11 mm used in the original scripts. With 2 mm of resin that gives
11 mm between layer radii instead of 14 mm, which allows 3-4 layers within the Form 4 envelope.

**Optimization.** The genetic algorithm (pymoo NSGA2, as in the MRI4ALL scripts) finds the layout type but results vary
strongly between random seeds, and it needs a second objective (maximize B0) to keep population diversity: runs with the
field pinned to a narrow band collapsed onto poor designs. The large gains came from local refinement of the best GA
design: minimax optimization of ring z positions, rotations and small radial shifts (SLSQP), plus moving single magnets
between rings (`refine_minimax.py`). In the one direct comparison (single runs, which also settled at different fields), optimizing
the as-built metric gave a design at least as good on both metrics as optimizing the paper metric, which overfit the sparse
check points; more seeds would be needed to call this general.

**Magnets and materials.**
- N42 1/4" cubes work; N52 gives the same ppm with ~12% more field but corner damage starts near 42 C (lower iHc).
  N52SH/N48SH would fix that. Thin plates (10 x 5 x 2 mm magnetized through 2 mm) demagnetize themselves.
- Forcing two repelling cubes together can exceed N42's coercivity at room temperature: use jigs during assembly.
- Br drops ~0.12 %/C (~1,200 ppm per degree): temperature stability matters more than the static homogeneity.
- Resin (Formlabs Tough 2000 or Precision Model) places magnets to ~0.05 mm; worth ~5-10% as built versus 0.1 mm.
  Slip-fit cube pockets with corner reliefs and lids between rings hold magnets without glue.

## Tools added

- `magsimulator.py` (additions only, original functions unchanged):
  `build_magnet_collection`, `compute_demag_fields`, `required_ihc`, `suggest_grade`, `define_sensor_points_on_cylinder`,
  `shim_residual`, `tolerance_homogeneity` (as-built homogeneity from linearized per-magnet error sensitivities,
  checked against a full Monte Carlo within ~10%).
- `check_demag.py`: reverse H field on every magnet of a ledger and the minimum iHc / grade needed.
  `python check_demag.py design.xlsx --Br 1320 --max-temp 40`
- `optimize_robust_aligned.py`: the MRI4ALL optimization script generalized. Arguments:
  `<region diameter> <generations> [max magnets] [min B0] [second objective: magnets|b0|none] [max length] [bore]
  [cylinder length, 0 = sphere] [max layers] [objective: asbuilt|paper] [seed] [max B0]`;
  environment: `POOL`, `POP`, `NSENS`, `TAG`, `ALIGN` (1 = shared ring positions), `LAYER_WALL` (mm, switches to
  sqrt(2)*a spacing). Saves the Pareto front every 100 generations to `results/`.
- `refine_design.py`: Powell refinement of ring z positions and rotations plus +-1 magnet moves.
- `refine_minimax.py`: minimax (SLSQP) refinement with optional radial shifts, per-ring Jacobian, magnet budget,
  and screened magnet moves between rings. Produced the recommended design:
  `python refine_minimax.py results/REFINED_from_FINAL_bore100_cyl40x40_649mag_3layer.xlsx out.xlsx --radial --rounds 4`
- `optimize_bore50mm_*.py`: earlier iterations for a 50 mm magnet bore, kept for reference.

## Holder CAD

Step-by-step instructions (coupons, generating, printing, assembly) are in [CAD.md](CAD.md); this section records the
design choices.


`make_ring_cad.py` turns a ledger into printable parts (STEP/STL) plus a top-view loading guide (PNG) per ring:

- **SLA rings**, one per ring z: magnet pockets open at the top face over a 2 mm floor (8.65 mm thick). Pockets are the
  cube plus a clearance per local axis (`--tol x y z`), with corner reliefs, an air hole through the floor and a small V
  notch in the wall at the N (magnetization) face. Where two rings are closer than 3 mm + ring thickness (slabs 0/1 and
  7/8, 9.35 mm apart) the upper ring reaches down and sits directly on the lower one. Slab 0 has a 3 mm floor; the top
  slab is closed by the end cap (not generated).
- **FDM spacers** in the other gaps: solid rings (bore to OD, with the bolt holes) whose bottom face caps the pockets of
  the ring below; the slicer adds walls and infill, and they print flat without supports. For the
  recommended design: three heights (6.06, 10.21, 17.22 mm), one file each, each printed twice (the +z and -z spacers
  are identical).
  `--solid` instead makes every SLA slab fill the gap below it (no spacers).
- **Bolts**: 8 M3 outside the magnets (one shifted 10 degrees so the stack only goes together one way) and 6 M2 between
  the bore and the inner magnets, placed where they clear the pockets of every slab; M3 does not fit inside (only
  ~2.3 mm of resin between the 100 mm bore and the inner pockets). SLA and FDM hole clearances are set separately
  (`--outer-tol`/`--inner-tol`, `--fdm-outer-tol`/`--fdm-inner-tol`). A V notch on the outside of every part marks +x (B0).
- `assembly.step` imports as `magnet_holder` > `sla_rings` (slabs) + `fdm_spacers` (one part per spacer height, placed
  at each position); `assembly_with_magnets.step` adds `magnets`, one cube part placed 649 times (with `--magnets`).

```bash
python make_test_coupon.py --out cad/coupon        # SLA: pocket clearances x rotations, thin walls, bolt holes
python make_fdm_coupon.py --out cad/fdm_coupon     # FDM: spacer heights, bolt holes
python make_ring_cad.py results/BEST_bore100_cyl40x40_649mag_427ppm.xlsx --out cad/best --tol 0.15 0.15 0.3 \
    --outer-tol 0.4 --inner-tol 0.4 --fdm-outer-tol 0.5 --fdm-inner-tol 0.5 --spacer-comp 0 --magnets
python check_bolt_loads.py results/BEST_bore100_cyl40x40_649mag_427ppm.xlsx
```

For the recommended design: OD 193.6 mm, stack 147.2 mm without end caps, 1.5 L of SLA resin (2.9 L with `--solid`)
and FDM spacers enclosing 1.4 L (plastic used depends on infill). The thinnest resin between neighbouring pockets is 1.02 mm (0.15 mm clearance, 0.4 mm corner
reliefs); these are pinch points where two pocket corners nearly meet, under 1.5 mm for only 0.1-0.3 mm along a face. The
SLA coupon's 0.6/0.8/1.0 mm pairs test them. The +z and -z SLA rings are not identical (the slab boundaries are not mirror
symmetric), so all 9 are printed from their own files. Print every ring with the same resin, wash and cure, and turn
every other ring 90 degrees about its own axis on the build plate so print ovality alternates through the stack (the
parts themselves always assemble in the keyed orientation).

**Spacer heights set the ring spacing**, which matters at the 0.05-0.1 mm level (ring z errors of 0.05-0.1 mm sigma
cost ~1,200-1,300 ppm in the warp study above). FDM loads are small (under 10 N per bolt), but typical FDM height accuracy
is ~0.1 mm: print the FDM coupon with the spacers' exact settings, measure the height samples and pass the mean error as
`--spacer-comp` (printed 0.1 mm tall -> `--spacer-comp -0.1`), then check every printed spacer with calipers. PETG or
PLA; the parts see room temperature and light clamping only.

### Forces, bolts and resin

`check_bolt_loads.py` computes the force on every magnet from the surface-charge model (exact for rigid magnets) and the
load each joint between slabs puts on the bolts (positions as placed by `make_ring_cad.py`):

- Per magnet: median 1.25 N, max 2.9 N.
- Separating force across a joint: 20-73 N (largest at the joints between slabs 1/2 and 7/8), 39.5 N on each end cap;
  shear up to 9.8 N between the two outermost slabs; torques about z under 0.06 N m.
- Shared by 8 M3 + 6 M2 bolts: at most 7.8 N per M3 and 3.1 N per M2 (with the current inner hole placement). Brass
  (CW614N/C360, taking a conservative 150 MPa yield) holds ~750 N (M3) and ~310 N (M2), a margin of ~100x. Brass is non-magnetic; avoid 18-8 stainless, whose
  cold-worked threads can be slightly magnetic. The stack is ~150 mm plus caps, so M3 brass threaded rod with nuts is
  the practical form. A light preload (snug, ~0.1-0.2 N m on M3 with washers) is far more than the 5-9 N per bolt
  needed to keep joints closed and avoids crushing or creeping the resin; friction from it also carries the shear.
- Rigid 4000 (tensile 69 MPa, flexural 105 MPa, glass-filled for low creep): pressing the worst-loaded magnet (2.9 N)
  on a 1.0 mm wall treated as a full-height cantilever gives ~10 MPa, 7x below tensile strength, and the real thin spots
  are corner pinch points backed by thicker resin, so actual stresses are much lower. Average contact pressure on a
  pocket face is ~0.07 MPa, far below the 0.45 MPa at which its 77 C heat-deflection temperature is rated.

## Reproducing the recommended design

```bash
source .venv/bin/activate
# 1. GA: 40x40 cylinder, 100 mm bore, >=45 mT, <=650 magnets, <=150 mm, 3 layers, as-built objective
POOL=11 python optimize_robust_aligned.py 40 2000 650 45 b0 150 100 40 3
# 2. Powell refinement of the best Pareto design
python refine_design.py results/FINAL_bore100_cyl40x40_649mag_3layer.xlsx results/REFINED_from_FINAL_bore100_cyl40x40_649mag_3layer.xlsx --sensors 300
# 3. minimax refinement with radial shifts and magnet moves
python refine_minimax.py results/REFINED_from_FINAL_bore100_cyl40x40_649mag_3layer.xlsx results/BEST.xlsx --radial --rounds 4
```

The GA step is stochastic; a different seed gives a different starting design.

## Results directory

- `BEST_*.xlsx`: recommended design. `FINAL_*.xlsx`: GA result before refinement (2,057 ppm).
  `REFINED_*`, `MM2_*`, `PAPEROPT_*`: intermediate and comparison designs.
- `*_summary.csv` / `*_pareto<N>.xlsx`: Pareto fronts saved by the optimizer (sim_ppm = paper, asbuilt_p90_ppm = as built).
- `*.log`: optimizer and refinement logs. Subfolders hold snapshots saved before runs were restarted.

## Known issues in the original code

- `extract_3Dfields` / `plot_3D_field`: the field grid is stored as B[y, x, z], so the panel titled "xz plane" shows the
  yz plane, "yz plane" shows the xz plane, and vertical axes are flipped. Values are correct.
- Needs magpylib 4 (mm, mT). magpylib 5 uses SI units and silently gives wrong fields. NumPy >= 2.4 breaks the scripts'
  `int()` of 1-element arrays; `requirements.txt` pins both.
- `generate_ring_of_magnets` calls `sys.exit` if a ring is too crowded, which kills a whole optimization run; keep the
  per-ring magnet bounds within `get_max_magnets_per_radius`.
- `addcopyfighandler` forces a Qt backend (PyQt6 is in requirements), and the scripts end in `plt.show()`, which blocks.

## Next steps

1. Measure every magnet (Br and magnetization angle) and add a mode to `refine_minimax.py` that assigns measured magnets to
   positions and rotations. This is the main route below ~3,000 ppm as built.
2. Print a test ring and a pocket-clearance coupon; confirm the 184 mm ring fits upright in PreForm with supports.
3. Plan passive or electrical shimming for the remaining error.

## Field-driven design (`field_design.py`)

A CoilGen-style alternative to the GA: target field (Bx = B0 on the 40 x 40 mm cylinder surface) -> sensitivity
matrix of every candidate ring (layer radius x z on a 1 mm grid) -> LP for the ideal continuous magnet density ->
MILP for whole ring slices with integer magnet counts (the initial design) -> gradient refinement of slice z, ring
radius and rotation (Adam, then SLSQP minimax, analytic dipole gradients with an exact-cuboid correction) plus screened
+-1 magnet moves. `--ga` refines the same initial design with a pymoo GA instead, for comparison. Same constraints as
`refine_minimax.py` (<= 650 magnets, B0 >= 45 mT, <= 150 mm, 100 mm bore, OD <= 200 mm, 9.35 mm slices, 3 layers).

```bash
python field_design.py results/field_design/FD_gd.xlsx                 # ~3 min (MILP 2 min, gradients <1 min)
python field_design.py out.xlsx --start <initial.xlsx> --ga --pool 11   # GA refinement for comparison
```

| Design | Magnets | B0 (mT) | ppm, ideal (20k pts) | As built p90 (ppm) | Run time |
|---|---|---|---|---|---|
| Previous BEST (GA + Powell + minimax) | 649 | 45.00 | 430 | 4,423 | hours |
| GA from scratch (`optimize_robust_aligned.py`, paper objective) | 634 | 47.58 | 3,669 | 6,886 | 10 min |
| LP + MILP initial placement | 648 | 45.10 | 64 | 4,716 | 2 min |
| MILP + GA refinement | 648 | 45.02 | 189 | 4,619 | 19 min |
| **MILP + gradient refinement** (`results/field_design/FD_gd.xlsx`) | 650 | 45.00 | **45** | 4,644 | **3 min** |

As built uses 1% Br, 0.05 mm, 1 degree, 400 builds, no shimming. The field-driven design is ~10x better on paper but
not better as built: below ~1,000 ppm magnet errors dominate (see above). Its smallest rings (13 and 16 magnets) are
valid only because the check points are dense enough to see their azimuthal ripple; always validate on >= 2,000 points.
