# 12 cm imaging-sphere Halbach magnet: requirements and design method

Script: [`field_design_12cm.py`](../field_design_12cm.py) (self-contained, ~390 lines). Results: `results/12cm/`.

## Requirements

| | |
|---|---|
| Imaging volume | 120 mm diameter sphere (DSV) |
| Homogeneity | < 50 ppm peak-to-peak of Bx over the sphere (ideal magnets, no shimming) |
| Bore (clear ID, plastic) | 200 mm (170 mm acceptable minimum); must fit RF coil + 3 gradient sets |
| Outer diameter | <= 325 mm including plastic |
| Magnet length (end to end) | <= 360 mm |
| Goal | Maximum B0 |

## Magnets

KJ Magnetics BX088SH: N42SH NdFeB block, 25.4 x 12.7 x 12.7 mm (1" x 1/2" x 1/2"), magnetized through a
12.7 mm face, Ni-Cu-Ni, 30.7 g, Br max 1.32 T, max operating temperature 150 C. Modeled with Br = 1.30 T; B0 scales
linearly with Br (every 1% of Br is 1% of B0, homogeneity unchanged).

Orientation: 25.4 mm length along the bore axis (z); the 12.7 x 12.7 mm cross-section in the ring plane;
magnetization in the ring plane at 2*theta for a block at azimuth theta (k = 1 Halbach), B0 along x.

## Layout rules (all command-line options)

| Rule | Value |
|---|---|
| Plastic wall at bore and at OD | 2 mm (`--bore-wall`, `--od-wall`) |
| Plastic between slices (axial) | 3 mm, so slices >= 28.4 mm apart (`--slice-wall`) |
| Plastic between ring layers (radial) | 2 mm, so rings in a slice >= 19.96 mm apart (block diagonal + 2 mm, `--layer-wall`) |
| Plastic between blocks in a ring | >= 1.5 mm at the corners (`--ring-wall`) |
| Ring structure | rings in mirror-symmetric slices (shared z planes at +-z, optional z = 0 slice); every block of a ring at the same z |
| Min blocks per ring | ripple rule: (60 mm / r)^(n-1) <= 2e-6 (`--ripple`), so no small near-bore rings whose lumpy field the linear model misses |

With a 200 mm bore and 325 mm OD, ring centers span r = 111.0-151.5 mm, exactly enough for three radial layers
(111.5 / 131.5 / 151.5 mm). With 3 mm walls only two layers fit and B0 drops to ~68 mT (coarse run).

## Method

1. Field table: exact Bx per magnet of a full ring (magpylib cuboids) for every candidate radius and slice z on a
   grid, on check points covering a quarter of the sphere surface (Bx is mirror symmetric in y and z; it is harmonic,
   so its extremes over the ball lie on the surface).
2. MILP (Gurobi, all cores): choose slices, rings and integer magnet counts that maximize B0 with
   |Bx - B0| <= eps*B0 on every check point. `eps` holds back a margin below the ppm cap (`--margin`; the dense check
   usually lands within ~1 ppm of the target, so 1-2 ppm is enough).
3. Gradient refinement (SLSQP, +-2 mm trust region per step): move slice z and ring radii to raise B0 under the same
   constraints. Point-dipole Jacobian (3 dipoles per block), exact field via a correction term.
4. Greedy +1 magnet per ring, re-refined, kept if B0 rises.
5. Final check on 20,000 points on the full sphere surface.

Run time ~2-3 min per design on 12 cores with Gurobi (free academic license, `~/gurobi.lic`). The field tables are
computed in parallel too.

B0 vs the ppm cap is a staircase rather than a slope: whole ring layers become fillable at particular ppm values, so a
few ppm near a step can change B0 by 50%. Sweep the cap when choosing specs.

## Designs

Each saved design lives in its own folder under [`designs/`](../designs/) with its block placements, specs, the
exact command and settings that produced it, and its gap check.

## Not yet done

- As-built tolerance (Br spread, magnetization angle, placement) and shimming: real builds land in the thousands of
  ppm before shimming, see SHIMMING.md.
- Demagnetization check (N42SH: iHc >= 1592 kA/m, large margin expected).
- Forces between rings / holder structure, CAD.
- Temperature: NdFeB Br drops ~0.12%/C, so B0 (and the Larmor frequency) drifts with temperature.
