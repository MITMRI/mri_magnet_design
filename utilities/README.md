# Utilities

Analysis scripts for a finished design. Each one reads a design's `placements.csv` (written by
`field_design_12cm.py`: one row per 12.7 x 12.7 x 25.4 mm N42SH block, columns `x, y, z, rot_z_deg`; magnetization
along the block's local x, rotated by `rot_z_deg` about the bore axis z). Run them from the repo root:

```
python utilities/<script>.py designs/<design>/placements.csv
```

All of the force scripts integrate the magnetic surface charge on each block's two pole faces (Gauss-Legendre
quadrature) in the exact cuboid fields from magpylib, assuming Br = 1.3 T, ideal blocks and no iron. The heavy ones are
parallelized and chunked to stay at about 4 GB of RAM; `--jobs` sets the number of worker processes.

| Script | What it answers | Writes | Runtime |
|---|---|---|---|
| `plot_design_3d.py` | 3D views of the blocks (iso, end, side, cutaway) | `magnets_*.png` | seconds |
| `check_gaps.py` | Smallest gap between blocks: within a ring, between rings of a slice, between slices | stdout | seconds |
| `tolerance_12cm.py` | Monte Carlo as-built homogeneity for position, magnetization-angle and Br errors | stdout | ~1 min |
| `forces_12cm.py` | Force and torque on every block in the finished magnet; net load per ring and per slice; force between the two halves | `forces.csv` + stdout | ~45 s |
| `assembly_12cm.py` | Force on each plate (and on each block) while stacking the plates in several orders, and while closing all plates together on long rods; rod length needed | `assembly_pairs.pkl` (cache) + stdout | ~4 min first run, seconds after |
| `insertion_12cm.py` | Axial force on seated blocks while a neighbour is pushed into its pocket, and on the block being inserted | stdout | ~1 min |
| `plot_slice_2d.py` | 2D cross-section of one plate (default: the one with the most blocks): block footprints, magnetization, bore, OD, DSV | `slice_z<z>.png` | seconds |
| `plot_shim_12cm.py` | Pictures of a shim solution: Bx over the DSV surface and central planes before/after (real cubes), the cubes ring by ring | `field_sphere_b<i>.png`, `field_planes_b<i>.png`, `shim_layout_b<i>.png` | ~1 min |
| `shim_12cm.py` | Passive shim insert (a stack of printed rings in the bore holding 1/8" cubes): placement from a field map or simulated as-built builds, as-built ppm with the cubes' own errors, loads on the cubes and the stack | `shim_<cube>mm/` | ~6 min for 3 builds |

Saved outputs live next to the design, e.g. `designs/12cm_bore200/forces.txt`, `assembly.txt`, `insertion.txt`
(stdout redirected), and `designs/12cm_bore200_84mT/tolerance*.txt`.

## Conventions

- `forces_12cm.py`: `Fr`, `Ft`, `Fz` are radial, tangential and axial components (N); `T`/`Tz` torque about the
  block centre (N m). Per-slice "outboard stack" is the axial load a holder must carry at that slice: the sum over it
  and every slice farther from the midplane.
- `assembly_12cm.py`: "push" > 0 means the incoming plate is pushed away from the stack (it must be pressed in). The
  signed per-block table is for building from the -z end; "back" = toward -z (toward the plates already built).
  The long-rod model treats plates as sliding without friction on the rods, with contact between touching plates;
  rod length = extra spread E + magnet span + 60 mm (10 mm end plates, 20 mm nut/washer/spare thread per end).
- `insertion_12cm.py`: pockets have a floor on the outward face (away from the midplane) and are open toward the
  midplane; plates are filled one at a time before stacking. "Pull-out" is force on a seated block toward its open
  face.
- `shim_12cm.py`: continuous optimum (least magnet for a target ppm, Gurobi barrier) -> cubes by walking each ring ->
  SLSQP refinement of every cube's ring angle and magnetization angle with exact fields. `--targets` sets the ppm
  targets tried; the cube count follows from the physics, at most `--max-cubes`. "real" = the cubes' own Br/angle
  errors (default +-3% / +-1.5 deg / +-0.1 mm), one installation; the main blocks get +-0.1 mm too (`--pos-tol`). `--parallel` builds run at once (~1.5-2 GB each).
  One shim ring per main plate (boundaries halfway between plates), cubes at each ring's +z face.
  `--map map.csv` (x, y, z mm, Bx mT): shim a measured map instead of simulated builds; a few hundred to 1,000 points
  on the DSV surface are enough (Bx is harmonic, so its extremes are on the surface).
- `tolerance_12cm.py`: errors are uniform within +-tol; the field is linearized in the errors and checked against
  exact magpylib on a few trials (`--check`).
