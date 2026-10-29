# MRI4ALL Open-Source Tools for Design of the Hallbach-Array Magnet

Please watch the <a href="https://www.youtube.com/embed/iKs5pwwoyoQ" target="_blank">video below</a> for an introduction to the tools provided in this repository.

[![Overview of the MRI4ALL Magnet Tools](https://img.youtube.com/vi/iKs5pwwoyoQ/0.jpg)](https://www.youtube.com/watch?v=iKs5pwwoyoQ)

## Setup

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The scripts use magpylib 4 units (mm, mT, kA/m); magpylib 5 switched to SI and will silently give wrong fields, so it is pinned in `requirements.txt`.

## Demagnetization check (required iHc)

`check_demag.py` computes the reverse H field every magnet in a ledger sees (from all other magnets plus its own self-demagnetization) and the minimum intrinsic coercivity iHc the magnets need:

```bash
python check_demag.py my_design.xlsx --max-temp 40 --safety 1.2 --plot
```

- `--max-temp`: hottest the magnets will ever get, including assembly and epoxy curing. iHc is derated at -0.6 %/°C (CY-Mag datasheet in `magnet_properties/`).
- `--Br`: remanence in mT (default 1270, as in the optimization scripts).
- `--out results.xlsx`: per-magnet reverse fields.

The required iHc only covers the finished array. Handling during assembly can be worse: two ½" cubes forced together face-to-face in repulsion see up to ~930 kA/m near their edges.

## Magnet holder CAD

See [CAD.md](CAD.md) for generating the printable ring holders (SLA rings + FDM spacers) from a design, the SLA and FDM
tolerance coupons to print first, and assembly:

```bash
python make_ring_cad.py results/BEST_bore100_cyl40x40_649mag_427ppm.xlsx --out cad/best --magnets
```

## Field-driven design (current method)

See [FIELD_DESIGN.md](FIELD_DESIGN.md): `field_design.py` designs the ring layout from a target field, CoilGen style
(sensitivity matrix -> LP ideal magnet distribution -> MILP ring placement -> gradient refinement), in about 3 minutes:

```bash
python field_design.py results/field_design/FD_gd.xlsx
python plot_magnets_3d.py views.png results/field_design/FD_gd.xlsx:"Field-driven"
```

## Compact magnet design study (this fork)

See [DESIGN_NOTES.md](DESIGN_NOTES.md) for the 100 mm bore / 40 x 40 mm region design built from 1/4" N42 cubes, the tools
added for it (`optimize_robust_aligned.py`, `refine_minimax.py`, as-built tolerance analysis), the ring-holder CAD generator (`make_ring_cad.py`,
`make_test_coupon.py`, `make_fdm_coupon.py`, `make_coupons_step.py`, `check_bolt_loads.py`), and what was learned.

See [OSII_COMPARISON.md](OSII_COMPARISON.md) for a comparison with the OSI² 30 cm "ROMA" magnet and a study of
magnet-magnet interaction (finite permeability) in our design.
