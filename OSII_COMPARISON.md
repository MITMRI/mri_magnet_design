# Comparison with the OSI² 30 cm "ROMA" magnet, and magnet-interaction study

This records a comparison between our recommended design (`results/BEST_bore100_cyl40x40_649mag_427ppm.xlsx`) and
the OSI² ONE magnet ([gitlab.com/osii/magnet/30cm-halbach-magnet](https://gitlab.com/osii/magnet/30cm-halbach-magnet),
`main` branch = v2.1), what explains the difference, and a follow-up study of magnet-magnet interaction (finite
permeability), which turned out to matter for our design too. All simulations are magpylib 4 (same models as
`magsimulator.py`) unless stated otherwise. The analysis scripts were throwaway and are not kept; the methods are
described below in enough detail to redo them (core code in the appendix).

## Summary

- **On paper (ideal magnets), OSII is ~10x more homogeneous than ours at the same relative size**: scaled to a 100 mm
  bore it gives 35 ppm on a 40 mm sphere vs our 306 ppm; 65 vs 417 ppm on a cylinder of diameter = length = 0.4 x bore.
- **As built, the gap mostly disappears**: with our usual error model both land at ~2,400-2,900 ppm (p90). The
  OSII magnet as actually built and measured (INRiM replica) was **24,600 ppm** over its 200 mm DSV against a ~400 ppm
  design.
- **Why theirs is better on paper**: their positions with pure Halbach angles would give 75 mT; the "rotation optimized"
  (ROMA) angles bring it down to 50 mT. They spend a third of the field, spread over 2,320 free angles, on homogeneity.
  Ours is pure Halbach at its maximum field (45 mT with 649 magnets), so it has no field to trade. Rotation-optimizing our
  layout does improve paper ppm if B0 is lowered (438 -> 174 ppm at 40 mT) but makes the as-built magnet slightly *worse*.
- **Their algorithm is not published.** The array was optimized at LUMC (Leiden) for homogeneity in a 200 mm DSV; earlier
  Leiden work used a genetic algorithm choosing {no magnet, Halbach magnet, flipped magnet} per slot. The v2.1 angles are
  continuous, so the method went beyond that.
- **Magnet-magnet interaction (finite recoil permeability) costs our design 110-1,900 ppm**, depending on the
  perpendicular susceptibility chi_perp (0.05-0.3), which we do not know. About 93% of it is a smooth second-order
  pattern, so it can be designed out (ring z / radius adjustment) or shimmed, but only if chi_perp is measured.
- Other findings: by our demag criterion the OSII N52 corners are past the knee at 20 C; their per-half axial force is
  ~650 N vs our ~20 N; magnetic torque rotating magnets in their pockets was the largest single cause of their as-built
  error, and our pocket clearance (0.15 mm on 6.35 mm, ~1.4 deg of play) is looser in angle than theirs.

## 1. The OSII design

From the repo (`README.md`, `src/BOM.csv`, `doc/*.md`, `CHANGELOG.md`, orientation files in
`src/magnet orientations/raw/`; the `.par`/`.FCStd` CAD files are Git LFS and were not downloaded):

| | OSII ONE v2.1 | Ours (BEST) |
|---|---|---|
| Magnets | 1936 x 12 mm cubes + 384 x 12x12x50 mm bars, N52 (Br 1.45 T in their file) | 649 x 6.35 mm cubes, N42 (Br 1.32 T) |
| Magnet volume / mass | 6,110 cm3, ~46 kg (whole magnet ~100 kg) | 166 cm3, 1.25 kg |
| B0 | 50.0 mT (our model; they quote ~50) | 45.0 mT |
| Bore | 322 mm | 100 mm |
| Layout | 8 ring pairs at a fixed 20 mm pitch (z = +-10 ... +-150), each 57 @ r 171 + 64 @ r 192; end rings at +-187 with three layers of 50 mm bars (57/64/71 @ r 171/192/213) | 9 slabs, irregular z (9.35-26 mm apart), per-ring radii, 6-39 magnets per layer, 3 layers |
| Optimized | each magnet's angle ("ROMA"), positions fixed | positions, counts, ring z; angles pure Halbach |
| Holder | milled 15/20 mm PP rings, PC lids screwed on, 12 M6 + 8 M4 brass rods, PA spacer sleeves | SLA rings + FDM spacers, M3/M2 brass rods |
| Shimming | 12 printed shim trays with 6 mm cubes | none yet |

The orientation file `analytic-18 rings.txt` has columns x, y, z, dx, dy, dz, Br, angle (mm, T, deg) for the z > 0
half (1,160 magnets); the per-ring files `ring 0..8.txt` carry the same positions and angles. The -z half is the mirror
image with the same in-plane magnetization. With magnetization at `angle` in the xy plane, our model gives
B0 = +50.0 mT along +x (their docs say "north in -x"; only the sign convention differs). Their uniform ring geometry
means one ring template and one bolt pattern for all rings.

**How far ROMA is from Halbach** (deviation of each angle from 2*theta):

- 50% within 15 deg, 23% at 15-45 deg, 20% at 45-135 deg, 7% reversed (> 135 deg); volume-weighted mean cos(dev) = 0.77.
- The deviations concentrate in a few rings: inner layer at z = 70 (median 68 deg, 15/57 reversed), z = 130 (both
  layers ~64-79 deg), z = 150 (median 106-130 deg, 25/57 and 22/64 reversed). The end rings (z = 187) and z = 90 are
  nearly pure Halbach.
- Same positions with pure Halbach angles: **75.2 mT**, but 5,243 / 12,673 / 28,140 ppm on DSV 100 / 150 / 200, against
  ROMA's 20 / 56 / 420 ppm.

## 2. Homogeneity on paper

Peak-to-peak / mean of Bx on the region surface (800 sensor points; field extremes of a source-free region lie on its
boundary).

OSII, spheres: DSV 40: 4, 100: 20, 140: 44, 150: 56, **200: 420**, 250: 20,005 ppm. (Their design target is a 200 mm DSV;
the INRiM paper's ideal-dipole value is 266 ppm on its own measurement grid.)

Same relative region size (D = region diameter, cylinders have length = D):

| D / holder bore | OSII cylinder | OSII sphere | Ours cylinder | Ours sphere |
|---|---|---|---|---|
| 0.3 | 19 | 18 | 143 | 86 |
| 0.4 | **65** (129 mm) | 35 | **417** (40 mm) | 306 |
| 0.5 | 660 | 74 | 3,837 | 1,052 |
| 0.6 | 5,639 | 291 | 16,347 | 3,508 |

Normalizing by the inner magnet-center radius instead (D = 0.70 r_in): OSII 37 ppm (120 mm), ours 443 ppm.

**Scaled to our size.** Scaling a magnet uniformly leaves B0 and ppm unchanged, so this is the same comparison done
explicitly: OSII scaled x0.311 (bore 322 -> 100 mm) gives **35 ppm on a 40 mm sphere**, 50.0 mT (45.5 mT with N42
Br 1.32); scaled x0.336 (inner magnet radius 171 -> 57.5 mm) gives 29 ppm. Ours: 306 ppm. The scaled OSII would need
~2,320 magnets of 3.7-4.0 mm (cubes) and 15.5 mm (bars).

## 3. Homogeneity as built (our error model)

Same model as `magsimulator.tolerance_homogeneity` (linearized per-magnet sensitivities, 1% Br sigma, 1 deg
magnetization-angle sigma, position sigma as listed, 400 random builds, linear shim terms removed), generalized to
non-cubic magnets. Median / 90th percentile ppm:

| Position sigma | OSII, 129 x 129 mm cylinder | OSII, DSV 200 | Ours, 40 x 40 mm cylinder |
|---|---|---|---|
| 0.05 mm | 1,864 / 2,421 | 4,717 / 6,038 | 2,173 / 2,905 |
| 0.10 mm | 1,884 / 2,434 | 4,824 / 6,156 | 2,396 / 3,100 |
| 0.20 mm | 1,971 / 2,540 | 5,144 / 6,643 | 3,082 / 4,053 |

Both magnets are dominated by magnet errors (Br and angle), as DESIGN_NOTES.md already found for ours. OSII's remaining
advantage comes from 3.6x more magnets (random errors average down) and 12 mm magnets, which make a given placement error
relatively small (insensitive to position sigma). The model does not include deterministic effects (next section), which
is why it predicts ~6,000 ppm on DSV 200 while the real build measured 24,600.

## 4. What happened when OSII was built

U. Zanovello et al., "Very-low-field MRI scanners: from the ideal to the real permanent magnet array", IEEE Trans.
Instrum. Meas. 2025, doi:10.1109/TIM.2025.3638912 ([arXiv 2509.11762](https://arxiv.org/pdf/2509.11762)). INRiM built a
replica of the v2.1 array and mapped it (4,224 points in the 200 mm DSV):

- Measured mean Bx 48.056 mT, **DIS1 (peak-to-peak) 24,600 +- 1,300 ppm**, DIS2 (integral metric) 4,400 ppm.
- Deterministic causes, added in sequence to the model:
  1. Finite permeability / magnet reaction (linear or nonlinear H-J curve, dipole model and FEM agree): ideal-dipole
     266 ppm -> **1,720 ppm**.
  2. Different measured characteristics of the 12 mm cubes and the 50 mm bars (both N52, same supplier): -> 6,960 ppm.
  3. Measured ring z positions (attraction/repulsion between rings shifts them despite the threaded rods): -> 6,520 ppm.
  4. Magnet rotation in the pockets driven by magnetic torque, +-0.8 / 1.0 / 1.2 deg with the sign of the computed
     torque (pockets 0.13-0.2 mm larger than the 12 mm magnet): -> 14,670 / 18,930 / **23,200 ppm**, matching the
     measurement.
- Stochastic effects (magnet characterization, orientation, positioning) add only a few percent on top.
- Their recommendations: avoid mixing magnet types, reduce pocket play, or mill pockets pre-rotated against the torque.

The OSI² certification review of v2.1.3 (approved, 14/18) contains no field measurements; it notes missing test
routines and a missing risk assessment. A related report says the built magnet was shimmed with passive trays
optimized with the Shimmer toolbox.

**The optimization algorithm** is not in the repo, the review or the paper. The paper says only that the orientations
"are the result of an optimization procedure designed to maximize the magnetic field homogeneity within a 200 mm DSV",
done at Leiden University Medical Center (frame engineered at PTB). Earlier Leiden work (O'Reilly, Teeuwisse, Webb,
[arXiv 2005.02834](https://arxiv.org/pdf/2005.02834)) optimized ring diameters and shimmed with a genetic algorithm
choosing, per grid slot, no magnet, a Halbach-oriented magnet, or the same magnet flipped 180 deg; the 7% reversed
magnets in v2.1 fit that lineage, but 43% sit at intermediate angles, so v2.1 used continuous angles. ISMRM 2022 #0063
(same authors) reports 280 ppm over 20 cm in simulation with higher-order Halbach shim inserts.

## 5. Why theirs is better on paper: experiments

**Field headroom.** Pure Halbach angles maximize B0 for given positions. OSII positions: 75 mT Halbach (68 mT at
Br 1.32), 50 mT ROMA, so ~1/3 of the field is spent on homogeneity. Ours: 45.0 mT is already the Halbach maximum for our
649 positions, and the requirement is B0 >= 45 mT.

**Best possible with free angles (convex bound).** The field is linear in each magnet's magnetization vector, so with
m_i = (m_x, m_y), |m_i| <= 1 (relaxation of unit magnets; 16-gon inscribed constraint), minimizing max |B - B0| on the
40 x 40 mm cylinder (500 points) at a fixed target B0 is a linear program (scipy HiGHS). Results (lower bounds; most
|m_i| came out at 1, so they are nearly attainable):

| Layout | B0 target | Bound on 40 x 40 mm cylinder |
|---|---|---|
| Ours | 45 mT | infeasible (needs every magnet at full Halbach; the 0.981 polygon radius gives 5,048) |
| Ours | 40 mT | 163 ppm |
| Ours | 35 mT | 77 ppm |
| Ours | 30 mT | 15 ppm |
| OSII scaled to 100 mm bore, Br 1.32 | 45 mT | 10 ppm |
| OSII scaled | 40 mT | 7 ppm |

So with free angles, our positions can only beat ~400 ppm by giving up field, and OSII's dense layout has far more
freedom at the same field.

**ROMA-style optimization of our positions.** Sequential LP on the angles (exact unit magnets): linearize B in the angle
changes, minimize max |B - B0| with a trust region on each step, continuation from Halbach while lowering the target B0.
Checked on 1,500 independent points; as built as in section 3 (0.05 mm).

| Target B0 | Paper ppm | Median deviation from Halbach | As built median / p90 |
|---|---|---|---|
| 45.0 (current Halbach) | 438 | 0 | 2,173 / 2,905 |
| 45 | 389 | 0.7 deg | 2,212 / 2,944 |
| 44 | 407 | 7.6 deg | |
| 43 | 309 | 10 deg | |
| 42 | 246 | 8 deg | 2,302 / 3,081 |
| 41 | 199 | 7 deg | |
| 40 | 174 | 7 deg | 2,456 / 3,319 |
| 38 | 143 | 11 deg | |
| 36 | 109 | 14 deg | |

(438 is our design in this model, which keeps cube bodies unrotated; the exact value is 425-427.) Paper ppm improves
steadily, as built gets slightly worse: magnets turned off Halbach partly cancel, so B0 drops while every magnet's error
stays the same size. **Rotation optimization does not help our as-built magnet.** It would be useful in combination with
measured magnets (assigning each measured magnet a position and rotation), where the angle is the knob that cancels
measured errors.

Failed attempts, for the record: L-BFGS on a smooth max from the Halbach start stalled at 390 ppm regardless of the B0
floor (local minimum); normalizing the LP solution to unit magnets and polishing gave 6,000-40,000 ppm (rounding destroys
the LP's cancellations). The sequential LP with continuation worked; below 36 mT the dense LP became numerically fragile
and was not pursued.

## 6. Other comparisons

**Demagnetization** (method of `check_demag.py`: reverse H from all other magnets plus self, at the center and on a
3^3 grid over 90% of each magnet):

| | Center, median / max | Grid max |
|---|---|---|
| OSII 12 mm cubes | 338 / 612 kA/m | 823 kA/m |
| OSII 50 mm bars | 503 / 677 kA/m | 866 kA/m |
| Ours | 344 / 425 kA/m | 621 kA/m |

Standard N52 has iHc ~876 kA/m (knee ~0.9 iHc = 790 kA/m), so by our criterion OSII's corners are past the knee at
20 C (localized corner damage, small effect on the field). Ours: 621 against 955 kA/m for N42, fine to ~66 C.

**Forces** (surface-charge model of `check_bolt_loads.py`, generalized to bars): OSII per magnet median 6.1 N, max
40.7 N (cubes), median 22.7 N, max 119 N (bars). The +z half is pushed away from the -z half with +649 N, while the end
ring is pulled inward with -626 N; the net axial force on everything beyond each ring varies between -626 and +649 N.
Ours: +19.6 N on the +z half. (Hence their M6 rods and thick PP rings.)

**Worth borrowing from OSII:**
1. Their QM process: check every magnet, pocket rotation < 1 deg, H9 pocket fit, a caliper measurement plan for pocket
   positions per ring. A similar per-slab checklist from our CAD would fit CAD.md.
2. Shim trays (12 printed trays with 6 mm cubes); the equivalent for us would use ~3 mm cubes.
3. Pocket play and torque: the dominant cause of their as-built error. Our 0.15 mm clearance on 6.35 mm cubes allows
   ~1.4 deg, more than their 1.2 deg. Tighten it, or pre-rotate each pocket against its computed torque.
4. Avoid mixing magnet types (they mix cubes and bars; we already use one type).

## 7. Magnet-magnet interaction in our design

Real NdFeB magnets have recoil permeability slightly above 1, so each magnet's magnetization responds to the field of
the others. magpylib treats magnets as rigid, so all our ppm numbers so far ignore this. The INRiM paper found it moves
OSII from 266 to 1,720 ppm before any manufacturing error.

**Model.** Each magnet is linear: M_i = Mr + chi . <H_others>, with chi_par along the easy axis and chi_perp across it,
H averaged over the magnet volume at 2x2x2 Gauss points, iterated to convergence (a few iterations, chi is small).
Self-demagnetization is omitted: for identical magnets it scales every magnet the same way in its own frame.
chi_par ~0.05 is well established (recoil permeability ~1.05); chi_perp is larger and poorly known, so it was swept.

**Validation on OSII:** rigid 430 ppm (DSV 200, 1,500 points) -> 1,500 ppm with chi = 0.05 both ways, -> 7,523 ppm with
chi_perp = 0.15. The paper's 266 -> 1,720 ppm is the same order, so the model is reasonable.

**Ours** (40 x 40 mm cylinder, 1,500 points; B0 changes by < 0.06 mT):

| chi_par / chi_perp | Total ppm | Change from interaction alone | Change after removing linear terms | Change after removing 2nd-order terms |
|---|---|---|---|---|
| rigid | 425 | - | - | - |
| 0.05 / 0.05 | 535 | 138 | 123 | 36 |
| 0.05 / 0.15 | 1,102 | 903 | 865 | 65 |
| 0.05 / 0.30 | 2,316 | 2,149 | 2,014 | 157 |

- The effect is significant: at chi_perp = 0.3 it adds as much as all our random magnet errors (~2,900 ppm as built).
- About 93% of it is a smooth second-order pattern. Gradient offsets (linear) do not remove it, but small ring z and
  radius changes can, and passive shimming targets exactly these terms.
- **Compensation plan:** treat the interaction correction dB as fixed for small design changes. Compute dB for the
  current design, refine positions/rotations (`refine_minimax.py`) so that B_rigid + dB is uniform, recompute dB, repeat
  (2-3 rounds expected). Verify the final design with the full iteration and with FEM.
- **Blocker:** the correction is only right if chi_perp is right (0.05 vs 0.3 is 900 vs 2,000 ppm of correction).
  Options: measure a few cubes (VSM or hysteresisgraph, easy axis and transverse); build a ring or two, map the field and
  fit chi_perp; or optimize the worst case over chi_perp in 0.05-0.3.

**FEM.** FEM with the recoil permeability (anisotropic, ideally the nonlinear curve) captures the interaction
automatically; the INRiM FEM agreed with their dipole-with-reaction model. Our PDE Toolbox run (`verify/`) has ~2% RMS
mesh error, 20-200x too coarse for 100-1,000 ppm effects. With more compute, either refine the mesh a lot, or run
mu_r = 1 and mu_r = 1.05 on the identical mesh and use the difference (mesh error largely cancels). Integral methods
(Radia, MagTense: only the magnets are discretized, no air mesh) are better suited to ppm-level fields; subdividing each
cube would refine the 8-point average used here.

## 8. Recommended next steps

1. Measure chi_perp (and chi_par, Br spread) on a sample of our cubes.
2. Add the interaction-correction loop to `refine_minimax.py`; until a measurement exists, use chi_perp = 0.15.
3. Compute the torque on every magnet and the resulting pocket rotation; tighten pocket clearance or pre-rotate pockets.
4. Cross-check the interaction effect with FEM (frozen-mesh difference) or Radia/MagTense.
5. Measured-magnet assignment (DESIGN_NOTES.md next step 1) and shim trays remain the route below ~3,000 ppm as built.

## Appendix: core code

Interaction iteration (magpylib 4; mm, mT, kA/m). `cs` is a list of `magpy.magnet.Cuboid` with magnetization along
local x:

```python
MU0 = 4e-7*np.pi*1e6   # mT per kA/m (1.2566)

def interact(cs, Br, chi_par, chi_perp, iters=6):
    g = np.array([[sx, sy, sz] for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)])/np.sqrt(3)/2
    pts = np.array([c.orientation.apply(g*np.array(c.dimension)) + c.position for c in cs])   # (n, 8, 3)
    n = len(cs)
    for c in cs: c.magnetization = (Br, 0, 0)
    for _ in range(iters):
        H = magpy.getH(cs, pts.reshape(-1, 3), sumup=False).reshape(n, n, 8, 3)
        Hav = (H.sum(0) - H[np.arange(n), np.arange(n)]).mean(1)          # field of the others, volume average
        for c, h in zip(cs, Hav):
            hl = c.orientation.inv().apply(h)
            c.magnetization = (Br + chi_par*MU0*hl[0], chi_perp*MU0*hl[1], chi_perp*MU0*hl[2])
    return cs
# field: magpy.getB(cs, points, sumup=True)   (with a list of sources, sumup defaults to False)
```

Free-angle bound: with `Bx[i, k]`, `By[i, k]` the x-field at sensor k from magnet i magnetized along global x / y, the
field is `mx @ Bx + my @ By`; minimize t subject to |field - B0| <= t and cos(a) mx_i + sin(a) my_i <= cos(pi/16) for
16 angles a (scipy `linprog`, HiGHS). Sequential LP: at angles phi, J = -sin(phi) Bx + cos(phi) By; minimize t subject
to |b + J^T d - B0| <= t, |d_i| <= step; accept and grow the step if max |b - B0| drops, else shrink it.

OSII geometry: load `src/magnet orientations/raw/analytic-18 rings.txt` (x, y, z, dx, dy, dz, Br, angle), mirror z
for the other half with the same angle, `Cuboid(magnetization=(Br*1000, 0, 0), dimension=(dx, dy, dz))` rotated by
`angle` about z.
