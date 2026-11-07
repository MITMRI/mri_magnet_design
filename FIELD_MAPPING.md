# Field mapping for shimming: harmonics, sensor density, accuracy

How densely the field has to be measured before shimming, and why a sparse map is enough. Simulated on `FD_gd`
(650 x 1/4" cubes) and `FDL_350` (349 x 12.7 mm long magnets), target 40 x 40 mm cylinder, as-built errors as in
SHIMMING.md (Br 0.8%, angle 0.7 deg, position +-0.25 mm). Both designs give nearly the same answers.

The analysis scripts were scratch work and are not in the repo; the numbers below are their output.

## Why a sparse map works (Laplace)

Inside the bore there are no magnets, so each field component obeys Laplace's equation. Two consequences:

- **The surface sets the inside.** Measuring a closed surface around the target fixes the field everywhere inside it;
  no 3D grid is needed.
- **Fine detail dies out.** Around the bore the field is a sum of sines, cos(m theta); a ripple of wavelength lambda
  shrinks like exp(-2 pi d / lambda) over a distance d from its source. The magnets are ~33 mm from the target surface,
  so short-wavelength detail barely reaches it.

So the field is fitted with a few smooth terms: solid spherical harmonics r^l Y_lm up to degree L, (L+1)^2 numbers.
The series converges inside the sphere that reaches the nearest magnet (r ~ 53 mm), which contains the 40 x 40 mm
cylinder (corners at 28.3 mm). A target 120 mm long does not fit: even the ideal field varies 58,000-75,000 ppm over
a 40 x 120 mm cylinder, and its ends approach the end rings.

## Harmonic content of the as-built field

`FD_gd`, ppm RMS on the R = 28.3 mm sphere (the sphere around the cylinder), average of 200 simulated builds:

| Degree l | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|---|---|---|
| as built | 740 | 618 | 461 | 333 | 218 | 140 | 53 | 19 | 6.5 |

- Each degree is ~0.6x the previous (0.55x for `FDL_350`), as expected for sources ~48-51 mm from the centre.
- Around the bore, on the r = 20 mm surface, azimuthal order m falls ~0.45x per step (798 ppm at m = 0, 3.8 at m = 8).
- Along z, wavelengths under ~30 mm hardly matter: keeping m <= 6 and wavelengths >= 30 mm leaves ~75 ppm p2p.
- Magnet errors dominate. The ideal designs' own terms are <= ~30 ppm, except a 127 ppm degree-12 ring ripple in
  `FDL_350` (its 12-magnet ring).
- The transverse components (By, Bz <= 0.2 mT) change |B| vs Bx by only 3-4 ppm, so mapping Bx is enough.

## Sensor layout and density

Fit error = fitted minus true field, peak-to-peak over the cylinder surface and interior, median of 120 builds, exact
sensors.

**Sphere of radius 28.3 mm** (an arc of sensors rotated about the bore axis; Gauss-Legendre latitudes):

| Grid (latitudes x rotations) | Points | Fit, `FD_gd` | Fit, `FDL_350` | Fit after subtracting the design field |
|---|---|---|---|---|
| 7 x 13 | 91 | 323 | 499 | 110-126 |
| 8 x 15 | 120 | 210 | 313 | 66-84 |
| 9 x 17 | 153 | 122 | 180 | 34-49 |
| **10 x 19** | **190** | **83** | **98** | **10-15** |
| 11 x 21 | 231 | 91 | 54 | 6-10 |

Recommended: **10 x 19 = 190 points**, fit to degree 10. Sensors along the arc at 13.0, 29.6, 47.2, 64.3, 81.4 deg
from the +z pole and mirrored (z = +-27.6, 24.6, 19.2, 12.3, 4.2 mm; 8.2-8.7 mm apart along the arc), rotated in
18.9 deg steps. Strict degree-10 sampling is 11 x 21 = 231 points (17.1 deg spacing, pi R / (L + 1/2) ~ 8.5 mm);
the general rule is angular spacing ~180 deg / (L + 1/2), independent of the sphere size.

"Subtracting the design field" means fitting measured minus simulated field and adding the simulation back; it only
helps as far as the real magnet matches the ledger.

**Cylinder 40 x 40 mm surface** (rings plus end caps) needs fewer points: 6 rings x 16 (8 mm apart) + caps = 114
points gives ~50-65 ppm; 5 x 12 + caps = 74 is the minimum (~100 ppm); 4 x 8 = 42 is too sparse (~430). The end caps
are required: without them a plain fit stays at ~700-800 ppm however dense the side.

Do not feed raw sparse points to `shim_design.py`: it then optimizes only at the sensors (one test: it believed
687 ppm, true 1,263). Fit first, then evaluate the fit on the dense check points.

## Accuracy matters more than density

Fit error on a 153-point sphere vs noise per reading (1 uT = 22 ppm at 45 mT):

| Noise | 0 | 1 uT | 2 uT | 5 uT | 10 uT |
|---|---|---|---|---|---|
| Fit error, ppm | 122 | 140 | 187 | 378 | 631 |

Sensor-to-sensor gain mismatch acts the same way: 0.01% -> ~340 ppm, 0.1% -> ~1,750, 1% -> ~7,700 (as large as the
field error being measured). Aim for ~1 uT per reading including mismatch and drift: one sensor moved point to
point, or an array calibrated to ~0.002%, plus a fixed reference sensor, since NdFeB drifts ~0.1%/C (a typical value,
not simulated) and slow drift during a rotating scan looks like field structure. If the arc rotates, use 3-axis
sensors or keep the sensing axis along x.

## Second-round maps

After the first sleeve is in, its cubes are ~20 mm from the sphere and add strong high-order field. Fitting the map
plainly then needs far more points (114-point cylinder: 176-304 ppm instead of ~60). Subtracting the modelled field of
the cubes that were placed brings it back to the first-round level (58-70 ppm at 114 points; with the 190-point sphere,
85-120 ppm including 1 uT noise). A check map after a second, inner layer of shims (cubes ~14 mm from the sphere) has
not been simulated.
