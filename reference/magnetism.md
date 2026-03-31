# Magnetism

*Part of the [GWT Complete Reference](../gwt_complete_reference.md). See also: [Standing Waves](standing_waves.md), [Foundation](foundation.md), [Nuclear](nuclear.md).*

---

## Magnetism from Springs

The lattice is a coupled spring network. Every spring has two responses to displacement:
- **Longitudinal**: compression/extension along the spring axis
- **Transverse**: twist perpendicular to the spring axis

Magnetism IS the transverse twist. Nothing more.

When a spring is displaced longitudinally, it pushes/pulls its neighbors = electric field.
When a spring is twisted transversely, it torques its neighbors = magnetic field.
Same spring, two modes of response. This is why E and B are always perpendicular
and always coupled — they are the two degrees of freedom of one spring.

---

## The T1g Irrep = Magnetic Field

On the d=3 cubic lattice, the gradient of a vector field decomposes under Oh as:

```
T1u x T1u = A1g(1) + Eg(2) + T1g(3) + T2g(3)

A1g = scalar (divergence)     → charge density
Eg  = symmetric traceless     → quadrupole (tidal forces)
T1g = antisymmetric           → CURL = magnetic field
T2g = symmetric off-diagonal  → shear strain
```

The magnetic field is the T1g (antisymmetric) part of the gradient tensor.
It has 3 independent components (Bx, By, Bz) because T1g is 3-dimensional at d=3.

**This is not an analogy.** The curl of the lattice displacement field IS what we
measure as the magnetic field. Maxwell's equations are the wave equations of the
d=3 cubic lattice decomposed into Oh irreps.

---

## Where Magnetic Dipoles Come From

### Kink circulation = nuclear magnetic moment

A kink wrapping a torus creates a circulating flow of field displacement.
This circulation IS a current loop. The loop has an area (the torus cross-section)
and a flow rate (the kink velocity around the loop).

```
Magnetic dipole moment = (flow rate) × (loop area)

For the proton:
  mu_p(bare) = (d^2-1)/d^2 × d = 8/3 nuclear magnetons

  (d^2-1)/d^2 = 8/9: only 8 of 9 Oh channels carry angular momentum
                       (A1g cannot — it is rotationally symmetric)
  d = 3: the kink wraps through 2d = 6 faces, normalized to d directions
```

The bare value 8/3 = 2.667 gets corrected by the pion cloud to 2.7928 (+0.03%).
See [Nuclear](nuclear.md) for the full derivation.

### Harmonic twist = mode magnetic moment

Each harmonic mode of the kink's field has an angular pattern (s, p, d, f).
The angular patterns with l >= 1 carry a twist — a circulating component in
the transverse plane. This twist IS the mode's magnetic moment.

```
l=0 (s, A1g):  no angular twist → no magnetic moment
l=1 (p, T1u):  twist in one plane → magnetic moment along the axis
l=2 (d, T2g+Eg): twist in two planes → larger magnetic moment
l=3 (f):       twist in three planes → larger still
```

The twist direction (clockwise vs counterclockwise) is what we call "spin up"
vs "spin down." Two modes with opposite twist cancel magnetically.

---

## Pairing and Cancellation

This is where the spring picture makes magnetism simple.

### Paired modes: no net twist

Two harmonic modes in the same angular pattern with opposite twist:

```
Mode A:  twisting clockwise in the x-y plane
Mode B:  twisting counterclockwise in the x-y plane
Sum:     zero net twist → zero magnetic moment
```

This is why most materials are not magnetic. In a bound configuration where
all modes are paired, every clockwise twist is cancelled by a counterclockwise one.
The springs twist back and forth but the NET twist is zero.

### Unpaired modes: net twist

When a mode has no partner with opposite twist, the net twist is nonzero.
The springs have a preferred twist direction. This creates a magnetic dipole.

```
One unpaired p-mode:  net twist in one plane → paramagnetic
Multiple unpaired d-modes: net twist in multiple planes → strongly paramagnetic
```

---

## The Three Types of Magnetic Response

### Diamagnetism — all modes paired

Every bound configuration has paired modes. When an external traveling wave
(applied B-field) arrives, it tries to twist all the springs. The paired modes
resist by generating opposing twist (Lenz's law in spring language):

```
External twist applied → springs resist → induced twist opposes applied twist
```

This is UNIVERSAL. Every bound configuration does this. The magnitude is small
because the response is second-order: the springs are already in equilibrium,
and the external twist must fight the existing pairing.

**GWT formula:** The diamagnetic susceptibility scales as the number of paired modes
times the square of the mode radius (how far the twist extends from the kink):

```
chi_dia ~ -(Z/d) × <r^2> × alpha^2

Z = number of paired modes
<r^2> = mean square radius of the harmonic clouds
alpha^2 = coupling strength (how much the springs respond)
```

The negative sign = opposing twist. Always present, always weak.

### Paramagnetism — unpaired modes present

When some modes are unpaired, their net twist can ALIGN with an external field.
This is energetically favorable: a twisted spring in a twisted field is lower energy
when they twist the same way.

```
External twist applied → unpaired mode twist aligns → net twist ADDS to applied field
```

This is stronger than diamagnetism because it is first-order: the unpaired twist
already exists, it just needs to rotate to align. The magnitude scales as:

```
chi_para ~ +N_unpaired × mu_mode^2 / (d × k_B × T)

N_unpaired = number of unpaired modes
mu_mode = magnetic moment per unpaired mode (from the twist amplitude)
T = temperature (thermal agitation randomizes the alignment)
```

The 1/T dependence (Curie's law) is natural: higher temperature = more thermal
spring vibration = harder to maintain alignment.

### Ferromagnetism — exchange coupling locks the twist

In certain configurations (Fe, Co, Ni — all with partially filled d-modes),
something stronger happens. The unpaired d-mode twists on neighboring kinks
don't just independently align with external fields — they LOCK to each other
through the lattice springs.

**The exchange mechanism in spring language:**

Two neighboring kinks each have unpaired d-modes. The springs between them
couple these modes. The coupling energy depends on the RELATIVE twist:

```
E_exchange = -J × cos(theta_12)

theta_12 = angle between twist directions of neighboring kinks
J > 0 (ferromagnetic): parallel twist is lower energy
J < 0 (antiferromagnetic): antiparallel twist is lower energy
```

J comes from the Oh tensor product of the d-mode coupling:

```
(T2g+Eg) x (T2g+Eg) contains A1g

The A1g component = the isotropic (scalar) coupling between neighboring twists.
Its sign depends on the orbital overlap:
  - Less than half-filled d-modes: J > 0 (Hund's rule → ferromagnetic)
  - More than half-filled: J < 0 or weak (antiferromagnetic or paramagnetic)
```

**Why Fe, Co, Ni are ferromagnetic:**

```
Fe: [Ar] 3d^6 4s^2 → 4 unpaired d-modes (less than half of 10 filled)
Co: [Ar] 3d^7 4s^2 → 3 unpaired d-modes
Ni: [Ar] 3d^8 4s^2 → 2 unpaired d-modes

All three have partially filled d-modes with J > 0.
The exchange coupling locks neighboring twists parallel.
Below the Curie temperature, thermal agitation cannot overcome J.
Result: spontaneous macroscopic alignment = permanent magnet.
```

**Why most elements are NOT ferromagnetic:**

- s and p modes: too delocalized. The springs between kinks are too long
  for strong exchange coupling. J is negligible.
- Full d-modes: all paired. No net twist to couple.
- Empty d-modes: no twist to couple.
- d^5 exactly (Mn, half-filled): J ~ 0 (balanced between ferro and antiferro).

Only the narrow window of partially-filled, localized d-modes gives strong
enough J with the right sign.

---

## Electromagnetic Waves = Coupled Spring Oscillations

A traveling wave on the lattice has both longitudinal and transverse components:

```
Longitudinal displacement oscillation = oscillating E-field
Transverse twist oscillation = oscillating B-field

They are 90 degrees out of phase (one is the time derivative of the other).
They propagate together at c = sqrt(k/eta) = 1 (impedance matched).
```

This IS a photon. The E-field is the springs compressing and extending.
The B-field is the springs twisting. Maxwell's equations describe how these
two modes couple on the d=3 lattice:

```
curl E = -dB/dt    →  longitudinal oscillation drives transverse twist
curl B = dE/dt     →  transverse twist drives longitudinal oscillation

div E = rho        →  net compression = charge
div B = 0          →  net twist = zero (twists always form loops)
```

The reason div B = 0 (no magnetic monopoles) is geometric: a twist on a spring
rotates BOTH ends. It cannot create a source or sink — only a loop. This is
the T1g irrep being traceless (antisymmetric matrices have zero trace).

---

## Magnetic Moments from Oh Channels

The magnetic moment of any configuration decomposes into Oh channels:

```
T1u x T1u = A1g(1) + Eg(2) + T1g(3) + T2g(3)

Total channels = 9 = d^2
Channels carrying angular momentum = 8 = d^2 - 1  (everything except A1g)
Channels carrying MAGNETIC moment = T1g = 3 = d    (the antisymmetric part)
```

**Key fractions:**
```
Magnetic / Total = T1g / (T1u x T1u) = 3/9 = 1/d     [magnetic coupling fraction]
Non-magnetic / Total = 1 - 1/d = (d-1)/d = 2/3         [non-magnetic fraction]
Symmetric / Antisymmetric = (Eg+T2g) / T1g = 5/3       [screening ratio]
```

The 1/d = 1/3 magnetic fraction appears throughout:
- Proton magnetic moment: corrections involve 1/d factors
- Electron g-2: the diamagnetic correction is -alpha/(2d-1) where 2d-1 = 5
  counts the symmetric modes that SCREEN the magnetic moment (Eg + T2g)
- Nuclear magnetic moments: bare values use (d^2-1)/d^2 = 8/9 (non-A1g fraction)

See [Nuclear](nuclear.md) for explicit calculations of mu_p, mu_n, g_A, and g-2.

---

## Summary: Magnetism in One Sentence

**Magnetism is the transverse twist mode of the lattice springs, and magnetic
phenomena are what happens when those twists align, cancel, or propagate.**

| Phenomenon | Spring language |
|-----------|----------------|
| Electric field | Longitudinal spring displacement |
| Magnetic field | Transverse spring twist (T1g) |
| Photon | Coupled longitudinal + transverse wave |
| Diamagnetism | Paired twists resist external twist |
| Paramagnetism | Unpaired twists align with external twist |
| Ferromagnetism | Exchange coupling locks neighboring twists parallel |
| Magnetic monopole | Impossible — twists form loops (T1g is traceless) |
| div B = 0 | Antisymmetric tensor has zero trace |
| No magnetic charge | Twist has no source/sink, only circulation |
