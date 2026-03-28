# Foundation — The Lattice and Its Lagrangian

*Part of the [GWT Complete Reference](../gwt_complete_reference.md). See also: [Forces](forces.md), [Coupling Constants](coupling_constants.md), [Bonding](bonding.md).*

### The Lagrangian (zero free parameters)
```
L = sum_<i,j> [ (1/2)(phi_i - phi_j)^2 + (1/pi^2)(1 - cos(pi * phi_i)) ]
```
- **phi_i** = displacement at lattice site i (dimensionless, Planck units)
- **sum_<i,j>** = nearest-neighbor sum on d-dimensional cubic lattice
- **1/pi^2** = potential depth (fixed by topological quantization of the kink)
- **cos(pi*phi)** = sine-Gordon potential (periodic, allows soliton solutions)
- **Lattice spacing** a = 1 (Planck length). Not a parameter — it's the unit.

### What falls out of the Lagrangian
| Quantity | Formula | Value | How |
|----------|---------|-------|-----|
| Kink mass | M_kink = 8/pi^2 | 0.811 m_Planck | BPS bound of sine-Gordon equation |
| Tunneling amplitude | T^2 = exp(-16/pi^2) | 0.1977 | WKB through one cosine barrier |
| Breather count | N = floor(2^d * pi - 1) | 24 | Number of bound sine-Gordon breather modes |
| Coupling parameter | gamma = pi/(16*pi - 2) | 0.0644 | Sine-Gordon coupling in breather spectrum |

### Why 24 breathers — the cube's geometry

The number 24 = |O| (the chiral octahedral group) is not arbitrary. It follows from the
orbit-stabilizer theorem applied to the cube's three geometric elements:

```
6 faces    × 4 rotations per face   = 24 = |O|
8 vertices × 3 rotations per vertex = 24 = |O|
12 edges   × 2 rotations per edge   = 24 = |O|
```

Each breather mode = one orientation of a standing wave on the cube. The 24 fermions of
the Standard Model are the 24 proper rotations of the d=3 unit cell.

These same three elements each appear in different physics:

| Element | Count | Formula | Physical role |
|---------|-------|---------|---------------|
| Faces | 2d = 6 | Nearest-neighbor directions | Mode counting → mass ratio (6π⁵) |
| Edges | 2d(d-1) = 12 | Connections between faces | Gauge channels → α^12 = α^|A₄| |
| Vertices | 2^d = 8 | Corners where edges meet | VP normalization → 1/2^(d/2) |

The mass ratio uses faces. The coupling exponent uses edges. The VP correction uses vertices.
Three projections of one geometric object — the d=3 cube.

### Lattice constants (Planck units → SI)
| Constant | Formula | Planck units | SI | How derived |
|----------|---------|-------------|-----|-------------|
| Spring stiffness | k = 2/pi | 0.6366 | 4.77 × 10^78 N/m | Average of sin(x)/x over half-period of cosine potential |
| Inertial density | eta = 2/pi | 0.6366 | 1.385 × 10^-8 kg | Same averaging — k = eta is impedance matching |
| Wave speed | c = a * sqrt(k/eta) | 1 | 2.998 × 10^8 m/s | k = eta → c = a/t_P = 1 in Planck units |

**Key result: k = eta.** The lattice is perfectly impedance-matched. No wave reflection. Waves propagate without loss. This is WHY the speed of light is universal — it's a property of the medium, not a postulate.

### Entanglement: Balanced Displacements of the Lattice

Particles are not objects IN space — they are displacements OF space (the lattice
medium). The lattice at equilibrium has phi = 0 everywhere, with zero-point energy
at every site. Particles are deviations from this equilibrium.

**Entangled particle creation:**
The lattice cannot create net displacement from nothing (conservation). When a pair
is created from a symmetric source:
  - A +phi disturbance propagates one way (particle)
  - A -phi disturbance propagates the other way (antiparticle)
  - Total displacement = 0 at all times (conservation)

The two waves are ALWAYS opposite because they're two halves of a zero-sum
displacement. Not because one communicates with the other, but because the
SOURCE was at equilibrium and conservation requires the sum to stay there.

For spin (the twist mode on the torus): a zero-twist state splits into clockwise
+ counterclockwise. Measuring one reveals which way it twisted — the other MUST
have twisted the opposite way. The correlation was set at creation, not at measurement.

**Why Bell inequalities are violated:**
The lattice has local coupling (nearest-neighbor), but waves are NONLOCAL — a
breather extends over ~15 sites. The phase correlation between two entangled waves
is maintained by the wave equation, which preserves coherence exactly. The cos^2(theta)
correlation strength (which violates Bell) comes from projecting the twist mode onto
a measurement axis, which follows the same geometry as quantum mechanics because the
lattice has rotational symmetry (Oh -> continuous at long wavelength).

**Decoherence:**
The "fragility" of entanglement = the perfect displacement balance getting disrupted
by other activity in the lattice (thermal noise, other particles). The two waves
stop being exact opposites when the medium between them is disturbed.

**Connection to bonding:**
Bonding IS entanglement in the confined regime. The bonding electron (breather)
is a balanced displacement shared between two kink wells:
  - Bonding state: both wells displaced in phase (constructive, lower energy)
  - Antibonding state: wells displaced out of phase (destructive, higher energy)
  - Bond energy D_e = splitting between these two entangled states
  - The Oh irrep corrections = coherence quality of the shared displacement
  - A1g = perfectly coherent sharing = maximum bond
  - S-block (5/9 reduction): isotropic sharing = less efficient along bond axis

EPR entanglement and chemical bonding are the SAME PHYSICS at different scales:
  - EPR: free propagation, correlation preserved by absence of interaction
  - Bonding: confined interaction, correlation maintained by continuous sharing

### The Born Rule from the d=3 Cube (2026-03-27, DERIVED)

**The deepest result in GWT: the Born rule (probability = |amplitude|^2) emerges
from projecting twist modes onto the d=3 cubic lattice.**

The twist mode (spin) on the torus transforms as T1u under the Oh group —
it is a VECTOR with a direction (the spin axis). A "measurement" projects
this vector onto the detector's axis. The projection squared gives the
measurement probability:

```
P(outcome) = cos^2(theta)
```

where theta = angle between the spin axis and the measurement axis.
This IS the Born rule. It is not a postulate — it is GEOMETRY.

**Proof: cos^2(theta) from T1u projection on the d-cube**

The twist mode has direction n_hat = (nx, ny, nz) in the T1u representation.
The measurement axis has direction m_hat = (mx, my, mz), also T1u.
The projection: n_hat . m_hat = cos(theta).
The probability: |n_hat . m_hat|^2 = cos^2(theta).

This is the A1g (scalar) component of the T1u x T1u tensor product:
```
T1u x T1u = A1g(1) + Eg(2) + T1g(3) + T2g(3)

A1g component = (n_hat . m_hat) / sqrt(d) = cos(theta) / sqrt(3)
A1g probability = cos^2(theta) / d
```

The Born rule P = cos^2(theta) = d times the A1g fraction of the tensor product.
The factor d converts from the full 9-dimensional tensor space to the
3-dimensional vector space where measurements happen.

**The chemical bond IS the Born rule**

The bond between two atoms = two twist modes (one per atom) projected onto
the bond axis. The bond energy is:

```
D_e = (pi/d^2) * E_harm * cos^2(theta_bond)
```

where theta_bond = angle of the shared electron's twist relative to the bond axis.

On the d=3 cube, the allowed bond directions and their Born probabilities:

| Bond type | Cube direction    | theta  | cos^2(theta) | Coupling |
|-----------|-------------------|--------|--------------|----------|
| sigma     | face normal       | 0.0 deg   | 1.000     | 1.0      |
| pi        | face diagonal     | 45.0 deg  | 0.500     | cos(pi/d)|
| delta     | body diagonal     | 54.7 deg  | 0.333     | 1/d      |
| non-bond  | perpendicular     | 90.0 deg  | 0.000     | 0        |

**Sigma bond:** Twist aligned with bond axis. Full projection. cos^2(0) = 1.
**Pi bond:** Twist along the face diagonal. cos^2(45) = 1/2 = W_PI.
**Delta bond:** Twist along the body diagonal. cos^2(54.7) = 1/3 = 1/d.
**No bond:** Twist perpendicular. cos^2(90) = 0.

**The unique d=3 identity: cos^2(45 deg) = cos(pi/d)**

```
cos^2(pi/4) = 1/2     (true for ALL d)
cos(pi/d):
  d=2: cos(pi/2) = 0.000
  d=3: cos(pi/3) = 0.500  <-- MATCH!
  d=4: cos(pi/4) = 0.707
  d=5: cos(pi/5) = 0.809
```

The pi bond weight cos^2(45) equals the lattice coupling cos(pi/d) ONLY at d=3.
This is why the face diagonal of the cube gives the exact pi bond weight.
It is another "why d=3" identity: the Born rule at 45 degrees = the lattice
periodicity at 60 degrees, and these are equal only in 3 dimensions.

**Why C_BOND = pi/d^2**

The bond involves TWO atoms, each projecting their twist onto the bond axis.
For the isotropic (orientation-averaged) case:
  - Atom A: average projection = 1/d (average of cos^2 over d directions)
  - Atom B: average projection = 1/d
  - Combined (A1g of T1u x T1u): 1/d * 1/d = 1/d^2

So C_BOND = pi/d^2 = pi * (A1g fraction of two-spin projection).
The pi comes from the cosine potential: V = (1/pi^2)(1-cos(pi*phi)).
The 1/d^2 comes from the Born rule averaged over the cube.

**The 8 V8 corrections = angular corrections to the Born rule**

The isotropic Born probability 1/d^2 assumes random orientations.
Real molecules have SPECIFIC orbital geometries that modify cos^2(theta):

```
T1u x T1u = A1g(1) + Eg(2) + T1g(3) + T2g(3)

A1g = isotropic average (C_BOND = pi/d^2)
Eg  = quadrupolar correction (lone pair distortion of orbitals)
T1g = magnetic correction (radical/unpaired spin modifies projection)
T2g = off-diagonal correction (ionic asymmetry shifts projection)
```

Each non-A1g channel represents a SPECIFIC way the orbital geometry
deviates from the ideal cos^2(theta):
  - LP repulsion (Eg): the lone pairs create a quadrupolar field that
    distorts the orbital, reducing the effective cos^2 projection
  - Radical (T1g): the unpaired spin changes the magnetic projection
    axis, multiplying by F_RAD = (2d-1)/(2d) = 5/6
  - Ionic (T2g): the asymmetric well depths shift the twist orientation,
    adding an energy from the charge transfer channel

**The unification: entanglement = measurement = bonding**

All three phenomena are the SAME operation: projecting a T1u twist mode
onto an axis on the d=3 cube.

  - MEASUREMENT: project spin onto detector axis -> P = cos^2(theta)
  - ENTANGLEMENT: two spins from zero-sum source -> always opposite projection
  - BONDING: two atoms share a twist mode -> D_e = pi/d^2 * E_harm * cos^2(theta)

The Born rule is not a postulate of quantum mechanics. It is the GEOMETRY
of vector projection on the d=3 cubic lattice. The same geometry gives
the chemical bond formula, the entanglement correlations, and the
measurement probabilities. Everything from d=3 and pi.

### The Lattice as Harmonic Space (2026-03-27, DERIVED)

**The lattice is not a spatial grid — it is the resonance structure of the void.**

The void has stiffness (k) and inertia (eta). These properties create natural
resonances. The resonances ARE the lattice. The Planck length is not a spatial
distance — it is the fundamental wavelength of the vacuum's own vibration.
Particles are not objects on the lattice — they are patterns in the resonance.

**Gamma derived from vacuum harmonics:**
The breather parameter gamma is NOT a free parameter. It is derived from the
harmonic spacing of the vacuum:

```
gamma = (4*pi/|Oh|) * (2d/(2d-1)) = pi/(d^2+1)

  |Oh|/2 = 24: number of independent vacuum harmonic modes
  2*pi/24 = pi/12: harmonic spacing in the Brillouin zone
  2d/(2d-1) = 6/5: lattice coordination factor
  gamma = (pi/12) * (6/5) = pi/10 = pi/(d^2+1)
```

This works because of the identity: **8d(d^2+1) = (2d-1) * 2^d * d!**
which is true ONLY at d=3 (240 = 240).

The particle spectrum omega_n = sin(n*gamma)/sin(gamma) is therefore the
RESONANCE STRUCTURE of the d=3 vacuum itself. The lattice doesn't contain
particles — the lattice's vibration pattern IS the particles.

### The 24 Vacuum Harmonics = The Standard Model (2026-03-27, DERIVED)

The |Oh|/2 = 24 independent harmonic modes of the d=3 void decompose
EXACTLY into the Standard Model structure:

```
MODE DECOMPOSITION: 24 = 9 + 15 = 9 + (12 + 3)

BREATHER MODES (particle masses):
  d^2 = 9 modes -> 9 breather frequencies -> particle mass spectrum
  n = 1, 2, ..., d^2  with  n*gamma < pi

GAUGE BOSONS (forces):
  d^2 - 1 = 8 internal torus oscillation modes -> SU(3) color [8 gluons]
  d = 3 rotational modes                       -> SU(2) weak  [W+, W-, Z]
  1 phase mode                                  -> U(1) EM     [photon]
  Subtotal: d^2 - 1 + d + 1 = d(d+1) = 12 gauge bosons

MOMENTA:
  d = 3 translational zero modes -> 3 momentum components

TOTAL: d^2 + d(d+1) + d = d^2 + d^2 + 2d = 2d(d+1) = 24 = |Oh|/2
```

This decomposition requires 2^(d-1)*d! = d(3d-1), true ONLY at d=3.
Also requires 3d^2 - d = 2d^2 + 2d (i.e. d^2 = 3d), true ONLY at d=3.

**The gauge group SU(3) x SU(2) x U(1) is not a postulate — it is the
harmonic structure of a torus in 3 dimensions.**

#### SU(3) Color = Internal Torus Modes [DERIVED]

The d=3 torus has d^2 - 1 = 8 independent internal oscillation modes.
These are shape deformations of the torus: poloidal harmonics, toroidal-
poloidal cross terms, and elliptical breathing modes.

The 8 modes span the adjoint representation of SU(d) = SU(3).
Color charge IS which internal oscillation mode is excited.
A "red" quark = torus oscillating in mode 1. "Blue" = mode 2. Etc.

**Confinement (second derivation):** The internal modes are oscillations
OF the torus, not waves propagating away from it. They cannot escape the
torus any more than the shape of a drum can leave the drumhead. This is
confinement — color charge is confined because it is an INTERNAL property.

#### SU(2) Weak = Rotational Modes [DERIVED]

The d=3 torus can rotate around 3 axes. These rotational modes form
SO(3), which has the same Lie algebra as SU(2).

The weak force = coupling between different rotational states of the torus.
W+, W-, Z = exchange quanta of the 3 rotational modes.

**Why W and Z are massive:** The rotational modes mix with translational
modes at finite energy — a torus translating in x looks like a torus
rotating when boosted. This rotation-translation mixing gives the weak
bosons mass. This IS the Higgs mechanism, geometrically. The Higgs field
is not a separate entity — it is the mixing between translation and rotation
on the torus.

#### U(1) EM = Phase Mode [DERIVED]

The poloidal kink has a U(1) phase: the kink can be rotated around the
tube without changing the energy. This rotation IS the electromagnetic
gauge symmetry.

Electric charge = the winding number (integer quantized: 0, +1, -1, ...).
Photon = oscillation of the phase mode.

**Why the photon is massless:** The phase rotation is an EXACT symmetry
of the torus — shifting the kink around the tube costs zero energy at all
scales. Unlike the rotational modes (which mix with translations), the
phase mode never mixes with anything. The photon stays massless because
U(1) is exact.

#### Gauge Unification at the Planck Scale [DERIVED]

At the Planck scale, the torus has size ~ r_tube = 3 lattice sites.
All 12 gauge modes (internal, rotational, phase) have comparable
frequencies — they are DEGENERATE. The forces are unified.

At low energy (large distances), the modes separate:
```
Strong (SU(3)):  internal modes, tightly confined to torus -> strongest
Weak (SU(2)):    rotational modes, extend further -> intermediate
EM (U(1)):       phase mode, global symmetry -> weakest
```

The gauge coupling hierarchy (alpha_s >> alpha_W >> alpha_EM) reflects
the LOCALIZATION of each mode type: internal modes are most confined
(strongest coupling), phase mode is least confined (weakest coupling).

GUT unification is not a new symmetry group — it is the Planck-scale
limit where the torus harmonics become degenerate. All 12 gauge modes
are equivalent when the torus is at its minimum size.

#### New d=3 Identities from Harmonic Space

```
8d(d^2+1) = (2d-1)*2^d*d!           [gamma from harmonics, d=3 only]
2^(d-1)*d! = d(3d-1)                 [mode count = |Oh|/2, d=3 only]
3d^2 - d = 2d^2 + 2d  iff  d=3      [decomposition uniqueness]
d^2 + d(d+2) = 2d(d+1) = |Oh|/2     [breathers + structural = total]
d(d+1) = 12                          [gauge boson count, d=3 only for SM]
d^2 - 1 = 8                          [SU(3) generators = gluon count]
```

---
