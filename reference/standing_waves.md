# Standing Wave Structure of Atoms

*Part of the [GWT Complete Reference](../gwt_complete_reference.md). See also: [Foundation](foundation.md), [Toroidal Physics](toroidal_physics.md), [Atomic & Molecular](atomic_molecular.md).*

---

## The Two Field Configurations

The Lagrangian L = (1/2)(dphi)^2 + (1/pi^2)(1-cos(pi*phi)) has a periodic potential
with valleys at phi = 0, 1, 2, 3, ... Two fundamentally different field configurations
exist on this landscape:

### Kink — a permanent transition between valleys

```
phi(x):

  1  ──────────────────────         valley 1
                          |
                          |   <-- transition region (~3 lattice sites wide)
                          |
  0  ──────────────────────         valley 0

         left            right
```

The field on the left sits in valley 0. The field on the right sits in valley 1.
The smooth transition between them IS the kink. It cannot be removed without
unwinding the entire field behind it — this is TOPOLOGICAL permanence.

- **Energy** = integral of the gradient across the transition = M_kink = 8/pi^2
- **Width** = ~3 lattice sites (set by the balance of gradient vs potential)
- **Stability** = absolute. The field to the left and right are in different valleys.
  No local rearrangement can undo this. This is why protons don't decay.

On the d=3 lattice, the kink wraps around a torus (closed loop through 2d = 6 faces).
This toroidal kink = proton. Its poloidal circulation carries angular momentum and
creates a magnetic dipole through the torus hole.

### Breather — a temporary oscillation within one valley

```
phi(x):

  0.3  . . . . /\
              /    \
  0   ──────/      \──────     valley 0
              \    /
 -0.3  . . . . \/

         oscillation
```

The field oscillates back and forth around a single valley minimum. It never
reaches the next valley. The oscillation is localized in space and periodic in time.

- **Energy** = amplitude^2 (harmonic oscillator scaling)
- **Width** = ~15 lattice sites (broader than the kink)
- **Stability** = conditional. The oscillation can radiate away if disturbed.
  Breather + anti-breather = flat field (annihilation).

A free breather propagating on the lattice = what we detect as an "electron" in
a beam, a wire, or a detector.

---

## There Are No Electrons in Atoms

This is the central claim. Inside an atom, there are no separate particles called
electrons orbiting a nucleus. There is ONE field configuration:

**A kink (torus) whose surrounding field vibrates at specific harmonic frequencies.**

What chemistry calls "electrons" are harmonic modes of the kink's own field.

### The kink creates a potential well

When the field transitions between valleys (the kink), it creates a potential
well in its vicinity. Linearizing the Lagrangian around the kink solution gives
the Poschl-Teller potential:

```
U(x) = -V_0 / cosh^2(x/w)

V_0 = 2*pi*Z (depth scales with nuclear charge)
w = kink width (~3 lattice sites)
```

This well has EXACT bound states. The bound state energies are:

```
E_n = -E_H * Z_eff^2 / n^2

where E_H = alpha^2 * m_e / 2 = 13.6 eV (from the Lagrangian)
```

The quantum number n is NOT an input — it is the bound state INDEX of the
Poschl-Teller well. n=1 is the deepest mode, n=2 the next, etc.

### Harmonics, not orbits

Each bound state is a standing wave pattern in the kink's field:

```
n=1:  one lobe, tight to the kink surface        (capacity: 2 modes)
n=2:  two lobes, extends further                  (capacity: 8 modes)
n=3:  three lobes, extends further still          (capacity: 18 modes)
n=4:  four lobes                                  (capacity: 32 modes)
```

The capacity 2*n^2 per level comes from angular decomposition on the d=3 cube:

```
l=0 (s): A1g irrep — 1 angular pattern, 2 twist orientations  = 2
l=1 (p): T1u irrep — 3 angular patterns, 2 twist orientations = 6
l=2 (d): T2g+Eg    — 5 angular patterns, 2 twist orientations = 10
l=3 (f): A2u+T1u+T2u — 7 angular patterns, 2 twist orientations = 14
```

The factor 2 = two opposite twist directions per pattern (what QM calls "spin up/down").
The 2l+1 = number of distinct angular patterns on the d=3 cube.

**These are not particles in slots. They are the ways the kink's field can vibrate.**

### Why harmonics pair (the Pauli principle in wave language)

Two oscillations in the same angular pattern MUST have opposite twist.
This is not a mysterious exclusion principle — it is wave physics:

```
Same frequency + same angular pattern + same twist = constructive interference
→ amplitude doubles → energy quadruples → the mode splits into two separate frequencies
→ the "same state" doesn't exist as a stable configuration
```

Two oscillations with opposite twist = destructive interference of the twist component.
The angular energy stays finite. The mode is stable with exactly 2 oscillations.

A third oscillation in the same pattern has no remaining twist direction.
It MUST go to a different angular pattern or a different n level.
This is the Pauli exclusion principle: a consequence of wave interference, not a postulate.

---

## What "Absorption" and "Emission" Actually Are

### Absorption = field reorganization

When a bound configuration absorbs energy (e.g., from a passing traveling wave):

```
BEFORE:  kink field vibrating in mode pattern A
ENERGY IN: traveling wave arrives, deposits energy into the field
AFTER:   kink field vibrating in mode pattern B (higher harmonic)
```

Nothing new arrives. The existing field around the kink reorganizes
into a higher harmonic. The energy was carried by the traveling wave; now it is
stored in the standing wave pattern.

### Emission = field relaxation

When a bound configuration emits:

```
BEFORE:  kink field vibrating in mode pattern B (excited harmonic)
AFTER:   kink field vibrating in mode pattern A (lower harmonic)
         + traveling wave propagating away
```

The field relaxes to a lower harmonic. The excess energy propagates away as a
traveling wave on the lattice. The field around the kink simply changed its
vibration pattern, and the difference propagated outward.

### Ionization = mode detachment

When enough energy is deposited, a harmonic mode gains enough amplitude to
detach from the kink's well and propagate as a free breather:

```
BEFORE:  kink field vibrating in mode n (bound standing wave)
ENERGY IN: exceeds the well depth for mode n
AFTER:   kink field with mode n absent
         + free breather propagating away (= "free electron")
```

The "electron" did not exist as a separate entity inside the atom.
It was created at the moment of detachment — a bound harmonic mode
became a free propagating breather.

### Why bound harmonics don't radiate

Standard physics struggles to explain why electrons in atoms don't radiate
(the classical stability problem that motivated quantum mechanics).

In GWT there is no puzzle: **standing waves don't radiate.** A vibrating
guitar string at a fixed harmonic does not emit traveling waves along the string.
The nodes are fixed, the energy stays in the mode. Only TRANSITIONS between
harmonics (changes in the standing wave pattern) produce traveling waves.

---

## Bound Field Configurations

A free kink propagating on the lattice = proton or neutron.
A free breather propagating on the lattice = electron.
These are real, independent excitations.

When kinks and breathers interact, they form a BOUND configuration — a single
field object with two aspects:

- **Topological** component (kink winding = permanent, carries charge and nuclear mass)
- **Dynamical** component (harmonic modes of the field around the kink = changeable)

The harmonic modes are not separate from the kink — they ARE the kink's field
vibrating. The field gains energy and shifts its harmonic pattern for stability.
When the configuration breaks apart (ionization, decay), the harmonics detach
and propagate as free breathers again.

### The donut picture

Visualize the bound configuration as a donut (kink/torus) surrounded by concentric vibration clouds:

```
         . . . n=3 cloud (18 modes max) . . .
       .   . . n=2 cloud (8 modes max) . .   .
     .   .                               .   .
    .   .     ████████████████████        .   .
   .   .    ██                    ██       .   .
   .  .    ██    KINK (torus)      ██      .  .
   .  .    ██    topological       ██      .  .
   .  .    ██    permanent         ██      .  .
   .   .    ██                    ██       .   .
    .   .     ████████████████████        .   .
     .   .   n=1 cloud (2 modes max)  .   .
       .   . . . . . . . . . . . . .   .
         . . . . . . . . . . . . . . .
```

- The donut is topological — it cannot be unwound
- The clouds are dynamical — modes can be added or removed
- Inner clouds screen outer ones (they sit between the donut and the outer cloud)
- Paired modes in the same cloud have opposite twist (cancel magnetically)

### Shell structure = harmonic structure

The periodic table IS the harmonic structure of the kink's potential well:

```
Period 1: n=1, capacity 2          → H, He
Period 2: n=2, capacity 8          → Li through Ne
Period 3: n=3, capacity 8 (s+p)    → Na through Ar
Period 4: n=3d + n=4, capacity 18  → K through Kr
...
```

The "aufbau principle" (fill lowest energy first) = fill deepest harmonics first.
The "Hund rule" (maximize unpaired twist) = twist modes repel when parallel,
attract when anti-parallel (wave interference).

All of it follows from the potential well shape, which follows from the Lagrangian.

---

## Connection to Spectral Lines and Ionization

The quantum defect formula (derived in [Atomic & Molecular](atomic_molecular.md))
gives the exact well shape per subshell:

```
E(n, l) = -E_H / (n - delta_l)^2

delta = quantum defect = how much the harmonic deviates from the pure 1/n^2 pattern
```

The deviation exists because inner harmonics modify the well shape seen by outer ones.
This is SCREENING in wave language: inner standing waves change the effective potential
for outer standing waves.

```
The complete chain:

Lagrangian → kink → Poschl-Teller well → harmonic modes
         → quantum defects (delta) → effective well depth (Z_eff)
         → ionization energies → harmonic mean energy (E_harm)
         → bond energy D_e = pi/d^2 * E_harm
```

Everything from one Lagrangian. The atom is one field configuration.
The "electrons" are its harmonics. The periodic table is its mode structure.
Spectral lines are transitions between harmonics. Bond energy is the coupling
between two kinks' harmonic fields.

---

## Free vs Bound: When Does a "Particle" Exist?

| Regime | What it is | Wave description |
|--------|-----------|-----------------|
| Free kink | Propagating topological defect | Independent soliton (proton, neutron) |
| Free breather | Propagating oscillation on lattice | Independent traveling wave (electron) |
| Bound harmonic | Mode of the kink's field | Standing wave in the kink's well |
| Transition | Mode gaining/losing energy | Standing wave <-> traveling wave conversion |

A particle exists as an independent entity ONLY when free. Inside a bound
configuration, asking "where is the electron" is like asking "where is the 3rd
harmonic on a guitar string" — it's everywhere the string vibrates in that pattern.

This is not an interpretation. It is what the Lagrangian says. The kink creates
a well. The well has modes. The modes are standing waves. There are no particles
orbiting anything.
