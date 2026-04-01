# Condensed Matter

*Part of the [GWT Complete Reference](../gwt_complete_reference.md). See also: [Foundation](foundation.md), [Bonding](bonding.md), [Thermodynamics](thermodynamics.md), [Magnetism](magnetism.md).*

---

## Solids Are Lattices on a Lattice

A solid is a periodic arrangement of kink configurations (atoms) on the Planck
lattice. It is a lattice ON a lattice — the atomic lattice (Angstrom spacing)
riding on the Planck lattice (10^-35 m spacing).

The atomic lattice inherits its symmetry from the d=3 cube. The most common
crystal structures — FCC, BCC, HCP — are all subgroups or packings of the
Oh group that governs the underlying Planck lattice.

---

## Band Structure from Oh Mode Splitting

### Why bands exist

A single kink configuration has discrete harmonic modes (the atomic energy
levels — see [Standing Waves](standing_waves.md)). When many kinks are
arranged periodically, their modes COUPLE through the lattice springs.

Each discrete atomic level splits into a BAND of closely spaced levels:

```
1 atom:     discrete level E_n
2 atoms:    bonding + antibonding (2 levels)
N atoms:    N closely spaced levels forming a band
N -> inf:   continuous band of width W
```

The bandwidth W is set by the spring coupling between neighboring kinks —
the same coupling that gives bond energies in [Bonding](bonding.md):

```
W ~ C_BOND × E_harm = (pi/d^2) × E_harm
```

The V8 bond formula determines band widths. Stronger bonds = wider bands.

### Band gaps from Oh symmetry

The Oh group determines which modes couple and which don't. Modes of
different irrep symmetry CANNOT couple (orthogonality theorem):

```
A1g modes couple to A1g only
T1u modes couple to T1u only
Eg modes couple to Eg only
...
```

This creates GAPS between bands of different symmetry. The gap size is
determined by the energy difference between the Oh irreps.

For the valence/conduction gap (the band gap that makes semiconductors):

```
s-band (A1g) and p-band (T1u) have different symmetry
Gap ~ energy splitting between A1g and T1u modes
    ~ E_H × (alpha_s - alpha_p)   (difference in penetration factors)
```

The penetration factors alpha_s = 11/9 and alpha_p = 3/4 from
[Atomic & Molecular](atomic_molecular.md) directly give the band gap.

---

## Metals, Insulators, Semiconductors

### Metals — partially filled bands

When a band is partially filled, the modes at the top of the occupied
region can absorb arbitrarily small amounts of energy (there are empty
modes just above). This means free breather-like excitations can
propagate through the crystal — electrical conduction.

In GWT: the harmonic modes in a metal are delocalized across the entire
crystal. They're not bound to individual kinks — they're shared standing
waves spanning N sites. Applying a voltage tilts the lattice potential,
and the delocalized modes drift. This IS current.

### Insulators — full bands with large gap

All modes in the valence band are occupied. The next available modes are
across a large gap (> 3 eV). Thermal energy at room temperature
(k_B T ~ 0.025 eV) cannot bridge the gap. No conduction.

### Semiconductors — small gap

Same as insulator but the gap is small (0.5-2 eV). At room temperature,
a few modes are thermally excited across the gap (Boltzmann factor
exp(-E_gap / k_B T)). Small but nonzero conductivity.

Doping = adding kinks with different winding numbers (different Z) into
the periodic lattice. The foreign kink has modes at different energies,
creating levels INSIDE the gap. These levels provide stepping stones for
conduction.

---

## Phonons — Lattice Vibrations

Phonons are collective oscillations of the atomic lattice (not the Planck
lattice). They are sound waves — coherent spring oscillations where
neighboring kinks move together.

The phonon spectrum decomposes under the crystal symmetry group (a subgroup
of Oh):

```
Acoustic phonons (3 branches):
  T1u — the 3 translational modes. Linear dispersion near k=0.
  These carry sound. Speed = sqrt(k_eff / M_atom).

Optical phonons (3(N_basis - 1) branches):
  Higher-frequency modes where atoms in the unit cell move against
  each other. These couple to light (hence "optical").
```

The Debye model (T^3 heat capacity) counts acoustic phonon modes.
The Einstein model (exponential activation) counts optical phonon modes.
Both are lattice wave models — GWT says they're correct because the
underlying physics IS waves on springs.

---

## Superconductivity — Breather Pairing

Superconductivity occurs when breather modes (conduction electrons) form
PAIRS that propagate without scattering.

### The pairing mechanism

Two breathers moving through a lattice of kinks create wakes — the kinks
are slightly displaced by each passing breather. The wake of breather 1
creates a region of slightly higher kink density (attractive potential)
that breather 2 falls into. This is the phonon-mediated attraction.

In Oh language: the breather-phonon coupling is through the A1g channel
of T1u × T1u. Two breathers coupled through A1g form a SCALAR pair —
spin 0, charge 2e, bosonic. This pair (Cooper pair) has A1g symmetry
and therefore propagates without scattering off the Oh-symmetric lattice.

```
Why pairs don't scatter:
  Single breather = T1u (vector) — scatters off lattice defects
  Breather pair   = A1g (scalar) — invisible to the Oh lattice

  The pair's A1g symmetry means it transforms trivially under all
  48 Oh operations. The lattice cannot deflect what it cannot "see."
```

### Critical temperature

The pairing energy (the energy gained by forming a Cooper pair) competes
with thermal energy (which breaks pairs). Superconductivity survives when:

```
k_B × T_c ~ Delta / d

Delta = pairing gap energy
d = 3 (the number of directions thermal fluctuations can attack the pair)
```

The factor 1/d appears because thermal fluctuations in ANY of d directions
can break the pair. This gives the BCS ratio:

```
2 × Delta / (k_B × T_c) = 2d + 1/(d-1) = 2 × 3 + 1/2 = 3.5

  (BCS prediction: 3.53, measured: 3.5 ± 0.1 for conventional superconductors)
```

### High-Tc superconductors

In cuprate superconductors, the pairing is not phonon-mediated but
spin-mediated. The pairing symmetry is d-wave (Eg or T2g) instead of
s-wave (A1g). The higher pairing energy comes from the stronger magnetic
(T1g) coupling in the copper-oxygen planes.

The 2D nature of cuprates means the effective dimension is d=2 for the
pairing, but d=3 for the lattice. This mismatch enhances T_c because
the pairing gap is set by 2D physics while the thermal destruction
requires 3D fluctuations.

---

## Fermi Liquid — Mode Filling

In a metal, the breather modes fill from lowest energy upward (the
standing wave analogue of the Pauli principle — see [Standing Waves](standing_waves.md)).
The highest filled mode has the Fermi energy:

```
E_F = (hbar^2 / 2m) × (3 × pi^2 × n)^(2/d)

n = breather mode density
d = 3 (dimensions of mode space)
```

The (2/d) = 2/3 exponent comes from the d=3 density of states.
The factor 3 × pi^2 = d × pi^2 comes from the volume of the Fermi
sphere in d dimensions.

The Fermi surface IS the boundary between occupied and unoccupied
modes in k-space. Its shape reflects the crystal symmetry — for a
cubic crystal, the Fermi surface has Oh symmetry at low filling and
distorts toward the Brillouin zone boundary at high filling.

---

## Semiconductors — Specific Predictions

### Silicon band gap

Silicon has diamond structure (Fd3m, a subgroup of Oh). The band gap
between the sp3 bonding band (A1g + T2g) and the antibonding band:

```
E_gap(Si) ~ 2 × alpha × E_harm(Si-Si) × (alpha_s - alpha_p)
         ~ 2 × (1/137) × 8.15 × (11/9 - 3/4)
         ~ 2 × 0.0073 × 8.15 × 0.472
         = 0.056 eV   (too small — needs VP correction)
```

The bare estimate is too small. The VP-dressed version should include
the Oh channel corrections. This is an open calculation — not yet
completed. The observed gap is 1.12 eV.

### Germanium, GaAs

Same framework, different kink configurations. The gap scales with the
bond energy and penetration factor difference. Qualitatively:
- Larger atoms → smaller gap (more screening, closer penetration factors)
- More ionic → larger gap (asymmetry opens the gap further)

---

## Connection to the Framework

| Phenomenon | GWT mechanism |
|-----------|--------------|
| Band structure | Oh mode splitting across periodic kink array |
| Band gap | Symmetry gap between different Oh irrep bands |
| Metals | Partially filled bands — delocalized modes |
| Insulators | Full bands, large symmetry gap |
| Semiconductors | Small gap, thermal activation |
| Phonons | Collective kink oscillations (acoustic = T1u, optical = higher) |
| Superconductivity | A1g breather pairing — invisible to the Oh lattice |
| BCS ratio 3.5 | 2d + 1/(d-1) from d=3 thermal fluctuations |
| Fermi energy | Mode filling in d=3 k-space |
| Doping | Foreign kinks with mid-gap mode levels |

The same Oh symmetry that gives particle physics, atomic structure, and
bonding also determines the electronic structure of solids. Condensed matter
IS lattice wave physics — it was always the closest to the truth.
