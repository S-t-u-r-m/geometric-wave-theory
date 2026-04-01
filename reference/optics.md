# Optics

*Part of the [GWT Complete Reference](../gwt_complete_reference.md). See also: [Foundation](foundation.md), [Relativity](relativity.md), [Standing Waves](standing_waves.md).*

---

## Light Is a Lattice Wave

A photon is a traveling wave on the lattice — a propagating oscillation of the
field phi that never reaches the cosine barrier. In vacuum, it propagates at:

```
c = a × sqrt(k/eta) = 1   (Planck units)

k   = 2/pi   spring stiffness
eta = 2/pi   inertial density
```

Impedance matched: k = eta. No reflection, no dispersion, no loss. This is why
light travels at c in vacuum — the medium is perfect.

---

## Refraction — Why Light Slows Down in Matter

A material is a region of lattice filled with kinks (nuclei) and their harmonic
modes (bound field vibrations). These modify the local spring properties:

### The mechanism

A kink sitting at a lattice site changes the local potential from the bare cosine
to the Poschl-Teller well. The springs near the kink are STIFFER (the field
gradient is larger in the transition region). But the kink also adds inertial
density — the field displacement carries energy.

The net effect on wave speed:

```
v = a × sqrt(k_eff / eta_eff)
```

In a material:
- k_eff increases (stiffer springs near kinks)
- eta_eff increases MORE (the harmonic modes add inertial density)
- Net: v < c

The ratio c/v = n (refractive index):

```
n = c/v = sqrt(eta_eff / eta_0) × sqrt(k_0 / k_eff)
```

Since eta_eff grows faster than k_eff, n > 1 for all normal materials.

### Why n depends on frequency (dispersion)

The harmonic modes around each kink have specific resonant frequencies
(the bound states of the Poschl-Teller well). When the traveling wave's
frequency approaches a mode resonance, the coupling between the wave and
the mode increases dramatically:

```
Far from resonance:  wave passes through, slight slowing (n ~ 1)
Near resonance:      wave strongly couples to the mode, large slowing (n >> 1)
At resonance:        wave is ABSORBED (converted to mode excitation)
Above resonance:     wave couples with opposite phase (n can drop below 1)
```

This is the origin of the Sellmeier equation and all chromatic dispersion.
The resonances are the spectral lines — the same transitions we predict
with the quantum defect formula.

### Snell's law from impedance mismatch

At a boundary between two media (different kink densities), the wave
encounters an impedance mismatch. The wave equation requires continuity
of phi and dphi/dx at the boundary. This gives:

```
n_1 × sin(theta_1) = n_2 × sin(theta_2)
```

This is Snell's law. It follows from the wave equation boundary conditions —
no additional physics needed. The angle changes because the wavefronts
must remain continuous across the boundary while the wave speed changes.

### Total internal reflection

When n_1 > n_2 (going from slower medium to faster medium) and the angle
exceeds the critical angle:

```
theta_c = arcsin(n_2 / n_1)
```

the wave cannot satisfy the boundary conditions with a transmitted wave.
All energy reflects. The wave on the far side becomes evanescent — an
exponentially decaying field that doesn't propagate. In GWT: the springs
on the far side can't oscillate fast enough to match the incoming wave's
transverse momentum.

---

## Reflection — Impedance Mismatch

In vacuum, k = eta (impedance matched, zero reflection). At a material
boundary, the impedance changes. The reflection coefficient:

```
R = |(Z_1 - Z_2) / (Z_1 + Z_2)|^2

Z = sqrt(k_eff × eta_eff)   (wave impedance)
```

This is identical to the formula for reflection of waves on a string, sound
in air, or electromagnetic waves at a boundary. Same physics — waves on
springs encountering a change in spring properties.

Metals have high reflectivity because their free breather modes (conduction
electrons) create a large impedance mismatch. The free modes oscillate in
response to the incoming wave, re-radiating it backward.

---

## Polarization — Twist Direction

The traveling wave on the d=3 lattice has a T1u vector character — it
oscillates along a specific direction transverse to propagation. This
direction is the polarization.

```
Propagation along z:
  x-polarized: phi oscillates along x
  y-polarized: phi oscillates along y
  circular: phi rotates in the x-y plane
```

The lattice has cubic symmetry (Oh), which at long wavelengths looks
isotropic. So unpolarized light has no preferred direction — all
transverse orientations are equivalent.

In a crystal, the Oh symmetry of the lattice is BROKEN by the crystal
structure (a different, macroscopic lattice superimposed on the Planck
lattice). Different polarizations see different effective spring constants
along different crystal axes. This is birefringence — the material has
different refractive indices for different polarization directions.

---

## Diffraction and Interference

Waves on the lattice obey superposition (the Lagrangian is linear for
small amplitudes). Two waves meeting at a point ADD their displacements:

```
phi_total = phi_1 + phi_2
```

Constructive interference: phi_1 and phi_2 in phase, amplitude doubles.
Destructive interference: phi_1 and phi_2 out of phase, amplitude cancels.

Diffraction occurs when a wave passes through an opening comparable to
its wavelength. The opening acts as a new source — the springs at the
edge re-radiate in all directions. The interference pattern of these
edge-radiated waves creates the diffraction pattern.

All of this is standard wave mechanics. GWT adds nothing new to diffraction
and interference — they work the same on any spring network. The Planck
lattice spacing is so far below optical wavelengths (~10^30 times smaller)
that the lattice discreteness is invisible. Light sees a continuous medium.

---

## The Photoelectric Effect — Mode Coupling Threshold

A traveling wave hitting a bound configuration (kink + harmonics) can
excite or detach a harmonic mode. The threshold condition:

```
hbar × omega >= IE   (the ionization energy of the outermost mode)
```

Below threshold: the traveling wave cannot excite any mode transition.
No energy is absorbed. The material is transparent at that frequency.

Above threshold: the wave couples to a mode, the mode detaches as a
free breather. The excess energy becomes the breather's kinetic energy.

This is the photoelectric effect. It appears "quantized" because the
harmonic modes are discrete (standing waves in a well have specific
energies). The traveling wave itself is continuous — but the modes it
couples to are not.

Einstein got the Nobel Prize for this. In GWT: it's a traveling wave
resonantly coupling to a standing wave. No photon "particles" needed.

---

## Cherenkov Radiation — Faster Than Light in a Medium

When a charged kink moves through a material faster than the local wave
speed v = c/n, it creates a shock cone — a Mach cone of lattice waves.

```
cos(theta) = v_wave / v_particle = 1 / (n × beta)
```

The kink outruns the waves it generates in the medium. The waves pile up
into a cone. This is Cherenkov radiation — the optical equivalent of a
sonic boom. Same physics, same formula, same springs.

---

## Connection to the Framework

All of optics reduces to: **waves on springs with varying stiffness and density.**

| Phenomenon | Spring interpretation |
|-----------|---------------------|
| Refraction | Wave speed changes with local spring properties |
| Reflection | Impedance mismatch at boundaries |
| Dispersion | Resonant coupling to bound modes (harmonic frequencies) |
| Polarization | Transverse oscillation direction (T1u vector) |
| Diffraction | Edge re-radiation + superposition |
| Photoelectric | Traveling wave coupling to standing wave threshold |
| Birefringence | Crystal breaks Oh symmetry, different axes ≠ stiffness |
| Cherenkov | Mach cone from exceeding local wave speed |

No new constants. No new physics. Just the same Lagrangian in regions
where kinks modify the local spring properties.
