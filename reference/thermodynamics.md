# Thermodynamics

*Part of the [GWT Complete Reference](../gwt_complete_reference.md). See also: [Foundation](foundation.md), [Lattice & Symmetry](lattice_and_symmetry.md), [Standing Waves](standing_waves.md).*

---

## Temperature Is Lattice Vibration

Temperature is not a fundamental quantity. It is the average kinetic energy
of the lattice oscillations — the random jiggling of the springs.

Every lattice site has zero-point energy (the springs are always vibrating).
Above zero-point, additional energy distributes across the available modes.
Temperature measures the AVERAGE energy per mode:

```
E_avg = (1/2) × k_B × T   per degree of freedom

k_B = Boltzmann constant = the conversion between energy and mode frequency
T   = temperature = the label for average mode energy
```

In GWT, each lattice site has d=3 spatial degrees of freedom. The total
thermal energy per site:

```
E_thermal = d × (1/2) × k_B × T = (3/2) × k_B × T
```

This IS the equipartition theorem. It follows directly from the lattice
having d=3 spring directions, each storing energy independently.

---

## Entropy Is Mode Counting

Entropy counts the number of ways energy can be distributed across modes.

A lattice with N sites and d directions has N × d oscillation modes.
The modes are labeled by wavevector k (which direction the wave propagates)
and frequency omega (how fast it oscillates). Each mode can hold energy
in integer multiples of hbar × omega.

```
S = k_B × ln(W)

W = number of ways to distribute total energy across all modes
```

For the d=3 lattice, the modes decompose under Oh symmetry:

```
T1u: d = 3 acoustic modes (long-wavelength sound)
A1g: 1 breathing mode
Eg:  2 shear modes
T2g: 3 additional shear modes
...
```

At low temperature, only the acoustic modes (T1u) are excited — these
have the lowest frequency. As temperature increases, higher-frequency
modes (Eg, T2g, etc.) become accessible. Each newly accessible mode
increases the entropy by k_B × ln(2) (one new binary choice).

---

## Heat Capacity from Mode Activation

Heat capacity = how much energy is needed to raise the temperature by 1 K.
This depends on how many modes can absorb energy at the current temperature.

### Low temperature: Debye T^3 law

At low T, only long-wavelength acoustic modes are excited. The number of
accessible modes scales as the volume of a sphere in k-space:

```
N_modes ~ (T / T_Debye)^d = (T / T_Debye)^3

C_V = (12/5) × pi^4 × N × k_B × (T / T_Debye)^3
```

The d=3 exponent in T^3 comes directly from the dimensionality of the
lattice. In d=2 it would be T^2 (and IS T^2 for 2D materials like graphene).

### High temperature: Dulong-Petit limit

At high T, all modes are excited. Each mode contributes k_B to the heat
capacity:

```
C_V = d × N × k_B = 3 × N × k_B
```

d modes per site, each contributing k_B. The factor d=3 gives the
Dulong-Petit value of 3R = 24.9 J/(mol·K). This is exact for the d=3
lattice.

### The Debye temperature

The Debye temperature T_D marks the crossover between T^3 (few modes)
and 3Nk_B (all modes). It is set by the maximum frequency of the lattice:

```
k_B × T_D = hbar × omega_max

omega_max = the highest frequency mode = the Brillouin zone boundary
```

For the Planck lattice: omega_max = c/a = the Planck frequency. So
T_D(Planck) = the Planck temperature ~ 10^32 K. In practice, no
material reaches this — atomic lattices have much lower T_D (100-1000 K)
because their spacing is ~10^25 times larger than the Planck length.

---

## The Boltzmann Distribution from Wave Statistics

The probability of a mode having energy E at temperature T:

```
P(E) = exp(-E / k_B T) / Z

Z = sum over all states of exp(-E_n / k_B T)   (partition function)
```

In GWT, this follows from the wave equation on the lattice. Each mode
is a harmonic oscillator with energy levels E_n = (n + 1/2) × hbar × omega.
The thermal population of each level is set by the Boltzmann factor,
which arises from maximizing the number of ways to distribute total
energy across modes (entropy maximization).

The Boltzmann distribution is NOT a postulate. It is the MOST PROBABLE
distribution for waves on a lattice with fixed total energy. Any other
distribution has fewer microstates and is therefore less likely.

---

## The Laws of Thermodynamics

### Zeroth law (thermal equilibrium)

Two lattice regions in contact exchange energy through their shared springs.
Energy flows from higher-amplitude oscillations to lower-amplitude until
the average amplitude (temperature) equalizes. Equilibrium = equal average
mode energy = same T.

### First law (energy conservation)

The total spring energy (kinetic + potential) is conserved. Energy can
transfer between regions or convert between modes, but the total is fixed.
This is a direct consequence of the Lagrangian being time-independent.

```
dU = delta_Q - delta_W

U = total lattice energy
Q = energy transferred via mode coupling (heat)
W = energy transferred via coherent displacement (work)
```

Heat = random mode-to-mode energy transfer. Work = organized, coherent
displacement of the lattice. Same energy, different organization.

### Second law (entropy increase)

Energy spreads across modes because there are MORE ways to distribute
energy evenly than unevenly. A single high-energy mode can distribute
its energy across d × N partner modes through spring coupling. The
reverse (all partners spontaneously sending energy to one mode) requires
exact phase alignment of d × N modes — probability ~ (1/N)^N.

Entropy increases because the lattice has d=3 directions at each site,
giving exponentially more spread-out configurations than concentrated ones.
This is not a law — it is STATISTICS of waves on a spring network.

### Third law (absolute zero)

At T=0, each mode is in its ground state (zero-point energy only).
There is exactly ONE configuration: all modes at minimum energy.
S = k_B × ln(1) = 0.

The zero-point energy remains: (1/2) × hbar × omega per mode. This
is the lattice's irreducible jiggle — the springs are never perfectly
still. Zero temperature means no THERMAL energy above zero-point,
not zero energy total.

---

## Phase Transitions from Mode Structure

Phase transitions occur when the mode structure changes qualitatively:

### Solid -> Liquid

In a solid, the kinks (atoms) are locked into a periodic arrangement.
Their modes are phonons — coherent lattice vibrations. When thermal
energy exceeds the binding energy between kinks (D_e), the kinks
delocalize. The mode structure shifts from discrete phonon bands
to a continuum. This is melting.

```
T_melt ~ D_e / (d × k_B)

D_e = bond energy (from V8 formula)
d = 3 degrees of freedom
```

### Liquid -> Gas

In a liquid, kinks are mobile but still interacting. When thermal
energy exceeds the interaction range, kinks separate completely.
The mode structure shifts from coupled oscillations to free propagation.
This is boiling.

### Magnetic transitions

At the Curie temperature, unpaired d-orbital modes (which create
ferromagnetism — see [Magnetism](magnetism.md)) lose their coherent
alignment. Thermal energy overwhelms the exchange coupling between
neighboring kink configurations. The T1g (magnetic) modes randomize.

---

## Connection to the Framework

| Concept | GWT interpretation |
|---------|-------------------|
| Temperature | Average spring vibration amplitude per mode |
| Entropy | Number of accessible Oh mode configurations |
| Heat capacity | Number of modes that can absorb energy at given T |
| Heat | Random mode-to-mode energy transfer |
| Work | Coherent lattice displacement |
| Phase transition | Qualitative change in available mode structure |
| Boltzmann distribution | Most probable wave energy distribution (statistics) |
| T^3 law | d=3 phase space volume for acoustic modes |
| Dulong-Petit (3Nk_B) | d=3 modes per lattice site |
| Absolute zero | All modes at zero-point, S = 0 |

No new constants. Temperature, entropy, and all thermodynamic laws
are emergent properties of waves on the d=3 spring network.
