# Radioactive Decay

*Part of the [GWT Complete Reference](../gwt_complete_reference.md). See also: [Foundation](foundation.md), [Nuclear](nuclear.md), [Toroidal Physics](toroidal_physics.md).*

---

## Decay Is Tunneling on the Cosine Potential

All radioactive decay in GWT reduces to one mechanism: a field configuration
tunneling through or restructuring on the cosine potential V = (1/pi^2)(1-cos(pi*phi)).

The potential has barriers of height 2/pi^2 between valleys. A configuration
trapped behind a barrier can tunnel through it (WKB), or the barrier can be
lowered by external conditions (other kinks nearby).

---

## Alpha Decay — Kink Cluster Tunneling

An alpha particle is a tightly bound cluster of 4 kink windings (2 protons +
2 neutrons in standard language, or winding number 2 + 2 uncharged windings
in GWT). This cluster sits inside a larger kink topology (the parent nucleus).

The parent nucleus has a potential well that confines the alpha cluster.
Beyond the well, the cluster would be free (lower energy). But between the
well and freedom, there is a Coulomb-like barrier — the combined electrostatic
repulsion from the remaining Z-2 windings.

### The tunneling rate

The WKB tunneling probability through a barrier:

```
T^2 = exp(-2 * integral sqrt(2m(V(r) - E)) dr)
```

For the cosine potential barrier, the base tunneling amplitude is:

```
T_0^2 = exp(-16/pi^2) = 0.1977
```

This is the SAME tunneling amplitude that gives alpha (the fine structure
constant). The barrier shape is universal — it's the cosine potential.

For alpha decay, the effective barrier includes the Coulomb repulsion
from the remaining nucleus:

```
V_barrier(r) = (Z-2) × 2 × alpha / r    (for r > R_nucleus)
```

The tunneling integral gives the Geiger-Nuttall law:

```
log(lambda) = a - b × (Z-2) / sqrt(E_alpha)

lambda = decay rate
E_alpha = kinetic energy of emitted alpha
a, b = constants from the barrier shape
```

In GWT: a and b come from the cosine potential + Coulomb tail. The
Geiger-Nuttall law IS the WKB integral on this combined potential.

### Why alpha clusters

Why does the nucleus emit alpha particles (4 windings) instead of
single kinks? Because 4 windings form an A1g-symmetric configuration
on the d=3 cube:

```
Alpha = 2 kink windings + 2 uncharged windings
      = 4 windings filling 4 of 2d=6 face directions
      = the most symmetric sub-configuration possible
```

The binding energy per winding is maximized for the alpha cluster
(B/A = 7.07 MeV, higher than any other light nucleus). This makes
it the lowest-energy "fragment" that can tunnel as a unit.

---

## Beta Decay — Torus Mode Conversion

Beta decay is NOT tunneling through a barrier. It is a mode conversion
on the torus — one type of winding converts to another.

### Beta-minus: neutron -> proton + electron + antineutrino

An uncharged kink winding converts to a charged winding. In torus language:
the poloidal phase goes from 0 (uncharged) to 1 (charged). This releases:
- A breather (electron) — the field oscillation created by the phase change
- A neutrino — the mode-basis excitation (see [Mixing & Neutrinos](mixing_neutrinos.md))

The energy released = n-p mass difference = m_e × 8/3 × (1-7*alpha) = 1.293 MeV.

### Beta-plus: proton -> neutron + positron + neutrino

The reverse: charged winding converts to uncharged. Requires energy input
(the proton is lighter than the neutron in free space), so it only occurs
inside nuclei where the binding energy difference provides the deficit.

### Weak interaction = torus rotation

The weak force mediates beta decay because it IS the rotational mode of
the torus (SU(2), see [Foundation](foundation.md)). Converting a charged
winding to an uncharged one requires rotating the torus — this is a W boson
exchange. The W mass (80.4 GeV) sets the rate: heavy mediator = slow process.

The beta decay rate scales as:

```
lambda ~ G_F^2 × E^5 / (30 * pi^3)

G_F = Fermi constant ~ alpha / (sin^2(theta_W) × M_W^2)
E = available energy (Q-value)
```

The E^5 dependence (Sargent's rule) comes from the 5-dimensional phase
space: 3 momentum directions + 2 spin states for the emitted leptons.
In GWT: 2d - 1 = 5, the same (2d-1) that appears in F_RAD and the
first p-shell quantum defect.

---

## Gamma Decay — Harmonic Relaxation

Gamma decay is the simplest: an excited nuclear configuration relaxes
to a lower energy state and emits the excess energy as a traveling wave.

This is IDENTICAL to atomic emission (see [Standing Waves](standing_waves.md))
— a harmonic mode drops to a lower level, and the energy difference
propagates away as a photon. The only difference is the energy scale:
nuclear modes are MeV (kink excitations), atomic modes are eV (breather
excitations).

```
Atomic emission:  breather mode transition,  E ~ eV
Gamma emission:   kink mode transition,      E ~ MeV
Same physics, different energy scale.
```

---

## Fission — Kink Topology Splitting

When a heavy nucleus splits, the single kink topology (high winding number)
divides into two separate topologies (lower winding numbers). The total
winding number is conserved (charge conservation).

The barrier to fission is the surface energy: the configuration must
deform through a saddle point where the surface area (spring energy) is
maximized before it splits into two separate pieces.

```
Barrier height ~ surface tension × deformation area
             ~ alpha × Z^2 / A^(1/3)   (fissility parameter)
```

Heavy nuclei (Z^2/A > 49) have barriers low enough for spontaneous fission.
The fissility threshold comes from the competition between:
- Coulomb repulsion (Z^2, tries to push apart)
- Surface tension (A^(2/3), tries to hold together)

The energy released per fission (~200 MeV) is the difference in binding
energy between the parent and the two daughter configurations. This comes
from the same Oh tensor product physics that gives B/A = 8.5 MeV in
[Nuclear](nuclear.md).

---

## Fusion — Kink Topology Merging

Fusion is the reverse: two separate kink topologies merge into one.
The barrier is the Coulomb repulsion between two charged configurations:

```
V_barrier = Z_1 × Z_2 × alpha / R_contact
```

For two protons: V ~ alpha / r_p ~ 0.6 MeV. The thermal energy needed
to overcome this (~10 keV in the sun) is much lower because quantum
tunneling provides the remaining factor — the same WKB amplitude T_0^2
from the cosine potential.

Fusion releases energy because light nuclei have lower B/A than the
iron peak (Z=26, where B/A is maximized). Merging two light configurations
into one heavier one releases the binding energy difference.

### Why iron is the endpoint

The B/A curve peaks at iron (A~56) because:
- Below iron: adding windings increases the A1g coupling (more binding)
- Above iron: the Coulomb repulsion (Z^2 term) overwhelms the A1g gain

The crossover happens at A ~ 2d × d! / (d-1) = 56 for d=3. Iron IS the
geometric sweet spot of the d=3 cube.

---

## Half-Lives — The Exponential Clock

All decay processes are probabilistic because tunneling is probabilistic.
The decay rate lambda is the tunneling probability per unit time:

```
N(t) = N_0 × exp(-lambda × t)
t_1/2 = ln(2) / lambda
```

The exponential decay law follows from the wave equation: at each oscillation
cycle, the field has probability T^2 of tunneling. Over many cycles, this
gives exponential survival.

Half-lives span 50+ orders of magnitude (10^-22 seconds for resonances
to 10^30 years for bismuth-209). This enormous range comes from the
EXPONENTIAL sensitivity of the tunneling probability to the barrier height
and width. Small changes in Z, A, or available energy produce huge changes
in T^2.

---

## Connection to the Framework

| Decay type | GWT mechanism | Key formula |
|-----------|--------------|-------------|
| Alpha | Kink cluster tunneling through Coulomb barrier | T^2 = exp(-16/pi^2) × Coulomb integral |
| Beta | Torus mode conversion (SU(2) rotation) | lambda ~ G_F^2 × E^5 |
| Gamma | Harmonic mode relaxation (same as atomic emission) | E = hbar × omega |
| Fission | Kink topology splitting at saddle point | Z^2/A fissility |
| Fusion | Kink topology merging through Coulomb barrier | Tunneling + thermal |

All five processes use the same potential, the same tunneling amplitude,
and the same topology. The cosine barrier sets the fundamental rate.
Everything else is geometry.
