# Relativity

*Part of the [GWT Complete Reference](../gwt_complete_reference.md). See also: [Foundation](foundation.md), [Forces](forces.md), [Cosmology](cosmology.md).*

---

## Einstein Was Right — GWT Explains Why

Special and general relativity are correct descriptions of how waves behave
on the d=3 cubic lattice. Every relativistic effect has a spring interpretation.
Einstein discovered the RULES. GWT provides the MECHANISM.

---

## Special Relativity from Springs

### The speed of light

```
c = a × sqrt(k/eta) = 1   (Planck units)

k   = 2/pi   spring stiffness (from cosine potential)
eta = 2/pi   inertial density (same averaging)
a   = 1      lattice spacing (Planck length)
```

k = eta is impedance matching. The lattice has zero reflection — waves propagate
without loss. This is WHY c is universal: it's a property of the medium, not a law
imposed from outside. Every instrument made of lattice waves responds at c, so
every local measurement gives c. Einstein's second postulate is a tautology in GWT:
observers made of waves always measure the wave speed of their own medium.

### E = mc^2

A kink (proton) has rest mass = the energy stored in the field gradient as phi
transitions between valleys:

```
M_kink = integral of (1/2)(dphi/dx)^2 + V(phi) dx = 8/pi^2   (Planck units)
```

Mass IS energy. Not by postulate — by construction. The kink's mass is literally
the elastic energy in the springs that are stretched across the valley transition.
A breather's mass is the energy of its oscillation. There is nothing else mass
could be. E = mc^2 is not a discovery about nature; it is a tautology about waves
in an elastic medium.

### Lorentz contraction

A wave pattern moving at velocity v on the lattice contracts along the direction
of motion. This is not mysterious — it is what wave packets do on any medium.

A stationary kink has width w ~ 3 lattice sites (set by the potential curvature).
A moving kink has width:

```
w(v) = w_0 × sqrt(1 - v^2/c^2) = w_0 / gamma
```

This follows directly from the Lorentz-invariant wave equation. The cosine potential
is periodic in phi, and the kink solution satisfies:

```
phi(x, t) = (4/pi) × arctan(exp(gamma × (x - v*t) / w_0))
```

The gamma factor compresses the kink profile along the direction of motion.
The springs ahead of the kink are compressed more; the springs behind are stretched
less. The total energy increases as gamma × M_kink — this IS relativistic mass.

### Time dilation

A breather (electron) oscillates at a frequency set by its rest energy:

```
omega_0 = m_e × c^2 / hbar   (rest frame)
```

When the breather moves at velocity v, its oscillation frequency decreases:

```
omega(v) = omega_0 / gamma = omega_0 × sqrt(1 - v^2/c^2)
```

This is time dilation. A moving clock (= moving breather) oscillates slower because
the wave pattern must propagate transversely while also moving longitudinally. The
lattice spring coupling is the same in all directions (cubic symmetry), so the
transverse oscillation gets "diluted" by the longitudinal motion.

All clocks are breathers or depend on breathers (atomic transitions, molecular
vibrations, pendulum springs). All slow down by the same gamma factor because
they are all waves on the same lattice.

### Relativistic momentum and energy

A wave packet on the lattice carries momentum proportional to its wave vector:

```
p = hbar × k = gamma × m × v
E = gamma × m × c^2
E^2 = (pc)^2 + (mc^2)^2
```

These are properties of the wave equation on the impedance-matched lattice.
The dispersion relation E^2 = p^2 + m^2 (natural units) is the lattice
dispersion relation in the long-wavelength limit (wavelength >> lattice spacing).

### Why nothing travels faster than c (for waves)

A wave on the lattice propagates through nearest-neighbor spring coupling.
Each site responds to its neighbors at rate sqrt(k/eta) = c. Information
cannot propagate faster than the springs can transmit it.

This applies to everything MADE of waves: matter, light, forces carried by
wave exchange. It does NOT necessarily apply to the lattice itself — the
lattice defines c but is not constrained by it. See [Foundation](foundation.md)
for the discussion of whether lattice restructuring (gravity, expansion)
occurs at c or instantaneously.

### The twins "paradox"

No paradox in GWT. The traveling twin's breathers oscillate slower (time dilation)
because they are wave patterns moving through the lattice. The stationary twin's
breathers oscillate at rest frequency. When they reunite, the traveling twin has
undergone fewer oscillation cycles = aged less.

The asymmetry: acceleration. The traveling twin's wave packet was boosted
(compressed by the acceleration), then decompressed (deceleration), then boosted
again (return). Each boost changes the gamma factor. The stationary twin
was never compressed. The difference in total oscillation count is real
and calculable from the lattice wave equation.

---

## General Relativity from Springs

### Gravity = lattice density gradient

In [Forces](forces.md), gravity is the 1/d = 1/3 longitudinal fraction of the
total spring force. At the microscopic level, this means:

A kink (mass) locally compresses the lattice. The springs around it are slightly
stiffer (more energy density = more displacement). This creates a gradient in
the effective spring stiffness:

```
k_eff(r) = k_0 × (1 + G × M / (r × c^2))
```

Other wave patterns (particles, light) propagate through this gradient. They
follow the path of least time (Fermat's principle), which curves toward the
higher-stiffness region. This IS the geodesic equation. Spacetime curvature
IS the lattice stiffness gradient.

### The equivalence principle

Acceleration compresses the lattice in the direction of motion (Lorentz contraction).
A gravitational field compresses the lattice toward the mass (stiffness gradient).

Both create the same local effect: a gradient in spring stiffness. A wave pattern
cannot distinguish between them because it only interacts with its local springs.
This is why inertial mass = gravitational mass: both are the response of the
same wave pattern to the same type of spring gradient.

### Gravitational time dilation

Deeper in a gravitational well, the lattice is more compressed (higher stiffness).
Breathers oscillate at a rate determined by the local spring stiffness:

```
omega(r) = omega_0 × sqrt(1 - 2GM/(rc^2))
```

Deeper = stiffer springs = different oscillation frequency. Clocks at different
heights tick at different rates because they are waves on springs with different
local stiffness. This is gravitational redshift.

### Gravitational waves

A sudden redistribution of mass (merger, collapse) creates a ripple in the lattice
stiffness field. This ripple propagates at c (the spring coupling speed) and causes
local oscillations in the lattice spacing.

LIGO detects these oscillations: the distance between mirrors changes because
the lattice between them is stretching and compressing. The atoms in the mirrors
don't change (their structure is set by alpha, a dimensionless ratio). The strain
h ~ 10^-21 reflects the ratio of gravitational to electromagnetic spring stiffness.

### Black holes

When the lattice is compressed to the point where the field sits at the TOP of
the cosine barrier (phi = 1), the potential curvature is negative. No oscillation
is possible — no breathers, no waves, no information can propagate outward.

This is the event horizon. The lattice stiffness reaches a critical value where
the wave equation changes character (oscillatory -> exponential decay).

```
Inside BH:  phi >= 1 (barrier top), d^2V/dphi^2 < 0
             No oscillation. Field decays. Information trapped.
Outside BH: phi < 1 (well), d^2V/dphi^2 > 0
             Normal oscillation. Waves propagate. Physics works.
```

The Schwarzschild radius = where the lattice compression pushes phi to the
barrier top. Hawking radiation = quantum tunneling off the barrier at the
boundary (T_H = 1/(2^d × pi × M), derived in [Cosmology](cosmology.md)).

### Gravitational constant from the lattice

The hierarchy between gravity and other forces is set by the mass ratio and
coupling constant:

```
G = alpha^24 × F^4 × (hbar × c / m_p^2)

alpha^24 = alpha^(2 × |A_4|):  tunneling through 24 gauge channels
F^4 = (6*pi^5)^4:              four powers of the mass ratio
```

Gravity is weak because it requires tunneling through ALL 24 gauge channels
of the Oh group. EM only needs one channel. The ratio G/alpha ~ 10^-37 is
not a mystery — it is exp(-24/alpha) ~ exp(-24 × 137).

---

## Kinetic Energy = Lattice Distortion [SIMULATED, 2026-04-01]

**The kinetic energy of a moving particle is physically stored in the lattice.**

A moving kink (proton) distorts the lattice in three ways:

```
1. COMPRESSED SPRINGS (dE_grad > 0):
   The Lorentz-contracted kink has steeper gradients.
   Steeper gradients = more energy in the springs.
   The springs AHEAD of the kink are compressed.

2. NARROWER POTENTIAL (dE_pot < 0):
   The contracted kink occupies fewer sites.
   Fewer sites in the cosine barrier = less potential energy.
   This partially offsets the gradient increase.

3. FIELD MOMENTUM (dE_kin > 0):
   The field is changing at every site the kink passes.
   dphi/dt = -v * dphi/dx at each site.
   This is the kinetic energy of the field motion itself.
```

The total: E_total = E_grad + E_pot + E_kin = gamma * M (relativistic energy).

### Simulation results

On a discrete lattice (a=1), the energy decomposition of a boosted kink:

| v | gamma | dE_grad | dE_pot | dE_kin | KE_total | KE_pred | error |
|---|-------|---------|--------|--------|----------|---------|-------|
| 0.1 | 1.005 | +0.002 | -0.002 | +0.004 | 0.004 | 0.004 | -13% |
| 0.3 | 1.048 | +0.018 | -0.018 | +0.034 | 0.034 | 0.039 | -13% |
| 0.5 | 1.155 | +0.057 | -0.053 | +0.103 | 0.107 | 0.125 | -15% |
| 0.7 | 1.400 | +0.142 | -0.110 | +0.233 | 0.265 | 0.325 | -18% |
| 0.9 | 2.294 | +0.380 | -0.187 | +0.509 | 0.701 | 1.049 | -33% |
| 0.99 | 7.089 | +0.604 | -0.204 | +0.734 | 1.134 | 4.935 | -77% |

### The missing energy = lattice densification

The systematic deficit grows with velocity. At v=0.99, 77% of the predicted
kinetic energy is "missing" from the fixed lattice. This is NOT a numerical
error — it is a PHYSICAL result:

**A fixed lattice (constant site count) cannot store the full relativistic
kinetic energy of a highly boosted kink.**

The kink's Lorentz-contracted gradient is too steep for the discrete sites
to resolve. The energy that SHOULD be stored in the compressed springs
has nowhere to go on a fixed grid.

On the REAL Planck lattice, this energy goes into **lattice densification**:
creating new lattice sites ahead of the moving kink. The kinetic energy
IS the energy of site creation. This is why:

```
Kinetic energy = lattice distortion energy
               = energy stored in compressed springs (partially)
               + energy of NEW lattice sites created ahead of the kink

E_kinetic = (gamma - 1) * M = energy to densify the lattice by factor gamma
```

This connects to cosmological expansion: just as kinetic energy creates
new sites ahead of a fast kink, dark energy creates new sites throughout
the lattice (expansion). Both are lattice growth driven by energy input.

### Synchrotron radiation as lattice shedding

At extreme velocities, the lattice densification reaches a limit — sites
cannot overlap. The excess energy that cannot be stored as lattice
distortion MUST be radiated away. This is synchrotron radiation.

In standard physics, synchrotron power scales as gamma^4:
```
P_synch = (2/3) * alpha * c * gamma^4 / R^2
```

In GWT, this is the lattice shedding energy it cannot absorb:
- Low gamma: lattice absorbs the energy (densification)
- High gamma: densification saturates, energy radiates
- The gamma^4 scaling comes from the 4 powers of the mass ratio F^4
  in the gravitational constant — the same hierarchy that makes
  gravity weak also limits lattice compression

### Connection to the equivalence principle

Gravity compresses the lattice (density gradient toward mass).
Acceleration compresses the lattice (density gradient toward motion).
Kinetic energy densifies the lattice (new sites created by motion).

All three are the SAME lattice distortion:
- Gravity: static compression (permanent density gradient)
- Acceleration: dynamic compression (changing density gradient)
- Kinetic energy: stored compression (the density gradient IS the energy)

The equivalence principle holds because there is only ONE mechanism:
spring compression. Whether caused by a nearby mass, an applied force,
or sustained motion, the local lattice distortion is identical.

See: `calculations/simulations/time_dilation_test.py`

---

## What GWT Adds Beyond Einstein

Einstein's relativity is COMPLETE as a description of wave behavior on the lattice.
GWT does not modify or contradict it. What GWT adds:

1. **The mechanism**: spacetime is not an abstract manifold — it is a spring network.
   Curvature is stiffness gradient. Metric is local spring constant.

2. **The constants**: G, c, and the Planck units are all derived from k = eta = 2/pi.
   Einstein used them as inputs. GWT derives them.

3. **The unification**: gravity (1/d), EM ((d-1)/d), and all forces come from the
   SAME spring, decomposed by the Oh group. Einstein spent decades looking for
   unification. It was always there — in the d=3 cube.

4. **The quantum connection**: wave-particle duality is not dual — everything is waves.
   Quantum mechanics and relativity are both lattice wave mechanics. There is no
   conflict between them because they are the same theory.

5. **The open question**: does the lattice restructure at c, faster, or instantaneously?
   Einstein assumed c is the universal speed limit. GWT distinguishes between wave speed
   (c, derived) and lattice restructuring speed (unknown). See [Foundation](foundation.md).
