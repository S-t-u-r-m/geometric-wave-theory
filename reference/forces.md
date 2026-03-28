# Forces — Hooke's Law Decomposition

*Part of the [GWT Complete Reference](../gwt_complete_reference.md). See also: [Foundation](foundation.md), [Bonding](bonding.md), [Cosmology](cosmology.md).*

### The master equation
```
F_total = -k * x     (Hooke's law, Planck units)
```
Every force in nature is this spring force decomposed into spatial directions.

### The 1/3 – 2/3 split
On a d-dimensional isotropic lattice, a radial displacement between two nodes creates:
- **Longitudinal component**: 1/d = 1/3 of the force (points radially toward/away)
- **Transverse component**: (d-1)/d = 2/3 of the force (perpendicular to radial)

| Component | Fraction | Physical identity | How |
|-----------|----------|-------------------|-----|
| Longitudinal | 1/d = 1/3 | **Gravity** | Inward pull through torus center. Always attractive (convergence-independent). Falls as 1/r^2 from geometry. |
| Transverse | (d-1)/d = 2/3 | **Dark energy** | Restoring pressure perpendicular to radial. Net outward push at cosmic scales. |

**Why gravity is weak**: It's 1/3 of the total force, and the 1/3 mostly cancels at distance. The residual falls as 1/r^2 — not because of a special law, but because that's how solid angles work in 3D.

### 1D energy principle

Energy is fundamentally 1D — confined in one cosine well. All 3D structure
(torus, quarks, confinement) is the lattice's RESPONSE to this confined energy.

**Physical picture:** Energy is structureless compressible pressure (like a gas)
trapped in a perfectly geometric container (the d=3 cosine lattice). The energy
has no intrinsic structure — it just pushes. ALL structure comes from the
container geometry (Oh symmetry, T1u⊗T1u channels, 8 modes, 24 orientations).
The container shapes the structureless pressure into particles.

The gas pushes on all walls but can't escape — every direction has another well.
The cosine potential sets a maximum energy density of 2/pi^2 per site. Below
this: the gas sloshes (breather = electron). At this limit: the gas punches
through to the next cell (kink = proton, a topological phase transition).
The container geometry determines what the deformation looks like
(torus, A1g ground state by Perron-Frobenius).

What this means for particles:
- **Mass** = how much fluid is trapped = how much pressure it exerts
- **Kink (proton)** = fluid has pushed over the hill permanently (topological, static)
- **Breather (electron)** = fluid sloshes back and forth in the well (oscillatory)
- **Photon** = the pressure wave propagating through the lattice, no trapped fluid (massless)
- **Torus** = the fluid's equilibrium pressure pattern (A1g = uniform, minimum stress)
- **Quarks** = the fluid pushes on 3 walls (3 axes), each wall gets a fraction of pressure
- **Confinement** = can't remove fluid from one wall only — pressure fills the whole well
- **Generations** = fluid restricted to fewer walls → higher pressure per wall → heavier
- **E = mc²** = trapped fluid energy (ℏω) relates to lattice propagation speed (c)

Evidence:
- Breathers are quasi-1D on 32³ lattice (modes 1-7 match 1D to <0.12%)
- Bond f_pi: breathers have d²+1 modes (d² spatial + 1 temporal oscillation);
  kinks have only d modes (no temporal DOF because they're static topology)
- VP_self = -0.7589 is universal — depends only on 1D breather profile sech²

The particle spectrum factors as: mass = f(1D mode) × g(Oh lattice response).
This is already what GWT computes: m_pi = (zero-mode energy) × (A1g fraction).

**Why d=3:** A 1D trapped fluid on a lattice needs d=3 for stability.
In d≤2: Derrick's theorem — the fluid radiates away, no stable trapping.
In d≥4: coupling too weak (alpha = 1/24753), fluid barely interacts with lattice.
Only d=3 gives stable trapping with the right interaction strength.

### Force channel separation at atomic scale [TESTED, CONFIRMED]

The 1/d longitudinal (gravity) channel does NOT contribute to atomic bonding.
Tested: adding any fraction of the gravity channel to C_bond = pi/d² makes
the bond model WORSE (7.5% → 32.1% at full strength, monotonic increase).

The two channels operate at DIFFERENT scales:
- Transverse (d-1)/d: EM coupling alpha → bonds at eV scale
- Longitudinal 1/d: gravitational coupling alpha_G → bonds at Planck scale

Same Oh fractions, same T1u⊗T1u decomposition, but different force carriers.
The 4/27 = m_pi/m_p connection applies to NUCLEAR bonding (MeV), not atomic.
