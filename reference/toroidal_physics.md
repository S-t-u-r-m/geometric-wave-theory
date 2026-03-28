# Toroidal Breather Physics

*Part of the [GWT Complete Reference](../gwt_complete_reference.md). See also: [Foundation](foundation.md), [Bonding](bonding.md), [Nuclear](nuclear.md).*

---

### What particles ARE
Fermions are **toroidal circulations** (smoke rings) on the lattice. Not radial pulsing (which would be spin-0). The lattice is impedance-matched (k = eta), so these vortex rings never dissipate.

### Why the torus — THEOREM (Perron-Frobenius) [DERIVED]

The A1g torus wrapping is the UNIQUE ground state of the kink Hamiltonian
on the d-cube. This is a theorem, not a conjecture.

**Proof:**

1. The Lagrangian L = (1/2)(dphi)^2 + (1/pi^2)(1-cos(pi*phi)) is Oh-symmetric.
   => Hamiltonian H commutes with all 48 Oh operations.
   => Eigenstates of H labeled by Oh irreps (A1g, T1u, Eg, T2g, ...).

2. The NN coupling between kinks on adjacent faces is ATTRACTIVE
   (adjacent kinks share an edge where the field gradient is reduced).
   => Off-diagonal elements of H are non-positive.
   => Transfer matrix T = exp(-beta*H) has ALL POSITIVE entries.

3. Perron-Frobenius theorem (for matrices with all positive entries):
   a) The largest eigenvalue of T is NON-DEGENERATE.
   b) The corresponding eigenvector has ALL POSITIVE components.
   c) This eigenvector = the ground state of H.

4. An all-positive vector on the 2^d cube vertices is TOTALLY SYMMETRIC
   under Oh = A1g irrep. This is the uniform (torus) wrapping.

**Therefore the A1g torus IS the ground state. QED.**

**Explicit spectrum on the d-cube (H = -J * adjacency):**
```
  k=0: E = -dJ = -3J,  mult = C(d,0) = 1,  irrep = A1g   [GROUND STATE]
  k=1: E = -(d-2)J = -1J, mult = C(d,1) = 3, irrep = T1u  [1st excited]
  k=2: E = +(d-2)J = +1J, mult = C(d,2) = 3, irrep = T2g  [2nd excited]
  k=3: E = +dJ = +3J,  mult = C(d,3) = 1,  irrep = A2u   [maximum]
```

Energy gap: E(T1u) - E(A1g) = 2J. Gap/ground = 2/d = 2/3.
The torus is robustly stable — the gap equals 2/3 of the ground state energy.

**Physical meaning of each state:**
  A1g (E = -3J): all 12 edges aligned = smooth torus = PROTON
  T1u (E = -1J): one axis flipped = vector excitation (3 states)
  T2g (E = +1J): two axes flipped = tensor excitation (3 states)
  A2u (E = +3J): all axes flipped = fully anti-aligned (1 state)

**Energy comparison of topological defects:**
```
Sphere (hedgehog, j_0):  E = 4*pi*a^2 * M_kink = 12.57 * a^2 * M
Torus (vortex ring):     E = 2*pi*a^2 * M_kink =  6.28 * a^2 * M
Knot:                    E > torus (longer path)
Double torus (genus 2):  E > 2*torus (more surface)
```
The torus has HALF the energy of the sphere — factor (d-1) = 2 lower.
The simple torus (genus 1) is the unique minimum among all closed topological defects.

**Lattice stabilization (Derrick's theorem defeated):**
- In continuum d>=2: Derrick's theorem forbids stable static scalar solitons
- On discrete lattice: theorem doesn't apply because:
  1. Minimum tube radius = lattice spacing a (can't shrink below one unit)
  2. Kink winding number is INTEGER (can't partially unwind)
  3. Energy barrier to unwind = M_kink * 2*pi*a (topological protection)

**Why 3 sub-components (only in d=3):**
- A torus always has exactly 3 independent motions (toroidal, poloidal, twist)
- In d=2: torus degenerates to circle (1 motion) — no quarks
- In d=4: 3 torus motions don't fill 4 lattice axes — mismatch
- Only in d=3: 3 torus motions = 3 lattice axes = perfect match
- This gives: 3 quarks, 3 colors, charge fractions 1/d and (d-1)/d

### Three torus motions = three quantum numbers
A torus in d=3 has exactly 3 independent motions:

| Motion | Physical quantity | Values |
|--------|-------------------|--------|
| Toroidal flow (around big ring) | Electric charge | quantized |
| Poloidal flow (through hole) | Color charge | 3 |
| Helical twist (around tube) | Spin | 0, 1/2, 1 |

### Quantum number mapping
| Torus property | Quantum number | Values |
|----------------|----------------|--------|
| Flow direction (which end converges) | Matter/antimatter | 2 |
| Internal twist direction (CW/CCW) | Chirality (L/R) | 2 |
| Internal twist count | Generation | 3 |
| External rotation | Spin | 2 |
| Toroidal winding rate | Electric charge | quantized |
| Poloidal winding rate | Color charge | 3 |
| Ring axis orientation | Orientation | 3 |

Total fermion states: 2 × 2 × 3 = 12 per type, × 2 for matter/antimatter = **24**.

### Quarks as sub-circulations
A toroidal vortex in 3D decomposes into 3 coupled flow loops (one per axis). Each loop is a quark.
```
1/d = 1/3 of flow along one axis  → charge 1/3 (down quark)
(d-1)/d = 2/3 across other two    → charge 2/3 (up quark)
```
Same 1/3 vs 2/3 split as gravity vs dark energy. Confinement is topological — can't separate one loop without destroying the torus.

### Stability
- Circulation on a discrete lattice is quantized: +1 or -1, no in-between
- Kelvin's circulation theorem: in a lossless medium (k=eta), circulation is conserved
- Winding number is an integer — can't change continuously. Flip costs 2mc^2.

### Annihilation
- Opposite flows ATTRACT (not repel), then cancel on contact
- Topology unravels, all energy radiates as photons: E = 2mc^2
- Like two opposite whirlpools merging — they don't bounce, they destroy

### Antimatter gravity prediction (confirmed)
- Gravity = inward pull through torus center (1/d longitudinal component)
- This pull is ALWAYS inward regardless of flow direction
- Both matter and antimatter gravitate identically
- ALPHA-g experiment at CERN (2023): antihydrogen falls down. **Confirmed.**

### Toroidal coupling modes (bonding decomposition)
Two toroidal breathers interact through 3 modes, each with a geometric weight:

| Mode | Weight | Bond formula parameter | Physical role |
|------|--------|------------------------|---------------|
| Toroidal | ~83% (pi^2/(pi^2+2)) | C = pi/d, f_pi = 9/10 | Bond energy (sigma/pi overlap) |
| Poloidal | ~17% (2/(pi^2+2)) | f_anti = 6/5, alpha = 7/10 | Equilibrium distance, bond angle |
| Twist | ~1% (alpha*kappa/k) | Pauli exclusion | Spin pairing |

**Dipole coupling model:**
- Toroidal-toroidal: E_TT = -m_T^2/r^3 (attractive, like parallel currents)
- Poloidal-poloidal: E_PP = +2*m_P^2/r^3 (repulsive, head-to-head dipoles)
- Cross-term: E_TP = 0 (perpendicular dipoles on-axis)
- R/a = pi gives m_T/m_P = pi, so toroidal dominates: pi^2/(pi^2+2) = 83.2%

**Connection to bond formula:**
- Sigma bond = toroidal overlap (main sin(phase) term)
- Pi bond = transverse toroidal overlap (reduced by f_pi = 9/10)
- Antibonding excess = poloidal repulsion: f_anti - 1 = 1/(2d-1) = 20%
- Bond angle = poloidal return flow symmetry: cos(theta) = -1/(d+1)
- Twist correction ~ 0.05 eV/bond (~1% of bond energy), controls spin pairing

The V8 formula already captures all three modes through its 6 parameters.
Residual ~2% errors are higher-multipole (3D geometry) corrections, not missing modes.

See: calculations/bonding/toroidal_coupling_modes.py for full calculation.
