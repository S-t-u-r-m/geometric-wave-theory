# Mixing Matrices & Neutrino Masses

*Part of the [GWT Complete Reference](../gwt_complete_reference.md). See also: [Mass Ratios](mass_ratios.md), [Nuclear](nuclear.md).*

## 9. MIXING MATRICES

### CKM Matrix (zero free parameters)

**Why sqrt (1/(d-1) = 1/2 power):**
Quarks are confined inside the proton, whose surface is (d-1) = 2 dimensional.
Mixing amplitudes go as mass^(1/(d-1)) = mass^(1/2) = sqrt.
Compare PMNS: leptons span the d=3 bulk, so they use 1/d = 1/3 (cube root).

**Derivation chain — three angles from quark mass ratios:**

```
N_top = d * 2^d + 1 = 25    [24 breathers + 1 kink = topological mode count]
N_br  = floor(2^d * pi - 1) = 24    [breather count from Lagrangian]
N_top/N_br = 25/24           [spectral completeness factor]
```

**Angle 1 — Cabibbo (both sectors in quadrature):**
```
sin^2(theta_12) = m_d/m_s + m_u/m_c
```
Both up and down sectors contribute on ORTHOGONAL planes of the 2D proton surface,
so their amplitudes add in quadrature (sin^2 = sum). No topological correction needed
because both sectors together span the full local flavor space.

**Angle 2 — nearest-neighbor (gen 2→3), partial topological correction:**
```
sin^2(theta_23) = (m_u/m_c) × (N_top/N_br)^(2/d)
```
Single sector (up only). Nearest-neighbor transition samples 1/d of the topological
space → exponent 2/d = 2/3 on (N_top/N_br) in sin^2.

**Angle 3 — full-range (gen 1→3), full topological correction:**
```
sin^2(theta_13) = (m_u/m_t) × (N_top/N_br)^2
```
Single sector (up only). The 1→3 transition spans ALL three generations, sampling the
complete topological space → full exponent 2 on (N_top/N_br) in sin^2.

**Exponent pattern in sin^2: 0, 2/d, 2  (= 0, 2/3, 2 for d=3)**
- theta_12: both sectors cancel the topological correction → 0
- theta_23: samples 1/d of the topology → 2/d
- theta_13: samples all of the topology → 2

**CP phase — antibonding geometry on proton surface:**
```
cos(delta_CKM) = 1/d + 2/(d+1)! = (2d-1)/(4d) = 5/12
```
Equivalently: cos(delta) = 1/(2*f_anti) where f_anti = 2d/(2d-1) = 6/5.
One axis (1/d = 1/3) plus one gauge gate (2/(d+1)! = 2/24 = 1/12).

**Validation of N_top/N_br corrections:**
- Without N_top/N_br on theta_13: V_ub off by -4.0%, Jarlskog off by -4.5%
- With (N_top/N_br)^2: V_ub error → -0.03%, Jarlskog → -0.5%
- Same factor with exponent 2/d fixes V_cb: -1.45% → -0.10%

| Element | Predicted | Observed | Error | sigma |
|---------|-----------|----------|-------|-------|
| V_us | 0.2242 | 0.2250 | -0.35% | 1.2 |
| V_cb | 0.04173 | 0.04182 | -0.10% | 0.1 |
| V_ub | 0.00369 | 0.00369 | -0.03% | 0.0 |
| delta_CKM | 65.38° | 65.5° | -0.2% | 0.0 |
| Jarlskog J | 3.06×10^-5 | 3.08×10^-5 | -0.5% | — |

All 9 matrix elements within 1.4 sigma. Mean error < 0.3%.

### PMNS Matrix (zero free parameters)

**Construction: single geometric rotation applied to tribimaximal mixing:**
```
U_PMNS = R(theta_corr, axis) × U_TBM

theta_corr = arcsin((m_e/m_mu)^(1/d))     [cube root — leptons span 3D bulk]
axis = (-1, sqrt(d), -(m_tau/m_p)^(1/d))  [normalized]
```
Leptons span 3D bulk → use 1/d = 1/3 power (cube root). Quarks on 2D surface → sqrt.

**Axis components (all from d=3 geometry):**
- **-1**: electron direction in TBM equilateral flavor triangle (tetrahedron vertex)
- **sqrt(3) = sqrt(d)**: muon vertex coordinate in the TBM degenerate subspace (lattice dimensionality)
- **-(m_tau/m_p)^(1/3)**: tau wrapping factor — tau sits inside the proton radius, so the proton "wraps around" it. The 1/d = 1/3 power is forced by 3D bulk overlap integrals. max(1, sigma_p/sigma_tau) applies only to tau (electron and muon are outside the proton).

**Rodrigues rotation formula (explicit construction, no external dependencies):**
```
K = skew-symmetric matrix from normalized axis vector
R = I + sin(theta_corr)*K + (1 - cos(theta_corr))*(K × K)
U_PMNS = R × U_TBM
```
This is the standard axis-angle rotation applied to tribimaximal mixing.
Every component is derived from d=3 geometry and GWT lepton masses.

**CP phase (lattice geometry):**
```
delta_PMNS = arccos(-1/d) = arccos(-1/3) = 109.47°
```
This is the tetrahedral dihedral angle — the angle between faces of a regular tetrahedron,
which is the fundamental angular unit of d=3 geometry.

| Angle | Predicted | Observed | Error |
|-------|-----------|----------|-------|
| theta_12 | 33.7° | 33.41° | +0.9% |
| theta_23 | 48.5° | 49.1° | -1.2% |
| theta_13 | 8.7° | 8.54° | +1.9% |
| delta_PMNS | 109.47° | ~230° | Poorly measured |

---

## 10. NEUTRINO MASSES

### Mass scale — GWT three-factor seesaw [DERIVED, 2026-03-28]

```
M_nu = m_e^3 / (d * m_p^2) = m_e / (d * F^2) = 50.5 meV
```

**Derivation: mode-basis perturbation theory on the d=3 lattice.**

The key difference from the standard seesaw: in GWT, the left-handed neutrino
(e_L) does NOT pre-exist as a field in the Lagrangian. The electron is breather
mode n=16 with definite chirality. Its opposite-chirality partner must be
CREATED by the kink interaction. This requires three factors of m_e/m_p
instead of two.

**The standard seesaw (2 factors):**
```
In standard physics, e_L is a pre-existing field in the Lagrangian.
You write a Dirac mass m_D and Majorana mass M_R:
  m_nu = m_D^2 / M_R       (2 factors of m_D)

If m_D = m_e, M_R = m_p:
  m_nu = m_e^2 / m_p = 278 eV    (too large by 10^4)
```

**The GWT seesaw (3 factors):**
```
In GWT, e_L is NOT a pre-existing mode. It must be CREATED.
Three factors of m_e/m_p, each with a distinct physical origin:

  Factor 1: m_e/m_p = e_R -> kink coupling
    The electron (breather) scatters into the kink (proton).
    Coupling = m_e (breather energy). Gap = m_p (kink mass).

  Factor 2: m_e/m_p = kink -> e_L coupling
    The kink produces the opposite-chirality state.
    Same coupling, same gap. This is the CHIRALITY FLIP.

  Factor 3: m_e/m_p = CREATION of the e_L mode
    In standard seesaw, e_L already exists — you just need to
    connect to it. In GWT, e_L is not any breather mode.
    It exists only as a virtual excitation INSIDE the kink
    interaction. Generating it costs one more mode-density gap.

  Axis averaging: 1/d (kink selects one of d equivalent axes)

  M_nu = m_e * (m_e/m_p) * (m_e/m_p) * (1/d)
       = m_e^3 / (d * m_p^2)
       = m_e / (d * F^2)
       = 50.5 meV
```

**Why e_L must be CREATED (the key GWT insight):**

In GWT, particles are specific excitation modes of the lattice. The electron
is breather mode n=16 — a definite mode with definite chirality. The 24
breather modes of the sine-Gordon spectrum account for all quarks and leptons.
The neutrino is NOT one of these 24 modes.

The neutrino exists only as a topologically-generated virtual state: when
the electron traverses the kink an ODD number of times, the intermediate
state has opposite chirality. This state cannot exist as a free propagating
mode — it is permanently confined to the kink interaction region.

The creation probability = m_e/m_p = one additional traversal of the
mode-density gap. This is the third factor that distinguishes GWT from
standard seesaw.

**Why the MODE basis, not position-space wavefunctions:**

The derivation works in the breather/kink mode basis, NOT in position-space.
Position-space Rayleigh-Schrödinger PT on the Pöschl-Teller wavefunctions
gives a divergent result because it double-counts: the potential that creates
the bound state is used again as the perturbation. In the mode basis, the
breather energy m_e and kink mass m_p are already the physical states.
The coupling between them IS m_e (the breather's energy = its coupling to
the kink that created it). The gap IS m_p (the kink mass = the next mode).

**Chirality = parity of kink traversals:**
```
The kink is a topological defect: phi goes 0 -> 2.
Each traversal FLIPS the field orientation.

  Even flips (0, 2, 4, ...): same chirality as electron (e_R)
  Odd flips  (1, 3, 5, ...): OPPOSITE chirality = neutrino (e_L)

Standard seesaw (2 traversals):
  e_R ->[kink]-> p ->[kink]-> e_R    (even = same chirality = VP correction)

GWT seesaw (3 traversals):
  e_R ->[kink]-> p ->[create e_L]-> e_L ->[kink]-> e_R
  The e_L creation is the THIRD traversal = odd = opposite chirality.
```

**Four neutrino properties from the mode-basis seesaw:**
1. **Left-handed:** e_L is the odd-traversal partner, opposite chirality
2. **Tiny mass:** THREE factors of m_e/m_p ≈ 5.4×10⁻⁴ gives (m_e/m_p)² ≈ 3×10⁻⁷
3. **Three flavors:** one per spatial axis (d = 3 equivalent kink orientations)
4. **Weak-only interaction:** the e_L exists only inside the kink interaction =
   the SU(d-1) rotational mode of the torus = the weak force

**Comparison to standard physics:**
```
Standard seesaw:  m_nu = m_D^2/M_R    (e_L pre-exists, 2 factors)
GWT seesaw:       m_nu = m_e^3/(d*m_p^2) (e_L created, 3 factors)

The extra factor m_e/m_p is WHY neutrinos are 10^4 lighter than
standard seesaw predicts (278 eV vs 50 meV). Standard seesaw needs
M_R >> m_p (a new high-energy scale) to get tiny neutrino masses.
GWT needs no new scale — the extra suppression comes from the
topological creation of e_L.
```

**Derivation status: [DERIVED]**
Every factor traces to the Lagrangian:
- m_e: breather bound state energy [DERIVED]
- m_p: F × m_e, F = 2d·π^(2d-1) [DERIVED]
- 1/d: axis averaging [STRUCTURAL]
- Third factor: e_L creation, does not pre-exist as a lattice mode [TOPOLOGICAL]

See: `calculations/core/neutrino_seesaw.py` for the full computation.

### Corrections to base mass [PATTERN — not yet formally derived]

The following corrections give excellent numerical results (0.1-2% on splittings)
but their derivation from the perturbation chain is incomplete. The factors
are geometrically motivated but the formal matrix element calculation connecting
them to the PT wavefunctions has not been done.

**Gauge gate correction:**
```
M_eff = M_nu * (1 + 1/(N_gauge * pi)) = M_nu * (1 + 1/(12*pi)) = 51.9 meV
```
1/(|A₄| × π): the PT chain passes through |A₄| = 12 gauge channels, each
contributing over one half-period (π) of the cosine potential. Physically
motivated but the 1/(|A₄|·π) product needs formal derivation from the
overlap integrals.

**Topological mode count:**
```
N_top = d * 2^d + 1 = |O| + 1 = 25      (proper cube rotations + identity)
N_eff = N_top * (1 + 1/(2*pi^2)) = 26.27   (V_0/2 = average potential perturbation)
```
d·2^d = 24 = |O| = chiral octahedral group. +1 for the vacuum. The correction
V₀/2 = 1/(2π²) comes from the Lagrangian potential depth. The connection between
N_eff and the splitting formulas below is pattern-matched, not derived.

### Mass splittings [PATTERN — 0.1% on ratio]
```
Delta_m^2_31 = (1 - 1/N_eff) * M_eff^2 = 2.586 × 10^-3 eV^2
  Observed: 2.534 × 10^-3 eV^2. Error: +2.1%

Delta_m^2_21 = (d/(4*N_eff)) * M_eff^2 = 7.68 × 10^-5 eV^2
  Observed: 7.53 × 10^-5 eV^2. Error: +2.0%

Ratio: Delta_m^2_31 / Delta_m^2_21 = 33.69
  Observed: 33.65. Error: +0.1%
```

### Individual masses
| State | Formula | Mass | Status |
|-------|---------|------|--------|
| nu_3 | M_eff | 51.9 meV | [PATTERN] |
| nu_2 | sqrt(m_1^2 + Delta_m^2_21) | 13.4 meV | [PATTERN] |
| nu_1 | M_eff / sqrt(N_eff) | 10.1 meV | [PATTERN] |
| Sum | | 75.4 meV (< 120 meV cosmo bound) | |

### Wave sizes (Compton wavelength)
```
lambda_C = hbar*c / m_nu

nu_3: hbar*c / 51.2 meV = 3.85 um
nu_2: hbar*c / 13.2 meV = 14.93 um
nu_1: hbar*c / 10.0 meV = 19.75 um
```
These are enormous by particle physics standards — comparable to biological cells.

### Ghostliness ratio
```
Ghostliness = lambda_C(nu_3) / r_weak
            = 3.85 um / (hbar*c / M_Z)
            = 3.85e-6 m / 2.16e-18 m
            = 1.78 × 10^12
```
The neutrino wave is nearly two trillion times larger than the range over which it can interact. This explains cross sections of order sigma ~ 10^-44 cm^2 — the wave simply does not "fit" into the interaction region.

### Lepton radii (classical radius)
```
r_rms = alpha * hbar*c / m

Electron: alpha * 197.3 / 0.511 = 2.818 fm   (obs: 2.8179 fm, exact)
Muon:     alpha * 197.3 / 105.7 = 0.01362 fm  (awaiting MUSE @ PSI)
Tau:      alpha * 197.3 / 1777  = 8.10e-4 fm   (awaiting)
```
The classical radius r = alpha * lambda_C is the distance at which the EM self-energy equals the rest mass. In GWT, this is the breather's electromagnetic interaction radius — the toroidal vortex has a sharp boundary at this scale. Free leptons (not confined) receive no VP correction.

---
