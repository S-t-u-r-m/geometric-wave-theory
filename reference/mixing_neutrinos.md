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

### Mass scale (third-order perturbative seesaw)
```
M_nu = m_e^3 / (d * m_p^2)

Using GWT-predicted m_p = 6*pi^5 * m_e:
M_nu = m_e^3 / (d * m_p^2) = m_e / (d * (6*pi^5)^2) ≈ 49.9 meV
```
Third-order perturbation: electron → proton → electron, averaged over d axes.

### Gauge gate correction (lattice-derived)
```
M_eff = M_nu * (1 + 1/(N_gauge * pi)) = M_nu * (1 + 1/(12*pi)) = 51.2 meV
```
1/(N_gauge * pi) = 1/(|A_4| * pi): one gauge gate contribution over one half-period of the cosine potential. Previously labeled "Wyler transverse sphere" — now pure lattice geometry.

### Topological mode count
```
N_top = d * 2^d + 1 = |O| + 1 = 25      (proper cube rotations + identity)
N_eff = N_top * (1 + 1/(2*pi^2)) = 26.27   (V_0/2 = average potential perturbation)
```
d*2^d = 24 = |O| = order of chiral octahedral group (proper rotations of the cube). +1 for the vacuum. The correction 1/(2*pi^2) = V_0/2 = half the Lagrangian potential depth. Previously labeled "D_IV(5) Shilov boundary" — now pure lattice quantities.

### Mass splittings
```
Delta_m^2_31 = (1 - 1/N_eff) * M_eff^2 = 2.523 × 10^-3 eV^2
  Observed: 2.534 × 10^-3 eV^2. Error: -0.4%

Delta_m^2_21 = (d/(4*N_eff)) * M_eff^2 = 7.49 × 10^-5 eV^2
  Observed: 7.53 × 10^-5 eV^2. Error: -0.5%

Ratio: Delta_m^2_31 / Delta_m^2_21 = 33.69
  Observed: 33.65. Error: +0.1%
```

### Individual masses
| State | Formula | Mass |
|-------|---------|------|
| nu_3 | M_eff | 51.2 meV |
| nu_2 | sqrt(m_1^2 + Delta_m^2_21) | 13.2 meV |
| nu_1 | M_eff / sqrt(N_eff) | 10.0 meV |
| Sum | | 74.4 meV (< 120 meV cosmological bound) |

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
