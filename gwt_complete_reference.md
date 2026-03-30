# GWT Complete Reference — Every Derivation, Equation, and Prediction

**Single input: d = 3 spatial dimensions. Everything else follows.**

---

## DERIVATION STATUS KEY

Each result is classified by its derivation completeness:

- **[DERIVED]** = follows from the Lagrangian through explicit calculation.
  Every step checkable. No freedom to change the result.
- **[PROVEN]** = confirmed by simulation (Hessian eigenvalue, time evolution).
  The result emerged from dynamics, not from assumption.
- **[STRUCTURAL]** = follows from d=3 by counting (group theory, topology).
  Exact, but depends on identifying the physical system with the lattice.
- **[HYPOTHESIS]** = physically motivated identification not yet rigorously
  derived. The formula works (matches observations) but the derivation
  chain has gaps that need filling.
- **[PATTERN]** = numerically matches but the full derivation connecting
  the factors is incomplete. Not fitting — the factors are derived
  individually — but the COMBINATION lacks an explicit dynamical mechanism.

---

## 0. MASTER EQUATION SHEET

Everything uses only: **d, pi, 2, factorials, and exp.**

```
INPUT: d = 3
       L = (1/2)(dphi)^2 + (1/pi^2)(1 - cos(pi*phi))   [Lagrangian]

COUPLING:
  alpha_bare = exp(-(2/d!) * (2^(2d+1)/pi^2 + ln(2d)))   = 1/137.042  [DERIVED: instanton on d-cube]
  alpha      = alpha_bare * (1 + alpha^2*(d^2-1)/d^2)     = 1/137.036  [DERIVED: VP_self series]
               Dressing: phi^4 PT on d=3 cube, non-A1g fraction = 8/9

MASSES:
  F_bare    = 2d * pi^(2d-1)                              = 1836.118   [DERIVED: phase space on half-BZ]
  F         = F_bare * (1 + alpha^2/2^(d/2))              = 1836.153   [DERIVED: VP law]
              VP correction: sum(Q^2)=1 for d=3 (quark charge theorem)  [STRUCTURAL]
              Confined proton: DFT norm 1/2^(d/2). Free electron: no VP.
  m_e       = m(16, 32) from breather spectrum               = 0.505 MeV  [DERIVED: sine-Gordon]
              Mode-counting shortcut: F * alpha^|A_4| * m_Pl  = 0.511 MeV  [DERIVED: equivalent]
              n=16=2N/d (2/3 harmonic), p=32=(d+1)*2^d (FORCED). All from d=3.
  m_p       = F * m_e                                      = 938.3 MeV  [follows from above]
  m_p/m_e   = F = 6*pi^5 * (1 + alpha^2/(2*sqrt(2)))     = 1836.153   [< 0.001 ppm]

BREATHERS:
  m(n,p)    = (2^(d+1)/pi^2) * sin(n*g) * exp(-2^(d+1)*p/pi^2) * m_Pl  [DERIVED: sine-Gordon]
  g         = pi / (2^(d+1)*pi - 2)                                      [DERIVED: exact SG]

LEPTONS:
  m_mu/m_e  = d/((d-1)*alpha)                            = 205.6        [DERIVED: unified generation factor * EM scale]
  Koide     = (d-1)/d                                    = 2/3          [DERIVED: cube C3 symmetry + cosine identities]
  m_tau     = from Koide + m_mu (constrained)            = 1777.0 MeV   [DERIVED: 0.006%]
  (Full parametrization below uses theta_0, M — equivalent but more detailed)

PROTON STRUCTURE:
  r_p(bare)  = (d+1) * hbar*c / m_p                      = 0.8412 fm    [DERIVED: zero-mode count, 0.02%]
  r_p(dress) = (d+1)*(1-alpha*(d^3-1)/(d^3*pi^2))*hbar*c/m_p = 0.84064 fm [DERIVED: toroidal self-energy]
  r_p(meas)  = (d+1)*(1-alpha/pi^2)*hbar*c/m_p            = 0.84062 fm   [sphere projection, 0.0001%]

MUON RADIUS:
  r_mu      = point-like (< 10^-20 fm)                                   [DERIVED: no topological zero modes]
              TESTABLE: MUSE at PSI (~2027). See papers/gwt_muse_prediction.md

CORRECTION HIERARCHY:
  1D breathers (e, mu, tau): no VP corrections (point-like, free on lattice)
  3D torus (p, n): VP corrections needed (toroidal self-energy)
  Full discussion: reference/nuclear.md, reference/mass_ratios.md

BARYONS:
  m_n - m_p = m_e * (d^2-1)/d * (1 - alpha*(2d+1))      = 1.293 MeV    [DERIVED: QCD-EM, 0.005%]

MESONS:
  m_pi      = m_p * (d+1)/d^3 = m_p * 4/27              = 139.0 MeV    [DERIVED: axial * A1g fraction, 0.4%]
  f_pi      = m_pi * (d-1)/d                              = 92.7 MeV     [DERIVED: transverse/weak fraction, 0.3%]
  m_K       = sqrt(m_pi^2 + (m_p/(d-1))^2)                = 489.3 MeV    [DERIVED: mass-shell pion+strange, 0.89%]
  m_rho     = m_p * sqrt((8/pi^2)^2 + (4/27)^2)          = 773.1 MeV    [DERIVED: mass-shell kink+pion, 0.28%]
  m_omega   = m_rho * (1 + alpha*(2d-1)/d)                = 782.5 MeV    [DERIVED: rho + 5/3 shape channels, 0.02%]

GENERATIONS (Koide parametrization):
  sqrt(m_n) = M * (1 + sqrt(2) * cos(theta_0 + 2n*pi/d))               [DERIVED: cube C3 symmetry]
  Koide K = (d-1)/d forces A = sqrt(2(d-2)) = sqrt(2) at d=3            [DERIVED: cosine sum identities]
  M, theta_0 fixed by m_e and m_mu -> m_tau is a PREDICTION             [DERIVED: 0.006%]

MIXING:
  cos(d_CKM) = 1/d + 2/(d+1)!                            = 5/12         [STRUCTURAL]

NEUTRINOS:
  M_nu      = m_e^3 / (d * m_p^2)                         = 50.5 meV    [DERIVED: 3-factor mode-basis seesaw]
  Gauge gate: 1+1/(|A_4|*pi), N_eff = (|O|+1)*(1+V_0/(d-1))                    [DERIVED]
  Splittings from Oh->D4h: Dm31 = (1-1/N_eff)*M^2, Dm21 = d/((d+1)*N_eff)*M^2  [DERIVED]
  Ratio = 33.69 (obs: 33.65, 0.1%). All 3 masses, zero free parameters.
  Full derivation: reference/mixing_neutrinos.md

ELECTRON g-2:
  a_e       = (alpha/(2*pi)) * (1 - alpha/(2d-1) - alpha^2/(2d+1))       [DERIVED: Oh channel decomposition]
            = 0.00115965182                              (obs: 0.00115965218, -0.31 ppm)

MUON g-2:
  a_mu      = a_e + alpha^2/(2*pi) * (m_mu/m_pi)^2 * d/(d-1) * F_Oh * F_D4h  [DERIVED]
              F_Oh = 169/198, F_D4h = 1 + alpha*11/10. Result: -0.85 ppm (2.5x better than SM)
              Full derivation + D4h table: reference/nuclear.md

BOND ENERGY:
  D_e       = pi/d^2 * E_H                               = 4.749 eV     [PROVEN: Hessian eigenvalues]
              Full derivation: reference/bonding.md

VP SELF-ENERGY CONSTANT:
  VP_self   = -0.7588963842629                                            [DERIVED: definite integral]
              Full derivation: reference/coupling_constants.md

BLACK HOLE THERMODYNAMICS:
  T_H       = 1/(2^d * pi * M)                              = 1/(8*pi*M)   [DERIVED: boundary cell tunneling]
  S_BH      = A / 2^(d-1)                                   = A/4          [DERIVED: surface channels]
              8 = 2^d (cube vertices), 4 = 2^(d-1) (surface channels), pi from Lagrangian.
              Full derivation: reference/cosmology.md

WHY 12:
  alpha^12  = alpha^((d+1)!/2)  =  alpha^|A_4|                           [STRUCTURAL]
  |Oh| = 48 -> |O| = 24 -> |A_4| = 12
  (d+1)!/2 = 2d(d-1) has UNIQUE solution d = 3
```

---

## DETAILED DERIVATIONS — By Topic

Each topic has its own file in [`reference/`](reference/). The master equation sheet above is the quick-lookup; the files below contain the full derivation chains, proofs, and discussion.

### [1. Foundation — The Lattice and Its Lagrangian](reference/foundation.md)
The Lagrangian, lattice constants, impedance matching (k = eta → c = 1), entanglement as balanced displacement, the Born rule from T1u projection on the d=3 cube, harmonic space interpretation, 24 vacuum harmonics = Standard Model, gauge group SU(3)×SU(2)×U(1) from torus modes.

### [2. Forces — Hooke's Law Decomposition](reference/forces.md)
Master equation F = -kx, the 1/d longitudinal (gravity) vs (d-1)/d transverse (EM/dark energy) split, 1D energy principle (energy is structureless pressure in a geometric container), force channel separation at atomic scale.

### [3–4. Structural Parameters & Coupling Constants](reference/coupling_constants.md)
Structural parameters forced by d=3 (generations, colors, gauge group, Koide, sin²θ_W). Full alpha derivation: instanton action, channel selection (8/9), Gray code prefactor, dressing via universal VP law. Alpha_s from Gibbs overshoot. Weak mixing angle (tree + one-loop + MS-bar). VP channel rule, VP geometric constant, VP_self series.

### [5–8. Mass Ratios, Fermion Spectrum & Generations](reference/mass_ratios.md)
Proton-electron mass ratio (Surface × Volume on half-BZ, VP correction, ΣQ²=1 theorem). Universal breather formula m(n,p). Complete fermion mass table (9 particles). Boson masses (W, Z, Higgs). Unified mode-counting. Koide formula = (d-1)/d, lepton mass chain, self-energy corrections.

### [9–10. Mixing Matrices & Neutrino Masses](reference/mixing_neutrinos.md)
CKM matrix from quark mass ratios (Cabibbo, θ₂₃, θ₁₃, CP phase, N_top/N_br corrections). PMNS matrix from tribimaximal + geometric rotation. Neutrino mass scale (seesaw), splittings, individual masses, wave sizes, ghostliness ratio, lepton radii.

### [11. Toroidal Breather Physics](reference/toroidal_physics.md)
What particles are (toroidal circulations). Perron-Frobenius proof of A1g ground state. Explicit d-cube spectrum. Energy comparison of defects. Lattice stabilization. Three torus motions = three quantum numbers. Quarks as sub-circulations. Stability, annihilation, antimatter gravity. Toroidal coupling modes.

### [12. Cosmological Parameters](reference/cosmology.md)
Dark energy Ω_Λ = (d-1)/d. Hubble constant H₀ = 64.5 (bare), with speculative corrections for CMB (67.0) and local (73.8). Cosmic age, cosmological constant, dark energy EOS (w = -1 exact), MOND acceleration, baryon asymmetry η_B = J × α² × d/2^d. **Black hole thermodynamics**: T_H = 1/(2^d·π·M), S = A/2^(d-1) — the 8 and 4 are cube geometry.

### [13a. Atomic & Molecular Physics](reference/atomic_molecular.md)
Hydrogen atom, H₂ bond energy, H₂O bond angle, wave channel geometry, atomic shell structure from the Lagrangian. **Quantum defect formula**: δ_p = 1/25 + first×5/6 + rest×8/9, giving per-subshell Z_eff and spectral lines. Chain: δ → Z_eff → IE → E_harm → D_e. H < 0.03%, Na -1.4%, Mg -0.2%, 11/12 under 10%.

### [13b–11b. Chemical & Nuclear Bonding](reference/bonding.md)
**The consolidated bond physics file.** V8 formula (8 corrections = 8 non-A1g channels). Oh-derived bond model (7.5% zero-parameter). Complete bond model rewrite: correct proton topology (poloidal winding), 3D ZPE mechanism, collective bonding (not mode-by-mode), scale separation (nuclear vs chemical), 5/9 universal reduction. Kink well physics (Morse from Pöschl-Teller). Born rule connection to bonding.

### [14. Nuclear Physics, Mesons & Electromagnetic Moments](reference/nuclear.md)
Proton charge radius (bare, dressed, measured). Proton cavity. Pion mass (direct: m_p × 4/27). Rho, kaon, omega, complete meson spectrum. Nuclear energy scales, deuteron binding, B/A saturation. n-p mass difference. Electron g-2 (Oh channels, 0.31 ppm). Muon g-2 (D4h restriction, -0.85 ppm). D4h character table and Oh→D4h branching rules. Gravitational constant. Rydberg, Bohr radius, fine structure, 21cm. Proton/neutron magnetic moments, g_A. QCD string tension. Confinement theorem.

### [15–17. Lattice Corrections, Symmetry & N-Body Oh Framework](reference/lattice_and_symmetry.md)
1D breather discreteness corrections. 8 stable breather modes = 8 Oh channels (3D confirmed). Why d=3 (triple coincidence). Summary scorecard (55+ predictions). N-body problem: Oh tensor products, closed-form A1g content, screening selection rules, key structural theorems (parity, orthogonality, half-fill exactness), sparsity, implementation.

---

## Source of Truth

This file + the [`reference/`](reference/) directory constitute the authoritative GWT reference.

Code implementation: [`calculations/core/gwt_lagrangian.py`](calculations/core/gwt_lagrangian.py)

### Key files
| File | Contents |
|------|----------|
| **Core** | |
| calculations/core/gwt_lagrangian.py | Master Lagrangian, all parameters |
| **Coupling constants** | |
| calculations/coupling/alpha_from_lattice.py | 7-step alpha derivation from lattice tunneling |
| calculations/coupling/alpha12_derivation.py | Why alpha^12: octahedral group derivation |
| calculations/coupling/alpha_s_formal.py | Formal alpha_s derivation from Lagrangian |
| calculations/coupling/springconstant.py | Hooke's law, c, k, eta from d=3 |
| **Masses & generations** | |
| calculations/masses/koide_final.py | Koide formula derivation |
| **Bonding** | |
| calculations/bonding/bond_v8_full.py | Bond energies V8 (24 molecules) |
| calculations/bonding/bond_3d_emerge.py | Morse well emergence from Hessian |
| calculations/bonding/bare_hessian_multimode.py | Multi-mode breather coupling |
| calculations/bonding/toroidal_coupling_modes.py | Three coupling modes for bonding |
| **Atomic** | |
| calculations/atomic/z_eff_v20.py | IE v20: Oh tensor product corrections (3.02%) |
| calculations/atomic/oh_nbody.py | N-body Oh tensor product analysis |
| **Simulations** | |
| calculations/simulations/breather_3d_kink.py | 3D breather stability proof |
| calculations/simulations/torus_poloidal_winding.py | Confinement from topology (GPU) |
| calculations/simulations/torus_confinement_test.py | Wrong torus comparison (7 tachyons) |
| calculations/simulations/eigenspectrum_proof.py | Breather eigenspectrum on lattice |
| **Vacuum polarization** | |
| calculations/vacuum_polarization/breather_vp_exact.py | VP from Hessian method |
| calculations/vacuum_polarization/breather_vp_closedform.py | VP high-precision + universality |
| **Cosmology** | |
| calculations/cosmology/kink_phase_baryogenesis.py | Baryon asymmetry from kink topology |
