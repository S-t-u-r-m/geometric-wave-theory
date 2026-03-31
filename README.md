# Geometric Wave Theory (GWT)

**One Lagrangian. One input (d = 3). Zero free parameters. 55+ predictions.**

GWT starts from one premise: space is a discrete elastic medium. The simplest Lagrangian on a d = 3 cubic lattice —

```
L = (1/2)(∂φ)² + (1/π²)(1 − cos(πφ))
```

— produces the Standard Model particle spectrum, all coupling constants, the cosmological parameters, and chemical bonding. No free parameters. Every prediction is a closed-form expression in d, π, and 2.

📖 **[Website](https://s-t-u-r-m.github.io/geometric-wave-theory/)** · 📄 **[Complete Reference](gwt_complete_reference.md)** · 🔬 **[Papers](papers/)**

## Quick Navigation

| Topic | Derivation | Calculations |
|-------|-----------|--------------|
| [Foundation & Lagrangian](reference/foundation.md) | Born rule, gauge group, harmonic space | — |
| [Forces & 1D energy](reference/forces.md) | Hooke's law, 1/3-2/3 split | — |
| [Coupling constants](reference/coupling_constants.md) | α, α_s, sin²θ_W, VP law, Gray codes | [calculations/coupling/](calculations/coupling/) |
| [Mass ratios & spectrum](reference/mass_ratios.md) | m_p/m_e, breathers, bosons, Koide | [calculations/masses/](calculations/masses/) |
| [Mixing & neutrinos](reference/mixing_neutrinos.md) | CKM, PMNS, neutrino seesaw | [calculations/core/](calculations/core/) |
| [Standing wave structure](reference/standing_waves.md) | Kinks, breathers, bound configurations, no particles inside | — |
| [Relativity](reference/relativity.md) | SR/GR from springs: c, E=mc^2, Lorentz, equivalence, BH | — |
| [Magnetism](reference/magnetism.md) | B-field = transverse twist (T1g), dia/para/ferromagnetism | — |
| [Toroidal physics](reference/toroidal_physics.md) | Quarks, confinement, torus structure | [calculations/simulations/](calculations/simulations/) |
| [Cosmology](reference/cosmology.md) | Ω_Λ, H₀, MOND, BH thermodynamics | [calculations/cosmology/](calculations/cosmology/) |
| [Atomic & molecular](reference/atomic_molecular.md) | H atom, Z_eff, quantum defects, spectra | [calculations/atomic/](calculations/atomic/) |
| [Bonding](reference/bonding.md) | V8 formula, Oh bonds, Morse, 3D ZPE | [calculations/bonding/](calculations/bonding/) |
| [Nuclear & moments](reference/nuclear.md) | Proton radius, mesons, g-2, D4h table | [calculations/core/](calculations/core/) |
| [Lattice & symmetry](reference/lattice_and_symmetry.md) | 8 modes, why d=3, Oh N-body, scorecard | [calculations/atomic/](calculations/atomic/) |
| [Plain-language guide](reference/gwt_simple_guide.md) | Non-technical overview of the theory | — |

---

## Headline Results

| Prediction | GWT Formula | GWT Value | Observed | Error |
|------------|-------------|-----------|----------|-------|
| [Proton-electron mass ratio](reference/mass_ratios.md#proton-electron-mass-ratio) | 6π⁵(1 + α²/2√2) | 1836.15267 | 1836.15267 | **< 0.001 ppm** |
| [Fine structure constant](reference/coupling_constants.md#fine-structure-constant-alpha) | instanton on d-cube | 1/137.036 | 1/137.036 | **0.66 ppm** |
| [Electron g-2](reference/nuclear.md#electron-g-2-three-terms-032-ppm) | (α/2π)(1 − α/5 − α²/7) | 0.00115965182 | 0.00115965218 | **0.32 ppm** |
| [Strong coupling α_s](reference/coupling_constants.md#strong-coupling-alpha_s-formal-derivation-from-the-lagrangian) | d²/(2^d π²) × VP | 0.11794 | 0.11790 | **0.030%** |
| [Proton radius](reference/nuclear.md#proton-charge-radius-derived-002) | (d+1)ℏc/m_p | 0.8412 fm | 0.8414 fm | **0.02%** |
| [Pion mass](reference/nuclear.md#pion-mass--direct-formula-derived-04) | m_p × 4/27 | 139.0 MeV | 139.6 MeV | **0.4%** |
| [Neutron-proton mass diff](reference/nuclear.md#neutron-proton-mass-difference-derived-0005) | m_e × 8/3 × (1−7α) | 1.293 MeV | 1.293 MeV | **0.005%** |
| [Dark energy Ω_Λ](reference/cosmology.md) | (d−1)/d | 0.667 | 0.685 | 2.7% |
| [Hubble constant (CMB)](reference/cosmology.md#h_0cmb--670-kms--the-torus-correction-derived-2026-03-29) | H₀ × d³/(d³−1) | 67.0 km/s/Mpc | 67.4 | **0.6%** |
| [Hubble constant (local)](reference/cosmology.md#h_0local--738-kms--the-mond-bias-derived-2026-03-29) | H₀ × (d/(d−1))^(1/d) | 73.8 km/s/Mpc | 73.0 | **1.1%** |
| [Gravitational constant](reference/nuclear.md#gravitational-constant-the-hierarchy-problem-solved) | F⁴ × α²⁴ | 5.903×10⁻³⁹ | 5.906×10⁻³⁹ | **0.05%** |
| [Water bond angle](reference/atomic_molecular.md#h2o-bond-angle) | arccos(−1/(d+1)) | 104.48° | 104.45° | **0.03%** |

---

## How It Works

**Particles are wave patterns on the lattice.** The cosine potential supports topological solitons (kinks = protons) and bound oscillations (breathers = electrons). The d = 3 cubic lattice has octahedral symmetry (Oh), which determines everything:

- **24 breather modes** = 24 fermions of the Standard Model (|O| = [chiral octahedral group](reference/foundation.md))
- **Gauge group SU(3)×SU(2)×U(1)** = internal, rotational, and phase modes of a [d = 3 torus](reference/toroidal_physics.md)
- **3 generations** = 3 spatial axes of the cube
- **Quark charges 1/3, 2/3** = flow fractions along one vs two axes of the torus
- **Confinement** = [Perron-Frobenius theorem](reference/toroidal_physics.md) (A1g ground state is unique, deviation costs energy linearly)
- **Born rule** P = cos²θ = [projection of T1u twist modes](reference/foundation.md) onto the d = 3 cube

### The Universal VP Law

One mechanism — φ⁴ perturbation theory on the cosine potential — [dresses all fundamental constants](reference/coupling_constants.md):

```
quantity_dressed = quantity_bare × (1 ± α² × (d²−1)/denominator)
```

| Constant | Denominator | Precision | Physics |
|----------|-------------|-----------|---------|
| [m_p/m_e](reference/mass_ratios.md#proton-electron-mass-ratio) | 2^(d/2) = 2√2 | < 0.001 ppm | Confined VP (DFT on cube vertices) |
| [1/α](reference/coupling_constants.md#fine-structure-constant-alpha) | d² = 9 | 0.66 ppm | Free VP (photon coupling dimensions) |
| [α_s](reference/coupling_constants.md#strong-coupling-alpha_s-formal-derivation-from-the-lagrangian) | d = 3 | 0.030% | Gluon VP (per color channel) |

Same 8 non-A1g channels from T1u⊗T1u = A1g + Eg + T1g + T2g. Different denominators because different confinement physics. No Feynman diagrams — just nonlinear springs.

---

## Full Scorecard

| Category | Count | Mean Error | Range |
|----------|-------|------------|-------|
| [Coupling constants](reference/coupling_constants.md) (α, α_s, sin²θ_W) | 3 | 0.07% | 0.0001% – 0.15% |
| [Mass ratios & VP law](reference/mass_ratios.md) | 6 | 0.13% | < 0.001 ppm – 0.61% |
| [Fermion masses](reference/mass_ratios.md) (breather spectrum) | 9 | 1.4% | 0.02% – 3.1% |
| [Boson masses](reference/mass_ratios.md) (W, Z, Higgs) | 4 | 0.01% | 0.00% – 0.03% |
| [Lepton masses](reference/mass_ratios.md) (Koide, zero params) | 3 | 0.04% | 0.007% – 0.11% |
| [Meson spectrum](reference/nuclear.md) (π, K, ρ, ω) | 5 | 0.6% | 0.02% – 1.5% |
| [CKM matrix elements](reference/mixing_neutrinos.md#ckm-matrix-zero-free-parameters) | 4 | 0.2% | 0.03% – 0.35% |
| [PMNS mixing angles](reference/mixing_neutrinos.md#pmns-matrix-zero-free-parameters) | 3 | 1.3% | 0.9% – 1.9% |
| [Neutrino mass splittings](reference/mixing_neutrinos.md#mass-splittings--from-oh--d4h-restriction-derived-2026-03-28) | 3 | 0.4% | 0.1% – 0.5% |
| [Nuclear moments](reference/nuclear.md#proton-magnetic-moment-with-pion-cloud-003) (μ_p, μ_n, g_A) | 4 | 0.1% | 0.03% – 0.20% |
| [Nuclear binding](reference/nuclear.md#nuclear-binding-energy-per-nucleon-derived-041) (B/A, deuteron) | 3 | 1.5% | 0.4% – 2.9% |
| [Atomic](reference/atomic_molecular.md) (Rydberg, fine structure, spectra) | 3 | 0.08% | 0.004% – 0.20% |
| [Molecular bonds & angles](reference/bonding.md) | 25 | 7.8% | 0.03% – 20% |
| [Cosmological](reference/cosmology.md) (Ω_Λ, H₀, η_B, MOND) | 8 | 2.7% | 0.3% – 9.1% |
| Structural (generations, colors, gauge group) | 8 | 0% | exact |

**Standard Model: 19 free parameters. GWT: 0.**

---

## Key Derivations

Every formula is closed-form. Compare:

| Quantity | Standard Physics | GWT |
|----------|-----------------|-----|
| [m_p/m_e](reference/mass_ratios.md) | Lattice QCD (supercomputer, years) | 2d × π^(2d−1) |
| [α_s(M_Z)](reference/coupling_constants.md) | 5-loop pQCD + lattice Monte Carlo | d²/(2^d × π²) |
| [Ω_Λ](reference/cosmology.md) | Measured, unexplained | (d−1)/d |
| [sin²θ_W](reference/coupling_constants.md#weak-mixing-angle-sin2theta_w) | Measured, unexplained | 15/64 |
| [η_B](reference/cosmology.md#baryon-asymmetry-derivation) (baryon asymmetry) | Electroweak baryogenesis (unsolved) | J × α² × d/2^d |
| [Bond angle of water](reference/atomic_molecular.md#h2o-bond-angle) | Numerical quantum chemistry | arccos(−1/(d+1)) |

### Proton Radius Puzzle — Resolved

GWT predicts r_p = [(d+1)ℏc/m_p](reference/nuclear.md#proton-charge-radius-derived-002) = 0.8412 fm, matching the muonic hydrogen measurement — not the old electronic value (0.875 fm). The puzzle was never a puzzle: the muonic measurement was correct. With toroidal VP dressing: r_p = 0.84062 fm (0.0001%).

### Why d = 3 Is Unique

Three algebraically independent expressions all equal 12 **only at d = 3**:

| Count | Formula | d=2 | d=3 | d=4 |
|-------|---------|-----|-----|-----|
| Gauge channels | 2d(d−1) | 4 | **12** | 24 |
| Even permutations of spacetime | (d+1)!/2 | 3 | **12** | 60 |
| Half the breather spectrum | ⌊2^d π−1⌋/2 | 5 | **12** | 24 |

The equation (d+1)!/2 = 2d(d−1) has unique solution d = 3. The "mysterious 3s" in physics (generations, colors, quarks, 1/3 charges) are all the same geometric fact. [Full derivation →](reference/lattice_and_symmetry.md)

---

## Falsifiable Predictions

| Prediction | GWT Value | How to Test |
|------------|-----------|-------------|
| [Proton radius](reference/nuclear.md) (μ-p = e-p) | 0.84062 fm | MUSE at PSI (~2027) |
| [Muon is point-like](papers/gwt_muse_prediction.md) | R_μ < 10⁻²⁰ fm | MUSE at PSI |
| [Dark energy EOS](reference/cosmology.md#dark-energy-equation-of-state) | w = −1 exactly, w_a = 0 | DESI, Euclid |
| [Neutrino mass sum](reference/mixing_neutrinos.md#individual-masses-derived) | 74.4 meV | KATRIN, cosmological bounds |
| [Hubble tension](reference/cosmology.md#testable-predictions-from-the-hubble-tension-resolution) | Ratio = 1.102 (permanent) | GW sirens, TRGB, future CMB |
| No new particles below Planck | — | LHC, future colliders |

---

## Repository Structure

```
gwt_complete_reference.md      Master equation sheet + index to all derivations

reference/                     Detailed derivations by topic (16 files)
  foundation.md                  Lagrangian, Born rule, harmonic space, gauge group
  forces.md                      Hooke's law, 1/3-2/3 split, 1D energy
  coupling_constants.md          α, α_s, sin²θ_W, VP law, Gray codes
  mass_ratios.md                 m_p/m_e, breather spectrum, bosons, Koide
  mixing_neutrinos.md            CKM, PMNS, neutrino masses
  standing_waves.md              Kinks, breathers, bound configurations, shell structure
  relativity.md                  SR/GR from springs: c, E=mc^2, Lorentz, equivalence, BH
  magnetism.md                   B = transverse twist, dia/para/ferromagnetism
  toroidal_physics.md            Torus structure, quarks, confinement
  cosmology.md                   Ω_Λ, H₀, MOND, baryon asymmetry, BH thermodynamics
  atomic_molecular.md            H atom, IE, Z_eff, quantum defects, spatial Z_eff
  bonding.md                     V8 formula, Oh bonds, Morse, 3D ZPE
  nuclear.md                     Proton radius, mesons, g-2, D4h table, moments
  lattice_and_symmetry.md        8 modes, why d=3, N-body Oh, scorecard
  toroidal_breathers.md          Breather dynamics on the torus
  gwt_simple_guide.md            Plain-language guide to the theory

calculations/                  Computational models (organized by topic)
  core/                          Master parameter registry, muon g-2, neutrino seesaw
  coupling/                      α, α_s, spring constant derivations
  masses/                        Koide, generation mass calculations
  bonding/                       Bond energy models and Hessian analysis
  atomic/                        Ionization energies, Oh N-body, spectrum predictor
  simulations/                   3D GPU lattice, breather dynamics, torus
  vacuum_polarization/           VP exact and closed-form computations
  cosmology/                     Baryogenesis, Hawking radiation

papers/                        Publications
website/                       Interactive presentation (GitHub Pages)
```

---

## Papers

- **Mass Ratio**: [gwt_mass_ratio.pdf](papers/gwt_mass_ratio.pdf) — Derives m_p/m_e = 6π⁵ from phase space on the half-BZ. Includes VP law, g-2, bond emergence, and 11 sections of supporting derivations.
- **Mixing Matrices**: [gwt_mixing_matrices.pdf](papers/gwt_mixing_matrices.pdf) — CKM and PMNS matrices from octahedral geometry with zero free parameters. All 9 CKM elements within 1.4σ.
- **MUSE Prediction**: [gwt_muse_prediction.md](papers/gwt_muse_prediction.md) — Falsifiable predictions for the MUSE experiment at PSI: proton radius, muon point-like, lepton universality. Includes derived muon g-2 (-0.85 ppm).
- **Lagrangian**: [gwt_lagrangian.pdf](papers/gwt_lagrangian.pdf) — The foundation: one Lagrangian on a d=3 cubic lattice produces all of physics with zero free parameters.

---

## Cite This Work

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19079880.svg)](https://zenodo.org/records/19079880)

```bibtex
@software{wollenberg2026gwt,
  author       = {Wollenberg, Jonathan},
  title        = {Geometric Wave Theory},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.19079880},
  url          = {https://zenodo.org/records/19079880}
}
```

## License

CC BY-SA 4.0 — Share and build on this work with attribution.

## Author

**Jonathan Wollenberg**
- ORCID: [0009-0009-5872-9076](https://orcid.org/0009-0009-5872-9076)
- GitHub: [S-t-u-r-m](https://github.com/S-t-u-r-m)
