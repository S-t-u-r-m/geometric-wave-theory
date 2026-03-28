# Geometric Wave Theory (GWT)

**One Lagrangian. One input (d = 3). Zero free parameters. 55+ predictions.**

GWT starts from one premise: space is a discrete elastic medium. The simplest Lagrangian on a d = 3 cubic lattice —

```
L = (1/2)(∂φ)² + (1/π²)(1 − cos(πφ))
```

— produces the Standard Model particle spectrum, all coupling constants, the cosmological parameters, and chemical bonding. No free parameters. Every prediction is a closed-form expression in d, π, and 2.

📖 **[Website](https://s-t-u-r-m.github.io/geometric-wave-theory/)** · 📄 **[Complete Reference](gwt_complete_reference.md)** · 🔬 **[Papers](papers/)**

---

## Headline Results

| Prediction | GWT Formula | GWT Value | Observed | Error |
|------------|-------------|-----------|----------|-------|
| Proton-electron mass ratio | 6π⁵(1 + α²/2√2) | 1836.15267 | 1836.15267 | **< 0.001 ppm** |
| Fine structure constant | instanton on d-cube | 1/137.036 | 1/137.036 | **0.66 ppm** |
| Electron g-2 | (α/2π)(1 − α/5 − α²/7) | 0.00115965182 | 0.00115965218 | **0.32 ppm** |
| Strong coupling α_s | d²/(2^d π²) × VP | 0.11794 | 0.11790 | **0.030%** |
| Proton radius | (d+1)ℏc/m_p | 0.8412 fm | 0.8414 fm | **0.02%** |
| Pion mass | m_p × 4/27 | 139.0 MeV | 139.6 MeV | **0.4%** |
| Neutron-proton mass diff | m_e × 8/3 × (1−7α) | 1.293 MeV | 1.293 MeV | **0.005%** |
| Dark energy fraction Ω_Λ | (d−1)/d | 0.667 | 0.685 | 2.7% |
| Hubble constant (CMB) | H₀ × d³/(d³−1) | 67.0 km/s/Mpc | 67.4 | **0.6%** |
| Hubble constant (local) | H₀ × (d/(d−1))^(1/d) | 73.8 km/s/Mpc | 73.0 | **1.1%** |
| Gravitational constant | F⁴ × α²⁴ | 5.903×10⁻³⁹ | 5.906×10⁻³⁹ | **0.05%** |
| Water bond angle | arccos(−1/(d+1)) | 104.48° | 104.45° | **0.03%** |

---

## How It Works

**Particles are wave patterns on the lattice.** The cosine potential supports topological solitons (kinks = protons) and bound oscillations (breathers = electrons). The d = 3 cubic lattice has octahedral symmetry (Oh), which determines everything:

- **24 breather modes** = 24 fermions of the Standard Model (|O| = chiral octahedral group)
- **Gauge group SU(3)×SU(2)×U(1)** = internal, rotational, and phase modes of a d = 3 torus
- **3 generations** = 3 spatial axes of the cube
- **Quark charges 1/3, 2/3** = flow fractions along one vs two axes of the torus
- **Confinement** = Perron-Frobenius theorem (A1g ground state is unique, deviation costs energy linearly)
- **Born rule** P = cos²θ = projection of T1u twist modes onto the d = 3 cube

### The Universal VP Law

One mechanism — φ⁴ perturbation theory on the cosine potential — dresses all fundamental constants:

```
quantity_dressed = quantity_bare × (1 ± α² × (d²−1)/denominator)
```

| Constant | Denominator | Precision | Physics |
|----------|-------------|-----------|---------|
| m_p/m_e | 2^(d/2) = 2√2 | < 0.001 ppm | Confined VP (DFT on cube vertices) |
| 1/α | d² = 9 | 0.66 ppm | Free VP (photon coupling dimensions) |
| α_s | d = 3 | 0.030% | Gluon VP (per color channel) |

Same 8 non-A1g channels from T1u⊗T1u = A1g + Eg + T1g + T2g. Different denominators because different confinement physics. No Feynman diagrams — just nonlinear springs.

---

## Full Scorecard

| Category | Count | Mean Error | Range |
|----------|-------|------------|-------|
| Coupling constants (α, α_s, sin²θ_W) | 3 | 0.07% | 0.0001% – 0.15% |
| Mass ratios & VP law | 6 | 0.13% | < 0.001 ppm – 0.61% |
| Fermion masses (breather spectrum) | 9 | 1.4% | 0.02% – 3.1% |
| Boson masses (W, Z, Higgs) | 4 | 0.01% | 0.00% – 0.03% |
| Lepton masses (Koide, zero params) | 3 | 0.04% | 0.007% – 0.11% |
| Meson spectrum (π, K, ρ, ω) | 5 | 0.6% | 0.02% – 1.5% |
| CKM matrix elements | 4 | 0.2% | 0.03% – 0.35% |
| PMNS mixing angles | 3 | 1.3% | 0.9% – 1.9% |
| Neutrino mass splittings | 3 | 0.4% | 0.1% – 0.5% |
| Nuclear moments (μ_p, μ_n, g_A) | 4 | 0.1% | 0.03% – 0.20% |
| Nuclear binding (B/A, deuteron) | 3 | 1.5% | 0.4% – 2.9% |
| Atomic (Rydberg, fine structure) | 3 | 0.08% | 0.004% – 0.20% |
| Molecular bonds & angles | 25 | 8.1% | 0.03% – 40% |
| Cosmological (Ω_Λ, H₀, η_B, MOND) | 8 | 2.7% | 0.3% – 9.1% |
| Structural (generations, colors, gauge group) | 8 | 0% | exact |

**Standard Model: 19 free parameters. GWT: 0.**

---

## Key Derivations

Every formula is closed-form. Compare:

| Quantity | Standard Physics | GWT |
|----------|-----------------|-----|
| m_p/m_e | Lattice QCD (supercomputer, years) | 2d × π^(2d−1) |
| α_s(M_Z) | 5-loop pQCD + lattice Monte Carlo | d²/(2^d × π²) |
| Ω_Λ | Measured, unexplained | (d−1)/d |
| sin²θ_W | Measured, unexplained | 15/64 |
| η_B (baryon asymmetry) | Electroweak baryogenesis (unsolved) | J × α² × d/2^d |
| Bond angle of water | Numerical quantum chemistry | arccos(−1/(d+1)) |

### Proton Radius Puzzle — Resolved

GWT predicts r_p = (d+1)ℏc/m_p = 0.8412 fm, matching the muonic hydrogen measurement — not the old electronic value (0.875 fm). The puzzle was never a puzzle: the muonic measurement was correct. With toroidal VP dressing: r_p = 0.84062 fm (0.0001%).

### Why d = 3 Is Unique

Three algebraically independent expressions all equal 12 **only at d = 3**:

| Count | Formula | d=2 | d=3 | d=4 |
|-------|---------|-----|-----|-----|
| Gauge channels | 2d(d−1) | 4 | **12** | 24 |
| Even permutations of spacetime | (d+1)!/2 | 3 | **12** | 60 |
| Half the breather spectrum | ⌊2^d π−1⌋/2 | 5 | **12** | 24 |

The equation (d+1)!/2 = 2d(d−1) has unique solution d = 3. The "mysterious 3s" in physics (generations, colors, quarks, 1/3 charges) are all the same geometric fact.

---

## Falsifiable Predictions

| Prediction | GWT Value | How to Test |
|------------|-----------|-------------|
| Proton radius (μ-p = e-p) | 0.84062 fm | MUSE at PSI (~2027) |
| Muon is point-like | R_μ < 10⁻²⁰ fm | MUSE at PSI |
| Dark energy EOS | w = −1 exactly, w_a = 0 | DESI, Euclid |
| Neutrino mass sum | 74.4 meV | KATRIN, cosmological bounds |
| No new particles below Planck | — | LHC, future colliders |

---

## Repository Structure

```
gwt_complete_reference.md      Master equation sheet + index to all derivations

reference/                     Detailed derivations by topic
  foundation.md                  Lagrangian, Born rule, harmonic space, gauge group
  forces.md                      Hooke's law, 1/3-2/3 split, 1D energy
  coupling_constants.md          α, α_s, sin²θ_W, VP law, Gray codes
  mass_ratios.md                 m_p/m_e, breather spectrum, bosons, Koide
  mixing_neutrinos.md            CKM, PMNS, neutrino masses
  toroidal_physics.md            Torus structure, quarks, confinement
  cosmology.md                   Ω_Λ, H₀, MOND, baryon asymmetry
  atomic_molecular.md            H atom, IE, Z_eff, shell structure
  bonding.md                     V8 formula, Oh bonds, Morse, 3D ZPE
  nuclear.md                     Proton radius, mesons, g-2, moments
  lattice_and_symmetry.md        8 modes, why d=3, N-body Oh, scorecard

calculations/                  Computational models (organized by topic)
  core/                          Master parameter registry
  coupling/                      α, α_s, spring constant derivations
  masses/                        Koide, generation mass calculations
  bonding/                       Bond energy models and Hessian analysis
  atomic/                        Ionization energies, Oh N-body framework
  simulations/                   3D GPU lattice, breather dynamics, torus
  vacuum_polarization/           VP exact and closed-form computations
  cosmology/                     Baryogenesis from kink topology

papers/                        Publications
  gwt_mass_ratio.md              Proton-electron mass ratio derivation
  gwt_mixing_matrices.md         CKM and PMNS from Oh geometry
  gwt_muse_prediction.md         Falsifiable MUSE predictions
  gwt_lagrangian.md              The Lagrangian and its consequences

website/                       Interactive presentation (GitHub Pages)
```

---

## Papers

- **Mass Ratio**: [gwt_mass_ratio.pdf](papers/gwt_mass_ratio.pdf) — Full derivation of m_p/m_e = 6π⁵ from phase space on the half-BZ
- **Mixing Matrices**: [gwt_mixing_matrices.pdf](papers/gwt_mixing_matrices.pdf) — CKM and PMNS from octahedral geometry
- **MUSE Prediction**: [gwt_muse_prediction.md](papers/gwt_muse_prediction.md) — Falsifiable predictions for the MUSE experiment at PSI
- **Lagrangian**: [gwt_lagrangian.pdf](papers/gwt_lagrangian.pdf) — The foundation: one Lagrangian, zero parameters

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
