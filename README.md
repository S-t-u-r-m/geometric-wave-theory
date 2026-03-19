# Geometric Wave Theory (GWT)

**One Lagrangian. One input (d=3). Zero free parameters. ~55 predictions.**

## What is GWT?

Geometric Wave Theory starts from one question: *what if space is a discrete elastic medium?*

The simplest possible Lagrangian on a d=3 cubic lattice:

```
L = (1/2)(dphi)^2 + (1/pi^2)(1 - cos(pi * phi))
```

produces:

- **Proton-electron mass ratio**: 6π⁵ × (1 + α²/2^(d/2)) = 1836.15267 (**< 0.001 ppm**)
- **Fine structure constant**: bare 1/137.042, dressed 1/137.036 (**0.66 ppm**)
- **Strong coupling**: α_s = 0.11794 (**0.030%**)
- **Electron g-2**: α/(2π) × (1 - α/5 - α²/7) = 0.00115965182 (**0.32 ppm**)
- **Gravitational constant**: α_G = F⁴ × α²⁴ (**0.05%**, hierarchy problem solved)
- **All 24 fermion masses** from the breather spectrum
- **CKM and PMNS mixing matrices** from octahedral group geometry
- **Ionization energies for 103 atoms** at 2.6% mean error
- **Bond energies** with ZPE correction (H₂ to 0.1%)

Everything derives from d=3. No observed values are used as inputs.

## The Universal VP Law (2026-03-18)

A single mechanism — second-order perturbation theory on the nonlinear spring — gives four fundamental constants:

```
V = (1/pi^2)(1 - cos(pi*phi))        Expand: phi^4 nonlinearity
    -> T1u scatters into T1u x T1u    9 channels on the d=3 cube
    -> A1g = secular (already in bare) 1 channel
    -> Non-A1g = correction            8 channels = (d^2-1) modes
    -> correction = alpha^2 × 8/denominator
```

| Constant | Formula | Precision | Denominator |
|----------|---------|-----------|-------------|
| m_p/m_e  | 6π⁵(1+α²/2^(d/2)) | **< 0.001 ppm** | 2^(d/2) (confined, DFT on cube) |
| 1/α      | bare × (1-α²×8/9) | **0.66 ppm** | d² (coupling dimensions) |
| α_s      | bare × (1+α_s²×8/3) | **0.030%** | d (color channels) |
| g-2      | α/(2π)(1-α/5-α²/7) | **0.32 ppm** | (2d-1), (2d+1) (directional modes) |

All four use the same 8 non-A1g channels from T1u⊗T1u. The denominator differs because the physics differs (confined vs free, colored vs colorless). This is textbook nonlinear wave perturbation theory — no Feynman diagrams, just springs.

## The Physics

Particles are not point objects. They are **localized wave oscillations** (breathers) trapped on a Planck-scale lattice:

- **Kink** = nucleus (topological soliton, baryon number = winding number)
- **Breather** = bound oscillation in the kink potential (what standard physics calls "electron")
- **Screening** = wave impedance (how breather modes block each other)
- **Bond** = breather shared between two kinks (resonant wave transfer)
- **Gravity** = 1/d of Hooke's law (longitudinal component of the lattice spring)

The cosine potential `(1/pi^2)(1 - cos(pi*phi))` has exactly **24 bound breather modes** — the 24 fermions of the Standard Model. The number 24 = |O| is the order of the chiral octahedral group, the symmetry group of the d=3 cube.

## Key Results

### Fundamental Constants
| Quantity | Formula | Predicted | Observed | Error |
|----------|---------|-----------|----------|-------|
| m_p/m_e | 6π⁵(1+α²/2^(d/2)) | 1836.15267 | 1836.15267 | < 0.001 ppm |
| 1/α (dressed) | bare × (1-α²×8/9) | 137.0359 | 137.0360 | 0.66 ppm |
| α_s (dressed) | bare × (1+α_s²×8/3) | 0.11794 | 0.11790 | 0.030% |
| g-2 | α/(2π)(1-α/5-α²/7) | 0.00115965182 | 0.00115965218 | 0.32 ppm |
| α_G | F⁴ × α²⁴ | 5.903×10⁻³⁹ | 5.906×10⁻³⁹ | 0.05% |
| α(M_Z) | lattice running | 1/127.1 | 1/127.9 | 0.61% |
| sin²θ_W | 15/64 - 3α/2 | 0.22343 | 0.22337 | 0.03% |

### Atomic Physics
| Quantity | Formula | Predicted | Observed | Error |
|----------|---------|-----------|----------|-------|
| E_H | α²m_e/2 | 13.602 eV | 13.598 eV | 0.03% |
| Rydberg | α²m_e/(4πℏc) | 10,972,730 m⁻¹ | 10,973,732 m⁻¹ | 0.009% |
| Bohr radius | ℏ/(m_e·c·α) | 0.52920 Å | 0.52918 Å | 0.004% |
| Fine structure (n=2) | α²E_H/16 | 10.947 GHz | 10.969 GHz | 0.20% |
| 21cm hyperfine | (16/3)R∞cα²(m_e/m_p)μ_p | 1420.9 MHz | 1420.4 MHz | 0.03% |
| H₂ bond (D₀) | πE_H/d² - ZPE | 4.481 eV | 4.478 eV | 0.1% |
| H₂O angle | arccos(-1/(d+1)) | 104.48° | 104.45° | 0.03% |

### Nuclear Moments (pion cloud = strong VP law)
| Quantity | Formula | Predicted | Observed | Error |
|----------|---------|-----------|----------|-------|
| μ_p | (8/3)(1+α_s²×11/3) | 2.7937 μ_N | 2.7928 | +0.03% |
| μ_n/μ_p | -(2/3)(1+α_s²×2) | -0.6840 | -0.6850 | -0.14% |
| μ_n | μ_p × ratio | -1.9109 μ_N | -1.9130 | -0.11% |
| g_A | (4/3)(1-α_s²×11/3) | 1.2698 | 1.2723 | -0.20% |

The "pion cloud" IS the strong VP law — same φ⁴ scattering, same Oh channels, just α_s instead of α.

### Ionization Energies (103 atoms, H through Lr)
| Model | Mean Error | Atoms < 5% | Method |
|-------|-----------|------------|--------|
| GWT v19 | 2.61% | 91/103 | Algebraic (instant) |
| GWT v20 (Oh) | 3.02% | 87/103 | Group theory (instant) |
| Hartree-Fock | 5-15% | — | Numerical (hours/atom) |
| DFT (B3LYP) | 2-5% | — | Numerical (minutes/atom) |

### Oh Tensor Product Framework (2026-03-18)

The N-body problem on the d=3 lattice has **exact closed-form solutions**:

```python
def a1g_T1u(n):  # p-electrons
    if n % 2 == 1: return 0          # Parity theorem: odd = zero!
    return (3**n + 15) // 24          # Closed form

def a1g_T2g(n):  # d-electrons (t2g)
    return (3**n + 6 + 9*(-1)**n) // 24

def a1g_Eg(n):   # d-electrons (eg)
    return (2**n + 2*(-1)**n) // 6
```

**Key theorems:**
1. **Parity**: Odd u-count → zero A1g (half-fill is exact, odd loops vanish)
2. **Orthogonality**: Pairwise A1g table = identity matrix (10×10, 90% zeros)
3. **Universality**: All dim-3 irreps share the same formula at even n
4. **f→d is three-body**: A1g(T2g⊗T1u⊗T1u) = 1 (mediated by T1u, not direct)

### Scorecard (~55 predictions)
| Category | Count | Mean error | Range |
|----------|-------|------------|-------|
| Structural (forced) | 8 | 0% | exact |
| Coupling constants | 3 | 0.07% | 0.0001% – 0.15% |
| Fundamental constants (VP law) | 6 | 0.13% | <0.001 ppm – 0.61% |
| Fermion masses | 9 | 1.4% | 0.02% – 3.1% |
| Boson masses | 4 | 0.01% | 0.00% – 0.03% |
| Generation masses (Koide) | 3 | 0.04% | 0.007% – 0.11% |
| CKM matrix | 4 | 1.2% | 0.2% – 4.0% |
| PMNS matrix | 3 | 1.3% | 0.9% – 1.9% |
| Neutrino masses | 3 | 2.3% | 0.1% – 2.4% |
| Atomic (fine structure, Rydberg) | 3 | 0.08% | 0.004% – 0.20% |
| Molecular (bonds, angles) | 3 | 0.6% | 0.03% – 1.7% |
| Cosmological | 3 | 5.3% | 2.7% – 9.1% |
| Nuclear moments (pion cloud) | 4 | 0.1% | 0.03% – 0.20% |

## Source of Truth

All formulas, derivations, and predictions are documented in:

**[`math/gwt_complete_reference.md`](math/gwt_complete_reference.md)** — the authoritative reference

Code implementation: **[`calculations/gwt_lagrangian.py`](calculations/gwt_lagrangian.py)**

## Repository Structure

```
math/                    Derivations and reference documents
  gwt_complete_reference.md   THE source of truth (all formulas)

calculations/            Computational models and simulations
  gwt_lagrangian.py          Master Lagrangian implementation
  z_eff_final.py             IE formula: Oh screening + alpha (3.07%)
  z_eff_v19.py               IE formula: v19 calibrated (2.61%)
  z_eff_v20.py               IE v20: Oh tensor product corrections (3.02%)
  gwt_bond_algorithm.py      Bond energy algorithm (Oh-derived)
  oh_nbody.py                N-body Oh tensor product analysis
  sim_3d_gpu.py              3D GPU simulation on cubic lattice

papers/                  Write-ups for publication
website/                 Public-facing presentation
```

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
- Zenodo: [10.5281/zenodo.19079880](https://zenodo.org/records/19079880)

## Acknowledgments

The lattice is cubic. The geometry is exact. The number is 3.
