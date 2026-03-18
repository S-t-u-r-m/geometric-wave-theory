# Geometric Wave Theory (GWT)

**One Lagrangian. One input (d=3). Zero free parameters. 200+ predictions.**

## What is GWT?

Geometric Wave Theory starts from one question: *what if space is a discrete elastic medium?*

The simplest possible Lagrangian on a d=3 cubic lattice:

```
L = (1/2)(dphi)^2 + (1/pi^2)(1 - cos(pi * phi))
```

produces:

- **Fine structure constant**: alpha = 1/137.042 (bare), 1/137.036 (Wyler dressed)
- **Proton-electron mass ratio**: 6 * pi^5 = 1836.12 (observed: 1836.15, 0.002%)
- **All 24 fermion masses** from 3 generations of the breather spectrum
- **CKM and PMNS mixing matrices** from octahedral group geometry
- **Ionization energies for 103 atoms** at 3% mean error, zero parameters
- **Bond energies** at 7.6% mean error, zero parameters

Everything derives from d=3. No observed values are used as inputs.

## The Physics

Particles are not point objects. They are **localized wave oscillations** (breathers) trapped on a Planck-scale lattice:

- **Kink** = nucleus (topological soliton, baryon number = winding number)
- **Breather** = bound oscillation in the kink potential (what standard physics calls "electron")
- **Screening** = wave impedance (how breather modes block each other)
- **Anti-screening** = destructive wave interference (d and f modes on cubic lattice)
- **Bond** = breather shared between two kinks (resonant wave transfer)

The cosine potential `(1/pi^2)(1 - cos(pi*phi))` has exactly **24 bound breather modes** — the 24 fermions of the Standard Model. The number 24 = |O| is the order of the chiral octahedral group, which IS the symmetry group of the d=3 cube.

## Key Results

### Fundamental Constants (from the Lagrangian)
| Quantity | Formula | Value | Observed | Error |
|----------|---------|-------|----------|-------|
| Speed of light | c = a * sqrt(k/eta) | 1 (Planck) | 2.998e8 m/s | exact |
| Fine structure | exp(-(2/d!)(2^7/pi^2 + ln(6))) | 1/137.042 | 1/137.036 | 0.005% |
| Proton/electron | 2d * pi^(2d-1) | 1836.12 | 1836.15 | 0.002% |
| Weak mixing angle | arctan(sqrt(d/(d+1))) | 28.97 deg | 28.74 deg | 0.8% |

### Ionization Energies (103 atoms, H through Lr)
| Model | Mean Error | Atoms < 5% | Method |
|-------|-----------|------------|--------|
| GWT formula (v19) | 2.61% | 91/103 | Algebraic (instant) |
| GWT Oh screening | 3.07% | 87/103 | Group theory (instant) |
| Hartree-Fock | 5-15% | — | Numerical (hours/atom) |
| DFT (B3LYP) | 2-5% | — | Numerical (minutes/atom) |
| Slater's rules | 15-30% | — | Empirical |

### The Derivation Chain
```
Lagrangian: L = (1/2)(dphi)^2 + (1/pi^2)(1-cos(pi*phi))
    |
    v  kink solution
Poschl-Teller potential -> bound states -> quantum number n
    |
    v  d=3 cubic lattice
Oh irrep decomposition -> angular channels l = 0,1,2,3 (s,p,d,f)
    |
    v  Oh Clebsch-Gordan coefficients
Screening matrix -> Z_net (effective nuclear charge)
    |
    v  breather mass ratio sin(g)/sin(2g) = w_pi = 1/2
Alpha exponent -> mode coupling corrections
    |
    v
E_ion = (Z_net^alpha / n)^2 * E_H    [zero free parameters]
```

### N-Body Problem (new, 2026-03-17)

The N-body problem for breather modes on the d=3 lattice is **exactly tractable**:

- Oh group has **10 irreps** — all couplings decompose into these
- Three-body correction: check A1g content of triple tensor product
- **Half-fill (p+p+p) has ZERO three-body** by Oh symmetry (provable!)
- N-body series is **finite** (max 24 modes), mostly zero, rapidly convergent

This reduces the N-body wave equation to a finite sum of group-theory table lookups.

## Source of Truth

All formulas, derivations, and predictions are documented in:

**[`math/gwt_complete_reference.md`](math/gwt_complete_reference.md)** — the authoritative reference

Code implementation: **[`calculations/gwt_lagrangian.py`](calculations/gwt_lagrangian.py)**

## Repository Structure

```
math/                    Derivations and reference documents
  gwt_complete_reference.md   THE source of truth (all formulas)
  gwt_lagrangian.py          Master Lagrangian implementation

calculations/            Computational models and simulations
  z_eff_final.py             IE formula: Oh screening + alpha (3.07%)
  z_eff_v19.py               IE formula: v19 calibrated (2.61%)
  oh_nbody.py                N-body Oh tensor product analysis
  wave_sim_v2.py             1D sine-Gordon simulator
  sim_3d_gpu.py              3D GPU simulation on cubic lattice
  bond_v7.py                 Bond energy formula (7.6%)

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
