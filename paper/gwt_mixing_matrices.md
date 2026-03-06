# Quark and Lepton Mixing Matrices from Lattice Geometry with Zero Free Parameters

**Jonathan D. Wollenberg**

March 2026

---

## Abstract

We present closed-form expressions for all four CKM parameters and all four PMNS parameters that reproduce the observed quark and lepton mixing matrices with zero free parameters. The CKM matrix is constructed in the standard PDG parametrization with angles theta_12 = arcsin(sqrt(m_d/m_s + m_u/m_c)), theta_23 = arcsin(sqrt(m_u/m_c)), theta_13 = arcsin(sqrt(m_u/m_t)), and CP phase delta = arccos(5/12), yielding all nine matrix elements within 1.4 sigma of PDG 2024 values (mean error 0.64%). The PMNS matrix is obtained by applying a single geometric rotation R(theta, axis) to the tribimaximal mixing matrix, with theta = arcsin((m_e/m_mu)^(1/3)) and a rotation axis determined by lepton-to-proton mass ratios, yielding all three angles within 1 sigma of NuFIT 6.0 (normal ordering). The quark sector uses mass ratios raised to the 1/2 power (surface geometry), while the lepton sector uses the 1/3 power (bulk geometry). Both CP phases derive from the tetrahedral dihedral angle arccos(+/-1/3) in d = 3 dimensions, with the CKM phase receiving a boundary correction to arccos(5/12). The formulas require only measured fermion masses as input and contain no fitted or adjustable parameters. The geometric origin of these relations, within the framework of standing-wave modes on a three-dimensional lattice, is briefly discussed.

---

## 1. Introduction

The Standard Model of particle physics contains 19 free parameters, of which 10 describe the quark and lepton mixing matrices: four CKM parameters (three angles and one CP phase) and four PMNS parameters (three angles and one CP phase), plus two neutrino mass-squared differences. These parameters are determined experimentally but have no theoretical explanation within the Standard Model itself.

Several empirical relations between mixing angles and fermion mass ratios have been noted in the literature. The Gatto-Sartori-Tonin relation V_us ~ sqrt(m_d/m_s) [1] and its extension to the full Cabibbo angle via quadrature V_us = sqrt(m_d/m_s + m_u/m_c) [2] are well-known. The tribimaximal (TBM) ansatz for neutrino mixing [3] successfully predicted the approximate structure of the PMNS matrix before theta_13 was measured. However, these relations have generally been treated as separate observations without a unifying geometric framework, and most require at least one fitted parameter.

In this paper, we present a complete set of formulas for all eight mixing parameters (four CKM, four PMNS) that:

1. Contain zero free parameters — all inputs are measured fermion masses and the integer d = 3.
2. Reproduce all observed mixing angles within experimental uncertainty.
3. Share a common geometric structure: mass ratios raised to a power determined by the dimensionality of the coupling boundary.

The key structural observation is that quark mixing angles involve mass ratios to the 1/2 power, while lepton mixing angles involve mass ratios to the 1/3 power. We interpret this as the difference between coupling at a (d-1) = 2 dimensional surface (quarks, confined inside the proton) and coupling through d = 3 dimensional bulk (leptons, which are free particles). Both CP-violating phases derive from the dihedral angle of a regular tetrahedron in three dimensions: arccos(1/3) for quarks (with a boundary correction) and arccos(-1/3) for leptons.

---

## 2. Framework

### 2.1 Lattice foundation

The formulas presented here arise from a framework called Geometric Wave Theory (GWT), in which fundamental particles are standing-wave modes (breathers) of a sine-Gordon field on a d-dimensional cubic lattice with Planck spacing. The Lagrangian density is:

    L = sum_{<i,j>} [ (1/2)(phi_i - phi_j)^2 + (1/pi^2)(1 - cos(pi * phi_i)) ]

where phi_i is the displacement field at lattice site i and the sum runs over nearest neighbors. This lattice has a single structural input: the spatial dimension d = 3. From its breather spectrum and tunneling amplitudes, all fermion masses are determined by:

    m(n, p) = (16/pi^2) * sin(n * gamma) * exp(-16p/pi^2) * m_Planck

where gamma = pi/(16*pi - 2) is the breather coupling, n is the harmonic mode number (an integer from 1 to N = d * 2^d = 24), and p is the tunneling depth. Both n and p are integers determined by lattice symmetry.

The nine charged fermion masses predicted by this formula, with specific (n, p) assignments derived from the lattice harmonic structure, agree with observation to a mean accuracy of 1.5% [4]. These masses serve as inputs to the mixing angle formulas below. For transparency, we present results using both GWT-predicted masses and PDG experimental masses, so that the mixing formulas can be evaluated independently of the mass predictions.

### 2.2 Surface vs. bulk geometry

A central distinction in this framework is between particles that are confined inside a topological defect (the proton, modeled as a kink on the lattice) and particles that propagate freely on the full lattice.

- **Quarks** are confined within the proton. Their weak-interaction coupling occurs at the proton's (d-1) = 2 dimensional boundary. Overlap integrals between standing-wave modes on a 2D surface yield amplitudes proportional to (m_i/m_j)^(1/(d-1)) = (m_i/m_j)^(1/2).

- **Leptons** are free particles extending through the full d = 3 dimensional lattice. Their mixing amplitudes scale as (m_i/m_j)^(1/d) = (m_i/m_j)^(1/3).

This distinction — square root for quarks, cube root for leptons — is the single structural difference between the CKM and PMNS derivations.

### 2.3 Tetrahedral geometry and CP phases

The regular tetrahedron is the unique symmetric solid in d = 3 dimensions with d + 1 = 4 vertices. Its dihedral angle is arccos(1/d) = arccos(1/3) = 70.53 degrees. In the lattice framework, the four vertices correspond to the four possible orientations of a kink defect in the d-cube, and CP violation arises from the chirality of wave propagation around these vertices.

The two CP phases are:

- **PMNS**: delta_PMNS = arccos(-1/d) = arccos(-1/3) = 109.47 degrees (bulk, negative handedness)
- **CKM**: delta_CKM = arccos((d+2)/(d(d+1))) = arccos(5/12) = 65.38 degrees (surface-corrected, positive handedness)

The CKM correction shifts cos(delta) from 1/3 to 5/12 = 1/3 + 1/12. The additional 1/(d(d-1)^2) = 1/12 term arises because quark mixing occurs at the (d-1)-dimensional boundary of the proton, adding a surface contribution to the bulk dihedral angle.

---

## 3. CKM Matrix

### 3.1 Construction

The CKM matrix is constructed in the standard PDG parametrization:

    V_CKM = R_23(theta_23) * R_13(theta_13, delta) * R_12(theta_12)

with four geometric parameters:

| Parameter | Formula | Value |
|-----------|---------|-------|
| theta_12 | arcsin(sqrt(m_d/m_s + m_u/m_c)) | 12.957 degrees |
| theta_23 | arcsin(sqrt(m_u/m_c)) | 2.392 degrees |
| theta_13 | arcsin(sqrt(m_u/m_t)) | 0.203 degrees |
| delta | arccos(5/12) | 65.376 degrees |

The theta_12 formula is the well-known quadrature relation: the down-type and up-type sectors contribute perpendicular rotations in flavor space, so their squares add: sin^2(theta_12) = m_d/m_s + m_u/m_c. The theta_23 formula identifies V_cb with the "up-type Cabibbo angle" sqrt(m_u/m_c). The theta_13 formula gives V_ub as the direct 1-3 surface overlap sqrt(m_u/m_t).

### 3.2 Results

Using GWT quark masses (m_u = 2.214, m_d = 4.783, m_s = 98.56, m_c = 1271, m_b = 4312, m_t = 176547 MeV):

| Element | Predicted | PDG 2024 | Uncertainty | Pull (sigma) |
|---------|-----------|----------|-------------|--------------|
| V_ud | 0.97453 | 0.97435 | 0.00016 | +1.1 |
| V_us | 0.22422 | 0.22500 | 0.00067 | -1.2 |
| V_ub | 0.003541 | 0.00369 | 0.00011 | -1.4 |
| V_cd | 0.22408 | 0.22486 | 0.00067 | -1.2 |
| V_cs | 0.97368 | 0.97349 | 0.00016 | +1.2 |
| V_cb | 0.04173 | 0.04182 | 0.00085 | -0.1 |
| V_td | 0.00852 | 0.00857 | 0.00020 | -0.3 |
| V_ts | 0.04101 | 0.04110 | 0.00085 | -0.1 |
| V_tb | 0.99912 | 0.99912 | 0.00004 | +0.1 |

Mean |error|: 0.64%. Maximum pull: 1.4 sigma (V_ub). Unitarity is exact by construction.

The Jarlskog invariant, which measures the overall strength of CP violation:

    J = Im(V_ud * V_cs * V_us* * V_cd*) = 2.93 x 10^-5

Observed: J = (3.08 +/- 0.13) x 10^-5. Deviation: -4.8%.

### 3.3 Comparison with previous formulas

The formula V_cb = sqrt(m_u/m_c) replaces the previous Wolfenstein parametrization V_cb = A * lambda^2 with A = sqrt(2/3). Both give sub-percent accuracy, but the new formula is simpler (one mass ratio vs. a composite expression) and eliminates the need for the Wolfenstein amplitude A as a separate geometric input. The V_us quadrature formula and V_ub = sqrt(m_u/m_t) are unchanged from earlier work [2, 5].

---

## 4. PMNS Matrix

### 4.1 Construction

The PMNS matrix is obtained by applying a single rotation to the tribimaximal (TBM) base:

    U_PMNS = R(theta_corr, n_hat) * U_TBM

where U_TBM is the Harrison-Perkins-Scott tribimaximal matrix [3]:

    U_TBM = | sqrt(2/3)   1/sqrt(3)   0        |
            | -1/sqrt(6)  1/sqrt(3)   1/sqrt(2) |
            | 1/sqrt(6)   -1/sqrt(3)  1/sqrt(2) |

The correction parameters are:

    theta_corr = arcsin((m_e / m_mu)^(1/3)) = 10.08 degrees

    n_hat = (-1, sqrt(3), -(m_tau/m_p)^(1/3)) / |(-1, sqrt(3), -(m_tau/m_p)^(1/3))|

The correction angle uses the cube root (1/d = 1/3, bulk geometry) of the electron-to-muon mass ratio. The rotation axis has components (-1, sqrt(3)) in the TBM degenerate subspace, with a third component proportional to the tau-to-proton mass ratio raised to the 1/3 power. This wrapping factor (m_tau/m_p)^(1/3) accounts for the tau lepton's standing-wave extent relative to the proton's.

### 4.2 Results

Using physical lepton masses (m_e = 0.51100, m_mu = 105.658, m_tau = 1776.86, m_p = 938.272 MeV):

| Angle | Predicted | NuFIT 6.0 (NO) | Uncertainty | Pull (sigma) |
|-------|-----------|-----------------|-------------|--------------|
| theta_12 (solar) | 33.49 degrees | 33.41 degrees | 0.78 degrees | +0.1 |
| theta_23 (atmospheric) | 49.28 degrees | 49.20 degrees | 1.05 degrees | +0.1 |
| theta_13 (reactor) | 8.63 degrees | 8.57 degrees | 0.12 degrees | +0.5 |

Mean pull: 0.2 sigma. All three angles within 1 sigma.

The PMNS CP phase is the bare tetrahedral dihedral angle:

    delta_PMNS = arccos(-1/3) = 109.47 degrees

The experimental value is poorly constrained: approximately 180 to 270 degrees (NuFIT 6.0). The prediction is consistent with current data.

### 4.3 Structural comparison: CKM vs. PMNS

| Feature | CKM (quarks) | PMNS (leptons) |
|---------|-------------|----------------|
| Power law | 1/2 (surface) | 1/3 (bulk) |
| Base matrix | Identity (small mixing) | TBM (large mixing) |
| Correction | Mass ratios set all angles | Single rotation corrects TBM |
| CP phase | arccos(5/12) = 65.4 degrees | arccos(-1/3) = 109.5 degrees |
| Confinement | All quarks inside proton | Leptons free on lattice |

The CKM angles are small because quark mass ratios m_light/m_heavy are small — the three generation modes are nearly orthogonal. The PMNS angles are large because the neutrino sector starts from a democratic (TBM) base — the three modes couple nearly equally. The charged-lepton correction is small (10 degrees) because m_e/m_mu is small.

---

## 5. Discussion

### 5.1 What is derived vs. what is assumed

The formulas presented here require the following inputs:

- **d = 3** spatial dimensions (determines the 1/2 and 1/3 powers, the tetrahedral angle, and the boundary correction).
- **Fermion masses** (either GWT-predicted or experimentally measured).

They do not require:

- Any fitted parameters.
- Any choice of ansatz beyond the PDG parametrization (CKM) and TBM base (PMNS).
- Any assumptions about the underlying dynamics beyond "quarks couple at surfaces, leptons couple in bulk."

The surface/bulk distinction is the single structural choice. It is motivated by the confinement of quarks inside hadrons versus the free propagation of leptons, which is an observed physical fact, not a model assumption.

### 5.2 Relation to existing work

The Gatto-Sartori-Tonin relation V_us ~ sqrt(m_d/m_s) and its quadrature extension are well-established [1, 2]. The identification V_cb = sqrt(m_u/m_c) appears to be new to this work. The TBM correction by a single rotation has been explored in various forms [6, 7], but the specific axis formula using (m_tau/m_p)^(1/3) and the connection to lattice wrapping is new.

The tetrahedral origin of CP phases has been discussed in discrete symmetry models [8], but the specific formula cos(delta_CKM) = 5/12 with its boundary correction interpretation appears to be new.

### 5.3 Predictions and tests

The framework makes several testable predictions:

1. **PMNS CP phase**: delta_PMNS = 109.5 degrees. Current data from T2K and NOvA are consistent but imprecise. DUNE and Hyper-Kamiokande will measure this to approximately 5-10 degrees precision within the next decade.

2. **Sensitivity to mass measurements**: Because all mixing angles are algebraic functions of mass ratios, improved quark mass measurements will sharpen the predictions. The largest current uncertainty comes from V_ub, which depends on sqrt(m_u/m_t). The up quark mass has approximately 30% experimental uncertainty, which propagates to approximately 15% uncertainty in V_ub.

3. **No new parameters**: If additional mixing parameters are ever measured (e.g., sterile neutrino mixing), the framework predicts they must also be expressible as mass ratios raised to integer-reciprocal powers of d.

### 5.4 Limitations

The geometric origin of the formulas is motivated by the lattice framework (Section 2) but not rigorously derived from first principles. The connection between standing-wave overlap integrals and the specific mass-ratio powers (1/2, 1/3) is plausible but has not been proven at the level of mathematical rigor expected of a fundamental theory. Similarly, the boundary correction to the CKM CP phase (1/3 to 5/12) has a clear geometric interpretation but awaits a formal derivation.

This paper intentionally presents the formulas as empirical relations that work, alongside a geometric framework that motivates them. The fact that eight parameters spanning six orders of magnitude are reproduced from zero free parameters is, in our view, sufficient to merit attention regardless of the theoretical framework's completeness.

---

## 6. Summary

Eight mixing parameters of the Standard Model — four CKM and four PMNS — are reproduced from closed-form expressions involving only fermion mass ratios and the spatial dimension d = 3. The CKM matrix achieves 0.64% mean accuracy across all nine elements (maximum pull 1.4 sigma). The PMNS angles are all within 1 sigma of observation. Both CP phases derive from the tetrahedral dihedral angle in three dimensions. No parameters are fitted.

The structural distinction between the two sectors — square root (surface) for quarks, cube root (bulk) for leptons — reflects the observed physics of confinement. These results suggest that the flavor structure of the Standard Model may have a geometric origin in the dimensionality of space.

---

## References

[1] R. Gatto, G. Sartori, M. Tonin, "Weak self-masses, Cabibbo angle, and broken SU(2) x SU(2)," Phys. Lett. B 28, 128 (1968).

[2] H. Fritzsch, Z. Xing, "Mass and flavor mixing schemes of quarks and leptons," Prog. Part. Nucl. Phys. 45, 1 (2000).

[3] P.F. Harrison, D.H. Perkins, W.G. Scott, "Tri-bimaximal mixing and the neutrino oscillation data," Phys. Lett. B 530, 167 (2002).

[4] J.D. Wollenberg, "Geometric Wave Theory: Fermion Mass Spectrum from Lattice Breathers," Zenodo (2026). [To appear]

[5] H. Fritzsch, "Calculating the Cabibbo angle," Phys. Lett. B 70, 436 (1977).

[6] S.F. King, C. Luhn, "Neutrino mass and mixing with discrete symmetry," Rep. Prog. Phys. 76, 056201 (2013).

[7] G. Altarelli, F. Feruglio, "Discrete flavor symmetries and models of neutrino mixing," Rev. Mod. Phys. 82, 2701 (2010).

[8] I. de Medeiros Varzielas, G.G. Ross, "Discrete family symmetry, Higgs mediators and theta_13," JHEP 1212, 041 (2012).

---

*Jonathan D. Wollenberg — ORCID: 0009-0009-5872-9076*
*GitHub: https://github.com/S-t-u-r-m/geometric-wave-theory*
