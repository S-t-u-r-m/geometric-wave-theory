# GWT Predictions for the MUSE Experiment

**Jonathan D. Wollenberg**
ORCID: [0009-0009-5872-9076](https://orcid.org/0009-0009-5872-9076)

March 26, 2026 — *Prediction documented prior to MUSE publication*

---

## Summary of Predictions

Geometric Wave Theory (GWT) makes the following **falsifiable predictions** for the MUSE experiment at PSI, which measures elastic muon-proton and electron-proton scattering at the same kinematics:

| Prediction | GWT Value | Falsified if |
|------------|-----------|-------------|
| Proton radius from μ-p scattering | 0.84062 ± 0.00002 fm | Differs from e-p by > 0.5% |
| Proton radius from e-p scattering | 0.84062 ± 0.00002 fm | Differs from μ-p by > 0.5% |
| Lepton universality | Exact (at tree level) | Any lepton-specific form factor |
| Muon charge radius | < 10⁻²⁰ fm (point-like) | R_μ > 0.001 fm |
| Two-photon exchange (μ vs e) | Identical structure | Lepton-mass-dependent TPE beyond QED |

**Core prediction: MUSE will find no discrepancy between μ-p and e-p extracted proton radii.** Both will yield r_p ≈ 0.841 fm, consistent with muonic hydrogen spectroscopy.

---

## 1. Why GWT Makes These Predictions

### 1.1 The Proton Radius

In GWT, the proton charge radius is derived from the zero-mode structure of a sine-Gordon kink on the d=3 cubic lattice:

$$r_p = (d+1)\left(1 - \alpha \cdot \frac{d^3-1}{d^3\pi^2}\right) \frac{\hbar c}{m_p} = 4 \times 0.99981 \times 0.21031 \text{ fm} = 0.84062 \text{ fm}$$

This is a property of the **proton** (the kink-antikink pair), not of the scattering probe. The proton looks the same to any lepton because the proton's structure is topological — it's determined by the lattice geometry, not by the probe particle.

### 1.2 The Muon Is Point-Like

In GWT, particles are classified by their topological structure on the lattice:

- **3D torus (proton, neutron):** Kink-antikink pair wrapping all 3 spatial dimensions. Has d+1 = 4 zero modes → measurable charge radius ~ 0.84 fm.
- **1D breather (electron):** Oscillation along one lattice direction. No topological winding → point-like. Charge radius at Planck scale (~10⁻³⁵ m).
- **Generation excitation (muon):** Same 1D breather structure as electron, with mass ratio m_μ/m_e = d/((d-1)α) = 205.6. Still 1D → still point-like.

The muon has **no topological zero modes** that would generate a classical charge radius. Its spatial extent is the breather width ~1/ε in Planck units, which is ~10⁻³⁵ m — seventeen orders of magnitude below MUSE's sensitivity (~10⁻¹⁸ m).

### 1.3 Lepton Universality

GWT predicts exact lepton universality at tree level because:

1. The electromagnetic coupling α is derived from the lattice tunneling rate, which is the same for all charged particles
2. The proton's form factor is a property of the proton (the kink topology), not the probe
3. The muon and electron differ only in their internal oscillation frequency, not in how they couple to the electromagnetic field

Any lepton universality violation in MUSE would require the probe particle's internal structure to affect the scattering — which in GWT would mean the muon has non-trivial topology. This would contradict the successful derivation of the muon mass as a generation excitation of the electron.

---

## 2. What MUSE Measures

MUSE (MUon Scattering Experiment) at the Paul Scherrer Institute uses the πM1 beamline to deliver simultaneous π, μ, e beams at momenta 115–210 MeV/c. It measures:

- Elastic μ⁺p and μ⁻p differential cross sections
- Elastic e⁺p and e⁻p differential cross sections (same detector, same kinematics)
- Ratio of μ-p to e-p cross sections (many systematics cancel)

The extracted observable is the proton charge form factor G_E(Q²) at low Q² ≈ 0.002–0.08 GeV², from which the charge radius is obtained via:

$$\langle r^2 \rangle = -6 \left.\frac{dG_E}{dQ^2}\right|_{Q^2=0}$$

### 2.1 Sensitivity to Muon Structure

If the muon has a charge radius R_μ, it modifies the measured cross section through the muon's own form factor:

$$F_\mu(Q^2) = 1 - \frac{Q^2 R_\mu^2}{6} + \ldots$$

This enters the cross section as a multiplicative factor, effectively shifting the extracted "proton radius" by:

$$\delta r_p^2 \approx R_\mu^2$$

For MUSE's target precision of ~1% on r_p (~0.008 fm), this means MUSE is sensitive to R_μ ≳ 0.008 fm. GWT predicts R_μ ~ 10⁻²⁰ fm, so **no muon structure will be visible**.

---

## 3. Falsification Criteria

GWT would be **challenged** if MUSE finds:

1. **r_p(μ-p) ≠ r_p(e-p)** at > 3σ significance
   - This would imply lepton non-universality or muon structure
   - GWT has no mechanism for this (both probe the same kink topology)

2. **r_p ≠ 0.841 fm** from either channel
   - The GWT value is 0.84062 fm with ~0.02% theoretical uncertainty
   - A measurement of 0.87 fm (the old CODATA value) would be problematic

3. **Anomalous Q²-dependence** in μ-p vs e-p ratio
   - A flat ratio = lepton universality (GWT prediction)
   - A Q²-dependent ratio = new physics beyond GWT

GWT would be **confirmed** if MUSE finds:

1. r_p(μ-p) = r_p(e-p) = 0.841 ± 0.01 fm (consistent with muonic H)
2. No lepton universality violation at any Q²
3. Muon consistent with point-like (no anomalous form factor)

---

## 4. Connection to the Proton Radius Puzzle

The proton radius puzzle (2010–2019) arose from a discrepancy between:
- Muonic hydrogen spectroscopy: r_p = 0.84087 ± 0.00039 fm
- Electronic hydrogen + e-p scattering: r_p = 0.8751 ± 0.0061 fm

GWT predicted the muonic value from the start (the formula gives 0.84062 fm). The resolution of the puzzle in favor of ~0.841 fm is consistent with GWT.

MUSE provides the **definitive test**: same detector, same kinematics, both leptons. If μ-p and e-p give the same r_p ≈ 0.841 fm, the puzzle is closed and GWT's prediction is confirmed.

---

## 5. The Deeper Question: Why Is the Muon Point-Like?

In the Standard Model, the muon is point-like by assumption (it's a fundamental fermion). GWT provides a **reason**: the muon is a 1D oscillation (breather) on the lattice, with no topological winding. Only 3D topological objects (kink-antikink tori) have measurable spatial extent.

This connects to the correction hierarchy discovered in GWT:
- **1D/2D modes (leptons):** No lattice distortion → no corrections to bare values → point-like
- **3D torus (proton):** Wraps all 3 dimensions → local lattice distortion → effective c reduced → correction terms (1 - α·26/(27π²)) → measurable radius

The muon mass formula m_μ/m_e = d/((d-1)α) = 205.6 (obs: 206.8, 0.6%) works **without** any spatial correction, confirming the muon has no 3D topological structure.

Supporting evidence: all three lepton masses (e, μ, τ) are predicted without corrections:
- m_μ/m_e = d/((d-1)α) = 205.6 (0.6%)
- Koide parameter = (d-1)/d = 2/3 (8.8 ppm)
- m_τ from Koide = 1777 MeV (0.006%)

---

## 6. Muon g-2: A Derived Prediction

**Update (March 28, 2026): The muon g-2 has been fully derived.**

The electron g-2 uses the full $O_h$ channel decomposition:

$$a_e = \frac{\alpha}{2\pi}\left(1 - \frac{\alpha}{2d-1} - \frac{\alpha^2}{2d+1}\right) = 0.00115965182 \quad (\text{obs: ...218, } {-0.31 \text{ ppm}})$$

The muon g-2 adds a hadronic VP correction with an NLO term from the $O_h \to D_{4h}$ subgroup restriction:

$$a_\mu = a_e + \frac{\alpha^2}{2\pi}\left(\frac{m_\mu}{m_\pi}\right)^2 \frac{d}{d-1} \cdot F_{O_h} \cdot F_{D_{4h}}$$

where:
- $F_{O_h} = \frac{169}{198} = \frac{(d^2+d+1)^2}{2d^2(d^2+d-1)}$ — bifundamental EM×QCD trace (Oh symmetry of 3D vacuum)
- $F_{D_{4h}} = 1 + \alpha \cdot \frac{d^2+d-1}{d^2+1} = 1 + \alpha \cdot \frac{11}{10}$ — NLO from D4h restriction

**Result:** $a_\mu = 0.00116591962$ (obs: $0.00116592061$, **−0.85 ppm**)

The SM predicts $0.00116591810$ (−2.15 ppm). **GWT is 2.5× closer to observation.**

**Derivation:** The muon lives on $d-1 = 2$ axes of the $d=3$ cube. Its local symmetry is $D_{4h}$ (the square), not $O_h$ (the cube). Restricting $T_{1u} \otimes T_{1u}$ from $O_h$ to $D_{4h}$:

$$T_{1u} \to A_{2u} + E_u \quad \text{in } D_{4h}$$

$$(A_{2u} + E_u) \otimes (A_{2u} + E_u) \text{ has } \textbf{2} \text{ } A_{1g} \text{ channels (not 1)}$$

The extra $A_{1g}$ comes from the $E_g$ irrep of $O_h$ splitting under $D_{4h}$: $E_g \to A_{1g} + B_{1g}$. The $d_{z^2}$ component looks scalar from the 2D muon frame. This channel is $\alpha$-suppressed and normalized by the ratio of QCD exchange paths to coupling modes: $(d^2+d-1)/(d^2+1) = 11/10$.

**Key identity:** $d^2+d-1 = d^2+2$ only at $d=3$ (requires $d-1=2$).

The D4h character table, branching rules, and full calculation are in `reference/nuclear.md` and `calculations/core/muon_g2_d4h.py`.

---

## 7. Timeline and Experimental Status

- **MUSE commissioning:** 2019–2022
- **Production data:** 2023–2025
- **First results:** Expected 2025–2026
- **This prediction documented:** March 26, 2026

All GWT predictions in this document were derived from the Lagrangian $\mathcal{L} = \frac{1}{2}(\partial\phi)^2 + \frac{1}{\pi^2}(1-\cos\pi\phi)$ on the d=3 cubic lattice with zero free parameters.

---

## References

[1] Antognini et al., Science 339, 417 (2013) — Muonic hydrogen r_p = 0.84087 fm
[2] Mohr et al., Rev. Mod. Phys. 84, 1527 (2012) — CODATA 2010 r_p = 0.8775 fm
[3] Brandt et al., Phys. Rev. Lett. (2026) — PRad-II r_p = 0.840615 fm (to be confirmed)
[4] Alexandrou et al., Phys. Rev. D 101, 114513 (2020) — Lattice QCD r_p
[5] Wollenberg, "Proton-Electron Mass Ratio from Sine-Gordon Lattice" (2026)
