# Conceptual Principles — Session 2026-06-03

*Working principles articulated during late-night exploration of the
baryon mass-splitting sector. These are HYPOTHESES that emerged from
trying to derive Σ-Λ from GWT first principles. They are PHYSICALLY
COHERENT and INTERNALLY CONSISTENT, but most are not yet derived
within the framework. Marked clearly below: derived / forced-but-unproven
/ fitted / open.*

**Status framework** (per critical review):
- **derived**: explicit derivation from primitives, math written out
- **forced-but-unproven**: structure follows from framework but mechanism not explicit
- **fitted**: matches data but exponent/coefficient chosen to fit
- **open**: hypothesis, awaiting derivation

---

## 1. The Σ-Λ exploration

### Geometric framing (Jon's insight)
A baryon is one torus on the d=3 lattice. The torus has two fundamental
rotations (toroidal around big ring, poloidal around small tube). These
are the same kind of motion at different locations on the torus. The
ratio R/a = π is already framework-derived (toroidal_physics.md line 137).

The torus has natural compression asymmetry: inner edge (r = R−a)
is more compressed (smaller area) than outer (r = R+a). For R/a = π:
`(R-a)/(R+a) = (π-1)/(π+1) ≈ 0.517`

### Candidate formula (FITTED)
```
Σ - Λ ≈ ω₀ × [(R-a)/(R+a)]² = ω₀ × [(π-1)/(π+1)]²
      = 78.5 MeV  vs observed 77 MeV (+2% with m_p anchor)
```

### Why fitted, not derived
1. **Anchor-dependent**: m_p anchor gives 78.5 MeV (matches), m_Λ anchor
   gives 93.4 MeV (21% off). Choosing m_p for a strange-containing
   splitting is unjustified.
2. **Exponent chosen by trial**: cubed gives 40 MeV, squared gives 78.5;
   squared chosen because it matched. The physical motivation (v² ∝ 1/r²)
   is plausible but the derivation isn't from a Hamiltonian.
3. **Doesn't generalize**: charm analog Σ_c-Λ_c = 169 MeV (formula gives
   191, 13% off); bottom Σ_b-Λ_b = 191 MeV (formula gives 470, 146% off).

### Status: OPEN candidate
The geometric framing is correct (compression asymmetry of torus). The
formula is suggestive. A real derivation requires Hamiltonian eigenvalue
analysis of two stable torus configurations. Σ-Λ remains an **open
target** until that derivation exists.

---

## 2. Σ and Λ as energy levels of one object (FORCED-BUT-UNPROVEN)

In GWT, "mass" is total energy of a lattice configuration. There is no
separate "rest mass" concept.

Σ and Λ are not separate "particles" — they are **two stable energy
eigenstates of the same fundamental torus structure**. The 77 MeV gap
is the energy difference between these eigenstates.

This is analogous to atomic energy levels:
- H ground state and H first excited state are not "two different atoms"
- They are eigenstates of the same H atom
- The 10.2 eV Lyman-α photon = energy gap between eigenstates

The observed decay **Σ⁰ → Λ + γ** with a 77 MeV photon is literally an
atomic-style level transition at the baryon scale.

**Status**: this framing is consistent with the framework but requires
explicit construction of the two eigenstates from breather configurations
to be fully derived.

---

## 3. Shells = breathers (OPEN unification hypothesis)

The same physics may describe shells at three scales:

| Scale | "Shell" identity | Implied structure |
|-------|------------------|-------------------|
| Atomic | Electron shells = breather modes around nucleus | Filling pattern 2,6,10,14 from breather mode counting |
| Nuclear | Proton/neutron shells = breather configurations around core | Magic numbers 2,8,20,28,50,82,126 from torus structure |
| Baryonic | Quark "modes" = breathers wrapping the torus | Σ vs Λ = different wrapping configurations |

**Implication**: One mechanism (breather configurations on tori) underlies
atomic shells, nuclear shells, and baryon internal structure.

**Status**: highly promising unification but requires deriving the magic
numbers and shell patterns from breather mode counting. Currently a
research program, not a derivation.

---

## 4. Protons aren't fundamental (OPEN reframing)

Proposition: "proton" and "neutron" are not different particles but two
configurations of the same kind of object. Their charges come from how
their breather shells are arranged, not from intrinsic properties.

**Supporting observations**:
- p-n mass difference is 1.293 MeV (0.14% — anomalously small if they
  were fundamentally different particles)
- They have similar magnetic moment magnitudes (|μ_p|=2.79, |μ_n|=1.91)
- Both bind via the same "strong" coupling
- Neutron decays naturally to proton (n → p + e⁻ + ν̄)

**Implications for nuclei**:
- "Neutrons" are filler breather configurations needed to maintain stable
  nuclear structures at larger radii
- The neutron-excess curve (N>Z for heavy nuclei) reflects breather
  filling requirements, not separate-particle physics
- Pb-208 has 82p, 126n because the lead-sized torus structure REQUIRES
  126 neutral breather slots to be filled

**What needs derivation**:
1. The 1.293 MeV mass split from breather configuration energetics
2. The magnetic moments from torus current geometry
3. The optimal N(Z) curve from torus-scaling and breather-filling rules

**Status**: bold reframing, internally consistent, requires substantial
derivation work to formalize.

---

## 5. Lifetime = cycle count (FORCED-BUT-UNPROVEN)

"Time" in GWT is emergent — the parameter tracking lattice state changes.
A particle's "lifetime" is therefore not a fundamental timescale but a
**cycle count**: number of breather oscillations before reconfiguration.

```
lifetime = (number of cycles) × (cycle time)
        ≈ (1 / probability per cycle) × (ℏ / m c²)
```

**Numerical estimates** (using cycle-count framing):
- Neutron lifetime: 880 s, breather freq ~1.4×10²⁴ Hz → ~10²⁷ cycles
- Muon lifetime: 2.2 μs, breather freq ~1.6×10²³ Hz → ~3.5×10¹⁷ cycles
- Tau lifetime: 2.9×10⁻¹³ s, freq ~2.7×10²⁴ Hz → ~7.8×10¹¹ cycles

Heavier particles decay in fewer cycles (more decay channels open, more
phase space, larger probability per cycle).

**Status**: this framing is consistent with time being emergent (already
in GWT framework). Predicting actual lifetimes from cycle counts requires
computing per-cycle reconfiguration probabilities from geometric overlaps
— not yet done.

---

## 6. Decay products = maximal stable standing wave under conservation (FORCED-BUT-UNPROVEN)

When a configuration reconfigures (decay), the released energy takes
the form of **the maximal stable standing wave configuration that
satisfies all relevant conservation laws**.

### Hierarchy of stable forms (by mass)
1. **Photon**: stable, no rest mass, no charge — universal energy carrier
2. **Neutrino**: stable, ~50 meV mass, no charge, lepton number ±1
3. **Electron**: stable, 0.511 MeV, charge ∓1, lepton number ±1
4. **Pion**: 140 MeV, decays but stable on hadronic timescales
5. **Kaon, nucleon**: heavier stable hadrons

### Predicted decay products
The framework predicts decay products by:
1. List conserved quantum numbers (charge, baryon, lepton, etc.)
2. Find the maximal-mass combination of stable particles that:
   - Total rest mass < initial particle rest mass
   - Conserves all quantum numbers
3. This is what's emitted

### Verification against observed decays
| Decay | Conservation | Predicted | Observed |
|-------|--------------|-----------|----------|
| n → ... | charge 0→+1 needs −1 | p + e⁻ + ν̄ | ✓ |
| μ⁻ → ... | charge −1, lepton +1 (no hadron threshold met) | e⁻ + ν̄_e + ν_μ | ✓ |
| τ⁻ → ... | charge −1, can emit pion (140 MeV < 1777) | π⁻ν, e⁻ν̄ν, μ⁻ν̄ν | ✓ |
| Σ⁰ → ... | charge 0, gap 77 MeV (below 2m_e pair threshold suppressed) | Λ + γ | ✓ |
| Δ → ... | charge varies, gap 293 MeV (above pion threshold) | N + π | ✓ |

### Predicted thresholds for decay channels
- Gap < ~1 MeV: only photon channel
- Gap 1-140 MeV: photon dominates (e⁺e⁻ pair suppressed)
- Gap >140 MeV: pion decay opens (often dominates if quantum numbers fit)
- Gap >494 MeV: kaon channels open
- Gap >938 MeV: nucleon channels open

**Status**: principle is consistent with all observed decays. Predicting
SPECIFIC branching ratios still requires computing matrix elements from
geometric overlaps.

---

## 7. Critical review actions (from gwt_baryon_null_test.py, gwt_splitting_ledger.py)

### Baryon spectrum null test (preregistered, awaiting full run)
- At loose tolerance (±15 MeV): matches are pure grid density (p ≈ 0.34)
- At tight tolerance (±2 MeV): marginal signal (p ≈ 0.02-0.04 uncorrected)
- Multiple-comparison correction will determine if signal survives
- Need rerun with: 4-star-only PDG list, J^PC compatibility, multiple-comparison correction

### Splitting ledger findings
| Splitting | Coupling type | Axis | Coefficient |
|-----------|---------------|------|-------------|
| n-p (charge) | α (EM, U(1)) | charge | (2d+1) = 7 |
| Δ-N (light spin) | (none, geometric) | rotational | -- |
| Σ*-Σ (strange spin) | α_s (SU(3)) | color | (d+1) = 4 |
| Σ-Λ (recouple) | ? | ? | OPEN |

**Result**: coupling TYPE and SIGN follow activated axis cleanly. BASE
factor and COEFFICIENT mapping are NOT unified — still separate per
transition class. The selection principle mapping axis → integer
coefficient is what's missing.

---

## 8. What needs to happen next

These principles are coherent and probably right, but they're not
derivations. To turn them into framework-level theorems:

**Tier 1 (closes the credibility gap)**:
- Walk back "derived/proved/SOLVED" language in framework docs where claims
  are actually fitted or forced-but-unproven
- Apply strict status vocabulary across all reference docs
- Add blind-application column to vacuum correction ledger

**Tier 2 (the real physics)**:
- Build Hamiltonian for two stable torus configurations and derive Σ-Λ
  eigenvalue gap (target observed 77 MeV)
- Derive 1.293 MeV n-p mass split from breather configurations
- Derive proton/neutron magnetic moments from torus current geometry
- Derive neutron-excess curve N(Z) from torus-filling rules
- Derive nuclear magic numbers from breather shell counting

**Tier 3 (validation)**:
- Rerun baryon null test rigorously with all corrections
- Predict tau g-2 from same framework structure
- Test "decay products = maximal stable wave" principle against
  comprehensive decay tables

---

## 9. Honest meta-comment

This document captures ~3 hours of late-night exploration during which
Jon articulated several profound physical pictures. The pictures are
internally consistent and probably right at the conceptual level, but
the framework's existing math doesn't yet derive them. The work of
turning these pictures into rigorous derivations is the next research
program — not tonight's work.

The Σ-Λ candidate at 2% (with the anchor caveat) is the closest we got
to a derivation. It's a candidate, not a result.

Author: Jonathan D. Wollenberg
Date: 2026-06-03 (overnight session)
Status: working document — principles articulated, derivations pending
