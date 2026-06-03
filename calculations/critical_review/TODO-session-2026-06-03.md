# TODO — Session 2026-06-03

Working notes from a critical-review session focused on the baryon resonance
spectrum, the vacuum-correction sector, and the Σ−Λ test. Organized by priority.
Tier 1 protects credibility (cheap, high-leverage). Tier 2 is the real theory
work. Tier 3 is validation and cleanup.

Status vocabulary to adopt repo-wide (see item 8): **derived** /
**forced-but-unproven** / **fitted** / **open**.

---

## Tier 1 — Framing & honesty (do first)

1. **Scoreboard split.** One table dividing every GWT comparison into:
   - *Clean observables* — particle masses, mass differences, magnetic moments,
     lifetimes, branching ratios (all directly measured; fair game, must match).
   - *Theory-laden extractions* — quark masses, running couplings, anything
     scheme-dependent (NOT clean tests; quarantine with an explicit note).
   Judge the theory only on the clean column. This stops the apples-to-oranges
   problem at the source and makes everything else defensible.

2. **Pull back "proved / exact to all orders" language** in
   `reference/baryon_resonances.md` and `reference/integrability_proof.md`.
   The DHN / dimensional-reduction sketch does not support "theorem." Downgrade
   to "suggestive; integrability conjectured." Overclaiming here is what makes a
   referee stop reading.

3. **Vacuum-correction ledger with the blind-application column.** One table:
   every correction, the breather/mode structure that *forces* it, and whether
   applying that same rule *blindly to every member of the class* (including the
   ones it makes worse) still improves the aggregate.
   - `(1 − α/24)` universal fermion correction — survives this test.
   - `(1 − πα)` on selected generations, and the bottom-quark extra `(1 − α)`
     patch — scrutinize these.
   The number of rows where you had to peek at the target = the real
   free-parameter count of the vacuum sector. (Source: `reference/vacuum_corrections.md`,
   Open Q#5 — "still no single derivation principle.")

---

## Tier 2 — The actual physics

4. **Derive the axis → integer-coefficient map.** The splitting ledger showed
   the coupling *type* and *sign* follow the activated axis cleanly (real result):
   charge → U(1) → α; strange-spin → SU(3) → α_s; light-spin → SU(2) → geometric.
   But the coefficients are still separate counts with separate labels
   (2d+1 = 7 for EM, d+1 = 4 for strong). ONE counting scheme that outputs both,
   indexed by the activated axis, is what turns "principle" into "rule."

5. **Build an internal-recoupling mechanism from the geometry** — no quark-model
   scaffolding. What does it cost, in winding/lattice terms, to change the
   relative coupling of two light sub-circulations at fixed L and fixed J?
   This is the gap Σ−Λ exposed.
   - NOTE: Σ−Λ was **never derived from GWT**. The 98 MeV was a hybrid (standard
     hyperfine spin algebra + GWT's ω₀ and cos(π/d)). GWT currently has no native
     scale for an internal recoupling. The angular ladder is for L ≥ 1; the
     strange spin-flip is for J=½→3/2; neither covers Σ vs Λ (both L=0, J=½).
   - Derive the scale, commit to a number, then test against Σ−Λ (obs 77.0 MeV)
     **and** N−Δ **and** the octet–decuplet pattern simultaneously.
   - Σ−Λ stays an **open target**, not a failed prediction, until this exists.
   - This is the single most decision-relevant piece of work: a clean,
     above-noise, parameter-free result on a measured observable the SM does not
     hand you.

6. **Fix bare-mass route consistency.** Electron (instanton/breather), muon
   (direct), tau (Koide) come from different derivations at precision *worse*
   than the ppm-level corrections being tested. Implied effective mode counts
   backed out from data: proton 23.9, electron 21.8, muon 18.5, tau 92.6 — the
   tau outlier confirms the routes are inconsistent. Until all masses come from
   one route at sub-10-ppm, the lepton-mass corrections cannot be tested (you'd
   be fitting derivation noise). The dimensional-correction idea (1D/2D/3D
   standing waves) is untestable on lepton masses for this reason — test it where
   corrections are percent-level instead.

---

## Tier 3 — Validation & cleanup

7. **Rerun the baryon null test properly.** Pre-register ONE tolerance (claimed
   intrinsic precision, ~2–3 MeV); use a 4-star-only PDG list; add J^PC
   compatibility as a match requirement; apply a multiple-comparison correction.
   Prior first-cut finding: at loose tolerance (±15 MeV) matches are pure grid
   density (p ≈ 0.34); at tight tolerance (±2 MeV) there is a marginal signal
   (p ≈ 0.02–0.04, uncorrected). The proper rerun decides if that signal is real.
   Script: `calculations/.../gwt_baryon_null_test.py` (see item 9).

8. **Organization pass — standardize status labels.** "SOLVED," "fits,"
   "validated," "derived," "candidate" are currently used loosely and sometimes
   interchangeably. Apply the strict vocabulary (derived / forced-but-unproven /
   fitted / open) everywhere so the real derivations are locatable.

9. **Commit the session scripts** so the machinery doesn't get lost again:
   - `gwt_baryon_null_test.py` — overfitting null test (random-ladder vs observed).
   - `gwt_splitting_ledger.py` — splitting decomposition + axis-rule unification test.

---

## One-line through-line

Items 1–3 keep the project from overclaiming. Items 4–5 are where a genuine win
is still available. **Item 5, tested against the clean-observable splittings, is
the one to chase** — a real, above-noise, zero-parameter prediction the Standard
Model doesn't give you for free. Don't reframe measured numbers to fit; meet them.
