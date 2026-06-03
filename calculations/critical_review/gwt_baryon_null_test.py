#!/usr/bin/env python3
"""
GWT baryon-ladder overfitting null test.

Question: does GWT's specific angular step omega_0 ~ 294 MeV match the observed
N*/Delta* spectrum better than a RANDOM ladder of identical combinatorial
structure and comparable grid density?

Method:
  - Build the GWT prediction set from one scale omega_0, anchored on m_p.
    (single angular modes, multi-quanta, two-mode mixed w/ -omega_0/(4pi) shift,
     P-wave, breathing). Everything scales with omega_0 except the fixed anchor.
  - Count how many observed PDG states each prediction set matches (within tol).
  - Null model: randomize omega_0 over a band, regenerate the SAME structure,
    record match counts. Control for grid density by also tracking #predictions.
  - p-value = fraction of random ladders (of comparable density) that match
    >= as many observed states as the real GWT ladder.

All assumptions are explicit and easy to change. This is a first cut.
"""
import numpy as np

rng = np.random.default_rng(20260603)

# ----------------------------------------------------------------------
# Fixed framework anchor (derived independently in GWT, not fit here)
# ----------------------------------------------------------------------
M_P      = 938.272          # proton anchor (MeV)
ALPHA_S  = 0.11794          # framework strong coupling
PI       = np.pi

WINDOW   = (1200.0, 3000.0) # mass window to test (MeV)
TOL_MEV  = 15.0             # doc's own confirmation criterion: +/-15 MeV
MAX_M    = 6                # highest angular harmonic considered
MAX_Q    = 4               # max total quanta in a configuration

# GWT angular step at d=3 (the ONE scale we will randomize in the null)
def omega0_gwt():
    d = 3
    return M_P * (2*d*PI) / (2*d*PI**2 + 1)   # ~293.8 MeV

# ----------------------------------------------------------------------
# Observed PDG N* and Delta* states, 1200-3000 MeV (3- and 4-star, approx).
# Pulled to represent real spectral density, NOT cherry-picked to GWT.
# ----------------------------------------------------------------------
OBSERVED = sorted(set([
    1232, 1440, 1520, 1535, 1600, 1620, 1650, 1675, 1680, 1700, 1710, 1720,
    1750, 1860, 1875, 1880, 1895, 1900, 1905, 1910, 1920, 1930, 1940, 1950,
    1990, 2000, 2040, 2060, 2100, 2120, 2150, 2190, 2200, 2220, 2250, 2300,
    2350, 2390, 2400, 2420, 2570, 2600, 2700, 2750, 2950,
]))
OBSERVED = [m for m in OBSERVED if WINDOW[0] <= m <= WINDOW[1]]

# ----------------------------------------------------------------------
# Build the full GWT prediction set from a given omega_0.
# Everything below scales linearly with omega_0 (verified against doc).
# ----------------------------------------------------------------------
def predictions(omega0):
    preds = set()

    # angular mode energies omega_m = omega0 * (m + (m-1)/pi)
    wm = {m: omega0 * (m + (m - 1) / PI) for m in range(1, MAX_M + 1)}

    # enumerate configs: choose counts k_m >= 0 over modes m=1..MAX_M
    # with total quanta sum(k_m) between 1 and MAX_Q
    def recurse(m, remaining, chosen):
        if m > MAX_M:
            q = sum(chosen.values())
            if q == 0:
                return
            energy = sum(k * wm[mm] for mm, k in chosen.items())
            n_distinct = sum(1 for k in chosen.values() if k > 0)
            shift = (omega0 / (4 * PI)) if n_distinct >= 2 else 0.0
            mass = M_P + energy - shift
            if WINDOW[0] <= mass <= WINDOW[1]:
                preds.add(round(mass, 1))
            return
        for k in range(0, remaining + 1):
            chosen[m] = k
            recurse(m + 1, remaining - k, chosen)
        chosen[m] = 0

    recurse(1, MAX_Q, {})

    # P-wave centroid: omega0 * (pi^2+1)/(pi+1)
    pw = M_P + omega0 * (PI**2 + 1) / (PI + 1)
    if WINDOW[0] <= pw <= WINDOW[1]:
        preds.add(round(pw, 1))

    # breathing (Roper): sqrt(3)*omega0*(1-alpha_s/8)
    br = M_P + np.sqrt(3) * omega0 * (1 - ALPHA_S / 8)
    if WINDOW[0] <= br <= WINDOW[1]:
        preds.add(round(br, 1))

    return np.array(sorted(preds))

# ----------------------------------------------------------------------
# How many observed states are matched (within tol) by a prediction set?
# Each observed state may be claimed by at most one prediction; each
# prediction may claim at most one observed state (greedy, nearest-first).
# ----------------------------------------------------------------------
def count_matches(preds, observed=OBSERVED, tol=TOL_MEV):
    obs = list(observed)
    used_pred = set()
    matched = 0
    # for each observed, take nearest unused prediction within tol
    for o in obs:
        best = None; bestd = tol + 1
        for i, p in enumerate(preds):
            if i in used_pred:
                continue
            d = abs(p - o)
            if d <= tol and d < bestd:
                best = i; bestd = d
        if best is not None:
            used_pred.add(best); matched += 1
    return matched

# ----------------------------------------------------------------------
# Run GWT
# ----------------------------------------------------------------------
w0 = omega0_gwt()
gwt_preds = predictions(w0)
gwt_match = count_matches(gwt_preds)

# analytic coverage baseline: fraction of window within tol of SOME prediction
def coverage(preds, tol=TOL_MEV):
    lo, hi = WINDOW
    intervals = sorted([(max(lo, p - tol), min(hi, p + tol)) for p in preds])
    merged = []
    for a, b in intervals:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    covered = sum(b - a for a, b in merged)
    return covered / (hi - lo)

cov = coverage(gwt_preds)
expected_random_matches = cov * len(OBSERVED)

print("="*64)
print("GWT BARYON LADDER  --  OVERFITTING NULL TEST")
print("="*64)
print(f"window               : {WINDOW[0]:.0f}-{WINDOW[1]:.0f} MeV")
print(f"observed N*/Delta*   : {len(OBSERVED)} states "
      f"(avg spacing {np.mean(np.diff(OBSERVED)):.0f} MeV)")
print(f"match tolerance      : +/-{TOL_MEV:.0f} MeV")
print(f"omega_0 (GWT, d=3)   : {w0:.2f} MeV")
print(f"# GWT predictions    : {len(gwt_preds)} distinct masses in window")
print(f"window coverage      : {cov*100:.1f}%  (fraction of window within tol "
      f"of a prediction)")
print("-"*64)
print(f"GWT states matched   : {gwt_match} / {len(OBSERVED)}")
print(f"expected by chance   : {expected_random_matches:.1f} "
      f"(coverage x #observed, if observed were random)")
print("-"*64)

# ----------------------------------------------------------------------
# Monte Carlo null: randomize omega_0, density-matched
# ----------------------------------------------------------------------
N_TRIALS = 20000
gwt_n_pred = len(gwt_preds)
null_matches = []
null_matches_density_matched = []
for _ in range(N_TRIALS):
    w = rng.uniform(0.5 * w0, 1.6 * w0)   # random step, wide band
    p = predictions(w)
    if len(p) == 0:
        continue
    mm = count_matches(p)
    null_matches.append(mm)
    # density-matched subset: only randoms with comparable # predictions
    if 0.7 * gwt_n_pred <= len(p) <= 1.3 * gwt_n_pred:
        null_matches_density_matched.append(mm)

null_matches = np.array(null_matches)
ndm = np.array(null_matches_density_matched)

p_all = np.mean(null_matches >= gwt_match)
p_dm  = np.mean(ndm >= gwt_match) if len(ndm) else float('nan')

print(f"MONTE CARLO NULL ({len(null_matches)} random ladders, "
      f"omega_0 ~ U[{0.5*w0:.0f},{1.6*w0:.0f}] MeV)")
print(f"  random match count : mean {null_matches.mean():.1f}, "
      f"median {np.median(null_matches):.0f}, "
      f"95th pct {np.percentile(null_matches,95):.0f}, "
      f"max {null_matches.max()}")
print(f"  p(random >= GWT)   : {p_all:.4f}   (all random ladders)")
print(f"  density-matched n  : {len(ndm)} ladders with ~same # predictions")
print(f"  p(random >= GWT)   : {p_dm:.4f}   (density-matched)")
print("="*64)
