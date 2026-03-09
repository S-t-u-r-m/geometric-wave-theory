"""
Quantitative fit: simulation zero crossings vs phase extension formula
======================================================================
Compare the GPU breather simulation's width-mismatch results against
the analytical phase extension prediction:

  R_cross(ratio) = R_cross_ref * (1/base^delta)

where base = 2*sqrt(ratio)/(1+ratio)  and delta is the fit parameter.

The GWT bond formula uses delta = d-1 = 2.
"""

import numpy as np
from scipy.optimize import curve_fit

# =====================================================================
# SIMULATION DATA (from breather_3d_mismatch_output.log)
# =====================================================================

# p-wave dE(++) data for each width ratio
sim_data = {
    1.00: [
        (1.5, +9.8398), (2.0, +6.9232), (2.5, +4.0025), (3.0, +1.3509),
        (3.5, -0.8922), (4.0, -2.8292), (4.5, -4.4751), (5.0, -5.6980),
        (5.5, -6.5239), (6.0, -7.1374), (6.5, -7.5545), (7.0, -7.6534),
        (7.5, -7.4985), (8.0, -7.2385), (8.5, -6.8130), (9.0, -6.0732),
        (9.5, -5.0823), (10.0, -3.9138),
    ],
    1.22: [
        (1.5, +9.6729), (2.0, +6.8057), (2.5, +3.9449), (3.0, +1.3597),
        (3.5, -0.8151), (4.0, -2.6854), (4.5, -4.2685), (5.0, -5.4334),
        (5.5, -6.2057), (6.0, -6.7692), (6.5, -7.1382), (7.0, -7.1893),
        (7.5, -6.9853), (8.0, -6.6737), (8.5, -6.1932), (9.0, -5.3948),
        (9.5, -4.3420), (10.0, -3.1094),
    ],
    1.42: [
        (1.5, +9.3437), (2.0, +6.5800), (2.5, +3.8449), (3.0, +1.3983),
        (3.5, -0.6336), (4.0, -2.3640), (4.5, -3.8148), (5.0, -4.8565),
        (5.5, -5.5141), (6.0, -5.9694), (6.5, -6.2336), (7.0, -6.1797),
        (7.5, -5.8676), (8.0, -5.4425), (8.5, -4.8414), (9.0, -3.9151),
        (9.5, -2.7281), (10.0, -1.3576),
    ],
    2.11: [
        (1.5, +8.0504), (2.0, +5.8585), (2.5, +3.8012), (3.0, +2.0949),
        (3.5, +0.8268), (4.0, -0.1457), (4.5, -0.8612), (5.0, -1.1950),
        (5.5, -1.1694), (6.0, -0.9572), (6.5, -0.5568), (7.0, +0.1733),
        (7.5, +1.1839), (8.0, +2.3350), (8.5, +3.6887), (9.0, +5.3863),
        (9.5, +7.3463), (10.0, +9.4658),
    ],
    3.00: [
        (1.5, +7.0438), (2.0, +6.0845), (2.5, +5.4257), (3.0, +5.2231),
        (3.5, +5.5211), (4.0, +6.1494), (4.5, +7.0549), (5.0, +8.3558),
        (5.5, +10.0244), (6.0, +11.8812), (6.5, +13.9198), (7.0, +16.2660),
        (7.5, +18.8399), (8.0, +21.4526), (8.5, +24.0977), (9.0, +26.8325),
        (9.5, +29.4815), (10.0, +31.8444),
    ],
}

# Also grab dE(+-) for ratio=2.11 to check anti-bonding crossing
sim_data_pm = {
    2.11: [
        (1.5, -8.8902), (2.0, -6.1944), (2.5, -3.7176), (3.0, -1.4028),
        (3.5, +0.7799), (4.0, +2.6784), (4.5, +4.2779), (5.0, +5.7641),
        (5.5, +7.1946), (6.0, +8.4673), (6.5, +9.6416), (7.0, +10.9056),
        (7.5, +12.2441), (8.0, +13.5358), (8.5, +14.8516), (9.0, +16.3359),
        (9.5, +17.9008), (10.0, +19.4231),
    ],
}


def find_zero_crossings(data):
    """Find zero crossings by linear interpolation."""
    crossings = []
    for i in range(1, len(data)):
        R_prev, dE_prev = data[i-1]
        R_curr, dE_curr = data[i]
        if dE_prev * dE_curr < 0:
            R_cross = R_prev + (R_curr - R_prev) * abs(dE_prev) / (abs(dE_prev) + abs(dE_curr))
            crossings.append(R_cross)
    return crossings


def base_factor(ratio):
    """Coherence base = 2*sqrt(ratio)/(1+ratio)"""
    return 2 * np.sqrt(ratio) / (1 + ratio)


# =====================================================================
# 1. EXTRACT ZERO CROSSINGS
# =====================================================================
print("=" * 70)
print("  ZERO CROSSING ANALYSIS")
print("=" * 70)

crossings = {}
for ratio in [1.0, 1.22, 1.42, 2.11, 3.0]:
    c = find_zero_crossings(sim_data[ratio])
    crossings[ratio] = c
    base = base_factor(ratio)
    c_str = ', '.join(f'{x:.3f}' for x in c)
    print(f"  ratio={ratio:.2f}  base={base:.4f}  crossings: [{c_str}]")

# Reference crossing
R0 = crossings[1.0][0]  # = 3.30
print(f"\n  Reference crossing R0 = {R0:.3f}")

# =====================================================================
# 2. FIT: R_cross = R0 * (1/base^delta)
# =====================================================================
print(f"\n{'='*70}")
print("  FIT: R_cross = R0 / base^delta")
print("=" * 70)

# Use ratios that have first crossings: 1.0, 1.22, 1.42, 2.11
fit_ratios = [1.0, 1.22, 1.42, 2.11]
fit_crossings = [crossings[r][0] for r in fit_ratios]
fit_bases = [base_factor(r) for r in fit_ratios]

# The crossing shift relative to reference
shifts = [c / R0 for c in fit_crossings]

print(f"\n  {'ratio':>6} {'base':>8} {'R_cross':>8} {'shift':>8} {'ln(shift)':>10} {'-ln(base)':>10}")
for ratio, base, R_c, shift in zip(fit_ratios, fit_bases, fit_crossings, shifts):
    ln_shift = np.log(shift) if shift > 0 else 0
    ln_base = -np.log(base) if base > 0 else 0
    print(f"  {ratio:6.2f} {base:8.4f} {R_c:8.3f} {shift:8.4f} {ln_shift:10.6f} {ln_base:10.6f}")

# Linear fit: ln(shift) = delta * (-ln(base))
# i.e., ln(R_cross/R0) = -delta * ln(base)
ln_shifts = np.array([np.log(s) for s in shifts])
neg_ln_bases = np.array([-np.log(b) for b in fit_bases])

# Exclude ratio=1.0 (trivially 0,0) for meaningful fit
mask = np.array(fit_ratios) != 1.0
if np.sum(mask) >= 2:
    # Weighted least squares through origin: ln(shift) = delta * (-ln(base))
    delta_fit = np.sum(ln_shifts[mask] * neg_ln_bases[mask]) / np.sum(neg_ln_bases[mask]**2)
    residuals = ln_shifts[mask] - delta_fit * neg_ln_bases[mask]
    rmse = np.sqrt(np.mean(residuals**2))
    print(f"\n  Linear fit (through origin): delta = {delta_fit:.3f}")
    print(f"  RMSE of log fit: {rmse:.6f}")
else:
    delta_fit = 2.0
    print(f"\n  Not enough points for fit, using delta = 2.0")

# =====================================================================
# 3. COMPARE PREDICTIONS FOR VARIOUS DELTA VALUES
# =====================================================================
print(f"\n{'='*70}")
print("  PREDICTED vs SIMULATED CROSSINGS for delta = 1, 2, 3, fit")
print("=" * 70)

print(f"\n  {'ratio':>6} {'sim_R':>8} | {'d=1':>8} {'err%':>6} | {'d=2':>8} {'err%':>6} | "
      f"{'d=3':>8} {'err%':>6} | {'d=fit':>8} {'err%':>6}")
print(f"  {'-'*80}")

for ratio in fit_ratios:
    base = base_factor(ratio)
    R_sim = crossings[ratio][0]

    preds = {}
    for delta in [1, 2, 3, delta_fit]:
        R_pred = R0 / base**delta
        err = 100 * (R_pred - R_sim) / R_sim
        preds[delta] = (R_pred, err)

    d1, d2, d3, df = preds[1], preds[2], preds[3], preds[delta_fit]
    print(f"  {ratio:6.2f} {R_sim:8.3f} | {d1[0]:8.3f} {d1[1]:+5.1f}% | {d2[0]:8.3f} {d2[1]:+5.1f}% | "
          f"{d3[0]:8.3f} {d3[1]:+5.1f}% | {df[0]:8.3f} {df[1]:+5.1f}%")


# =====================================================================
# 4. BONDING WELL DEPTH ANALYSIS
# =====================================================================
print(f"\n{'='*70}")
print("  BONDING WELL DEPTH ANALYSIS")
print("=" * 70)

print(f"\n  Peak bonding = minimum dE(++) in bonding region (R > first crossing)")

for ratio in [1.0, 1.22, 1.42, 2.11]:
    data = sim_data[ratio]
    first_cross = crossings[ratio][0]
    # Find minimum dE(++) after first crossing
    bonding_pts = [(R, dE) for R, dE in data if R > first_cross]
    if bonding_pts:
        R_min, dE_min = min(bonding_pts, key=lambda x: x[1])
        base = base_factor(ratio)
        print(f"  ratio={ratio:.2f}: peak bonding dE={dE_min:.4f} at R={R_min:.1f}  "
              f"(base={base:.4f})")

# Amplitude suppression: compare peak depths
ref_depth = min(dE for R, dE in sim_data[1.0] if R > crossings[1.0][0])
print(f"\n  Reference depth (ratio=1.0): {ref_depth:.4f}")
print(f"\n  {'ratio':>6} {'depth':>8} {'suppress':>10} {'base^2':>8} {'base^3':>8} {'base^5':>8}")

for ratio in [1.0, 1.22, 1.42, 2.11]:
    data = sim_data[ratio]
    first_cross = crossings[ratio][0]
    bonding_pts = [(R, dE) for R, dE in data if R > first_cross]
    if bonding_pts:
        depth = min(dE for R, dE in bonding_pts)
        suppress = depth / ref_depth if ref_depth != 0 else 0
        base = base_factor(ratio)
        print(f"  {ratio:6.2f} {depth:8.4f} {suppress:10.4f} {base**2:8.4f} {base**3:8.4f} {base**5:8.4f}")


# =====================================================================
# 5. SECOND CROSSING (ratio=2.11) — period doubling?
# =====================================================================
print(f"\n{'='*70}")
print("  SECOND CROSSING ANALYSIS (ratio=2.11)")
print("=" * 70)

c211 = crossings[2.11]
if len(c211) >= 2:
    half_period = c211[1] - c211[0]
    full_period = 2 * half_period
    print(f"  First crossing:  R = {c211[0]:.3f}")
    print(f"  Second crossing: R = {c211[1]:.3f}")
    print(f"  Half-period: {half_period:.3f}")
    print(f"  Full period: {full_period:.3f}")

    # Compare to homonuclear: ratio=1.0 has crossing at 3.30 but no second crossing
    # in R=[1.5, 10] — the bonding well extends beyond R=10
    # So the mismatch CREATES a finite bonding region

    # The bonding well width for ratio=2.11
    well_width = c211[1] - c211[0]
    print(f"\n  Bonding well width: {well_width:.3f}")
    print(f"  (homonuclear well extends beyond R=10, effectively infinite in our box)")


# =====================================================================
# 6. OVERALL FIT QUALITY — RESCALED CURVES
# =====================================================================
print(f"\n{'='*70}")
print("  RESCALED CURVE COMPARISON")
print(f"  If phase extension is correct, dE(R) ~ f(R * base^delta)")
print("=" * 70)

# For each ratio, rescale R -> R * base^delta_fit and compare to reference
ref_R = np.array([R for R, _ in sim_data[1.0]])
ref_dE = np.array([dE for _, dE in sim_data[1.0]])

print(f"\n  Using delta = {delta_fit:.3f}")
print(f"\n  For each ratio, we rescale R -> R * base^delta and interpolate")
print(f"  to compare against the reference (ratio=1.0) curve.\n")

for ratio in [1.22, 1.42, 2.11]:
    base = base_factor(ratio)
    data = sim_data[ratio]
    R_raw = np.array([R for R, _ in data])
    dE_raw = np.array([dE for _, dE in data])

    # Rescale R
    R_rescaled = R_raw * base**delta_fit

    # Interpolate reference at rescaled R points
    # (only where within reference range)
    mask_interp = (R_rescaled >= ref_R[0]) & (R_rescaled <= ref_R[-1])
    if np.sum(mask_interp) > 0:
        ref_interp = np.interp(R_rescaled[mask_interp], ref_R, ref_dE)
        # Also need amplitude rescaling?
        # Compare raw dE vs reference at rescaled R
        dE_at_rescaled = dE_raw[mask_interp]

        # RMS difference (no amplitude correction)
        rms_no_amp = np.sqrt(np.mean((dE_at_rescaled - ref_interp)**2))

        # Best amplitude scale factor
        amp_scale = np.sum(dE_at_rescaled * ref_interp) / np.sum(ref_interp**2)
        rms_with_amp = np.sqrt(np.mean((dE_at_rescaled - amp_scale * ref_interp)**2))

        print(f"  ratio={ratio:.2f} (base={base:.4f}):")
        print(f"    Phase-only rescaling: RMS = {rms_no_amp:.4f}")
        print(f"    Phase + amplitude (scale={amp_scale:.4f}): RMS = {rms_with_amp:.4f}")
        print(f"    Amplitude suppression: {amp_scale:.4f}")


# =====================================================================
# 7. BEST-FIT DELTA WITH AMPLITUDE
# =====================================================================
print(f"\n{'='*70}")
print("  SCAN: best delta (phase) and amplitude suppression")
print("=" * 70)

print(f"\n  {'delta':>6} {'amp_1.22':>10} {'amp_1.42':>10} {'amp_2.11':>10} "
      f"{'rms_1.22':>10} {'rms_1.42':>10} {'rms_2.11':>10} {'rms_total':>10}")

best_rms = 1e10
best_delta = 0

for delta_test in np.arange(0.0, 6.1, 0.25):
    rms_vals = []
    amp_vals = []

    for ratio in [1.22, 1.42, 2.11]:
        base = base_factor(ratio)
        data = sim_data[ratio]
        R_raw = np.array([R for R, _ in data])
        dE_raw = np.array([dE for _, dE in data])

        R_rescaled = R_raw * base**delta_test
        mask_interp = (R_rescaled >= ref_R[0]) & (R_rescaled <= ref_R[-1])

        if np.sum(mask_interp) >= 3:
            ref_interp = np.interp(R_rescaled[mask_interp], ref_R, ref_dE)
            dE_at = dE_raw[mask_interp]
            amp = np.sum(dE_at * ref_interp) / np.sum(ref_interp**2)
            rms = np.sqrt(np.mean((dE_at - amp * ref_interp)**2))
        else:
            amp = 1.0
            rms = 999.0

        amp_vals.append(amp)
        rms_vals.append(rms)

    rms_total = np.sqrt(np.mean(np.array(rms_vals)**2))

    if rms_total < best_rms:
        best_rms = rms_total
        best_delta = delta_test

    if delta_test in [0, 1, 2, 3, 4, 5] or abs(delta_test - delta_fit) < 0.13:
        print(f"  {delta_test:6.2f} {amp_vals[0]:10.4f} {amp_vals[1]:10.4f} {amp_vals[2]:10.4f} "
              f"{rms_vals[0]:10.4f} {rms_vals[1]:10.4f} {rms_vals[2]:10.4f} {rms_total:10.4f}")

print(f"\n  Best delta = {best_delta:.2f} (RMS = {best_rms:.4f})")

# Show all deltas near the best
print(f"\n  Fine scan around best delta:")
for delta_test in np.arange(max(0, best_delta - 1), best_delta + 1.1, 0.25):
    rms_vals = []
    amp_vals = []

    for ratio in [1.22, 1.42, 2.11]:
        base = base_factor(ratio)
        data = sim_data[ratio]
        R_raw = np.array([R for R, _ in data])
        dE_raw = np.array([dE for _, dE in data])

        R_rescaled = R_raw * base**delta_test
        mask_interp = (R_rescaled >= ref_R[0]) & (R_rescaled <= ref_R[-1])

        if np.sum(mask_interp) >= 3:
            ref_interp = np.interp(R_rescaled[mask_interp], ref_R, ref_dE)
            dE_at = dE_raw[mask_interp]
            amp = np.sum(dE_at * ref_interp) / np.sum(ref_interp**2)
            rms = np.sqrt(np.mean((dE_at - amp * ref_interp)**2))
        else:
            amp = 1.0
            rms = 999.0

        amp_vals.append(amp)
        rms_vals.append(rms)

    rms_total = np.sqrt(np.mean(np.array(rms_vals)**2))
    print(f"  {delta_test:6.2f} {amp_vals[0]:10.4f} {amp_vals[1]:10.4f} {amp_vals[2]:10.4f} "
          f"{rms_vals[0]:10.4f} {rms_vals[1]:10.4f} {rms_vals[2]:10.4f} {rms_total:10.4f}")


# =====================================================================
# 8. SUMMARY
# =====================================================================
print(f"\n{'='*70}")
print("  SUMMARY")
print("=" * 70)

print(f"""
  Phase extension hypothesis: R_eff = R / base^delta
    base = 2*sqrt(Z_ratio) / (1 + Z_ratio)

  From zero-crossing fit:  delta = {delta_fit:.3f}
  GWT formula uses:        delta = d-1 = 2

  Zero crossing predictions (delta={delta_fit:.3f}):
""")

for ratio in fit_ratios:
    base = base_factor(ratio)
    R_sim = crossings[ratio][0]
    R_pred = R0 / base**delta_fit
    err = 100 * (R_pred - R_sim) / R_sim
    print(f"    ratio={ratio:.2f}: sim={R_sim:.3f}  pred={R_pred:.3f}  err={err:+.1f}%")

print(f"""
  Key findings:
  1. Width mismatch SHIFTS zero crossings (phase effect confirmed)
  2. Width mismatch SUPPRESSES bonding depth (amplitude effect)
  3. Large mismatch (ratio=2.11) creates FINITE bonding well
  4. Extreme mismatch (ratio=3.0) eliminates bonding entirely
""")
