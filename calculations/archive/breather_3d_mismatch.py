"""
3D Breather Width-Mismatch Simulation — GPU Accelerated
=========================================================
Key question: When two breathers have the SAME mode but DIFFERENT
spatial widths, how is the interaction energy suppressed?

This directly tests the phase extension hypothesis:
  - Equal-width breathers = homonuclear (N2-like)
  - Mismatched-width breathers = heteronuclear (BF-like, CO-like)

The breather width is controlled by a scaling factor w_scale:
  profile(r) = (4/pi) * arctan(eta / (w * cosh(eta * r / w_scale)))

w_scale > 1 = wider breather (lower Z_eff, like B in BF)
w_scale < 1 = narrower breather (higher Z_eff, like F in BF)
w_scale = 1 = standard (reference)

For BF analog: w_scale ratio = Z_F/Z_B = 5.10/2.42 = 2.11
For CO analog: w_scale ratio = Z_O/Z_C = 4.45/3.14 = 1.42
For CN analog: w_scale ratio = Z_N/Z_C = 3.83/3.14 = 1.22

We measure dE_int(R) for each width ratio and compare.
"""

import cupy as cp
import numpy as np
import time as timer
import sys

pi = float(np.pi)
V_0 = 1.0 / pi**2
N_br = 24
gamma_br = pi / (N_br + 1)

# Progress tracking
total_runs = 0
completed_runs = 0
start_time = None


def progress(desc=""):
    global completed_runs
    completed_runs += 1
    elapsed = timer.time() - start_time
    rate = completed_runs / elapsed if elapsed > 0 else 0
    remaining = (total_runs - completed_runs) / rate if rate > 0 else 0
    pct = 100 * completed_runs / total_runs
    mins_left = remaining / 60
    print(f"  [PROGRESS] {completed_runs}/{total_runs} ({pct:.0f}%) | "
          f"Elapsed: {elapsed/60:.1f}min | ETA: {mins_left:.1f}min | {desc}",
          flush=True)


def run_mismatch_gpu(N_grid, L, n_mode, R_sep, w_scale1=1.0, w_scale2=1.0,
                     sign2=+1, l_ang=0, n_steps=3000, dt_factor=0.15):
    """Run 3D breather interaction with width mismatch.

    Both breathers use the SAME mode number n_mode, but different
    spatial widths controlled by w_scale1 and w_scale2.

    w_scale > 1 means wider (more diffuse, like lower Z_eff)
    w_scale < 1 means narrower (more compact, like higher Z_eff)
    """
    dx = L / N_grid
    dt = dt_factor * dx

    coords = cp.linspace(-L/2, L/2, N_grid, dtype=cp.float32)
    X, Y, Z = cp.meshgrid(coords, coords, coords, indexing='ij')

    w = float(np.cos(n_mode * gamma_br))
    eta = float(np.sin(n_mode * gamma_br))

    def breather_3d(X, Y, Z, center_x, amplitude, w_scale, l_ang_val):
        rx = X - center_x
        r = cp.sqrt(rx**2 + Y**2 + Z**2)
        r = cp.maximum(r, dx/2)
        # Width scaling: r -> r/w_scale makes the breather wider
        r_scaled = r / w_scale
        profile = amplitude * (4.0/pi) * cp.arctan(eta / (w * cp.cosh(eta * r_scaled)))
        if l_ang_val == 1:
            cos_theta = rx / r
            profile = profile * cos_theta
        return profile

    def compute_energy(phi):
        gx = cp.diff(phi, axis=0)
        gy = cp.diff(phi, axis=1)
        gz = cp.diff(phi, axis=2)
        E_grad = 0.5 * (cp.sum(gx**2) + cp.sum(gy**2) + cp.sum(gz**2)) / dx * dx**2
        E_pot = V_0 * cp.sum(1 - cp.cos(pi * phi)) * dx**3
        return float(E_grad + E_pot)

    def compute_force(phi):
        force = cp.zeros_like(phi)
        force[1:-1,:,:] += phi[:-2,:,:] + phi[2:,:,:] - 2*phi[1:-1,:,:]
        force[:,1:-1,:] += phi[:,:-2,:] + phi[:,2:,:] - 2*phi[:,1:-1,:]
        force[:,:,1:-1] += phi[:,:,:-2] + phi[:,:,2:] - 2*phi[:,:,1:-1]
        force /= dx**2
        force -= V_0 * pi * cp.sin(pi * phi)
        return force

    def evolve_and_measure(phi_init, n_steps):
        phi = phi_init.copy()
        dphi = cp.zeros_like(phi)
        force = compute_force(phi)
        dphi += 0.5 * dt * force
        E_samples = []
        for step in range(n_steps):
            phi += dt * dphi
            force = compute_force(phi)
            dphi += dt * force
            if step % 50 == 0 and step > n_steps//4:
                dphi_half = dphi - 0.5*dt*force
                E_kin = 0.5 * float(cp.sum(dphi_half**2)) * dx**3
                E_pg = compute_energy(phi)
                E_samples.append(E_kin + E_pg)
        return np.mean(E_samples) if E_samples else compute_energy(phi_init)

    # Single breather 1 (with w_scale1)
    phi1 = breather_3d(X, Y, Z, 0.0, 1.0, w_scale1, l_ang)
    E1 = evolve_and_measure(phi1, n_steps)

    # Single breather 2 (with w_scale2)
    phi2 = breather_3d(X, Y, Z, 0.0, sign2, w_scale2, l_ang)
    E2 = evolve_and_measure(phi2, n_steps)

    # Two-breather system
    phi_two = (breather_3d(X, Y, Z, -R_sep/2, 1.0, w_scale1, l_ang) +
               breather_3d(X, Y, Z, +R_sep/2, sign2, w_scale2, l_ang))
    E12 = evolve_and_measure(phi_two, n_steps)

    E_int = E12 - E1 - E2
    cp.cuda.Stream.null.synchronize()
    return E1, E2, E12, E_int


# =============================================================================
# SIMULATION PLAN
# =============================================================================
# Width ratios to test (Z_ratio analogs):
#   1.0  = homonuclear (N2 reference)
#   1.22 = CN analog
#   1.42 = CO analog
#   2.11 = BF analog
#   3.0  = extreme mismatch (control)
#
# For each ratio, scan R from 1.5 to 10 with step 0.5
# Mode n=3 p-wave (the oscillatory case from previous sim)
# Both +- and ++ phases
#
# Width convention: breather1 has w_scale=sqrt(ratio), breather2 has w_scale=1/sqrt(ratio)
# This keeps the geometric mean width = 1 (same as homonuclear reference)

width_ratios = [1.0, 1.22, 1.42, 2.11, 3.0]
R_values = np.arange(1.5, 10.1, 0.5)  # 18 points
n_mode = 3

# Count runs: for each ratio, each R, we do +- and ++ = 2 runs per R
# Plus the reference single-breather measurements
# Each ratio × 18 R × 2 phases × 3 evolutions = 108 runs per ratio
# But single breather E1, E2 only need to be computed once per ratio
# Simplified: each (ratio, R, sign) = 1 run call = 3 evolutions
runs_per_ratio = len(R_values) * 2
total_runs = len(width_ratios) * runs_per_ratio

print("=" * 80)
print("  3D BREATHER WIDTH-MISMATCH SIMULATION")
print("  Testing: does spatial mismatch suppress interaction?")
print("=" * 80)

# Warm up GPU
_ = cp.zeros(10)
cp.cuda.Stream.null.synchronize()

N = 60
L = 20.0
n_steps = 3000

print(f"  Grid: {N}x{N}x{N} = {N**3} points")
print(f"  Box: L={L}, dx={L/N:.3f}")
print(f"  Steps: {n_steps}")
print(f"  Mode: n={n_mode}, p-wave (l=1)")
print(f"  Width ratios: {width_ratios}")
print(f"  R range: {R_values[0]:.1f} to {R_values[-1]:.1f}, step 0.5")
print(f"  Total runs: {total_runs}")
print(f"  Each run = 3 evolutions x {n_steps} steps on {N}^3 grid")
sys.stdout.flush()

start_time = timer.time()

# =============================================================================
# SECTION 1: p-wave (l=1) with width mismatch — dense R scan
# =============================================================================
print(f"\n{'='*80}")
print(f"  p-WAVE WIDTH MISMATCH SCAN")
print(f"  mode n={n_mode}, l=1")
print(f"{'='*80}", flush=True)

all_results = {}

for w_ratio in width_ratios:
    # Width scaling: geometric mean preserved at 1.0
    # w1 = sqrt(ratio), w2 = 1/sqrt(ratio)
    # This means breather1 is wider, breather2 is narrower
    w1 = np.sqrt(w_ratio)
    w2 = 1.0 / np.sqrt(w_ratio)

    # The Z-analog coherence factor
    # base = 2*sqrt(Z1*Z2)/(Z1+Z2) = 2*sqrt(w_ratio)/(1+w_ratio)
    # (since Z ~ 1/w_scale, higher Z = narrower)
    base = 2*np.sqrt(w_ratio)/(1+w_ratio)

    print(f"\n  --- Width ratio = {w_ratio:.2f} (w1={w1:.3f}, w2={w2:.3f}, "
          f"base={base:.4f}, base^2={base**2:.4f}, base^5={base**5:.4f}) ---")
    print(f"  {'R':>6} {'dE(+-)':>10} {'dE(++)':>10}", flush=True)

    results = []
    for R in R_values:
        t0 = timer.time()

        E1, E2, E12_pm, dE_pm = run_mismatch_gpu(
            N, L, n_mode, R, w_scale1=w1, w_scale2=w2,
            sign2=-1, l_ang=1, n_steps=n_steps)
        progress(f"ratio={w_ratio:.2f} R={R:.1f} +-")

        _, _, E12_pp, dE_pp = run_mismatch_gpu(
            N, L, n_mode, R, w_scale1=w1, w_scale2=w2,
            sign2=+1, l_ang=1, n_steps=n_steps)
        progress(f"ratio={w_ratio:.2f} R={R:.1f} ++")

        elapsed = timer.time() - t0
        print(f"  {R:6.1f} {dE_pm:+10.4f} {dE_pp:+10.4f}  ({elapsed:.1f}s)", flush=True)

        results.append((R, dE_pm, dE_pp))

    all_results[w_ratio] = results

# =============================================================================
# SECTION 2: s-wave comparison (should show monotonic for all ratios)
# =============================================================================
print(f"\n{'='*80}")
print(f"  s-WAVE WIDTH MISMATCH (control — should be monotonic)")
print(f"  mode n={n_mode}, l=0")
print(f"{'='*80}", flush=True)

# Only test the extreme ratios for s-wave (to save time)
# Update total_runs for progress tracking
s_wave_ratios = [1.0, 2.11]
s_wave_R = np.arange(2.0, 10.1, 1.0)  # coarser
extra_runs = len(s_wave_ratios) * len(s_wave_R) * 2
total_runs += extra_runs

for w_ratio in s_wave_ratios:
    w1 = np.sqrt(w_ratio)
    w2 = 1.0 / np.sqrt(w_ratio)
    base = 2*np.sqrt(w_ratio)/(1+w_ratio)

    print(f"\n  --- s-wave, ratio={w_ratio:.2f} (base={base:.4f}) ---")
    print(f"  {'R':>6} {'dE(+-)':>10} {'dE(++)':>10}", flush=True)

    for R in s_wave_R:
        _, _, _, dE_pm = run_mismatch_gpu(
            N, L, n_mode, R, w_scale1=w1, w_scale2=w2,
            sign2=-1, l_ang=0, n_steps=n_steps)
        progress(f"s-wave ratio={w_ratio:.2f} R={R:.1f} +-")

        _, _, _, dE_pp = run_mismatch_gpu(
            N, L, n_mode, R, w_scale1=w1, w_scale2=w2,
            sign2=+1, l_ang=0, n_steps=n_steps)
        progress(f"s-wave ratio={w_ratio:.2f} R={R:.1f} ++")

        print(f"  {R:6.1f} {dE_pm:+10.4f} {dE_pp:+10.4f}", flush=True)


# =============================================================================
# ANALYSIS: Compare interaction strength vs width ratio
# =============================================================================
print(f"\n{'='*80}")
print(f"  ANALYSIS: Interaction suppression vs width ratio")
print(f"{'='*80}")

# For each ratio, find the peak |dE| and compare to homonuclear
ref_results = all_results[1.0]
ref_peak = max(abs(dE_pp) for _, _, dE_pp in ref_results)

print(f"\n  Homonuclear (ratio=1.0) peak |dE(++)| = {ref_peak:.4f}")
print(f"\n  {'ratio':>6} {'peak_dE':>10} {'suppression':>12} {'base^2':>8} {'base^5':>8}")
print(f"  {'-'*50}")

for w_ratio in width_ratios:
    results = all_results[w_ratio]
    peak = max(abs(dE_pp) for _, _, dE_pp in results)
    suppression = peak / ref_peak if ref_peak > 0 else 0
    base = 2*np.sqrt(w_ratio)/(1+w_ratio)
    print(f"  {w_ratio:6.2f} {peak:10.4f} {suppression:12.4f} {base**2:8.4f} {base**5:8.4f}")


# =============================================================================
# ANALYSIS: Phase shift detection
# =============================================================================
print(f"\n{'='*80}")
print(f"  ANALYSIS: Does width mismatch shift the oscillation period?")
print(f"{'='*80}")

for w_ratio in width_ratios:
    results = all_results[w_ratio]
    # Find zero crossings of dE(++) (sign changes)
    crossings = []
    for i in range(1, len(results)):
        R_prev, _, dE_prev = results[i-1]
        R_curr, _, dE_curr = results[i]
        if dE_prev * dE_curr < 0:
            # Linear interpolation for crossing
            R_cross = R_prev + (R_curr - R_prev) * abs(dE_prev) / (abs(dE_prev) + abs(dE_curr))
            crossings.append(R_cross)

    base = 2*np.sqrt(w_ratio)/(1+w_ratio)
    cross_str = ', '.join(f'{c:.2f}' for c in crossings)
    print(f"  ratio={w_ratio:.2f}: zero crossings at R = [{cross_str}]")
    if len(crossings) >= 2:
        period = 2 * (crossings[1] - crossings[0])
        print(f"         half-period = {crossings[1]-crossings[0]:.2f}, "
              f"full period ~ {period:.2f}")


total_elapsed = timer.time() - start_time
print(f"\n{'='*80}")
print(f"  DONE — Total time: {total_elapsed/60:.1f} minutes")
print(f"{'='*80}")
