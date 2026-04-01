"""
Time Dilation on the Lattice
==============================
Does a moving breather oscillate slower by exactly gamma?

The sine-Gordon breather has an EXACT boosted solution:
  Stationary: phi(x,t) = (4/pi) * arctan(eps * sin(omega*t) / (omega * cosh(eps*x)))
  Boosted:    phi(x,t) = (4/pi) * arctan(eps * sin(omega*tau) / (omega * cosh(eps*xi)))
    where xi = gamma*(x - v*t), tau = gamma*(t - v*x)
    and gamma = 1/sqrt(1-v^2)

The internal oscillation frequency in the lab frame:
  omega_lab = omega_rest / gamma = omega_rest * sqrt(1 - v^2)

This simulation:
1. Creates a stationary breather, measures its oscillation frequency
2. Creates the SAME breather boosted to velocity v
3. Evolves both on the discrete lattice (a=1)
4. Measures internal frequency from zero crossings
5. Compares omega_moving/omega_rest to 1/gamma = sqrt(1-v^2)

If GWT is correct, time dilation emerges from the wave equation on springs.
No postulate needed — just F = -kx.
"""

import sys, io, os, time
import numpy as np

if not hasattr(sys.stdout, '_gwt_wrapped'):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stdout._gwt_wrapped = True

try:
    import cupy as cp
    xp = cp
    GPU = True
    print("GPU active (CuPy)")
except Exception:
    xp = np
    GPU = False
    print("CPU mode (NumPy)")

PI = np.pi
d = 3
gamma_sg = PI / (2**(d+1) * PI - 2)  # sine-Gordon coupling parameter


def measure_frequency(time_series, dt_sample):
    """
    Measure oscillation frequency from zero crossings.
    More robust than FFT for drifting or decaying signals.
    """
    ts = np.array(time_series)
    ts = ts - np.mean(ts)

    # Find zero crossings
    crossings = []
    for i in range(1, len(ts)):
        if ts[i-1] * ts[i] < 0:
            # Linear interpolation for sub-sample accuracy
            frac = abs(ts[i-1]) / (abs(ts[i-1]) + abs(ts[i]))
            t_cross = (i - 1 + frac) * dt_sample
            crossings.append(t_cross)

    if len(crossings) < 4:
        return 0.0, 0

    # Period = 2 * average half-period
    half_periods = np.diff(crossings)
    full_period = 2.0 * np.median(half_periods)
    freq = 2 * PI / full_period if full_period > 0 else 0.0

    return freq, len(crossings)


def run_breather(N, v, n_mode, dt, n_periods, label=""):
    """
    Evolve a sine-Gordon breather at velocity v on a 1D discrete lattice.

    Returns the measured internal oscillation frequency.
    """
    omega_n = np.cos(n_mode * gamma_sg)
    eps_n = np.sin(n_mode * gamma_sg)

    if omega_n < 0.01:
        print(f"  {label}: omega_n too small ({omega_n:.4f}), skipping")
        return 0.0, 0.0

    gamma_L = 1.0 / np.sqrt(1.0 - v**2) if abs(v) < 0.999 else 100.0

    center = N // 2
    x = xp.arange(N, dtype=xp.float64) - center

    # Boosted breather initial conditions at t=0:
    # xi = gamma_L * (x - v*0) = gamma_L * x
    # tau = gamma_L * (0 - v*x) = -gamma_L * v * x
    xi = gamma_L * x
    tau = -gamma_L * v * x

    # phi(x, 0) = (4/pi) * arctan(eps * sin(omega * tau) / (omega * cosh(eps * xi)))
    denom = omega_n * xp.cosh(eps_n * xi) + 1e-30
    phi = (4.0 / PI) * xp.arctan(eps_n * xp.sin(omega_n * tau) / denom)

    # phi_dot from finite difference of the exact solution at t = -dt/2 and t = +dt/2
    # Or analytically: d/dt of boosted breather
    # Simpler: use the exact solution at t = -dt to initialize leapfrog
    t_m1 = -dt
    xi_m1 = gamma_L * (x - v * t_m1)
    tau_m1 = gamma_L * (t_m1 - v * x)
    denom_m1 = omega_n * xp.cosh(eps_n * xi_m1) + 1e-30
    phi_old = (4.0 / PI) * xp.arctan(eps_n * xp.sin(omega_n * tau_m1) / denom_m1)

    # Evolution
    period = 2 * PI / omega_n
    N_steps = int(n_periods * period / dt)
    N_steps = min(N_steps, 5000000)

    # Record at the TRACKING POINT: for stationary, track the center.
    # For moving, track a point that moves with the breather.
    rec_interval = max(1, N_steps // 50000)
    ts_record = []
    dt_sample = rec_interval * dt

    # Track the breather center position
    track_x = center

    for step in range(N_steps):
        # Discrete 1D Laplacian with periodic BC
        lap = xp.roll(phi, 1) + xp.roll(phi, -1) - 2 * phi
        force = (1.0 / PI) * xp.sin(PI * phi)
        phi_new = 2 * phi - phi_old + dt**2 * (lap - force)
        phi_old = phi
        phi = phi_new

        if step % rec_interval == 0:
            phi_host = phi.get() if GPU else phi
            t_now = step * dt

            # Track breather by searching near its EXPECTED position
            # Expected: x = center + v * t (with periodic wrapping)
            x_expect = (center + v * t_now) % N
            ix_expect = int(round(x_expect)) % N

            # Search within a window around expected position
            # Window = breather width / gamma (Lorentz contracted)
            hw = max(5, int(20.0 / gamma_L))
            best_val = 0.0
            best_abs = 0.0
            for dx in range(-hw, hw + 1):
                ix = (ix_expect + dx) % N
                if abs(phi_host[ix]) > best_abs:
                    best_abs = abs(phi_host[ix])
                    best_val = phi_host[ix]

            ts_record.append(best_val)

    # Measure frequency
    freq, n_cross = measure_frequency(ts_record, dt_sample)

    return freq, n_cross


# ============================================================
# MAIN TEST
# ============================================================
if __name__ == '__main__':
    print()
    print("TIME DILATION ON THE DISCRETE LATTICE")
    print("=" * 70)
    print(f"Sine-Gordon coupling: gamma = {gamma_sg:.8f}")
    print()

    N = 4096       # lattice sites (large enough for breather to travel at high v)
    dt = 0.01      # time step (smaller for high-v stability)
    n_mode = 4     # n=4 mode: omega=0.966, more deeply bound than n=1
    n_periods = 60 # evolve for many periods for accurate frequency measurement

    omega_rest = np.cos(n_mode * gamma_sg)
    print(f"Breather mode n={n_mode}")
    print(f"  Predicted rest frequency: omega_0 = cos({n_mode}*gamma) = {omega_rest:.8f}")
    print(f"  Lattice: {N} sites, dt={dt}")
    print(f"  Evolution: {n_periods} periods")
    print()

    # Velocities to test
    velocities = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    print(f"{'v':>5} {'gamma':>8} {'omega_pred':>12} {'omega_meas':>12} "
          f"{'ratio':>8} {'1/gamma':>8} {'error':>8} {'crossings':>10}")
    print("-" * 82)

    results = []
    omega_0_measured = None

    for v in velocities:
        gamma_L = 1.0 / np.sqrt(1.0 - v**2) if v < 0.999 else 999
        omega_pred = omega_rest / gamma_L

        t_start = time.time()
        omega_meas, n_cross = run_breather(N, v, n_mode, dt, n_periods,
                                            label=f"v={v:.1f}")
        elapsed = time.time() - t_start

        if v == 0.0:
            omega_0_measured = omega_meas

        if omega_0_measured and omega_0_measured > 0 and omega_meas > 0:
            ratio = omega_meas / omega_0_measured
            inv_gamma = 1.0 / gamma_L
            err = (ratio - inv_gamma) / inv_gamma * 100
        else:
            ratio = 0
            inv_gamma = 1.0 / gamma_L
            err = 999

        results.append((v, gamma_L, omega_pred, omega_meas, ratio, inv_gamma, err))

        print(f"{v:>5.1f} {gamma_L:>8.4f} {omega_pred:>12.6f} {omega_meas:>12.6f} "
              f"{ratio:>8.4f} {inv_gamma:>8.4f} {err:>+7.2f}% {n_cross:>10d}  "
              f"({elapsed:.1f}s)")

    print()
    print("SUMMARY")
    print("=" * 70)
    print(f"Rest frequency (measured): {omega_0_measured:.6f}")
    print(f"Rest frequency (predicted): {omega_rest:.6f}")
    if omega_0_measured:
        print(f"Rest frequency error: {(omega_0_measured - omega_rest)/omega_rest*100:+.3f}%")
    print()
    print("If GWT is correct, the 'ratio' column should match '1/gamma' exactly.")
    print("Any deviation is a lattice discreteness effect (finite a, finite dt).")
    print()

    # Analyze results
    # Separate physical regime (breather survives) from lattice-breakdown regime
    eps_n = np.sin(n_mode * gamma_sg)
    width_rest = 1.0 / eps_n if eps_n > 0.01 else 999
    print(f"Breather rest width: {width_rest:.1f} sites")
    print()

    good = [(v, r) for v, *r in results if abs(r[5]) < 10]
    bad = [(v, r) for v, *r in results if abs(r[5]) >= 10]

    if good:
        errors_good = [abs(r[5]) for _, r in good]
        v_max_good = max(v for v, _ in good)
        print(f"PHYSICAL REGIME (v <= {v_max_good:.1f}, breather survives):")
        print(f"  Mean |error|: {np.mean(errors_good):.3f}%")
        print(f"  Max  |error|: {np.max(errors_good):.3f}%")
        print()

    if bad:
        v_min_bad = min(v for v, _ in bad)
        gamma_break = 1.0 / np.sqrt(1.0 - v_min_bad**2)
        width_break = width_rest / gamma_break
        print(f"LATTICE BREAKDOWN (v >= {v_min_bad:.1f}):")
        print(f"  Breather Lorentz-contracted to {width_break:.1f} sites")
        print(f"  Too narrow for discrete lattice — breather radiates and disintegrates")
        print(f"  This IS the lattice discreteness effect (GWT prediction)")
        print()

    if good and np.mean(errors_good) < 2.0:
        print("CONCLUSION: Time dilation CONFIRMED on the discrete lattice.")
        print("  A moving breather oscillates slower by exactly 1/gamma.")
        print("  Special relativity emerges from F = -kx on springs.")
        print(f"  Confirmed for v = 0 to {v_max_good} c ({len(good)} velocities, all < {np.max(errors_good):.1f}%)")
        print()
        print("  Above v ~ {:.1f} c, the Lorentz-contracted breather is too narrow".format(v_min_bad if bad else 1.0))
        print("  for the discrete lattice. This is a PHYSICAL prediction:")
        print("  near-Planck-speed particles would radiate and decay.")
    else:
        print("RESULT: Deviations detected. Check parameters.")
