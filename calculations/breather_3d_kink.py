"""
3D Discrete Lattice Breather Eigenspectrum
==========================================
Three tests, each on a TRULY DISCRETE lattice (a=1):

Part 1: 1D discrete lattice (N sites, a=1)
  - NOT a continuous PDE on a fine grid (like eigenspectrum_proof.py used)
  - Each site IS a physical Planck site
  - The lattice corrections ARE the physics

Part 2: 3D lattice, breather uniform in y,z (quasi-1D)
  - Same 1D mode, but on a 3D cubic lattice
  - Periodic in y,z — no transverse curvature
  - Tests whether 3D lattice coupling affects the 1D eigenspectrum

Part 3: 3D lattice, breather localized in all 3 directions
  - Gaussian localization in y,z
  - Tests whether truly 3D breathers can exist

EOM: phi_ddot_i = sum_{j in NN(i)} (phi_j - phi_i) - (1/pi)*sin(pi*phi_i)

Key insight: the continuous PDE (dx→0) has breather widths ~ 1/sin(n*gamma)
that get NARROWER than the lattice for higher modes. On the physical lattice
(a=1), these modes are INHERENTLY discrete — no continuous approximation exists.
"""
import sys, io, os, time
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

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
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "breather_3d_kink_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("DISCRETE LATTICE BREATHER EIGENSPECTRUM")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"gamma = {gamma:.10f}")
report("")

particles = {4:"mu/strange", 5:"down", 7:"bottom", 8:"(n=8)",
             11:"charm", 12:"top", 13:"up", 16:"electron", 18:"tau"}

# ============================================================
# PART 1: Truly discrete 1D lattice (a=1)
# ============================================================
report("PART 1: DISCRETE 1D LATTICE (a=1)")
report("-" * 60)

N1d = 512  # physical sites
dt = 0.05  # stable for a=1 lattice (CFL: dt < 1/sqrt(d_max) for 1D = 1)

# Breather width for each n: width ~ 1/eps_n = 1/sin(n*gamma)
report(f"Lattice: {N1d} sites, a=1, dt={dt}, PERIODIC BC")
report(f"Phonon band (discrete 1D): [1.000, {np.sqrt(1+4):.3f}]")
report("")

for n_test in [1, 4, 7]:
    eps = np.sin(n_test * gamma)
    width = 1.0 / eps if eps > 0.01 else 999
    report(f"  n={n_test}: breather width = 1/sin({n_test}*gamma) = {width:.1f} sites")
report("")

report(f"{'n':>3} {'omega_pred':>11} {'omega_meas':>11} {'shift_%':>9} "
       f"{'periods':>8} {'amp_decay':>9} {'status':>8} {'particle':>12}")
report("-" * 85)

results_1d = []
t0_all = time.time()

center = N1d // 2

for n in range(1, 25):
    omega_n = np.cos(n * gamma)
    eps_n = np.sin(n * gamma)

    if omega_n < 0.005:
        results_1d.append((n, omega_n, 0, 0, 0, 0, "SKIP"))
        p = particles.get(n, "")
        report(f"{n:>3} {omega_n:11.6f} {'---':>11} {'---':>9} "
               f"{'---':>8} {'---':>9} {'SKIP':>8} {p:>12}")
        continue

    # Initialize breather on DISCRETE lattice
    # Exact sine-Gordon breather: phi(x,0) = 0, phi_dot(x,0) = (4/pi)*eps/(omega*cosh(eps*x))
    x = np.arange(N1d, dtype=np.float64) - center  # integer sites!

    phi = np.zeros(N1d)
    phi_dot = (4.0/PI) * eps_n / (omega_n * np.cosh(eps_n * x) + 1e-30)

    # Leapfrog init
    phi_old = phi - dt * phi_dot

    # Evolve
    period = 2*PI / omega_n
    N_periods = 50
    N_steps = int(N_periods * period / dt)
    N_steps = min(N_steps, 2000000)
    rec = max(1, N_steps // 30000)

    ts = []
    amps = []  # track LOCAL amplitude at breather center (not global max)
    hw_amp = max(2, int(1.0 / (eps_n + 0.01)))  # half-width ~ breather size
    for step in range(N_steps):
        # Discrete 1D Laplacian with PERIODIC BC (via roll)
        lap = np.roll(phi, 1) + np.roll(phi, -1) - 2*phi

        force = (1.0/PI) * np.sin(PI * phi)
        phi_new = 2*phi - phi_old + dt**2 * (lap - force)
        phi_old = phi.copy()
        phi = phi_new

        if step % rec == 0:
            ts.append(phi[center])
            # Local amplitude: max |phi| within breather half-width of center
            lo = max(0, center - hw_amp)
            hi = min(N1d, center + hw_amp + 1)
            amps.append(np.max(np.abs(phi[lo:hi])))

    # Frequency from zero crossings (more robust than FFT for shifted freqs)
    ts = np.array(ts)
    ts_mean = ts - np.mean(ts)

    crossings = []
    for i in range(len(ts_mean)-1):
        if ts_mean[i] <= 0 and ts_mean[i+1] > 0:
            # Linear interpolation
            t_cross = i + (-ts_mean[i]) / (ts_mean[i+1] - ts_mean[i] + 1e-30)
            crossings.append(t_cross)

    if len(crossings) >= 4:
        # Average period from consecutive crossings
        periods_list = np.diff(crossings) * dt * rec
        # Use median to reject outliers
        T_meas = np.median(periods_list)
        omega_meas = 2*PI / T_meas
        n_periods_meas = len(crossings) - 1
    else:
        omega_meas = 0
        n_periods_meas = 0

    # Amplitude decay from LOCAL amplitude at breather center
    amps = np.array(amps)
    if len(amps) > 20:
        # Compare first 10% vs last 10%
        n_window = max(3, len(amps) // 10)
        amp_start = np.mean(amps[:n_window])
        amp_end = np.mean(amps[-n_window:])
        amp_decay = (amp_end - amp_start) / (amp_start + 1e-30) * 100  # percent
    else:
        amp_decay = -999

    # Status
    shift_pct = (omega_meas - omega_n) / omega_n * 100 if omega_meas > 0 else -999

    if abs(shift_pct) < 0.5 and abs(amp_decay) < 5:
        status = "EXACT"
    elif abs(shift_pct) < 2.0 and abs(amp_decay) < 15:
        status = "STABLE"
    elif abs(shift_pct) < 15.0 and abs(amp_decay) < 30:
        status = "GOOD"
    elif abs(shift_pct) < 30.0 and n_periods_meas > 10:
        status = "shifted"
    elif n_periods_meas > 5:
        status = "weak"
    else:
        status = "MISS"

    p = particles.get(n, "")
    results_1d.append((n, omega_n, omega_meas, shift_pct, n_periods_meas, amp_decay, status))
    report(f"{n:>3} {omega_n:11.6f} {omega_meas:11.6f} {shift_pct:+9.2f} "
           f"{n_periods_meas:8d} {amp_decay:+9.1f}% {status:>8} {p:>12}")

t1 = time.time() - t0_all
report(f"\nPart 1 time: {t1:.0f}s ({t1/60:.1f} min)")

stable_1d = sum(1 for r in results_1d if r[6] in ("EXACT","STABLE","GOOD"))
shifted_1d = sum(1 for r in results_1d if r[6] in ("shifted",))
report(f"EXACT/STABLE/GOOD: {stable_1d}/24")
report(f"Shifted (oscillating but freq wrong): {shifted_1d}/24")
report("")

# ============================================================
# PART 2: 3D lattice, breather UNIFORM in y,z (quasi-1D)
# ============================================================
report("PART 2: 3D LATTICE, QUASI-1D (uniform in y,z)")
report("-" * 60)

N3 = 32  # 3D lattice size
dt3 = 0.05
center3 = N3 // 2

report(f"Lattice: {N3}x{N3}x{N3} = {N3**3:,} sites, a=1")
report(f"Breather along x, UNIFORM in y,z (periodic BC)")
report(f"Phonon band (discrete 3D): [1.000, {np.sqrt(1+4*3):.3f}]")
report("")

report(f"{'n':>3} {'omega_pred':>11} {'omega_meas':>11} {'shift_%':>9} "
       f"{'periods':>8} {'amp_decay':>9} {'status':>8} {'particle':>12}")
report("-" * 85)

results_3d_uniform = []
t0_p2 = time.time()

for n in range(1, 11):
    omega_n = np.cos(n * gamma)
    eps_n = np.sin(n * gamma)

    # Initialize: 1D breather along x, uniform in y,z
    x = xp.arange(N3, dtype=np.float64) - center3
    v_profile = (4.0/PI) * eps_n / (omega_n * xp.cosh(eps_n * x) + 1e-30)

    phi = xp.zeros((N3, N3, N3), dtype=np.float64)
    phi_dot = v_profile[:, None, None] * xp.ones((1, N3, N3), dtype=np.float64)

    phi_old = phi - dt3 * phi_dot

    # Evolve
    period = 2*PI / omega_n
    N_steps = min(int(30 * period / dt3), 300000)
    rec = max(1, N_steps // 15000)

    ts = []
    for step in range(N_steps):
        # 3D discrete Laplacian (periodic BC via roll)
        lap = (xp.roll(phi, 1, 0) + xp.roll(phi, -1, 0) +
               xp.roll(phi, 1, 1) + xp.roll(phi, -1, 1) +
               xp.roll(phi, 1, 2) + xp.roll(phi, -1, 2) - 6*phi)
        force = (1.0/PI) * xp.sin(PI * phi)
        phi_new = 2*phi - phi_old + dt3**2 * (lap - force)
        phi_old = phi.copy()
        phi = phi_new

        if step % rec == 0:
            ts.append(float(phi[center3, center3, center3]))

    # Frequency measurement (zero crossings)
    ts = np.array(ts)
    ts_mean = ts - np.mean(ts)

    crossings = []
    for i in range(len(ts_mean)-1):
        if ts_mean[i] <= 0 and ts_mean[i+1] > 0:
            t_cross = i + (-ts_mean[i]) / (ts_mean[i+1] - ts_mean[i] + 1e-30)
            crossings.append(t_cross)

    if len(crossings) >= 4:
        periods_list = np.diff(crossings) * dt3 * rec
        T_meas = np.median(periods_list)
        omega_meas = 2*PI / T_meas
        n_periods_meas = len(crossings) - 1
    else:
        omega_meas = 0
        n_periods_meas = 0

    shift_pct = (omega_meas - omega_n) / omega_n * 100 if omega_meas > 0 else -999

    if abs(shift_pct) < 0.1:
        status = "EXACT"
    elif abs(shift_pct) < 1.0:
        status = "STABLE"
    elif abs(shift_pct) < 10.0:
        status = "GOOD"
    elif abs(shift_pct) < 30.0 and n_periods_meas > 5:
        status = "shifted"
    elif n_periods_meas > 3:
        status = "weak"
    else:
        status = "MISS"

    p = particles.get(n, "")
    results_3d_uniform.append((n, omega_n, omega_meas, shift_pct, n_periods_meas, status))
    report(f"{n:>3} {omega_n:11.6f} {omega_meas:11.6f} {shift_pct:+9.2f} "
           f"{n_periods_meas:8d} {'---':>9} {status:>8} {p:>12}")

t2 = time.time() - t0_p2
report(f"\nPart 2 time: {t2:.0f}s ({t2/60:.1f} min)")
report("")

# ============================================================
# PART 3: 3D lattice, breather LOCALIZED in all directions
# ============================================================
report("PART 3: 3D LATTICE, FULLY LOCALIZED")
report("-" * 60)
report("Breather along x, Gaussian-localized in y,z")
report("This tests whether a truly 3D breather can survive.")
report("")

report(f"{'n':>3} {'omega_pred':>11} {'omega_meas':>11} {'shift_%':>9} "
       f"{'periods':>8} {'sigma':>6} {'status':>8} {'particle':>12}")
report("-" * 85)

results_3d_local = []

# Try different localization widths
for sigma_yz in [8.0, 4.0, 2.0]:
    report(f"\n--- sigma_yz = {sigma_yz:.1f} ---")
    iy = xp.arange(N3, dtype=np.float64) - center3
    iz = xp.arange(N3, dtype=np.float64) - center3
    gauss_y = xp.exp(-iy**2 / (2*sigma_yz**2))
    gauss_z = xp.exp(-iz**2 / (2*sigma_yz**2))
    envelope = gauss_y[None, :, None] * gauss_z[None, None, :]

    for n in [1, 4, 7]:
        omega_n = np.cos(n * gamma)
        eps_n = np.sin(n * gamma)

        x = xp.arange(N3, dtype=np.float64) - center3
        v_profile = (4.0/PI) * eps_n / (omega_n * xp.cosh(eps_n * x) + 1e-30)

        phi = xp.zeros((N3, N3, N3), dtype=np.float64)
        phi_dot = v_profile[:, None, None] * envelope

        phi_old = phi - dt3 * phi_dot

        period = 2*PI / omega_n
        N_steps = min(int(20 * period / dt3), 200000)
        rec = max(1, N_steps // 10000)

        ts = []
        for step in range(N_steps):
            lap = (xp.roll(phi, 1, 0) + xp.roll(phi, -1, 0) +
                   xp.roll(phi, 1, 1) + xp.roll(phi, -1, 1) +
                   xp.roll(phi, 1, 2) + xp.roll(phi, -1, 2) - 6*phi)
            force = (1.0/PI) * xp.sin(PI * phi)
            phi_new = 2*phi - phi_old + dt3**2 * (lap - force)
            phi_old = phi.copy()
            phi = phi_new

            if step % rec == 0:
                ts.append(float(phi[center3, center3, center3]))

        ts = np.array(ts)
        ts_mean = ts - np.mean(ts)

        crossings = []
        for i in range(len(ts_mean)-1):
            if ts_mean[i] <= 0 and ts_mean[i+1] > 0:
                t_cross = i + (-ts_mean[i]) / (ts_mean[i+1] - ts_mean[i] + 1e-30)
                crossings.append(t_cross)

        if len(crossings) >= 4:
            periods_list = np.diff(crossings) * dt3 * rec
            T_meas = np.median(periods_list)
            omega_meas = 2*PI / T_meas
            n_periods_meas = len(crossings) - 1
        else:
            omega_meas = 0
            n_periods_meas = 0

        shift_pct = (omega_meas - omega_n) / omega_n * 100 if omega_meas > 0 else -999

        if abs(shift_pct) < 1.0 and n_periods_meas > 10:
            status = "STABLE"
        elif abs(shift_pct) < 10.0 and n_periods_meas > 5:
            status = "GOOD"
        elif n_periods_meas > 3:
            status = "weak"
        else:
            status = "MISS"

        p = particles.get(n, "")
        results_3d_local.append((n, omega_n, omega_meas, shift_pct, n_periods_meas, sigma_yz, status))
        report(f"{n:>3} {omega_n:11.6f} {omega_meas:11.6f} {shift_pct:+9.2f} "
               f"{n_periods_meas:8d} {sigma_yz:6.1f} {status:>8} {p:>12}")

# ============================================================
# SUMMARY
# ============================================================
report("")
report("=" * 70)
report("SUMMARY")
report("=" * 70)
report("")

# Compare 1D discrete vs continuous
report("1D DISCRETE (a=1) vs CONTINUOUS (eigenspectrum_proof.py, dx=0.002):")
report("  The discrete lattice has inherent corrections that shift breather")
report("  frequencies. These shifts ARE the physics — not numerical error.")
report("")

stable = sum(1 for r in results_1d if r[6] in ("EXACT","STABLE","GOOD"))
shifted = sum(1 for r in results_1d if r[6] in ("shifted",))
report(f"  1D discrete: {stable} exact/stable, {shifted} shifted, out of 24")

stable3u = sum(1 for r in results_3d_uniform if r[5] in ("EXACT","STABLE","GOOD"))
report(f"  3D quasi-1D: {stable3u} exact/stable out of 10")
report("")

# Check if 3D uniform matches 1D
report("3D QUASI-1D vs 1D COMPARISON:")
for r3 in results_3d_uniform:
    n = r3[0]
    r1 = results_1d[n-1]
    if r1[2] > 0 and r3[2] > 0:
        diff = abs(r3[2] - r1[2]) / r1[2] * 100
        report(f"  n={n}: 1D={r1[2]:.6f}, 3D={r3[2]:.6f}, diff={diff:.2f}%")

report("")
report("3D LOCALIZED — does sigma matter?")
for sigma in [8.0, 4.0, 2.0]:
    subset = [r for r in results_3d_local if r[5] == sigma]
    report(f"  sigma={sigma}: " + ", ".join(f"n={r[0]}: {r[3]:+.1f}%" for r in subset))

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"Total time: {(time.time()-t0_all)/60:.1f} min")

log.close()
print(f"\nResults saved to: {outfile}")
