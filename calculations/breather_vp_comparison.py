"""
VP Self-Energy: Scalar vs Vector Breather
==========================================
Measure the DIFFERENCE between:
  - Scalar: V = (1/pi^2)(1 - cos(pi*phi_x))  — no VP, single component
  - Vector: V = (1/pi^2)(1 - cos(pi*|phi|))   — VP from cross-component coupling

The VP correction should appear as a small shift in:
  1. Breather frequency (omega)
  2. Breather energy (E)
  3. The shift should scale as alpha^2 × Oh_fraction

Both on the same 3D discrete lattice (32^3, a=1, periodic BC).
Same initial conditions. Same integrator. Only the potential differs.
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
alpha = np.exp(-(2/6) * (2**7/PI**2 + np.log(6)))

outfile = os.path.join(os.path.dirname(__file__), "breather_vp_comparison_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("VP SELF-ENERGY: SCALAR vs VECTOR POTENTIAL")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"alpha = {alpha:.6f}, alpha^2 = {alpha**2:.6e}")
report(f"gamma = {gamma:.10f}")
report("")

N = 32
dt = 0.04
center = N // 2

# ============================================================
# MODEL A: SCALAR — V = (1/pi^2)(1 - cos(pi*phi_x))
# Each component has its OWN independent cosine potential.
# No cross-component coupling. No VP.
# ============================================================
def acc_scalar(phi):
    """Scalar model: each component has independent cosine potential.
    phi shape: (N,N,N,3)
    Force_a = laplacian(phi_a) - (1/pi)*sin(pi*phi_a)
    """
    lap = xp.zeros_like(phi)
    for axis in range(3):
        lap += xp.roll(phi, 1, axis=axis) + xp.roll(phi, -1, axis=axis)
    lap -= 6 * phi
    # Independent cosine per component — NO cross-coupling
    force = (1.0/PI) * xp.sin(PI * phi)
    return lap - force

def energy_scalar(phi, phi_dot):
    KE = 0.5 * xp.sum(phi_dot**2)
    # On-site: independent cosine per component
    PE_site = xp.sum((1.0/PI**2) * (1 - xp.cos(PI * phi)))
    PE_coup = 0.0
    for axis in range(3):
        dphi = xp.roll(phi, -1, axis=axis) - phi
        PE_coup += 0.5 * xp.sum(dphi**2)
    return float(KE + PE_site + PE_coup)

# ============================================================
# MODEL B: VECTOR — V = (1/pi^2)(1 - cos(pi*|phi|))
# Magnitude-dependent cosine. Cross-component coupling through |phi|.
# VP emerges from phi_x^2 * phi_y^2 terms in the phi^4 expansion.
# ============================================================
def acc_vector(phi):
    """Vector model: cosine depends on |phi| = sqrt(phi_x^2+phi_y^2+phi_z^2).
    phi shape: (N,N,N,3)
    Force_a = laplacian(phi_a) - (phi_a/|phi|) * (1/pi)*sin(pi*|phi|)
    """
    lap = xp.zeros_like(phi)
    for axis in range(3):
        lap += xp.roll(phi, 1, axis=axis) + xp.roll(phi, -1, axis=axis)
    lap -= 6 * phi
    # Magnitude-dependent force
    mag = xp.sqrt(xp.sum(phi**2, axis=-1, keepdims=True) + 1e-30)
    force_mag = xp.where(mag > 1e-10,
                          (1.0/PI) * xp.sin(PI * mag) / mag,
                          xp.ones_like(mag))
    force = phi * force_mag
    return lap - force

def energy_vector(phi, phi_dot):
    KE = 0.5 * xp.sum(phi_dot**2)
    mag = xp.sqrt(xp.sum(phi**2, axis=-1) + 1e-30)
    PE_site = xp.sum((1.0/PI**2) * (1 - xp.cos(PI * mag)))
    PE_coup = 0.0
    for axis in range(3):
        dphi = xp.roll(phi, -1, axis=axis) - phi
        PE_coup += 0.5 * xp.sum(dphi**2)
    return float(KE + PE_site + PE_coup)

# ============================================================
# RUN BOTH MODELS WITH IDENTICAL INITIAL CONDITIONS
# ============================================================
def run_model(acc_fn, energy_fn, model_name, n_mode, pol_axis=0):
    """Run a single breather and measure frequency + energy precisely."""
    omega_n = np.cos(n_mode * gamma)
    eps_n = np.sin(n_mode * gamma)

    # Initialize: breather along x-axis, polarized in pol_axis
    phi = xp.zeros((N, N, N, 3), dtype=np.float64)
    phi_dot = xp.zeros_like(phi)

    x = xp.arange(N, dtype=np.float64) - center
    v_profile = (4.0/PI) * eps_n / (omega_n * xp.cosh(eps_n * x) + 1e-30)
    phi_dot[:, :, :, pol_axis] = v_profile[:, None, None]

    phi_old = phi - dt * phi_dot

    # Evolve for many periods — need precision
    period = 2*PI / omega_n
    N_periods = 60
    N_steps = int(N_periods * period / dt)
    settle = int(5 * period / dt)
    rec = max(1, (N_steps - settle) // 20000)

    ts = []
    energies = []

    t0 = time.time()
    for step in range(N_steps):
        a = acc_fn(phi)
        phi_new = 2*phi - phi_old + dt**2 * a
        phi_old = phi.copy()
        phi = phi_new

        if step >= settle and step % rec == 0:
            ts.append(float(phi[center, center, center, pol_axis]))
            v_approx = (phi - phi_old) / dt
            energies.append(energy_fn(phi, v_approx))

    elapsed = time.time() - t0

    # Precise frequency from zero crossings
    ts = np.array(ts)
    ts_mean = ts - np.mean(ts)
    crossings = []
    for i in range(len(ts_mean)-1):
        if ts_mean[i] <= 0 and ts_mean[i+1] > 0:
            t_cross = i + (-ts_mean[i]) / (ts_mean[i+1] - ts_mean[i] + 1e-30)
            crossings.append(t_cross)

    if len(crossings) >= 10:
        # Use ALL crossing intervals for maximum precision
        periods_list = np.diff(crossings) * dt * rec
        T_meas = np.mean(periods_list)  # mean for precision
        T_std = np.std(periods_list)
        omega_meas = 2*PI / T_meas
        omega_err = 2*PI * T_std / T_meas**2  # error propagation
        n_periods_meas = len(crossings) - 1
    else:
        omega_meas = 0
        omega_err = 999
        n_periods_meas = 0

    E_avg = np.mean(energies)
    E_std = np.std(energies)

    # Cross-polarization energy
    phi_np = phi.get() if GPU else phi
    other = [a for a in range(3) if a != pol_axis]
    max_main = np.max(np.abs(phi_np[:,:,:,pol_axis]))
    max_cross = max(np.max(np.abs(phi_np[:,:,:,a])) for a in other)

    return {
        'model': model_name,
        'n': n_mode,
        'omega_pred': omega_n,
        'omega_meas': omega_meas,
        'omega_err': omega_err,
        'E_avg': E_avg,
        'E_std': E_std,
        'n_periods': n_periods_meas,
        'leakage': max_cross / (max_main + 1e-30),
        'elapsed': elapsed
    }

# ============================================================
# COMPARISON: Multiple breather modes
# ============================================================
report("COMPARISON: Scalar vs Vector potential")
report("Same lattice, same IC, same integrator. Only the potential differs.")
report("")
report(f"{'n':>3} {'model':>8} {'omega_meas':>12} {'omega_err':>10} "
       f"{'E_avg':>12} {'E_std':>10} {'periods':>8} {'leak':>8}")
report("-" * 80)

all_results = []

for n_mode in [1, 2, 3, 4, 5, 6, 7]:
    r_scalar = run_model(acc_scalar, energy_scalar, "SCALAR", n_mode)
    r_vector = run_model(acc_vector, energy_vector, "VECTOR", n_mode)
    all_results.append((r_scalar, r_vector))

    for r in [r_scalar, r_vector]:
        report(f"{r['n']:>3} {r['model']:>8} {r['omega_meas']:12.8f} {r['omega_err']:10.6f} "
               f"{r['E_avg']:12.6f} {r['E_std']:10.6f} {r['n_periods']:8d} "
               f"{r['leakage']:8.4f}")

report("")

# ============================================================
# VP EXTRACTION
# ============================================================
report("VP EXTRACTION: omega_vector - omega_scalar")
report("-" * 60)
report("")
report(f"{'n':>3} {'omega_scalar':>14} {'omega_vector':>14} {'delta_omega':>14} "
       f"{'delta/omega':>12} {'alpha^2':>10}")
report("-" * 70)

for rs, rv in all_results:
    n = rs['n']
    delta = rv['omega_meas'] - rs['omega_meas']
    frac = delta / rs['omega_meas'] if rs['omega_meas'] > 0 else 0

    report(f"{n:>3} {rs['omega_meas']:14.8f} {rv['omega_meas']:14.8f} {delta:+14.8f} "
           f"{frac:+12.2e} {alpha**2:10.2e}")

report("")
report("If delta/omega ~ alpha^2 × Oh_fraction, VP is confirmed.")
report(f"alpha^2 = {alpha**2:.4e}")
report(f"alpha^2 × 8/9 = {alpha**2 * 8/9:.4e}")
report(f"alpha^2 × 8/d^2 = {alpha**2 * 8/9:.4e}")
report("")

# Energy comparison
report("ENERGY COMPARISON:")
report(f"{'n':>3} {'E_scalar':>14} {'E_vector':>14} {'delta_E':>14} {'delta_E/E':>12}")
report("-" * 60)

for rs, rv in all_results:
    n = rs['n']
    dE = rv['E_avg'] - rs['E_avg']
    dE_frac = dE / rs['E_avg'] if rs['E_avg'] != 0 else 0
    report(f"{n:>3} {rs['E_avg']:14.6f} {rv['E_avg']:14.6f} {dE:+14.6f} {dE_frac:+12.4e}")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))

log.close()
print(f"\nResults saved to: {outfile}")
