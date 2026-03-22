"""
VP from Quantum Vacuum: ZPE-Seeded Breather Simulation
=======================================================
The VP correction is a QUANTUM effect — it comes from zero-point fluctuations
in the vacuum scattering off the breather through the phi^4 nonlinearity.

Setup:
  - 3D discrete lattice (32^3), vector field (3 components)
  - n=1 breather polarized along x (classical excitation)
  - y and z components seeded with zero-point energy (quantum vacuum)
  - ZPE amplitude per mode: sqrt(hbar*omega / 2) ~ 1/sqrt(2*N^3) in natural units

Compare:
  A) Scalar cos(pi*phi_x): no VP possible (no cross-component coupling)
  B) Vector cos(pi*|phi|) WITHOUT ZPE: no VP (y,z stay at zero)
  C) Vector cos(pi*|phi|) WITH ZPE: VP emerges from phi^4 scattering

The VP should appear as:
  - Frequency shift: delta_omega / omega ~ alpha^2 * Oh_fraction
  - Energy shift: delta_E / E ~ same

We run multiple realizations (different random ZPE seeds) and average
to extract the systematic VP shift from the stochastic noise.
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

outfile = os.path.join(os.path.dirname(__file__), "breather_vp_zpe_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("VP FROM QUANTUM VACUUM: ZPE-SEEDED SIMULATION")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"alpha = {alpha:.6f}, alpha^2 = {alpha**2:.6e}")
report("")

N = 32
dt = 0.04
center = N // 2

# ============================================================
# ZPE INITIALIZATION
# ============================================================
# On a lattice with N^3 sites, the ZPE per mode is hbar*omega/2.
# In natural units (hbar=1), the phonon frequency omega ~ 1 (mass gap).
# ZPE per mode = 1/2.
# Total ZPE distributed over N^3 sites and 3 components:
#   amplitude per site per component ~ sqrt(1/(2*N^3))
#
# But this is the FULL ZPE. We can scale it to test:
#   zpe_scale = 1.0 = physical ZPE
#   zpe_scale < 1 = reduced vacuum (check linearity)
#   zpe_scale = 0 = no vacuum (classical, should give zero VP)

def make_zpe(shape, zpe_scale, rng):
    """Generate zero-point energy fluctuations.

    Uses random phases in Fourier space to create a thermalized vacuum state.
    Each Fourier mode gets energy hbar*omega_k/2 with random phase.
    """
    Nx, Ny, Nz, n_comp = shape

    phi_zpe = xp.zeros(shape, dtype=np.float64)
    v_zpe = xp.zeros(shape, dtype=np.float64)

    if zpe_scale == 0:
        return phi_zpe, v_zpe

    for comp in range(n_comp):
        # Random phases in Fourier space
        if GPU:
            phases = xp.asarray(rng.uniform(0, 2*PI, (Nx, Ny, Nz)))
            amplitudes = xp.asarray(rng.standard_normal((Nx, Ny, Nz)))
        else:
            phases = rng.uniform(0, 2*PI, (Nx, Ny, Nz))
            amplitudes = rng.standard_normal((Nx, Ny, Nz))

        # Dispersion relation: omega_k^2 = 1 + 4*sum sin^2(k_a/2)
        kx = xp.fft.fftfreq(Nx) * 2 * PI
        ky = xp.fft.fftfreq(Ny) * 2 * PI
        kz = xp.fft.fftfreq(Nz) * 2 * PI
        KX, KY, KZ = xp.meshgrid(kx, ky, kz, indexing='ij')
        omega_k = xp.sqrt(1.0 + 4*(xp.sin(KX/2)**2 + xp.sin(KY/2)**2 + xp.sin(KZ/2)**2))

        # ZPE amplitude per mode: sqrt(1/(2*omega_k * N^3))
        amp_k = zpe_scale * xp.sqrt(1.0 / (2 * omega_k * Nx * Ny * Nz))

        # phi(x) = sum_k amp_k * cos(k.x + phase_k) -> real-space via iFFT
        # v(x) = sum_k amp_k * omega_k * sin(k.x + phase_k)
        fk_phi = amp_k * xp.exp(1j * phases) * amplitudes
        fk_vel = amp_k * omega_k * xp.exp(1j * (phases + PI/2)) * amplitudes

        phi_zpe[:,:,:,comp] = xp.real(xp.fft.ifftn(fk_phi)) * Nx * Ny * Nz
        v_zpe[:,:,:,comp] = xp.real(xp.fft.ifftn(fk_vel)) * Nx * Ny * Nz

    return phi_zpe, v_zpe

# ============================================================
# MODELS
# ============================================================
def acc_scalar(phi):
    lap = xp.zeros_like(phi)
    for axis in range(3):
        lap += xp.roll(phi, 1, axis=axis) + xp.roll(phi, -1, axis=axis)
    lap -= 6 * phi
    force = (1.0/PI) * xp.sin(PI * phi)
    return lap - force

def acc_vector(phi):
    lap = xp.zeros_like(phi)
    for axis in range(3):
        lap += xp.roll(phi, 1, axis=axis) + xp.roll(phi, -1, axis=axis)
    lap -= 6 * phi
    mag = xp.sqrt(xp.sum(phi**2, axis=-1, keepdims=True) + 1e-30)
    force_mag = xp.where(mag > 1e-10,
                          (1.0/PI) * xp.sin(PI * mag) / mag,
                          xp.ones_like(mag))
    force = phi * force_mag
    return lap - force

def measure_frequency(phi_series, dt_rec):
    """Extract frequency from time series via zero crossings."""
    ts = np.array(phi_series)
    ts_mean = ts - np.mean(ts)
    crossings = []
    for i in range(len(ts_mean)-1):
        if ts_mean[i] <= 0 and ts_mean[i+1] > 0:
            t_cross = i + (-ts_mean[i]) / (ts_mean[i+1] - ts_mean[i] + 1e-30)
            crossings.append(t_cross)
    if len(crossings) >= 6:
        periods = np.diff(crossings) * dt_rec
        return 2*PI / np.mean(periods), 2*PI * np.std(periods) / np.mean(periods)**2
    return 0, 999

def run_breather(acc_fn, phi_init, v_init, n_periods=50):
    """Evolve and measure breather frequency precisely."""
    omega_1 = np.cos(gamma)
    period = 2*PI / omega_1

    phi = phi_init.copy()
    phi_old = phi - dt * v_init

    N_steps = int(n_periods * period / dt)
    settle = int(5 * period / dt)
    rec = max(1, (N_steps - settle) // 15000)

    ts = []
    for step in range(N_steps):
        a = acc_fn(phi)
        phi_new = 2*phi - phi_old + dt**2 * a
        phi_old = phi.copy()
        phi = phi_new
        if step >= settle and step % rec == 0:
            ts.append(float(phi[center, center, center, 0]))  # x-component

    omega, omega_err = measure_frequency(ts, dt * rec)
    return omega, omega_err

# ============================================================
# EXPERIMENT: Multiple ZPE realizations
# ============================================================
omega_1 = np.cos(gamma)
eps_1 = np.sin(gamma)

report(f"Breather: n=1, omega_pred = {omega_1:.8f}")
report(f"Lattice: {N}^3, dt={dt}")
report("")

# Breather initial condition (x-polarized, along x-axis)
x = xp.arange(N, dtype=np.float64) - center
v_breather = xp.zeros((N, N, N, 3), dtype=np.float64)
v_breather[:,:,:,0] = ((4.0/PI) * eps_1 / (omega_1 * xp.cosh(eps_1 * x) + 1e-30))[:, None, None]

# ZPE scales to test
zpe_scales = [0.0, 0.5, 1.0, 2.0, 4.0]
N_realizations = 8  # average over random seeds

report(f"ZPE scales: {zpe_scales}")
report(f"Realizations per scale: {N_realizations}")
report("")

# Reference: scalar model, no ZPE (the "bare" frequency)
report("REFERENCE: Scalar, no ZPE")
phi_bare = xp.zeros((N, N, N, 3), dtype=np.float64)
omega_bare, err_bare = run_breather(acc_scalar, phi_bare, v_breather, n_periods=60)
report(f"  omega_bare = {omega_bare:.8f} +/- {err_bare:.6f}")
report("")

# Main comparison
report(f"{'scale':>6} {'model':>8} {'omega_avg':>14} {'omega_std':>12} "
       f"{'delta':>14} {'delta/omega':>12}")
report("-" * 72)

all_data = []

for zpe_scale in zpe_scales:
    for model_name, acc_fn in [("SCALAR", acc_scalar), ("VECTOR", acc_vector)]:
        omegas = []
        for seed in range(N_realizations):
            rng = np.random.default_rng(seed + 42)

            # Initialize: breather + ZPE
            phi_zpe, v_zpe = make_zpe((N, N, N, 3), zpe_scale, rng)

            # For scalar model: ZPE only in x-component (consistent)
            # For vector model: ZPE in all 3 components
            if model_name == "SCALAR":
                phi_zpe[:,:,:,1] = 0
                phi_zpe[:,:,:,2] = 0
                v_zpe[:,:,:,1] = 0
                v_zpe[:,:,:,2] = 0

            phi_init = phi_zpe  # ZPE position fluctuations
            v_init = v_breather + v_zpe  # breather + ZPE velocity

            omega, _ = run_breather(acc_fn, phi_init, v_init, n_periods=50)
            if omega > 0:
                omegas.append(omega)

        if len(omegas) >= 3:
            omega_avg = np.mean(omegas)
            omega_std = np.std(omegas) / np.sqrt(len(omegas))  # standard error
            delta = omega_avg - omega_bare
            delta_frac = delta / omega_bare
        else:
            omega_avg = 0
            omega_std = 0
            delta = 0
            delta_frac = 0

        all_data.append((zpe_scale, model_name, omega_avg, omega_std, delta, delta_frac))
        report(f"{zpe_scale:6.1f} {model_name:>8} {omega_avg:14.8f} {omega_std:12.8f} "
               f"{delta:+14.8f} {delta_frac:+12.4e}")

    report("")

# ============================================================
# VP EXTRACTION
# ============================================================
report("VP EXTRACTION")
report("-" * 60)
report("")

report("The VP = (vector+ZPE) - (scalar+ZPE) at each ZPE scale.")
report("This isolates the cross-component coupling effect.")
report("")

report(f"{'scale':>6} {'scalar':>14} {'vector':>14} {'VP_shift':>14} "
       f"{'VP/omega':>12} {'alpha^2':>10}")
report("-" * 72)

for zpe_scale in zpe_scales:
    scalar_data = [d for d in all_data if d[0]==zpe_scale and d[1]=="SCALAR"]
    vector_data = [d for d in all_data if d[0]==zpe_scale and d[1]=="VECTOR"]
    if scalar_data and vector_data:
        s = scalar_data[0]
        v = vector_data[0]
        vp_shift = v[2] - s[2]  # omega_vector - omega_scalar
        vp_frac = vp_shift / omega_bare if omega_bare > 0 else 0
        report(f"{zpe_scale:6.1f} {s[2]:14.8f} {v[2]:14.8f} {vp_shift:+14.8f} "
               f"{vp_frac:+12.4e} {alpha**2:10.2e}")

report("")
report(f"Expected VP: delta_omega/omega ~ alpha^2 * Oh_fraction")
report(f"  alpha^2 = {alpha**2:.4e}")
report(f"  8/9 (free photon) = {8/9:.4f}")
report(f"  1/2^(d/2) (confined) = {1/2**(d/2):.4f}")
report(f"  alpha^2 * 8/9 = {alpha**2 * 8/9:.4e}")
report(f"  alpha^2 / 2sqrt(2) = {alpha**2 / (2*np.sqrt(2)):.4e}")
report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))

log.close()
print(f"\nResults saved to: {outfile}")
