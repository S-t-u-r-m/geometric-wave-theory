"""
Vector Field Breather Simulation — VP from Dynamics
====================================================
The CORRECT GWT lattice has a 3-component displacement vector at each site.
The VP emerges from the nonlinear coupling between components.

Lagrangian:
  L = sum_<i,j> (1/2)|phi_i - phi_j|^2  +  sum_i (1/pi^2)(1 - cos(pi*|phi_i|))

The on-site potential V(|phi|) = (1/pi^2)(1 - cos(pi*|phi|)) depends on the
MAGNITUDE of the vector. Expanding to 4th order:
  V ~ |phi|^2/2 - pi^2 |phi|^4/24 + ...
  |phi|^4 = (phi_x^2 + phi_y^2 + phi_z^2)^2
          = phi_x^4 + phi_y^4 + phi_z^4 + 2(phi_x^2*phi_y^2 + ...)

The cross terms phi_x^2 * phi_y^2 ARE the T1u x T1u coupling that creates VP.
A scalar field can't generate these. Only the vector field can.

EOM for component a at site i:
  phi_ddot_i^a = sum_{j in NN} (phi_j^a - phi_i^a)
                 - (phi_i^a / |phi_i|) * (1/pi) * sin(pi * |phi_i|)

Test plan:
  Part 1: Single breather stability (T1u_x mode on 3D lattice)
  Part 2: Two breathers — parallel (both x) vs orthogonal (x and y)
  Part 3: Extract Oh channel decomposition from the orientation dependence
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

outfile = os.path.join(os.path.dirname(__file__), "breather_vector_vp_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("VECTOR FIELD BREATHER — VP FROM DYNAMICS")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"gamma = {gamma:.10f}")
report("")

# ============================================================
# VECTOR FIELD ENGINE
# ============================================================
N = 32  # lattice size per axis
dt = 0.04  # needs to be smaller for vector field (3 DOF per site)

def vec_magnitude(phi):
    """Compute |phi| at each site. phi shape: (N,N,N,3)"""
    return xp.sqrt(xp.sum(phi**2, axis=-1, keepdims=True) + 1e-30)

def vec_acceleration(phi):
    """EOM: phi_ddot = laplacian(phi) - (phi/|phi|) * (1/pi) * sin(pi*|phi|)

    phi shape: (N,N,N,3)
    Returns: acceleration with same shape
    """
    # Laplacian: sum of (phi_j - phi_i) over 6 NN, for each component
    lap = xp.zeros_like(phi)
    for axis in range(3):
        lap += xp.roll(phi, 1, axis=axis) + xp.roll(phi, -1, axis=axis)
    lap -= 6 * phi  # subtract 2d * phi_i

    # On-site force: -dV/dphi_a = -(phi_a/|phi|) * (1/pi) * sin(pi*|phi|)
    mag = vec_magnitude(phi)  # (N,N,N,1)
    # At |phi| = 0: force = -phi_a (linear limit). Handle smoothly:
    force_mag = xp.where(mag > 1e-10,
                          (1.0/PI) * xp.sin(PI * mag) / mag,
                          xp.ones_like(mag))  # lim_{r->0} sin(pi*r)/(pi*r) = 1 -> force = phi
    on_site_force = phi * force_mag  # (N,N,N,3)

    return lap - on_site_force

def vec_energy(phi, phi_dot):
    """Total energy of vector field on lattice."""
    # Kinetic
    KE = 0.5 * xp.sum(phi_dot**2)
    # On-site potential: V = (1/pi^2)(1 - cos(pi*|phi|))
    mag = vec_magnitude(phi)[:,:,:,0]  # (N,N,N)
    PE_site = xp.sum((1.0/PI**2) * (1 - xp.cos(PI * mag)))
    # NN coupling: (1/2) * sum |phi_j - phi_i|^2, positive direction only
    PE_coup = 0.0
    for axis in range(3):
        dphi = xp.roll(phi, -1, axis=axis) - phi  # (N,N,N,3)
        PE_coup += 0.5 * xp.sum(dphi**2)
    return float(KE + PE_site + PE_coup)

# ============================================================
# PART 1: Single vector breather stability
# ============================================================
report("PART 1: SINGLE VECTOR BREATHER (T1u modes)")
report("-" * 60)

center = N // 2
omega_1 = np.cos(gamma)
eps_1 = np.sin(gamma)

report(f"Lattice: {N}^3 = {N**3:,} sites x 3 components = {N**3*3:,} DOF")
report(f"n=1 breather: omega = {omega_1:.6f}, eps = {eps_1:.6f}")
report(f"Breather width: {1/eps_1:.1f} sites")
report("")

# Test each polarization
for pol_name, pol_axis in [("T1u_x", 0), ("T1u_y", 1), ("T1u_z", 2)]:
    # Breather along x-axis, polarized in pol_axis direction
    # Velocity profile: v_a(x) = (4/pi)*eps/(omega*cosh(eps*x)) * delta_{a, pol_axis}
    # Uniform in y,z (quasi-1D)
    phi = xp.zeros((N, N, N, 3), dtype=np.float64)
    phi_dot = xp.zeros_like(phi)

    x = xp.arange(N, dtype=np.float64) - center
    v_profile = (4.0/PI) * eps_1 / (omega_1 * xp.cosh(eps_1 * x) + 1e-30)

    # Set velocity along pol_axis, uniform in y,z
    phi_dot[:, :, :, pol_axis] = v_profile[:, None, None]

    phi_old = phi - dt * phi_dot

    # Evolve for 20 periods
    period = 2*PI / omega_1
    N_steps = int(20 * period / dt)
    rec = max(1, N_steps // 5000)

    E0 = vec_energy(phi, phi_dot)
    ts = []
    energies = []

    t0 = time.time()
    for step in range(N_steps):
        acc = vec_acceleration(phi)
        phi_new = 2*phi - phi_old + dt**2 * acc
        phi_old = phi.copy()
        phi = phi_new

        if step % rec == 0:
            ts.append(float(phi[center, center, center, pol_axis]))
            v_approx = (phi - phi_old) / dt
            energies.append(vec_energy(phi, v_approx))

    elapsed = time.time() - t0
    E_final = energies[-1]

    # Frequency
    ts = np.array(ts)
    ts_mean = ts - np.mean(ts)
    crossings = []
    for i in range(len(ts_mean)-1):
        if ts_mean[i] <= 0 and ts_mean[i+1] > 0:
            t_cross = i + (-ts_mean[i]) / (ts_mean[i+1] - ts_mean[i] + 1e-30)
            crossings.append(t_cross)
    if len(crossings) >= 3:
        T_meas = np.median(np.diff(crossings)) * dt * rec
        omega_meas = 2*PI / T_meas
    else:
        omega_meas = 0

    # Check cross-polarization leakage
    phi_np = phi.get() if GPU else phi
    other_axes = [a for a in range(3) if a != pol_axis]
    max_main = np.max(np.abs(phi_np[:,:,:,pol_axis]))
    max_cross = max(np.max(np.abs(phi_np[:,:,:,a])) for a in other_axes)
    leakage = max_cross / (max_main + 1e-30) * 100

    shift = (omega_meas - omega_1) / omega_1 * 100 if omega_meas > 0 else -999
    dE = (E_final - E0) / E0 * 100

    report(f"  {pol_name}: omega={omega_meas:.6f} (shift {shift:+.2f}%), "
           f"dE={dE:+.2f}%, leakage={leakage:.1f}%, ({elapsed:.0f}s)")

report("")

# ============================================================
# PART 2: Two breathers — parallel vs orthogonal
# ============================================================
report("PART 2: TWO-BREATHER INTERACTIONS")
report("-" * 60)
report("Test: do parallel (x,x) and orthogonal (x,y) breathers interact differently?")
report("The difference IS the VP — the Oh channel decomposition.")
report("")

# First, get single breather energy
phi_single = xp.zeros((N, N, N, 3), dtype=np.float64)
v_single = xp.zeros_like(phi_single)
x = xp.arange(N, dtype=np.float64) - center
v_profile = (4.0/PI) * eps_1 / (omega_1 * xp.cosh(eps_1 * x) + 1e-30)
v_single[:, :, :, 0] = v_profile[:, None, None]  # T1u_x
phi_old_s = phi_single - dt * v_single

# Evolve and measure energy
period = 2*PI / omega_1
N_steps_ref = int(25 * period / dt)
settle = int(5 * period / dt)
rec_ref = max(1, (N_steps_ref - settle) // 3000)

E_singles = []
for step in range(N_steps_ref):
    acc = vec_acceleration(phi_single)
    phi_new = 2*phi_single - phi_old_s + dt**2 * acc
    phi_old_s = phi_single.copy()
    phi_single = phi_new
    if step >= settle and step % rec_ref == 0:
        v_approx = (phi_single - phi_old_s) / dt
        E_singles.append(vec_energy(phi_single, v_approx))

E_single_avg = np.mean(E_singles)
E_single_std = np.std(E_singles)
report(f"Single breather (T1u_x): E = {E_single_avg:.6f} +/- {E_single_std:.6f}")
report("")

# Scan separations for parallel and orthogonal configurations
R_values = [4, 6, 8, 10, 12, 14, 16, 20, 24]

report(f"{'R':>4} {'E_parallel':>12} {'E_ortho':>12} {'V_par':>12} {'V_orth':>12} "
       f"{'ratio':>8} {'VP_frac':>8}")
report("-" * 78)

results_2b = []

for R_sep in R_values:
    pos_A = center - R_sep // 2
    pos_B = center + R_sep // 2

    x = xp.arange(N, dtype=np.float64)
    v_A = (4.0/PI) * eps_1 / (omega_1 * xp.cosh(eps_1 * (x - pos_A)) + 1e-30)
    v_B = (4.0/PI) * eps_1 / (omega_1 * xp.cosh(eps_1 * (x - pos_B)) + 1e-30)

    # === PARALLEL: both T1u_x ===
    phi_par = xp.zeros((N, N, N, 3), dtype=np.float64)
    v_par = xp.zeros_like(phi_par)
    v_par[:, :, :, 0] = (v_A - v_B)[:, None, None]  # opposite phase (bonding)
    phi_old_p = phi_par - dt * v_par

    E_pars = []
    for step in range(N_steps_ref):
        acc = vec_acceleration(phi_par)
        phi_new = 2*phi_par - phi_old_p + dt**2 * acc
        phi_old_p = phi_par.copy()
        phi_par = phi_new
        if step >= settle and step % rec_ref == 0:
            v_approx = (phi_par - phi_old_p) / dt
            E_pars.append(vec_energy(phi_par, v_approx))

    E_par_avg = np.mean(E_pars)

    # === ORTHOGONAL: A is T1u_x, B is T1u_y ===
    phi_ort = xp.zeros((N, N, N, 3), dtype=np.float64)
    v_ort = xp.zeros_like(phi_ort)
    v_ort[:, :, :, 0] = v_A[:, None, None]   # A polarized along x
    v_ort[:, :, :, 1] = -v_B[:, None, None]   # B polarized along y (opposite phase)
    phi_old_o = phi_ort - dt * v_ort

    E_orts = []
    for step in range(N_steps_ref):
        acc = vec_acceleration(phi_ort)
        phi_new = 2*phi_ort - phi_old_o + dt**2 * acc
        phi_old_o = phi_ort.copy()
        phi_ort = phi_new
        if step >= settle and step % rec_ref == 0:
            v_approx = (phi_ort - phi_old_o) / dt
            E_orts.append(vec_energy(phi_ort, v_approx))

    E_ort_avg = np.mean(E_orts)

    # Interaction energies
    V_par = E_par_avg - 2 * E_single_avg
    V_ort = E_ort_avg - 2 * E_single_avg

    # The ratio V_ort/V_par should reveal the Oh decomposition
    # If parallel = full A1g coupling, orthogonal = reduced coupling,
    # then the VP fraction = 1 - V_ort/V_par
    ratio = V_ort / V_par if abs(V_par) > 1e-10 else 0
    vp_frac = 1 - ratio if abs(V_par) > 1e-10 else 0

    results_2b.append((R_sep, E_par_avg, E_ort_avg, V_par, V_ort, ratio, vp_frac))

    report(f"{R_sep:4d} {E_par_avg:12.6f} {E_ort_avg:12.6f} {V_par:+12.6f} {V_ort:+12.6f} "
           f"{ratio:8.4f} {vp_frac:8.4f}")

report("")

# ============================================================
# PART 3: Analysis — Oh channel extraction
# ============================================================
report("PART 3: Oh CHANNEL ANALYSIS")
report("-" * 60)
report("")

# GWT predictions:
# T1u x T1u = A1g(1) + Eg(2) + T1g(3) + T2g(3) = 9 dims
# A1g fraction = 1/9 (scalar coupling, same orientation)
# Non-A1g fraction = 8/9 (VP correction)
#
# Parallel breathers: FULL coupling (all 9 channels)
# Orthogonal breathers: only channels that mix x and y
#   Of T1u_x x T1u_y:
#   A1g: 0 (scalar can't mix x,y)
#   Eg: 1 (the d_xy component)
#   T1g: 1 (the antisymmetric xy combination)
#   T2g: 1 (the off-diagonal symmetric)
#   Total: 3/9 of full coupling
#
# So V_ort/V_par should be ~ 3/9 = 1/3 if VP channels dominate
# Or V_ort/V_par ~ 0 if only A1g couples (pure scalar limit)

report("GWT predictions for T1u x T1u coupling:")
report(f"  Parallel (x,x): all 9 channels active")
report(f"  Orthogonal (x,y): 3/9 channels active (Eg_xy + T1g_xy + T2g_xy)")
report(f"  Predicted ratio V_ort/V_par = 3/9 = {3/9:.4f} (if nonlinear coupling works)")
report(f"  Predicted ratio V_ort/V_par = 0.0 (if only linear/scalar coupling)")
report("")

# Average ratio over separations
ratios = [r[5] for r in results_2b if abs(r[3]) > 1e-8]
if ratios:
    avg_ratio = np.mean(ratios)
    report(f"Measured average V_ort/V_par = {avg_ratio:.4f}")
    report(f"  If ≈ 0.33: VP channels active, T1u⊗T1u confirmed")
    report(f"  If ≈ 0.00: scalar limit, no cross-component coupling")
    report(f"  If ≈ 1.00: fully isotropic, no Oh structure")
report("")

# VP fraction
vp_fracs = [r[6] for r in results_2b if abs(r[3]) > 1e-8]
if vp_fracs:
    avg_vp = np.mean(vp_fracs)
    report(f"VP fraction (1 - V_ort/V_par) = {avg_vp:.4f}")
    report(f"  GWT prediction: (d^2-1)/d^2 = 8/9 = {8/9:.4f}")
report("")

report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"Total time: {(time.time()-time.time()):.1f} min")

log.close()
print(f"\nResults saved to: {outfile}")
