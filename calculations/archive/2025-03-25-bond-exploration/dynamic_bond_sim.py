"""
Dynamic Bond Simulation — Full Sine-Gordon Time Evolution
==========================================================
Instead of static Hessian perturbation theory, EVOLVE the actual field.

Two kink-antikink pairs ("protons") on a 1D lattice, each with breather
modes 1-7 excited at their natural frequencies. The sine-Gordon equation
is evolved with leapfrog integration. The bond curve V(R) emerges from
the time-averaged total energy vs separation.

EOM: phi_ddot_i = Σ_{j∈NN(i)} (phi_j - phi_i) - (1/π) sin(π phi_i)

This is the FULL nonlinear dynamics — no perturbation theory, no
approximations beyond the discretization. The breathers interact
through the actual field equation, including all orders of coupling.

Measurement:
  - Total energy E(t) = kinetic + gradient + potential
  - Time-averaged <E>(R) for each separation R
  - V(R) = <E>(R) - <E>(∞) → the bond curve
  - Also track: energy in each proton's region, energy exchange rate

All parameters from GWT Lagrangian (d=3, zero free parameters).
"""
import sys, io, os, time
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "dynamic_bond_sim_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("DYNAMIC BOND SIMULATION — FULL SINE-GORDON EVOLUTION")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"d = {d}, gamma = {gamma:.10f}")
report("")

# ============================================================
# LATTICE SETUP
# ============================================================
N = 1024       # lattice sites (big enough for two protons + buffer)
dt = 0.05      # stable for a=1 (CFL: dt < 1)
kink_width = 3

x = np.arange(N, dtype=np.float64)
center = N // 2

report(f"Lattice: {N} sites, a=1, dt={dt}, periodic BC")
report(f"Kink width: {kink_width}")
report("")

# ============================================================
# INITIAL CONDITIONS
# ============================================================
def kink_profile(x, x_center):
    return (4.0/PI) * np.arctan(np.exp(x - x_center))

def antikink_profile(x, x_center):
    return 2.0 - (4.0/PI) * np.arctan(np.exp(x - x_center))

def make_proton(x, center, width):
    """Static kink-antikink pair (proton at rest)."""
    return kink_profile(x, center - width/2) + antikink_profile(x, center + width/2) - 2.0

def add_breather_kick(phi_dot, x, center, n, amplitude=1.0):
    """Add breather mode n to the velocity field.

    The sine-Gordon breather starts with:
      phi(x, 0) = 0  (we add to existing static background)
      phi_dot(x, 0) = (4/π) × eps_n / (omega_n × cosh(eps_n × (x - center)))

    The amplitude parameter scales the kick (1.0 = full breather).
    """
    eps_n = np.sin(n * gamma)
    omega_n = np.cos(n * gamma)
    kick = amplitude * (4.0/PI) * eps_n / (omega_n * np.cosh(eps_n * (x - center) + 1e-30))
    phi_dot += kick
    return phi_dot

# ============================================================
# ENERGY MEASUREMENT
# ============================================================
def total_energy(phi, phi_old, dt):
    """Total energy on the 1D lattice.

    E = Σ_i [ (1/2)(phi_dot_i)^2 + (1/2)(phi_i - phi_{i+1})^2 + (1/π²)(1 - cos(π phi_i)) ]

    phi_dot ≈ (phi - phi_old) / dt  (leapfrog midpoint velocity)
    """
    phi_dot = (phi - phi_old) / dt

    # Kinetic
    KE = 0.5 * np.sum(phi_dot**2)

    # Gradient (nearest-neighbor)
    grad = np.roll(phi, -1) - phi  # phi_{i+1} - phi_i
    GE = 0.5 * np.sum(grad**2)

    # Potential
    PE = (1.0/PI**2) * np.sum(1.0 - np.cos(PI * phi))

    return KE + GE + PE, KE, GE, PE

def region_energy(phi, phi_old, dt, lo, hi):
    """Energy in a subregion [lo, hi)."""
    phi_dot = (phi - phi_old) / dt

    KE = 0.5 * np.sum(phi_dot[lo:hi]**2)

    # Gradient within region (use periodic for edge)
    phi_ext = np.concatenate([phi[lo:hi], [phi[hi % N]]])
    grad = np.diff(phi_ext)
    GE = 0.5 * np.sum(grad**2)

    PE = (1.0/PI**2) * np.sum(1.0 - np.cos(PI * phi[lo:hi]))

    return KE + GE + PE

# ============================================================
# RUN SIMULATION FOR A GIVEN SEPARATION R
# ============================================================
def run_bond_sim(R, breather_modes, amplitude=0.3, n_periods=200, report_fn=None):
    """Evolve two protons at separation R with given breather modes.

    Returns: time-averaged total energy, energy fluctuation, diagnostics.

    amplitude: scale factor for breather kicks (< 1 to avoid blowup
               from overlapping modes; 0.3 = moderate excitation)
    """
    pos_A = center - R // 2
    pos_B = center + R // 2

    # Static background: two protons
    phi = make_proton(x, pos_A, kink_width) + make_proton(x, pos_B, kink_width)

    # Initial velocity: add breather kicks to both wells
    phi_dot = np.zeros(N)
    for n in breather_modes:
        phi_dot = add_breather_kick(phi_dot, x, pos_A, n, amplitude)
        phi_dot = add_breather_kick(phi_dot, x, pos_B, n, amplitude)

    # Leapfrog initialization: phi_old = phi - dt * phi_dot
    phi_old = phi - dt * phi_dot

    # Evolution parameters
    # Use the slowest breather (n=1) to define the period
    omega_1 = np.cos(gamma)
    T_1 = 2 * PI / omega_1
    N_steps = int(n_periods * T_1 / dt)

    # Thermalization: skip the first 20 periods to let transients die
    n_therm = int(20 * T_1 / dt)

    # Measurement interval
    rec_interval = max(1, int(T_1 / dt / 10))  # ~10 samples per period

    energies = []
    energies_A = []
    energies_B = []
    energies_mid = []

    # Region boundaries
    mid = (pos_A + pos_B) // 2
    region_A = (max(0, pos_A - 20), mid)
    region_B = (mid, min(N, pos_B + 20))

    for step in range(N_steps):
        # Discrete 1D Laplacian (periodic BC via roll)
        lap = np.roll(phi, 1) + np.roll(phi, -1) - 2 * phi

        # Force: -dV/dphi = -(1/π)sin(πφ)
        force = (1.0/PI) * np.sin(PI * phi)

        # Leapfrog update
        phi_new = 2 * phi - phi_old + dt**2 * (lap - force)
        phi_old = phi.copy()
        phi = phi_new

        # Record energy after thermalization
        if step >= n_therm and step % rec_interval == 0:
            E_tot, KE, GE, PE = total_energy(phi, phi_old, dt)
            energies.append(E_tot)

            E_A = region_energy(phi, phi_old, dt, region_A[0], region_A[1])
            E_B = region_energy(phi, phi_old, dt, region_B[0], region_B[1])
            energies_A.append(E_A)
            energies_B.append(E_B)

    energies = np.array(energies)
    energies_A = np.array(energies_A)
    energies_B = np.array(energies_B)

    E_mean = np.mean(energies)
    E_std = np.std(energies)

    return {
        'E_mean': E_mean,
        'E_std': E_std,
        'E_A_mean': np.mean(energies_A),
        'E_B_mean': np.mean(energies_B),
        'n_samples': len(energies),
    }

# ============================================================
# PART 1: SINGLE PROTON REFERENCE
# ============================================================
report("PART 1: SINGLE PROTON REFERENCE ENERGY")
report("-" * 55)

# Single proton with all 7 breather modes
# Use a large R (well separated) as the "isolated" reference
R_ref = 400  # far enough that protons don't interact

amplitudes_to_test = [0.1, 0.2, 0.3, 0.5]
report("Testing different breather excitation amplitudes:")
report(f"{'amp':>6} {'E_mean':>12} {'E_std':>10} {'samples':>8}")
report("-" * 40)

for amp in amplitudes_to_test:
    result = run_bond_sim(R_ref, list(range(1, 8)), amplitude=amp, n_periods=100)
    report(f"  {amp:6.2f} {result['E_mean']:12.4f} {result['E_std']:10.4f} "
           f"{result['n_samples']:8d}")

report("")

# Choose a good amplitude
amp_use = 0.2  # moderate — enough to excite all modes, not so much they blow up
report(f"Using amplitude = {amp_use}")
report("")

# Reference: isolated protons (large R)
ref_result = run_bond_sim(R_ref, list(range(1, 8)), amplitude=amp_use, n_periods=300)
E_ref = ref_result['E_mean']
report(f"Reference (R={R_ref}): <E> = {E_ref:.6f} ± {ref_result['E_std']:.6f}")
report("")

# ============================================================
# PART 2: BOND CURVE — V(R) FROM DYNAMICS
# ============================================================
report("PART 2: DYNAMIC BOND CURVE")
report("-" * 55)
report("V(R) = <E(R)> - <E(∞)>")
report("")

R_values = [8, 10, 12, 14, 16, 20, 25, 30, 40, 50, 60, 80]

report(f"{'R':>4} {'<E>':>12} {'<E>-E_ref':>12} {'E_std':>10} "
       f"{'E_A':>10} {'E_B':>10} {'samples':>8}")
report("-" * 72)

bond_curve = {}

for R in R_values:
    t0 = time.time()
    result = run_bond_sim(R, list(range(1, 8)), amplitude=amp_use, n_periods=300)
    elapsed = time.time() - t0

    V_R = result['E_mean'] - E_ref
    bond_curve[R] = result

    report(f"{R:4d} {result['E_mean']:12.6f} {V_R:+12.6f} {result['E_std']:10.4f} "
           f"{result['E_A_mean']:10.4f} {result['E_B_mean']:10.4f} "
           f"{result['n_samples']:8d}  [{elapsed:.1f}s]")

report("")

# ============================================================
# PART 3: BARE BOND CURVE (no breathers, for comparison)
# ============================================================
report("PART 3: BARE BOND CURVE (static protons, no breather excitation)")
report("-" * 55)

# With no breathers, the kink-antikink pairs just sit there
# The energy difference comes purely from the static field overlap
report(f"{'R':>4} {'E_static':>12} {'V_static':>12}")
report("-" * 30)

for R in R_values:
    pos_A = center - R // 2
    pos_B = center + R // 2
    phi_static = make_proton(x, pos_A, kink_width) + make_proton(x, pos_B, kink_width)

    # Static energy
    grad = np.roll(phi_static, -1) - phi_static
    GE = 0.5 * np.sum(grad**2)
    PE = (1.0/PI**2) * np.sum(1.0 - np.cos(PI * phi_static))
    E_static = GE + PE

    # Reference (isolated)
    phi_ref = make_proton(x, center - R_ref//2, kink_width) + make_proton(x, center + R_ref//2, kink_width)
    grad_ref = np.roll(phi_ref, -1) - phi_ref
    E_ref_static = 0.5 * np.sum(grad_ref**2) + (1.0/PI**2) * np.sum(1.0 - np.cos(PI * phi_ref))

    V_static = E_static - E_ref_static
    report(f"{R:4d} {E_static:12.6f} {V_static:+12.6f}")

report("")

# ============================================================
# PART 4: MODE-BY-MODE CONTRIBUTION
# ============================================================
report("PART 4: WHICH MODES CONTRIBUTE TO BONDING?")
report("-" * 55)
report("Run with subsets of breather modes to see individual contributions.")
report("")

R_test = 10
mode_sets = {
    'none': [],
    'n=1 only': [1],
    'n=1-3 (wide)': [1, 2, 3],
    'n=1-4': [1, 2, 3, 4],
    'n=5-7 (narrow)': [5, 6, 7],
    'n=1-7 (all)': [1, 2, 3, 4, 5, 6, 7],
}

report(f"R = {R_test}:")
report(f"{'modes':>20} {'<E>':>12} {'V(R)':>12} {'E_std':>10}")
report("-" * 58)

# References for each mode set
ref_results = {}
for name, mode_list in mode_sets.items():
    ref_r = run_bond_sim(R_ref, mode_list, amplitude=amp_use, n_periods=200)
    ref_results[name] = ref_r['E_mean']

for name, mode_list in mode_sets.items():
    result = run_bond_sim(R_test, mode_list, amplitude=amp_use, n_periods=200)
    V = result['E_mean'] - ref_results[name]
    report(f"{name:>20} {result['E_mean']:12.6f} {V:+12.6f} {result['E_std']:10.4f}")

report("")

# ============================================================
# PART 5: ENERGY EXCHANGE RATE
# ============================================================
report("PART 5: ENERGY FLOW BETWEEN PROTONS")
report("-" * 55)
report("Track how energy oscillates between the two proton regions.")
report("Fast oscillation = strong coupling, slow = weak coupling.")
report("")

# Run at R=10 with detailed energy recording
R_flow = 10
pos_A = center - R_flow // 2
pos_B = center + R_flow // 2

phi = make_proton(x, pos_A, kink_width) + make_proton(x, pos_B, kink_width)
phi_dot = np.zeros(N)
for n in range(1, 8):
    phi_dot = add_breather_kick(phi_dot, x, pos_A, n, amp_use)
    phi_dot = add_breather_kick(phi_dot, x, pos_B, n, amp_use)

phi_old = phi - dt * phi_dot

omega_1 = np.cos(gamma)
T_1 = 2 * PI / omega_1
N_steps = int(100 * T_1 / dt)
n_therm = int(10 * T_1 / dt)

mid = (pos_A + pos_B) // 2
E_A_series = []
E_B_series = []
t_series = []

rec = max(1, int(T_1 / dt / 20))

for step in range(N_steps):
    lap = np.roll(phi, 1) + np.roll(phi, -1) - 2 * phi
    force = (1.0/PI) * np.sin(PI * phi)
    phi_new = 2 * phi - phi_old + dt**2 * (lap - force)
    phi_old = phi.copy()
    phi = phi_new

    if step >= n_therm and step % rec == 0:
        E_A = region_energy(phi, phi_old, dt, max(0, pos_A-20), mid)
        E_B = region_energy(phi, phi_old, dt, mid, min(N, pos_B+20))
        E_A_series.append(E_A)
        E_B_series.append(E_B)
        t_series.append(step * dt)

E_A_series = np.array(E_A_series)
E_B_series = np.array(E_B_series)

# Energy imbalance: E_A - E_B should oscillate if energy flows between protons
imbalance = E_A_series - E_B_series
report(f"Energy imbalance (E_A - E_B) at R={R_flow}:")
report(f"  Mean: {np.mean(imbalance):.6f}")
report(f"  Std:  {np.std(imbalance):.6f}")
report(f"  Max:  {np.max(np.abs(imbalance)):.6f}")

# Frequency of energy oscillation (from FFT of imbalance)
if len(imbalance) > 20:
    fft_imb = np.abs(np.fft.rfft(imbalance - np.mean(imbalance)))
    freqs = np.fft.rfftfreq(len(imbalance), d=(t_series[1] - t_series[0]))
    peak_idx = np.argmax(fft_imb[1:]) + 1  # skip DC
    f_exchange = freqs[peak_idx]
    report(f"  Exchange frequency: {f_exchange:.6f}")
    report(f"  Exchange period: {1/f_exchange:.2f} Planck times" if f_exchange > 0 else "  No clear oscillation")
    report(f"  Compare omega_1 = {omega_1:.6f}")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report("Dynamic bond curve V(R) from full sine-Gordon evolution:")
report(f"{'R':>4} {'V(R)':>12}")
report("-" * 20)
for R in R_values:
    V = bond_curve[R]['E_mean'] - E_ref
    report(f"{R:4d} {V:+12.6f}")

# Check for Morse-like shape
V_arr = np.array([bond_curve[R]['E_mean'] - E_ref for R in R_values])
R_arr = np.array(R_values, dtype=float)
i_min = np.argmin(V_arr)

report("")
if V_arr[i_min] < -0.001:
    report(f"MORSE WELL DETECTED:")
    report(f"  R_eq = {R_arr[i_min]:.0f}")
    report(f"  D_e = {-V_arr[i_min]:.6f}")
    report(f"  Bond energy has EMERGED from the dynamics!")
elif V_arr[i_min] < 0:
    report(f"Weak attraction detected at R={R_arr[i_min]:.0f}, V={V_arr[i_min]:.6f}")
    report(f"May need longer evolution or different amplitude.")
else:
    report("No clear bonding detected.")
    report("The breather modes may need different amplitude or more evolution time.")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
