"""
Direct Minimization SCF — 7-Parameter Energy Optimization
============================================================
Instead of iterating mode-by-mode (which oscillates), minimize the
total energy DIRECTLY as a function of 7 mode amplitudes.

E_total(a_1, ..., a_7) = lowest eigenvalue of Hessian(background + Σ a_n × pert_n)

This is a 7-dimensional optimization — small enough for Nelder-Mead.
No convergence issues, no oscillation, finds the actual minimum.

For bonding:
  E_single = min_{a} E(single proton, a)
  E_double(R) = min_{a} E(two protons at R, a)  [symmetric: same a for both]
  V(R) = E_double(R) - 2 × E_single
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.optimize import minimize

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "scf_direct_minimize_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("DIRECT MINIMIZATION SCF — 7-PARAMETER BOND MODEL")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"d = {d}, gamma = {gamma:.10f}")
report("")

# ============================================================
# SETUP
# ============================================================
N = 512
x = np.arange(N, dtype=np.float64)
center = N // 2
kink_width = 3
n_modes = 7
modes = list(range(1, n_modes + 1))

def make_proton(x, pos, kw=3):
    return (4.0/PI) * (np.arctan(np.exp(x - pos + kw/2))
                      - np.arctan(np.exp(x - pos - kw/2)))

def build_hessian(phi_bg, extra_diag=None):
    diag = 2.0 + np.cos(PI * phi_bg)
    if extra_diag is not None:
        diag = diag + extra_diag
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()

# Pre-compute breather profiles at each position (reusable)
def precompute_profiles(positions):
    """Pre-compute -(π²/6) × φ_n²(x) for each mode at each position."""
    profiles = {}
    for n in modes:
        eps_n = np.sin(n * gamma)
        prof = np.zeros(N)
        for pos in positions:
            phi = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * (x - pos) + 1e-30))
            prof += -(PI**2 / 6.0) * phi**2
        profiles[n] = prof
    return profiles

# ============================================================
# THE ENERGY FUNCTION
# ============================================================
eval_count = [0]

def energy_function(a_vec, phi_bg, profiles):
    """Total energy as a function of 7 mode amplitudes.

    a_vec: array of 7 amplitudes (one per mode)
    Returns: lowest eigenvalue of Hessian with perturbation Σ a_n × profile_n
    """
    eval_count[0] += 1

    # Build total perturbation
    pert = np.zeros(N)
    for i, n in enumerate(modes):
        pert += a_vec[i] * profiles[n]

    # Check: perturbation shouldn't make diagonal negative (causes eigsh failure)
    diag_test = 2.0 + np.cos(PI * phi_bg) + pert
    if np.min(diag_test) < -10:
        return 100.0  # penalty for insane parameters

    try:
        H = build_hessian(phi_bg, extra_diag=pert)
        ev, _ = eigsh(H, k=1, which='SM', maxiter=3000)
        return ev[0]
    except Exception:
        return 100.0  # penalty

# ============================================================
# PART 1: SINGLE PROTON — FIND OPTIMAL MODE AMPLITUDES
# ============================================================
report("PART 1: SINGLE PROTON — OPTIMAL MODE AMPLITUDES")
report("-" * 55)

phi_single = make_proton(x, center, kink_width)
profiles_single = precompute_profiles([center])

# Start from no perturbation (a = 0)
a0 = np.zeros(n_modes)
eval_count[0] = 0
t0 = time.time()

result_s = minimize(
    energy_function, a0,
    args=(phi_single, profiles_single),
    method='Nelder-Mead',
    options={'maxiter': 5000, 'xatol': 1e-6, 'fatol': 1e-10, 'adaptive': True}
)

t_single = time.time() - t0
E_single = result_s.fun
a_single = result_s.x

report(f"Optimization: {eval_count[0]} evaluations in {t_single:.1f}s")
report(f"Converged: {result_s.success}")
report(f"E_single = {E_single:.8f}")
report("")

# Bare reference
H_bare = build_hessian(phi_single)
ev_bare, _ = eigsh(H_bare, k=1, which='SM')
E_bare = ev_bare[0]

report(f"Bare (no modes): {E_bare:.8f}")
report(f"Optimized:       {E_single:.8f}")
report(f"SCF correction:  {E_single - E_bare:+.8f}")
report("")

report("Optimal mode amplitudes (single proton):")
report(f"{'n':>3} {'eps_n':>8} {'a_n':>10} {'a_n/a_1':>10}")
report("-" * 35)
for i, n in enumerate(modes):
    eps_n = np.sin(n * gamma)
    ratio = a_single[i] / a_single[0] if abs(a_single[0]) > 1e-10 else 0
    report(f"  {n:3d} {eps_n:8.4f} {a_single[i]:+10.4f} {ratio:+10.4f}")
report("")

# ============================================================
# PART 2: TWO PROTONS — BOND CURVE
# ============================================================
report("PART 2: BOND CURVE — V(R) FROM DIRECT MINIMIZATION")
report("-" * 55)
report("Symmetric: same mode amplitudes at both protons.")
report("")

R_values = [6, 7, 8, 9, 10, 12, 14, 16, 20, 24, 30, 40]

report(f"{'R':>4} {'E_double':>12} {'V(R)':>12} {'V_bare':>12} "
       f"{'evals':>6} {'time':>6}")
report("-" * 58)

bond_curve = {}
a_bonded = {}

for R in R_values:
    pos_A = center - R // 2
    pos_B = center + R // 2
    phi_double = make_proton(x, pos_A, kink_width) + make_proton(x, pos_B, kink_width)
    profiles_double = precompute_profiles([pos_A, pos_B])

    # Bare reference
    H_bare_d = build_hessian(phi_double)
    ev_bare_d, _ = eigsh(H_bare_d, k=1, which='SM')
    V_bare = ev_bare_d[0] - 2 * E_bare

    # Optimize starting from single-proton solution
    eval_count[0] = 0
    t0 = time.time()

    result_d = minimize(
        energy_function, a_single,  # start from single-proton optimum
        args=(phi_double, profiles_double),
        method='Nelder-Mead',
        options={'maxiter': 5000, 'xatol': 1e-6, 'fatol': 1e-10, 'adaptive': True}
    )

    t_bond = time.time() - t0
    E_double = result_d.fun
    V_scf = E_double - 2 * E_single

    bond_curve[R] = {'E': E_double, 'V': V_scf, 'V_bare': V_bare, 'a': result_d.x.copy()}
    a_bonded[R] = result_d.x.copy()

    report(f"{R:4d} {E_double:12.6f} {V_scf:+12.8f} {V_bare:+12.8f} "
           f"{eval_count[0]:6d} {t_bond:5.1f}s")

report("")

# ============================================================
# PART 3: MORSE WELL ANALYSIS
# ============================================================
report("PART 3: MORSE WELL ANALYSIS")
report("-" * 55)

R_arr = np.array(sorted(bond_curve.keys()), dtype=float)
V_arr = np.array([bond_curve[int(R)]['V'] for R in R_arr])
V_bare_arr = np.array([bond_curve[int(R)]['V_bare'] for R in R_arr])

i_min = np.argmin(V_arr)
i_min_bare = np.argmin(V_bare_arr)

report(f"SCF: R_eq = {R_arr[i_min]:.0f}, D_e = {-V_arr[i_min]:.8f}")
report(f"Bare: R_eq = {R_arr[i_min_bare]:.0f}, D_e = {-V_bare_arr[i_min_bare]:.8f}")
report(f"V(R_max) = {V_arr[-1]:.8f} (should → 0)")
report("")

if V_arr[i_min] < -1e-6:
    report("*** MORSE WELL DETECTED ***")
    mask = (V_arr < -1e-8) & (R_arr > R_arr[i_min])
    if np.sum(mask) > 2:
        coeffs = np.polyfit(R_arr[mask], np.log(-V_arr[mask]), 1)
        report(f"  Morse decay: {-coeffs[0]:.4f}")
report("")

# ============================================================
# PART 4: HOW MODE WEIGHTS CHANGE WITH R
# ============================================================
report("PART 4: MODE AMPLITUDE SHIFTS DURING BONDING")
report("-" * 55)
report("How do the 7 mode amplitudes change from isolated to bonded?")
report("")

header = f"{'R':>4}"
for n in modes:
    header += f"  {'da_'+str(n):>8}"
report(header)
report("-" * (4 + n_modes * 10))

for R in R_values:
    da = a_bonded[R] - a_single
    line = f"{R:4d}"
    for i in range(n_modes):
        line += f"  {da[i]:+8.4f}"
    report(line)

report("")

# Which modes shift most?
report("AVERAGE MODE SHIFTS (across all R):")
for i, n in enumerate(modes):
    shifts = [a_bonded[R][i] - a_single[i] for R in R_values]
    avg = np.mean(shifts)
    std = np.std(shifts)
    report(f"  n={n}: avg shift = {avg:+.4f} ± {std:.4f}")

report("")

# ============================================================
# PART 5: THE RATIO STRUCTURE
# ============================================================
report("PART 5: MODE AMPLITUDE RATIOS")
report("-" * 55)
report("Are the bonded-state ratios simple functions of n?")
report("")

report("Single proton ratios a_n/a_1:")
for i, n in enumerate(modes):
    r = a_single[i] / a_single[0] if abs(a_single[0]) > 1e-10 else 0
    report(f"  n={n}: {r:+.4f}")

report("")
report(f"At R=10 (bonded) ratios a_n/a_1:")
if 10 in a_bonded:
    for i, n in enumerate(modes):
        r = a_bonded[10][i] / a_bonded[10][0] if abs(a_bonded[10][0]) > 1e-10 else 0
        report(f"  n={n}: {r:+.4f}")

report("")

# ============================================================
# PART 6: ASYMMETRIC BONDING (different mode weights per proton)
# ============================================================
report("PART 6: ASYMMETRIC BONDING")
report("-" * 55)
report("Allow different mode amplitudes at each proton (14 parameters).")
report("Does breaking symmetry lower the energy further?")
report("")

R_asym = 10
pos_A = center - R_asym // 2
pos_B = center + R_asym // 2
phi_double = make_proton(x, pos_A, kink_width) + make_proton(x, pos_B, kink_width)

# Now each proton has its own profiles
profiles_A = {}
profiles_B = {}
for n in modes:
    eps_n = np.sin(n * gamma)
    phi_A = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * (x - pos_A) + 1e-30))
    phi_B = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * (x - pos_B) + 1e-30))
    profiles_A[n] = -(PI**2 / 6.0) * phi_A**2
    profiles_B[n] = -(PI**2 / 6.0) * phi_B**2

def energy_asym(a_vec_14, phi_bg):
    """Energy with 14 parameters: 7 per proton."""
    pert = np.zeros(N)
    for i, n in enumerate(modes):
        pert += a_vec_14[i] * profiles_A[n]          # proton A
        pert += a_vec_14[i + n_modes] * profiles_B[n]  # proton B

    diag_test = 2.0 + np.cos(PI * phi_bg) + pert
    if np.min(diag_test) < -10:
        return 100.0

    try:
        H = build_hessian(phi_bg, extra_diag=pert)
        ev, _ = eigsh(H, k=1, which='SM', maxiter=3000)
        return ev[0]
    except Exception:
        return 100.0

# Start symmetric: a_A = a_B = a_single
a0_asym = np.concatenate([a_single, a_single])
eval_count[0] = 0
t0 = time.time()

result_asym = minimize(
    energy_asym, a0_asym,
    args=(phi_double,),
    method='Nelder-Mead',
    options={'maxiter': 10000, 'xatol': 1e-6, 'fatol': 1e-10, 'adaptive': True}
)

t_asym = time.time() - t0
E_asym = result_asym.fun
a_asym = result_asym.x

E_sym = bond_curve[R_asym]['E'] if R_asym in bond_curve else 999

report(f"R = {R_asym}")
report(f"Symmetric (7 params):   E = {E_sym:.8f}")
report(f"Asymmetric (14 params): E = {E_asym:.8f}")
report(f"Symmetry breaking gain: {E_asym - E_sym:+.8f}")
report(f"({eval_count[0]} evaluations in {t_asym:.1f}s)")
report("")

if abs(E_asym - E_sym) > 1e-6:
    report("Asymmetric weights (proton A vs B):")
    report(f"{'n':>3} {'a_A':>10} {'a_B':>10} {'diff':>10}")
    report("-" * 35)
    for i, n in enumerate(modes):
        a_A = a_asym[i]
        a_B = a_asym[i + n_modes]
        report(f"  {n:3d} {a_A:+10.4f} {a_B:+10.4f} {a_A - a_B:+10.4f}")
    report("")
    report("Nonzero diff = symmetry breaking = ionic character!")
else:
    report("No significant symmetry breaking — bond is purely covalent at this R.")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report(f"Single proton: E = {E_single:.8f} (bare: {E_bare:.8f})")
report(f"SCF correction: {E_single - E_bare:+.8f}")
report("")

report("Bond curve V(R):")
for R in R_values:
    v = bond_curve[R]['V']
    marker = " <-- min" if R == R_arr[i_min] else ""
    report(f"  R={R:3d}: V = {v:+.8f}{marker}")

report("")
if V_arr[i_min] < 0:
    report(f"D_e = {-V_arr[i_min]:.8f} at R_eq = {R_arr[i_min]:.0f}")
else:
    report("No bonding detected.")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
