"""
SCF Continuation Method — Adiabatic Bond Curve
================================================
The mode amplitudes reorganize at each R — that IS the chemistry.
But the optimizer must follow the SAME basin continuously.

Method: start at large R (= two isolated protons), decrease R step
by step, using each solution as the starting point for the next.
The physical solution deforms smoothly.

R=100 → R=80 → R=60 → R=40 → R=30 → R=24 → R=20 → ... → R=4

Uses EXACT lattice potential: cos(π × φ_total), no approximation.
CPU scipy (faster than GPU for 512-site 1D problem).
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

outfile = os.path.join(os.path.dirname(__file__), "scf_continuation_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("SCF CONTINUATION METHOD — ADIABATIC BOND CURVE")
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

def kink_field(x, pos, kw=3):
    return (4.0/PI) * (np.arctan(np.exp(x - pos + kw/2))
                      - np.arctan(np.exp(x - pos - kw/2)))

def breather_field(x, pos, eps_n):
    return (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * (x - pos) + 1e-30))

def build_hessian(phi_total):
    diag = 2.0 + np.cos(PI * phi_total)
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()

def get_E0(phi_total):
    H = build_hessian(phi_total)
    ev, _ = eigsh(H, k=1, which='SM', maxiter=3000)
    return ev[0]

def precompute_breather_fields(positions):
    fields = {}
    for n in modes:
        eps_n = np.sin(n * gamma)
        f = np.zeros(N)
        for pos in positions:
            f += breather_field(x, pos, eps_n)
        fields[n] = f
    return fields

def total_field(phi_kink, bf, a_vec):
    phi = phi_kink.copy()
    for i, n in enumerate(modes):
        phi += a_vec[i] * bf[n]
    return phi

def energy_fn(a_vec, phi_kink, bf):
    phi = total_field(phi_kink, bf, a_vec)
    return get_E0(phi)

# ============================================================
# PART 1: SINGLE PROTON REFERENCE
# ============================================================
report("PART 1: SINGLE PROTON")
report("-" * 55)

phi_kink_s = kink_field(x, center, kink_width)
bf_s = precompute_breather_fields([center])

E_bare_s = get_E0(phi_kink_s)
report(f"Bare: {E_bare_s:.8f}")

# Optimize single proton
a0 = np.zeros(n_modes)
result_s = minimize(energy_fn, a0, args=(phi_kink_s, bf_s),
                    method='Powell', bounds=[(-1,1)]*n_modes,
                    options={'maxiter': 3000, 'ftol': 1e-10})

E_single = result_s.fun
a_single = result_s.x.copy()

report(f"Optimized: {E_single:.8f} (correction: {E_single - E_bare_s:+.8f})")
report(f"Amplitudes: {' '.join(f'{a:+.4f}' for a in a_single)}")
report("")

# ============================================================
# PART 2: CONTINUATION BOND CURVE — LARGE R TO SMALL R
# ============================================================
report("PART 2: CONTINUATION BOND CURVE")
report("-" * 55)
report("Starting at R=100, decreasing to R=4.")
report("Each R starts from the previous R's solution.")
report("")

# Dense R scan, from large to small
R_scan = [100, 80, 60, 50, 40, 35, 30, 27, 24, 22, 20, 18,
          16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4]

# Start from single-proton solution (= large R limit)
a_current = a_single.copy()

report(f"{'R':>4} {'E_double':>12} {'V(R)':>12} {'V_bare':>12} "
       f"{'|da|':>8} {'evals':>5} {'t':>5}")
report("-" * 62)

bond_curve = {}
a_history = {}

for R in R_scan:
    pos_A = center - R // 2
    pos_B = center + R // 2

    phi_kink_d = kink_field(x, pos_A, kink_width) + kink_field(x, pos_B, kink_width)
    bf_d = precompute_breather_fields([pos_A, pos_B])

    # Bare
    E_bare_d = get_E0(phi_kink_d)
    V_bare = E_bare_d - 2 * E_bare_s

    # Optimize from PREVIOUS solution (continuation!)
    t0 = time.time()
    result_d = minimize(energy_fn, a_current, args=(phi_kink_d, bf_d),
                        method='Powell', bounds=[(-1,1)]*n_modes,
                        options={'maxiter': 2000, 'ftol': 1e-10})
    t_opt = time.time() - t0

    E_double = result_d.fun
    a_new = result_d.x.copy()
    V_scf = E_double - 2 * E_single

    # How much did the amplitudes change from previous step?
    da = np.sqrt(np.sum((a_new - a_current)**2))

    bond_curve[R] = {'E': E_double, 'V': V_scf, 'V_bare': V_bare, 'a': a_new}
    a_history[R] = a_new.copy()

    report(f"{R:4d} {E_double:12.8f} {V_scf:+12.8f} {V_bare:+12.8f} "
           f"{da:8.5f} {result_d.nfev:5d} {t_opt:5.1f}s")

    # Update for next step
    a_current = a_new.copy()

report("")

# ============================================================
# PART 3: ALSO SCAN SMALL R TO LARGE R (check hysteresis)
# ============================================================
report("PART 3: REVERSE SCAN (R=4 → R=100) — HYSTERESIS CHECK")
report("-" * 55)
report("If the forward and reverse curves match, the solution is unique.")
report("If they differ, there are multiple basins (hysteresis).")
report("")

# Start from the R=4 solution
a_current_rev = a_history[R_scan[-1]].copy()

report(f"{'R':>4} {'V_fwd':>12} {'V_rev':>12} {'diff':>12}")
report("-" * 42)

bond_curve_rev = {}

for R in reversed(R_scan):
    pos_A = center - R // 2
    pos_B = center + R // 2

    phi_kink_d = kink_field(x, pos_A, kink_width) + kink_field(x, pos_B, kink_width)
    bf_d = precompute_breather_fields([pos_A, pos_B])

    result_rev = minimize(energy_fn, a_current_rev, args=(phi_kink_d, bf_d),
                          method='Powell', bounds=[(-1,1)]*n_modes,
                          options={'maxiter': 2000, 'ftol': 1e-10})

    V_rev = result_rev.fun - 2 * E_single
    V_fwd = bond_curve[R]['V']
    diff = V_rev - V_fwd

    bond_curve_rev[R] = V_rev
    a_current_rev = result_rev.x.copy()

    report(f"{R:4d} {V_fwd:+12.8f} {V_rev:+12.8f} {diff:+12.8f}")

report("")

# ============================================================
# PART 4: MODE AMPLITUDE EVOLUTION
# ============================================================
report("PART 4: HOW MODES EVOLVE WITH R (forward scan)")
report("-" * 55)

header = f"{'R':>4}"
for n in modes:
    header += f"  {'a'+str(n):>8}"
report(header)
report("-" * (4 + n_modes * 10))

for R in R_scan:
    line = f"{R:4d}"
    for i in range(n_modes):
        line += f"  {a_history[R][i]:+8.5f}"
    report(line)

report(f"\n ref")
line = "    "
for i in range(n_modes):
    line += f"  {a_single[i]:+8.5f}"
report(line)
report("")

# ============================================================
# PART 5: MORSE WELL ANALYSIS
# ============================================================
report("PART 5: MORSE WELL ANALYSIS")
report("-" * 55)

R_arr = np.array(R_scan, dtype=float)
V_arr = np.array([bond_curve[R]['V'] for R in R_scan])
V_bare_arr = np.array([bond_curve[R]['V_bare'] for R in R_scan])

i_min = np.argmin(V_arr)
i_min_bare = np.argmin(V_bare_arr)

report(f"Forward scan:")
report(f"  SCF:  R_eq = {R_arr[i_min]:.0f}, V_min = {V_arr[i_min]:+.8f}")
report(f"  Bare: R_eq = {R_arr[i_min_bare]:.0f}, V_min = {V_bare_arr[i_min_bare]:+.8f}")
report(f"  V(R=100) = {V_arr[0]:+.8f} (should → 0)")
report("")

if V_arr[i_min] < -1e-6:
    report(f"*** MORSE WELL: D_e = {-V_arr[i_min]:.8f} at R = {R_arr[i_min]:.0f} ***")

    mask = (V_arr < -1e-8) & (R_arr < R_arr[i_min])  # repulsive wall (small R)
    mask2 = (V_arr < -1e-8) & (R_arr > R_arr[i_min])  # attractive tail (large R)
    if np.sum(mask2) > 2:
        coeffs = np.polyfit(R_arr[mask2], np.log(-V_arr[mask2]), 1)
        report(f"  Morse tail decay: {-coeffs[0]:.4f}")
elif np.min(V_arr) < 0:
    report(f"Weak attraction: V_min = {np.min(V_arr):+.8f} at R = {R_arr[np.argmin(V_arr)]:.0f}")
else:
    report("No bonding detected.")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report(f"Single proton: E = {E_single:.8f}")
report(f"Method: continuation (adiabatic following from R=100 to R=4)")
report("")

# Check smoothness
dV = np.diff(V_arr)
smooth = np.all(np.abs(np.diff(dV)) < 0.5 * np.max(np.abs(dV)))
report(f"Bond curve smoothness: {'SMOOTH' if smooth else 'NOT SMOOTH'}")
report(f"V(R=100) = {V_arr[0]:+.8f} (should ≈ 0)")
report(f"V(R=4) = {V_arr[-1]:+.8f}")
report("")

# Hysteresis
V_rev_arr = np.array([bond_curve_rev[R] for R in R_scan])
max_hysteresis = np.max(np.abs(V_arr - V_rev_arr))
report(f"Max hysteresis (fwd vs rev): {max_hysteresis:.8f}")
if max_hysteresis < 0.01:
    report("→ Solution is UNIQUE (no hysteresis)")
else:
    report("→ HYSTERESIS detected — multiple basins exist")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
