"""
SCF Mode Progression — 1 Mode, 2 Modes, 3 Modes, ...
=======================================================
Start with 1 breather mode and get a clean bond curve.
Then add modes one at a time to see when/if the landscape breaks.

Uses exact lattice potential: cos(π × φ_total)
Continuation method: R from large to small.
At each mode count, the previous mode count's solution seeds the next.
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

outfile = os.path.join(os.path.dirname(__file__), "scf_mode_progression_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("SCF MODE PROGRESSION — 1, 2, 3, ... MODES")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("")

# ============================================================
# SETUP
# ============================================================
N = 512
x = np.arange(N, dtype=np.float64)
center = N // 2
kw = 3

def kink_field(x, pos):
    return (4.0/PI) * (np.arctan(np.exp(x - pos + kw/2))
                      - np.arctan(np.exp(x - pos - kw/2)))

def breather_field(x, pos, n):
    eps_n = np.sin(n * gamma)
    return (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * (x - pos) + 1e-30))

def get_E0(phi_total):
    diag = 2.0 + np.cos(PI * phi_total)
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    H = H.tocsr()
    ev, _ = eigsh(H, k=1, which='SM', maxiter=3000)
    return ev[0]

# Bare references
phi_single = kink_field(x, center)
E_bare_single = get_E0(phi_single)
report(f"Bare single proton: E0 = {E_bare_single:.8f}")
report("")

# R values for bond curve (continuation: large to small)
R_scan = [80, 60, 40, 30, 24, 20, 18, 16, 14, 12, 11, 10, 9, 8, 7, 6, 5, 4]

# ============================================================
# RUN FOR EACH MODE COUNT
# ============================================================
all_results = {}

for n_modes in [1, 2, 3, 4, 5, 7]:
    mode_list = list(range(1, n_modes + 1))
    report(f"{'='*60}")
    report(f"N_MODES = {n_modes}: modes {mode_list}")
    report(f"{'='*60}")

    # Pre-compute breather fields for single proton
    bf_fields_single = {}
    for n in mode_list:
        bf_fields_single[n] = breather_field(x, center, n)

    def energy_single(a_vec):
        phi = phi_single.copy()
        for i, n in enumerate(mode_list):
            phi += a_vec[i] * bf_fields_single[n]
        return get_E0(phi)

    # Optimize single proton
    a0 = np.zeros(n_modes)
    res_s = minimize(energy_single, a0, method='Powell',
                     bounds=[(-1, 1)] * n_modes,
                     options={'maxiter': 2000, 'ftol': 1e-10})
    E_single_opt = res_s.fun
    a_single_opt = res_s.x.copy()

    report(f"Single proton: E = {E_single_opt:.8f} "
           f"(correction: {E_single_opt - E_bare_single:+.8f})")
    report(f"Amplitudes: {' '.join(f'{a:+.5f}' for a in a_single_opt)}")
    report("")

    # Bond curve with continuation
    report(f"{'R':>4} {'V_scf':>12} {'V_bare':>12} {'|da|':>8} {'evals':>5}")
    report("-" * 48)

    a_current = a_single_opt.copy()
    bond_data = {}

    for R in R_scan:
        pos_A = center - R // 2
        pos_B = center + R // 2
        phi_kink_d = kink_field(x, pos_A) + kink_field(x, pos_B)

        # Precompute breather fields for both protons
        bf_double = {}
        for n in mode_list:
            bf_double[n] = breather_field(x, pos_A, n) + breather_field(x, pos_B, n)

        # Bare
        E_bare_d = get_E0(phi_kink_d)
        V_bare = E_bare_d - 2 * E_bare_single

        def energy_double(a_vec):
            phi = phi_kink_d.copy()
            for i, n in enumerate(mode_list):
                phi += a_vec[i] * bf_double[n]
            return get_E0(phi)

        res_d = minimize(energy_double, a_current, method='Powell',
                         bounds=[(-1, 1)] * n_modes,
                         options={'maxiter': 2000, 'ftol': 1e-10})

        E_double = res_d.fun
        V_scf = E_double - 2 * E_single_opt
        da = np.sqrt(np.sum((res_d.x - a_current)**2))

        bond_data[R] = {'V': V_scf, 'V_bare': V_bare, 'a': res_d.x.copy()}
        a_current = res_d.x.copy()

        report(f"{R:4d} {V_scf:+12.8f} {V_bare:+12.8f} {da:8.5f} {res_d.nfev:5d}")

    report("")

    # Analysis
    R_arr = np.array(R_scan, dtype=float)
    V_arr = np.array([bond_data[R]['V'] for R in R_scan])

    i_min = np.argmin(V_arr)
    V_at_inf = V_arr[0]  # R=80

    report(f"Results for {n_modes} modes:")
    report(f"  V(R=80) = {V_at_inf:+.8f} (should ≈ 0)")
    report(f"  V_min = {V_arr[i_min]:+.8f} at R = {R_arr[i_min]:.0f}")

    if V_arr[i_min] < -1e-6:
        report(f"  D_e = {-V_arr[i_min]:.8f}")
        report(f"  *** BONDING DETECTED ***")
    elif V_arr[i_min] < 0:
        report(f"  Weak attraction: {V_arr[i_min]:+.8f}")
    else:
        report(f"  No bonding (all V > 0)")

    # Smoothness
    jumps = np.max(np.abs(np.diff(V_arr)))
    report(f"  Max step-to-step jump: {jumps:.6f}")
    report(f"  Smooth: {'YES' if jumps < 0.05 else 'NO'}")

    all_results[n_modes] = {
        'E_single': E_single_opt, 'a_single': a_single_opt,
        'bond_data': bond_data, 'V_at_inf': V_at_inf,
        'V_min': V_arr[i_min], 'R_eq': R_arr[i_min],
    }
    report("")

# ============================================================
# COMPARISON ACROSS MODE COUNTS
# ============================================================
report("COMPARISON ACROSS MODE COUNTS")
report("=" * 70)
report("")

report(f"{'n_modes':>7} {'E_single':>12} {'V(R=80)':>12} {'V_min':>12} "
       f"{'R_eq':>5} {'D_e':>12} {'smooth':>7}")
report("-" * 70)

for nm in sorted(all_results.keys()):
    r = all_results[nm]
    De = -r['V_min'] if r['V_min'] < 0 else 0
    bd = r['bond_data']
    V_arr = np.array([bd[R]['V'] for R in R_scan])
    jumps = np.max(np.abs(np.diff(V_arr)))
    smooth = "YES" if jumps < 0.05 else "NO"

    report(f"  {nm:5d} {r['E_single']:12.8f} {r['V_at_inf']:+12.8f} "
           f"{r['V_min']:+12.8f} {r['R_eq']:5.0f} {De:12.8f} {smooth:>7}")

report("")

# Full bond curves side by side
report("BOND CURVES V(R) — ALL MODE COUNTS:")
header = f"{'R':>4}"
for nm in sorted(all_results.keys()):
    header += f"  {'V('+str(nm)+'m)':>12}"
header += f"  {'V_bare':>12}"
report(header)
report("-" * (4 + (len(all_results) + 1) * 14))

for R in R_scan:
    line = f"{R:4d}"
    for nm in sorted(all_results.keys()):
        v = all_results[nm]['bond_data'][R]['V']
        line += f"  {v:+12.8f}"
    vb = all_results[1]['bond_data'][R]['V_bare']  # bare is same for all
    line += f"  {vb:+12.8f}"
    report(line)

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
