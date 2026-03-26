"""
Bare Hessian Multi-Mode Bond — The Right Framework
=====================================================
The Hessian eigenvalues ARE the breather mode energies.
No need to add breather fields to the background — the eigenstates
are already the breathers.

Method:
  1. Single kink well → find ALL bound eigenvalues (= breather modes)
  2. Two kink wells at R → each eigenvalue splits into bonding/antibonding
  3. Bond energy for each mode = bonding eigenvalue - single eigenvalue
  4. Total bond energy = sum over OCCUPIED modes
  5. σ/π ratio = ratio of mode splittings (NOT angular geometry)

This is what bond_3d_emerge.py did for 1 mode. Here we track ALL modes.
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "bare_hessian_multimode_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("BARE HESSIAN MULTI-MODE BOND MODEL")
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
kw = 3

def kink_field(x, pos):
    return (4.0/PI) * (np.arctan(np.exp(x - pos + kw/2))
                      - np.arctan(np.exp(x - pos - kw/2)))

def build_hessian(phi):
    diag = 2.0 + np.cos(PI * phi)
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()

# ============================================================
# PART 1: SINGLE WELL — ALL BOUND STATES
# ============================================================
report("PART 1: SINGLE KINK WELL — ALL BREATHER MODES")
report("-" * 55)

phi_single = kink_field(x, center)
H_single = build_hessian(phi_single)
n_eig = 20
evals_s, evecs_s = eigsh(H_single, k=n_eig, which='SM')
idx = np.argsort(evals_s)
evals_s = evals_s[idx]

n_bound = np.sum(evals_s < 1.0)
report(f"Bound states (below mass gap omega^2=1): {n_bound}")
report("")

report(f"{'mode':>4} {'omega^2':>12} {'omega':>10} {'GWT omega':>10} {'err%':>8} {'status':>8}")
report("-" * 58)

single_evals = []
for i in range(min(n_eig, 10)):
    ev = evals_s[i]
    omega = np.sqrt(abs(ev))
    gwt = np.cos((i+1) * gamma) if i < 24 else 0
    err = abs(omega - gwt) / gwt * 100 if gwt > 0 else 999
    status = "BOUND" if ev < 1.0 else "band"
    single_evals.append(ev)
    report(f"  {i:4d} {ev:12.8f} {omega:10.6f} {gwt:10.6f} {err:8.2f}% {status:>8}")

report("")

# ============================================================
# PART 2: TWO WELLS — ALL MODE SPLITTINGS VS R
# ============================================================
report("PART 2: ALL MODE SPLITTINGS VS SEPARATION R")
report("-" * 55)
report("Each single-well eigenvalue splits into bonding + antibonding.")
report("Splitting = tunnel coupling strength for that mode.")
report("")

R_scan = list(range(4, 52, 2))
n_track = min(n_bound, 7)  # track up to 7 bound modes

# Store: splitting[mode_i][R] = antibonding - bonding eigenvalue
splittings = {i: {} for i in range(n_track)}
bond_shifts = {i: {} for i in range(n_track)}

# Header
header = f"{'R':>4}"
for i in range(n_track):
    header += f"  {'spl_'+str(i):>11}"
report(header)
report("-" * (4 + n_track * 13))

for R in R_scan:
    pos_A = center - R // 2
    pos_B = center + R // 2
    phi_double = kink_field(x, pos_A) + kink_field(x, pos_B)

    H_double = build_hessian(phi_double)
    # Need 2× as many eigenvalues (each mode doubles)
    evals_d, _ = eigsh(H_double, k=min(2 * n_track + 4, N - 2), which='SM')
    evals_d = np.sort(evals_d)

    # Pair up eigenvalues: modes 0,1 are bonding/antibonding of mode 0, etc.
    line = f"{R:4d}"
    for i in range(n_track):
        if 2*i + 1 < len(evals_d):
            e_bond = evals_d[2*i]
            e_anti = evals_d[2*i + 1]
            split = e_anti - e_bond
            shift = e_bond - single_evals[i]  # bonding shift

            splittings[i][R] = split
            bond_shifts[i][R] = shift
            line += f"  {split:11.8f}"
        else:
            line += f"  {'---':>11}"

    report(line)

report("")

# ============================================================
# PART 3: BOND SHIFT (V(R)) FOR EACH MODE
# ============================================================
report("PART 3: BONDING SHIFT V_n(R) = E_bond_n(R) - E_single_n")
report("-" * 55)
report("Negative = attractive (bonding). This is the energy gain per mode.")
report("")

header = f"{'R':>4}"
for i in range(n_track):
    header += f"  {'V_'+str(i):>11}"
report(header)
report("-" * (4 + n_track * 13))

for R in R_scan:
    line = f"{R:4d}"
    for i in range(n_track):
        v = bond_shifts[i].get(R, 0)
        line += f"  {v:+11.8f}"
    report(line)

report("")

# ============================================================
# PART 4: TOTAL BOND ENERGY FOR DIFFERENT OCCUPATIONS
# ============================================================
report("PART 4: TOTAL BOND ENERGY BY OCCUPATION")
report("-" * 55)
report("V_total(R) = sum of V_n(R) over occupied modes.")
report("")
report("Occupations (like electron configurations):")
report("  1 mode:  H-like (1 electron)")
report("  2 modes: He-like (2 electrons)")
report("  3 modes: Li/Be-like")
report("  4 modes: C-like (4 valence)")
report("  All:     all bound modes")
report("")

occupations = {
    '1 mode': [0],
    '2 modes': [0, 1],
    'all': list(range(n_track)),
}
if n_track >= 3:
    occupations['3 modes'] = [0, 1, 2]

header = f"{'R':>4}"
for name in occupations:
    header += f"  {name:>12}"
report(header)
report("-" * (4 + len(occupations) * 14))

total_curves = {}
for name in occupations:
    total_curves[name] = {}

for R in R_scan:
    line = f"{R:4d}"
    for name, occ in occupations.items():
        V_total = sum(bond_shifts[i].get(R, 0) for i in occ)
        total_curves[name][R] = V_total
        line += f"  {V_total:+12.8f}"
    report(line)

report("")

# ============================================================
# PART 5: MORSE WELL ANALYSIS FOR EACH OCCUPATION
# ============================================================
report("PART 5: MORSE WELL ANALYSIS")
report("-" * 55)

for name, occ in occupations.items():
    V_arr = np.array([total_curves[name][R] for R in R_scan])
    R_arr = np.array(R_scan, dtype=float)

    i_min = np.argmin(V_arr)
    D_e = -V_arr[i_min] if V_arr[i_min] < 0 else 0
    R_eq = R_arr[i_min]

    # Decay rate
    a_morse = 0
    if D_e > 1e-8:
        mask = (V_arr < -1e-10) & (R_arr > R_eq)
        if np.sum(mask) > 2:
            coeffs = np.polyfit(R_arr[mask], np.log(-V_arr[mask]), 1)
            a_morse = -coeffs[0]

    report(f"  {name:>12}: R_eq={R_eq:5.0f}, D_e={D_e:.8f}, Morse_a={a_morse:.4f}")

report("")

# ============================================================
# PART 6: SPLITTING RATIOS — THE REAL σ/π
# ============================================================
report("PART 6: SPLITTING RATIOS — σ/π FROM EIGENVALUE STRUCTURE")
report("-" * 55)
report("The ratio split_1/split_0 is the REAL pi/sigma coupling ratio.")
report("This should relate to W_PI = cos(pi/d) = 0.5")
report("")

report(f"{'R':>4} {'split_0':>12} {'split_1':>12} {'split_2':>12} "
       f"{'s1/s0':>8} {'s2/s0':>8}")
report("-" * 58)

for R in R_scan:
    s0 = splittings[0].get(R, 0)
    s1 = splittings[1].get(R, 0) if 1 < n_track else 0
    s2 = splittings[2].get(R, 0) if 2 < n_track else 0

    r10 = s1 / s0 if s0 > 1e-12 else 0
    r20 = s2 / s0 if s0 > 1e-12 else 0

    if s0 > 1e-10:  # only show where splitting is measurable
        report(f"{R:4d} {s0:12.8f} {s1:12.8f} {s2:12.8f} {r10:8.4f} {r20:8.4f}")

report("")

# Decay rates for each mode
report("SPLITTING DECAY RATES:")
report(f"{'mode':>4} {'decay_rate':>11} {'ratio_to_m0':>12} {'eps_n':>10} {'width':>8}")
report("-" * 50)

decay_rates = {}
for i in range(n_track):
    R_arr_i = np.array([R for R in R_scan if splittings[i].get(R, 0) > 1e-12])
    S_arr_i = np.array([splittings[i][R] for R in R_arr_i])

    if len(R_arr_i) > 3:
        coeffs = np.polyfit(R_arr_i, np.log(S_arr_i + 1e-30), 1)
        decay = -coeffs[0]
    else:
        decay = 0

    decay_rates[i] = decay
    eps_n = np.sin((i+1) * gamma)
    ratio = decay / decay_rates[0] if decay_rates[0] > 0 else 0

    report(f"  {i:4d} {decay:11.6f} {ratio:12.4f} {eps_n:10.6f} {1/eps_n:8.2f}")

report("")

# The σ/π ratio should be related to the decay rate ratio
if decay_rates[0] > 0 and len(decay_rates) > 1 and decay_rates[1] > 0:
    report("KEY RATIO:")
    report(f"  Mode 0 (σ) decay rate: {decay_rates[0]:.6f}")
    report(f"  Mode 1 (π) decay rate: {decay_rates.get(1, 0):.6f}")
    report(f"  Ratio mode1/mode0: {decay_rates.get(1,0)/decay_rates[0]:.6f}")
    report(f"  V8 W_PI = cos(π/d) = {np.cos(PI/d):.6f}")
    report("")

    # At equilibrium R, the splitting ratio gives W_PI(R_eq)
    R_eq_idx = np.argmin([total_curves['1 mode'].get(R, 999) for R in R_scan])
    R_eq = R_scan[R_eq_idx]
    s0_eq = splittings[0].get(R_eq, 0)
    s1_eq = splittings[1].get(R_eq, 0) if 1 < n_track else 0
    if s0_eq > 1e-12:
        report(f"  At R_eq={R_eq}: split_1/split_0 = {s1_eq/s0_eq:.6f}")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report(f"Single kink well: {n_bound} bound states")
report("")

for name in occupations:
    V_arr = np.array([total_curves[name][R] for R in R_scan])
    V_min = np.min(V_arr)
    R_eq = R_scan[np.argmin(V_arr)]
    De = -V_min if V_min < 0 else 0
    report(f"  {name:>12}: D_e = {De:.6f} at R = {R_eq}")

report("")
report("The eigenvalue splittings ARE the bond — no perturbation needed.")
report("Each mode's splitting = its contribution to the bond energy.")
report("The σ/π ratio comes from the DECAY RATE ratio of the splittings.")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
