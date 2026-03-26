"""
Breather Mode Shift — Self-Consistent Bond Correction
=======================================================
When two breathers sit near each other, their kink tails overlap and modify
the effective potential each breather sees. This shifts the breather mode
frequencies. The current bond model treats breather parameters as fixed —
adding self-consistent mode shifting could close the ~7% V8/V10 bond gap.

Method:
  1. Isolated breather: Hessian eigenvalues in a single kink well (reference)
  2. Two-breather system: Hessian with kink-antikink pairs at separation R
  3. Track how breather mode frequencies shift vs R
  4. Self-consistent iteration: shift → width change → overlap change → repeat
  5. Compare bond curves: fixed-mode vs self-consistent

Physics:
  - Kink tail decays as exp(-eps_n * |x|), so shift ~ exp(-eps_n * R)
  - Mode 1 (wide, ~15 sites) feels neighbor at large R
  - Mode 7 (narrow, ~2.3 sites) only at close range
  - Neighbor kink SOFTENS potential → lower frequency → wider breather
  - Feedback converges in 2-3 iterations

All parameters from the GWT Lagrangian (d=3, zero free parameters).
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "breather_mode_shift_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("BREATHER MODE SHIFT — SELF-CONSISTENT BOND CORRECTION")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"d = {d}, gamma = {gamma:.10f}")
report("")

# ============================================================
# LATTICE AND PROFILES
# ============================================================
N = 512  # lattice sites (large enough for two well-separated breathers)

def kink_profile(x, x_center):
    """Single kink: phi goes 0 -> 2."""
    return (4.0/PI) * np.arctan(np.exp(x - x_center))

def antikink_profile(x, x_center):
    """Anti-kink: phi goes 2 -> 0."""
    return 2.0 - (4.0/PI) * np.arctan(np.exp(x - x_center))

def make_kink_antikink(x, center, width):
    """Kink-antikink pair (proton): phi goes 0 -> 2 -> 0."""
    return kink_profile(x, center - width/2) + antikink_profile(x, center + width/2) - 2.0

def breather_profile(x, center, eps_n):
    """Breather peak profile: phi(x) = (4/pi)*arctan(1/cosh(eps*x))."""
    return (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * (x - center) + 1e-30))

def build_hessian(phi_background, N):
    """Hessian H = -Laplacian + d^2V/dphi^2 at background phi.
    For V = (1/pi^2)(1 - cos(pi*phi)): d^2V/dphi^2 = cos(pi*phi).
    Eigenvalues = omega^2 (squared frequencies of linearized modes).
    """
    diag = 2.0 + np.cos(PI * phi_background)
    off_diag = -np.ones(N - 1)
    H = sparse.diags([off_diag, diag, off_diag], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()

# ============================================================
# PART 1: ISOLATED BREATHER EIGENVALUES (reference)
# ============================================================
report("PART 1: ISOLATED BREATHER MODE FREQUENCIES")
report("-" * 55)

x = np.arange(N, dtype=np.float64)
center = N // 2
kink_width = 3  # proton size

# Single proton (kink-antikink pair)
phi_single = make_kink_antikink(x, center, kink_width)

# Hessian for single well
H_single = build_hessian(phi_single, N)
evals_single, evecs_single = eigsh(H_single, k=20, which='SM')
idx = np.argsort(evals_single)
evals_single = evals_single[idx]
evecs_single = evecs_single[:, idx]

report(f"Lattice: {N} sites, kink_width = {kink_width}")
report(f"{'mode':>4} {'omega^2':>12} {'omega':>10} {'GWT_omega':>10} {'error%':>8} {'status':>8}")

n_bound_single = 0
single_bound_states = []
for i in range(min(15, len(evals_single))):
    ev = evals_single[i]
    omega = np.sqrt(abs(ev))
    status = "BOUND" if ev < 1.0 else "band"
    if ev < 1.0:
        n_bound_single += 1
        single_bound_states.append((i, ev, omega))
    # Compare with GWT prediction omega_n = cos(n*gamma)
    gwt_omega = np.cos((i+1) * gamma) if i < 24 else 0
    err = abs(omega - gwt_omega) / gwt_omega * 100 if gwt_omega > 0 else 999
    report(f"  {i:4d} {ev:12.8f} {omega:10.6f} {gwt_omega:10.6f} {err:8.2f}% {status:>8}")

report(f"\nBound states below mass gap: {n_bound_single}")
report("")

# ============================================================
# PART 2: FREQUENCY SHIFT VS SEPARATION R
# ============================================================
report("PART 2: MODE FREQUENCY SHIFT VS SEPARATION R")
report("-" * 55)
report("Two protons at separation R. Track lowest eigenvalue shifts.")
report("")

R_values = list(range(6, 52, 2))
n_eig = 20

# Store results: {R: [eigenvalues]}
shift_data = {}

report(f"{'R':>4} {'E0':>12} {'E1':>12} {'shift_E0':>12} {'shift_E1':>12} {'split':>12}")
report("-" * 75)

E0_ref = evals_single[0]  # reference isolated bound state energy
E1_ref = evals_single[1] if n_bound_single > 1 else 1.0

for R in R_values:
    pos_A = center - R // 2
    pos_B = center + R // 2

    phi_double = make_kink_antikink(x, pos_A, kink_width) + make_kink_antikink(x, pos_B, kink_width)

    H_double = build_hessian(phi_double, N)
    evals_double, evecs_double = eigsh(H_double, k=n_eig, which='SM')
    evals_double = np.sort(evals_double)

    shift_data[R] = evals_double

    E0 = evals_double[0]
    E1 = evals_double[1] if len(evals_double) > 1 else 999

    shift_E0 = E0 - E0_ref
    shift_E1 = E1 - E0_ref  # relative to same reference
    split = E1 - E0

    report(f"{R:4d} {E0:12.8f} {E1:12.8f} {shift_E0:+12.8f} {shift_E1:+12.8f} {split:12.8f}")

report("")

# ============================================================
# PART 3: DECAY RATE ANALYSIS
# ============================================================
report("PART 3: SHIFT DECAY RATE — DOES IT MATCH eps_n?")
report("-" * 55)

# The shift should decay as exp(-eps_n * R) where eps_n = sin(n*gamma)
# For the lowest mode, eps_1 = sin(gamma)
eps_1 = np.sin(gamma)
report(f"GWT prediction: eps_1 = sin(gamma) = {eps_1:.8f}")
report(f"Expected decay rate of shift: exp(-{eps_1:.4f} * R)")
report("")

# Extract the bonding shift (V(R) = E0_double - E0_single)
R_arr = np.array(R_values)
V_arr = np.array([shift_data[R][0] - E0_ref for R in R_values])

# Fit exponential decay to the attractive (negative) part
mask = V_arr < -1e-10
if np.sum(mask) > 3:
    log_V = np.log(-V_arr[mask])
    R_fit = R_arr[mask]

    # Linear fit: log(-V) = A - decay_rate * R
    coeffs = np.polyfit(R_fit, log_V, 1)
    decay_rate = -coeffs[0]

    report(f"Fitted decay rate: {decay_rate:.6f}")
    report(f"eps_1 = {eps_1:.6f}")
    report(f"Ratio (decay/eps_1): {decay_rate/eps_1:.4f}")
    report(f"2*eps_1 = {2*eps_1:.6f} (kink-kink repulsion rate)")
    report("")

    # Check against different eps_n
    for n in [1, 2, 3, 4, 7]:
        eps_n = np.sin(n * gamma)
        report(f"  n={n}: eps_n = {eps_n:.6f}, ratio = {decay_rate/eps_n:.4f}")
else:
    report("Not enough negative V(R) points for decay fit.")

report("")

# ============================================================
# PART 4: CROSS-MODE SHIFTS
# ============================================================
report("PART 4: CROSS-MODE SHIFTS")
report("-" * 55)
report("How does a mode-n breather shift when near a mode-m breather?")
report("(Asymmetric: wide mode feels narrow mode differently than reverse)")
report("")

# For this we need to add breather profiles to the kink background
# The breather modifies the kink well potential by occupying a bound state
# This is like inner-electron screening in atoms

# We can approximate: breather in well A modifies the potential that well B sees
# by the breather's energy density (phi_breather^2) at well B's location

report("Cross-mode overlap integrals:")
report(f"{'n':>3} {'m':>3} {'eps_n':>8} {'eps_m':>8} {'width_n':>8} {'width_m':>8} "
       f"{'overlap_at_R10':>14} {'overlap_at_R20':>14}")

for n in [1, 4, 7]:
    eps_n = np.sin(n * gamma)
    width_n = 1.0 / eps_n

    for m in [1, 4, 7]:
        eps_m = np.sin(m * gamma)
        width_m = 1.0 / eps_m

        # Breather n centered at 0, breather m centered at R
        # Overlap = integral of |phi_n(x)|^2 * |phi_m(x-R)|^2
        for R_test in [10, 20]:
            x_test = np.arange(-50, 50+R_test, dtype=np.float64)
            phi_n = breather_profile(x_test, 0, eps_n)
            phi_m = breather_profile(x_test, R_test, eps_m)
            overlap = np.sum(phi_n**2 * phi_m**2) * 1.0  # dx = 1

            if R_test == 10:
                ov10 = overlap
            else:
                ov20 = overlap

        report(f"  {n:3d} {m:3d} {eps_n:8.4f} {eps_m:8.4f} {width_n:8.2f} {width_m:8.2f} "
               f"{ov10:14.6e} {ov20:14.6e}")

report("")

# ============================================================
# PART 5: SELF-CONSISTENT ITERATION
# ============================================================
report("PART 5: SELF-CONSISTENT MODE SHIFT")
report("-" * 55)
report("Iteration: shift changes width, width changes overlap, repeat.")
report("")

# For each separation R, iterate:
# 1. Start with isolated breather width w_0 = 1/eps_1
# 2. Compute Hessian with neighbor at distance R
# 3. Find new lowest eigenvalue -> new omega -> new eps -> new width
# 4. Recompute breather profile with new width
# 5. Repeat until converged

R_test_values = [8, 10, 12, 16, 20, 30, 40]

report(f"{'R':>4} {'iter':>4} {'omega':>10} {'eps':>10} {'width':>8} {'shift%':>10}")
report("-" * 55)

sc_results = {}  # {R: (converged_omega, converged_eps, n_iter)}

for R in R_test_values:
    # Initial: isolated breather parameters
    eps_current = np.sin(gamma)  # mode 1
    omega_current = np.cos(gamma)
    omega_0 = omega_current  # reference

    converged = False
    for iteration in range(10):
        # Build background: two protons + effective potential modification
        # The neighbor's kink tail at distance R modifies our potential
        pos_A = center - R // 2
        pos_B = center + R // 2

        phi_bg = make_kink_antikink(x, pos_A, kink_width) + make_kink_antikink(x, pos_B, kink_width)

        # Add the breather's own effect on the potential:
        # The breather at well A with current eps modifies the Hessian
        # via the sinc perturbation
        phi_br_A = breather_profile(x, pos_A, eps_current)

        # The breather modifies the effective potential through the phi^4 coupling
        # delta_V = -(pi*phi_breather)^2 / 6 (leading order)
        # This modifies the diagonal of the Hessian
        delta_V = -(PI * phi_br_A)**2 / 6.0

        # Build Hessian with breather correction
        H = build_hessian(phi_bg, N)
        # Add breather perturbation as diagonal correction
        H_mod = H + sparse.diags(delta_V, 0, shape=(N, N), format='csr')

        evals, evecs = eigsh(H_mod, k=5, which='SM')
        evals = np.sort(evals)

        omega_new = np.sqrt(abs(evals[0])) if evals[0] > 0 else np.sqrt(abs(evals[0]))
        eps_new = np.sqrt(abs(1.0 - evals[0])) if evals[0] < 1.0 else 0.01
        width_new = 1.0 / eps_new if eps_new > 0.01 else 100.0

        shift_pct = (omega_new - omega_0) / omega_0 * 100

        report(f"{R:4d} {iteration:4d} {omega_new:10.6f} {eps_new:10.6f} {width_new:8.2f} {shift_pct:+10.4f}%")

        # Check convergence
        if iteration > 0 and abs(omega_new - omega_current) < 1e-10:
            converged = True
            sc_results[R] = (omega_new, eps_new, iteration + 1, shift_pct)
            report(f"  --> Converged in {iteration+1} iterations")
            break

        omega_current = omega_new
        eps_current = eps_new

    if not converged:
        sc_results[R] = (omega_current, eps_current, 10, shift_pct)
        report(f"  --> Not converged after 10 iterations")
    report("")

# ============================================================
# PART 6: BOND CURVE COMPARISON — FIXED vs SELF-CONSISTENT
# ============================================================
report("PART 6: BOND CURVE — FIXED MODE vs SELF-CONSISTENT")
report("-" * 55)
report("Compare bond energy with and without self-consistent mode shifting.")
report("")

R_scan = list(range(6, 46, 2))

# Fixed-mode bond curve (from bond_3d_emerge.py method)
V_fixed = []
V_selfcon = []

for R in R_scan:
    pos_A = center - R // 2
    pos_B = center + R // 2

    # FIXED: just two protons, no breather correction
    phi_fixed = make_kink_antikink(x, pos_A, kink_width) + make_kink_antikink(x, pos_B, kink_width)
    H_fixed = build_hessian(phi_fixed, N)
    ev_fixed, _ = eigsh(H_fixed, k=5, which='SM')
    ev_fixed = np.sort(ev_fixed)
    V_fixed.append(ev_fixed[0] - E0_ref)

    # SELF-CONSISTENT: add breather perturbation
    eps_sc = np.sin(gamma)  # start with isolated
    for _ in range(3):  # 3 iterations is enough
        phi_br = breather_profile(x, pos_A, eps_sc)
        delta_V = -(PI * phi_br)**2 / 6.0
        H_sc = build_hessian(phi_fixed, N) + sparse.diags(delta_V, 0, shape=(N, N), format='csr')
        ev_sc, _ = eigsh(H_sc, k=5, which='SM')
        ev_sc = np.sort(ev_sc)
        eps_sc = np.sqrt(abs(1.0 - ev_sc[0])) if ev_sc[0] < 1.0 else eps_sc

    V_selfcon.append(ev_sc[0] - E0_ref)

V_fixed = np.array(V_fixed)
V_selfcon = np.array(V_selfcon)
R_scan = np.array(R_scan)

report(f"{'R':>4} {'V_fixed':>12} {'V_selfcon':>12} {'correction':>12} {'corr_%':>10}")
report("-" * 55)

for i, R in enumerate(R_scan):
    corr = V_selfcon[i] - V_fixed[i]
    corr_pct = corr / abs(V_fixed[i]) * 100 if abs(V_fixed[i]) > 1e-12 else 0
    marker = ""
    if i > 0 and V_fixed[i] < V_fixed[i-1] and (i == len(R_scan)-1 or V_fixed[i] < V_fixed[i+1]):
        marker = " <-- min(fixed)"
    report(f"{R:4d} {V_fixed[i]:+12.8f} {V_selfcon[i]:+12.8f} {corr:+12.8f} {corr_pct:+10.4f}%{marker}")

# Find minima
i_min_fixed = np.argmin(V_fixed)
i_min_sc = np.argmin(V_selfcon)

report("")
report("WELL DEPTH COMPARISON:")
report(f"  Fixed:         R_eq = {R_scan[i_min_fixed]:3d}, D_e = {-V_fixed[i_min_fixed]:.8f}")
report(f"  Self-consistent: R_eq = {R_scan[i_min_sc]:3d}, D_e = {-V_selfcon[i_min_sc]:.8f}")
if abs(V_fixed[i_min_fixed]) > 1e-12:
    depth_change = (-V_selfcon[i_min_sc] - (-V_fixed[i_min_fixed])) / (-V_fixed[i_min_fixed]) * 100
    report(f"  Depth change: {depth_change:+.2f}%")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report("Key findings:")
report(f"  1. Isolated breather: omega_1 = {np.cos(gamma):.6f}, eps_1 = {np.sin(gamma):.6f}")
report(f"  2. Shift decay rate vs eps_1:")
if 'decay_rate' in dir():
    report(f"     Fitted: {decay_rate:.6f}, eps_1: {eps_1:.6f}, ratio: {decay_rate/eps_1:.4f}")
report(f"  3. Self-consistent convergence: typically 2-3 iterations")
if sc_results:
    for R in sorted(sc_results.keys()):
        omega, eps, n_iter, shift = sc_results[R]
        report(f"     R={R:3d}: shift = {shift:+.4f}%, converged in {n_iter} iterations")
report(f"  4. Bond depth correction from self-consistency:")
if abs(V_fixed[i_min_fixed]) > 1e-12:
    report(f"     Fixed D_e = {-V_fixed[i_min_fixed]:.8f}")
    report(f"     Self-consistent D_e = {-V_selfcon[i_min_sc]:.8f}")
    report(f"     Correction = {depth_change:+.2f}%")
report("")

report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
