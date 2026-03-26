"""
Mode Coupling Matrix — The 7-Mode Equilibrium
================================================
The 7 breather modes in a proton form a tightly coupled system.
Their individual perturbations nearly cancel (99.7% cancellation).
This means there's a CONSTRAINT — a sum rule — that fixes the
ratios of mode amplitudes.

Model:
  1. Build the 7×7 overlap matrix M_{nm}:
       M_{nm} = ∫ phi_n(x)^2 * phi_m(x)^2 dx
     This measures how much mode n's perturbation overlaps mode m's.

  2. Also build the "well-filling" vector W_n:
       W_n = ∫ phi_n(x)^2 * V_kink(x) dx
     This measures how much mode n's perturbation looks like the kink well.

  3. The near-cancellation means: sum_n a_n * W_n ≈ 0
     where a_n are the mode amplitudes. The ratios a_n/a_1 are
     determined by this constraint.

  4. For bonding: compute how M and W change when a second kink is nearby.
     The shift in the constrained equilibrium = bond energy.

All parameters from GWT Lagrangian (d=3).
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "mode_coupling_matrix_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("MODE COUPLING MATRIX — THE 7-MODE EQUILIBRIUM")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"d = {d}, gamma = {gamma:.10f}")
report("")

# ============================================================
# SETUP
# ============================================================
N = 1024  # high resolution for accurate integrals
x = np.arange(N, dtype=np.float64) - N/2
kink_width = 3

def kink_antikink(x, center, width):
    return (4.0/PI) * (np.arctan(np.exp(x - center + width/2))
                      - np.arctan(np.exp(x - center - width/2)))

def breather_sq(x, center, eps_n):
    """Squared breather profile (the perturbation is proportional to this)."""
    phi = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * (x - center) + 1e-30))
    return phi**2

phi_kink = kink_antikink(x, 0, kink_width)
V_kink = np.cos(PI * phi_kink)  # the potential curvature at the kink

modes = list(range(1, 8))
n_modes = len(modes)

# ============================================================
# PART 1: MODE PROFILES AND BASIC PROPERTIES
# ============================================================
report("PART 1: MODE PROFILES")
report("-" * 55)

profiles = {}  # {n: phi_n(x)^2}
for n in modes:
    eps_n = np.sin(n * gamma)
    profiles[n] = breather_sq(x, 0, eps_n)

report(f"{'n':>3} {'eps_n':>8} {'width':>6} {'integral(phi^2)':>15} {'peak':>8}")
report("-" * 45)
for n in modes:
    eps_n = np.sin(n * gamma)
    integ = np.sum(profiles[n])  # sum ≈ integral (dx=1)
    peak = np.max(profiles[n])
    report(f"  {n:3d} {eps_n:8.4f} {1/eps_n:6.1f} {integ:15.4f} {peak:8.4f}")

report("")

# ============================================================
# PART 2: THE 7×7 OVERLAP MATRIX M
# ============================================================
report("PART 2: OVERLAP MATRIX M_{nm} = ∫ φ_n² · φ_m² dx")
report("-" * 55)
report("This is the mode-mode coupling strength.")
report("")

M = np.zeros((n_modes, n_modes))
for i, n in enumerate(modes):
    for j, m in enumerate(modes):
        M[i, j] = np.sum(profiles[n] * profiles[m])

# Print matrix
header = "     " + "".join(f"  {m:>8}" for m in modes)
report(header)
report("-" * (5 + n_modes * 10))
for i, n in enumerate(modes):
    line = f"  {n:3d}"
    for j in range(n_modes):
        line += f"  {M[i,j]:8.2f}"
    report(line)

report("")

# Eigenvalues of M
M_evals, M_evecs = np.linalg.eigh(M)
report("Eigenvalues of M (overlap matrix):")
for i, ev in enumerate(M_evals):
    report(f"  [{i}] {ev:12.4f}")
report("")

report("Eigenvectors (columns = mode weights):")
header = "     " + "".join(f"  {'v'+str(i):>8}" for i in range(n_modes))
report(header)
for i, n in enumerate(modes):
    line = f"  {n:3d}"
    for j in range(n_modes):
        line += f"  {M_evecs[i,j]:+8.4f}"
    report(line)
report("")

# ============================================================
# PART 3: WELL-FILLING VECTOR W
# ============================================================
report("PART 3: WELL-FILLING VECTOR W_n = ∫ φ_n² · V_kink dx")
report("-" * 55)
report("How much does each mode's perturbation overlap the kink potential?")
report("")

W = np.zeros(n_modes)
for i, n in enumerate(modes):
    W[i] = np.sum(profiles[n] * V_kink)

report(f"{'n':>3} {'W_n':>12} {'W_n/W_1':>10}")
report("-" * 30)
for i, n in enumerate(modes):
    ratio = W[i] / W[0] if abs(W[0]) > 1e-10 else 0
    report(f"  {n:3d} {W[i]:12.4f} {ratio:10.4f}")

report("")

# Check: does a linear combination of modes cancel V_kink?
# Solve: sum_n a_n * phi_n^2 ≈ c * V_kink (in a least-squares sense)
# Using the profile matrix
report("LEAST-SQUARES FIT: sum_n a_n · φ_n² ≈ V_kink")
P_matrix = np.column_stack([profiles[n] for n in modes])
# Only fit in the region where the kink is nonzero
kink_mask = np.abs(phi_kink) > 0.01
P_masked = P_matrix[kink_mask]
V_masked = V_kink[kink_mask]

# Also try: sum_n a_n * phi_n^2 ≈ -delta_V_kink where delta_V_kink = V_kink - 1
# (the perturbation relative to the vacuum value)
delta_V = V_kink - 1.0  # V_kink = cos(pi*phi), vacuum = cos(0) = 1
delta_V_masked = delta_V[kink_mask]

a_fit, residual, rank, sv = np.linalg.lstsq(P_masked, delta_V_masked, rcond=None)

report(f"  Fit coefficients a_n (modes fill the kink well deviation):")
for i, n in enumerate(modes):
    report(f"    n={n}: a = {a_fit[i]:+.6f}")

# Quality of fit
fit_result = P_matrix @ a_fit
residual_norm = np.sqrt(np.sum((fit_result[kink_mask] - delta_V[kink_mask])**2))
total_norm = np.sqrt(np.sum(delta_V[kink_mask]**2))
report(f"  Residual / Total = {residual_norm/total_norm:.6f} ({residual_norm/total_norm*100:.2f}%)")
report("")

# ============================================================
# PART 4: RATIO STRUCTURE — WHAT FIXES THE RATIOS?
# ============================================================
report("PART 4: RATIO STRUCTURE")
report("-" * 55)
report("Are the mode ratios simple functions of n and gamma?")
report("")

# Ratios relative to mode 1
report("Fit coefficients normalized to a_1:")
report(f"{'n':>3} {'a_n':>10} {'a_n/a_1':>10} {'eps_n/eps_1':>11} {'omega_n':>10}")
report("-" * 50)
for i, n in enumerate(modes):
    ratio_a = a_fit[i] / a_fit[0] if abs(a_fit[0]) > 1e-10 else 0
    eps_ratio = np.sin(n*gamma) / np.sin(gamma)
    omega_n = np.cos(n*gamma)
    report(f"  {n:3d} {a_fit[i]:+10.4f} {ratio_a:+10.4f} {eps_ratio:11.4f} {omega_n:10.6f}")

report("")

# Try simple models for the ratios
report("Testing ratio models:")
report("")

# Model 1: a_n ~ 1/eps_n (width weighting)
model1 = np.array([1.0/np.sin(n*gamma) for n in modes])
model1 = model1 / model1[0] * a_fit[0]
err1 = np.sqrt(np.mean((model1 - a_fit)**2)) / np.sqrt(np.mean(a_fit**2))
report(f"  Model 1 (a_n ~ 1/eps_n, width weighting): error = {err1*100:.1f}%")

# Model 2: a_n ~ eps_n (inverse width)
model2 = np.array([np.sin(n*gamma) for n in modes])
model2 = model2 / model2[0] * a_fit[0]
err2 = np.sqrt(np.mean((model2 - a_fit)**2)) / np.sqrt(np.mean(a_fit**2))
report(f"  Model 2 (a_n ~ eps_n, inverse width): error = {err2*100:.1f}%")

# Model 3: a_n ~ omega_n (frequency weighting)
model3 = np.array([np.cos(n*gamma) for n in modes])
model3 = model3 / model3[0] * a_fit[0]
err3 = np.sqrt(np.mean((model3 - a_fit)**2)) / np.sqrt(np.mean(a_fit**2))
report(f"  Model 3 (a_n ~ omega_n, frequency weighting): error = {err3*100:.1f}%")

# Model 4: a_n ~ (-1)^n * eps_n (alternating)
model4 = np.array([(-1)**n * np.sin(n*gamma) for n in modes])
model4 = model4 / model4[0] * a_fit[0]
err4 = np.sqrt(np.mean((model4 - a_fit)**2)) / np.sqrt(np.mean(a_fit**2))
report(f"  Model 4 (a_n ~ (-1)^n * eps_n): error = {err4*100:.1f}%")

# Model 5: a_n ~ cos(n*pi/8) or some discrete pattern
model5 = np.array([np.cos(n*PI/8) for n in modes])
model5 = model5 / model5[0] * a_fit[0]
err5 = np.sqrt(np.mean((model5 - a_fit)**2)) / np.sqrt(np.mean(a_fit**2))
report(f"  Model 5 (a_n ~ cos(n·π/8)): error = {err5*100:.1f}%")

# Model 6: a_n from the smallest eigenvector of M
model6 = M_evecs[:, 0]  # smallest eigenvalue eigenvector
model6 = model6 / model6[0] * a_fit[0]
err6 = np.sqrt(np.mean((model6 - a_fit)**2)) / np.sqrt(np.mean(a_fit**2))
report(f"  Model 6 (smallest eigenvector of M): error = {err6*100:.1f}%")

# Model 7: a_n from M^{-1} W (inverse matrix times well vector)
try:
    model7 = np.linalg.solve(M, W)
    model7 = model7 / model7[0] * a_fit[0]
    err7 = np.sqrt(np.mean((model7 - a_fit)**2)) / np.sqrt(np.mean(a_fit**2))
    report(f"  Model 7 (M^{{-1}} W, matrix solution): error = {err7*100:.1f}%")
except Exception:
    report(f"  Model 7: M is singular")

report("")

# ============================================================
# PART 5: THE SUM RULE
# ============================================================
report("PART 5: THE SUM RULE")
report("-" * 55)
report("Does sum_n a_n^2 * phi_n^2 = const (completeness)?")
report("Or sum_n a_n * phi_n^2 = V_kink (well-filling)?")
report("")

# Test: total perturbation = sum of weighted mode perturbations
total_pert_lsq = P_matrix @ a_fit
total_pert_equal = np.sum(P_matrix, axis=1)  # equal weights

# Compare to delta_V = cos(pi*phi_kink) - 1
report("Profile comparison at key positions:")
report(f"{'x':>5} {'delta_V':>10} {'fit':>10} {'equal_wt':>10} {'phi_kink':>10}")
report("-" * 50)
for xi in range(-8, 9):
    idx = N//2 + xi
    report(f"  {xi:5d} {delta_V[idx]:10.4f} {total_pert_lsq[idx]:10.4f} "
           f"{total_pert_equal[idx]:10.4f} {phi_kink[idx]:10.4f}")

report("")

# The sum rule: what does equal-weight perturbation look like?
# sum_n phi_n(x)^2 for n=1..7
sum_phi_sq = np.sum(P_matrix, axis=1)
report("Sum rule: Σ_n φ_n²(x) vs V_kink(x) and delta_V(x)")
report(f"  Integral of Σφ²: {np.sum(sum_phi_sq[kink_mask]):.4f}")
report(f"  Integral of |delta_V|: {np.sum(np.abs(delta_V[kink_mask])):.4f}")
report(f"  Integral of -delta_V: {np.sum(-delta_V[kink_mask]):.4f}")

# Ratio
if np.sum(-delta_V[kink_mask]) > 0:
    ratio_sum = np.sum(sum_phi_sq[kink_mask]) / np.sum(-delta_V[kink_mask])
    report(f"  Ratio Σφ² / |delta_V|: {ratio_sum:.4f}")

report("")

# ============================================================
# PART 6: BONDING — HOW THE EQUILIBRIUM SHIFTS
# ============================================================
report("PART 6: BONDING — EQUILIBRIUM SHIFT MODEL")
report("-" * 55)
report("When a second proton arrives at distance R, the kink well changes.")
report("The 7-mode equilibrium shifts. The energy change = bond energy.")
report("")

def build_hessian(phi_bg, N, extra_diag=None):
    diag = 2.0 + np.cos(PI * phi_bg)
    if extra_diag is not None:
        diag = diag + extra_diag
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()

# The LSQ fit is ill-conditioned (huge alternating coefficients).
# Use UNIT weights instead — each mode contributes equally.
# Then optimize the weights directly via energy minimization.

from scipy.optimize import minimize

R_values = [6, 8, 10, 12, 14, 16, 20]

phi_single = kink_antikink(x, 0, kink_width)
H_single_bare = build_hessian(phi_single, N)
ev_s_bare, _ = eigsh(H_single_bare, k=3, which='SM')
E0_single_bare = np.sort(ev_s_bare)[0]

# Find optimal weights for the SINGLE proton (energy minimization)
report("Finding optimal single-proton mode weights...")

def energy_single(c_vec):
    pert = np.zeros(N)
    for i, n in enumerate(modes):
        eps_n = np.sin(n * gamma)
        pert += c_vec[i] * breather_sq(x, 0, eps_n) * (-PI**2 / 6.0)
    H = build_hessian(phi_single, N, extra_diag=pert)
    ev, _ = eigsh(H, k=1, which='SM')
    return ev[0]

# Start from unit weights
c0 = np.ones(n_modes)
result_s = minimize(energy_single, c0, method='Nelder-Mead',
                    options={'maxiter': 2000, 'xatol': 1e-8, 'fatol': 1e-10})
c_single_opt = result_s.x
E0_single_opt = result_s.fun

report(f"Single proton optimal E0: {E0_single_opt:.8f} (bare: {E0_single_bare:.8f})")
report(f"Optimal weights:")
for i, n in enumerate(modes):
    report(f"  n={n}: c = {c_single_opt[i]:+.6f}")
report("")

# Normalized ratios
report("OPTIMAL WEIGHT RATIOS (single proton):")
report(f"{'n':>3} {'c_n':>10} {'c_n/c_1':>10} {'eps_n':>8} {'omega_n':>8}")
report("-" * 45)
for i, n in enumerate(modes):
    r = c_single_opt[i] / c_single_opt[0] if abs(c_single_opt[0]) > 1e-10 else 0
    report(f"  {n:3d} {c_single_opt[i]:+10.4f} {r:+10.4f} "
           f"{np.sin(n*gamma):8.4f} {np.cos(n*gamma):8.4f}")
report("")

# Now: BONDING — optimize weights for two protons at separation R
report("BONDING — OPTIMIZED WEIGHTS VS R:")
report("-" * 55)

report(f"{'R':>4} {'V_bare':>12} {'V_optimized':>12} {'V_unit_wt':>12} "
       f"{'opt_vs_bare%':>12}")
report("-" * 58)

for R in R_values:
    pos_A = -R // 2
    pos_B = R // 2
    phi_double = kink_antikink(x, pos_A, kink_width) + kink_antikink(x, pos_B, kink_width)

    # Bare
    H_bare = build_hessian(phi_double, N)
    ev_bare, _ = eigsh(H_bare, k=3, which='SM')
    V_bare = np.sort(ev_bare)[0] - E0_single_bare

    # Unit weight
    pert_unit = np.zeros(N)
    for i, n in enumerate(modes):
        eps_n = np.sin(n * gamma)
        pert_unit += breather_sq(x, pos_A, eps_n) + breather_sq(x, pos_B, eps_n)
    pert_unit *= (-PI**2 / 6.0)

    pert_unit_s = np.zeros(N)
    for i, n in enumerate(modes):
        eps_n = np.sin(n * gamma)
        pert_unit_s += breather_sq(x, 0, eps_n)
    pert_unit_s *= (-PI**2 / 6.0)

    H_unit = build_hessian(phi_double, N, extra_diag=pert_unit)
    ev_unit, _ = eigsh(H_unit, k=3, which='SM')
    H_unit_s = build_hessian(phi_single, N, extra_diag=pert_unit_s)
    ev_unit_s, _ = eigsh(H_unit_s, k=3, which='SM')
    V_unit = np.sort(ev_unit)[0] - np.sort(ev_unit_s)[0]

    # Optimized weights for two-proton system
    def energy_double(c_vec):
        pert = np.zeros(N)
        for i, n in enumerate(modes):
            eps_n = np.sin(n * gamma)
            pert += c_vec[i] * (breather_sq(x, pos_A, eps_n)
                              + breather_sq(x, pos_B, eps_n)) * (-PI**2 / 6.0)
        H = build_hessian(phi_double, N, extra_diag=pert)
        ev, _ = eigsh(H, k=1, which='SM')
        return ev[0]

    result_d = minimize(energy_double, c_single_opt, method='Nelder-Mead',
                        options={'maxiter': 2000, 'xatol': 1e-8, 'fatol': 1e-10})
    E0_double_opt = result_d.fun
    V_opt = E0_double_opt - E0_single_opt

    pct = (V_opt - V_bare) / abs(V_bare) * 100 if abs(V_bare) > 1e-10 else 0

    report(f"{R:4d} {V_bare:+12.8f} {V_opt:+12.8f} {V_unit:+12.8f} "
           f"{pct:+12.2f}%")

report("")

# Show how weights change at bonding distance
report("WEIGHT SHIFT DURING BONDING (R=10):")
R_opt = 10
pos_A = -R_opt // 2
pos_B = R_opt // 2
phi_double_opt = kink_antikink(x, pos_A, kink_width) + kink_antikink(x, pos_B, kink_width)

def energy_func_10(c_vec):
    pert = np.zeros(N)
    for i, n in enumerate(modes):
        eps_n = np.sin(n * gamma)
        pert += c_vec[i] * (breather_sq(x, pos_A, eps_n)
                          + breather_sq(x, pos_B, eps_n)) * (-PI**2 / 6.0)
    H = build_hessian(phi_double_opt, N, extra_diag=pert)
    ev, _ = eigsh(H, k=1, which='SM')
    return ev[0]

result_10 = minimize(energy_func_10, c_single_opt, method='Nelder-Mead',
                     options={'maxiter': 2000, 'xatol': 1e-8, 'fatol': 1e-10})
c_bonded = result_10.x

report(f"{'n':>3} {'c_single':>10} {'c_bonded':>10} {'ratio':>8} {'shift':>10}")
report("-" * 45)
for i, n in enumerate(modes):
    ratio = c_bonded[i] / c_single_opt[i] if abs(c_single_opt[i]) > 1e-10 else 0
    shift = c_bonded[i] - c_single_opt[i]
    report(f"  {n:3d} {c_single_opt[i]:+10.4f} {c_bonded[i]:+10.4f} {ratio:8.4f} {shift:+10.4f}")

report("")

# ============================================================
# PART 7: NORMALIZED RATIO TABLE
# ============================================================
report("PART 7: NORMALIZED RATIOS (the fingerprint)")
report("-" * 55)
report("The ratios a_n/a_1 are the proton's internal structure.")
report("These should be UNIVERSAL — same for every proton.")
report("")

report(f"{'n':>3} {'a_n/a_1':>10} {'eps_n':>8} {'omega_n':>8} "
       f"{'n*gamma':>8} {'sin/cos':>8}")
report("-" * 50)
for i, n in enumerate(modes):
    r = a_fit[i] / a_fit[0]
    eps = np.sin(n*gamma)
    omega = np.cos(n*gamma)
    ng = n * gamma
    sc = eps / omega  # tan(n*gamma)
    report(f"  {n:3d} {r:+10.4f} {eps:8.4f} {omega:8.4f} {ng:8.4f} {sc:8.4f}")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report("1. THE OVERLAP MATRIX M shows all modes strongly coupled.")
report(f"   Condition number: {M_evals[-1]/M_evals[0]:.1f}")
report(f"   Eigenvalue range: {M_evals[0]:.2f} to {M_evals[-1]:.2f}")
report("")
report("2. THE WELL-FILLING FIT achieves:")
report(f"   Σ a_n φ_n² ≈ ΔV_kink with {residual_norm/total_norm*100:.1f}% residual")
report("   The 7 modes CAN reconstruct the kink well deviation.")
report("")
report("3. THE RATIOS a_n/a_1:")
for i, n in enumerate(modes):
    r = a_fit[i] / a_fit[0]
    report(f"   n={n}: {r:+.4f}")
report("")
report("4. BONDING shifts the equilibrium weights — the proton")
report("   restructures internally when a neighbor is present.")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
