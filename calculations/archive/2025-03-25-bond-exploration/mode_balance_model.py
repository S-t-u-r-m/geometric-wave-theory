"""
Mode Balance Model — Constrained 1-2 Parameter Bond Model
============================================================
The 7 breather modes nearly cancel (99.7% cancellation with equal weights).
The well-filling vector W shows:
  Modes 1-4: SOFTEN the kink well (W > 0)
  Mode 5:    NEUTRAL (W ≈ 0)
  Modes 6-7: STIFFEN the kink well (W < 0)

The sum rule Σφ²/|ΔV| = 7 fixes the total perturbation strength.
The BALANCE between softening and stiffening determines the eigenvalue.

Constrained model:
  c_n = 1 + α × W_n / |W|   (one parameter: tilt along W direction)

When bonding, the neighbor shifts the effective W_n for each mode
(wider modes feel it more), changing the optimal α.

This reduces 7 free parameters to 1 (or 2 with a second direction).
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "mode_balance_model_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("MODE BALANCE MODEL — CONSTRAINED BOND MODEL")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"d = {d}, gamma = {gamma:.10f}")
report("")

# ============================================================
# SETUP
# ============================================================
N = 512
x = np.arange(N, dtype=np.float64) - N/2
kink_width = 3
modes = list(range(1, 8))
n_modes = len(modes)

def kink_antikink(x, center, width):
    return (4.0/PI) * (np.arctan(np.exp(x - center + width/2))
                      - np.arctan(np.exp(x - center - width/2)))

def breather_sq(x, center, eps_n):
    phi = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * (x - center) + 1e-30))
    return phi**2

def build_hessian(phi_bg, extra_diag=None):
    diag = 2.0 + np.cos(PI * phi_bg)
    if extra_diag is not None:
        diag = diag + extra_diag
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()

phi_kink = kink_antikink(x, 0, kink_width)
V_kink = np.cos(PI * phi_kink)
delta_V = V_kink - 1.0  # deviation from vacuum

# ============================================================
# BUILD THE BASIS
# ============================================================
# Profiles centered at origin (single proton)
profiles_0 = {}
for n in modes:
    eps_n = np.sin(n * gamma)
    profiles_0[n] = breather_sq(x, 0, eps_n)

# Well-filling vector W
W = np.zeros(n_modes)
for i, n in enumerate(modes):
    W[i] = np.sum(profiles_0[n] * delta_V)

W_norm = np.sqrt(np.sum(W**2))
W_hat = W / W_norm  # unit direction

report("WELL-FILLING DIRECTION W:")
report(f"{'n':>3} {'W_n':>10} {'W_hat_n':>10} {'group':>10}")
report("-" * 38)
for i, n in enumerate(modes):
    group = "SOFTEN" if W[i] > 0.1 else "STIFFEN" if W[i] < -0.1 else "NEUTRAL"
    report(f"  {n:3d} {W[i]:10.4f} {W_hat[i]:10.4f} {group:>10}")
report(f"  |W| = {W_norm:.4f}")
report("")

# Overlap matrix — get the second direction from its eigenvectors
M = np.zeros((n_modes, n_modes))
for i, n in enumerate(modes):
    for j, m in enumerate(modes):
        M[i, j] = np.sum(profiles_0[n] * profiles_0[m])

M_evals, M_evecs = np.linalg.eigh(M)

# The dominant eigenvector of M (largest eigenvalue)
v_dom = M_evecs[:, -1]  # all same sign — this is the "total size" direction
# Second direction: orthogonal to W_hat, next most important
v2 = M_evecs[:, -2]
# Orthogonalize v2 against W_hat
v2 = v2 - np.dot(v2, W_hat) * W_hat
v2 = v2 / (np.linalg.norm(v2) + 1e-30)

report("SECOND DIRECTION (orthogonal to W, from M eigenvector):")
for i, n in enumerate(modes):
    report(f"  n={n}: v2 = {v2[i]:+.4f}")
report("")

# ============================================================
# PART 1: 1-PARAMETER MODEL — SCAN α
# ============================================================
report("PART 1: 1-PARAMETER MODEL (tilt along W)")
report("-" * 55)
report("c_n = 1 + α × W_hat_n")
report("Scan α to find the optimal single-proton eigenvalue.")
report("")

# The perturbation for weight vector c is:
# delta_diag(x) = -(π²/6) × Σ_n c_n × φ_n²(x)

def compute_E0(c_vec, phi_bg, profiles_dict):
    """Compute lowest eigenvalue with given mode weights."""
    pert = np.zeros(N)
    for i, n in enumerate(modes):
        pert += c_vec[i] * profiles_dict[n]
    pert *= (-PI**2 / 6.0)
    try:
        H = build_hessian(phi_bg, extra_diag=pert)
        ev, _ = eigsh(H, k=1, which='SM', maxiter=5000)
        return ev[0]
    except Exception:
        return 999.0

# Scan alpha for single proton
alpha_values = np.linspace(-0.5, 0.5, 51)
E0_vs_alpha = []

report(f"{'alpha':>8} {'E0':>12} {'c_1':>8} {'c_4':>8} {'c_7':>8}")
report("-" * 45)

for alpha in alpha_values:
    c = 1.0 + alpha * W_hat
    E0 = compute_E0(c, phi_kink, profiles_0)
    E0_vs_alpha.append(E0)
    if abs(alpha) < 0.01 or abs(alpha - 0.1) < 0.01 or abs(alpha + 0.1) < 0.01 or \
       abs(alpha - 0.3) < 0.01 or abs(alpha + 0.3) < 0.01 or abs(alpha - 0.5) < 0.01:
        report(f"  {alpha:+8.3f} {E0:12.6f} {c[0]:8.4f} {c[3]:8.4f} {c[6]:8.4f}")

E0_vs_alpha = np.array(E0_vs_alpha)
i_min = np.argmin(E0_vs_alpha)
alpha_opt = alpha_values[i_min]
E0_opt = E0_vs_alpha[i_min]

report("")
report(f"OPTIMAL: α = {alpha_opt:.3f}, E0 = {E0_opt:.8f}")

c_opt_single = 1.0 + alpha_opt * W_hat
report(f"Optimal weights c_n:")
for i, n in enumerate(modes):
    report(f"  n={n}: c = {c_opt_single[i]:.4f}")
report("")

# ============================================================
# PART 2: BONDING — HOW α SHIFTS WITH R
# ============================================================
report("PART 2: BONDING — OPTIMAL α VS SEPARATION R")
report("-" * 55)
report("At each R, scan α to find the optimal two-proton eigenvalue.")
report("The bond energy = E0(R, α_opt(R)) - E0(∞, α_opt(∞))")
report("")

R_values = [6, 7, 8, 9, 10, 12, 14, 16, 20, 24, 30]

# Reference: single proton with optimal α
E0_ref = E0_opt

# For the bare (no breather) comparison
H_bare_single = build_hessian(phi_kink)
ev_bare_s, _ = eigsh(H_bare_single, k=1, which='SM')
E0_bare_single = ev_bare_s[0]

report(f"{'R':>4} {'alpha_opt':>10} {'E0_double':>12} {'V(R)':>12} "
       f"{'V_bare':>12} {'ratio':>8}")
report("-" * 65)

bond_results = {}

for R in R_values:
    pos_A = -R / 2.0
    pos_B = R / 2.0

    phi_double = kink_antikink(x, pos_A, kink_width) + kink_antikink(x, pos_B, kink_width)

    # Build profiles for both wells
    profiles_AB = {}
    for n in modes:
        eps_n = np.sin(n * gamma)
        profiles_AB[n] = breather_sq(x, pos_A, eps_n) + breather_sq(x, pos_B, eps_n)

    # Bare (no perturbation)
    H_bare = build_hessian(phi_double)
    ev_bare, _ = eigsh(H_bare, k=1, which='SM')
    V_bare = ev_bare[0] - E0_bare_single

    # Scan alpha
    best_E0 = 999
    best_alpha = 0
    for alpha in np.linspace(-0.5, 0.5, 101):
        c = 1.0 + alpha * W_hat
        E0 = compute_E0(c, phi_double, profiles_AB)
        if E0 < best_E0:
            best_E0 = E0
            best_alpha = alpha

    # Fine-tune around best
    for alpha in np.linspace(best_alpha - 0.02, best_alpha + 0.02, 41):
        c = 1.0 + alpha * W_hat
        E0 = compute_E0(c, phi_double, profiles_AB)
        if E0 < best_E0:
            best_E0 = E0
            best_alpha = alpha

    V_opt = best_E0 - 2 * E0_ref  # two isolated protons
    ratio = V_opt / V_bare if abs(V_bare) > 1e-10 else 0

    bond_results[R] = {
        'alpha': best_alpha, 'E0': best_E0,
        'V': V_opt, 'V_bare': V_bare
    }

    report(f"{R:4.0f} {best_alpha:+10.4f} {best_E0:12.6f} {V_opt:+12.8f} "
           f"{V_bare:+12.8f} {ratio:8.3f}")

report("")

# ============================================================
# PART 3: WEIGHT SHIFT DURING BONDING
# ============================================================
report("PART 3: HOW EACH MODE'S WEIGHT CHANGES DURING BONDING")
report("-" * 55)
report(f"Single proton: α = {alpha_opt:.4f}")
report("")

report(f"{'R':>4} {'Δα':>8} {'Δc_1':>8} {'Δc_2':>8} {'Δc_3':>8} "
       f"{'Δc_4':>8} {'Δc_5':>8} {'Δc_6':>8} {'Δc_7':>8}")
report("-" * 75)

for R in R_values:
    alpha_R = bond_results[R]['alpha']
    d_alpha = alpha_R - alpha_opt
    c_R = 1.0 + alpha_R * W_hat
    dc = c_R - c_opt_single
    line = f"{R:4.0f} {d_alpha:+8.4f}"
    for i in range(n_modes):
        line += f" {dc[i]:+8.4f}"
    report(line)

report("")

# ============================================================
# PART 4: 2-PARAMETER MODEL
# ============================================================
report("PART 4: 2-PARAMETER MODEL (α along W, β along v2)")
report("-" * 55)
report("c_n = 1 + α × W_hat_n + β × v2_n")
report("")

# Scan (α, β) grid for single proton
alpha_range = np.linspace(-0.3, 0.3, 31)
beta_range = np.linspace(-0.3, 0.3, 31)

best_E0_2d = 999
best_ab = (0, 0)

for alpha in alpha_range:
    for beta in beta_range:
        c = 1.0 + alpha * W_hat + beta * v2
        E0 = compute_E0(c, phi_kink, profiles_0)
        if E0 < best_E0_2d:
            best_E0_2d = E0
            best_ab = (alpha, beta)

report(f"Single proton: α={best_ab[0]:.3f}, β={best_ab[1]:.3f}, "
       f"E0={best_E0_2d:.8f}")
report(f"  (1-param gave E0={E0_opt:.8f}, improvement: "
       f"{(best_E0_2d - E0_opt):.8f})")
report("")

# 2-param bond curve
report("2-PARAMETER BOND CURVE:")
report(f"{'R':>4} {'α_opt':>8} {'β_opt':>8} {'V_2param':>12} "
       f"{'V_1param':>12} {'V_bare':>12}")
report("-" * 60)

E0_ref_2d = best_E0_2d

for R in R_values:
    pos_A = -R / 2.0
    pos_B = R / 2.0
    phi_double = kink_antikink(x, pos_A, kink_width) + kink_antikink(x, pos_B, kink_width)

    profiles_AB = {}
    for n in modes:
        eps_n = np.sin(n * gamma)
        profiles_AB[n] = breather_sq(x, pos_A, eps_n) + breather_sq(x, pos_B, eps_n)

    best_E0_R = 999
    best_ab_R = (0, 0)
    for alpha in np.linspace(-0.3, 0.3, 31):
        for beta in np.linspace(-0.3, 0.3, 31):
            c = 1.0 + alpha * W_hat + beta * v2
            E0 = compute_E0(c, phi_double, profiles_AB)
            if E0 < best_E0_R:
                best_E0_R = E0
                best_ab_R = (alpha, beta)

    V_2p = best_E0_R - 2 * E0_ref_2d
    V_1p = bond_results[R]['V']
    V_bare = bond_results[R]['V_bare']

    report(f"{R:4.0f} {best_ab_R[0]:+8.3f} {best_ab_R[1]:+8.3f} "
           f"{V_2p:+12.8f} {V_1p:+12.8f} {V_bare:+12.8f}")

report("")

# ============================================================
# PART 5: PHYSICAL INTERPRETATION — THE SOFTENING/STIFFENING BALANCE
# ============================================================
report("PART 5: PHYSICAL INTERPRETATION")
report("-" * 55)
report("")
report("The W vector divides modes into:")
report("  SOFTENING (W>0): modes 1-4  (wide, low n)")
report("  NEUTRAL (W≈0):   mode 5")
report("  STIFFENING (W<0): modes 6-7  (narrow, high n)")
report("")
report("Equilibrium: softening and stiffening nearly cancel.")
report("Bonding shifts α → more/less softening depending on R.")
report("")

# Compute: at each R, what fraction of the change is in
# softening vs stiffening modes?
report("SOFTENING vs STIFFENING BALANCE:")
report(f"{'R':>4} {'Δα':>8} {'Δ(soft)':>10} {'Δ(stiff)':>10} {'ratio':>8}")
report("-" * 45)

for R in R_values:
    alpha_R = bond_results[R]['alpha']
    d_alpha = alpha_R - alpha_opt
    dc = d_alpha * W_hat

    # Softening modes (1-4)
    d_soft = sum(dc[i] for i in range(4))
    # Stiffening modes (6-7)
    d_stiff = sum(dc[i] for i in range(5, 7))

    ratio = d_soft / d_stiff if abs(d_stiff) > 1e-10 else 0
    report(f"{R:4.0f} {d_alpha:+8.4f} {d_soft:+10.6f} {d_stiff:+10.6f} {ratio:+8.3f}")

report("")

# ============================================================
# PART 6: WHAT DETERMINES THE RATIO Σφ²/|ΔV| = 7?
# ============================================================
report("PART 6: THE SUM RULE")
report("-" * 55)
report("")

# Check: is the ratio exactly 7 (number of modes)?
sum_phi_sq = sum(np.sum(profiles_0[n]) for n in modes)
sum_abs_dV = np.sum(np.abs(delta_V))
ratio_7 = sum_phi_sq / sum_abs_dV
report(f"Σ integral(φ_n²) = {sum_phi_sq:.4f}")
report(f"integral(|ΔV|) = {sum_abs_dV:.4f}")
report(f"Ratio = {ratio_7:.4f}")
report("")

# Each mode's integral of phi_n^2
report("Individual integrals:")
report(f"{'n':>3} {'∫φ_n²':>10} {'∫φ_n²/∫|ΔV|':>13} {'eps_n':>8} {'1/eps_n':>8}")
report("-" * 45)
for n in modes:
    integ = np.sum(profiles_0[n])
    eps_n = np.sin(n * gamma)
    report(f"  {n:3d} {integ:10.4f} {integ/sum_abs_dV:13.4f} {eps_n:8.4f} {1/eps_n:8.2f}")

report("")
report("Note: ∫φ_n² ∝ 1/eps_n (breather width).")
report("The sum rule emerges because the TOTAL width of all 7 modes")
report("spans the kink well exactly.")
report("")

# Is sum of 1/eps_n related to anything?
sum_inv_eps = sum(1.0/np.sin(n*gamma) for n in modes)
report(f"Σ 1/eps_n = {sum_inv_eps:.4f}")
report(f"This should relate to the total available phase space.")

# Sum of widths in units of kink_width
report(f"Σ width / kink_width = {sum_inv_eps / kink_width:.4f}")
report(f"N_modes × avg_width / kink_width = {sum_inv_eps / (n_modes * kink_width) * n_modes:.4f}")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report("1. The 7 breather modes form a CONSTRAINED EQUILIBRIUM.")
report("   Modes 1-4 soften, modes 6-7 stiffen, the balance nearly cancels.")
report("")
report(f"2. Single parameter α along the W direction captures the physics.")
report(f"   Optimal single proton: α = {alpha_opt:.4f}")
report("")
report("3. During bonding, α shifts because the neighbor modifies")
report("   the effective well. Wider modes (1-2) feel it first.")
report("")
report("4. The bond energy comes from the SHIFT in equilibrium,")
report("   not from individual mode contributions.")

if bond_results:
    V_8 = bond_results.get(8, {}).get('V', 0)
    V_bare_8 = bond_results.get(8, {}).get('V_bare', 0)
    if V_bare_8 != 0:
        report(f"\n5. At R=8: V_constrained/V_bare = {V_8/V_bare_8:.3f}")
        report(f"   Bond energy is modified by the mode balance correction.")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
