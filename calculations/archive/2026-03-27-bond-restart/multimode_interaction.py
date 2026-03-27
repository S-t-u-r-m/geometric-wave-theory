"""
Multi-Mode Breather Interaction Matrix
========================================
A real proton has 7 stable breather modes (n=1..7) coexisting in the same
kink well. When two protons bond, ALL modes interact.

Three questions:
  1. INTRA-PROTON: Do the 7 breathers shift each other's frequencies?
     (Are they independent or coupled within a single proton?)

  2. INTER-PROTON, SAME MODE: When two protons approach, how does each
     mode's splitting depend on its width (1/eps_n)?

  3. INTER-PROTON, CROSS MODE: Does breather n on proton A couple to
     breather m on proton B? Is this additive or non-linear?

Method:
  - 1D lattice with two kink-antikink pairs at separation R
  - Populate each well with breather perturbations (modes 1-7)
  - The breather at mode n modifies the Hessian diagonal by:
      delta_V_n(x) = -(pi * phi_n(x))^2 / 6  (leading phi^4 coupling)
  - Compare eigenspectra:
      (a) Empty wells (bare kink-antikink, no breathers)
      (b) Single breather mode in each well
      (c) All 7 modes in each well simultaneously
  - If (c) = sum of (b)'s, modes are independent. If not, they're coupled.

All parameters from GWT Lagrangian (d=3, zero free parameters).
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "multimode_interaction_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("MULTI-MODE BREATHER INTERACTION MATRIX")
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

def kink_profile(x, x_center):
    return (4.0/PI) * np.arctan(np.exp(x - x_center))

def antikink_profile(x, x_center):
    return 2.0 - (4.0/PI) * np.arctan(np.exp(x - x_center))

def make_kink_antikink(x, center, width):
    return kink_profile(x, center - width/2) + antikink_profile(x, center + width/2) - 2.0

def breather_profile(x, center, eps_n):
    return (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * (x - center) + 1e-30))

def build_hessian(phi_background, N, extra_diag=None):
    """Hessian with optional extra diagonal perturbation."""
    diag = 2.0 + np.cos(PI * phi_background)
    if extra_diag is not None:
        diag = diag + extra_diag
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()

def breather_perturbation(x, center, eps_n):
    """Leading-order phi^4 coupling: delta_V = -(pi*phi_breather)^2 / 6."""
    phi = breather_profile(x, center, eps_n)
    return -(PI * phi)**2 / 6.0

# ============================================================
# BREATHER MODE TABLE
# ============================================================
report("BREATHER MODE TABLE")
report("-" * 55)
report(f"{'n':>3} {'omega_n':>10} {'eps_n':>10} {'width':>8} {'phi_max':>8}")
report("-" * 45)

modes = list(range(1, 8))  # n = 1..7
for n in modes:
    omega = np.cos(n * gamma)
    eps = np.sin(n * gamma)
    width = 1.0 / eps
    phi_max = (4.0/PI) * np.arctan(1.0)  # at center, same for all modes
    report(f"  {n:3d} {omega:10.6f} {eps:10.6f} {width:8.2f} {phi_max:8.4f}")

report("")

# ============================================================
# PART 1: SINGLE PROTON — INTRA-PROTON COUPLING
# ============================================================
report("PART 1: INTRA-PROTON MODE COUPLING")
report("-" * 55)
report("Does mode n shift when modes 1..n-1 are also present?")
report("")

# Base: single kink-antikink, no breathers
phi_base = make_kink_antikink(x, center, kink_width)
H_bare = build_hessian(phi_base, N)
ev_bare, _ = eigsh(H_bare, k=10, which='SM')
ev_bare = np.sort(ev_bare)
E0_bare = ev_bare[0]

report(f"Bare kink well: E0 = {E0_bare:.8f}")
report("")

# Add breathers one at a time, track how E0 shifts
report("CUMULATIVE MODE ADDITION:")
report(f"{'modes_present':>20} {'E0':>12} {'shift_from_bare':>15} "
       f"{'shift_from_prev':>15} {'notes':>15}")
report("-" * 80)

report(f"{'none':>20} {E0_bare:12.8f} {0:+15.8f} {0:+15.8f} {'bare':>15}")

# Track eigenvalue with each mode added
cumulative_perturbation = np.zeros(N)
prev_E0 = E0_bare
individual_shifts = {}

for n in modes:
    eps_n = np.sin(n * gamma)

    # Individual mode shift (mode n alone)
    pert_n = breather_perturbation(x, center, eps_n)
    H_n = build_hessian(phi_base, N, extra_diag=pert_n)
    ev_n, _ = eigsh(H_n, k=5, which='SM')
    ev_n = np.sort(ev_n)
    individual_shifts[n] = ev_n[0] - E0_bare

    # Cumulative: add mode n to all previous
    cumulative_perturbation += pert_n
    H_cum = build_hessian(phi_base, N, extra_diag=cumulative_perturbation.copy())
    ev_cum, _ = eigsh(H_cum, k=5, which='SM')
    ev_cum = np.sort(ev_cum)

    shift_bare = ev_cum[0] - E0_bare
    shift_prev = ev_cum[0] - prev_E0
    prev_E0 = ev_cum[0]

    modes_str = f"1..{n}"
    report(f"{modes_str:>20} {ev_cum[0]:12.8f} {shift_bare:+15.8f} "
           f"{shift_prev:+15.8f}")

report("")

# Compare: sum of individual shifts vs actual cumulative shift
sum_individual = sum(individual_shifts[n] for n in modes)
H_all = build_hessian(phi_base, N, extra_diag=cumulative_perturbation)
ev_all, _ = eigsh(H_all, k=5, which='SM')
ev_all = np.sort(ev_all)
actual_shift = ev_all[0] - E0_bare

report("ADDITIVITY TEST:")
report(f"  Sum of individual mode shifts:  {sum_individual:+.8f}")
report(f"  Actual cumulative shift (all 7): {actual_shift:+.8f}")
report(f"  Difference (non-additive part): {actual_shift - sum_individual:+.8f}")
report(f"  Non-additive fraction: {(actual_shift - sum_individual)/actual_shift * 100:.2f}%")
report("")

# Individual shifts table
report("INDIVIDUAL MODE SHIFTS (each mode alone):")
report(f"{'n':>3} {'eps_n':>10} {'width':>8} {'shift':>15} {'shift/E0':>12}")
report("-" * 55)
for n in modes:
    eps = np.sin(n * gamma)
    width = 1.0/eps
    shift = individual_shifts[n]
    report(f"  {n:3d} {eps:10.6f} {width:8.2f} {shift:+15.8f} {shift/E0_bare:+12.6f}")

report("")

# ============================================================
# PART 2: INTER-PROTON — EACH MODE'S BONDING RANGE
# ============================================================
report("PART 2: INTER-PROTON — MODE-BY-MODE BONDING RANGE")
report("-" * 55)
report("Two protons at separation R. How does each mode's splitting")
report("depend on R? Wide modes (small n) should bond at larger R.")
report("")

R_values = [6, 8, 10, 12, 14, 16, 20, 24, 30]

# For each mode, add that mode's breather perturbation to BOTH wells
# and track the bonding-antibonding splitting

header = f"{'R':>4} " + "".join(f"  {'split_n'+str(n):>11}" for n in modes)
report(header)
report("-" * (4 + 7 * 13))

split_data = {}  # {n: [splits at each R]}
for n in modes:
    split_data[n] = []

for R in R_values:
    pos_A = center - R // 2
    pos_B = center + R // 2

    # Two bare kink-antikink pairs
    phi_double = make_kink_antikink(x, pos_A, kink_width) + make_kink_antikink(x, pos_B, kink_width)

    line = f"{R:4d} "
    for n in modes:
        eps_n = np.sin(n * gamma)

        # Add mode-n breather to BOTH wells
        pert = breather_perturbation(x, pos_A, eps_n) + breather_perturbation(x, pos_B, eps_n)

        H = build_hessian(phi_double, N, extra_diag=pert)
        ev, _ = eigsh(H, k=4, which='SM')
        ev = np.sort(ev)

        split = ev[1] - ev[0] if len(ev) >= 2 else 0
        split_data[n].append(split)
        line += f"  {split:11.8f}"

    report(line)

report("")

# ============================================================
# PART 3: SPLITTING DECAY RATES — DOES WIDTH MATTER?
# ============================================================
report("PART 3: SPLITTING DECAY RATES BY MODE")
report("-" * 55)
report("Fit each mode's splitting vs R to extract decay rate.")
report("Prediction: wider mode (smaller n, larger 1/eps_n) decays slower.")
report("")

R_arr = np.array(R_values, dtype=float)

report(f"{'n':>3} {'eps_n':>8} {'width':>6} {'decay_rate':>11} "
       f"{'ratio/eps':>10} {'split_R8':>11} {'split_R20':>11}")
report("-" * 65)

for n in modes:
    eps_n = np.sin(n * gamma)
    splits = np.array(split_data[n])

    # Fit log(split) vs R where split > 0
    mask = splits > 1e-12
    if np.sum(mask) > 2:
        log_s = np.log(splits[mask])
        R_fit = R_arr[mask]
        try:
            coeffs = np.polyfit(R_fit, log_s, 1)
            decay = -coeffs[0]
        except Exception:
            decay = 0
    else:
        decay = 0

    split_8 = split_data[n][R_values.index(8)] if 8 in R_values else 0
    split_20 = split_data[n][R_values.index(20)] if 20 in R_values else 0

    report(f"  {n:3d} {eps_n:8.4f} {1/eps_n:6.1f} {decay:11.6f} "
           f"{decay/eps_n:10.2f} {split_8:11.8f} {split_20:11.8f}")

report("")

# ============================================================
# PART 4: FULL MULTI-MODE BONDING — ALL 7 IN BOTH WELLS
# ============================================================
report("PART 4: FULL MULTI-MODE BOND CURVE")
report("-" * 55)
report("All 7 breather modes in BOTH kink wells simultaneously.")
report("Compare to: bare wells, and sum-of-individual-modes prediction.")
report("")

report(f"{'R':>4} {'V_bare':>12} {'V_all7':>12} {'V_sum_indiv':>12} "
       f"{'nonadditive':>12} {'nonadd_%':>10}")
report("-" * 68)

# First get single-well reference with all 7 modes
pert_all_single = np.zeros(N)
for n in modes:
    eps_n = np.sin(n * gamma)
    pert_all_single += breather_perturbation(x, center, eps_n)

H_ref = build_hessian(phi_base, N, extra_diag=pert_all_single)
ev_ref, _ = eigsh(H_ref, k=5, which='SM')
E0_ref_all = np.sort(ev_ref)[0]

for R in R_values:
    pos_A = center - R // 2
    pos_B = center + R // 2
    phi_double = make_kink_antikink(x, pos_A, kink_width) + make_kink_antikink(x, pos_B, kink_width)

    # (a) Bare: no breather perturbation
    H_bare2 = build_hessian(phi_double, N)
    ev_bare2, _ = eigsh(H_bare2, k=4, which='SM')
    ev_bare2 = np.sort(ev_bare2)
    V_bare = ev_bare2[0] - E0_bare

    # (b) All 7 modes in both wells
    pert_all = np.zeros(N)
    for n in modes:
        eps_n = np.sin(n * gamma)
        pert_all += breather_perturbation(x, pos_A, eps_n)
        pert_all += breather_perturbation(x, pos_B, eps_n)

    H_all7 = build_hessian(phi_double, N, extra_diag=pert_all)
    ev_all7, _ = eigsh(H_all7, k=4, which='SM')
    ev_all7 = np.sort(ev_all7)
    V_all7 = ev_all7[0] - E0_ref_all

    # (c) Sum of individual mode shifts
    V_sum = 0
    for n in modes:
        eps_n = np.sin(n * gamma)
        pert_n = breather_perturbation(x, pos_A, eps_n) + breather_perturbation(x, pos_B, eps_n)
        H_n = build_hessian(phi_double, N, extra_diag=pert_n)
        ev_n, _ = eigsh(H_n, k=4, which='SM')
        ev_n = np.sort(ev_n)
        V_sum += (ev_n[0] - E0_bare) - individual_shifts[n]  # subtract self-energy

    nonadd = V_all7 - V_sum if abs(V_sum) > 1e-12 else 0
    nonadd_pct = nonadd / V_all7 * 100 if abs(V_all7) > 1e-12 else 0

    report(f"{R:4d} {V_bare:+12.8f} {V_all7:+12.8f} {V_sum:+12.8f} "
           f"{nonadd:+12.8f} {nonadd_pct:+10.2f}%")

report("")

# ============================================================
# PART 5: CROSS-MODE COUPLING MATRIX
# ============================================================
report("PART 5: CROSS-MODE COUPLING MATRIX AT R=10")
report("-" * 55)
report("How does adding mode m to well B change mode n's splitting?")
report("Rows = mode in well A, Columns = mode added to well B.")
report("Entry = change in splitting when cross-mode is present.")
report("")

R_test = 10
pos_A = center - R_test // 2
pos_B = center + R_test // 2
phi_double = make_kink_antikink(x, pos_A, kink_width) + make_kink_antikink(x, pos_B, kink_width)

# Baseline: mode n in well A only, nothing extra in well B
baseline_splits = {}
for n in modes:
    eps_n = np.sin(n * gamma)
    pert = breather_perturbation(x, pos_A, eps_n) + breather_perturbation(x, pos_B, eps_n)
    H = build_hessian(phi_double, N, extra_diag=pert)
    ev, _ = eigsh(H, k=4, which='SM')
    ev = np.sort(ev)
    baseline_splits[n] = ev[1] - ev[0]

# Cross-mode: mode n in both wells + mode m in both wells
header = f"{'n\\m':>5}"
for m in modes:
    header += f"  {m:>9}"
report(header)
report("-" * (5 + 7 * 11))

for n in modes:
    eps_n = np.sin(n * gamma)
    pert_n = breather_perturbation(x, pos_A, eps_n) + breather_perturbation(x, pos_B, eps_n)
    line = f"  {n:3d}"

    for m in modes:
        eps_m = np.sin(m * gamma)
        pert_m = breather_perturbation(x, pos_A, eps_m) + breather_perturbation(x, pos_B, eps_m)

        if m == n:
            # Same mode — just report baseline
            line += f"  {baseline_splits[n]:9.2e}"
        else:
            # Both modes present
            pert_nm = pert_n + pert_m
            H = build_hessian(phi_double, N, extra_diag=pert_nm)
            ev, _ = eigsh(H, k=4, which='SM')
            ev = np.sort(ev)
            split_nm = ev[1] - ev[0]

            # Change from baseline (mode n alone)
            delta = split_nm - baseline_splits[n]
            line += f"  {delta:+9.2e}"

    report(line)

report("")
report("Diagonal = baseline split for that mode (both wells)")
report("Off-diagonal = CHANGE in split when cross-mode is added")
report("If all off-diagonal ≈ 0, modes are independent.")
report("")

# ============================================================
# PART 6: OCCUPATION-WEIGHTED BOND POTENTIAL
# ============================================================
report("PART 6: OCCUPATION-WEIGHTED BOND POTENTIAL")
report("-" * 55)
report("In a real atom, not all 7 modes are equally occupied.")
report("Weight each mode by its occupation (from electron config).")
report("Example: hydrogen = 1 electron in mode n=1.")
report("         Carbon = modes 1-4 occupied (2+2+2 electrons).")
report("")

# Compute bond curve for different occupation patterns
configs = {
    'H (n=1 only)': [1],
    'He (n=1 paired)': [1],
    'C (n=1,2,3,4)': [1, 2, 3, 4],
    'N (n=1..5)': [1, 2, 3, 4, 5],
    'O (n=1..6)': [1, 2, 3, 4, 5, 6],
    'F (n=1..7)': [1, 2, 3, 4, 5, 6, 7],
}

header = f"{'R':>4}" + "".join(f"  {name:>14}" for name in configs)
report(header)
report("-" * (4 + len(configs) * 16))

for R in R_values:
    pos_A = center - R // 2
    pos_B = center + R // 2
    phi_double = make_kink_antikink(x, pos_A, kink_width) + make_kink_antikink(x, pos_B, kink_width)

    line = f"{R:4d}"
    for name, occ_modes in configs.items():
        # Build perturbation for occupied modes only
        pert = np.zeros(N)
        pert_single = np.zeros(N)  # single-well reference
        for n in occ_modes:
            eps_n = np.sin(n * gamma)
            pert += breather_perturbation(x, pos_A, eps_n) + breather_perturbation(x, pos_B, eps_n)
            pert_single += breather_perturbation(x, center, eps_n)

        # Two-well system
        H = build_hessian(phi_double, N, extra_diag=pert)
        ev, _ = eigsh(H, k=4, which='SM')
        ev = np.sort(ev)

        # Single-well reference with same modes
        H_s = build_hessian(phi_base, N, extra_diag=pert_single)
        ev_s, _ = eigsh(H_s, k=4, which='SM')
        ev_s = np.sort(ev_s)

        V = ev[0] - ev_s[0]
        line += f"  {V:+14.8f}"

    report(line)

report("")

# ============================================================
# PART 7: MODE COMPETITION — ENERGY PARTITION
# ============================================================
report("PART 7: MODE COMPETITION AT BONDING DISTANCE")
report("-" * 55)
report("At the bonding distance, which modes contribute most to D_e?")
report("This tells us the effective bond order decomposition.")
report("")

R_bond = 8  # typical bonding distance
pos_A = center - R_bond // 2
pos_B = center + R_bond // 2
phi_double = make_kink_antikink(x, pos_A, kink_width) + make_kink_antikink(x, pos_B, kink_width)

# Full 7-mode bond energy
pert_full = np.zeros(N)
pert_full_single = np.zeros(N)
for n in modes:
    eps_n = np.sin(n * gamma)
    pert_full += breather_perturbation(x, pos_A, eps_n) + breather_perturbation(x, pos_B, eps_n)
    pert_full_single += breather_perturbation(x, center, eps_n)

H_full = build_hessian(phi_double, N, extra_diag=pert_full)
ev_full, _ = eigsh(H_full, k=4, which='SM')
ev_full = np.sort(ev_full)

H_full_s = build_hessian(phi_base, N, extra_diag=pert_full_single)
ev_full_s, _ = eigsh(H_full_s, k=4, which='SM')
ev_full_s = np.sort(ev_full_s)

V_full = ev_full[0] - ev_full_s[0]
report(f"Full 7-mode bond energy at R={R_bond}: V = {V_full:+.8f}")
report("")

# Remove one mode at a time — see how much D_e changes
report(f"{'removed':>10} {'V_without':>14} {'delta_V':>14} {'contrib_%':>10}")
report("-" * 55)

for n_remove in modes:
    pert_minus = np.zeros(N)
    pert_minus_s = np.zeros(N)
    for n in modes:
        if n == n_remove:
            continue
        eps_n = np.sin(n * gamma)
        pert_minus += breather_perturbation(x, pos_A, eps_n) + breather_perturbation(x, pos_B, eps_n)
        pert_minus_s += breather_perturbation(x, center, eps_n)

    H_m = build_hessian(phi_double, N, extra_diag=pert_minus)
    ev_m, _ = eigsh(H_m, k=4, which='SM')
    ev_m = np.sort(ev_m)

    H_m_s = build_hessian(phi_base, N, extra_diag=pert_minus_s)
    ev_m_s, _ = eigsh(H_m_s, k=4, which='SM')
    ev_m_s = np.sort(ev_m_s)

    V_without = ev_m[0] - ev_m_s[0]
    delta_V = V_full - V_without
    contrib_pct = delta_V / V_full * 100 if abs(V_full) > 1e-12 else 0

    report(f"  n={n_remove:3d}   {V_without:+14.8f} {delta_V:+14.8f} {contrib_pct:+10.2f}%")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report("1. INTRA-PROTON COUPLING:")
report(f"   Sum of individual shifts: {sum_individual:+.8f}")
report(f"   Actual all-7 shift:       {actual_shift:+.8f}")
na_frac = (actual_shift - sum_individual)/actual_shift * 100
report(f"   Non-additive fraction:    {na_frac:.1f}%")
if abs(na_frac) < 5:
    report("   → Modes are approximately INDEPENDENT within a single proton.")
elif abs(na_frac) < 20:
    report("   → Modes have WEAK coupling within a single proton.")
else:
    report("   → Modes are STRONGLY coupled within a single proton!")
report("")

report("2. INTER-PROTON DECAY RATES:")
report("   Wide modes (small n) should bond at larger R.")
for n in modes:
    eps_n = np.sin(n * gamma)
    splits = np.array(split_data[n])
    mask = splits > 1e-12
    if np.sum(mask) > 2:
        coeffs = np.polyfit(R_arr[mask], np.log(splits[mask]), 1)
        decay = -coeffs[0]
        report(f"   n={n}: width={1/eps_n:.1f} sites, decay_rate={decay:.4f}, "
               f"ratio/eps={decay/eps_n:.1f}")

report("")
report("3. CROSS-MODE COUPLING: see matrix above.")
report("   If off-diagonal entries << diagonal, modes bond independently.")
report("")
report("4. KEY QUESTION ANSWERED: are the 7 breather modes independent")
report("   channels for bonding, or do they form a coupled system?")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
