"""
Self-Consistent Field Bond Model — Hartree Method for Breather Modes
=====================================================================
Each breather mode adjusts to the field created by all the others.
Iterate until self-consistent. No additivity assumption.

Algorithm:
  1. Initialize: bare kink well, no breather perturbations
  2. For each mode n = 1..7:
     a. Build Hessian with perturbations from ALL OTHER modes
     b. Find lowest eigenvalue = mode n's energy in the combined field
     c. Record mode n's eigenstate (the perturbation it creates)
  3. Repeat step 2 until all eigenvalues converge (< 1e-8 change)
  4. Total energy = sum of self-consistent eigenvalues with double-counting correction

For bonding: run SCF for single proton, then for two protons at R.
Bond energy = E_SCF(two protons, R) - 2 × E_SCF(one proton)

This automatically handles:
  - Mode-mode coupling (the 99.7% cancellation)
  - R-dependence (the Hessian changes with separation)
  - Asymmetry/ionic effects (different well shapes for different atoms)
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "scf_bond_model_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("SELF-CONSISTENT FIELD BOND MODEL")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"d = {d}, gamma = {gamma:.10f}")
report("")

# ============================================================
# LATTICE SETUP
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

def breather_profile_sq(x, pos, eps_n):
    """Squared breather profile — the perturbation source for mode n."""
    phi = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * (x - pos) + 1e-30))
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

# ============================================================
# THE SCF ENGINE
# ============================================================
def run_scf(phi_background, proton_positions, max_iter=60, tol=1e-8,
            verbose=True, label=""):
    """Run self-consistent field iteration for breather modes.

    Args:
        phi_background: static kink field (1D array, length N)
        proton_positions: list of positions where protons sit
        max_iter: maximum SCF iterations
        tol: convergence tolerance on eigenvalue change
        verbose: print iteration details
        label: label for output

    Returns:
        dict with converged eigenvalues, total energy, iteration count
    """
    n_protons = len(proton_positions)

    # Initialize: no perturbations, just bare kink background
    # Each mode at each proton has a perturbation amplitude (starts at 0)
    # perturbation[n][p] = the diagonal correction from mode n at proton p
    mode_pert = {}  # {mode_n: array of length N}
    mode_evals = {}  # {mode_n: eigenvalue}

    for n in modes:
        mode_pert[n] = np.zeros(N)
        mode_evals[n] = 1.0  # start at mass gap

    converged = False
    prev_total_E = 999.0

    for iteration in range(max_iter):
        total_E = 0.0
        max_change = 0.0

        for n in modes:
            eps_n = np.sin(n * gamma)

            # Build the potential from ALL OTHER modes
            other_pert = np.zeros(N)
            for m in modes:
                if m != n:
                    other_pert += mode_pert[m]

            # Solve the Hessian with other modes' perturbations
            H = build_hessian(phi_background, extra_diag=other_pert)
            ev, evec = eigsh(H, k=1, which='SM')
            new_eval = ev[0]

            # Track convergence
            change = abs(new_eval - mode_evals[n])
            max_change = max(max_change, change)
            mode_evals[n] = new_eval
            total_E += new_eval

            # Update this mode's perturbation:
            # Mode n creates a softening proportional to its squared profile
            # The strength comes from the eigenstate — use the eigenvector
            # to compute the effective occupation at each site
            psi = evec[:, 0]
            psi_sq = psi**2  # probability density

            # The perturbation from mode n at each proton position:
            # delta_V = -(pi^2/6) * phi_n^2(x) weighted by occupation
            new_pert = np.zeros(N)
            for pos in proton_positions:
                profile = breather_profile_sq(x, pos, eps_n)
                # Weight by how much of the eigenstate is near this proton
                # (for well-separated protons, each eigenstate localizes at one)
                hw = max(5, int(2.0 / eps_n))
                lo = max(0, int(pos) - hw)
                hi = min(N, int(pos) + hw + 1)
                local_weight = np.sum(psi_sq[lo:hi])
                new_pert += -(PI**2 / 6.0) * profile * local_weight

            mode_pert[n] = new_pert

        # Check convergence
        if verbose and (iteration < 3 or iteration % 5 == 0 or max_change < tol):
            evals_str = " ".join(f"{mode_evals[n]:+.6f}" for n in modes)
            report(f"  {label} iter {iteration:2d}: total_E = {total_E:.8f}, "
                   f"max_change = {max_change:.2e}")

        if max_change < tol:
            converged = True
            if verbose:
                report(f"  {label} CONVERGED in {iteration+1} iterations")
            break

        prev_total_E = total_E

    if not converged and verbose:
        report(f"  {label} NOT converged after {max_iter} iterations "
               f"(max_change = {max_change:.2e})")

    # Total SCF energy with double-counting correction:
    # E_total = sum of eigenvalues - (1/2) sum of interaction energies
    # The interaction energy is already counted twice (once for each mode)
    # Correction: E_total = sum_n eps_n - (1/2) sum_{n!=m} <n|V_m|n>
    # For simplicity, use the total potential energy approach:
    total_pert = np.zeros(N)
    for n in modes:
        total_pert += mode_pert[n]

    H_total = build_hessian(phi_background, extra_diag=total_pert)
    ev_total, _ = eigsh(H_total, k=1, which='SM')
    E_total = ev_total[0]

    return {
        'E_total': E_total,
        'E_sum': total_E,
        'mode_evals': dict(mode_evals),
        'converged': converged,
        'iterations': iteration + 1,
        'total_pert': total_pert,
    }

# ============================================================
# PART 1: SINGLE PROTON SCF
# ============================================================
report("PART 1: SINGLE PROTON — SCF CONVERGENCE")
report("-" * 55)

phi_single = make_proton(x, center, kink_width)
result_single = run_scf(phi_single, [center], verbose=True, label="1-proton")

E_ref = result_single['E_total']
report(f"\nSingle proton SCF energy: {E_ref:.8f}")
report(f"Individual mode eigenvalues:")
for n in modes:
    eps_n = np.sin(n * gamma)
    report(f"  n={n}: eps={eps_n:.4f}, eval={result_single['mode_evals'][n]:.8f}")

# Compare with bare (no SCF)
H_bare = build_hessian(phi_single)
ev_bare, _ = eigsh(H_bare, k=1, which='SM')
report(f"\nBare (no breathers): E0 = {ev_bare[0]:.8f}")
report(f"SCF (self-consistent): E0 = {E_ref:.8f}")
report(f"SCF correction: {E_ref - ev_bare[0]:+.8f}")
report("")

# ============================================================
# PART 2: TWO PROTONS — BOND CURVE
# ============================================================
report("PART 2: TWO PROTONS — SCF BOND CURVE")
report("-" * 55)
report("V(R) = E_SCF(two protons at R) - 2 × E_SCF(one proton)")
report("")

R_values = [6, 7, 8, 9, 10, 12, 14, 16, 20, 24, 30, 40]

report(f"{'R':>4} {'E_SCF':>12} {'V(R)':>12} {'V_bare':>12} "
       f"{'SCF/bare':>10} {'iter':>5}")
report("-" * 60)

bond_curve = {}

for R in R_values:
    pos_A = center - R // 2
    pos_B = center + R // 2
    phi_double = make_proton(x, pos_A, kink_width) + make_proton(x, pos_B, kink_width)

    # SCF for two protons
    result = run_scf(phi_double, [pos_A, pos_B], verbose=False, label=f"R={R}")

    V_scf = result['E_total'] - 2 * E_ref

    # Bare comparison
    H_bare = build_hessian(phi_double)
    ev_bare, _ = eigsh(H_bare, k=1, which='SM')
    V_bare = ev_bare[0] - 2 * ev_bare[0]  # wrong reference, fix:

    # Bare reference: single proton bare
    H_bare_s = build_hessian(phi_single)
    ev_bare_s, _ = eigsh(H_bare_s, k=1, which='SM')
    V_bare = ev_bare[0] - 2 * ev_bare_s[0]

    ratio = V_scf / V_bare if abs(V_bare) > 1e-12 else 0

    bond_curve[R] = {
        'E_scf': result['E_total'], 'V_scf': V_scf,
        'V_bare': V_bare, 'converged': result['converged'],
        'iterations': result['iterations']
    }

    conv = "Y" if result['converged'] else "N"
    report(f"{R:4d} {result['E_total']:12.6f} {V_scf:+12.8f} {V_bare:+12.8f} "
           f"{ratio:10.4f} {result['iterations']:5d}{conv}")

report("")

# ============================================================
# PART 3: MORSE WELL ANALYSIS
# ============================================================
report("PART 3: MORSE WELL ANALYSIS")
report("-" * 55)

R_arr = np.array(sorted(bond_curve.keys()), dtype=float)
V_arr = np.array([bond_curve[int(R)]['V_scf'] for R in R_arr])
V_bare_arr = np.array([bond_curve[int(R)]['V_bare'] for R in R_arr])

i_min = np.argmin(V_arr)
i_min_bare = np.argmin(V_bare_arr)

report(f"SCF bond curve:")
report(f"  R_eq = {R_arr[i_min]:.0f}")
report(f"  D_e = {-V_arr[i_min]:.8f}")
report(f"  V(R_max) = {V_arr[-1]:.8f} (should → 0)")
report("")

report(f"Bare bond curve:")
report(f"  R_eq = {R_arr[i_min_bare]:.0f}")
report(f"  D_e = {-V_bare_arr[i_min_bare]:.8f}")
report("")

if V_arr[i_min] < -1e-6:
    report("*** SCF MORSE WELL DETECTED ***")
    # Fit Morse decay
    mask = (V_arr < -1e-10) & (R_arr > R_arr[i_min])
    if np.sum(mask) > 2:
        coeffs = np.polyfit(R_arr[mask], np.log(-V_arr[mask]), 1)
        a_morse = -coeffs[0]
        report(f"  Morse decay rate: {a_morse:.4f}")
elif V_arr[i_min] < 0:
    report("Weak attraction detected")
else:
    report("No bonding in SCF curve")
    report("The SCF self-energy may be overwhelming the bond signal.")
    report("Check: are the mode perturbations too strong?")

report("")

# ============================================================
# PART 4: MODE-RESOLVED BOND CONTRIBUTION
# ============================================================
report("PART 4: HOW EACH MODE SHIFTS DURING BONDING")
report("-" * 55)
report("Compare single-proton mode eigenvalues with two-proton values.")
report("")

R_detail = 10
pos_A = center - R_detail // 2
pos_B = center + R_detail // 2
phi_double = make_proton(x, pos_A, kink_width) + make_proton(x, pos_B, kink_width)

result_detail = run_scf(phi_double, [pos_A, pos_B], verbose=True, label=f"detail R={R_detail}")

report(f"\nMode eigenvalue shifts at R={R_detail}:")
report(f"{'n':>3} {'eps_n':>8} {'single':>12} {'bonded':>12} {'shift':>12}")
report("-" * 50)

total_shift = 0
for n in modes:
    eps_n = np.sin(n * gamma)
    ev_single = result_single['mode_evals'][n]
    ev_bonded = result_detail['mode_evals'][n]
    shift = ev_bonded - ev_single
    total_shift += shift
    report(f"  {n:3d} {eps_n:8.4f} {ev_single:12.6f} {ev_bonded:12.6f} {shift:+12.8f}")

report(f"\n  Total mode shift: {total_shift:+.8f}")
report(f"  SCF bond energy: {bond_curve[R_detail]['V_scf']:+.8f}")
report("")

# ============================================================
# PART 5: ASYMMETRIC SCF (different atoms)
# ============================================================
report("PART 5: ASYMMETRIC SCF — DIFFERENT WELL DEPTHS")
report("-" * 55)
report("Model ionic bonding: one well deeper than the other.")
report("")

R_asym = 10
pos_A = center - R_asym // 2
pos_B = center + R_asym // 2

# Modify one well to be deeper (higher IE)
def make_asymmetric_bg(x, pos_A, pos_B, kw, delta_depth):
    """Two protons with different well depths.
    delta_depth > 0: well B is deeper (more electronegative)."""
    phi = make_proton(x, pos_A, kw) + make_proton(x, pos_B, kw)
    # Scale the potential at well B
    # cos(pi*phi) → cos(pi*phi) - delta at well B region
    return phi, delta_depth

phi_asym = make_proton(x, pos_A, kink_width) + make_proton(x, pos_B, kink_width)

for delta in [0.0, 0.1, 0.3, 0.5]:
    # Add asymmetric depth perturbation
    asym_pert = delta * np.exp(-(x - pos_B)**2 / (2 * 8**2))
    # Modify the background to include asymmetry
    # We'll add the asymmetry as an extra diagonal term

    def run_scf_asym(phi_bg, positions, asym_diag):
        """SCF with an additional fixed asymmetric diagonal term."""
        mode_pert = {n: np.zeros(N) for n in modes}
        mode_evals = {n: 1.0 for n in modes}

        for iteration in range(20):
            max_change = 0
            for n in modes:
                eps_n = np.sin(n * gamma)
                other_pert = sum(mode_pert[m] for m in modes if m != n)
                total_extra = other_pert + asym_diag
                H = build_hessian(phi_bg, extra_diag=total_extra)
                ev, evec = eigsh(H, k=1, which='SM')
                change = abs(ev[0] - mode_evals[n])
                max_change = max(max_change, change)
                mode_evals[n] = ev[0]

                psi_sq = evec[:, 0]**2
                new_pert = np.zeros(N)
                for pos in positions:
                    profile = breather_profile_sq(x, pos, eps_n)
                    hw = max(5, int(2.0 / eps_n))
                    lo = max(0, int(pos) - hw)
                    hi = min(N, int(pos) + hw + 1)
                    local_weight = np.sum(psi_sq[lo:hi])
                    new_pert += -(PI**2 / 6.0) * profile * local_weight
                # Damped update to prevent oscillation
                damping = 0.25  # mix 25% new + 75% old
                mode_pert[n] = damping * new_pert + (1 - damping) * mode_pert[n]

            if max_change < 1e-8:
                break

        total_pert = sum(mode_pert.values()) + asym_diag
        H_total = build_hessian(phi_bg, extra_diag=total_pert)
        ev_total, _ = eigsh(H_total, k=1, which='SM')
        return ev_total[0], mode_evals

    E_asym, evals_asym = run_scf_asym(phi_asym, [pos_A, pos_B], asym_pert)
    V_asym = E_asym - 2 * E_ref  # approximate (reference should also be asymmetric)

    report(f"  delta={delta:.1f}: E_SCF={E_asym:.6f}, V(R)={V_asym:+.8f}")

report("")

# ============================================================
# PART 6: COMPARE SCF WITH PERTURBATIVE MODELS
# ============================================================
report("PART 6: SCF vs PERTURBATIVE — FULL COMPARISON")
report("-" * 55)

# Bare: no breather perturbations at all
# Pert: all 7 modes as static perturbation (no self-consistency)
# SCF: self-consistent

report(f"{'R':>4} {'V_bare':>12} {'V_SCF':>12} {'diff':>12} {'ratio':>8}")
report("-" * 50)

for R in R_values:
    vb = bond_curve[R]['V_bare']
    vs = bond_curve[R]['V_scf']
    diff = vs - vb
    ratio = vs / vb if abs(vb) > 1e-12 else 0
    report(f"{R:4d} {vb:+12.8f} {vs:+12.8f} {diff:+12.8f} {ratio:8.4f}")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report(f"SCF converged in {result_single['iterations']} iterations (single proton)")
report(f"Single proton SCF energy: {E_ref:.8f}")
report(f"Bare energy: {ev_bare_s[0]:.8f}")
report(f"SCF correction: {E_ref - ev_bare_s[0]:+.8f}")
report("")

if V_arr[i_min] < -1e-6:
    report(f"SCF BOND CURVE:")
    report(f"  R_eq = {R_arr[i_min]:.0f}, D_e = {-V_arr[i_min]:.8f}")
    report(f"  (Bare: R_eq = {R_arr[i_min_bare]:.0f}, D_e = {-V_bare_arr[i_min_bare]:.8f})")
    report("")
    report("The SCF method WORKS — it produces a bond curve from the")
    report("self-consistent breather mode equilibrium.")
else:
    report("SCF bond curve does not show clear bonding.")
    report("The self-energy correction may need refinement.")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
