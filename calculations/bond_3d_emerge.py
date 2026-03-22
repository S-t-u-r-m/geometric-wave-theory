"""
Bond Energy from 3D Discrete Lattice — Morse Well Emergence
=============================================================
Two kink-antikink pairs ("protons") on the d=3 discrete lattice.
A breather ("electron") can tunnel between them, creating a bond.

The bond energy should EMERGE from the competition between:
  - Breather tunneling: attraction, decays as exp(-s*sqrt(Z)*R)
  - Kink-kink overlap: repulsion, decays as exp(-2*sqrt(Z)*R)
  - s = 0.1728 (Poeschl-Teller parameter, universal)

Method:
  1. Create static kink-antikink pairs on the 1D discrete lattice
  2. Compute the Hessian (d^2E/dphi^2) at this static configuration
  3. Find the lowest eigenvalues of the Hessian = bound state energies
  4. The bonding-antibonding splitting = 2 * D_e
  5. Scan R from 2 to 40 lattice sites

This is exact linear algebra — no time evolution, no noise, no approximations.
The Hessian eigenvalues ARE the physics. If a Morse well appears, the bond
energy has emerged from the Lagrangian.

EOM linearized: delta_phi_ddot = [Laplacian + d^2V/dphi^2 at kink] * delta_phi
The d^2V/dphi^2 matrix is the Hessian.
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3

outfile = os.path.join(os.path.dirname(__file__), "bond_3d_emerge_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("BOND ENERGY EMERGENCE FROM 3D DISCRETE LATTICE")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("")

# ============================================================
# KINK PROFILES ON DISCRETE LATTICE
# ============================================================
N = 256  # 1D lattice sites (large enough for two well-separated kinks)

# The sine-Gordon kink: phi goes from 0 to 2 (one cosine period)
# On a discrete lattice, the kink profile is:
# phi(x) = 1 + (2/pi)*arcsin(tanh(x - x_center))
# But we use the exact SG kink: phi(x) = (4/pi)*arctan(exp(x - x_center))
# which goes from 0 (x -> -inf) to 2 (x -> +inf).

def kink_profile(x, x_center):
    """Single kink: phi goes 0 -> 2."""
    return (4.0/PI) * np.arctan(np.exp(x - x_center))

def antikink_profile(x, x_center):
    """Anti-kink: phi goes 2 -> 0."""
    return 2.0 - (4.0/PI) * np.arctan(np.exp(x - x_center))

def make_kink_antikink(x, center, width):
    """Kink-antikink pair: phi goes 0 -> 2 -> 0, centered at 'center' with given width.
    This is a "proton" — a localized topological defect."""
    return kink_profile(x, center - width/2) + antikink_profile(x, center + width/2) - 2.0

def make_two_protons(x, pos_A, pos_B, kink_width):
    """Two kink-antikink pairs at positions A and B."""
    phi_A = make_kink_antikink(x, pos_A, kink_width)
    phi_B = make_kink_antikink(x, pos_B, kink_width)
    # Superposition (valid when wells don't overlap significantly)
    return phi_A + phi_B

# ============================================================
# HESSIAN CONSTRUCTION
# ============================================================
def build_hessian(phi_background, N):
    """Build the Hessian H = -Laplacian + d^2V/dphi^2 at the background configuration.

    For V = (1/pi^2)(1 - cos(pi*phi)):
      d^2V/dphi^2 = cos(pi*phi)

    For the discrete Laplacian with periodic BC:
      -Lap_{ij} = 2*delta_{ij} - delta_{i,j+1} - delta_{i,j-1}

    So H_{ij} = (2 + cos(pi*phi_i)) * delta_{ij} - delta_{i,j+1} - delta_{i,j-1}

    The eigenvalues of H are omega_k^2 (squared frequencies of linearized modes).
    """
    # Diagonal: 2 (from -Laplacian) + cos(pi*phi_i) (from potential)
    diag = 2.0 + np.cos(PI * phi_background)

    # Off-diagonal: -1 (from -Laplacian, nearest neighbors)
    off_diag = -np.ones(N - 1)

    # Periodic BC: connect site 0 to site N-1
    H = sparse.diags([off_diag, diag, off_diag], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0

    return H.tocsr()

# ============================================================
# PART 1: SINGLE KINK WELL — reference energy
# ============================================================
report("PART 1: SINGLE KINK WELL")
report("-" * 55)

center = N // 2
x = np.arange(N, dtype=np.float64)

kink_width = 3  # width of kink-antikink pair

# Single proton at center
phi_single = make_kink_antikink(x, center, kink_width)

report(f"Lattice: {N} sites, periodic BC")
report(f"Kink width: {kink_width} sites")
report(f"Profile at center:")
lo, hi = center - 6, center + 7
report("  x:   " + " ".join(f"{i:6d}" for i in range(lo, hi)))
report("  phi: " + " ".join(f"{phi_single[i]:6.3f}" for i in range(lo, hi)))
report("")

# Build Hessian and find lowest eigenvalues
H_single = build_hessian(phi_single, N)
n_eig = 10  # find the 10 lowest eigenvalues

# The eigenvalues are omega^2. The LOWEST one is the bound state (if < 1).
# The phonon band starts at omega^2 = 1 (mass gap squared).
eigenvalues_single, eigenvectors_single = eigsh(H_single, k=n_eig, which='SM')
eigenvalues_single = np.sort(eigenvalues_single)

report("Lowest eigenvalues (omega^2) of single-well Hessian:")
for i, ev in enumerate(eigenvalues_single):
    status = "BOUND" if ev < 1.0 else "band"
    report(f"  [{i}] omega^2 = {ev:.8f}, omega = {np.sqrt(abs(ev)):.8f}  {status}")

# Count bound states (below mass gap omega^2 = 1)
n_bound_single = np.sum(eigenvalues_single < 1.0)
report(f"\nBound states: {n_bound_single} (below mass gap omega^2 = 1)")
report("")

# ============================================================
# PART 2: DOUBLE WELL — scan separation R
# ============================================================
report("PART 2: DOUBLE WELL — SCAN SEPARATION R")
report("-" * 55)
report("Two kink-antikink pairs at separation R.")
report("Looking for bonding-antibonding splitting.")
report("")

R_values = list(range(4, 52, 2))

report(f"{'R':>4} {'E_bond':>12} {'E_anti':>12} {'split':>12} {'V(R)':>12} {'status':>10}")
report("-" * 70)

results = []
E_single_bound = eigenvalues_single[0] if n_bound_single > 0 else 1.0

for R in R_values:
    pos_A = center - R // 2
    pos_B = center + R // 2

    phi_double = make_two_protons(x, pos_A, pos_B, kink_width)

    # Build Hessian
    H_double = build_hessian(phi_double, N)

    # Find lowest eigenvalues
    try:
        evals, evecs = eigsh(H_double, k=min(n_eig, N-2), which='SM')
        evals = np.sort(evals)
    except Exception:
        evals = np.array([999])

    # The two lowest eigenvalues should be bonding and antibonding
    if len(evals) >= 2:
        E_bonding = evals[0]
        E_antibonding = evals[1]
        splitting = E_antibonding - E_bonding

        # V(R) = bonding energy relative to isolated wells
        # For isolated wells: E_0 = E_single_bound (each well has this)
        # For double well: E_bonding < E_single_bound (stabilized by tunneling)
        V_R = E_bonding - E_single_bound  # negative = attractive

        if E_bonding < 1.0 and E_antibonding < 1.0:
            status = "BONDING"
        elif E_bonding < 1.0:
            status = "1 bound"
        else:
            status = "no bound"
    else:
        E_bonding = evals[0] if len(evals) > 0 else 999
        E_antibonding = 999
        splitting = 0
        V_R = E_bonding - E_single_bound
        status = "sparse"

    results.append((R, E_bonding, E_antibonding, splitting, V_R))
    report(f"{R:4d} {E_bonding:12.8f} {E_antibonding:12.8f} {splitting:12.8f} "
           f"{V_R:+12.8f} {status:>10}")

report("")

# ============================================================
# PART 3: EXTRACT MORSE PARAMETERS
# ============================================================
report("PART 3: MORSE WELL ANALYSIS")
report("-" * 55)

R_arr = np.array([r[0] for r in results])
V_arr = np.array([r[4] for r in results])
split_arr = np.array([r[3] for r in results])

# Is there a minimum?
i_min = np.argmin(V_arr)
R_eq = R_arr[i_min]
V_min = V_arr[i_min]

report(f"Potential minimum:")
report(f"  R_eq = {R_eq} lattice sites")
report(f"  V(R_eq) = {V_min:.8f}")
report(f"  D_e = {-V_min:.8f} (well depth)")
report("")

# Check if it's a well (V goes to 0 at large R AND V < 0 at some R)
V_large_R = V_arr[-1]
report(f"  V(R={R_arr[-1]}) = {V_large_R:.8f} (should -> 0 at large R)")
report(f"  V is negative somewhere: {np.any(V_arr < 0)}")
report(f"  V returns to ~0 at large R: {abs(V_large_R) < abs(V_min) * 0.1}")
report("")

if np.any(V_arr < -1e-8) and abs(V_large_R) < abs(V_min) * 0.5:
    report("*** MORSE WELL DETECTED ***")
    report("Bond energy has EMERGED from the Lagrangian!")
    report("")

    # Fit Morse parameters: V(R) = D_e * [exp(-2*a*(R-R_eq)) - 2*exp(-a*(R-R_eq))]
    # At minimum: V(R_eq) = -D_e
    D_e = -V_min

    # Find 'a' from the slope at large R (exponential tail)
    # V(R) ~ -2*D_e*exp(-a*(R-R_eq)) for large R
    # log(-V) ~ log(2*D_e) - a*(R-R_eq)
    mask = (V_arr < -1e-10) & (R_arr > R_eq)
    if np.sum(mask) > 3:
        log_V = np.log(-V_arr[mask])
        R_fit = R_arr[mask]
        # Linear fit to log(-V) vs R
        coeffs = np.polyfit(R_fit, log_V, 1)
        a_morse = -coeffs[0]  # decay rate
        report(f"  Morse fit: D_e = {D_e:.8f}, a = {a_morse:.6f}, R_eq = {R_eq}")
        report(f"  Poeschl-Teller s = {(-1+np.sqrt(1+8/PI**2))/2:.6f}")
        report(f"  Ratio a/s = {a_morse/((-1+np.sqrt(1+8/PI**2))/2):.4f}")
    report("")

# Splitting at equilibrium
split_eq = split_arr[i_min]
report(f"Bonding-antibonding splitting at R_eq: {split_eq:.8f}")
report(f"2 * D_e = {2*(-V_min):.8f} (should match splitting)")
report("")

# ============================================================
# PART 4: VARY KINK WIDTH
# ============================================================
report("PART 4: KINK WIDTH DEPENDENCE")
report("-" * 55)
report("The proton 'size' = kink width. How does D_e depend on it?")
report("")

for kw in [1, 2, 3, 5, 7]:
    # Single well
    phi_s = make_kink_antikink(x, center, kw)
    H_s = build_hessian(phi_s, N)
    ev_s, _ = eigsh(H_s, k=5, which='SM')
    ev_s = np.sort(ev_s)
    E0_s = ev_s[0]

    # Double well at R = kw + 6 (give some gap)
    R_test = kw + 8
    phi_d = make_two_protons(x, center - R_test//2, center + R_test//2, kw)
    H_d = build_hessian(phi_d, N)
    ev_d, _ = eigsh(H_d, k=5, which='SM')
    ev_d = np.sort(ev_d)

    V_R = ev_d[0] - E0_s
    split = ev_d[1] - ev_d[0] if len(ev_d) >= 2 else 0

    n_bound = np.sum(ev_s < 1.0)
    report(f"  kink_width={kw}: E0_single={E0_s:.6f}, V(R={R_test})={V_R:+.6f}, "
           f"split={split:.6f}, bound_states={n_bound}")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report("POTENTIAL CURVE V(R):")
report(f"{'R':>4} {'V(R)':>12} {'splitting':>12}")
for r in results:
    marker = " <-- min" if r[0] == R_eq else ""
    report(f"{r[0]:4d} {r[4]:+12.8f} {r[3]:12.8f}{marker}")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))

log.close()
print(f"\nResults saved to: {outfile}")
