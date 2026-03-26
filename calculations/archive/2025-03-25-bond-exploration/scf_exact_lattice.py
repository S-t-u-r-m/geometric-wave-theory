"""
SCF with Exact Lattice Potential — No Approximations
======================================================
Previous attempts used -(π²/6)φ² as the breather perturbation.
This is a Taylor expansion that overshoots the real lattice potential.

The FIX: use the EXACT cosine potential of the lattice.

The field at each site is:
  φ_total(x) = φ_kink(x) + Σ_n a_n × φ_breather_n(x)

The Hessian diagonal is EXACTLY:
  H_ii = 2 + cos(π × φ_total(x_i))

No approximation. The cosine IS the lattice potential.
This automatically enforces:
  - Boundedness (|cos| ≤ 1, so diagonal ∈ [1, 3])
  - Periodicity (the cosine repeats every Δφ = 2)
  - The exact barrier height and well shape
  - All higher-order coupling terms

The optimizer finds amplitudes a_n ∈ [0, 1] that minimize
the lowest eigenvalue of this exact Hessian.

GPU: CuPy sparse eigsh for speed.
"""
import sys, io, os, time
import numpy as np
from scipy.optimize import minimize

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

# Try GPU, fall back to CPU
try:
    import cupy as cp
    import cupyx.scipy.sparse as csp
    import cupyx.scipy.sparse.linalg as csla
    from scipy import sparse as sp_sparse
    GPU = True
    print("GPU active (CuPy)")
except Exception:
    from scipy import sparse as sp_sparse
    from scipy.sparse.linalg import eigsh as cpu_eigsh
    GPU = False
    print("CPU mode (SciPy)")

outfile = os.path.join(os.path.dirname(__file__), "scf_exact_lattice_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("SCF WITH EXACT LATTICE POTENTIAL")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"d = {d}, gamma = {gamma:.10f}")
report(f"GPU: {GPU}")
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
    """Static kink-antikink field."""
    return (4.0/PI) * (np.arctan(np.exp(x - pos + kw/2))
                      - np.arctan(np.exp(x - pos - kw/2)))

def breather_field(x, pos, eps_n):
    """Breather peak profile (the field displacement, not squared)."""
    return (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * (x - pos) + 1e-30))

def build_hessian_exact(phi_total):
    """Build Hessian from EXACT lattice potential.
    H_ii = 2 + cos(π × φ_total_i)
    H_{i,i±1} = -1 (nearest neighbor, periodic BC)
    """
    diag = 2.0 + np.cos(PI * phi_total)
    off = -np.ones(N - 1)
    H = sp_sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()

def get_E0(phi_total):
    """Lowest eigenvalue of the exact Hessian."""
    H_cpu = build_hessian_exact(phi_total)
    if GPU:
        H_gpu = csp.csr_matrix(H_cpu)
        ev, _ = csla.eigsh(H_gpu, k=1, which='SA')
        return float(ev[0])
    else:
        ev, _ = cpu_eigsh(H_cpu, k=1, which='SM', maxiter=3000)
        return ev[0]

# Pre-compute breather field profiles
def precompute_breather_fields(positions):
    """Pre-compute breather field profiles for all modes at given positions."""
    fields = {}
    for n in modes:
        eps_n = np.sin(n * gamma)
        f = np.zeros(N)
        for pos in positions:
            f += breather_field(x, pos, eps_n)
        fields[n] = f
    return fields

def total_field(phi_kink, breather_fields, a_vec):
    """Total field = kink + Σ a_n × breather_n. No approximation."""
    phi = phi_kink.copy()
    for i, n in enumerate(modes):
        phi += a_vec[i] * breather_fields[n]
    return phi

# ============================================================
# ENERGY FUNCTION (exact lattice)
# ============================================================
eval_count = [0]

def energy_exact(a_vec, phi_kink, breather_fields):
    """Total energy from exact lattice Hessian."""
    eval_count[0] += 1
    phi = total_field(phi_kink, breather_fields, a_vec)
    return get_E0(phi)

# ============================================================
# PART 1: SINGLE PROTON
# ============================================================
report("PART 1: SINGLE PROTON — EXACT LATTICE OPTIMIZATION")
report("-" * 55)

phi_kink_single = kink_field(x, center, kink_width)
bf_single = precompute_breather_fields([center])

# Bare energy (no breathers)
E_bare = get_E0(phi_kink_single)
report(f"Bare (no breathers): E0 = {E_bare:.8f}")

# Optimize with bounded amplitudes [0, 1]
# Physical constraint: breather can't exceed one full oscillation
bounds = [(-1.0, 1.0)] * n_modes

a0 = np.zeros(n_modes)
eval_count[0] = 0
t0 = time.time()

# Use Powell (supports bounds) for efficiency
result_s = minimize(
    energy_exact, a0,
    args=(phi_kink_single, bf_single),
    method='Powell',
    bounds=bounds,
    options={'maxiter': 3000, 'ftol': 1e-10}
)

t_single = time.time() - t0
E_single = result_s.fun
a_single = result_s.x

report(f"Optimized: E0 = {E_single:.8f}")
report(f"SCF correction: {E_single - E_bare:+.8f}")
report(f"({eval_count[0]} evaluations in {t_single:.1f}s)")
report("")

report("Optimal mode amplitudes:")
report(f"{'n':>3} {'eps_n':>8} {'a_n':>10} {'field_max':>10}")
report("-" * 35)
for i, n in enumerate(modes):
    eps_n = np.sin(n * gamma)
    fmax = np.max(np.abs(a_single[i] * bf_single[n]))
    report(f"  {n:3d} {eps_n:8.4f} {a_single[i]:+10.6f} {fmax:10.6f}")

# Show the total field modification
phi_total_single = total_field(phi_kink_single, bf_single, a_single)
report(f"\nTotal field range: [{phi_total_single.min():.4f}, {phi_total_single.max():.4f}]")
report(f"Kink field range: [{phi_kink_single.min():.4f}, {phi_kink_single.max():.4f}]")
report(f"Max field modification: {np.max(np.abs(phi_total_single - phi_kink_single)):.6f}")
report("")

# ============================================================
# PART 2: BOND CURVE
# ============================================================
report("PART 2: BOND CURVE — V(R)")
report("-" * 55)

R_values = [4, 6, 7, 8, 9, 10, 12, 14, 16, 20, 24, 30, 40, 60, 80]

report(f"{'R':>4} {'E_double':>12} {'V(R)':>12} {'V_bare':>12} "
       f"{'ratio':>8} {'evals':>6} {'time':>6}")
report("-" * 65)

bond_curve = {}
a_bonded = {}

for R in R_values:
    pos_A = center - R // 2
    pos_B = center + R // 2

    phi_kink_double = kink_field(x, pos_A, kink_width) + kink_field(x, pos_B, kink_width)
    bf_double = precompute_breather_fields([pos_A, pos_B])

    # Bare
    E_bare_d = get_E0(phi_kink_double)
    V_bare = E_bare_d - 2 * E_bare

    # Optimize starting from single-proton solution
    eval_count[0] = 0
    t0 = time.time()

    result_d = minimize(
        energy_exact, a_single,
        args=(phi_kink_double, bf_double),
        method='Powell',
        bounds=bounds,
        options={'maxiter': 3000, 'ftol': 1e-10}
    )

    t_bond = time.time() - t0
    E_double = result_d.fun
    V_scf = E_double - 2 * E_single

    ratio = V_scf / V_bare if abs(V_bare) > 1e-12 else 0

    bond_curve[R] = {'E': E_double, 'V': V_scf, 'V_bare': V_bare}
    a_bonded[R] = result_d.x.copy()

    report(f"{R:4d} {E_double:12.8f} {V_scf:+12.8f} {V_bare:+12.8f} "
           f"{ratio:8.4f} {eval_count[0]:6d} {t_bond:5.1f}s")

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

report(f"SCF: R_eq = {R_arr[i_min]:.0f}, V_min = {V_arr[i_min]:+.8f}")
report(f"Bare: R_eq = {R_arr[i_min_bare]:.0f}, V_min = {V_bare_arr[i_min_bare]:+.8f}")
report(f"V(R→∞) = {V_arr[-1]:+.8f} (should → 0)")
report("")

if V_arr[i_min] < -1e-6:
    D_e = -V_arr[i_min]
    report(f"*** MORSE WELL DETECTED ***")
    report(f"  D_e = {D_e:.8f}")
    report(f"  R_eq = {R_arr[i_min]:.0f}")

    # Decay rate
    mask = (V_arr < -1e-8) & (R_arr > R_arr[i_min])
    if np.sum(mask) > 2:
        coeffs = np.polyfit(R_arr[mask], np.log(-V_arr[mask]), 1)
        report(f"  Morse decay: {-coeffs[0]:.4f}")
    report("")

    # Compare with bare
    D_e_bare = -V_bare_arr[i_min_bare]
    report(f"  D_e(SCF) / D_e(bare) = {D_e / D_e_bare:.4f}" if D_e_bare > 0 else "")
else:
    report("No bonding detected in SCF curve.")
    report("Checking if V(R) approaches 0 at large R...")
    for R in [40, 60, 80]:
        if R in bond_curve:
            report(f"  V({R}) = {bond_curve[R]['V']:+.8f}")

report("")

# ============================================================
# PART 4: MODE WEIGHT EVOLUTION WITH R
# ============================================================
report("PART 4: MODE AMPLITUDES VS R")
report("-" * 55)

header = f"{'R':>4}"
for n in modes:
    header += f"  {'a_'+str(n):>8}"
report(header)
report("-" * (4 + n_modes * 10))

for R in R_values:
    line = f"{R:4d}"
    for i in range(n_modes):
        line += f"  {a_bonded[R][i]:+8.5f}"
    report(line)

report("")
report("Single proton reference:")
line = " ref"
for i in range(n_modes):
    line += f"  {a_single[i]:+8.5f}"
report(line)

report("")

# How stable are the weights? Do they converge to the single-proton values?
report("CONVERGENCE: |a(R) - a(single)| for each mode")
report(f"{'R':>4}", end="")
for n in modes:
    print(f"  {'d'+str(n):>8}", end="")
print()
log.write(f"{'R':>4}" + "".join(f"  {'d'+str(n):>8}" for n in modes) + "\n")
report("-" * (4 + n_modes * 10))

for R in R_values:
    line = f"{R:4d}"
    for i in range(n_modes):
        da = abs(a_bonded[R][i] - a_single[i])
        line += f"  {da:8.5f}"
    report(line)

report("")

# ============================================================
# PART 5: VERIFY THE TOTAL FIELD IS PHYSICAL
# ============================================================
report("PART 5: TOTAL FIELD VERIFICATION")
report("-" * 55)
report("The total field should stay within physical bounds.")
report("")

for R in [8, 10, 14, 40]:
    if R in bond_curve:
        pos_A = center - R // 2
        pos_B = center + R // 2
        phi_kink = kink_field(x, pos_A, kink_width) + kink_field(x, pos_B, kink_width)
        bf = precompute_breather_fields([pos_A, pos_B])
        phi_total = total_field(phi_kink, bf, a_bonded[R])

        report(f"  R={R}: phi_total ∈ [{phi_total.min():.4f}, {phi_total.max():.4f}], "
               f"kink ∈ [{phi_kink.min():.4f}, {phi_kink.max():.4f}], "
               f"max_mod = {np.max(np.abs(phi_total - phi_kink)):.4f}")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report("EXACT LATTICE POTENTIAL — no Taylor expansion.")
report("Hessian diagonal = 2 + cos(π × (φ_kink + Σ aₙφₙ))")
report("Amplitudes bounded: |aₙ| ≤ 1")
report("")
report(f"Single proton: E = {E_single:.8f} (bare: {E_bare:.8f})")
report(f"SCF correction: {E_single - E_bare:+.8f}")
report("")

report("Bond curve V(R):")
for R in R_values:
    v = bond_curve[R]['V']
    vb = bond_curve[R]['V_bare']
    marker = " <-- min" if R == R_arr[i_min] else ""
    report(f"  R={R:3d}: V_scf={v:+.8f}, V_bare={vb:+.8f}{marker}")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
