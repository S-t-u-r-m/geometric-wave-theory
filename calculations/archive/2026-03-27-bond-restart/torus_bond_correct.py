"""
Two-Torus Bonding — Correct Poloidal Winding Topology
=======================================================
Two topologically correct tori (poloidal winding) on the 3D lattice.
Track how the 40 bound modes split as the tori approach.

This is the FIRST bond calculation with the correct proton model:
  - No tachyons (all eigenvalues positive)
  - Topological stability (can't annihilate)
  - Real confinement

The two tori are coplanar (both rings in xy-plane), offset along x.
Separation R = distance between ring centers.

GPU: RTX 4070 Ti. N=64, each eigsh ~2-3s.
"""
import sys, io, os, time
import numpy as np
import cupy as cp
import cupyx.scipy.sparse as csp
import cupyx.scipy.sparse.linalg as csla
from scipy import sparse

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3

outfile = os.path.join(os.path.dirname(__file__), "torus_bond_correct_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("TWO-TORUS BONDING — CORRECT POLOIDAL WINDING")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("")

N = 64
Ntot = N**3
R_maj = 8
r_tube = 3
kink_K = r_tube  # kink sharpness

report(f"Lattice: {N}^3 = {Ntot:,} sites")
report(f"Torus: R_major = {R_maj}, r_tube = {r_tube}")
report("")

ix = np.arange(N, dtype=np.float64) - N/2

def make_poloidal_torus(X, Y, Z, center_x, center_y):
    """Build a single torus with poloidal winding at (center_x, center_y, 0)."""
    Xc = X - center_x
    Yc = Y - center_y
    rho_xy = np.sqrt(Xc**2 + Yc**2 + 1e-30)
    rho_tube = np.sqrt((rho_xy - R_maj)**2 + Z**2)
    phi_pol = np.arctan2(Z, rho_xy - R_maj)

    # Kink profile: 0→2 around poloidal direction
    field = (4.0/PI) * np.arctan(np.exp(kink_K * phi_pol))

    # Cutoff: snap to vacuum outside tube region
    cutoff = r_tube + 3
    mask_far_upper = (rho_tube >= cutoff) & (phi_pol >= 0)
    mask_far_lower = (rho_tube >= cutoff) & (phi_pol < 0)
    field[mask_far_upper] = 2.0
    field[mask_far_lower] = 0.0

    return field

def build_hessian_gpu(phi_flat):
    diag = 2.0 * d + np.cos(PI * phi_flat)
    idx = np.arange(Ntot, dtype=np.int32)
    ix_a = idx % N
    iy_a = (idx // N) % N
    iz_a = idx // (N * N)
    neighbors = []
    for dx, dy, dz in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
        jx = (ix_a + dx) % N
        jy = (iy_a + dy) % N
        jz = (iz_a + dz) % N
        j = (jz * N * N + jy * N + jx).astype(np.int32)
        neighbors.append(j)
    all_rows = np.concatenate([idx] * 7)
    all_cols = np.concatenate([idx] + neighbors)
    all_vals = np.concatenate([diag] + [-np.ones(Ntot, dtype=np.float64)] * 6)
    H_cpu = sparse.coo_matrix((all_vals, (all_rows, all_cols)),
                               shape=(Ntot, Ntot)).tocsr()
    return csp.csr_matrix(H_cpu)

# ============================================================
# SINGLE TORUS REFERENCE
# ============================================================
report("SINGLE TORUS REFERENCE (poloidal winding)")
report("-" * 55)

X3, Y3, Z3 = np.meshgrid(ix, ix, ix, indexing='ij')
phi_single = make_poloidal_torus(X3, Y3, Z3, 0, 0)

H_gpu = build_hessian_gpu(phi_single.ravel())
n_eig = 20
evals_s, evecs_s = csla.eigsh(H_gpu, k=n_eig, which='SA')
cp.cuda.Stream.null.synchronize()
evals_s = cp.asnumpy(evals_s)
evals_s = np.sort(evals_s)

n_neg = np.sum(evals_s < 0)
report(f"Negative eigenvalues: {n_neg} (should be 0)")
report(f"Lowest 10 eigenvalues:")
for i in range(min(10, len(evals_s))):
    report(f"  [{i}] {evals_s[i]:+.6f}")
report("")

# ============================================================
# TWO TORI — SCAN SEPARATION
# ============================================================
report("TWO TORI — EIGENVALUE SPLITTING VS SEPARATION")
report("-" * 55)
report("Two tori in the xy-plane, offset along x-axis.")
report("R = distance between ring centers.")
report(f"Tube-to-tube contact at R = 2*R_maj = {2*R_maj}")
report(f"Tubes overlap when R < 2*R_maj + 2*r_tube = {2*R_maj + 2*r_tube}")
report("")

# R must be > 2*R_maj for the rings not to intersect
# Tubes start touching at R ≈ 2*R_maj + 2*r_tube = 22
# Tubes strongly overlap at R ≈ 2*R_maj = 16
# For N=64: max R before periodic image ≈ 32 - 2*R_maj = 16... tight.
# We need R from about 18 to 30.

# Actually, let me reconsider. The tori are in the xy-plane.
# Torus A centered at (-R/2, 0), torus B at (+R/2, 0).
# The rightmost point of A's ring is at (-R/2 + R_maj, 0) = (R_maj - R/2, 0)
# The leftmost point of B's ring is at (R/2 - R_maj, 0)
# They meet when R_maj - R/2 = R/2 - R_maj → R = 2*R_maj = 16
# Tubes touch when the gap between ring edges = 2*r_tube
# Gap = R - 2*R_maj, tubes touch when gap ≈ 2*r_tube = 6
# So tubes touch at R ≈ 22.
# For N=64 (box = 64 sites), max R ≈ 30 to avoid periodic images.

R_values = [30, 28, 26, 24, 22, 21, 20, 19, 18, 17]

report(f"R range: {R_values[0]} to {R_values[-1]}")
report(f"Ring gap = R - 2*R_maj. Tubes touch at gap ≈ 2*r_tube = {2*r_tube}")
report("")

header = f"{'R':>3} {'gap':>4} {'n_neg':>5}"
for i in range(6):
    header += f" {'ev'+str(i):>10}"
report(header)
report("-" * (12 + 6*11))

all_results = {}

for R in R_values:
    t0 = time.time()

    # Place two tori along x-axis
    phi_A = make_poloidal_torus(X3, Y3, Z3, -R/2, 0)
    phi_B = make_poloidal_torus(X3, Y3, Z3, +R/2, 0)

    # Superposition (valid when tubes don't overlap much)
    # Both fields go 0→2 poloidally. In the overlap region,
    # the fields add, giving φ > 2. But cos(πφ) treats φ=2 as vacuum.
    # For well-separated tori, the overlap is minimal.
    phi_double = phi_A + phi_B

    # Clip to [0, 4] for safety (shouldn't matter far from overlap)
    phi_double = np.clip(phi_double, 0, 4)

    H_gpu = build_hessian_gpu(phi_double.ravel())

    evals_d, evecs_d = csla.eigsh(H_gpu, k=n_eig, which='SA')
    cp.cuda.Stream.null.synchronize()
    evals_d = cp.asnumpy(evals_d)
    evals_d = np.sort(evals_d)

    gap = R - 2*R_maj
    n_neg_d = np.sum(evals_d < -0.01)

    all_results[R] = evals_d.copy()

    line = f"{R:3d} {gap:4d} {n_neg_d:5d}"
    for i in range(min(6, len(evals_d))):
        line += f" {evals_d[i]:+10.6f}"
    elapsed = time.time() - t0
    report(f"{line}  [{elapsed:.1f}s]")

report("")

# ============================================================
# SPLITTING ANALYSIS
# ============================================================
report("SPLITTING ANALYSIS")
report("-" * 55)
report("At large R: each single-torus eigenvalue appears twice.")
report("As R decreases: pairs split into bonding + antibonding.")
report("")

# Pair up eigenvalues at largest R
ev_ref = all_results[R_values[0]]
pairs = []
i = 0
while i < len(ev_ref) - 1:
    if abs(ev_ref[i+1] - ev_ref[i]) < 0.02:
        pairs.append((i, i+1))
        i += 2
    else:
        pairs.append((i, None))
        i += 1

report(f"Identified {len(pairs)} pairs at R={R_values[0]}:")
for p, (i, j) in enumerate(pairs[:6]):
    if j is not None:
        avg = (ev_ref[i] + ev_ref[j]) / 2
        split = ev_ref[j] - ev_ref[i]
        report(f"  Pair {p}: avg={avg:+.6f}, split={split:.2e}")
    else:
        report(f"  Pair {p}: unpaired, ev={ev_ref[i]:+.6f}")

report("")

# Track splittings vs R
report("PAIR SPLITTINGS VS R:")
header2 = f"{'R':>3} {'gap':>4}"
for p in range(min(5, len(pairs))):
    header2 += f" {'split'+str(p):>12}"
report(header2)
report("-" * (7 + 5*13))

for R in R_values:
    ev = all_results[R]
    line = f"{R:3d} {R-2*R_maj:4d}"
    for p, (i, j) in enumerate(pairs[:5]):
        if j is not None and j < len(ev):
            split = ev[j] - ev[i]
            line += f" {split:12.8f}"
        else:
            line += f" {'---':>12}"
    report(line)

report("")

# ============================================================
# BOND ENERGY (eigenvalue shift)
# ============================================================
report("BOND ENERGY — EIGENVALUE SHIFTS")
report("-" * 55)
report("V_n(R) = bonding eigenvalue of pair n - single torus eigenvalue")
report("")

header3 = f"{'R':>3} {'gap':>4}"
for p in range(min(5, len(pairs))):
    header3 += f" {'V'+str(p):>12}"
header3 += f" {'V_total':>12}"
report(header3)
report("-" * (7 + 6*13))

for R in R_values:
    ev = all_results[R]
    line = f"{R:3d} {R-2*R_maj:4d}"
    V_total = 0
    for p, (i, j) in enumerate(pairs[:5]):
        V = ev[i] - evals_s[p] if i < len(ev) and p < len(evals_s) else 0
        V_total += V
        line += f" {V:+12.8f}"
    line += f" {V_total:+12.8f}"
    report(line)

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report("TWO POLOIDAL-WINDING TORI — CORRECT PROTON MODEL")
report(f"  Single torus: {n_neg} tachyons (should be 0)")
report(f"  Lowest single eigenvalue: {evals_s[0]:+.6f}")
report("")

# Any bonding detected?
V_min = 999
R_min = 0
for R in R_values:
    ev = all_results[R]
    V = ev[0] - evals_s[0]
    if V < V_min:
        V_min = V
        R_min = R

if V_min < -1e-6:
    report(f"  BONDING DETECTED: V_min = {V_min:+.8f} at R = {R_min}")
    report(f"  D_e = {-V_min:.8f}")
elif V_min < 0:
    report(f"  Weak attraction: V_min = {V_min:+.8f} at R = {R_min}")
else:
    report(f"  No bonding detected (V_min = {V_min:+.8f})")

# Any tachyons in double-torus?
any_tach = any(np.sum(all_results[R] < -0.01) > 0 for R in R_values)
if any_tach:
    for R in R_values:
        nt = np.sum(all_results[R] < -0.01)
        if nt > 0:
            report(f"  WARNING: {nt} tachyon(s) at R={R}")
else:
    report("  No tachyons at any R — topology holds for bonding too")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
