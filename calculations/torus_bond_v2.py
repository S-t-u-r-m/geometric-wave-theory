"""
Two-Torus Bonding V2 — Nearest-Torus Field Assignment
=======================================================
V1 failed because φ_A + φ_B double-counts in the overlap region.
Fix: assign each lattice site to its NEAREST torus.

For well-separated tori: each site is clearly closer to one torus.
At the boundary: smooth blending using distance-weighted average.
This creates a "domain wall" between the two tori = the bond.

Also try: modular field (φ mod 2) to handle the periodicity.
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

outfile = os.path.join(os.path.dirname(__file__), "torus_bond_v2_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("TWO-TORUS BONDING V2 — NEAREST-TORUS ASSIGNMENT")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("")

N = 64
Ntot = N**3
R_maj = 8
r_tube = 3
kink_K = r_tube

ix = np.arange(N, dtype=np.float64) - N/2
X, Y, Z = np.meshgrid(ix, ix, ix, indexing='ij')

def torus_field_and_distance(X, Y, Z, cx, cy):
    """Compute poloidal winding field AND distance to tube center."""
    Xc = X - cx
    Yc = Y - cy
    rho_xy = np.sqrt(Xc**2 + Yc**2 + 1e-30)
    rho_tube = np.sqrt((rho_xy - R_maj)**2 + Z**2)
    phi_pol = np.arctan2(Z, rho_xy - R_maj)

    field = (4.0/PI) * np.arctan(np.exp(kink_K * phi_pol))

    cutoff = r_tube + 3
    mask_far_upper = (rho_tube >= cutoff) & (phi_pol >= 0)
    mask_far_lower = (rho_tube >= cutoff) & (phi_pol < 0)
    field[mask_far_upper] = 2.0
    field[mask_far_lower] = 0.0

    return field, rho_tube

def build_two_torus_field(R):
    """Build field for two tori using nearest-torus assignment."""
    phi_A, rho_A = torus_field_and_distance(X, Y, Z, -R/2, 0)
    phi_B, rho_B = torus_field_and_distance(X, Y, Z, +R/2, 0)

    # Method 1: Nearest torus (sharp boundary)
    use_A = rho_A <= rho_B
    phi_nearest = np.where(use_A, phi_A, phi_B)

    # Method 2: Distance-weighted blend (smooth boundary)
    w_A = 1.0 / (rho_A + 0.1)
    w_B = 1.0 / (rho_B + 0.1)
    phi_blend = (w_A * phi_A + w_B * phi_B) / (w_A + w_B)

    # Method 3: Modular (take sum mod 2, preserving winding)
    phi_sum = phi_A + phi_B
    phi_mod = phi_sum - 2.0 * np.floor(phi_sum / 2.0)

    return phi_nearest, phi_blend, phi_mod, phi_A, phi_B

def build_hessian_gpu(phi_flat):
    diag = 2.0 * d + np.cos(PI * phi_flat)
    idx = np.arange(Ntot, dtype=np.int32)
    ix_a = idx % N; iy_a = (idx // N) % N; iz_a = idx // (N * N)
    neighbors = []
    for dx, dy, dz in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
        jx = (ix_a + dx) % N; jy = (iy_a + dy) % N; jz = (iz_a + dz) % N
        neighbors.append((jz * N * N + jy * N + jx).astype(np.int32))
    all_rows = np.concatenate([idx] * 7)
    all_cols = np.concatenate([idx] + neighbors)
    all_vals = np.concatenate([diag] + [-np.ones(Ntot, dtype=np.float64)] * 6)
    H_cpu = sparse.coo_matrix((all_vals, (all_rows, all_cols)),
                               shape=(Ntot, Ntot)).tocsr()
    return csp.csr_matrix(H_cpu)

def get_evals(phi_3d, n_eig=20):
    H = build_hessian_gpu(phi_3d.ravel())
    ev, _ = csla.eigsh(H, k=n_eig, which='SA')
    cp.cuda.Stream.null.synchronize()
    ev = cp.asnumpy(ev)
    return np.sort(ev)

# ============================================================
# SINGLE TORUS REFERENCE
# ============================================================
report("SINGLE TORUS REFERENCE")
report("-" * 40)
phi_s, _ = torus_field_and_distance(X, Y, Z, 0, 0)
evals_s = get_evals(phi_s)
report(f"Lowest 6: {' '.join(f'{e:+.6f}' for e in evals_s[:6])}")
report(f"Tachyons: {np.sum(evals_s < -0.01)}")
report("")

# ============================================================
# TEST ALL THREE METHODS AT ONE R
# ============================================================
R_test = 24
report(f"TESTING 3 METHODS AT R={R_test} (gap={R_test-2*R_maj})")
report("-" * 55)

phi_near, phi_blend, phi_mod, _, _ = build_two_torus_field(R_test)

for name, phi in [("nearest", phi_near), ("blend", phi_blend), ("modular", phi_mod)]:
    ev = get_evals(phi)
    n_neg = np.sum(ev < -0.01)
    V0 = ev[0] - evals_s[0]
    report(f"  {name:>8}: lowest={ev[0]:+.6f}, tachyons={n_neg}, "
           f"V0={V0:+.6f}, range=[{phi.min():.2f},{phi.max():.2f}]")

report("")

# Pick the best method (fewest tachyons, cleanest bonding)
# Use nearest-torus for the full scan
report("Using NEAREST-TORUS method for full scan.")
report("")

# ============================================================
# FULL R SCAN — NEAREST TORUS
# ============================================================
report("FULL R SCAN — NEAREST TORUS")
report("-" * 55)

R_values = list(range(17, 32))

report(f"{'R':>3} {'gap':>4} {'n_neg':>5} {'ev0':>10} {'ev1':>10} "
       f"{'V0':>10} {'split0':>10}")
report("-" * 60)

bond_data = {}

for R in R_values:
    phi_near, _, _, _, _ = build_two_torus_field(R)
    ev = get_evals(phi_near)

    gap = R - 2*R_maj
    n_neg = np.sum(ev < -0.01)
    V0 = ev[0] - evals_s[0]
    split0 = ev[1] - ev[0] if len(ev) > 1 else 0

    bond_data[R] = {'ev': ev, 'V0': V0, 'n_neg': n_neg}

    report(f"{R:3d} {gap:4d} {n_neg:5d} {ev[0]:+10.6f} {ev[1]:+10.6f} "
           f"{V0:+10.6f} {split0:10.6f}")

report("")

# ============================================================
# BOND CURVE
# ============================================================
report("BOND CURVE V(R)")
report("-" * 55)

R_arr = np.array(R_values, dtype=float)
V_arr = np.array([bond_data[R]['V0'] for R in R_values])

i_min = np.argmin(V_arr)

report(f"V(R) for lowest mode:")
for R in R_values:
    V = bond_data[R]['V0']
    marker = " <-- min" if R == R_arr[i_min] else ""
    report(f"  R={R:3d} (gap={R-2*R_maj:2d}): V = {V:+.8f}{marker}")

report("")

if V_arr[i_min] < -1e-6:
    report(f"*** MORSE WELL: D_e = {-V_arr[i_min]:.8f} at R = {R_arr[i_min]:.0f} ***")

    # Decay rate
    mask = (V_arr < -1e-8) & (R_arr > R_arr[i_min])
    if np.sum(mask) > 2:
        coeffs = np.polyfit(R_arr[mask], np.log(-V_arr[mask]), 1)
        report(f"  Morse decay rate: {-coeffs[0]:.4f}")
else:
    report("No bonding detected.")

# Any tachyons?
any_tach = any(bond_data[R]['n_neg'] > 0 for R in R_values)
if any_tach:
    report("WARNING: tachyons at some R values")
    for R in R_values:
        if bond_data[R]['n_neg'] > 0:
            report(f"  R={R}: {bond_data[R]['n_neg']} tachyon(s)")
else:
    report("No tachyons at any R — clean bonding!")

# ============================================================
# MULTI-MODE BOND SHIFTS
# ============================================================
report("")
report("MULTI-MODE BOND SHIFTS V_n(R)")
report("-" * 55)
n_pairs = 5

header = f"{'R':>3}"
for p in range(n_pairs):
    header += f" {'V'+str(p):>10}"
header += f" {'V_total':>10}"
report(header)
report("-" * (3 + (n_pairs+1)*11))

for R in R_values:
    ev = bond_data[R]['ev']
    line = f"{R:3d}"
    V_total = 0
    for p in range(n_pairs):
        V = ev[p] - evals_s[min(p, len(evals_s)-1)]
        V_total += V
        line += f" {V:+10.6f}"
    line += f" {V_total:+10.6f}"
    report(line)

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
