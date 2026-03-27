"""
Sigma vs Pi Bond from 2D Hessian — The W_PI Test
===================================================
The σ/π bond strength ratio should be W_PI = cos(π/d) = 0.5.

In 1D, there's no transverse direction, so we can't distinguish σ from π.
In 2D, we can: σ has no node perpendicular to the bond, π has one node.

Setup:
  - 2D lattice (Nx × Nz), bond along x, transverse along z
  - Two torus cross-sections at x = ±R/2
  - Each cross-section: the REAL profile from the 3D torus, sliced in xz
  - The 2D Hessian: -Laplacian_2d + cos(π·φ)
  - Eigenvalues: σ bond = lowest (symmetric in z), π bond = next (antisymmetric)
  - The ratio of bonding-antibonding splittings: split_π / split_σ = W_PI?

The torus cross-section in the xz-plane (at y=0) shows the tube as two
bumps at z = ±0 (the tube top and bottom). When two tori approach in x,
the σ mode connects the bumps along x (head-on), while the π mode
connects them with a z-node (side-by-side).
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "sigma_pi_2d_hessian_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("SIGMA vs PI BOND FROM 2D HESSIAN")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"d = {d}, gamma = {gamma:.10f}")
report(f"V8 W_PI = cos(π/d) = {np.cos(PI/d):.6f}")
report("")

# ============================================================
# BUILD THE 3D TORUS AND EXTRACT 2D CROSS-SECTION
# ============================================================
N3 = 64
R_maj = 8
kink_width = 3

ix = np.arange(N3, dtype=np.float64) - N3/2
X3, Y3, Z3 = np.meshgrid(ix, ix, ix, indexing='ij')

def make_torus(X, Y, Z, center_x, center_y, z_offset, R_major, kw):
    """Torus with ring in the xy-plane, centered at (center_x, center_y, z_offset)."""
    rho_xy = np.sqrt((X - center_x)**2 + (Y - center_y)**2 + 1e-30)
    rho_tube = np.sqrt((rho_xy - R_major)**2 + (Z - z_offset)**2)
    return (4.0/PI) * (np.arctan(np.exp(rho_tube + kw/2.0))
                      - np.arctan(np.exp(rho_tube - kw/2.0)))

# Single torus at origin, ring in xy-plane
phi_torus = make_torus(X3, Y3, Z3, 0, 0, 0, R_maj, kink_width)

# Extract the xz-plane cross-section (y=0)
mid_y = N3 // 2
phi_xz = phi_torus[:, mid_y, :]  # shape (N3, N3), indexed by (ix, iz)

report("Torus cross-section in xz-plane (y=0):")
report(f"  Shape: {phi_xz.shape}")
report(f"  Peak: {phi_xz.max():.4f}")
report(f"  The tube appears as two bumps at x ≈ ±R_maj = ±{R_maj}")
report("")

# Show the cross-section
mid_x = N3 // 2
mid_z = N3 // 2
report("Cross-section profile (z=0 row, should show two bumps):")
row = phi_xz[:, mid_z]
sig = np.abs(row) > 0.01
if np.any(sig):
    i_lo, i_hi = np.where(sig)[0][0], np.where(sig)[0][-1]
    report("  x:   " + " ".join(f"{ix[i]:5.0f}" for i in range(max(0,i_lo-1), min(N3,i_hi+2))))
    report("  phi: " + " ".join(f"{row[i]:5.3f}" for i in range(max(0,i_lo-1), min(N3,i_hi+2))))
report("")

# Cross-section at x=R_maj (through one bump, varying z):
ix_bump = int(R_maj + N3/2)
col = phi_xz[ix_bump, :]
report(f"Cross-section through tube (x=+{R_maj}, varying z):")
sig_z = np.abs(col) > 0.01
if np.any(sig_z):
    j_lo, j_hi = np.where(sig_z)[0][0], np.where(sig_z)[0][-1]
    report("  z:   " + " ".join(f"{ix[j]:5.0f}" for j in range(max(0,j_lo-1), min(N3,j_hi+2))))
    report("  phi: " + " ".join(f"{col[j]:5.3f}" for j in range(max(0,j_lo-1), min(N3,j_hi+2))))
report("")

# ============================================================
# 2D HESSIAN CONSTRUCTION
# ============================================================
Nx = 256
Nz = 64
Ntot = Nx * Nz

report(f"2D lattice: {Nx} × {Nz} = {Ntot} sites")
report("")

def build_2d_hessian(phi_2d, Nx, Nz):
    """Build sparse 2D Hessian: -Laplacian_2d + cos(π·φ).
    phi_2d: shape (Nx, Nz), flattened to 1D with index = ix * Nz + iz.
    """
    Ntot = Nx * Nz
    phi_flat = phi_2d.ravel()

    # Diagonal: 4 (from 2D Laplacian) + cos(π·φ)
    diag = 4.0 + np.cos(PI * phi_flat)

    # Build COO
    idx = np.arange(Ntot, dtype=np.int32)
    ix_arr = idx // Nz
    iz_arr = idx % Nz

    rows_list = [idx]
    cols_list = [idx]
    vals_list = [diag]

    # x neighbors (periodic)
    for dx in [1, -1]:
        jx = (ix_arr + dx) % Nx
        j = jx * Nz + iz_arr
        rows_list.append(idx)
        cols_list.append(j)
        vals_list.append(-np.ones(Ntot))

    # z neighbors (periodic)
    for dz in [1, -1]:
        jz = (iz_arr + dz) % Nz
        j = ix_arr * Nz + jz
        rows_list.append(idx)
        cols_list.append(j)
        vals_list.append(-np.ones(Ntot))

    all_rows = np.concatenate(rows_list)
    all_cols = np.concatenate(cols_list)
    all_vals = np.concatenate(vals_list)

    H = sparse.coo_matrix((all_vals, (all_rows, all_cols)),
                          shape=(Ntot, Ntot)).tocsr()
    return H

def make_torus_xz_2d(x_arr, z_arr, center_x, R_major, kw):
    """2D torus cross-section in the xz-plane.
    The torus ring is in the xy-plane, so at y=0 we see the tube
    at x = center_x ± R_major, extending in z.
    """
    X2, Z2 = np.meshgrid(x_arr, z_arr, indexing='ij')

    # Distance from ring (at radius R_major from center_x in the x-direction)
    # In the xz-plane at y=0: the ring intersects at x = center_x ± R_major
    # The tube radius = kink profile

    # Distance from the nearest point on the ring:
    # The ring is a circle of radius R_major in the xy-plane.
    # At y=0, the two closest points are at (center_x ± R_major, 0, 0).
    # For each point (x, 0, z), the distance to the ring is:
    rho_from_center = np.abs(X2 - center_x)  # distance from torus center in x
    rho_tube = np.sqrt((rho_from_center - R_major)**2 + Z2**2)

    phi = (4.0/PI) * (np.arctan(np.exp(rho_tube + kw/2.0))
                      - np.arctan(np.exp(rho_tube - kw/2.0)))
    return phi

# ============================================================
# PART 1: SINGLE TORUS — REFERENCE EIGENSPECTRUM
# ============================================================
report("PART 1: SINGLE TORUS CROSS-SECTION — 2D EIGENSPECTRUM")
report("-" * 55)

x_2d = np.arange(Nx, dtype=np.float64) - Nx/2
z_2d = np.arange(Nz, dtype=np.float64) - Nz/2

phi_single = make_torus_xz_2d(x_2d, z_2d, 0, R_maj, kink_width)

report(f"Single torus xz cross-section: peak = {phi_single.max():.4f}")

H_single = build_2d_hessian(phi_single, Nx, Nz)
n_eig = 12
evals_s, evecs_s = eigsh(H_single, k=n_eig, which='SM')
idx_sort = np.argsort(evals_s)
evals_s = evals_s[idx_sort]
evecs_s = evecs_s[:, idx_sort]

report(f"{'mode':>4} {'omega^2':>12} {'status':>8}")
for i in range(n_eig):
    status = "BOUND" if evals_s[i] < 1.0 else "band"
    report(f"  {i:4d} {evals_s[i]:12.6f} {status:>8}")

E0_ref = evals_s[0]
report(f"\nReference: E0 = {E0_ref:.8f}")
report("")

# Classify modes by z-symmetry
report("MODE SYMMETRY IN z (σ vs π character):")
for i in range(min(8, n_eig)):
    evec = evecs_s[:, i].reshape(Nx, Nz)
    # Check z-symmetry: is the mode symmetric (σ) or antisymmetric (π) in z?
    mid_z_2d = Nz // 2
    # Compare evec(x, z) with evec(x, -z)
    evec_plus = evec[:, mid_z_2d:]  # z >= 0
    evec_minus = evec[:, mid_z_2d::-1][:, :evec_plus.shape[1]]  # z <= 0, reversed
    if evec_plus.shape == evec_minus.shape:
        sym = np.sum(evec_plus * evec_minus) / (np.sum(evec_plus**2) + 1e-30)
        z_type = "σ-type" if sym > 0.5 else "π-type" if sym < -0.5 else "mixed"
    else:
        z_type = "?"

    # Also check: where is the mode localized?
    # Is it at the bump positions (x ≈ ±R_maj)?
    evec_sq = evec**2
    x_mean = np.sum(x_2d[:, None] * evec_sq) / np.sum(evec_sq)

    report(f"  Mode {i}: ω² = {evals_s[i]:.6f}, z-sym = {sym:+.3f} ({z_type}), "
           f"<x> = {x_mean:.1f}")

report("")

# ============================================================
# PART 2: TWO TORUS CROSS-SECTIONS — SCAN R
# ============================================================
report("PART 2: TWO TORI — σ AND π SPLITTING VS R")
report("-" * 55)
report("Two torus cross-sections separated by R along x-axis.")
report("Track σ (z-symmetric) and π (z-antisymmetric) splittings.")
report("")

# The two tori are placed at x = ±R_sep/2
# Their rings are at distance R_sep apart, so the tube surfaces
# are at distance R_sep - 2*R_maj apart (when R_sep > 2*R_maj)
# For R_sep < 2*R_maj, the rings overlap.

# We need R_sep such that the tubes are close:
# Tube-to-tube distance ≈ R_sep - 2*R_maj (for the inner edges)
# For bonding: want tube-to-tube ≈ 2-10 sites
# So R_sep ≈ 2*R_maj + 2 to 2*R_maj + 10 = 18 to 26

R_values = list(range(14, 36, 2))

report(f"{'R':>4} {'E0_bond':>12} {'E0_anti':>12} {'E1_bond':>12} {'E1_anti':>12} "
       f"{'σ_split':>10} {'π_split':>10} {'π/σ':>8}")
report("-" * 85)

results = {}

for R_sep in R_values:
    # Two torus cross-sections
    phi_A = make_torus_xz_2d(x_2d, z_2d, -R_sep/2, R_maj, kink_width)
    phi_B = make_torus_xz_2d(x_2d, z_2d, +R_sep/2, R_maj, kink_width)
    phi_double = phi_A + phi_B

    H_double = build_2d_hessian(phi_double, Nx, Nz)
    evals_d, evecs_d = eigsh(H_double, k=12, which='SM')
    idx_sort = np.argsort(evals_d)
    evals_d = evals_d[idx_sort]
    evecs_d = evecs_d[:, idx_sort]

    # Classify each mode by z-symmetry
    sigma_evals = []
    pi_evals = []

    for i in range(min(8, len(evals_d))):
        evec = evecs_d[:, i].reshape(Nx, Nz)
        mid_z_2d = Nz // 2
        evec_plus = evec[:, mid_z_2d:]
        evec_minus = evec[:, mid_z_2d::-1][:, :evec_plus.shape[1]]
        if evec_plus.shape == evec_minus.shape:
            sym = np.sum(evec_plus * evec_minus) / (np.sum(evec_plus**2) + 1e-30)
        else:
            sym = 0
        if sym > 0.3:
            sigma_evals.append(evals_d[i])
        elif sym < -0.3:
            pi_evals.append(evals_d[i])

    # σ splitting = difference between first two σ-type modes
    sigma_split = sigma_evals[1] - sigma_evals[0] if len(sigma_evals) >= 2 else 0

    # π splitting = difference between first two π-type modes
    pi_split = pi_evals[1] - pi_evals[0] if len(pi_evals) >= 2 else 0

    ratio = pi_split / sigma_split if sigma_split > 1e-12 else 0

    results[R_sep] = {
        'sigma_split': sigma_split, 'pi_split': pi_split, 'ratio': ratio,
        'sigma_evals': sigma_evals[:4], 'pi_evals': pi_evals[:4],
        'all_evals': evals_d[:8].tolist()
    }

    e0b = sigma_evals[0] if sigma_evals else 0
    e0a = sigma_evals[1] if len(sigma_evals) > 1 else 0
    e1b = pi_evals[0] if pi_evals else 0
    e1a = pi_evals[1] if len(pi_evals) > 1 else 0

    report(f"{R_sep:4d} {e0b:12.6f} {e0a:12.6f} {e1b:12.6f} {e1a:12.6f} "
           f"{sigma_split:10.6f} {pi_split:10.6f} {ratio:8.4f}")

report("")

# ============================================================
# PART 3: THE KEY RATIO — π/σ SPLITTING
# ============================================================
report("PART 3: THE π/σ SPLITTING RATIO")
report("-" * 55)
report(f"V8 prediction: W_PI = cos(π/d) = {np.cos(PI/d):.6f}")
report("")

R_arr = np.array(sorted(results.keys()))
ratio_arr = np.array([results[R]['ratio'] for R in R_arr])
sigma_arr = np.array([results[R]['sigma_split'] for R in R_arr])
pi_arr_data = np.array([results[R]['pi_split'] for R in R_arr])

# Where both splittings are significant
mask = (sigma_arr > 1e-6) & (pi_arr_data > 1e-6)
if np.sum(mask) > 0:
    avg_ratio = np.mean(ratio_arr[mask])
    report(f"Average π/σ ratio (where both > 1e-6): {avg_ratio:.6f}")
    report(f"V8 W_PI = cos(π/d) = {np.cos(PI/d):.6f}")
    report(f"Ratio of ratios: {avg_ratio / np.cos(PI/d):.4f}")
    report("")

    # At each R
    report(f"{'R':>4} {'σ_split':>12} {'π_split':>12} {'π/σ':>10} {'vs W_PI':>10}")
    report("-" * 50)
    for R in R_arr:
        if results[R]['sigma_split'] > 1e-8 and results[R]['pi_split'] > 1e-8:
            r = results[R]['ratio']
            vs = r / np.cos(PI/d)
            report(f"{R:4.0f} {results[R]['sigma_split']:12.8f} "
                   f"{results[R]['pi_split']:12.8f} {r:10.6f} {vs:10.4f}")

report("")

# ============================================================
# PART 4: ALTERNATIVE — USE EIGENVALUE SHIFTS INSTEAD
# ============================================================
report("PART 4: BONDING SHIFTS (eigenvalue - reference)")
report("-" * 55)
report("The bonding SHIFT (not splitting) for σ and π modes.")
report("")

report(f"{'R':>4} {'σ_shift':>12} {'π_shift':>12} {'π/σ shift':>12}")
report("-" * 42)

for R in R_arr:
    sigma_ev = results[R]['sigma_evals']
    pi_ev = results[R]['pi_evals']

    sigma_shift = sigma_ev[0] - E0_ref if sigma_ev else 0
    pi_shift = pi_ev[0] - E0_ref if pi_ev else 0

    ratio_shift = pi_shift / sigma_shift if abs(sigma_shift) > 1e-12 else 0

    report(f"{R:4.0f} {sigma_shift:+12.8f} {pi_shift:+12.8f} {ratio_shift:+12.6f}")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report(f"V8 W_PI = cos(π/d) = {np.cos(PI/d):.6f}")
report("")

if np.sum(mask) > 0:
    report(f"Measured π/σ splitting ratio: {avg_ratio:.6f}")
    report(f"Agreement with W_PI: {abs(avg_ratio - np.cos(PI/d))/np.cos(PI/d)*100:.1f}%")
else:
    report("Could not measure π/σ ratio (insufficient data)")

report("")
report("If the ratio matches, the W_PI = cos(π/d) weight in the V8 bond")
report("model comes from the TRANSVERSE MODE STRUCTURE of the torus — ")
report("the ratio of z-antisymmetric (π) to z-symmetric (σ) tunnel coupling.")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
