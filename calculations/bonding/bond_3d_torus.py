"""
3D Two-Torus Bond Simulation
==============================
Places two poloidal-winding tori on a 3D lattice, computes the
Hessian eigenvalues, and extracts the bond energy from the ZPE
difference: D_e = ZPE(two tori at R) - 2*ZPE(single torus).

Uses the CORRECT topology: poloidal winding (0 tachyons).
GPU accelerated (RTX 4070 Ti, CuPy).

This is the full 3D calculation, not 1D approximation.
"""
import sys, io, os, time
import numpy as np

try:
    import cupy as cp
    import cupyx.scipy.sparse as csp
    import cupyx.scipy.sparse.linalg as csla
    HAS_GPU = True
except ImportError:
    HAS_GPU = False
    print("No CuPy — falling back to CPU (slower)")

from scipy import sparse
from scipy.sparse.linalg import eigsh

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3

print("=" * 70)
print("3D TWO-TORUS BOND SIMULATION")
print("=" * 70)
print(f"GPU: {'RTX 4070 Ti (CuPy)' if HAS_GPU else 'CPU fallback'}")
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print()


# ============================================================
# TORUS FIELD BUILDER
# ============================================================

def build_poloidal_torus(X, Y, Z, R_maj, r_tube, center_x=0, center_y=0, center_z=0):
    """
    Build a single torus with poloidal winding.
    Field winds 0 -> 2 around the tube cross-section.

    The torus ring lies in the XY plane, centered at (center_x, center_y, center_z).
    """
    # Shift coordinates to torus center
    Xs = X - center_x
    Ys = Y - center_y
    Zs = Z - center_z

    # Torus coordinates
    rho_xy = np.sqrt(Xs**2 + Ys**2 + 1e-30)  # distance from z-axis
    rho_tube = np.sqrt((rho_xy - R_maj)**2 + Zs**2)  # distance from tube center
    phi_pol = np.arctan2(Zs, rho_xy - R_maj)  # poloidal angle

    # Kink profile: winds 0->2 around the poloidal direction
    K = r_tube  # sharpness = tube radius
    phi_kink = (4.0/PI) * np.arctan(np.exp(K * phi_pol))

    # Envelope: kink exists near tube, decays to vacuum outside
    cutoff = r_tube + 3
    mask_near = rho_tube < cutoff

    phi = np.zeros_like(X)
    phi[mask_near] = phi_kink[mask_near]

    # Far field: snap to nearest vacuum
    mask_far_upper = (~mask_near) & (phi_pol >= 0)
    mask_far_lower = (~mask_near) & (phi_pol < 0)
    phi[mask_far_upper] = 2.0
    phi[mask_far_lower] = 0.0

    return phi


def build_two_tori(N, R_maj, r_tube, gap):
    """
    Build two tori separated by 'gap' along the x-axis.
    Both rings in the XY plane, offset along x.

    Torus A centered at x = -(R_maj + gap/2)
    Torus B centered at x = +(R_maj + gap/2)

    Returns: phi field, coordinate arrays
    """
    ix = np.arange(N, dtype=np.float64) - N/2
    X, Y, Z = np.meshgrid(ix, ix, ix, indexing='ij')

    # Center positions along x-axis
    cx_A = -(R_maj + gap/2)
    cx_B = +(R_maj + gap/2)

    phi_A = build_poloidal_torus(X, Y, Z, R_maj, r_tube, center_x=cx_A)
    phi_B = build_poloidal_torus(X, Y, Z, R_maj, r_tube, center_x=cx_B)

    # Combine: in the gap region, both should be near vacuum (0 or 2)
    # Use the one with the stronger signal (closer to its torus)
    rho_A = np.sqrt((np.sqrt((X - cx_A)**2 + Y**2 + 1e-30) - R_maj)**2 + Z**2)
    rho_B = np.sqrt((np.sqrt((X - cx_B)**2 + Y**2 + 1e-30) - R_maj)**2 + Z**2)

    # Each site belongs to whichever torus is closer
    phi = np.where(rho_A <= rho_B, phi_A, phi_B)

    return phi, X, Y, Z


# ============================================================
# HESSIAN BUILDER (GPU or CPU)
# ============================================================

def build_hessian_3d(phi_flat, N):
    """Build the 3D Hessian: H = -Laplacian + cos(pi*phi), periodic BC."""
    Ntot = N**3
    diag = 2.0 * d + np.cos(PI * phi_flat)

    idx = np.arange(Ntot, dtype=np.int32)
    ix_arr = idx % N
    iy_arr = (idx // N) % N
    iz_arr = idx // (N * N)

    rows_list = [idx]
    cols_list = [idx]
    vals_list = [diag]

    for dx, dy, dz in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
        jx = (ix_arr + dx) % N
        jy = (iy_arr + dy) % N
        jz = (iz_arr + dz) % N
        j = (jz * N * N + jy * N + jx).astype(np.int32)
        rows_list.append(idx)
        cols_list.append(j)
        vals_list.append(-np.ones(Ntot, dtype=np.float64))

    all_rows = np.concatenate(rows_list)
    all_cols = np.concatenate(cols_list)
    all_vals = np.concatenate(vals_list)

    H = sparse.coo_matrix((all_vals, (all_rows, all_cols)), shape=(Ntot, Ntot))
    return H.tocsr()


def get_lowest_eigenvalues(H, k=10):
    """Get k lowest eigenvalues using GPU if available."""
    if HAS_GPU:
        H_csr = H.tocsr()
        H_gpu = csp.csr_matrix(
            (cp.array(H_csr.data), cp.array(H_csr.indices), cp.array(H_csr.indptr)),
            shape=H_csr.shape
        )
        evals = csla.eigsh(H_gpu, k=k, which='SA', return_eigenvectors=False)
        return cp.asnumpy(np.sort(evals))
    else:
        evals = eigsh(H, k=k, which='SM', return_eigenvectors=False)
        return np.sort(evals)


def compute_zpe(evals, mass_gap=0.95):
    """ZPE from bound modes (eigenvalue below mass gap)."""
    bound = evals[evals < mass_gap]
    # For tachyonic modes (should be 0 for correct topology):
    n_tachyon = np.sum(bound < -0.01)
    zpe = np.sum(np.sqrt(np.abs(bound))) / 2
    return zpe, len(bound), n_tachyon


# ============================================================
# MAIN COMPUTATION
# ============================================================

R_maj = 6
r_tube = 3
N = 64  # MUST use same grid for single and double torus

print(f"Torus parameters: R_maj={R_maj}, r_tube={r_tube}")
print(f"Grid: {N}^3 = {N**3:,} sites")
print()

# --- Step 1: Single torus reference ---
print("STEP 1: Single torus reference ZPE")
print("-" * 50)

ix = np.arange(N, dtype=np.float64) - N/2
X, Y, Z = np.meshgrid(ix, ix, ix, indexing='ij')
phi_single = build_poloidal_torus(X, Y, Z, R_maj, r_tube)

n_active = np.sum((phi_single > 0.1) & (phi_single < 1.9))
print(f"  Active sites (0.1 < phi < 1.9): {n_active:,}")
print(f"  Field range: [{phi_single.min():.3f}, {phi_single.max():.3f}]")

t0 = time.time()
H_single = build_hessian_3d(phi_single.ravel(), N)
evals_single = get_lowest_eigenvalues(H_single, k=20)
t1 = time.time()

zpe_single, n_bound_single, n_tach_single = compute_zpe(evals_single)
print(f"  Eigenvalues ({t1-t0:.1f}s):")
for i, ev in enumerate(evals_single[:10]):
    status = "TACHYON!" if ev < -0.01 else ("bound" if ev < 0.95 else "phonon")
    print(f"    [{i}] {ev:+.6f}  {status}")
print(f"  Bound modes: {n_bound_single}, Tachyons: {n_tach_single}")
print(f"  ZPE(single) = {zpe_single:.6f}")
print()

# --- Step 2: Two-torus bond scan ---
print("STEP 2: Two-torus bond energy scan")
print("-" * 50)

# Gap = distance between nearest tube surfaces
# Total center-to-center = 2*R_maj + gap
gaps = [6, 8, 10, 12, 14, 16, 18, 20]

N2 = N  # SAME grid for both single and double (consistent ZPE baseline)
print(f"  Using N={N2} for two-torus calculation (same as single)")
print()

zpe_ref = 2 * zpe_single  # reference: two non-interacting tori
print(f"  ZPE(reference) = 2 * {zpe_single:.6f} = {zpe_ref:.6f}")
print()

print(f"  {'gap':>4} {'ZPE_AB':>10} {'V(gap)':>10} {'n_bound':>7} {'n_tach':>6} {'time':>6}")
print(f"  {'-'*4} {'-'*10} {'-'*10} {'-'*7} {'-'*6} {'-'*6}")

results = []
for gap in gaps:
    # Check if tori fit in the grid
    extent = 2 * (R_maj + r_tube + 3) + gap
    if extent > N2 - 4:
        print(f"  {gap:4d}  SKIP (extent {extent:.0f} > grid {N2})")
        continue

    t0 = time.time()
    phi_two, _, _, _ = build_two_tori(N2, R_maj, r_tube, gap)
    H_two = build_hessian_3d(phi_two.ravel(), N2)
    evals_two = get_lowest_eigenvalues(H_two, k=20)
    t1 = time.time()

    zpe_two, n_bound_two, n_tach_two = compute_zpe(evals_two)
    V = zpe_two - zpe_ref

    results.append((gap, V, n_bound_two, n_tach_two))

    bar = '*' * max(0, int(-V * 500)) if V < 0 else ''
    print(f"  {gap:4d} {zpe_two:10.6f} {V:+10.6f} {n_bound_two:7d} {n_tach_two:6d} {t1-t0:5.1f}s {bar}")

print()

# --- Step 3: Analysis ---
print("STEP 3: Bond energy extraction")
print("-" * 50)

if results:
    gaps_arr = np.array([r[0] for r in results])
    V_arr = np.array([r[1] for r in results])

    # Method 1: raw subtraction from 2*ZPE_single
    i_min = np.argmin(V_arr)
    D_e_raw = -V_arr[i_min]
    gap_eq = gaps_arr[i_min]

    # Method 2: use ASYMPTOTE as reference (more reliable)
    # The largest gap gives V_inf (non-interacting limit on same grid)
    V_inf = V_arr[-1]  # largest gap
    D_e_asym = -(V_arr[i_min] - V_inf)

    print(f"  Equilibrium gap: {gap_eq} sites")
    print(f"  V(minimum):   {V_arr[i_min]:+.6f}")
    print(f"  V(asymptote): {V_inf:+.6f}")
    print(f"  D_e (raw, from 2*ZPE_single): {D_e_raw:.6f}")
    print(f"  D_e (from asymptote):         {D_e_asym:.6f}")
    print()

    # Print V(R) relative to asymptote
    print(f"  Bond curve (relative to asymptote):")
    for gap, V, nb, nt in results:
        D = -(V - V_inf)
        bar = '*' * max(0, int(D * 500))
        print(f"    gap={gap:2d}: D = {D:+.6f} {bar}")
    print()

    # Convert to eV
    s_PT = (-1 + np.sqrt(1 + 8/PI**2)) / 2
    alpha = np.exp(-(2/6) * (2**7/PI**2 + np.log(6)))
    E_H = alpha**2 * 0.51100e6 / 2
    M_kink = 8/PI**2

    D_e = D_e_asym  # use the asymptote method

    D_e_eV_v1 = D_e * E_H / (2*s_PT)  # s-cancellation conversion
    D_e_eV_v2 = D_e / M_kink * E_H     # M_kink conversion
    D_e_eV_v3 = PI/d**2 * E_H           # formula prediction

    print(f"  Conversion to eV:")
    print(f"    D_e * E_H/(2s) = {D_e:.6f} * {E_H/(2*s_PT):.2f} = {D_e_eV_v1:.3f} eV")
    print(f"    D_e/M_kink*E_H = {D_e/M_kink:.6f} * {E_H:.2f} = {D_e_eV_v2:.3f} eV")
    print(f"    Formula pi/d^2*E_H = {D_e_eV_v3:.3f} eV")
    print(f"    Observed H2:   4.478 eV")
    print()
    print(f"  D_e/M_kink = {D_e/M_kink:.6f}")
    print(f"  pi/d^2     = {PI/d**2:.6f}")
    print(f"  Ratio: {D_e/M_kink / (PI/d**2):.4f}")
else:
    print("  No results computed.")

print()
print(f"Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}")
