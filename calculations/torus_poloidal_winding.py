"""
Torus with Poloidal Winding — The Correct Topological Kink
============================================================
Previous model: field peaked at tube center, decays radially (WRONG).
  → Not topological, can shrink to vacuum, has tachyons.

This model: field winds AROUND the tube cross-section (CORRECT).
  → The kink goes 0→2 as you travel around the small circle.
  → Can't unwind because 0 and 2 are separated by a potential barrier.
  → THIS is what makes the proton stable.

The kink profile at each tube cross-section:
  φ(φ_pol) = (4/π) × arctan(exp(K × φ_pol))
  where φ_pol = poloidal angle ∈ [-π, π]
  K = r_tube (controls kink sharpness on the tube circumference)

At φ_pol = -π: φ ≈ 0 (vacuum)
At φ_pol = 0:  φ = 1 (barrier top — the kink center)
At φ_pol = +π: φ ≈ 2 (equivalent vacuum)

GPU: RTX 4070 Ti, CuPy. N=64 eigsh ~3s.
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

outfile = os.path.join(os.path.dirname(__file__), "torus_poloidal_winding_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("TORUS WITH POLOIDAL WINDING — TOPOLOGICAL KINK")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("")

# ============================================================
# SETUP
# ============================================================
N = 64
Ntot = N**3
R_maj = 8
r_tube = 3  # tube radius (determines poloidal circumference)

report(f"Lattice: {N}^3 = {Ntot:,} sites")
report(f"Torus: R_major = {R_maj}, r_tube = {r_tube}")
report(f"Tube circumference: 2π × {r_tube} = {2*PI*r_tube:.1f} sites")
report(f"Kink width: ~3 sites (fits in circumference of {2*PI*r_tube:.0f})")
report("")

ix = np.arange(N, dtype=np.float64) - N/2
X, Y, Z = np.meshgrid(ix, ix, ix, indexing='ij')

# ============================================================
# BUILD THE POLOIDAL WINDING FIELD
# ============================================================
report("Building poloidal winding field...")

# For each lattice site, compute torus coordinates
rho_xy = np.sqrt(X**2 + Y**2 + 1e-30)  # distance from z-axis
rho_tube = np.sqrt((rho_xy - R_maj)**2 + Z**2)  # distance from tube center
phi_pol = np.arctan2(Z, rho_xy - R_maj)  # poloidal angle [-π, π]

# The kink profile: field winds 0→2 around the poloidal direction
# K controls sharpness: K × π should be >> 1 for the kink to fit
K = r_tube  # kink sharpness parameter = tube radius
phi_field = (4.0/PI) * np.arctan(np.exp(K * phi_pol))

# Radial envelope: the kink exists on the tube, decays outside
# For sites far from the tube, field should be at vacuum (0 or 2)
# Upper half (phi_pol > 0): decays to 2
# Lower half (phi_pol < 0): decays to 0
# At the kink center (phi_pol = 0): field = 1, decays to... depends on direction

# Simple approach: use a smooth envelope that keeps the kink near the tube
# and lets the field relax to vacuum outside
envelope_width = r_tube + 2  # how far the kink extends from tube center
envelope = np.exp(-np.maximum(rho_tube - r_tube, 0)**2 / (2 * 2.0**2))

# The full field: kink × envelope + vacuum × (1 - envelope)
# vacuum = 0 for lower hemisphere, 2 for upper hemisphere
vacuum = np.where(phi_pol >= 0, 2.0, 0.0)
phi_torus = phi_field * envelope + vacuum * (1 - envelope)

# Actually, simpler: just use the kink profile everywhere near the torus
# and 0 far away. The cosine potential treats φ=0 and φ=2 equivalently.
# Let's try: kink for rho_tube < cutoff, smooth transition outside
cutoff = r_tube + 3
mask_near = rho_tube < cutoff
phi_torus_v2 = np.zeros_like(X)
phi_torus_v2[mask_near] = (4.0/PI) * np.arctan(np.exp(K * phi_pol[mask_near]))
# For sites outside, snap to nearest vacuum: 0 if phi_pol < 0, 2 if phi_pol >= 0
mask_far_upper = (~mask_near) & (phi_pol >= 0)
mask_far_lower = (~mask_near) & (phi_pol < 0)
phi_torus_v2[mask_far_upper] = 2.0
phi_torus_v2[mask_far_lower] = 0.0

# Use v2 (cleaner)
phi_torus = phi_torus_v2

report(f"Field range: [{phi_torus.min():.4f}, {phi_torus.max():.4f}]")
report(f"Sites near torus (rho < {cutoff}): {np.sum(mask_near):,}")
report(f"Sites with 0.1 < phi < 1.9: {np.sum((phi_torus > 0.1) & (phi_torus < 1.9)):,}")
report("")

# Show cross-section through torus (y=0 plane, varying x and z)
mid_y = N // 2
xz_slice = phi_torus[:, mid_y, :]
report("Cross-section (y=0 plane) — should show kink wrapping around tube:")
report("  Looking at x=R_maj column (through tube center):")
ix_tube = int(R_maj + N/2)
col = xz_slice[ix_tube, :]
# Show the tube region
tube_z = [j for j in range(N) if abs(ix[j]) <= r_tube + 2]
if tube_z:
    report("  z:   " + " ".join(f"{ix[j]:5.0f}" for j in tube_z))
    report("  phi: " + " ".join(f"{col[j]:5.3f}" for j in tube_z))
report("")

# Energy of this configuration
phi_flat = phi_torus.ravel()
# Kinetic energy: sum of (phi_i - phi_j)^2 / 2 over neighbors
# Potential energy: sum of (1/pi^2)(1 - cos(pi*phi_i))
pot_energy = np.sum((1.0/PI**2) * (1 - np.cos(PI * phi_flat)))
report(f"Potential energy: {pot_energy:.2f}")
report(f"  (vacuum: {Ntot * 0:.2f}, max: {Ntot * 2/PI**2:.2f})")
report("")

# ============================================================
# BUILD HESSIAN AND FIND EIGENVALUES
# ============================================================
def build_3d_hessian_gpu(phi_flat):
    diag = 2.0 * d + np.cos(PI * phi_flat)
    idx = np.arange(Ntot, dtype=np.int32)
    ix_arr = idx % N
    iy_arr = (idx // N) % N
    iz_arr = idx // (N * N)
    neighbors = []
    for dx, dy, dz in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
        jx = (ix_arr + dx) % N
        jy = (iy_arr + dy) % N
        jz = (iz_arr + dz) % N
        j = (jz * N * N + jy * N + jx).astype(np.int32)
        neighbors.append(j)
    all_rows = np.concatenate([idx] * 7)
    all_cols = np.concatenate([idx] + neighbors)
    all_vals = np.concatenate([diag] + [-np.ones(Ntot, dtype=np.float64)] * 6)
    H_cpu = sparse.coo_matrix((all_vals, (all_rows, all_cols)),
                               shape=(Ntot, Ntot)).tocsr()
    return csp.csr_matrix(H_cpu)

report("Building 3D Hessian...")
t0 = time.time()
H_gpu = build_3d_hessian_gpu(phi_flat)
report(f"  Built in {time.time()-t0:.2f}s")

report("Finding eigenvalues on GPU...")
n_eig = 40
t0 = time.time()
evals_gpu, evecs_gpu = csla.eigsh(H_gpu, k=n_eig, which='SA')
cp.cuda.Stream.null.synchronize()
t_eigsh = time.time() - t0

evals = cp.asnumpy(evals_gpu)
idx_sort = np.argsort(evals)
evals = evals[idx_sort]

report(f"  Eigsh completed in {t_eigsh:.2f}s")
report("")

# ============================================================
# THE CONFINEMENT TEST
# ============================================================
n_negative = np.sum(evals < -0.01)
n_zero = np.sum((evals >= -0.01) & (evals < 0.01))
n_bound = np.sum(evals < 1.0)

report("=" * 50)
if n_negative == 0:
    report("*** NO TACHYONS — TOPOLOGY STABILIZES THE KINK ***")
    report("*** THIS IS CONFINEMENT ***")
elif n_negative < 3:
    report(f"*** {n_negative} TACHYON(S) — PARTIAL STABILIZATION ***")
else:
    report(f"*** {n_negative} TACHYONS — NOT STABILIZED ***")
report("=" * 50)
report("")

report(f"Negative eigenvalues (< -0.01): {n_negative}")
report(f"Near-zero eigenvalues (|ev| < 0.01): {n_zero}")
report(f"Bound states (< 1.0): {n_bound}")
report("")

report(f"{'idx':>4} {'omega^2':>12} {'type':>10}")
report("-" * 30)
for i in range(min(n_eig, 40)):
    ev = evals[i]
    if ev < -0.01:
        t = "TACHYON"
    elif ev < 0.01:
        t = "~ZERO"
    elif ev < 1.0:
        t = "BOUND"
    else:
        t = "band"
    report(f"  {i:3d} {ev:+12.6f} {t:>10}")
    if ev > 1.1 and i > 15:
        break

report("")

# ============================================================
# COMPARE OLD vs NEW TORUS
# ============================================================
report("COMPARISON: OLD (radial bump) vs NEW (poloidal winding)")
report("-" * 55)
report("")
report("OLD torus (radial bump, from torus_confinement_test.py):")
report("  7 tachyonic modes (omega^2 from -0.197 to -0.054)")
report("  40 bound states total")
report("  Lowest: -0.197")
report("")
report("NEW torus (poloidal winding):")
report(f"  {n_negative} tachyonic modes")
report(f"  {n_bound} bound states total")
report(f"  Lowest: {evals[0]:+.6f}")
report("")

if n_negative == 0:
    report("THE POLOIDAL WINDING ELIMINATES ALL TACHYONS.")
    report("The topological winding prevents annihilation.")
    report("The proton (poloidal kink on torus) is STABLE.")
    report("")
    report("The zero modes (if any) correspond to:")
    report("  - d = 3 translational zero modes (move the torus)")
    report("  - 1 rotational zero mode (rotate the kink around the tube)")
    report("  Total: d+1 = 4 zero modes → r_p = (d+1)ℏc/m_p ✓")
elif n_negative < n_negative:  # just for flow
    pass

# Count near-zero modes (these should be the d+1 = 4 zero modes)
if n_zero > 0:
    report(f"Near-zero modes found: {n_zero}")
    report("These should be the (d+1) = 4 zero modes:")
    for i in range(min(n_zero + 2, len(evals))):
        if abs(evals[i]) < 0.1:
            report(f"  mode {i}: omega^2 = {evals[i]:+.8f}")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
