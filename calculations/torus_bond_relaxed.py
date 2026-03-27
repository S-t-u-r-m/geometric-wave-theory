"""
Two-Torus Bonding V3 — Energy Relaxation
==========================================
The initial two-torus profile is approximate. The Hessian depends on
being at a TRUE minimum, not just a guess.

Fix: start from the nearest-torus profile and RELAX the field to
minimize the lattice energy using gradient descent.

E = Σ_{<ij>} (φ_i - φ_j)²/2 + Σ_i (1/π²)(1 - cos(πφ_i))
∂E/∂φ_i = Σ_{j∈NN} (φ_i - φ_j) + (1/π)sin(πφ_i)

Update: φ_i → φ_i - dt × ∂E/∂φ_i (gradient descent)

After relaxation, compute the Hessian eigenvalues of the TRUE minimum.
GPU for both relaxation and eigsh.
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

outfile = os.path.join(os.path.dirname(__file__), "torus_bond_relaxed_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("TWO-TORUS BONDING V3 — ENERGY RELAXATION")
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
    Xc = X - cx; Yc = Y - cy
    rho_xy = np.sqrt(Xc**2 + Yc**2 + 1e-30)
    rho_tube = np.sqrt((rho_xy - R_maj)**2 + Z**2)
    phi_pol = np.arctan2(Z, rho_xy - R_maj)
    field = (4.0/PI) * np.arctan(np.exp(kink_K * phi_pol))
    cutoff = r_tube + 3
    field[(rho_tube >= cutoff) & (phi_pol >= 0)] = 2.0
    field[(rho_tube >= cutoff) & (phi_pol < 0)] = 0.0
    return field, rho_tube

def compute_gradient_gpu(phi_gpu, N):
    """Compute ∂E/∂φ on GPU using CuPy.
    grad_i = Σ_{j∈NN}(φ_i - φ_j) + (1/π)sin(πφ_i)
    """
    phi = phi_gpu.reshape(N, N, N)

    # Laplacian: sum of (phi_i - phi_j) for 6 neighbors
    lap = 6.0 * phi
    lap -= cp.roll(phi, 1, axis=0) + cp.roll(phi, -1, axis=0)
    lap -= cp.roll(phi, 1, axis=1) + cp.roll(phi, -1, axis=1)
    lap -= cp.roll(phi, 1, axis=2) + cp.roll(phi, -1, axis=2)

    # Potential gradient: (1/π)sin(πφ)
    pot_grad = (1.0/PI) * cp.sin(PI * phi)

    grad = lap + pot_grad
    return grad.ravel()

def compute_energy_gpu(phi_gpu, N):
    """Total lattice energy on GPU."""
    phi = phi_gpu.reshape(N, N, N)

    # Kinetic: Σ (φ_i - φ_j)²/2
    kin = 0.0
    for ax in range(3):
        diff = phi - cp.roll(phi, 1, axis=ax)
        kin += 0.5 * cp.sum(diff**2)

    # Potential: Σ (1/π²)(1 - cos(πφ))
    pot = cp.sum((1.0/PI**2) * (1 - cp.cos(PI * phi)))

    return float(kin + pot)

def relax_field(phi_init, N, max_iter=5000, dt=0.02, tol=1e-6, report_fn=None):
    """Topology-preserving gradient descent.

    Each site is clamped to its initial potential basin:
      Sites starting at φ < 1: stay in [0, 1] (lower vacuum basin)
      Sites starting at φ > 1: stay in [1, 2] (upper vacuum basin)
      Sites near φ = 1 (within 0.05): frozen (they ARE the kink)

    This prevents the gradient descent from unwinding the topology
    while allowing the field to relax within each basin.
    """
    phi_flat = phi_init.ravel().astype(np.float64)

    # Classify each site by its initial basin
    basin_upper = phi_flat > 1.05  # in upper vacuum, clamp to [1, 2]
    basin_lower = phi_flat < 0.95  # in lower vacuum, clamp to [0, 1]
    basin_kink = ~basin_upper & ~basin_lower  # near the kink, freeze

    phi_gpu = cp.asarray(phi_flat)
    basin_kink_gpu = cp.asarray(basin_kink)

    E_prev = compute_energy_gpu(phi_gpu, N)
    if report_fn:
        n_upper = int(cp.sum(cp.asarray(basin_upper)))
        n_lower = int(cp.sum(cp.asarray(basin_lower)))
        n_kink = int(cp.sum(basin_kink_gpu))
        report_fn(f"  Initial energy: {E_prev:.4f}")
        report_fn(f"  Basin counts: upper={n_upper}, lower={n_lower}, kink(frozen)={n_kink}")

    for iteration in range(max_iter):
        grad = compute_gradient_gpu(phi_gpu, N)

        # Zero gradient at frozen kink sites
        grad[basin_kink_gpu] = 0.0

        phi_gpu -= dt * grad

        # Clamp to basins: upper sites stay in [1, 2], lower in [0, 1]
        phi_3d = phi_gpu.reshape(N, N, N)
        phi_gpu = cp.clip(phi_gpu, 0.0, 2.0)
        # Don't need per-basin clip since the frozen kink sites prevent crossing

        if iteration % 200 == 0 or iteration == max_iter - 1:
            E = compute_energy_gpu(phi_gpu, N)
            grad_norm = float(cp.sqrt(cp.sum(grad**2)))
            if report_fn and (iteration % 1000 == 0 or iteration < 3):
                report_fn(f"  iter {iteration:5d}: E={E:.4f}, |grad|={grad_norm:.4f}")

            if abs(E - E_prev) < tol and grad_norm < 0.1:
                if report_fn:
                    report_fn(f"  Converged at iter {iteration}: E={E:.4f}")
                break
            E_prev = E

    phi_out = cp.asnumpy(phi_gpu).reshape(N, N, N)
    return phi_out, E_prev

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
    return np.sort(cp.asnumpy(ev))

# ============================================================
# SINGLE TORUS — RELAX AND GET REFERENCE
# ============================================================
report("SINGLE TORUS — RELAX TO TRUE MINIMUM")
report("-" * 55)

phi_init_s, _ = torus_field_and_distance(X, Y, Z, 0, 0)
phi_relaxed_s, E_s = relax_field(phi_init_s, N, max_iter=3000, report_fn=report)

report(f"Relaxed energy: {E_s:.4f}")
report(f"Field range: [{phi_relaxed_s.min():.4f}, {phi_relaxed_s.max():.4f}]")

evals_s = get_evals(phi_relaxed_s)
n_neg_s = np.sum(evals_s < -0.01)
report(f"Tachyons after relaxation: {n_neg_s}")
report(f"Lowest 8: {' '.join(f'{e:+.6f}' for e in evals_s[:8])}")
report("")

# ============================================================
# TWO TORI — RELAX AT EACH R
# ============================================================
report("TWO TORI — RELAX AND COMPUTE EIGENVALUES")
report("-" * 55)

R_values = [30, 28, 26, 24, 22, 20, 19, 18, 17]

report(f"{'R':>3} {'gap':>4} {'E_relax':>10} {'n_neg':>5} {'ev0':>10} "
       f"{'V0':>10} {'split0':>10}")
report("-" * 60)

bond_data = {}

for R in R_values:
    t0 = time.time()

    # Initial guess: nearest-torus assignment
    phi_A, rho_A = torus_field_and_distance(X, Y, Z, -R/2, 0)
    phi_B, rho_B = torus_field_and_distance(X, Y, Z, +R/2, 0)
    phi_init = np.where(rho_A <= rho_B, phi_A, phi_B)

    # Relax to true minimum
    phi_relax, E_relax = relax_field(phi_init, N, max_iter=3000)

    # Eigenvalues of relaxed configuration
    ev = get_evals(phi_relax)
    gap = R - 2*R_maj
    n_neg = np.sum(ev < -0.01)
    V0 = ev[0] - evals_s[0]
    split0 = ev[1] - ev[0]

    bond_data[R] = {'ev': ev, 'V0': V0, 'n_neg': n_neg, 'E': E_relax}
    elapsed = time.time() - t0

    report(f"{R:3d} {gap:4d} {E_relax:10.2f} {n_neg:5d} {ev[0]:+10.6f} "
           f"{V0:+10.6f} {split0:10.6f}  [{elapsed:.1f}s]")

report("")

# ============================================================
# BOND CURVE
# ============================================================
report("BOND CURVE V(R) — FROM RELAXED CONFIGURATIONS")
report("-" * 55)

for R in R_values:
    V = bond_data[R]['V0']
    marker = ""
    if V == min(bond_data[r]['V0'] for r in R_values):
        marker = " <-- min"
    report(f"  R={R:3d} (gap={R-2*R_maj:2d}): V = {V:+.8f}{marker}")

V_min = min(bond_data[R]['V0'] for R in R_values)
R_min = [R for R in R_values if bond_data[R]['V0'] == V_min][0]

report("")
if V_min < -1e-6:
    report(f"*** BONDING: D_e = {-V_min:.8f} at R = {R_min} ***")
else:
    report("No bonding detected.")

# Tachyon check
any_tach = any(bond_data[R]['n_neg'] > 0 for R in R_values)
if not any_tach:
    report("NO TACHYONS at any R — clean topological bonding!")
else:
    report("Tachyons present at some R:")
    for R in R_values:
        if bond_data[R]['n_neg'] > 0:
            report(f"  R={R}: {bond_data[R]['n_neg']} tachyon(s)")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
