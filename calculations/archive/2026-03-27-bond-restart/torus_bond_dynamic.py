"""
Two-Torus Bonding — Damped Time Evolution (The Real Physics)
==============================================================
The sine-Gordon wave equation AUTOMATICALLY preserves topology.
Unlike gradient descent, the dynamics can't unwind a kink.

Method:
  1. Start from nearest-torus initial condition (has both windings)
  2. Evolve: φ̈ = Σ(φ_j - φ_i) - (1/π)sin(πφ)
  3. Add damping: φ̇ *= (1 - γ) each step (cool to T=0)
  4. Wait until kinetic energy → 0
  5. The settled field = true two-torus minimum
  6. Compute Hessian eigenvalues of the equilibrium

GPU: all operations are array rolls and element-wise — maximally parallel.
"""
import sys, io, os, time
import numpy as np
import cupy as cp
import cupyx.scipy.sparse as csp
import cupyx.scipy.sparse.linalg as csla
from scipy import sparse

try:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
except:
    pass
sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout, 'reconfigure') else None

PI = np.pi
d = 3

outfile = os.path.join(os.path.dirname(__file__), "torus_bond_dynamic_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("TWO-TORUS BONDING — DAMPED TIME EVOLUTION")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("")

N = 64
Ntot = N**3
R_maj = 8
r_tube = 3
kink_K = r_tube
dt = 0.15       # time step (larger for faster convergence, CFL limit ~0.577)
damping = 0.01   # stronger damping for faster settling

report(f"Lattice: {N}^3 = {Ntot:,} sites")
report(f"Torus: R_major = {R_maj}, r_tube = {r_tube}")
report(f"dt = {dt}, damping = {damping}")
report("")

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

def compute_acceleration_gpu(phi):
    """φ̈ = Laplacian(φ) - (1/π)sin(πφ) on 3D periodic lattice."""
    # Laplacian: sum of (φ_j - φ_i) for 6 neighbors = -6φ + Σφ_j
    lap = -6.0 * phi
    lap += cp.roll(phi, 1, axis=0) + cp.roll(phi, -1, axis=0)
    lap += cp.roll(phi, 1, axis=1) + cp.roll(phi, -1, axis=1)
    lap += cp.roll(phi, 1, axis=2) + cp.roll(phi, -1, axis=2)
    # Potential force: -(1/π)sin(πφ)
    pot = -(1.0/PI) * cp.sin(PI * phi)
    return lap + pot

def compute_energies_gpu(phi, phi_dot):
    """Kinetic, potential, and total energy."""
    # Kinetic: Σ φ̇²/2
    KE = 0.5 * float(cp.sum(phi_dot**2))
    # Potential (gradient part): Σ (φ_i - φ_j)²/2
    PE_grad = 0.0
    for ax in range(3):
        diff = phi - cp.roll(phi, 1, axis=ax)
        PE_grad += 0.5 * float(cp.sum(diff**2))
    # Potential (on-site): Σ (1/π²)(1 - cos(πφ))
    PE_site = float(cp.sum((1.0/PI**2) * (1 - cp.cos(PI * phi))))
    return KE, PE_grad + PE_site, KE + PE_grad + PE_site

def count_winding(phi_np, cx, cy, R_maj, N):
    """Count the poloidal winding number around a torus.
    Sample φ around the tube cross-section at one toroidal angle."""
    # Sample at θ_tor = 0 (positive x-axis), going around the tube
    mid = N // 2
    ix_ring = int(cx + R_maj + mid)  # x-position on the ring
    iy_ring = int(cy + mid)
    if ix_ring >= N: ix_ring = ix_ring % N

    # Extract φ along z at this (x,y) position = poloidal cut
    col = phi_np[ix_ring, iy_ring, :]
    # The winding = (φ at top - φ at bottom) / 2
    # For a winding of 1: should be ≈ 2
    z_arr = np.arange(N) - N/2
    top = col[N//2 + r_tube]
    bot = col[N//2 - r_tube]
    winding = (top - bot) / 2.0
    return winding

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

# ============================================================
# SINGLE TORUS — SETTLE TO EQUILIBRIUM
# ============================================================
report("SINGLE TORUS — DAMPED EVOLUTION TO EQUILIBRIUM")
report("-" * 55)

phi_init_s, _ = torus_field_and_distance(X, Y, Z, 0, 0)
phi_gpu = cp.asarray(phi_init_s)
phi_dot_gpu = cp.zeros_like(phi_gpu)  # start at rest

n_steps = 2000
report(f"Evolving for {n_steps} steps with dt={dt}, damping={damping}")

t0 = time.time()
for step in range(n_steps):
    accel = compute_acceleration_gpu(phi_gpu)
    phi_dot_gpu += dt * accel
    phi_dot_gpu *= (1 - damping)  # damping
    phi_gpu += dt * phi_dot_gpu

    if step % 500 == 0 or step == n_steps - 1:
        KE, PE, TE = compute_energies_gpu(phi_gpu, phi_dot_gpu)
        phi_np = cp.asnumpy(phi_gpu)
        w = count_winding(phi_np, 0, 0, R_maj, N)
        report(f"  step {step:5d}: KE={KE:.2f}, PE={PE:.2f}, TE={TE:.2f}, "
               f"winding={w:.3f}")

t_evolve = time.time() - t0
report(f"  Evolution time: {t_evolve:.1f}s")
report("")

phi_settled_s = cp.asnumpy(phi_gpu)
report(f"Settled field range: [{phi_settled_s.min():.4f}, {phi_settled_s.max():.4f}]")

# Winding check
w_s = count_winding(phi_settled_s, 0, 0, R_maj, N)
report(f"Winding number: {w_s:.3f} (should be ≈ 1.0)")

# Eigenvalues — k=6 (enough for first 3 pairs), sigma=0 for shift-invert
report("Computing eigenvalues...")
H_s = build_hessian_gpu(phi_settled_s.ravel())
evals_s_gpu, _ = csla.eigsh(H_s, k=6, which='SA')
cp.cuda.Stream.null.synchronize()
evals_s = np.sort(cp.asnumpy(evals_s_gpu))
n_neg_s = np.sum(evals_s < -0.01)
report(f"Tachyons: {n_neg_s}")
report(f"Lowest 8: {' '.join(f'{e:+.6f}' for e in evals_s[:8])}")
report("")

# ============================================================
# TWO TORI — SETTLE AND COMPUTE BOND
# ============================================================
report("TWO TORI — DAMPED EVOLUTION AT EACH R")
report("-" * 55)

R_values = [28, 26, 24, 22, 20, 19, 18, 17]

report(f"{'R':>3} {'gap':>4} {'KE_final':>10} {'PE_final':>10} {'n_neg':>5} "
       f"{'ev0':>10} {'V0':>10} {'w_A':>6} {'w_B':>6}")
report("-" * 72)

bond_data = {}

for R in R_values:
    # Initial condition: nearest-torus
    phi_A, rho_A = torus_field_and_distance(X, Y, Z, -R/2, 0)
    phi_B, rho_B = torus_field_and_distance(X, Y, Z, +R/2, 0)
    phi_init = np.where(rho_A <= rho_B, phi_A, phi_B)

    phi_gpu = cp.asarray(phi_init)
    phi_dot_gpu = cp.zeros_like(phi_gpu)

    # Evolve
    for step in range(n_steps):
        accel = compute_acceleration_gpu(phi_gpu)
        phi_dot_gpu += dt * accel
        phi_dot_gpu *= (1 - damping)
        phi_gpu += dt * phi_dot_gpu

    phi_settled = cp.asnumpy(phi_gpu)
    KE, PE, TE = compute_energies_gpu(phi_gpu, phi_dot_gpu)

    # Check windings
    w_A = count_winding(phi_settled, -R/2, 0, R_maj, N)
    w_B = count_winding(phi_settled, +R/2, 0, R_maj, N)

    # Eigenvalues
    ev_gpu = csla.eigsh(build_hessian_gpu(phi_settled.ravel()),
                        k=6, which='SA')
    cp.cuda.Stream.null.synchronize()
    ev = np.sort(cp.asnumpy(ev_gpu[0]))
    n_neg = np.sum(ev < -0.01)
    V0 = ev[0] - evals_s[0]

    bond_data[R] = {'ev': ev, 'V0': V0, 'n_neg': n_neg,
                     'KE': KE, 'PE': PE, 'w_A': w_A, 'w_B': w_B}

    report(f"{R:3d} {R-2*R_maj:4d} {KE:10.2f} {PE:10.2f} {n_neg:5d} "
           f"{ev[0]:+10.6f} {V0:+10.6f} {w_A:6.3f} {w_B:6.3f}")

report("")

# ============================================================
# BOND CURVE
# ============================================================
report("BOND CURVE")
report("-" * 40)

for R in R_values:
    V = bond_data[R]['V0']
    w = f"w={bond_data[R]['w_A']:.2f},{bond_data[R]['w_B']:.2f}"
    marker = " <-- min" if V == min(bond_data[r]['V0'] for r in R_values) else ""
    report(f"  R={R:3d}: V={V:+.8f}  {w}{marker}")

V_min = min(bond_data[R]['V0'] for R in R_values)
report("")
if V_min < -1e-6:
    report(f"BONDING: D_e = {-V_min:.8f}")
else:
    report("No clear bonding from V0 shifts.")

# Check: did both windings survive?
report("")
report("WINDING SURVIVAL:")
for R in R_values:
    wA = bond_data[R]['w_A']
    wB = bond_data[R]['w_B']
    status = "BOTH" if abs(wA) > 0.5 and abs(wB) > 0.5 else \
             "ONE LOST" if abs(wA) > 0.5 or abs(wB) > 0.5 else "BOTH LOST"
    report(f"  R={R}: w_A={wA:.3f}, w_B={wB:.3f} — {status}")

# Tachyon check
report("")
any_tach = any(bond_data[R]['n_neg'] > 0 for R in R_values)
if not any_tach:
    report("NO TACHYONS — clean topological bonding!")
else:
    for R in R_values:
        if bond_data[R]['n_neg'] > 0:
            report(f"  R={R}: {bond_data[R]['n_neg']} tachyon(s)")

# Splittings
report("")
report("EIGENVALUE PAIR SPLITTINGS:")
report(f"{'R':>3} {'split0':>10} {'split1':>10} {'split2':>10}")
report("-" * 36)
for R in R_values:
    ev = bond_data[R]['ev']
    s0 = ev[1]-ev[0] if len(ev)>1 else 0
    s1 = ev[3]-ev[2] if len(ev)>3 else 0
    s2 = ev[5]-ev[4] if len(ev)>5 else 0
    report(f"{R:3d} {s0:10.6f} {s1:10.6f} {s2:10.6f}")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
