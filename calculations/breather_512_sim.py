"""
Breather 512^3 Large-Box Simulation
====================================

512^3 grid, L=60 box — large enough to see true asymptotic behavior.
Sparse R scan first to check if |sin(phase)| and exp(-eta*R) decay emerge.

Mode n=5: eta=0.5878, width~1.7, half-wavelength pi/eta=5.34
Box half-width = 30, so R up to ~40 should be clean of boundary effects.
"""

import cupy as cp
import numpy as np
import time as timer
import gc

pi = float(np.pi)
V_0 = 1.0 / pi**2
N_br = 24
gamma_br = pi / (N_br + 1)


def run_pair(N_grid, L, n_mode, R_sep, node1=0, node2=0,
             l1=0, l2=0, n_steps=3000, dt_factor=0.15):
    """Compute interaction energy for breather pair. Memory-optimized."""
    dx = L / N_grid
    dt = dt_factor * dx

    w = float(np.cos(n_mode * gamma_br))
    eta = float(np.sin(n_mode * gamma_br))

    # Build coordinate grid
    coords = cp.linspace(-L/2, L/2, N_grid, dtype=cp.float32)
    X, Y, Z = cp.meshgrid(coords, coords, coords, indexing='ij')
    del coords

    def make_breather(center_z, n_nodes, l_ang, sign=1.0):
        rz = Z - center_z
        r = cp.sqrt(X**2 + Y**2 + rz**2)
        r = cp.maximum(r, dx/2)

        profile = sign * (4.0/pi) * cp.arctan(
            cp.float32(eta) / (cp.float32(w) * cp.cosh(cp.float32(eta) * r))
        )

        if n_nodes > 0:
            r_scale = 3.0 / eta
            profile *= cp.cos(cp.float32(n_nodes) * cp.float32(pi) * r / cp.float32(r_scale))

        if l_ang == 1:
            profile *= (rz / r)

        del rz, r
        return profile

    def compute_energy(phi):
        gx = cp.diff(phi, axis=0)
        E_grad = 0.5 * cp.sum(gx**2)
        del gx
        gy = cp.diff(phi, axis=1)
        E_grad += 0.5 * cp.sum(gy**2)
        del gy
        gz = cp.diff(phi, axis=2)
        E_grad += 0.5 * cp.sum(gz**2)
        del gz
        E_grad = float(E_grad) / dx * dx**2
        E_pot = float(V_0 * cp.sum(1 - cp.cos(pi * phi))) * dx**3
        return E_grad + E_pot

    def compute_force(phi):
        force = cp.zeros_like(phi)
        force[1:-1,:,:] += phi[:-2,:,:] + phi[2:,:,:] - 2*phi[1:-1,:,:]
        force[:,1:-1,:] += phi[:,:-2,:] + phi[:,2:,:] - 2*phi[:,1:-1,:]
        force[:,:,1:-1] += phi[:,:,:-2] + phi[:,:,2:] - 2*phi[:,:,1:-1]
        force /= dx**2
        force -= V_0 * pi * cp.sin(pi * phi)
        return force

    def evolve_and_measure(phi_init):
        phi = phi_init.copy()
        dphi = cp.zeros_like(phi)
        force = compute_force(phi)
        dphi += 0.5 * dt * force
        E_samples = []
        for step in range(n_steps):
            phi += dt * dphi
            force = compute_force(phi)
            dphi += dt * force
            if step % 50 == 0 and step > n_steps // 4:
                dphi_half = dphi - 0.5 * dt * force
                E_kin = 0.5 * float(cp.sum(dphi_half**2)) * dx**3
                del dphi_half
                E_pot = compute_energy(phi)
                E_samples.append(E_kin + E_pot)
        del phi, dphi, force
        return np.mean(E_samples) if E_samples else compute_energy(phi_init)

    # Build and measure each config, freeing as we go
    phi1 = make_breather(-R_sep/2, node1, l1)
    E1 = evolve_and_measure(phi1)
    del phi1
    cp.get_default_memory_pool().free_all_blocks()

    phi2 = make_breather(+R_sep/2, node2, l2)
    E2 = evolve_and_measure(phi2)
    del phi2
    cp.get_default_memory_pool().free_all_blocks()

    phi_both = make_breather(-R_sep/2, node1, l1) + make_breather(+R_sep/2, node2, l2)
    E_both = evolve_and_measure(phi_both)
    del phi_both
    cp.get_default_memory_pool().free_all_blocks()

    # Free the meshgrid
    del X, Y, Z
    cp.get_default_memory_pool().free_all_blocks()

    E_int = E_both - E1 - E2
    return E_int


# =============================================================================
print("=" * 80)
print("  BREATHER 512^3 LARGE-BOX SIMULATION")
print("=" * 80)

# Warm up GPU
_ = cp.zeros(10)
cp.cuda.Stream.null.synchronize()

# Check available memory
mem_free = cp.cuda.Device().mem_info[0] / 1e9
mem_total = cp.cuda.Device().mem_info[1] / 1e9
print(f"\n  GPU memory: {mem_free:.1f} GB free / {mem_total:.1f} GB total")

N = 512
L = 60.0
n_steps = 3000
n_mode = 5

w_val = np.cos(n_mode * gamma_br)
eta_val = np.sin(n_mode * gamma_br)

print(f"\n  Mode n={n_mode}: w={w_val:.4f}, eta={eta_val:.4f}, width~{1/eta_val:.2f}")
print(f"  Grid: {N}^3 = {N**3/1e6:.0f}M cells, L={L}, dx={L/N:.4f}")
print(f"  Breather width in grid points: {1/eta_val / (L/N):.1f}")
print(f"  Box half-width: {L/2:.0f}  (breathers safe up to R~{L-2/eta_val:.0f})")
print(f"  Expected oscillation: period = pi/eta = {pi/eta_val:.2f}")
print(f"  Expected decay: exp(-eta*R), half-life R = ln(2)/eta = {np.log(2)/eta_val:.2f}")

est_mem = 12 * N**3 * 4 / 1e9
print(f"\n  Estimated peak VRAM: ~{est_mem:.1f} GB")
if est_mem > mem_free:
    print(f"  WARNING: May be tight! Consider N=384 if OOM.")

# =============================================================================
# PHASE 1: Sparse R scan — s-wave, nodeless, mode n=5
# 10 points from R=2 to R=30, check for oscillation + decay
# =============================================================================
print(f"\n{'='*80}")
print(f"  PHASE 1: Sparse R scan (s-wave, 0 nodes, mode n={n_mode})")
print(f"  Looking for: oscillation with period ~{pi/eta_val:.1f}, decay as exp(-{eta_val:.3f}*R)")
print(f"{'='*80}")

R_sparse = [2.0, 4.0, 6.0, 8.0, 10.0, 14.0, 18.0, 22.0, 26.0, 30.0]

print(f"\n  {'R':>5}  {'E_int':>12}  {'exp(-eta*R)':>12}  {'E/decay':>12}  {'|sin(eta*R)|':>12}")

results = []
for R in R_sparse:
    t0 = timer.time()
    E_int = run_pair(N, L, n_mode, float(R), node1=0, node2=0, n_steps=n_steps)
    elapsed = timer.time() - t0

    decay = np.exp(-eta_val * R)
    s1 = abs(np.sin(eta_val * R))
    E_over_decay = E_int / decay if decay > 1e-15 else float('inf')

    results.append((R, E_int, decay, E_over_decay, s1))
    print(f"  {R:5.1f}  {E_int:12.6f}  {decay:12.6f}  {E_over_decay:12.4f}  {s1:12.4f}  ({elapsed:.0f}s)", flush=True)


# =============================================================================
# Quick analysis after phase 1
# =============================================================================
print(f"\n{'='*80}")
print(f"  PHASE 1 ANALYSIS")
print(f"{'='*80}")

# Check if E_int decays
print(f"\n  Decay check (does E_int decrease with R?):")
for i in range(1, len(results)):
    R_prev, E_prev = results[i-1][0], abs(results[i-1][1])
    R_curr, E_curr = results[i][0], abs(results[i][1])
    ratio = E_curr / E_prev if E_prev > 1e-15 else 0
    expected = np.exp(-eta_val * (R_curr - R_prev))
    print(f"    R={R_prev:.0f}->{R_curr:.0f}: |E| ratio = {ratio:.4f}, expected exp(-eta*dR) = {expected:.4f}")

# Check for sign changes
print(f"\n  Sign changes (zero crossings):")
n_crossings = 0
for i in range(1, len(results)):
    if results[i-1][1] * results[i][1] < 0:
        R_prev, R_curr = results[i-1][0], results[i][0]
        print(f"    Between R={R_prev:.1f} and R={R_curr:.1f}")
        n_crossings += 1
if n_crossings == 0:
    print(f"    NONE -- E_int does not change sign")

# Check E/decay stability (if |sin| works, E/decay should oscillate around a constant)
E_over_decay_vals = [r[3] for r in results if abs(r[3]) < 1e10]
if len(E_over_decay_vals) > 2:
    print(f"\n  E/decay range: {min(E_over_decay_vals):.4f} to {max(E_over_decay_vals):.4f}")
    print(f"  If exp(-eta*R) is the right decay, these should be roughly constant")
    print(f"  Spread factor: {max(E_over_decay_vals)/min(E_over_decay_vals):.1f}x"
          if min(E_over_decay_vals) > 0 else "  (has negative values)")


# =============================================================================
# PHASE 2: If phase 1 shows decay, do dense scan near first expected zero
# Zero of sin(eta*R) at R = pi/eta = 5.34
# =============================================================================
print(f"\n{'='*80}")
print(f"  PHASE 2: Dense scan near first expected zero (R ~ {pi/eta_val:.1f})")
print(f"  6 points from R=3.5 to R=7.5")
print(f"{'='*80}")

R_dense = [3.5, 4.5, 5.0, 5.5, 6.5, 7.5]

print(f"\n  {'R':>5}  {'E_int':>12}  {'exp(-eta*R)':>12}  {'E/decay':>12}  {'|sin(eta*R)|':>12}")

for R in R_dense:
    t0 = timer.time()
    E_int = run_pair(N, L, n_mode, float(R), node1=0, node2=0, n_steps=n_steps)
    elapsed = timer.time() - t0

    decay = np.exp(-eta_val * R)
    s1 = abs(np.sin(eta_val * R))
    E_over_decay = E_int / decay if decay > 1e-15 else float('inf')

    results.append((R, E_int, decay, E_over_decay, s1))
    print(f"  {R:5.1f}  {E_int:12.6f}  {decay:12.6f}  {E_over_decay:12.4f}  {s1:12.4f}  ({elapsed:.0f}s)", flush=True)


# Sort all results by R for final table
results.sort(key=lambda x: x[0])

print(f"\n{'='*80}")
print(f"  COMBINED RESULTS (sorted by R)")
print(f"{'='*80}")
print(f"\n  {'R':>5}  {'E_int':>12}  {'decay':>12}  {'E/decay':>12}  {'|sin|':>8}")
for R, E, d, Ed, s in results:
    print(f"  {R:5.1f}  {E:12.6f}  {d:12.6f}  {Ed:12.4f}  {s:8.4f}")


cp.cuda.Stream.null.synchronize()
print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
