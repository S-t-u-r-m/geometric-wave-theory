"""
Breather Node Interaction — GPU
==================================

Tests how radial nodes affect breather-breather interaction energy.
Uses exact DHN sine-Gordon breather profiles (not Gaussians).

The GWT breather profile is:
  phi(r) = (4/pi) * arctan[ eta / (w * cosh(eta * r)) ]
where w = cos(n*gamma), eta = sin(n*gamma), gamma = pi/(N+1).

This is nodeless. To model atoms with radial nodes (2s, 3s),
we multiply by a node function:
  phi_noded(r) = phi(r) * cos(k_node * pi * r / r_scale)
where k_node = number of nodes.

The node scale r_scale is set so nodes fall within the breather width.
For breather of width ~1/eta, r_scale = 3/eta puts the first node
at r ≈ 1.5/eta (roughly half the breather extent).

Key questions:
1. Does |sin(phase)| emerge for nodeless breather pairs?
2. What suppression do nodes cause vs the nodeless baseline?
3. Is asymmetric (noded+nodeless) different from symmetric (noded+noded)?
4. What is the actual suppression factor? (V6 claims S/n_lobes)
"""

import cupy as cp
import numpy as np
import time as timer

pi = float(np.pi)
V_0 = 1.0 / pi**2
N_br = 24
gamma_br = pi / (N_br + 1)


def run_breather_interaction(N_grid, L, n_mode, R_sep, node1=0, node2=0,
                              l1=0, l2=0, n_steps=3000, dt_factor=0.15):
    """
    Compute interaction energy for two breather-like waves.

    n_mode: DHN mode number (sets internal frequency/width)
    node1, node2: number of radial nodes (0=1s-like, 1=2s-like, 2=3s-like)
    l1, l2: angular momentum (0=s, 1=p)
    """
    dx = L / N_grid
    dt = dt_factor * dx
    coords = cp.linspace(-L/2, L/2, N_grid, dtype=cp.float32)
    X, Y, Z = cp.meshgrid(coords, coords, coords, indexing='ij')

    w = float(np.cos(n_mode * gamma_br))
    eta = float(np.sin(n_mode * gamma_br))

    def make_breather(center_z, n_nodes, l_ang, sign=1.0):
        """Create a breather field with optional radial nodes."""
        rz = Z - center_z
        r = cp.sqrt(X**2 + Y**2 + rz**2)
        r = cp.maximum(r, dx/2)

        # Base breather profile (nodeless)
        profile = sign * (4.0/pi) * cp.arctan(
            cp.float32(eta) / (cp.float32(w) * cp.cosh(cp.float32(eta) * r))
        )

        # Add radial nodes
        if n_nodes > 0:
            r_scale = 3.0 / eta  # node spacing relative to breather width
            profile *= cp.cos(cp.float32(n_nodes) * cp.float32(pi) * r / cp.float32(r_scale))

        # Angular structure
        if l_ang == 1:
            cos_theta = rz / r
            profile *= cos_theta

        return profile

    def compute_energy(phi):
        gx = cp.diff(phi, axis=0)
        gy = cp.diff(phi, axis=1)
        gz = cp.diff(phi, axis=2)
        E_grad = 0.5 * (cp.sum(gx**2) + cp.sum(gy**2) + cp.sum(gz**2)) / dx * dx**2
        E_pot = V_0 * cp.sum(1 - cp.cos(pi * phi)) * dx**3
        return float(E_grad + E_pot)

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
                E_pot = compute_energy(phi)
                E_samples.append(E_kin + E_pot)
        return np.mean(E_samples) if E_samples else compute_energy(phi_init)

    phi1 = make_breather(-R_sep/2, node1, l1, sign=1.0)
    E1 = evolve_and_measure(phi1)

    phi2 = make_breather(+R_sep/2, node2, l2, sign=1.0)
    E2 = evolve_and_measure(phi2)

    phi_both = make_breather(-R_sep/2, node1, l1, sign=1.0) + \
               make_breather(+R_sep/2, node2, l2, sign=1.0)
    E_both = evolve_and_measure(phi_both)

    E_int = E_both - E1 - E2
    cp.cuda.Stream.null.synchronize()
    return E_int, E1, E2


# =============================================================================
print("=" * 80)
print("  BREATHER NODE INTERACTION SIMULATION")
print("  Exact DHN breather profiles with radial nodes")
print("=" * 80)

_ = cp.zeros(10)
cp.cuda.Stream.null.synchronize()

N = 60
L = 20.0
n_steps = 3000

# Use mode n=1 (widest breather, best resolved on grid)
n_mode = 1
w_val = np.cos(n_mode * gamma_br)
eta_val = np.sin(n_mode * gamma_br)
print(f"\n  Mode n={n_mode}: w={w_val:.4f}, eta={eta_val:.4f}, width~{1/eta_val:.2f}")
print(f"  Grid: {N}^3, L={L}, dx={L/N:.3f}")
print(f"  Steps: {n_steps}")


# =============================================================================
# TEST 1: Nodeless baseline — 0-node + 0-node vs R
# This should show whether |sin(phase)| emerges for breathers
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 1: Nodeless s-wave breather pair (0+0 nodes) vs R")
print(f"  Baseline for all node comparisons")
print(f"{'='*80}")

R_values = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0]

print(f"\n  {'R':>5}  {'E_int':>10}  {'E1':>10}  {'E2':>10}")

results_00 = {}
for R in R_values:
    t0 = timer.time()
    E_int, E1, E2 = run_breather_interaction(N, L, n_mode, R, node1=0, node2=0,
                                               n_steps=n_steps)
    dt = timer.time() - t0
    results_00[R] = E_int
    print(f"  {R:5.1f}  {E_int:10.4f}  {E1:10.4f}  {E2:10.4f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 2: 1-node + 0-node (ASYMMETRIC) — models LiH (2s+1s)
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 2: 1-node + 0-node (asymmetric, like LiH 2s+1s) vs R")
print(f"{'='*80}")

print(f"\n  {'R':>5}  {'E_int':>10}  {'E_00':>10}  {'ratio':>8}")

results_10 = {}
for R in R_values:
    t0 = timer.time()
    E_int, _, _ = run_breather_interaction(N, L, n_mode, R, node1=1, node2=0,
                                            n_steps=n_steps)
    dt = timer.time() - t0
    E_00 = results_00[R]
    ratio = E_int / E_00 if abs(E_00) > 1e-6 else float('inf')
    results_10[R] = E_int
    print(f"  {R:5.1f}  {E_int:10.4f}  {E_00:10.4f}  {ratio:8.4f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 3: 1-node + 1-node (SYMMETRIC) — models Li2 (2s+2s)
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 3: 1-node + 1-node (symmetric, like Li2 2s+2s) vs R")
print(f"{'='*80}")

print(f"\n  {'R':>5}  {'E_int':>10}  {'E_00':>10}  {'ratio':>8}")

results_11 = {}
for R in R_values:
    t0 = timer.time()
    E_int, _, _ = run_breather_interaction(N, L, n_mode, R, node1=1, node2=1,
                                            n_steps=n_steps)
    dt = timer.time() - t0
    E_00 = results_00[R]
    ratio = E_int / E_00 if abs(E_00) > 1e-6 else float('inf')
    results_11[R] = E_int
    print(f"  {R:5.1f}  {E_int:10.4f}  {E_00:10.4f}  {ratio:8.4f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 4: 2-node + 0-node (ASYMMETRIC) — models NaH (3s+1s)
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 4: 2-node + 0-node (asymmetric, like NaH 3s+1s) vs R")
print(f"{'='*80}")

print(f"\n  {'R':>5}  {'E_int':>10}  {'E_00':>10}  {'ratio':>8}")

results_20 = {}
for R in R_values:
    t0 = timer.time()
    E_int, _, _ = run_breather_interaction(N, L, n_mode, R, node1=2, node2=0,
                                            n_steps=n_steps)
    dt = timer.time() - t0
    E_00 = results_00[R]
    ratio = E_int / E_00 if abs(E_00) > 1e-6 else float('inf')
    results_20[R] = E_int
    print(f"  {R:5.1f}  {E_int:10.4f}  {E_00:10.4f}  {ratio:8.4f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 5: 2-node + 2-node (SYMMETRIC) — models Na2 (3s+3s)
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 5: 2-node + 2-node (symmetric, like Na2 3s+3s) vs R")
print(f"{'='*80}")

print(f"\n  {'R':>5}  {'E_int':>10}  {'E_00':>10}  {'ratio':>8}")

results_22 = {}
for R in R_values:
    t0 = timer.time()
    E_int, _, _ = run_breather_interaction(N, L, n_mode, R, node1=2, node2=2,
                                            n_steps=n_steps)
    dt = timer.time() - t0
    E_00 = results_00[R]
    ratio = E_int / E_00 if abs(E_00) > 1e-6 else float('inf')
    results_22[R] = E_int
    print(f"  {R:5.1f}  {E_int:10.4f}  {E_00:10.4f}  {ratio:8.4f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 6: Summary comparison at multiple R
# =============================================================================
print(f"\n{'='*80}")
print(f"  SUMMARY: Node suppression ratios (E_noded / E_nodeless)")
print(f"  sym = both same nodes, asym = one noded one not")
print(f"{'='*80}")

print(f"\n  {'R':>5}  {'0+0':>10}  {'1+0 asym':>10}  {'1+1 sym':>10}  {'2+0 asym':>10}  {'2+2 sym':>10}")

for R in R_values:
    E00 = results_00[R]
    r10 = results_10[R] / E00 if abs(E00) > 1e-6 else 0
    r11 = results_11[R] / E00 if abs(E00) > 1e-6 else 0
    r20 = results_20[R] / E00 if abs(E00) > 1e-6 else 0
    r22 = results_22[R] / E00 if abs(E00) > 1e-6 else 0
    print(f"  {R:5.1f}  {E00:10.4f}  {r10:10.4f}  {r11:10.4f}  {r20:10.4f}  {r22:10.4f}")


# =============================================================================
# TEST 7: Higher mode (n=3) — more compact breather
# Compare nodeless vs 1-node at a few R values
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 7: Mode n=3 (more compact breather)")
n_mode_3 = 3
w3 = np.cos(n_mode_3 * gamma_br)
eta3 = np.sin(n_mode_3 * gamma_br)
print(f"  w={w3:.4f}, eta={eta3:.4f}, width~{1/eta3:.2f}")
print(f"{'='*80}")

R_test = [3.0, 4.0, 5.0, 6.0, 8.0]
print(f"\n  {'R':>5}  {'E_00':>10}  {'E_10':>10}  {'ratio_10':>10}  {'E_11':>10}  {'ratio_11':>10}")

for R in R_test:
    t0 = timer.time()
    E00, _, _ = run_breather_interaction(N, L, n_mode_3, R, node1=0, node2=0, n_steps=n_steps)
    E10, _, _ = run_breather_interaction(N, L, n_mode_3, R, node1=1, node2=0, n_steps=n_steps)
    E11, _, _ = run_breather_interaction(N, L, n_mode_3, R, node1=1, node2=1, n_steps=n_steps)
    dt = timer.time() - t0

    r10 = E10 / E00 if abs(E00) > 1e-6 else 0
    r11 = E11 / E00 if abs(E00) > 1e-6 else 0
    print(f"  {R:5.1f}  {E00:10.4f}  {E10:10.4f}  {r10:10.4f}  {E11:10.4f}  {r11:10.4f}  ({dt:.1f}s)", flush=True)


cp.cuda.Stream.null.synchronize()
print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
