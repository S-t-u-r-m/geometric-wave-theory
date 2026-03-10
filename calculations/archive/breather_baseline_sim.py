"""
Breather Baseline — Does |sin(phase)| emerge?
================================================

Dense R scan with compact breather (mode n=5, width ~1.7).
Tests the FOUNDATION of the V6 formula: E_int ~ |sin(kR)|.

The breather has internal spatial frequency eta = sin(n*gamma).
Two breathers separated by R should interact with oscillatory
dependence on R due to the standing wave structure.

If |sin(phase)| works, E_int/envelope should oscillate sinusoidally.
If not, the entire V6 formula needs rethinking.

Also tests: do nodes change the oscillation, or just suppress it?
"""

import cupy as cp
import numpy as np
import time as timer

pi = float(np.pi)
V_0 = 1.0 / pi**2
N_br = 24
gamma_br = pi / (N_br + 1)


def run_pair(N_grid, L, n_mode, R_sep, node1=0, node2=0,
             l1=0, l2=0, n_steps=3000, dt_factor=0.15):
    """Compute interaction energy for breather pair."""
    dx = L / N_grid
    dt = dt_factor * dx
    coords = cp.linspace(-L/2, L/2, N_grid, dtype=cp.float32)
    X, Y, Z = cp.meshgrid(coords, coords, coords, indexing='ij')

    w = float(np.cos(n_mode * gamma_br))
    eta = float(np.sin(n_mode * gamma_br))

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

    phi1 = make_breather(-R_sep/2, node1, l1)
    E1 = evolve_and_measure(phi1)

    phi2 = make_breather(+R_sep/2, node2, l2)
    E2 = evolve_and_measure(phi2)

    phi_both = make_breather(-R_sep/2, node1, l1) + make_breather(+R_sep/2, node2, l2)
    E_both = evolve_and_measure(phi_both)

    E_int = E_both - E1 - E2
    cp.cuda.Stream.null.synchronize()
    return E_int


# =============================================================================
print("=" * 80)
print("  BREATHER BASELINE: Does |sin(phase)| emerge?")
print("=" * 80)

_ = cp.zeros(10)
cp.cuda.Stream.null.synchronize()

N = 64
L = 24.0
n_steps = 3000

# Mode n=5: compact enough to resolve, wide enough for good grid sampling
n_mode = 5
w_val = np.cos(n_mode * gamma_br)
eta_val = np.sin(n_mode * gamma_br)
print(f"\n  Mode n={n_mode}: w={w_val:.4f}, eta={eta_val:.4f}, width~{1/eta_val:.2f}")
print(f"  Grid: {N}^3, L={L}, dx={L/N:.3f}")
print(f"  Breather half-wavelength: pi/eta = {pi/eta_val:.2f}")
print(f"  Expect oscillation period in R: ~pi/eta = {pi/eta_val:.2f}")

# Also check mode n=3 for comparison
n_mode_3 = 3
w3 = np.cos(n_mode_3 * gamma_br)
eta3 = np.sin(n_mode_3 * gamma_br)
print(f"\n  Mode n={n_mode_3}: w={w3:.4f}, eta={eta3:.4f}, width~{1/eta3:.2f}")
print(f"  Breather half-wavelength: pi/eta = {pi/eta3:.2f}")


# =============================================================================
# TEST 1: Dense R scan, s-wave, nodeless, mode n=5
# If |sin(phase)| works, E_int should oscillate with period ~pi/eta
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 1: Dense R scan — mode n={n_mode} s-wave, 0 nodes")
print(f"  Oscillation period if |sin(eta*R)|: T = pi/eta = {pi/eta_val:.2f}")
print(f"  Zeros at R = k*pi/eta = {pi/eta_val:.2f}, {2*pi/eta_val:.2f}, ...")
print(f"{'='*80}")

# Dense scan from R=1 to R=14
R_dense = np.arange(1.0, 14.5, 0.5)

print(f"\n  {'R':>5}  {'E_int':>10}  {'|sin(eta*R)|':>12}  {'|sin(2etaR)|':>12}  {'E/decay':>10}")

results_n5 = []
for R in R_dense:
    t0 = timer.time()
    E_int = run_pair(N, L, n_mode, float(R), node1=0, node2=0, n_steps=n_steps)
    dt = timer.time() - t0

    # Candidate phases
    s1 = abs(np.sin(eta_val * R))          # single eta
    s2 = abs(np.sin(2 * eta_val * R))      # double eta (R/n1 + R/n2 = 2*eta*R for same mode)

    # Remove exponential decay to see oscillation
    decay = np.exp(-eta_val * R)
    E_over_decay = E_int / decay if decay > 1e-10 else float('inf')

    results_n5.append((R, E_int, s1, s2, decay, E_over_decay))
    print(f"  {R:5.1f}  {E_int:10.4f}  {s1:12.4f}  {s2:12.4f}  {E_over_decay:10.2f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 2: Dense R scan, s-wave, nodeless, mode n=3
# Wider breather — check if same oscillation pattern
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 2: Dense R scan — mode n={n_mode_3} s-wave, 0 nodes")
print(f"  Oscillation period if |sin(eta*R)|: T = pi/eta = {pi/eta3:.2f}")
print(f"{'='*80}")

R_dense_3 = np.arange(1.0, 14.5, 0.5)

print(f"\n  {'R':>5}  {'E_int':>10}  {'|sin(eta*R)|':>12}  {'|sin(2etaR)|':>12}  {'E/decay':>10}")

results_n3 = []
for R in R_dense_3:
    t0 = timer.time()
    E_int = run_pair(N, L, n_mode_3, float(R), node1=0, node2=0, n_steps=n_steps)
    dt = timer.time() - t0

    s1 = abs(np.sin(eta3 * R))
    s2 = abs(np.sin(2 * eta3 * R))
    decay = np.exp(-eta3 * R)
    E_over_decay = E_int / decay if decay > 1e-10 else float('inf')

    results_n3.append((R, E_int, s1, s2, decay, E_over_decay))
    print(f"  {R:5.1f}  {E_int:10.4f}  {s1:12.4f}  {s2:12.4f}  {E_over_decay:10.2f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 3: p-wave (l=1) dense scan, mode n=5
# Previous 3D GPU sim showed p-wave oscillates — verify with dense scan
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 3: Dense R scan — mode n={n_mode} P-WAVE (l=1), 0 nodes")
print(f"  p-wave should show clearer oscillation (cos_theta factor)")
print(f"{'='*80}")

print(f"\n  {'R':>5}  {'E_int':>10}  {'|sin(eta*R)|':>12}  {'|sin(2etaR)|':>12}  {'E/decay':>10}")

results_p5 = []
for R in R_dense:
    t0 = timer.time()
    E_int = run_pair(N, L, n_mode, float(R), node1=0, node2=0,
                     l1=1, l2=1, n_steps=n_steps)
    dt = timer.time() - t0

    s1 = abs(np.sin(eta_val * R))
    s2 = abs(np.sin(2 * eta_val * R))
    decay = np.exp(-eta_val * R)
    E_over_decay = E_int / decay if decay > 1e-10 else float('inf')

    results_p5.append((R, E_int, s1, s2, decay, E_over_decay))
    print(f"  {R:5.1f}  {E_int:10.4f}  {s1:12.4f}  {s2:12.4f}  {E_over_decay:10.2f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 4: 1-node vs 0-node at select R values, mode n=5
# Quick check of node suppression with compact breather
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 4: Node suppression — mode n={n_mode}, select R values")
print(f"{'='*80}")

R_node_test = [3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
print(f"\n  {'R':>5}  {'E_00':>10}  {'E_10':>10}  {'E_11':>10}  {'r_10':>8}  {'r_11':>8}")

for R in R_node_test:
    t0 = timer.time()
    E00 = run_pair(N, L, n_mode, R, node1=0, node2=0, n_steps=n_steps)
    E10 = run_pair(N, L, n_mode, R, node1=1, node2=0, n_steps=n_steps)
    E11 = run_pair(N, L, n_mode, R, node1=1, node2=1, n_steps=n_steps)
    dt = timer.time() - t0

    r10 = E10 / E00 if abs(E00) > 1e-6 else 0
    r11 = E11 / E00 if abs(E00) > 1e-6 else 0
    print(f"  {R:5.1f}  {E00:10.4f}  {E10:10.4f}  {E11:10.4f}  {r10:8.4f}  {r11:8.4f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# ANALYSIS
# =============================================================================
print(f"\n{'='*80}")
print(f"  ANALYSIS: Looking for oscillation in E_int(R)")
print(f"{'='*80}")

print(f"\n  Mode n={n_mode} (eta={eta_val:.4f}):")
print(f"  Zero crossings of E_int (sign changes):")
for i in range(1, len(results_n5)):
    R_prev, E_prev = results_n5[i-1][0], results_n5[i-1][1]
    R_curr, E_curr = results_n5[i][0], results_n5[i][1]
    if E_prev * E_curr < 0:
        # Linear interpolation for zero crossing
        R_zero = R_prev + (R_curr - R_prev) * abs(E_prev) / (abs(E_prev) + abs(E_curr))
        print(f"    R ~ {R_zero:.2f}  (between {R_prev:.1f} and {R_curr:.1f})")
        print(f"    Expected from sin(eta*R)=0: R = k*pi/eta, nearest = {round(R_zero*eta_val/pi)*pi/eta_val:.2f}")
        print(f"    Expected from sin(2*eta*R)=0: R = k*pi/(2*eta), nearest = {round(R_zero*2*eta_val/pi)*pi/(2*eta_val):.2f}")

print(f"\n  Mode n={n_mode_3} (eta={eta3:.4f}):")
print(f"  Zero crossings of E_int (sign changes):")
for i in range(1, len(results_n3)):
    R_prev, E_prev = results_n3[i-1][0], results_n3[i-1][1]
    R_curr, E_curr = results_n3[i][0], results_n3[i][1]
    if E_prev * E_curr < 0:
        R_zero = R_prev + (R_curr - R_prev) * abs(E_prev) / (abs(E_prev) + abs(E_curr))
        print(f"    R ~ {R_zero:.2f}  (between {R_prev:.1f} and {R_curr:.1f})")
        print(f"    Expected from sin(eta*R)=0: R = k*pi/eta, nearest = {round(R_zero*eta3/pi)*pi/eta3:.2f}")

print(f"\n  p-wave mode n={n_mode}:")
print(f"  Zero crossings:")
for i in range(1, len(results_p5)):
    R_prev, E_prev = results_p5[i-1][0], results_p5[i-1][1]
    R_curr, E_curr = results_p5[i][0], results_p5[i][1]
    if E_prev * E_curr < 0:
        R_zero = R_prev + (R_curr - R_prev) * abs(E_prev) / (abs(E_prev) + abs(E_curr))
        print(f"    R ~ {R_zero:.2f}  (between {R_prev:.1f} and {R_curr:.1f})")


cp.cuda.Stream.null.synchronize()
print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
