"""
3D Breather Interaction — GPU Accelerated (CuPy)
=================================================
Key question: Does angular interference in 3D produce
oscillatory sin(kR)-like interaction energy?

Uses GPU for ~50-100x speedup over CPU NumPy.
Includes progress tracking with ETA.
"""

import cupy as cp
import numpy as np
import time as timer
import sys

pi = float(np.pi)
V_0 = 1.0 / pi**2
N_br = 24
gamma_br = pi / (N_br + 1)

# Progress tracking
total_runs = 0
completed_runs = 0
start_time = None


def progress(desc=""):
    global completed_runs
    completed_runs += 1
    elapsed = timer.time() - start_time
    rate = completed_runs / elapsed if elapsed > 0 else 0
    remaining = (total_runs - completed_runs) / rate if rate > 0 else 0
    pct = 100 * completed_runs / total_runs
    mins_left = remaining / 60
    print(f"  [PROGRESS] {completed_runs}/{total_runs} ({pct:.0f}%) | "
          f"Elapsed: {elapsed/60:.1f}min | ETA: {mins_left:.1f}min | {desc}",
          flush=True)


def run_3d_gpu(N_grid, L, n1, n2, R_sep, sign2=+1,
               l1=0, l2=0, n_steps=2000, dt_factor=0.15):
    """Run 3D breather interaction on GPU."""
    dx = L / N_grid
    dt = dt_factor * dx

    coords = cp.linspace(-L/2, L/2, N_grid, dtype=cp.float32)
    X, Y, Z = cp.meshgrid(coords, coords, coords, indexing='ij')

    def breather_3d(X, Y, Z, n, center_x, amplitude, l_ang):
        w = float(np.cos(n * gamma_br))
        eta = float(np.sin(n * gamma_br))
        rx = X - center_x
        r = cp.sqrt(rx**2 + Y**2 + Z**2)
        r = cp.maximum(r, dx/2)
        profile = amplitude * (4.0/pi) * cp.arctan(eta / (w * cp.cosh(eta * r)))
        if l_ang == 1:
            cos_theta = rx / r
            profile = profile * cos_theta
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

    def evolve_and_measure(phi_init, n_steps):
        phi = phi_init.copy()
        dphi = cp.zeros_like(phi)
        force = compute_force(phi)
        dphi += 0.5 * dt * force
        E_samples = []
        for step in range(n_steps):
            phi += dt * dphi
            force = compute_force(phi)
            dphi += dt * force
            if step % 50 == 0 and step > n_steps//4:
                dphi_half = dphi - 0.5*dt*force
                E_kin = 0.5 * float(cp.sum(dphi_half**2)) * dx**3
                E_pg = compute_energy(phi)
                E_samples.append(E_kin + E_pg)
        return np.mean(E_samples) if E_samples else compute_energy(phi_init)

    # Single breather 1
    phi1 = breather_3d(X, Y, Z, n1, 0.0, 1.0, l1)
    E1 = evolve_and_measure(phi1, n_steps)

    # Single breather 2
    phi2 = breather_3d(X, Y, Z, n2, 0.0, sign2, l2)
    E2 = evolve_and_measure(phi2, n_steps)

    # Two-breather system
    phi_two = (breather_3d(X, Y, Z, n1, -R_sep/2, 1.0, l1) +
               breather_3d(X, Y, Z, n2, +R_sep/2, sign2, l2))
    E12 = evolve_and_measure(phi_two, n_steps)

    E_int = E12 - E1 - E2
    cp.cuda.Stream.null.synchronize()
    return E1, E2, E12, E_int


# =============================================================================
# Count total runs for progress tracking
# =============================================================================
# Section 1: s-wave - 3 pairs x 6 R x 2 phases = 36
# Section 2: p-wave - 2 pairs x 6 R x 2 phases = 24
# Section 3: mixed  - 9 R x 2 phases = 18
# Section 4: dense s-wave (1,1) - 23 R x 2 = 46
# Section 5: dense p-wave (3,3) - 23 R x 2 = 46
# Section 6: dense mixed 1s+3p  - 23 R x 2 = 46
total_runs = 36 + 24 + 18 + 46 + 46 + 46  # = 216

# =============================================================================
# MAIN
# =============================================================================
print("=" * 80)
print("  3D BREATHER INTERACTION — GPU ACCELERATED")
print("=" * 80)

# Warm up GPU
_ = cp.zeros(10)
cp.cuda.Stream.null.synchronize()

N = 60
L = 20.0
n_steps = 3000

print(f"  Grid: {N}x{N}x{N} = {N**3} points")
print(f"  Box: L={L}, dx={L/N:.3f}")
print(f"  Steps: {n_steps}")
print(f"  Total simulation runs: {total_runs}")
print(f"  Each run = 3 evolutions x {n_steps} steps on {N}^3 grid")
sys.stdout.flush()

start_time = timer.time()

# --- s-wave (l=0) breathers ---
print(f"\n{'='*80}")
print(f"  s-WAVE (l=0) BREATHER INTERACTION")
print(f"{'='*80}", flush=True)

for n1, n2, label in [(1,1,"mode 1+1"), (3,3,"mode 3+3"), (1,3,"mode 1+3")]:
    print(f"\n  --- ({n1},{n2}) {label}, s-wave ---")
    print(f"  {'R':>6} {'E1':>10} {'E2':>10} {'E12':>10} {'dE(+-)':>10} {'dE(++)':>10}", flush=True)

    for R in [2.0, 3.0, 4.0, 5.0, 7.0, 10.0]:
        t0 = timer.time()
        E1, E2, E12_pm, dE_pm = run_3d_gpu(N, L, n1, n2, R, sign2=-1,
                                             l1=0, l2=0, n_steps=n_steps)
        progress(f"({n1},{n2}) R={R} +-")
        _, _, E12_pp, dE_pp = run_3d_gpu(N, L, n1, n2, R, sign2=+1,
                                          l1=0, l2=0, n_steps=n_steps)
        progress(f"({n1},{n2}) R={R} ++")
        elapsed = timer.time() - t0
        print(f"  {R:6.1f} {E1:10.4f} {E2:10.4f} {E12_pm:10.4f} {dE_pm:+10.4f} {dE_pp:+10.4f}  ({elapsed:.1f}s)", flush=True)


# --- p-wave (l=1) breathers ---
print(f"\n{'='*80}")
print(f"  p-WAVE (l=1) BREATHER INTERACTION")
print(f"{'='*80}", flush=True)

for n1, n2, label in [(3,3,"mode 3+3 p-wave"), (1,3,"mode 1s+3p")]:
    l1_val = 1 if n1 >= 3 else 0
    l2_val = 1 if n2 >= 3 else 0
    print(f"\n  --- ({n1},{n2}) {label}, l=({l1_val},{l2_val}) ---")
    print(f"  {'R':>6} {'E1':>10} {'E2':>10} {'E12':>10} {'dE(+-)':>10} {'dE(++)':>10}", flush=True)

    for R in [2.0, 3.0, 4.0, 5.0, 7.0, 10.0]:
        t0 = timer.time()
        E1, E2, E12_pm, dE_pm = run_3d_gpu(N, L, n1, n2, R, sign2=-1,
                                             l1=l1_val, l2=l2_val, n_steps=n_steps)
        progress(f"p-wave ({n1},{n2}) R={R} +-")
        _, _, E12_pp, dE_pp = run_3d_gpu(N, L, n1, n2, R, sign2=+1,
                                          l1=l1_val, l2=l2_val, n_steps=n_steps)
        progress(f"p-wave ({n1},{n2}) R={R} ++")
        elapsed = timer.time() - t0
        print(f"  {R:6.1f} {E1:10.4f} {E2:10.4f} {E12_pm:10.4f} {dE_pm:+10.4f} {dE_pp:+10.4f}  ({elapsed:.1f}s)", flush=True)


# --- Mixed s+p ---
print(f"\n{'='*80}")
print(f"  MIXED s+p INTERACTION")
print(f"{'='*80}", flush=True)

print(f"\n  --- mode 1(s) + mode 3(p) ---")
print(f"  {'R':>6} {'dE(+-)':>10} {'dE(++)':>10}", flush=True)

for R in [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]:
    t0 = timer.time()
    _, _, _, dE_pm = run_3d_gpu(N, L, 1, 3, R, sign2=-1, l1=0, l2=1, n_steps=n_steps)
    progress(f"mixed s+p R={R} +-")
    _, _, _, dE_pp = run_3d_gpu(N, L, 1, 3, R, sign2=+1, l1=0, l2=1, n_steps=n_steps)
    progress(f"mixed s+p R={R} ++")
    elapsed = timer.time() - t0
    print(f"  {R:6.1f} {dE_pm:+10.4f} {dE_pp:+10.4f}  ({elapsed:.1f}s)", flush=True)

# --- DENSE R SCAN for oscillation detection ---
print(f"\n{'='*80}")
print(f"  DENSE R SCAN — Looking for oscillation")
print(f"{'='*80}", flush=True)

# s-wave (1,1) dense scan
print(f"\n  --- Dense scan: mode (1,1) s-wave ---")
print(f"  {'R':>6} {'dE(+-)':>10} {'dE(++)':>10}", flush=True)
for R in np.arange(1.0, 12.1, 0.5):
    _, _, _, dE_pm = run_3d_gpu(N, L, 1, 1, R, sign2=-1, l1=0, l2=0, n_steps=n_steps)
    progress(f"dense s(1,1) R={R}")
    _, _, _, dE_pp = run_3d_gpu(N, L, 1, 1, R, sign2=+1, l1=0, l2=0, n_steps=n_steps)
    progress(f"dense s(1,1) R={R}")
    print(f"  {R:6.1f} {dE_pm:+10.4f} {dE_pp:+10.4f}", flush=True)

# p-wave (3,3) dense scan
print(f"\n  --- Dense scan: mode (3,3) p-wave ---")
print(f"  {'R':>6} {'dE(+-)':>10} {'dE(++)':>10}", flush=True)
for R in np.arange(1.0, 12.1, 0.5):
    _, _, _, dE_pm = run_3d_gpu(N, L, 3, 3, R, sign2=-1, l1=1, l2=1, n_steps=n_steps)
    progress(f"dense p(3,3) R={R}")
    _, _, _, dE_pp = run_3d_gpu(N, L, 3, 3, R, sign2=+1, l1=1, l2=1, n_steps=n_steps)
    progress(f"dense p(3,3) R={R}")
    print(f"  {R:6.1f} {dE_pm:+10.4f} {dE_pp:+10.4f}", flush=True)

# Mixed s+p dense scan
print(f"\n  --- Dense scan: mode 1(s) + 3(p) mixed ---")
print(f"  {'R':>6} {'dE(+-)':>10} {'dE(++)':>10}", flush=True)
for R in np.arange(1.0, 12.1, 0.5):
    _, _, _, dE_pm = run_3d_gpu(N, L, 1, 3, R, sign2=-1, l1=0, l2=1, n_steps=n_steps)
    progress(f"dense mixed R={R}")
    _, _, _, dE_pp = run_3d_gpu(N, L, 1, 3, R, sign2=+1, l1=0, l2=1, n_steps=n_steps)
    progress(f"dense mixed R={R}")
    print(f"  {R:6.1f} {dE_pm:+10.4f} {dE_pp:+10.4f}", flush=True)

total_elapsed = timer.time() - start_time
print(f"\n{'='*80}")
print(f"  DONE — Total time: {total_elapsed/60:.1f} minutes")
print(f"{'='*80}")
