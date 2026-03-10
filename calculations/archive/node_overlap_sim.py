"""
GWT Node Overlap Simulation — GPU
====================================

Tests how radial nodes affect wave-wave interaction energy.

Key questions for V6 corrections 1 & 7:
1. Does E_int scale as |sin(phase)| / n_lobes when one wave has nodes?
2. Does ASYMMETRIC node mismatch (noded + nodeless) differ from
   SYMMETRIC (both noded or both nodeless)?
3. What is the actual suppression exponent? V6 uses:
   - Symmetric: S /= n_lobes
   - Asymmetric: S /= n_lobes^(4/3)
   Is 4/3 right, or should it be 1?

Molecules affected:
  Correction 1 (symmetric nodes): Li2(2s+2s), Na2(3s+3s)
  Correction 7 (asymmetric nodes): LiH(2s+1s), NaH(3s+1s), HCl(1s+3p), LiF(2s+2p), NaCl(3s+3p)

Simulation approach:
  - 3D sine-Gordon field on GPU (same physics as other sims)
  - Wave with radial nodes: envelope *= cos(n_nodes * pi * r / (3*width))
  - Compare: 1s+1s, 2s+1s, 2s+2s, 3s+1s, 3s+3s at various R
  - Extract E_int and compare to |sin(phase)| predictions
"""

import cupy as cp
import numpy as np
import time as timer

pi = float(np.pi)
k_long = 1.0
k_trans = 0.5


def run_interaction(N, L, R, wave1, wave2, n_steps=2500, dt_factor=0.12):
    """Compute interaction energy for two waves separated by R along z."""
    dx = L / N
    dt = dt_factor * dx
    coords = cp.linspace(-L/2, L/2, N, dtype=cp.float32)
    X, Y, Z = cp.meshgrid(coords, coords, coords, indexing='ij')

    def make_wave(wave, center_z):
        field = cp.zeros((*X.shape, 3), dtype=cp.float32)
        rz = Z - center_z
        r = cp.sqrt(X**2 + Y**2 + rz**2)
        r = cp.maximum(r, dx/2)

        n = wave['n']
        l = wave['l']
        width = wave['width']
        amp = wave['amplitude']

        envelope = amp * cp.exp(-r**2 / (2 * width**2))
        n_nodes = n - l - 1
        if n_nodes > 0:
            envelope *= cp.cos(n_nodes * pi * r / (3 * width))

        # Channel mapping: sigma(0)->z(2), px(1)->x(0), py(2)->y(1)
        ch_to_comp = {0: 2, 1: 0, 2: 1}
        channels = wave.get('channels', {0: 1.0})

        for ch_idx, ch_amp in channels.items():
            comp = ch_to_comp[ch_idx]
            if l == 0:
                field[:,:,:,comp] += ch_amp * envelope
            elif l == 1:
                if ch_idx == 0:
                    field[:,:,:,comp] += ch_amp * envelope * (rz / r)
                elif ch_idx == 1:
                    field[:,:,:,comp] += ch_amp * envelope * (X / r)
                elif ch_idx == 2:
                    field[:,:,:,comp] += ch_amp * envelope * (Y / r)
        return field

    def compute_energy(phi):
        E = cp.float32(0.0)
        V_0 = 1.0 / pi**2
        for axis in range(3):
            dphi = cp.diff(phi, axis=axis)
            E += 0.5 * k_long * cp.sum(dphi[:,:,:,axis]**2) * dx
            for comp in range(3):
                if comp != axis:
                    E += 0.5 * k_trans * cp.sum(dphi[:,:,:,comp]**2) * dx
        for comp in range(3):
            E += V_0 * cp.sum(1 - cp.cos(pi * phi[:,:,:,comp])) * dx**3
        return float(E)

    def compute_force(phi):
        force = cp.zeros_like(phi)
        V_0 = 1.0 / pi**2
        for axis in range(3):
            pp = cp.roll(phi, -1, axis=axis)
            pm = cp.roll(phi, 1, axis=axis)
            lap = (pp + pm - 2*phi) / dx**2
            force[:,:,:,axis] += k_long * lap[:,:,:,axis]
            for comp in range(3):
                if comp != axis:
                    force[:,:,:,comp] += k_trans * lap[:,:,:,comp]
        for comp in range(3):
            force[:,:,:,comp] -= V_0 * pi * cp.sin(pi * phi[:,:,:,comp])
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

    phi1 = make_wave(wave1, -R/2)
    E1 = evolve_and_measure(phi1)

    phi2 = make_wave(wave2, +R/2)
    E2 = evolve_and_measure(phi2)

    phi_both = make_wave(wave1, -R/2) + make_wave(wave2, +R/2)
    E_both = evolve_and_measure(phi_both)

    E_int = E_both - E1 - E2
    cp.cuda.Stream.null.synchronize()
    return E_int


# =============================================================================
print("=" * 80)
print("  NODE OVERLAP SIMULATION")
print("  How do radial nodes affect interaction energy?")
print("=" * 80)

_ = cp.zeros(10)
cp.cuda.Stream.null.synchronize()

N = 48
L = 16.0
n_steps = 2500
amp = 0.3  # moderate amplitude to stay in linear regime
width = 2.0


# =============================================================================
# TEST 1: 1s+1s baseline vs R
# No nodes — this is the reference for all comparisons
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 1: 1s+1s (no nodes) — baseline E_int vs R")
print(f"{'='*80}")

w_1s = {'amplitude': amp, 'width': width, 'n': 1, 'l': 0, 'channels': {0: 1.0}}

R_values = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0]

print(f"\n  {'R':>4}  {'E_int':>10}  {'|sin(2R)|':>10}  {'E/|sin|':>10}")

results_1s1s = {}
for R in R_values:
    t0 = timer.time()
    E_int = run_interaction(N, L, R, w_1s, w_1s, n_steps)
    dt = timer.time() - t0
    phase = 2*R  # R/n1 + R/n2 = R/1 + R/1 = 2R
    S = abs(np.sin(phase))
    ratio = E_int / S if abs(S) > 1e-4 else float('inf')
    results_1s1s[R] = (E_int, S, ratio)
    print(f"  {R:4.1f}  {E_int:10.4f}  {S:10.4f}  {ratio:10.4f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 2: 2s+1s (ASYMMETRIC nodes) vs R
# Models LiH (2s on Li, 1s on H)
# V6 says: S /= n_lobes^(4/3) for asymmetric
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 2: 2s+1s (asymmetric nodes) — E_int vs R")
print(f"  Models: LiH (Li 2s + H 1s)")
print(f"  V6 correction 7: S /= n_lobes^(4/3)")
print(f"{'='*80}")

w_2s = {'amplitude': amp, 'width': width, 'n': 2, 'l': 0, 'channels': {0: 1.0}}

# Phase for 2s+1s: R/n1^b1 + R/n2^b2
# With beta=0.95: b_2s = 1+0.95 = 1.95, b_1s = 1
# phase = R/1^1 + R/2^1.95 = R + R/3.81 ≈ 1.262*R
# But for raw simulation, n_2s=2 gives extra node in the wave

print(f"\n  {'R':>4}  {'E_2s1s':>10}  {'E_1s1s':>10}  {'ratio':>8}  {'|sin(ph)|':>10}  {'E/|sin|':>10}")

results_2s1s = {}
for R in R_values:
    t0 = timer.time()
    E_int = run_interaction(N, L, R, w_2s, w_1s, n_steps)
    dt = timer.time() - t0

    E_1s1s = results_1s1s[R][0]
    ratio = E_int / E_1s1s if abs(E_1s1s) > 1e-6 else 0

    # Phase for 2s+1s with sim parameters: R/1 + R/2 = 1.5*R
    # (no beta correction in sim — both have same width)
    phase = R + R/2  # = 1.5*R
    S = abs(np.sin(phase))
    E_over_S = E_int / S if abs(S) > 1e-4 else float('inf')

    results_2s1s[R] = (E_int, S, E_over_S, ratio)
    print(f"  {R:4.1f}  {E_int:10.4f}  {E_1s1s:10.4f}  {ratio:8.4f}  {S:10.4f}  {E_over_S:10.4f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 3: 2s+2s (SYMMETRIC nodes) vs R
# Models Li2 (2s on both)
# V6 correction 1: S /= n_lobes (divides by integer)
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 3: 2s+2s (symmetric nodes) — E_int vs R")
print(f"  Models: Li2 (2s + 2s)")
print(f"  V6 correction 1: S /= n_lobes")
print(f"{'='*80}")

print(f"\n  {'R':>4}  {'E_2s2s':>10}  {'E_1s1s':>10}  {'ratio':>8}  {'|sin(ph)|':>10}  {'E/|sin|':>10}")

results_2s2s = {}
for R in R_values:
    t0 = timer.time()
    E_int = run_interaction(N, L, R, w_2s, w_2s, n_steps)
    dt = timer.time() - t0

    E_1s1s = results_1s1s[R][0]
    ratio = E_int / E_1s1s if abs(E_1s1s) > 1e-6 else 0

    # Phase for 2s+2s: R/2 + R/2 = R
    phase = R
    S = abs(np.sin(phase))
    E_over_S = E_int / S if abs(S) > 1e-4 else float('inf')

    results_2s2s[R] = (E_int, S, E_over_S, ratio)
    print(f"  {R:4.1f}  {E_int:10.4f}  {E_1s1s:10.4f}  {ratio:8.4f}  {S:10.4f}  {E_over_S:10.4f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 4: 3s+1s (MORE nodes, asymmetric) vs R
# Models NaH (3s on Na, 1s on H)
# This should show stronger node suppression
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 4: 3s+1s (2 nodes, asymmetric) — E_int vs R")
print(f"  Models: NaH (Na 3s + H 1s)")
print(f"  3s wave has 2 radial nodes (n-l-1 = 3-0-1 = 2)")
print(f"{'='*80}")

w_3s = {'amplitude': amp, 'width': width, 'n': 3, 'l': 0, 'channels': {0: 1.0}}

print(f"\n  {'R':>4}  {'E_3s1s':>10}  {'E_1s1s':>10}  {'ratio':>8}  {'|sin(ph)|':>10}  {'E/|sin|':>10}")

results_3s1s = {}
for R in R_values:
    t0 = timer.time()
    E_int = run_interaction(N, L, R, w_3s, w_1s, n_steps)
    dt = timer.time() - t0

    E_1s1s = results_1s1s[R][0]
    ratio = E_int / E_1s1s if abs(E_1s1s) > 1e-6 else 0

    # Phase for 3s+1s: R/1 + R/3 = 4R/3
    phase = R + R/3
    S = abs(np.sin(phase))
    E_over_S = E_int / S if abs(S) > 1e-4 else float('inf')

    results_3s1s[R] = (E_int, S, E_over_S, ratio)
    print(f"  {R:4.1f}  {E_int:10.4f}  {E_1s1s:10.4f}  {ratio:8.4f}  {S:10.4f}  {E_over_S:10.4f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 5: 3s+3s (symmetric nodes) vs R
# Models Na2 (3s on both)
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 5: 3s+3s (symmetric, 2 nodes each) — E_int vs R")
print(f"  Models: Na2 (3s + 3s)")
print(f"{'='*80}")

print(f"\n  {'R':>4}  {'E_3s3s':>10}  {'E_1s1s':>10}  {'ratio':>8}  {'|sin(ph)|':>10}  {'E/|sin|':>10}")

results_3s3s = {}
for R in R_values:
    t0 = timer.time()
    E_int = run_interaction(N, L, R, w_3s, w_3s, n_steps)
    dt = timer.time() - t0

    E_1s1s = results_1s1s[R][0]
    ratio = E_int / E_1s1s if abs(E_1s1s) > 1e-6 else 0

    # Phase for 3s+3s: R/3 + R/3 = 2R/3
    phase = 2*R/3
    S = abs(np.sin(phase))
    E_over_S = E_int / S if abs(S) > 1e-4 else float('inf')

    results_3s3s[R] = (E_int, S, E_over_S, ratio)
    print(f"  {R:4.1f}  {E_int:10.4f}  {E_1s1s:10.4f}  {ratio:8.4f}  {S:10.4f}  {E_over_S:10.4f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 6: Direct node suppression measurement
# At a fixed R where |sin(phase)| is large, compare all combinations
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 6: Node suppression at R=4.0 — all combinations")
print(f"  Directly measures how nodes suppress E_int")
print(f"{'='*80}")

R_fix = 4.0
combos = [
    ('1s+1s', w_1s, w_1s, 0, 0),
    ('2s+1s', w_2s, w_1s, 1, 0),
    ('1s+2s', w_1s, w_2s, 0, 1),
    ('2s+2s', w_2s, w_2s, 1, 1),
    ('3s+1s', w_3s, w_1s, 2, 0),
    ('1s+3s', w_1s, w_3s, 0, 2),
    ('3s+3s', w_3s, w_3s, 2, 2),
    ('3s+2s', w_3s, w_2s, 2, 1),
    ('2s+3s', w_2s, w_3s, 1, 2),
]

print(f"\n  R = {R_fix}")
print(f"  {'combo':>7}  {'nodes_L':>7}  {'nodes_R':>7}  {'E_int':>10}  {'vs 1s1s':>8}  {'type':>10}")

E_ref = None
for label, w1, w2, n1, n2 in combos:
    t0 = timer.time()
    E_int = run_interaction(N, L, R_fix, w1, w2, n_steps)
    dt = timer.time() - t0

    if E_ref is None:
        E_ref = E_int

    ratio = E_int / E_ref if abs(E_ref) > 1e-6 else 0
    if n1 == n2:
        ntype = 'symmetric'
    elif n1 == 0 or n2 == 0:
        ntype = 'asymmetric'
    else:
        ntype = 'both_noded'

    print(f"  {label:>7}  {n1:>7}  {n2:>7}  {E_int:10.4f}  {ratio:8.4f}  {ntype:>10}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 7: Overlap floor test
# At R values where |sin(phase)| -> 0, does E_int also -> 0?
# Or does it bottom out at some minimum? (Tests correction 2)
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 7: Does E_int -> 0 when |sin(phase)| -> 0?")
print(f"  Tests whether overlap floor (correction 2) is real")
print(f"  Using 1s+1s at R values near sin zeros (R = k*pi/2)")
print(f"{'='*80}")

# sin(2R) = 0 when R = k*pi/2
# Test R near pi/2 ≈ 1.571, pi ≈ 3.14, 3pi/2 ≈ 4.71
R_fine = [1.0, 1.2, 1.4, 1.5, 1.571, 1.65, 1.8, 2.0,
          2.5, 3.0, 3.14, 3.3, 3.5, 4.0, 4.5, 4.71, 5.0, 5.5, 6.0]

print(f"\n  {'R':>5}  {'E_int':>10}  {'|sin(2R)|':>10}  {'E/|sin|':>10}  {'near_zero':>10}")

for R in R_fine:
    t0 = timer.time()
    E_int = run_interaction(N, L, R, w_1s, w_1s, n_steps)
    dt = timer.time() - t0

    phase = 2*R
    S = abs(np.sin(phase))
    E_over_S = E_int / S if abs(S) > 0.01 else float('inf')
    near_zero = '*' if S < 0.1 else ''

    print(f"  {R:5.3f}  {E_int:10.4f}  {S:10.4f}  {E_over_S:10.4f}  {near_zero:>10}  ({dt:.1f}s)", flush=True)


cp.cuda.Stream.null.synchronize()
print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
