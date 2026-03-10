"""
GWT sp Bond Energy Sweep — GPU
================================

Directly measures the covalent interaction energy for s(1s) + p(2p) bonds
as a function of separation R, and compares to V6 formula predictions.

Goal: identify whether the V6 D_cov formula (C_bond * E_scale * |sin(phase)|)
accurately captures the wave physics, or if there's a systematic deviation
for sp bonds that explains the H-X overshoot (OH +9.5%, NH +5.0%).

Key tests:
1. E_int(R) curve shape — does |sin(R/n1 + R/n2)| match?
2. Amplitude scaling — is C_bond*E_scale correct?
3. Does the p-wave angular structure (cos theta) modify the coupling?
4. Compare s+s vs s+p vs p+p interaction at same parameters
"""

import cupy as cp
import numpy as np
import time as timer

pi = float(np.pi)
k_long = 1.0
k_trans = 0.5

def run_bond_sim(N_grid, L, R_sep, wave1, wave2, n_steps=2000, dt_factor=0.12):
    """Run 3-component 3D bond simulation on GPU."""
    dx = L / N_grid
    dt = dt_factor * dx

    coords = cp.linspace(-L/2, L/2, N_grid, dtype=cp.float32)
    X, Y, Z = cp.meshgrid(coords, coords, coords, indexing='ij')

    def make_wave(X, Y, Z, wave, center_z):
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

        channels = wave.get('channels', {0: 1.0})
        ch_to_comp = {0: 2, 1: 0, 2: 1}

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
            if step % 50 == 0 and step > n_steps // 4:
                dphi_half = dphi - 0.5 * dt * force
                E_kin = 0.5 * float(cp.sum(dphi_half**2)) * dx**3
                E_pot = compute_energy(phi)
                E_samples.append(E_kin + E_pot)
        return np.mean(E_samples) if E_samples else compute_energy(phi_init)

    phi1 = make_wave(X, Y, Z, wave1, -R_sep/2)
    E1 = evolve_and_measure(phi1, n_steps)

    phi2 = make_wave(X, Y, Z, wave2, +R_sep/2)
    E2 = evolve_and_measure(phi2, n_steps)

    phi_both = make_wave(X, Y, Z, wave1, -R_sep/2) + make_wave(X, Y, Z, wave2, +R_sep/2)
    E_both = evolve_and_measure(phi_both, n_steps)

    E_int = E_both - E1 - E2
    cp.cuda.Stream.null.synchronize()
    return E1, E2, E_both, E_int


# =============================================================================
print("=" * 80)
print("  SP BOND ENERGY SWEEP")
print("  Comparing s+s, s+p, p+p interaction vs R")
print("=" * 80)

_ = cp.zeros(10)
cp.cuda.Stream.null.synchronize()

N = 48
L = 16.0
n_steps = 2500
width = 2.0
amp = 0.5

# =============================================================================
# TEST 1: s+p sigma bond vs R — the key curve
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 1: s(1s) + p(2p,sigma) bond energy vs R")
print(f"  V6 formula: D_cov = C_bond * E_scale * |sin(R/n1 + R/n2)|")
print(f"{'='*80}")

w_s = {'amplitude': amp, 'width': width, 'n': 1, 'l': 0, 'channels': {0: 1.0}}
w_p = {'amplitude': amp, 'width': width, 'n': 2, 'l': 1, 'channels': {0: 1.0}}

R_values = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0]

print(f"\n  {'R':>4}  {'E_int(sim)':>10}  {'|sin(ph)|':>10}  {'phase':>7}  {'ratio':>8}")

sp_results = []
for R in R_values:
    t0 = timer.time()
    _, _, _, E_int = run_bond_sim(N, L, R, w_s, w_p, n_steps)
    dt = timer.time() - t0

    # V6 phase: R/n1 + R/n2 (for 1s + 2p, no node corrections)
    phase = R/1.0 + R/2.0  # = 1.5*R
    S = abs(np.sin(phase))

    ratio = E_int / S if abs(S) > 1e-4 else float('inf')
    sp_results.append((R, E_int, phase, S, ratio))
    print(f"  {R:4.1f}  {E_int:10.4f}  {S:10.4f}  {phase:7.3f}  {ratio:8.4f}  ({dt:.1f}s)",
          flush=True)


# =============================================================================
# TEST 2: s+s bond vs R — baseline comparison
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 2: s(1s) + s(1s) bond energy vs R")
print(f"{'='*80}")

w_s2 = {'amplitude': amp, 'width': width, 'n': 1, 'l': 0, 'channels': {0: 1.0}}

print(f"\n  {'R':>4}  {'E_int(sim)':>10}  {'|sin(ph)|':>10}  {'phase':>7}  {'ratio':>8}")

ss_results = []
for R in R_values:
    t0 = timer.time()
    _, _, _, E_int = run_bond_sim(N, L, R, w_s2, w_s2, n_steps)
    dt = timer.time() - t0

    phase = R/1.0 + R/1.0  # = 2*R
    S = abs(np.sin(phase))

    ratio = E_int / S if abs(S) > 1e-4 else float('inf')
    ss_results.append((R, E_int, phase, S, ratio))
    print(f"  {R:4.1f}  {E_int:10.4f}  {S:10.4f}  {phase:7.3f}  {ratio:8.4f}  ({dt:.1f}s)",
          flush=True)


# =============================================================================
# TEST 3: p+p sigma bond vs R
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 3: p(2p,sigma) + p(2p,sigma) bond energy vs R")
print(f"{'='*80}")

print(f"\n  {'R':>4}  {'E_int(sim)':>10}  {'|sin(ph)|':>10}  {'phase':>7}  {'ratio':>8}")

pp_results = []
for R in R_values:
    t0 = timer.time()
    _, _, _, E_int = run_bond_sim(N, L, R, w_p, w_p, n_steps)
    dt = timer.time() - t0

    phase = R/2.0 + R/2.0  # = R
    S = abs(np.sin(phase))

    ratio = E_int / S if abs(S) > 1e-4 else float('inf')
    pp_results.append((R, E_int, phase, S, ratio))
    print(f"  {R:4.1f}  {E_int:10.4f}  {S:10.4f}  {phase:7.3f}  {ratio:8.4f}  ({dt:.1f}s)",
          flush=True)


# =============================================================================
# ANALYSIS
# =============================================================================
print(f"\n{'='*80}")
print(f"  ANALYSIS: Compare E_int/|sin(phase)| ratios")
print(f"  If V6 D_cov formula is correct, this ratio should be CONSTANT")
print(f"  (= C_bond * E_scale for each bond type)")
print(f"{'='*80}")

print(f"\n  {'R':>4}  {'ss ratio':>10}  {'sp ratio':>10}  {'pp ratio':>10}  {'sp/ss':>8}")
for i, R in enumerate(R_values):
    ss_r = ss_results[i][4]
    sp_r = sp_results[i][4]
    pp_r = pp_results[i][4]
    sp_ss = sp_r / ss_r if abs(ss_r) > 1e-4 else 0
    print(f"  {R:4.1f}  {ss_r:10.4f}  {sp_r:10.4f}  {pp_r:10.4f}  {sp_ss:8.4f}")


# =============================================================================
# TEST 4: Width dependence — does wave width affect sp vs ss ratio?
# The p-wave has angular structure (cos theta) that concentrates amplitude
# along the bond axis. Does this affect the coupling differently at small R?
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 4: sp/ss ratio at different wave widths")
print(f"  Does the p-wave angular structure modify coupling?")
print(f"{'='*80}")

R_test = 4.0
print(f"\n  R = {R_test}")
print(f"  {'width':>6}  {'E_ss':>10}  {'E_sp':>10}  {'sp/ss':>8}")

for w in [1.5, 2.0, 2.5, 3.0]:
    w_s_t = {'amplitude': amp, 'width': w, 'n': 1, 'l': 0, 'channels': {0: 1.0}}
    w_p_t = {'amplitude': amp, 'width': w, 'n': 2, 'l': 1, 'channels': {0: 1.0}}
    t0 = timer.time()
    _, _, _, E_ss = run_bond_sim(N, L, R_test, w_s_t, w_s_t, n_steps)
    _, _, _, E_sp = run_bond_sim(N, L, R_test, w_s_t, w_p_t, n_steps)
    dt = timer.time() - t0
    ratio = E_sp / E_ss if abs(E_ss) > 1e-6 else 0
    print(f"  {w:6.1f}  {E_ss:10.4f}  {E_sp:10.4f}  {ratio:8.4f}  ({dt:.1f}s)",
          flush=True)


print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
