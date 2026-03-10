"""
GWT Bond Channel Simulation Part 2 — Pi Repulsion Test
========================================================

Test 4 showed p-wave pi channels are REPULSIVE when interacting with another wave.
This test measures: for s+p bonds (like H-X), how much do non-bonding pi channels
on the p-atom ADD repulsive energy?

If this is the source of the OH/NH/CH overshoot, then:
- More non-bonding pi electrons → more repulsion → lower D_e
- The V6 formula ignores this → overpredicts
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
        ch_to_comp = {0: 2, 1: 0, 2: 1}  # sigma→z, pi_x→x, pi_y→y

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
print("  TEST A: s(1s) + p(2p) with FIXED sigma, varying pi channels")
print("  This measures pi REPULSION in sp bonds")
print("=" * 80)

_ = cp.zeros(10)
cp.cuda.Stream.null.synchronize()

N = 48
L = 16.0
n_steps = 2500
width = 2.0
amp = 0.5

w_H = {'amplitude': amp, 'width': width, 'n': 1, 'l': 0, 'channels': {0: 1.0}}

# Key test: p-wave with sigma=1.0 FIXED, add pi channels
# If pi adds repulsion, E_int should become less negative (weaker bond)
print(f"\n  R=4, H(1s,sigma) + X(2p), sigma FIXED at 1.0, adding pi channels")
print(f"  {'Config':25s}  {'E_int':>10s}  {'change':>8s}")

E_ref = None
for label, channels in [
    ("sig=1.0 only",          {0: 1.0}),
    ("sig=1.0 + pi_x=0.3",   {0: 1.0, 1: 0.3}),
    ("sig=1.0 + pi_x=0.5",   {0: 1.0, 1: 0.5}),
    ("sig=1.0 + 2pi=0.3",    {0: 1.0, 1: 0.3, 2: 0.3}),
    ("sig=1.0 + 2pi=0.5",    {0: 1.0, 1: 0.5, 2: 0.5}),
    ("sig=1.0 + 2pi=0.7",    {0: 1.0, 1: 0.7, 2: 0.7}),
    ("sig=1.0 + 2pi=1.0",    {0: 1.0, 1: 1.0, 2: 1.0}),
]:
    w_X = {'amplitude': amp, 'width': width, 'n': 2, 'l': 1, 'channels': channels}
    t0 = timer.time()
    _, _, _, E_int = run_bond_sim(N, L, 4.0, w_H, w_X, n_steps)
    dt = timer.time() - t0

    if E_ref is None:
        E_ref = E_int
        vs = "  (ref)"
    else:
        pct = (E_int - E_ref) / abs(E_ref) * 100 if abs(E_ref) > 1e-10 else 0
        vs = f"{pct:+7.1f}%"

    print(f"  {label:25s}  {E_int:10.4f}  {vs}  ({dt:.1f}s)", flush=True)


# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST B: Same test at multiple distances")
print(f"  Does pi repulsion scale differently with R than sigma bonding?")
print(f"{'='*80}")

print(f"\n  {'R':>4s}  {'sig_only':>10s}  {'sig+2pi':>10s}  {'delta':>10s}  {'delta%':>8s}")

for R in [3.0, 4.0, 5.0, 6.0, 8.0]:
    w_sig = {'amplitude': amp, 'width': width, 'n': 2, 'l': 1, 'channels': {0: 1.0}}
    w_full = {'amplitude': amp, 'width': width, 'n': 2, 'l': 1,
              'channels': {0: 1.0, 1: 0.7, 2: 0.7}}

    t0 = timer.time()
    _, _, _, E_sig = run_bond_sim(N, L, R, w_H, w_sig, n_steps)
    _, _, _, E_full = run_bond_sim(N, L, R, w_H, w_full, n_steps)
    dt = timer.time() - t0

    delta = E_full - E_sig
    pct = delta / abs(E_sig) * 100 if abs(E_sig) > 1e-10 else 0
    print(f"  {R:4.1f}  {E_sig:10.4f}  {E_full:10.4f}  {delta:10.4f}  {pct:+7.1f}%  ({dt:.1f}s)",
          flush=True)


# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST C: p+p bond — decompose all channel combinations")
print(f"  Which channels bond, which repel?")
print(f"{'='*80}")

print(f"\n  R=4, p(2p) + p(2p)")
print(f"  {'Atom1':>12s}  {'Atom2':>12s}  {'E_int':>10s}  {'nature':>10s}")

combos = [
    ("sig", "sig",   {0: 1.0}, {0: 1.0}),
    ("sig", "pi_x",  {0: 1.0}, {1: 1.0}),
    ("sig", "pi_y",  {0: 1.0}, {2: 1.0}),
    ("pi_x", "pi_x", {1: 1.0}, {1: 1.0}),
    ("pi_x", "pi_y", {1: 1.0}, {2: 1.0}),
    ("pi_y", "pi_y", {2: 1.0}, {2: 1.0}),
]

for l1, l2, ch1, ch2 in combos:
    w1 = {'amplitude': amp, 'width': width, 'n': 2, 'l': 1, 'channels': ch1}
    w2 = {'amplitude': amp, 'width': width, 'n': 2, 'l': 1, 'channels': ch2}
    t0 = timer.time()
    _, _, _, E_int = run_bond_sim(N, L, 4.0, w1, w2, n_steps)
    dt = timer.time() - t0

    nature = "BONDING" if E_int < -0.01 else ("REPULSIVE" if E_int > 0.01 else "~ZERO")
    print(f"  {l1:>12s}  {l2:>12s}  {E_int:10.4f}  {nature:>10s}  ({dt:.1f}s)", flush=True)


# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST D: s(1s) + p(2p) — cross-channel interactions")
print(f"  Does H(sigma) interact with X(pi)?")
print(f"{'='*80}")

print(f"\n  R=4, H(1s) channel vs X(2p) channel")
print(f"  {'H_ch':>6s}  {'X_ch':>6s}  {'E_int':>10s}  {'nature':>10s}")

for h_ch, x_ch, h_channels, x_channels in [
    ("sig", "sig",  {0: 1.0}, {0: 1.0}),
    ("sig", "pi_x", {0: 1.0}, {1: 1.0}),
    ("sig", "pi_y", {0: 1.0}, {2: 1.0}),
]:
    w1 = {'amplitude': amp, 'width': width, 'n': 1, 'l': 0, 'channels': h_channels}
    w2 = {'amplitude': amp, 'width': width, 'n': 2, 'l': 1, 'channels': x_channels}
    t0 = timer.time()
    _, _, _, E_int = run_bond_sim(N, L, 4.0, w1, w2, n_steps)
    dt = timer.time() - t0

    nature = "BONDING" if E_int < -0.01 else ("REPULSIVE" if E_int > 0.01 else "~ZERO")
    print(f"  {h_ch:>6s}  {x_ch:>6s}  {E_int:10.4f}  {nature:>10s}  ({dt:.1f}s)", flush=True)


print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
