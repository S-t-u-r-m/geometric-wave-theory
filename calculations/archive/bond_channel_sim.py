"""
GWT Bond Channel Simulation — GPU
===================================

Extends breather_3d_gpu.py to test:
1. Does filling non-bonding channels on a p-wave atom reduce the sigma bond?
2. Does the sigma/pi interaction ratio match k/kappa = 2?
3. What correction naturally emerges from the wave physics?

Uses the EXISTING 3D sine-Gordon Hamiltonian with multi-component displacement.

Key change from breather_3d_gpu: displacement is now a 3-component VECTOR field
(sigma + 2 pi channels), with k_long for sigma coupling and k_trans = k/2 for pi.
"""

import cupy as cp
import numpy as np
import time as timer
import sys

pi = float(np.pi)

# Lattice coupling constants (from GWT Hamiltonian)
k_long = 1.0         # longitudinal = sigma
k_trans = 0.5         # transverse = pi = k/2 (isotropy)

def run_bond_sim(N_grid, L, R_sep, wave1, wave2, n_steps=2000, dt_factor=0.12):
    """Run 3-component 3D bond simulation on GPU.

    wave1, wave2: dicts with keys:
        'amplitude': float
        'width': float (Gaussian width in lattice units)
        'n': int (principal quantum number)
        'l': int (angular: 0=s, 1=p)
        'channels': dict {0: amp, 1: amp, 2: amp}
                    0=sigma(z), 1=pi_x, 2=pi_y
    """
    dx = L / N_grid
    dt = dt_factor * dx

    coords = cp.linspace(-L/2, L/2, N_grid, dtype=cp.float32)
    X, Y, Z = cp.meshgrid(coords, coords, coords, indexing='ij')

    def make_wave(X, Y, Z, wave, center_z):
        """Create a 3-component wave field."""
        # shape: (N, N, N, 3) for 3 displacement components
        field = cp.zeros((*X.shape, 3), dtype=cp.float32)

        rz = Z - center_z
        r = cp.sqrt(X**2 + Y**2 + rz**2)
        r = cp.maximum(r, dx/2)

        n = wave['n']
        l = wave['l']
        width = wave['width']
        amp = wave['amplitude']

        # Radial envelope
        envelope = amp * cp.exp(-r**2 / (2 * width**2))

        # Radial nodes
        n_nodes = n - l - 1
        if n_nodes > 0:
            envelope *= cp.cos(n_nodes * pi * r / (3 * width))

        channels = wave.get('channels', {0: 1.0})

        # Channel mapping: bond axis is Z, so
        #   channel 0 = sigma → component 2 (u_z, longitudinal to bond)
        #   channel 1 = pi_x  → component 0 (u_x, transverse)
        #   channel 2 = pi_y  → component 1 (u_y, transverse)
        ch_to_comp = {0: 2, 1: 0, 2: 1}

        for ch_idx, ch_amp in channels.items():
            comp = ch_to_comp[ch_idx]
            if l == 0:
                # s-wave: spherically symmetric
                field[:,:,:,comp] += ch_amp * envelope
            elif l == 1:
                if ch_idx == 0:
                    # p_z (sigma) = cos(theta) * envelope → u_z
                    field[:,:,:,comp] += ch_amp * envelope * (rz / r)
                elif ch_idx == 1:
                    # p_x (pi_x) = sin(theta)*cos(phi) * envelope → u_x
                    field[:,:,:,comp] += ch_amp * envelope * (X / r)
                elif ch_idx == 2:
                    # p_y (pi_y) = sin(theta)*sin(phi) * envelope → u_y
                    field[:,:,:,comp] += ch_amp * envelope * (Y / r)

        return field

    def compute_energy(phi):
        """Total energy of 3-component field with k_long/k_trans coupling."""
        E = cp.float32(0.0)

        # Gradient energy along each spatial axis
        for axis in range(3):
            dphi = cp.diff(phi, axis=axis)  # shape: (..., 3) with one axis shortened

            # Longitudinal: component matching the spatial axis
            E += 0.5 * k_long * cp.sum(dphi[:,:,:,axis]**2) * dx

            # Transverse: other components
            for comp in range(3):
                if comp != axis:
                    E += 0.5 * k_trans * cp.sum(dphi[:,:,:,comp]**2) * dx

        # On-site potential: sine-Gordon for each component
        V_0 = 1.0 / pi**2
        for comp in range(3):
            E += V_0 * cp.sum(1 - cp.cos(pi * phi[:,:,:,comp])) * dx**3

        return float(E)

    def compute_force(phi):
        """Force from 3-component Hamiltonian."""
        force = cp.zeros_like(phi)

        for axis in range(3):
            # Laplacian contribution along this axis
            pp = cp.roll(phi, -1, axis=axis)
            pm = cp.roll(phi, 1, axis=axis)
            lap = (pp + pm - 2*phi) / dx**2

            # Longitudinal component
            force[:,:,:,axis] += k_long * lap[:,:,:,axis]

            # Transverse components
            for comp in range(3):
                if comp != axis:
                    force[:,:,:,comp] += k_trans * lap[:,:,:,comp]

        # On-site force
        V_0 = 1.0 / pi**2
        for comp in range(3):
            force[:,:,:,comp] -= V_0 * pi * cp.sin(pi * phi[:,:,:,comp])

        return force

    def evolve_and_measure(phi_init, n_steps):
        """Evolve with velocity Verlet and measure average energy."""
        phi = phi_init.copy()
        dphi = cp.zeros_like(phi)

        # Initial half-step
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

    # --- Compute interaction energy ---

    # Single wave 1
    phi1 = make_wave(X, Y, Z, wave1, -R_sep/2)
    E1 = evolve_and_measure(phi1, n_steps)

    # Single wave 2
    phi2 = make_wave(X, Y, Z, wave2, +R_sep/2)
    E2 = evolve_and_measure(phi2, n_steps)

    # Combined system
    phi_both = make_wave(X, Y, Z, wave1, -R_sep/2) + make_wave(X, Y, Z, wave2, +R_sep/2)
    E_both = evolve_and_measure(phi_both, n_steps)

    E_int = E_both - E1 - E2

    cp.cuda.Stream.null.synchronize()
    return E1, E2, E_both, E_int


# =============================================================================
# MAIN SIMULATION
# =============================================================================
print("=" * 80)
print("  GWT BOND CHANNEL SIMULATION — GPU")
print("  3-component waves with k_long/k_trans coupling")
print(f"  k_long = {k_long}, k_trans = {k_trans}, ratio = {k_trans/k_long}")
print("=" * 80)

# Warm up GPU
_ = cp.zeros(10)
cp.cuda.Stream.null.synchronize()

N = 48
L = 16.0
n_steps = 2000
width = 2.0
amp = 0.5

print(f"  Grid: {N}^3 = {N**3}, Box: L={L}, Steps: {n_steps}")

# =============================================================================
# TEST 1: Sigma vs Pi bond energy ratio
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 1: Sigma vs Pi interaction energy ratio")
print(f"  Pure sigma bond vs pure pi bond at same R")
print(f"  Expected ratio: k_long/k_trans = {k_long/k_trans:.1f}")
print(f"{'='*80}")

s_wave = {'amplitude': amp, 'width': width, 'n': 1, 'l': 0}

print(f"\n  {'R':>4}  {'E_sigma':>10}  {'E_pi':>10}  {'ratio':>8}")
for R in [3.0, 4.0, 5.0, 6.0, 8.0]:
    t0 = timer.time()

    # Pure sigma bond: both waves in channel 0
    w_sig = {**s_wave, 'channels': {0: 1.0}}
    _, _, _, E_sig = run_bond_sim(N, L, R, w_sig, w_sig, n_steps)

    # Pure pi bond: both waves in channel 1
    w_pi = {**s_wave, 'channels': {1: 1.0}}
    _, _, _, E_pi = run_bond_sim(N, L, R, w_pi, w_pi, n_steps)

    ratio = E_sig / E_pi if abs(E_pi) > 1e-10 else float('inf')
    dt = timer.time() - t0
    print(f"  {R:4.1f}  {E_sig:10.4f}  {E_pi:10.4f}  {ratio:8.3f}  ({dt:.1f}s)",
          flush=True)


# =============================================================================
# TEST 2: Does non-bonding channel occupancy affect bond energy?
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 2: Non-bonding channel effect on sigma bond")
print(f"  H(1s,sigma) + X(2p,varying channels)")
print(f"  If lone-pair screening exists, bond energy should decrease")
print(f"  as more non-bonding pi channels are filled")
print(f"{'='*80}")

R_bond = 4.0
w_H = {'amplitude': amp, 'width': width, 'n': 1, 'l': 0, 'channels': {0: 1.0}}

configs = [
    ("B-like: sig only",     {'amplitude': amp, 'width': width, 'n': 2, 'l': 1,
                               'channels': {0: 1.0}}),
    ("C-like: sig + 1pi",    {'amplitude': amp, 'width': width, 'n': 2, 'l': 1,
                               'channels': {0: 0.8, 1: 0.6}}),
    ("N-like: sig + 2pi",    {'amplitude': amp, 'width': width, 'n': 2, 'l': 1,
                               'channels': {0: 0.7, 1: 0.5, 2: 0.5}}),
    ("O-like: sig + 2pi(f)", {'amplitude': amp, 'width': width, 'n': 2, 'l': 1,
                               'channels': {0: 0.6, 1: 0.7, 2: 0.5}}),
    ("F-like: sig + 2pi(ff)",{'amplitude': amp, 'width': width, 'n': 2, 'l': 1,
                               'channels': {0: 0.5, 1: 0.7, 2: 0.7}}),
]

print(f"\n  R = {R_bond}")
print(f"  {'Config':22s}  {'E_int':>10s}  {'vs_ref':>8s}  {'channels':>20s}")

E_ref = None
for label, w_X in configs:
    t0 = timer.time()
    _, _, _, E_int = run_bond_sim(N, L, R_bond, w_H, w_X, n_steps)
    dt = timer.time() - t0

    if E_ref is None:
        E_ref = E_int
        vs = "  (ref)"
    else:
        pct = (E_int - E_ref) / abs(E_ref) * 100 if abs(E_ref) > 1e-10 else 0
        vs = f"{pct:+7.1f}%"

    ch_str = str(w_X['channels'])
    print(f"  {label:22s}  {E_int:10.4f}  {vs:>8s}  {ch_str:>20s}  ({dt:.1f}s)",
          flush=True)


# =============================================================================
# TEST 3: Same test but keep sigma amplitude FIXED
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 3: Fixed sigma amplitude, add pi channels")
print(f"  This isolates the screening effect from amplitude dilution")
print(f"{'='*80}")

configs2 = [
    ("sig=1.0, pi=0",     {0: 1.0}),
    ("sig=1.0, pi_x=0.5", {0: 1.0, 1: 0.5}),
    ("sig=1.0, pi=0.5",   {0: 1.0, 1: 0.5, 2: 0.5}),
    ("sig=1.0, pi=0.7",   {0: 1.0, 1: 0.7, 2: 0.7}),
    ("sig=1.0, pi=1.0",   {0: 1.0, 1: 1.0, 2: 1.0}),
]

print(f"\n  R = {R_bond}, sigma amplitude fixed at 1.0")
print(f"  {'Config':22s}  {'E_int':>10s}  {'vs_ref':>8s}")

E_ref2 = None
for label, channels in configs2:
    w_X = {'amplitude': amp, 'width': width, 'n': 2, 'l': 1, 'channels': channels}
    t0 = timer.time()
    _, _, _, E_int = run_bond_sim(N, L, R_bond, w_H, w_X, n_steps)
    dt = timer.time() - t0

    if E_ref2 is None:
        E_ref2 = E_int
        vs = "  (ref)"
    else:
        pct = (E_int - E_ref2) / abs(E_ref2) * 100 if abs(E_ref2) > 1e-10 else 0
        vs = f"{pct:+7.1f}%"

    print(f"  {label:22s}  {E_int:10.4f}  {vs:>8s}  ({dt:.1f}s)", flush=True)


# =============================================================================
# TEST 4: pp bond — sigma vs pi contributions
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 4: p+p bond decomposition")
print(f"  Separate sigma and pi channel contributions")
print(f"{'='*80}")

pp_configs = [
    ("sigma only",       {0: 1.0},            {0: 1.0}),
    ("pi_x only",        {1: 1.0},            {1: 1.0}),
    ("2 pi",             {1: 0.7, 2: 0.7},    {1: 0.7, 2: 0.7}),
    ("sig + 2pi (N2)",   {0: 0.6, 1: 0.6, 2: 0.6}, {0: 0.6, 1: 0.6, 2: 0.6}),
    ("sig+2pi+pi*(O2)",  {0: 0.6, 1: 0.8, 2: 0.6}, {0: 0.6, 1: 0.8, 2: 0.6}),
]

print(f"\n  R = {R_bond}, p-wave (n=2, l=1)")
print(f"  {'Config':22s}  {'E_int':>10s}")

for label, ch1, ch2 in pp_configs:
    w1 = {'amplitude': amp, 'width': width, 'n': 2, 'l': 1, 'channels': ch1}
    w2 = {'amplitude': amp, 'width': width, 'n': 2, 'l': 1, 'channels': ch2}
    t0 = timer.time()
    _, _, _, E_int = run_bond_sim(N, L, R_bond, w1, w2, n_steps)
    dt = timer.time() - t0
    print(f"  {label:22s}  {E_int:10.4f}  ({dt:.1f}s)", flush=True)


# =============================================================================
# SUMMARY
# =============================================================================
total = timer.time() - timer.time()  # placeholder
print(f"\n{'='*80}")
print(f"  SIMULATION COMPLETE")
print(f"{'='*80}")
print(f"""
  KEY QUESTIONS ANSWERED:
  1. Does sigma/pi ratio = k_long/k_trans = 2?
  2. Does non-bonding channel occupancy screen the sigma bond?
  3. Is the screening proportional to number of filled channels?
  4. What correction form matches the simulation?
""")
