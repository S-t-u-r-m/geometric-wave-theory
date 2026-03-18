"""
Direct Screening Measurement from Sine-Gordon Dynamics
=======================================================
For each atom: add core breather modes one by one to the kink.
Measure how each mode changes the valence mode's binding energy.
The change IS the screening — measured, not assumed.

Compare to v19's screening weights to see where the formula
is right and where it's approximate.
"""

import numpy as np
from math import factorial
import time

d = 3
gamma_sg = np.pi / (2**(d+1) * np.pi - 2)
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6
V_0 = 1.0 / np.pi**2
w_pi = np.cos(np.pi / d)

# Grid
Nx = 1500
x_max = 25.0
dx = 2 * x_max / Nx
x = np.linspace(-x_max, x_max, Nx)
dt = 0.4 * dx
N_steps = 12000

_Z_well = 1.0

def set_Z(Z):
    global _Z_well
    _Z_well = Z

def kink(x_grid, Z):
    return (4.0 / np.pi) * np.arctan(np.exp(np.sqrt(Z) * x_grid))

def add_breather(phi_bg, x_grid, omega, x0=0.0, amplitude=1.0):
    """Add breather mode at position x0 with frequency omega."""
    eps = np.sqrt(max(1.0 - omega**2, 1e-12))
    num = eps * amplitude
    den = omega * np.cosh(eps * (x_grid - x0))
    return phi_bg + (4.0 / np.pi) * np.arctan(num / (den + 1e-30))

def evolve(phi_init, n_steps=None):
    """Evolve and return time-averaged energy."""
    if n_steps is None:
        n_steps = N_steps
    phi = phi_init.copy()
    phi_old = phi.copy()
    energies = []
    for step in range(n_steps):
        lap = np.zeros_like(phi)
        lap[1:-1] = (phi[2:] - 2*phi[1:-1] + phi[:-2]) / dx**2
        force = _Z_well * (1.0 / np.pi) * np.sin(np.pi * phi)
        phi_new = 2*phi - phi_old + dt**2 * (lap - force)
        phi_new[0] = phi_init[0]
        phi_new[-1] = phi_init[-1]
        phi_old = phi.copy()
        phi = phi_new
        if step > n_steps // 2 and step % 50 == 0:
            KE = 0.5 * np.sum(((phi - phi_old)/dt)**2) * dx
            GE = 0.5 * np.sum(((phi[1:] - phi[:-1])/dx)**2) * dx
            PE = np.sum(_Z_well * V_0 * (1 - np.cos(np.pi * phi))) * dx
            energies.append(KE + GE + PE)
    return np.mean(energies) if energies else 0.0


def measure_screening(Z, core_modes, val_mode):
    """Measure screening: how core modes change the valence binding.

    core_modes: list of (omega, x_position, label)
    val_mode: (omega, x_position, label) for the valence mode

    Returns:
        E_bare: valence binding without core (in bare kink potential)
        E_screened: valence binding with core present
        S_eff: effective screening (Z_net = Z - S_eff)
    """
    set_Z(Z)
    phi_k = kink(x, Z)

    # Step 1: kink alone
    E_k = evolve(phi_k)

    # Step 2: kink + valence only
    omega_v, x_v, _ = val_mode
    phi_kv = add_breather(phi_k, x, omega_v, x0=x_v)
    E_kv = evolve(phi_kv)
    E_bare = E_kv - E_k  # valence binding in bare potential

    # Step 3: kink + all core modes
    phi_kc = phi_k.copy()
    for omega_c, x_c, _ in core_modes:
        phi_kc = add_breather(phi_kc, x, omega_c, x0=x_c)
    E_kc = evolve(phi_kc)

    # Step 4: kink + core + valence
    phi_kcv = add_breather(phi_kc, x, omega_v, x0=x_v)
    E_kcv = evolve(phi_kcv)
    E_screened = E_kcv - E_kc  # valence binding with core screening

    # Effective screening: how much Z_net is reduced
    # E ~ Z_net^(2*alpha), so E_screened/E_bare = (Z_net/Z)^(2*alpha)
    # For small screening: S_eff ≈ Z * (1 - E_screened/E_bare)
    ratio = E_screened / E_bare if abs(E_bare) > 1e-10 else 1.0

    return E_bare, E_screened, ratio


# === TEST ATOMS ===
print("Direct Screening Measurement from Sine-Gordon Dynamics")
print("=" * 65)
print(f"Grid: {Nx}, N_steps: {N_steps}")
print()

# Breather frequencies for different angular channels:
# l=0 (s): omega_s = cos(1*gamma)   [fundamental]
# l=1 (p): omega_p = cos(2*gamma)   [2nd harmonic]
# l=2 (d): omega_d = cos(3*gamma)   [3rd harmonic]
omega_s = np.cos(1 * gamma_sg)
omega_p = np.cos(2 * gamma_sg)
omega_d = np.cos(3 * gamma_sg)

print(f"Mode frequencies: omega_s={omega_s:.6f}, omega_p={omega_p:.6f}, omega_d={omega_d:.6f}")
print()

# Test 1: He-like (Z=2, 1s2)
# Core: one s-mode. Valence: second s-mode (same channel pairing)
print("--- He-like: Z=2, core=1s, val=1s (s-pair) ---")
t0 = time.time()
E_bare, E_scr, ratio = measure_screening(
    Z=2,
    core_modes=[(omega_s, 0.0, '1s')],
    val_mode=(omega_s, 0.0, '1s_val')
)
print(f"  E_bare={E_bare:.4f}  E_screened={E_scr:.4f}  ratio={ratio:.4f}")
print(f"  Screening reduces binding by {(1-ratio)*100:.1f}%")
print(f"  v19: S_core = 0.5 (one s-channel at w_pi)")
print(f"  Time: {time.time()-t0:.1f}s")
print(flush=True)

# Test 2: Li-like (Z=3, 1s2 + 2s1)
# Core: two s-modes (1s pair). Valence: s-mode at higher energy (2s)
# Use lower omega for 2s (less tightly bound)
omega_2s = np.cos(1 * gamma_sg) * 0.95  # slightly lower freq for outer shell
print("\n--- Li-like: Z=3, core=1s2, val=2s ---")
t0 = time.time()
E_bare, E_scr, ratio = measure_screening(
    Z=3,
    core_modes=[(omega_s, 0.0, '1s'), (omega_s, 0.0, '1s')],
    val_mode=(omega_2s, 0.0, '2s')
)
print(f"  E_bare={E_bare:.4f}  E_screened={E_scr:.4f}  ratio={ratio:.4f}")
print(f"  Screening reduces binding by {(1-ratio)*100:.1f}%")
print(f"  v19: S_core = 0.5 (one s-channel at w_pi)")
print(f"  Time: {time.time()-t0:.1f}s")
print(flush=True)

# Test 3: B-like (Z=5, 1s2 2s2 + 2p1)
# Core: two s-pairs (1s2 + 2s2). Valence: p-mode
print("\n--- B-like: Z=5, core=1s2+2s2, val=2p ---")
t0 = time.time()
E_bare, E_scr, ratio = measure_screening(
    Z=5,
    core_modes=[(omega_s, 0.0, '1s'), (omega_s, 0.0, '1s'),
                (omega_2s, 0.0, '2s'), (omega_2s, 0.0, '2s')],
    val_mode=(omega_p, 0.0, '2p')
)
print(f"  E_bare={E_bare:.4f}  E_screened={E_scr:.4f}  ratio={ratio:.4f}")
print(f"  Screening reduces binding by {(1-ratio)*100:.1f}%")
print(f"  v19: S_core = 2.5 (1s: 0.5 + 2s: 0.5 + 2p_self: 1.5)")
print(f"  Time: {time.time()-t0:.1f}s")
print(flush=True)

# Test 4: Incremental screening — add modes one at a time
print("\n--- Incremental screening at Z=10 (Ne-like) ---")
print("Adding modes one at a time, measuring cumulative screening")
print()

modes_to_add = [
    (omega_s, 0.0, '1s_1'),
    (omega_s, 0.0, '1s_2'),
    (omega_2s, 0.0, '2s_1'),
    (omega_2s, 0.0, '2s_2'),
    (omega_p, -0.5, '2p_1'),
    (omega_p, 0.0, '2p_2'),
    (omega_p, 0.5, '2p_3'),
]

# Valence: outermost p-mode
val = (omega_p, 1.0, '2p_val')

set_Z(10)
phi_k = kink(x, 10)
E_k = evolve(phi_k)

# Bare valence binding
phi_kv = add_breather(phi_k, x, val[0], x0=val[1])
E_kv = evolve(phi_kv)
E_bare_val = E_kv - E_k

print(f"  {'Core':>15} {'E_val_bind':>10} {'ratio':>8} {'delta_S':>8}")
print(f"  {'-'*45}")
print(f"  {'(bare)':>15} {E_bare_val:10.4f} {1.000:8.4f} {'---':>8}")

phi_core = phi_k.copy()
E_core = E_k
prev_ratio = 1.0

for omega_c, x_c, label in modes_to_add:
    phi_core = add_breather(phi_core, x, omega_c, x0=x_c)
    E_core = evolve(phi_core)

    phi_cv = add_breather(phi_core, x, val[0], x0=val[1])
    E_cv = evolve(phi_cv)
    E_val = E_cv - E_core

    ratio = E_val / E_bare_val
    delta = ratio - prev_ratio

    print(f"  +{label:>14} {E_val:10.4f} {ratio:8.4f} {delta:+8.4f}", flush=True)
    prev_ratio = ratio

print()
print("Each delta shows how much that mode changed the valence binding.")
print("This IS the per-mode screening — directly from wave dynamics.")
