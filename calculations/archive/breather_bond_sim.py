"""
Breather Bonding Simulation — Energy Minimization
===================================================
Properly compute two-breather interaction by MINIMIZING the total
lattice energy, not just superposing static profiles.

The bond energy = E(two breathers relaxed) - 2*E(single breather relaxed)

If E(coupled) < 2*E(single), the breathers ATTRACT (bonding).
If E(coupled) > 2*E(single), they REPEL (antibonding).

Uses scipy.optimize to find the minimum-energy field configuration.
"""

import numpy as np
from scipy.optimize import minimize
import time

pi = np.pi

# Lattice parameters (GWT)
V_0 = 1.0 / pi**2
N_br = 24
gamma = pi / (N_br + 1)


def breather_profile(x, n, center=0.0, amplitude=1.0):
    """Static breather profile at maximum displacement."""
    w_n = np.cos(n * gamma)
    eta = np.sin(n * gamma)
    eta = max(eta, 1e-10)
    return amplitude * (4.0/pi) * np.arctan(eta / np.cosh(eta * (x - center)))


def lattice_energy(phi, dx):
    """Total energy of field configuration."""
    grad = np.diff(phi)
    E_grad = 0.5 * np.sum(grad**2) / dx
    E_pot = V_0 * np.sum(1 - np.cos(pi * phi)) * dx
    return E_grad + E_pot


def lattice_energy_gradient(phi, dx):
    """Gradient of energy w.r.t. phi (for optimization)."""
    N = len(phi)
    grad = np.zeros(N)

    # Gradient energy contribution: d/dphi_i [(phi_i - phi_{i+1})^2 + (phi_{i-1} - phi_i)^2]
    # = 2*phi_i - phi_{i+1} - phi_{i-1}  (interior points)
    grad[1:-1] = (2*phi[1:-1] - phi[:-2] - phi[2:]) / dx
    grad[0] = (phi[0] - phi[1]) / dx
    grad[-1] = (phi[-1] - phi[-2]) / dx

    # Potential energy contribution
    grad += V_0 * pi * np.sin(pi * phi) * dx

    return grad


def minimize_energy(phi_init, dx, method='L-BFGS-B'):
    """Find minimum energy field configuration."""
    result = minimize(
        lambda phi: lattice_energy(phi, dx),
        phi_init,
        jac=lambda phi: lattice_energy_gradient(phi, dx),
        method=method,
        options={'maxiter': 5000, 'ftol': 1e-15, 'gtol': 1e-10}
    )
    return result.x, result.fun


# =============================================================================
# SIMULATION SETUP
# =============================================================================

N_lattice = 1000
L = 60.0  # total lattice length
x = np.linspace(-L/2, L/2, N_lattice)
dx = x[1] - x[0]


def single_breather_energy_minimized(n):
    """Minimized energy of a single breather."""
    phi_init = np.array([breather_profile(xi, n) for xi in x])
    phi_min, E_min = minimize_energy(phi_init, dx)
    return E_min, phi_min


def two_breather_interaction(n1, n2, R, signs=(1, 1)):
    """Compute interaction energy of two breathers at separation R.

    signs = (s1, s2): relative signs of the breather amplitudes.
    (+,+) = in-phase (bonding for bosonic breathers)
    (+,-) = out-of-phase (antibonding)
    """
    # Initial guess: superposition of two breather profiles
    phi_init = np.zeros(N_lattice)
    for i, xi in enumerate(x):
        phi_init[i] = (signs[0] * breather_profile(xi, n1, center=-R/2) +
                       signs[1] * breather_profile(xi, n2, center=+R/2))

    # Minimize
    phi_min, E_coupled = minimize_energy(phi_init, dx)

    return E_coupled, phi_min


# =============================================================================
# MAIN COMPUTATION
# =============================================================================

print("=" * 80)
print("  BREATHER BONDING SIMULATION")
print("  Energy minimization on sine-Gordon lattice")
print("=" * 80)
print(f"  Lattice: N={N_lattice}, L={L}, dx={dx:.4f}")
print(f"  V_0 = 1/pi^2 = {V_0:.6f}")

# Step 1: Single breather energies
print(f"\n--- Single breather energies (minimized) ---")
single_energies = {}
for n in range(1, 8):
    E, phi = single_breather_energy_minimized(n)
    w = 1.0 / np.sin(n * gamma)
    print(f"  mode {n}: E = {E:.8f}, width = {w:.2f}")
    single_energies[n] = E


# Step 2: Two-breather interaction vs R
print(f"\n{'='*80}")
print(f"  TWO-BREATHER INTERACTION ENERGY vs SEPARATION")
print(f"  E_int = E(coupled,minimized) - E_1 - E_2")
print(f"{'='*80}")

mode_pairs = [
    (1, 1, "H-H"),
    (1, 3, "H-2p (CH,BH,NH,HF,OH)"),
    (3, 3, "2p-2p (N2,O2,F2,B2)"),
    (1, 2, "H-Li (LiH)"),
    (2, 3, "Li-2p (LiF)"),
    (1, 4, "H-Na (NaH)"),
]

R_range = np.arange(1.0, 12.1, 0.5)

for n1, n2, label in mode_pairs:
    E1 = single_energies[n1]
    E2 = single_energies[n2]

    print(f"\n  --- ({n1},{n2}): {label} ---")
    print(f"  E_single: {E1:.6f} + {E2:.6f} = {E1+E2:.6f}")
    print(f"  {'R':>6} {'E_coupled':>12} {'E_int(++)':>12} {'E_int(+-)':>12} {'bonding':>8}")

    for R in R_range:
        # In-phase (bonding attempt)
        E_pp, _ = two_breather_interaction(n1, n2, R, signs=(+1, +1))
        E_int_pp = E_pp - E1 - E2

        # Out-of-phase (antibonding attempt)
        E_pm, _ = two_breather_interaction(n1, n2, R, signs=(+1, -1))
        E_int_pm = E_pm - E1 - E2

        # Bonding = lower energy
        bond = "BOND" if E_int_pp < -0.001 or E_int_pm < -0.001 else ""

        print(f"  {R:6.1f} {E_pp:12.6f} {E_int_pp:+12.6f} {E_int_pm:+12.6f} {bond:>8}")


# Step 3: Focus on key molecules
print(f"\n{'='*80}")
print(f"  COMPARISON TO REAL BOND ENERGIES")
print(f"  (Need to find the right energy scale conversion)")
print(f"{'='*80}")

molecules = [
    ('H2',   1.401,  4.745, 1, 1),
    ('CH',   2.116,  3.47,  1, 3),
    ('BH',   2.329,  3.42,  1, 3),
    ('NH',   1.958,  3.57,  1, 3),
    ('HF',   1.733,  5.869, 1, 3),
    ('OH',   1.834,  4.392, 1, 3),
    ('N2',   2.074,  9.759, 3, 3),
    ('LiH',  3.015,  2.515, 1, 2),
    ('NaH',  3.566,  1.97,  1, 4),
    ('Li2',  5.051,  1.056, 2, 2),
    ('LiF',  2.955,  5.939, 2, 3),
    ('BF',   2.386,  7.81,  3, 3),
]

print(f"\n{'Mol':<7} {'R':>6} {'De_exp':>7} {'E_int++':>10} {'E_int+-':>10} {'best':>10} {'scale':>8}")
print("-" * 65)

for name, R, De_exp, n1, n2 in molecules:
    E1 = single_energies[n1]
    E2 = single_energies[n2]

    E_pp, _ = two_breather_interaction(n1, n2, R, signs=(+1, +1))
    E_pm, _ = two_breather_interaction(n1, n2, R, signs=(+1, -1))

    E_int_pp = E_pp - E1 - E2
    E_int_pm = E_pm - E1 - E2

    # Use whichever gives lower energy (bonding)
    E_best = min(E_int_pp, E_int_pm)
    sign_best = "++" if E_int_pp < E_int_pm else "+-"

    # Scale factor needed to match De_exp
    scale = De_exp / abs(E_best) if abs(E_best) > 1e-8 else float('inf')

    print(f"{name:<7} {R:6.3f} {De_exp:7.3f} {E_int_pp:+10.4f} {E_int_pm:+10.4f} "
          f"{E_best:+10.4f} {scale:8.1f}")
