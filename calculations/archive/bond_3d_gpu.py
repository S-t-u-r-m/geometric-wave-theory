#!/usr/bin/env python3
"""
3D Bond Simulation — Two Breathers on the d=3 Cubic Lattice (GPU)
==================================================================
The REAL simulation. No 1D approximations. All forces from ONE Lagrangian:
  L = sum_<ij> [(1/2)(phi_i - phi_j)^2] + sum_i [(1/pi^2)(1 - cos(pi*phi_i))]

Two atoms = two kinks (nuclear cores) + breather modes (valence electrons).
Place them at separation R along z-axis. Evolve. Measure energy.

In 3D on the cubic lattice:
- Oh symmetry emerges naturally (it's a cube!)
- Angular channels (sigma, pi, delta) are REAL geometric projections
- LP repulsion comes from the 3D field trying to go over the barrier
- The cosine potential + 3D gradients give the COMPLETE interaction

The Oh lookup table precomputes WHICH modes to place on each atom.
The GPU computes HOW they interact — no formulas needed.

Uses CuPy (CUDA) on RTX 4070 Ti.
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial
import time

try:
    import cupy as cp
    GPU = True
    xp = cp
    print("GPU: CuPy + CUDA active")
except ImportError:
    cp = np
    xp = np
    GPU = False
    print("GPU: CuPy not available, using NumPy (will be slow)")

# ============================================================
# CONSTANTS FROM d=3
# ============================================================
d = 3
V_0 = 1.0 / np.pi**2
gamma_sg = np.pi / (2**(d+1) * np.pi - 2)
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6  # 13.604 eV

# Breather frequencies for each angular mode
omega = [float(np.cos(n * gamma_sg)) for n in range(1, 5)]
# omega[0] = s-mode, omega[1] = p-mode, omega[2] = d-mode, omega[3] = f-mode

# Amplitudes by occupancy
A_HALF = 0.10   # single electron: small perturbation on the kink
A_FULL = 0.20   # paired electrons: double amplitude


# ============================================================
# 3D GRID (cubic lattice)
# ============================================================
N = 64          # 64^3 grid
BOX = 10.0      # box size in Bohr (±BOX)
dx = 2 * BOX / N
dt = 0.25 * dx  # CFL condition: dt < dx/sqrt(3) ≈ 0.58*dx

print(f"Grid: {N}^3 = {N**3:,} points, box = ±{BOX} Bohr, dx = {dx:.3f}")
print()

# 1D coordinate arrays
x1d = np.linspace(-BOX, BOX, N, endpoint=False)


def make_grid(center_z=0.0):
    """Create 3D coordinate arrays centered on a point along z-axis."""
    X, Y, Z = np.meshgrid(x1d, x1d, x1d - center_z, indexing='ij')
    R = np.sqrt(X**2 + Y**2 + Z**2) + 1e-10
    return X, Y, Z, R


# ============================================================
# FIELD OPERATIONS (GPU)
# ============================================================
def laplacian_3d(phi):
    """6-point stencil Laplacian on periodic cubic lattice."""
    return (
        xp.roll(phi, 1, 0) + xp.roll(phi, -1, 0) +
        xp.roll(phi, 1, 1) + xp.roll(phi, -1, 1) +
        xp.roll(phi, 1, 2) + xp.roll(phi, -1, 2) -
        6 * phi
    ) / dx**2


def total_energy(phi):
    """Total Lagrangian energy of 3D field."""
    # Gradient energy: (1/2) sum of squared differences
    GE = 0.5 * (
        (xp.roll(phi, 1, 0) - phi)**2 +
        (xp.roll(phi, 1, 1) - phi)**2 +
        (xp.roll(phi, 1, 2) - phi)**2
    )
    # Potential energy: (1/pi^2)(1 - cos(pi*phi))
    PE = V_0 * (1.0 - xp.cos(np.pi * phi))
    return float(xp.sum(GE + PE) * dx**3)


def evolve(phi_init, n_steps=3000, damping=0.01):
    """
    Evolve 3D sine-Gordon with light damping → find equilibrium.

    Damped wave equation:
    phi_tt = laplacian(phi) - (1/pi)*sin(pi*phi) - damping * phi_t

    Light damping lets the system ring down to its ground state
    without artificially constraining the field.
    """
    phi = xp.asarray(phi_init.copy(), dtype=np.float64)
    phi_old = phi.copy()

    energies = []
    for step in range(n_steps):
        lap = laplacian_3d(phi)
        force = (1.0 / np.pi) * xp.sin(np.pi * phi)
        phi_new = (2 - damping*dt) * phi - (1 - damping*dt) * phi_old + dt**2 * (lap - force)
        phi_old = phi.copy()
        phi = phi_new

        # Sample energy in second half (after transient)
        if step > n_steps * 2 // 3 and step % 50 == 0:
            E = total_energy(phi)
            energies.append(E)

    E_avg = np.mean(energies) if energies else total_energy(phi)
    return phi, E_avg


# ============================================================
# ATOM FIELD BUILDERS
# ============================================================
def kink_field(R_arr, Z_eff=1.0):
    """Spherical kink (nuclear core). phi = (4/pi)*arctan(exp(-sqrt(Z)*r))"""
    R_g = xp.asarray(R_arr)
    return (4.0 / np.pi) * xp.arctan(xp.exp(-float(np.sqrt(Z_eff)) * R_g))


def breather_mode(X, Y, Z, R, l, m, amp, omega_n):
    """
    Breather with angular momentum (l,m) on the cubic lattice.
    Uses real spherical harmonics for the angular part.
    """
    X_g = xp.asarray(X)
    Y_g = xp.asarray(Y)
    Z_g = xp.asarray(Z)
    R_g = xp.asarray(R)
    eps = float(np.sqrt(max(1.0 - omega_n**2, 1e-12)))

    # Angular part (real spherical harmonics)
    if l == 0:
        ang = xp.ones_like(R_g)
    elif l == 1:
        if m == 0:   ang = Z_g / R_g     # p_z (sigma)
        elif m == 1: ang = X_g / R_g     # p_x (pi)
        elif m == -1: ang = Y_g / R_g    # p_y (pi)
        else: ang = xp.ones_like(R_g)
    else:
        ang = xp.ones_like(R_g)

    # Radial envelope
    radial = amp / (omega_n * xp.cosh(eps * R_g) + 1e-10)

    return (4.0 / np.pi) * xp.arctan(eps * ang * radial)


def build_atom_field(sym, center_z, Z_eff=1.0):
    """
    Build the complete 3D field for an atom at position center_z on the z-axis.

    Returns the field as a GPU array.
    """
    X, Y, Z, R = make_grid(center_z)

    # Nuclear core (kink)
    phi = kink_field(R, Z_eff)

    # Add valence breather modes based on atom type
    atom_modes = {
        'H':  [(0, 0, A_HALF, omega[0])],                          # 1s half
        'Li': [(0, 0, A_HALF, omega[0])],                          # 2s half
        'C':  [(1, 0, A_HALF, omega[1]), (1, 1, A_HALF, omega[1])],   # p_z + p_x half
        'N':  [(1, 0, A_HALF, omega[1]), (1, 1, A_HALF, omega[1]),
               (1, -1, A_HALF, omega[1])],                            # p_z + p_x + p_y all half
        'O':  [(1, 0, A_FULL, omega[1]), (1, 1, A_HALF, omega[1]),
               (1, -1, A_HALF, omega[1])],                            # p_z full, p_x,p_y half
        'F':  [(1, 0, A_HALF, omega[1]), (1, 1, A_FULL, omega[1]),
               (1, -1, A_FULL, omega[1])],                            # p_z half (bond), p_x,p_y full (LP)
        'Cl': [(1, 0, A_HALF, omega[1]), (1, 1, A_FULL, omega[1]),
               (1, -1, A_FULL, omega[1])],
    }

    modes = atom_modes.get(sym, [(0, 0, A_HALF, omega[0])])
    for l, m, amp, om in modes:
        phi = phi + breather_mode(X, Y, Z, R, l, m, amp, om)

    return phi


# ============================================================
# BOND SIMULATION
# ============================================================
def simulate_bond_3d(sym_a, sym_b, R_bond, Z_a=1.0, Z_b=1.0, n_steps=4000):
    """
    Place two atoms at ±R_bond/2 along z-axis.
    Evolve the combined field.
    Measure: E_combined - E_A_isolated - E_B_isolated = V_interaction
    """
    t0 = time.time()

    # Build isolated atom fields and measure their energies
    phi_a = build_atom_field(sym_a, -R_bond/2, Z_a)
    _, E_a = evolve(phi_a, n_steps=n_steps // 2, damping=0.02)

    phi_b = build_atom_field(sym_b, +R_bond/2, Z_b)
    _, E_b = evolve(phi_b, n_steps=n_steps // 2, damping=0.02)

    # Build combined field (superposition as initial condition, then evolve)
    # The NONLINEAR evolution will find the true equilibrium
    phi_combined = build_atom_field(sym_a, -R_bond/2, Z_a) + build_atom_field(sym_b, +R_bond/2, Z_b)
    # Subtract vacuum (double-counted background)
    # Actually for kink fields, the vacuum is phi=0, so superposition is fine as initial condition

    _, E_ab = evolve(phi_combined, n_steps=n_steps, damping=0.02)

    V_int = E_ab - E_a - E_b
    elapsed = time.time() - t0

    return V_int, E_ab, E_a, E_b, elapsed


# ============================================================
# RUN: H2 first
# ============================================================
print("=" * 70)
print("  3D BOND SIMULATION ON CUBIC LATTICE (GPU)")
print("  Two kinks + breather modes, full nonlinear evolution")
print("  All forces from L = (1/2)(dphi)^2 + (1/pi^2)(1-cos(pi*phi))")
print("=" * 70)
print()

# H2: scan R to find equilibrium
print("=== H2: Scanning bond distance ===")
R_values = [1.0, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0]

print(f"  {'R':>5} {'V_int':>10} {'E_ab':>10} {'E_a':>10} {'E_b':>10} {'time':>6}")
print("  " + "-" * 50)

results = []
for R in R_values:
    V, E_ab, E_a, E_b, t = simulate_bond_3d('H', 'H', R, Z_a=1.0, Z_b=1.0, n_steps=3000)
    results.append((R, V))
    print(f"  {R:5.1f} {V:+10.4f} {E_ab:10.4f} {E_a:10.4f} {E_b:10.4f} {t:5.1f}s")

# Find minimum
R_arr = [r[0] for r in results]
V_arr = [r[1] for r in results]
i_min = np.argmin(V_arr)
print(f"\n  Minimum: R = {R_arr[i_min]:.1f} Bohr, V = {V_arr[i_min]:.4f}")
print(f"  Observed: R_eq = 1.401 Bohr, D_e = 4.748 eV")
print()

# If H2 works, try F2 (the LP test)
print("=== F2: The LP repulsion test ===")
R_f2_values = [1.2, 1.4, 1.6, 2.0, 2.5, 3.0]
print(f"  {'R':>5} {'V_int':>10} {'time':>6}")
for R in R_f2_values:
    V, _, _, _, t = simulate_bond_3d('F', 'F', R, Z_a=9.0, Z_b=9.0, n_steps=3000)
    print(f"  {R:5.1f} {V:+10.4f} {t:5.1f}s")
