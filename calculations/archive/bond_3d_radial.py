#!/usr/bin/env python3
"""
3D Radial Curve Precomputation — Oh Lookup Engine Backend
==========================================================
Compute 3 key radial interaction curves on the 3D cubic lattice:
  1. SIGMA: half+half p_z modes (bonding along z-axis)
  2. PI:    half+half p_x modes (bonding perpendicular to z)
  3. LP:    full+full p_x modes (lone pair overlap perpendicular)

These 3 curves + the Oh angular lookup table = complete bond predictions.

Method:
  - Place two breather modes at ±R/2 along z (no kinks — isolate the breather interaction)
  - Evolve on 3D cubic lattice with damping → find equilibrium
  - Measure interaction energy V(R) = E_combined - E_A - E_B
  - The 3D geometry naturally gives sigma ≠ pi ≠ LP

Higher resolution: 96^3 grid, smaller box focused on bond region.
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial
import time

try:
    import cupy as cp
    xp = cp
    GPU = True
    print("GPU: CuPy + CUDA active")
except ImportError:
    xp = np
    GPU = False
    print("GPU: CuPy not available, using NumPy")

# ============================================================
# CONSTANTS
# ============================================================
d = 3
V_0 = 1.0 / np.pi**2
gamma_sg = np.pi / (2**(d+1) * np.pi - 2)
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6

# Breather parameters
omega_p = float(np.cos(2 * gamma_sg))  # p-mode frequency
eps_p = float(np.sqrt(max(1.0 - omega_p**2, 1e-12)))

A_HALF = 0.08    # half-filled mode amplitude
A_FULL = 0.16    # full mode amplitude (2x)

# ============================================================
# 3D GRID — higher resolution, smaller box
# ============================================================
N = 96   # 96^3 ≈ 885K points
BOX = 6.0  # ±6 Bohr — focused on bond region
dx = 2 * BOX / N
dt = 0.25 * dx

print(f"Grid: {N}^3 = {N**3:,} points, box = ±{BOX}, dx = {dx:.4f}")
print(f"dt = {dt:.5f}, omega_p = {omega_p:.4f}, eps_p = {eps_p:.4f}")
print()

# Coordinate arrays
x1d = xp.linspace(-BOX, BOX, N, endpoint=False, dtype=np.float64)
# Meshgrid on GPU
X, Y, Z = xp.meshgrid(x1d, x1d, x1d, indexing='ij')


# ============================================================
# FIELD OPERATIONS
# ============================================================
def laplacian(phi):
    return (
        xp.roll(phi, 1, 0) + xp.roll(phi, -1, 0) +
        xp.roll(phi, 1, 1) + xp.roll(phi, -1, 1) +
        xp.roll(phi, 1, 2) + xp.roll(phi, -1, 2) -
        6 * phi
    ) / dx**2


def total_energy(phi):
    GE = 0.5 * (
        (xp.roll(phi, 1, 0) - phi)**2 +
        (xp.roll(phi, 1, 1) - phi)**2 +
        (xp.roll(phi, 1, 2) - phi)**2
    )
    PE = V_0 * (1.0 - xp.cos(np.pi * phi))
    return float(xp.sum(GE + PE) * dx**3)


def grad_energy(phi):
    """Gradient energy only."""
    GE = 0.5 * (
        (xp.roll(phi, 1, 0) - phi)**2 +
        (xp.roll(phi, 1, 1) - phi)**2 +
        (xp.roll(phi, 1, 2) - phi)**2
    )
    return float(xp.sum(GE) * dx**3)


def pot_energy(phi):
    """Cosine potential energy only."""
    PE = V_0 * (1.0 - xp.cos(np.pi * phi))
    return float(xp.sum(PE) * dx**3)


def evolve(phi_init, n_steps=4000, damping=0.02):
    """Damped evolution to equilibrium."""
    phi = phi_init.copy()
    phi_old = phi.copy()

    energies = []
    for step in range(n_steps):
        lap = laplacian(phi)
        force = (1.0 / np.pi) * xp.sin(np.pi * phi)
        phi_new = (2 - damping*dt) * phi - (1 - damping*dt) * phi_old + dt**2 * (lap - force)
        phi_old = phi.copy()
        phi = phi_new

        if step > n_steps * 3 // 4 and step % 50 == 0:
            energies.append(total_energy(phi))

    E_avg = np.mean(energies) if energies else total_energy(phi)
    return phi, E_avg


# ============================================================
# BREATHER MODE BUILDERS (no kink — pure valence modes)
# ============================================================
def p_z_mode(center_z, amp):
    """p_z breather: angular pattern = z/r, centered at center_z on z-axis."""
    Z_shifted = Z - center_z
    R = xp.sqrt(X**2 + Y**2 + Z_shifted**2) + 1e-10
    ang = Z_shifted / R  # p_z angular pattern
    radial = amp / (omega_p * xp.cosh(eps_p * R) + 1e-10)
    return (4.0 / np.pi) * xp.arctan(eps_p * ang * radial)


def p_x_mode(center_z, amp):
    """p_x breather: angular pattern = x/r, centered at center_z on z-axis."""
    Z_shifted = Z - center_z
    R = xp.sqrt(X**2 + Y**2 + Z_shifted**2) + 1e-10
    ang = X / R  # p_x angular pattern (perpendicular to bond)
    radial = amp / (omega_p * xp.cosh(eps_p * R) + 1e-10)
    return (4.0 / np.pi) * xp.arctan(eps_p * ang * radial)


def p_y_mode(center_z, amp):
    """p_y breather: angular pattern = y/r."""
    Z_shifted = Z - center_z
    R = xp.sqrt(X**2 + Y**2 + Z_shifted**2) + 1e-10
    ang = Y / R
    radial = amp / (omega_p * xp.cosh(eps_p * R) + 1e-10)
    return (4.0 / np.pi) * xp.arctan(eps_p * ang * radial)


# ============================================================
# COMPUTE INTERACTION ENERGY FOR A MODE PAIR
# ============================================================
def compute_V(R_bond, mode_func, amp_a, amp_b, n_evolve=3000):
    """
    Interaction energy of two breather modes at separation R_bond.

    mode_func: function(center_z, amp) -> 3D field
    amp_a, amp_b: amplitudes for atom A and B
    """
    # Build isolated fields
    phi_a = mode_func(-R_bond/2, amp_a)
    phi_b = mode_func(+R_bond/2, amp_b)

    # Isolated energies (quick evolve)
    _, E_a = evolve(phi_a, n_steps=n_evolve//2, damping=0.03)
    _, E_b = evolve(phi_b, n_steps=n_evolve//2, damping=0.03)

    # Combined: start from superposition, evolve nonlinearly
    phi_init = phi_a + phi_b
    _, E_ab = evolve(phi_init, n_steps=n_evolve, damping=0.02)

    return E_ab - E_a - E_b


# ============================================================
# PRECOMPUTE THE 3 KEY RADIAL CURVES
# ============================================================
R_scan = np.array([0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0])

print("=" * 70)
print("  PRECOMPUTING 3 RADIAL CURVES ON 3D CUBIC LATTICE")
print("  These curves + Oh lookup = complete bond predictions")
print("=" * 70)
print()

# --- CURVE 1: SIGMA (half+half p_z along bond axis) ---
print("CURVE 1: SIGMA — half+half p_z (bonding along z)")
print(f"  {'R':>5} {'V_sigma':>10} {'time':>6}")
V_sigma = []
for R in R_scan:
    t0 = time.time()
    V = compute_V(R, p_z_mode, A_HALF, A_HALF, n_evolve=3000)
    elapsed = time.time() - t0
    V_sigma.append(V)
    print(f"  {R:5.2f} {V:+10.6f} {elapsed:5.1f}s")

print()

# --- CURVE 2: PI (half+half p_x perpendicular to bond) ---
print("CURVE 2: PI — half+half p_x (bonding perpendicular)")
print(f"  {'R':>5} {'V_pi':>10} {'time':>6}")
V_pi = []
for R in R_scan:
    t0 = time.time()
    V = compute_V(R, p_x_mode, A_HALF, A_HALF, n_evolve=3000)
    elapsed = time.time() - t0
    V_pi.append(V)
    print(f"  {R:5.2f} {V:+10.6f} {elapsed:5.1f}s")

print()

# --- CURVE 3: LP (full+full p_x — lone pair overlap) ---
print("CURVE 3: LP — full+full p_x (lone pair repulsion?)")
print(f"  {'R':>5} {'V_LP':>10} {'time':>6}")
V_lp = []
for R in R_scan:
    t0 = time.time()
    V = compute_V(R, p_x_mode, A_FULL, A_FULL, n_evolve=3000)
    elapsed = time.time() - t0
    V_lp.append(V)
    print(f"  {R:5.2f} {V:+10.6f} {elapsed:5.1f}s")

print()

# ============================================================
# ANALYSIS
# ============================================================
V_sigma = np.array(V_sigma)
V_pi = np.array(V_pi)
V_lp = np.array(V_lp)

print("=" * 70)
print("  COMPARISON OF 3 CURVES")
print("=" * 70)
print(f"  {'R':>5} {'V_sigma':>10} {'V_pi':>10} {'V_LP':>10} {'LP/sig':>8} {'pi/sig':>8}")
for i, R in enumerate(R_scan):
    ratio_lp = V_lp[i] / V_sigma[i] if abs(V_sigma[i]) > 1e-8 else 0
    ratio_pi = V_pi[i] / V_sigma[i] if abs(V_sigma[i]) > 1e-8 else 0
    print(f"  {R:5.2f} {V_sigma[i]:+10.6f} {V_pi[i]:+10.6f} {V_lp[i]:+10.6f} {ratio_lp:+8.3f} {ratio_pi:+8.3f}")

print()
print("KEY QUESTIONS:")
print("  1. Is V_sigma attractive (negative)? → bonding works")
print("  2. Is V_LP repulsive (positive)? → LP repulsion from 3D geometry")
print("  3. Is V_pi weaker than V_sigma? → angular dependence correct")
print("  4. Ratio V_pi/V_sigma ≈ 0.5 (= w_pi = cos(pi/3))? → Oh prediction")

# Find sigma minimum
i_min = np.argmin(V_sigma)
print(f"\n  Sigma minimum: R = {R_scan[i_min]:.2f}, V = {V_sigma[i_min]:.6f}")
print(f"  (Need to calibrate to D_e(H2) = 4.748 eV)")
