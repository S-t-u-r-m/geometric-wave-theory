#!/usr/bin/env python3
"""
Relaxed Bond Simulation — Nonlinear Field on d=3 Lattice
=========================================================
No linear superposition. For each separation R:
1. Initialize field as phi_A + phi_B (two breathers)
2. RELAX the field by minimizing the Lagrangian energy
3. The cosine potential naturally handles bonding vs repulsion:
   - phi stays in well (< 2): bonding (energy lowering)
   - phi pushed over barrier (> 2): lattice resists → repulsion

The field equation (static sine-Gordon):
  d²phi/dx² = (1/pi) sin(pi*phi)

Relaxation: gradient descent on E = integral[(1/2)(dphi)² + (1/pi²)(1-cos(pi*phi))]

Oh angular coupling applied AFTER the radial interaction is computed.
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial

# ============================================================
# CONSTANTS
# ============================================================
d = 3
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6
C_bond = np.pi / d
w_ch = [np.cos(k * np.pi / d) for k in range(d + 1)]

# Breather amplitudes
A_HALF = 2.0 / np.pi   # ~0.637 — single breather, stays in well
A_FULL = 4.0 / np.pi   # ~1.273 — double breather, near barrier top

# Grid
N = 2000
X_MAX = 12.0
DX = 2 * X_MAX / N
x = np.linspace(-X_MAX, X_MAX, N)


# ============================================================
# FIELD OPERATIONS
# ============================================================
def sech_profile(x, center, L, amp):
    """Breather spatial profile: amp * sech((x-center)/L)"""
    arg = np.clip((x - center) / L, -40, 40)
    return amp / np.cosh(arg)


def energy(phi):
    """Total Lagrangian energy of field configuration."""
    dphi = np.gradient(phi, DX)
    kinetic = 0.5 * dphi**2
    potential = (1.0 / np.pi**2) * (1.0 - np.cos(np.pi * phi))
    return np.sum(kinetic + potential) * DX


def force_field(phi):
    """
    Functional derivative dE/dphi = -d²phi/dx² + (1/pi)*sin(pi*phi)
    This is what drives the relaxation.
    """
    # Laplacian via finite differences
    d2phi = np.zeros_like(phi)
    d2phi[1:-1] = (phi[2:] - 2*phi[1:-1] + phi[:-2]) / DX**2
    d2phi[0] = d2phi[1]
    d2phi[-1] = d2phi[-2]

    return -d2phi + (1.0/np.pi) * np.sin(np.pi * phi)


def relax_field(phi_init, n_steps=2000, dt=0.0005, boundary_width=30):
    """
    Relax the field to minimize energy.
    Uses gradient descent with clipping for stability.
    """
    phi = phi_init.copy()
    phi_boundary = phi_init.copy()

    # Pin boundaries
    mask = np.ones(N)
    mask[:boundary_width] = 0
    mask[-boundary_width:] = 0

    E_prev = energy(phi)
    for step in range(n_steps):
        dEdphi = force_field(phi)

        # Clip gradient to prevent blowup
        max_grad = 10.0
        dEdphi = np.clip(dEdphi, -max_grad, max_grad)

        phi_new = phi - dt * dEdphi * mask

        # Restore boundary values
        phi_new[:boundary_width] = phi_boundary[:boundary_width]
        phi_new[-boundary_width:] = phi_boundary[-boundary_width:]

        # Check for NaN
        if np.any(np.isnan(phi_new)):
            break

        # Accept step
        phi = phi_new

    return phi


def interaction_energy(R, L_a, L_b, amp_a, amp_b, do_relax=True):
    """
    Compute interaction energy between two breather modes at separation R.

    1. Build initial field: superposition of two breathers
    2. Relax to minimum energy (nonlinear)
    3. Subtract isolated breather energies

    The cosine potential naturally handles:
    - Bonding: field merges smoothly, energy drops
    - Repulsion: field hits barrier, energy rises after relaxation
    """
    # Individual breather fields
    phi_a = sech_profile(x, -R/2, L_a, amp_a)
    phi_b = sech_profile(x, +R/2, L_b, amp_b)

    # Isolated energies (exact, no relaxation needed for single breathers)
    E_a = energy(phi_a)
    E_b = energy(phi_b)

    # Combined field: initial superposition
    phi_init = phi_a + phi_b

    # Relax the combined field
    if do_relax:
        phi_relaxed = relax_field(phi_init, n_steps=300, dt=0.003)
    else:
        phi_relaxed = phi_init

    E_ab = energy(phi_relaxed)

    # Interaction energy
    V_int = E_ab - E_a - E_b

    # Also report max field amplitude (diagnostic)
    phi_max_init = np.max(phi_init)
    phi_max_relax = np.max(phi_relaxed)

    return V_int, phi_max_init, phi_max_relax


# ============================================================
# DIAGNOSTICS: Compare bonding vs LP overlap
# ============================================================
print("=" * 75)
print("  Relaxed Bond Simulation — Nonlinear Sine-Gordon on d=3 Lattice")
print("  L = (1/2)(dphi)^2 + (1/pi^2)(1-cos(pi*phi))")
print("  Cosine well boundary at phi = 2 (next minimum)")
print("=" * 75)
print()

L_test = 0.49  # H-like breather width
R_test = 1.4   # H2 equilibrium

print(f"Mode overlap at R = {R_test} Bohr, L = {L_test}:")
print(f"  {'Config':25s} {'V_raw':>8} {'V_relax':>8} {'phi_init':>8} {'phi_relax':>8} {'Bonding?':>10}")

for label, amp_a, amp_b in [
    ("half + half (bond)",  A_HALF, A_HALF),
    ("full + half (LP+bd)", A_FULL, A_HALF),
    ("full + full (LP+LP)", A_FULL, A_FULL),
]:
    V_raw, _, _ = interaction_energy(R_test, L_test, L_test, amp_a, amp_b, do_relax=False)
    V_rel, phi_i, phi_r = interaction_energy(R_test, L_test, L_test, amp_a, amp_b, do_relax=True)
    bonding = "BOND" if V_rel < -0.001 else ("REPEL" if V_rel > 0.001 else "~zero")
    print(f"  {label:25s} {V_raw:+8.4f} {V_rel:+8.4f} {phi_i:8.3f} {phi_r:8.3f} {bonding:>10}")

print()

# ============================================================
# H2 CALIBRATION
# ============================================================
print("=== H2 Potential Curve (half + half, L=0.49) ===")
R_scan = np.linspace(0.4, 5.0, 100)
V_scan = np.zeros(len(R_scan))
phi_max_scan = np.zeros(len(R_scan))

for i, R in enumerate(R_scan):
    V_scan[i], _, phi_max_scan[i] = interaction_energy(R, L_test, L_test, A_HALF, A_HALF, do_relax=True)

i_min = np.argmin(V_scan)
V_min = V_scan[i_min]
R_eq = R_scan[i_min]
E_SCALE = -4.748 / V_min if V_min != 0 else 1.0

print(f"  R_eq = {R_eq:.3f} Bohr (obs: 1.401)")
print(f"  V_min(raw) = {V_min:.6f}")
print(f"  E_SCALE = {E_SCALE:.3f} eV/unit")
print(f"  D_e = {-V_min * E_SCALE:.3f} eV")
print()

# Show curve
print("  R      V(eV)    phi_max")
for i in range(0, len(R_scan), 5):
    print(f"  {R_scan[i]:5.2f}  {V_scan[i]*E_SCALE:+8.3f}  {phi_max_scan[i]:.3f}")

print()

# ============================================================
# FULL+FULL (LP) POTENTIAL CURVE
# ============================================================
print("=== LP Potential Curve (full + full, L=0.49) ===")
V_lp_scan = np.zeros(len(R_scan))
for i, R in enumerate(R_scan):
    V_lp_scan[i], _, _ = interaction_energy(R, L_test, L_test, A_FULL, A_FULL, do_relax=True)

i_min_lp = np.argmin(V_lp_scan)  # or argmax for repulsion
print(f"  V at R=1.4: {V_lp_scan[np.argmin(np.abs(R_scan-1.4))]*E_SCALE:+.3f} eV")
print(f"  V at R=2.0: {V_lp_scan[np.argmin(np.abs(R_scan-2.0))]*E_SCALE:+.3f} eV")
print(f"  V at R=3.0: {V_lp_scan[np.argmin(np.abs(R_scan-3.0))]*E_SCALE:+.3f} eV")
print()
print("  R      V_bond(eV)  V_LP(eV)   Ratio")
for i in range(0, len(R_scan), 5):
    vb = V_scan[i] * E_SCALE
    vlp = V_lp_scan[i] * E_SCALE
    ratio = vlp / vb if abs(vb) > 0.001 else 0
    print(f"  {R_scan[i]:5.2f}  {vb:+9.3f}  {vlp:+9.3f}  {ratio:+6.2f}")

print()

# ============================================================
# QUICK BOND TEST: F2 (the hard case)
# ============================================================
print("=== F2 Test: sigma(half+half) + 2×LP(full+full) ===")
L_F = 1.0  # n=2, L = n/2

# Sigma (bonding)
R_f2 = np.linspace(0.5, 6.0, 80)
V_bond_f2 = np.zeros(len(R_f2))
V_lp_f2 = np.zeros(len(R_f2))

for i, R in enumerate(R_f2):
    V_bond_f2[i], _, _ = interaction_energy(R, L_F, L_F, A_HALF, A_HALF, do_relax=True)
    V_lp_f2[i], _, _ = interaction_energy(R, L_F, L_F, A_FULL, A_FULL, do_relax=True)

# F2: 1 sigma bond + 2 LP facing
# Oh weights: sigma = 1, LP pi-channel = 0.5 each
E_harm_F2 = atoms_IE_F = 17.423
E_ratio_F2 = E_harm_F2 / E_H

V_total_f2 = (V_bond_f2 * w_ch[0] + 2 * V_lp_f2 * w_ch[1]) * E_SCALE * E_ratio_F2

i_min_f2 = np.argmin(V_total_f2)
print(f"  D_e(F2) = {-V_total_f2[i_min_f2]:.3f} eV (obs: 1.602)")
print(f"  R_eq(F2) = {R_f2[i_min_f2]:.3f} Bohr (obs: 1.412)")
print()
print("  R      V_bond    V_LP     V_total")
for i in range(0, len(R_f2), 5):
    vb = V_bond_f2[i] * E_SCALE * E_ratio_F2 * w_ch[0]
    vlp = 2 * V_lp_f2[i] * E_SCALE * E_ratio_F2 * w_ch[1]
    print(f"  {R_f2[i]:5.2f}  {vb:+8.3f}  {vlp:+8.3f}  {vb+vlp:+8.3f}")
