#!/usr/bin/env python3
"""
Breather-Breather Potential from the Lagrangian
================================================
Numerically compute the interaction between two sine-Gordon breathers
on the d=3 cubic lattice.

L = (1/2)(dphi/dx)^2 + (1/pi^2)(1 - cos(pi*phi))

The breather is a bound kink-antikink pair:
  phi_breather(x) = (4/pi) * arctan[sin(w*t) / cosh(x/L)]

At t=0 (maximum amplitude):
  phi(x) = 0 (the breather passes through zero)

At t = pi/(2w) (peak amplitude):
  phi(x) = (4/pi) * arctan[1/cosh(x/L)] = (2/pi) * am(x/L)

The peak spatial profile is approximately sech(x/L) scaled by 2/pi.

For two breathers separated by R, we compute the total field and energy.
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial

# Constants from d=3
d = 3
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6  # 13.604 eV
C_bond = np.pi / d
w_pi = np.cos(np.pi / d)  # 0.5

# Spatial grid for numerical integration
N_grid = 4000
x_max = 20.0  # Bohr
dx = 2 * x_max / N_grid
x = np.linspace(-x_max, x_max, N_grid)


def breather_profile(x, center, L):
    """
    Sine-Gordon breather at peak amplitude, centered at 'center'.

    phi(x) = (4/pi) * arctan(1 / cosh((x - center) / L))

    This is the exact static profile at the moment of maximum amplitude.
    L = breather width parameter ~ n * a_0 for principal quantum number n.
    """
    arg = (x - center) / L
    # Avoid overflow in cosh
    arg = np.clip(arg, -50, 50)
    return (4.0 / np.pi) * np.arctan(1.0 / np.cosh(arg))


def lagrangian_energy(phi, dphi_dx):
    """
    Energy density from the Lagrangian:
    E = (1/2)(dphi/dx)^2 + (1/pi^2)(1 - cos(pi*phi))
    """
    kinetic = 0.5 * dphi_dx**2
    potential = (1.0 / np.pi**2) * (1.0 - np.cos(np.pi * phi))
    return kinetic + potential


def compute_interaction(R, L_a, L_b, E_scale):
    """
    Compute the interaction energy between two breathers at separation R.

    Method:
    1. Place breather A at -R/2, breather B at +R/2
    2. Compute total field phi_total = phi_A + phi_B
    3. Compute energy of combined system: E_AB
    4. Compute energy of isolated A and B: E_A + E_B
    5. Interaction = E_AB - E_A - E_B

    This is the EXACT interaction from the Lagrangian.
    """
    # Individual breather profiles
    phi_a = breather_profile(x, -R/2, L_a)
    phi_b = breather_profile(x, +R/2, L_b)

    # Combined field (superposition)
    phi_total = phi_a + phi_b

    # Gradient of combined field
    dphi_total = np.gradient(phi_total, dx)

    # Energy of combined system
    E_density_AB = lagrangian_energy(phi_total, dphi_total)
    E_AB = np.trapezoid(E_density_AB, x)

    # Energy of isolated breathers (far apart)
    dphi_a = np.gradient(phi_a, dx)
    dphi_b = np.gradient(phi_b, dx)
    E_A = np.trapezoid(lagrangian_energy(phi_a, dphi_a), x)
    E_B = np.trapezoid(lagrangian_energy(phi_b, dphi_b), x)

    # Interaction energy (in Lagrangian units)
    V_int = E_AB - E_A - E_B

    # Convert to eV: scale by E_scale
    # The Lagrangian energy in natural units relates to physical energy by E_scale
    V_eV = V_int * E_scale

    return V_eV, E_AB * E_scale, E_A * E_scale, E_B * E_scale


def scan_potential(L_a, L_b, E_scale, R_min=0.2, R_max=6.0, N_R=300):
    """Scan the potential V(R) for a range of separations."""
    R_arr = np.linspace(R_min, R_max, N_R)
    V_arr = np.zeros(N_R)

    for i, R in enumerate(R_arr):
        V_arr[i], _, _, _ = compute_interaction(R, L_a, L_b, E_scale)

    return R_arr, V_arr


# ============================================================
# H2 VALIDATION
# ============================================================
print("=" * 70)
print("  Breather-Breather Potential from the Lagrangian")
print("  L = (1/2)(dphi/dx)^2 + (1/pi^2)(1-cos(pi*phi))")
print("=" * 70)
print()

# H2: two n=1 breathers
# Breather width L = 1 Bohr (for n=1, L ~ n*a_0 but in Bohr a_0 = 1)
# Energy scale: for sigma bond, project onto 1/d of coupling
# E_scale converts Lagrangian units to eV via the channel coupling
# From H2: D_e = pi*E_H/d^2 at sin(2R) = 1/d
# The Lagrangian energy unit = C_bond * E_H / normalization

print("=== H2: Two n=1 breathers ===")
print(f"  E_H = {E_H:.3f} eV")
print(f"  C_bond = pi/d = {C_bond:.5f}")
print()

# Scan with different L and E_scale to calibrate
L_H = 1.0  # Breather width for H (n=1), in Bohr

# The energy scale maps Lagrangian units → eV
# We need to find it by matching H2: D_e = 4.748 eV at R_eq = 1.401 Bohr

# First, scan the raw Lagrangian potential (E_scale = 1)
R_arr, V_raw = scan_potential(L_H, L_H, 1.0, R_min=0.3, R_max=5.0, N_R=200)

# Find minimum
i_min = np.argmin(V_raw)
V_min_raw = V_raw[i_min]
R_eq_raw = R_arr[i_min]

print(f"  Raw Lagrangian potential:")
print(f"    R_eq = {R_eq_raw:.3f} Bohr")
print(f"    V_min = {V_min_raw:.6f} (Lagrangian units)")
print()

# Calibrate E_scale so that V_min × E_scale = -D_e
D_e_H2 = 4.748  # eV (equilibrium dissociation energy)
E_scale_calibrated = -D_e_H2 / V_min_raw if V_min_raw != 0 else 1.0

print(f"  Calibrated E_scale = {E_scale_calibrated:.3f} eV per Lagrangian unit")
print(f"  -> V_min = {V_min_raw * E_scale_calibrated:.3f} eV (should be -4.748)")
print(f"  -> R_eq = {R_eq_raw:.3f} Bohr (should be 1.401)")
print()

# Show the potential curve
print("  R(Bohr)   V(eV)     V(raw)")
for i in range(0, len(R_arr), 10):
    if abs(i - i_min) < 25 or i % 30 == 0:
        print(f"  {R_arr[i]:6.3f}   {V_raw[i]*E_scale_calibrated:+8.3f}   {V_raw[i]:+10.6f}")

print()

# Now scan more finely around the minimum
R_fine, V_fine = scan_potential(L_H, L_H, E_scale_calibrated, R_min=0.5, R_max=4.0, N_R=500)
i_min_f = np.argmin(V_fine)
print(f"  Fine scan: R_eq = {R_fine[i_min_f]:.4f} Bohr, D_e = {-V_fine[i_min_f]:.4f} eV")

# ============================================================
# KEY QUESTION: Does the shape match?
# ============================================================
print()
print("=== KEY PHYSICS: Does the Lagrangian give the right potential SHAPE? ===")
print()

# If the shape is right, then:
# 1. R_eq should scale with n (period)
# 2. D_e should scale with E_ion
# 3. The Oh lookup table handles the angular part
# 4. Only E_scale needs to be set once — then ALL bonds predicted

# Test: vary L (breather size) and see how R_eq and V_min change
print("  Effect of breather width L on equilibrium:")
print(f"  {'L':>6} {'R_eq':>8} {'V_min(raw)':>12} {'D_e(eV)':>10}")
for L_test in [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]:
    R_t, V_t = scan_potential(L_test, L_test, 1.0, R_min=0.2, R_max=8.0, N_R=200)
    i_m = np.argmin(V_t)
    D_e_t = -V_t[i_m] * E_scale_calibrated
    print(f"  {L_test:6.2f} {R_t[i_m]:8.3f} {V_t[i_m]:12.6f} {D_e_t:10.3f}")
