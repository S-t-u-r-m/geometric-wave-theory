#!/usr/bin/env python3
"""
3D Bond Simulation v2 — With Kink Wells + Proper Amplitudes
=============================================================
The bonding mechanism: breather (electron) resonantly couples between
two kink wells (nuclei). The breather lowers its energy by spreading
across both wells — this IS the chemical bond.

Strategy:
  1. Create double-well potential: two kinks at ±R/2
  2. Place breather mode in the double well
  3. Compare energy to breather in single well
  4. Difference = bond energy

This isolates the BONDING interaction from the kink-kink repulsion.
The kink-kink repulsion is handled separately (it's just nuclear repulsion).

Run THREE configurations per R:
  A: kink_L + kink_R (nuclear repulsion baseline)
  B: kink_L + breather_L + kink_R (single-atom electron)
  C: kink_L + kink_R + breather_spread (bonded electron)

Bond energy from electron = E_C - E_B (energy gain from spreading)
Total bond = sum over all electron modes × Oh angular weights
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
    print("No GPU")

d = 3
V_0 = 1.0 / np.pi**2
gamma_sg = np.pi / (2**(d+1) * np.pi - 2)
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6

# Grid
N = 96
BOX = 8.0
dx = 2 * BOX / N
dt = 0.2 * dx

print(f"Grid: {N}^3, box=±{BOX}, dx={dx:.4f}")

x1d = xp.linspace(-BOX, BOX, N, endpoint=False, dtype=np.float64)
X, Y, Z = xp.meshgrid(x1d, x1d, x1d, indexing='ij')


def laplacian(phi):
    return (
        xp.roll(phi, 1, 0) + xp.roll(phi, -1, 0) +
        xp.roll(phi, 1, 1) + xp.roll(phi, -1, 1) +
        xp.roll(phi, 1, 2) + xp.roll(phi, -1, 2) -
        6 * phi
    ) / dx**2


def total_energy(phi):
    GE = 0.5 * ((xp.roll(phi,1,0)-phi)**2 + (xp.roll(phi,1,1)-phi)**2 + (xp.roll(phi,1,2)-phi)**2)
    PE = V_0 * (1.0 - xp.cos(np.pi * phi))
    return float(xp.sum(GE + PE) * dx**3)


def evolve(phi_init, n_steps=5000, damping=0.015):
    phi = phi_init.copy()
    phi_old = phi.copy()
    energies = []
    for step in range(n_steps):
        lap = laplacian(phi)
        force = (1.0/np.pi) * xp.sin(np.pi * phi)
        phi_new = (2-damping*dt)*phi - (1-damping*dt)*phi_old + dt**2*(lap - force)
        phi_old = phi.copy()
        phi = phi_new
        if step > n_steps*3//4 and step % 100 == 0:
            energies.append(total_energy(phi))
    return phi, np.mean(energies) if energies else total_energy(phi)


def kink(center_z, Z_eff=1.0):
    """Spherical kink at center_z along z-axis."""
    R = xp.sqrt(X**2 + Y**2 + (Z - center_z)**2) + 1e-10
    return (4.0/np.pi) * xp.arctan(xp.exp(-float(np.sqrt(Z_eff)) * R))


def breather_pz(center_z, amp=0.1):
    """p_z breather mode (sigma bonding, along z)."""
    Zs = Z - center_z
    R = xp.sqrt(X**2 + Y**2 + Zs**2) + 1e-10
    omega = float(np.cos(2*gamma_sg))
    eps = float(np.sqrt(1-omega**2))
    ang = Zs / R
    rad = amp / (omega * xp.cosh(eps * R) + 1e-10)
    return (4.0/np.pi) * xp.arctan(eps * ang * rad)


def breather_px(center_z, amp=0.1):
    """p_x breather mode (pi bonding or LP, perpendicular)."""
    Zs = Z - center_z
    R = xp.sqrt(X**2 + Y**2 + Zs**2) + 1e-10
    omega = float(np.cos(2*gamma_sg))
    eps = float(np.sqrt(1-omega**2))
    ang = X / R
    rad = amp / (omega * xp.cosh(eps * R) + 1e-10)
    return (4.0/np.pi) * xp.arctan(eps * ang * rad)


# ============================================================
# APPROACH: measure the full atom-atom interaction directly
# but with MUCH better statistics (longer evolution, more samples)
# ============================================================

def bond_energy(sym_a, sym_b, R_bond, n_steps=6000):
    """
    Full bond energy at separation R.

    Config A: atom_a isolated (kink + breathers)
    Config B: atom_b isolated
    Config C: both atoms together (kink + kink + all breathers)

    V(R) = E_C - E_A - E_B
    """
    # Atom configs: (Z_eff, [(mode_func, amplitude), ...])
    atom_config = {
        'H':  (1.0,  [(breather_pz, 0.5)]),
        'F':  (9.0,  [(breather_pz, 0.5), (breather_px, 1.0)]),  # pz half, px full (LP)
        'N':  (7.0,  [(breather_pz, 0.5), (breather_px, 0.5)]),
        'O':  (8.0,  [(breather_pz, 1.0), (breather_px, 0.5)]),  # pz full (LP), px half
        'C':  (6.0,  [(breather_pz, 0.5), (breather_px, 0.5)]),
    }

    Za, modes_a = atom_config.get(sym_a, (1.0, [(breather_pz, 0.5)]))
    Zb, modes_b = atom_config.get(sym_b, (1.0, [(breather_pz, 0.5)]))

    # Build isolated atom A
    phi_a = kink(-R_bond/2, Za)  # Use same position as in combined
    for mode_fn, amp in modes_a:
        phi_a = phi_a + mode_fn(-R_bond/2, amp)

    # Build isolated atom B (far away, at +R_bond/2 but we compute separately)
    phi_b = kink(+R_bond/2, Zb)
    for mode_fn, amp in modes_b:
        phi_b = phi_b + mode_fn(+R_bond/2, amp)

    # Build combined
    phi_c = kink(-R_bond/2, Za) + kink(+R_bond/2, Zb)
    for mode_fn, amp in modes_a:
        phi_c = phi_c + mode_fn(-R_bond/2, amp)
    for mode_fn, amp in modes_b:
        phi_c = phi_c + mode_fn(+R_bond/2, amp)

    t0 = time.time()

    # Evolve all three
    _, E_a = evolve(phi_a, n_steps=n_steps)
    _, E_b = evolve(phi_b, n_steps=n_steps)
    _, E_c = evolve(phi_c, n_steps=n_steps)

    V = E_c - E_a - E_b
    elapsed = time.time() - t0

    return V, E_c, E_a, E_b, elapsed


# ============================================================
# RUN
# ============================================================
print()
print("=" * 70)
print("  3D BOND SIMULATION v2")
print("  Kink wells + breathers, full nonlinear, GPU")
print("=" * 70)
print()

# H2 scan
print("=== H2 (Z=1+1, sigma bond) ===")
print(f"  {'R':>5} {'V_int':>10} {'E_AB':>10} {'E_A':>10} {'E_B':>10} {'time':>6}")
h2_results = []
for R in [1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 4.0, 5.0, 6.0]:
    V, Eab, Ea, Eb, t = bond_energy('H', 'H', R, n_steps=5000)
    h2_results.append((R, V))
    print(f"  {R:5.1f} {V:+10.5f} {Eab:10.5f} {Ea:10.5f} {Eb:10.5f} {t:5.1f}s")

# Find minimum
Rs = [r[0] for r in h2_results]
Vs = [r[1] for r in h2_results]
i_min = np.argmin(Vs)
print(f"\n  Min: R={Rs[i_min]:.1f}, V={Vs[i_min]:.5f}")
print(f"  V spans: {min(Vs):.5f} to {max(Vs):.5f}")
print(f"  Range (well depth): {max(Vs)-min(Vs):.5f}")

# If there's a minimum, estimate E_SCALE
if Vs[i_min] < 0:
    E_SC = -4.748 / Vs[i_min]
    print(f"  E_SCALE = {E_SC:.1f} eV/unit")
    print(f"  Predicted D_e = {-Vs[i_min]*E_SC:.3f} eV")

print()

# F2 if H2 looks reasonable
print("=== F2 (Z=9+9, sigma + LP) ===")
for R in [1.4, 2.0, 3.0, 4.0]:
    V, Eab, Ea, Eb, t = bond_energy('F', 'F', R, n_steps=5000)
    print(f"  R={R:4.1f} V={V:+10.5f} E_AB={Eab:10.4f} E_A={Ea:10.4f} E_B={Eb:10.4f} ({t:.0f}s)")
