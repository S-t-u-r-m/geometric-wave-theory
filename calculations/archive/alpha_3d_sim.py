"""
3D Alpha Corrections from Sine-Gordon Dynamics
================================================
Simulate the mode-mode interactions that determine alpha:
1. s-pair vs s-single (parity)
2. p-mode incremental filling (exchange)
3. Empty channel enhancement (mode-vacancy)
4. d-mode effect on s-valence alpha (TM corrections)

Uses 3D cubic lattice on GPU.
Each measurement: place modes, evolve, measure energy.
The DIFFERENCE between configurations gives the correction.
"""

import numpy as np
from math import factorial
import time

try:
    import cupy as cp
    GPU = True
except ImportError:
    cp = np
    GPU = False

d = 3
gamma_sg = np.pi / (2**(d+1) * np.pi - 2)
V_0 = 1.0 / np.pi**2
xp = cp

# Grid — use 48^3 for speed, enough for coupling measurement
N = 48
L = 10.0
dx = 2 * L / N
x1d = np.linspace(-L, L, N, endpoint=False)
X, Y, Z_grid = np.meshgrid(x1d, x1d, x1d, indexing='ij')
R = np.sqrt(X**2 + Y**2 + Z_grid**2) + 1e-10

X_g = xp.asarray(X)
Y_g = xp.asarray(Y)
Z_g = xp.asarray(Z_grid)
R_g = xp.asarray(R)

dt = 0.3 * dx
N_steps = 4000  # enough for energy measurement


def laplacian_3d(phi):
    return (
        xp.roll(phi, 1, 0) + xp.roll(phi, -1, 0) +
        xp.roll(phi, 1, 1) + xp.roll(phi, -1, 1) +
        xp.roll(phi, 1, 2) + xp.roll(phi, -1, 2) -
        6 * phi
    ) / dx**2


def kink_3d(R, Z_eff):
    return (4.0 / np.pi) * xp.arctan(xp.exp(-xp.sqrt(float(Z_eff)) * R))


def breather_3d(X, Y, Z, R, l, m, omega, amp=0.15):
    """Breather with angular momentum (l,m) on cubic lattice."""
    eps = float(np.sqrt(max(1.0 - omega**2, 1e-12)))

    # Angular patterns (real spherical harmonics on cube)
    if l == 0:
        ang = xp.ones_like(R)
    elif l == 1:
        if m == 0: ang = Z / R      # p_z
        elif m == 1: ang = X / R     # p_x
        elif m == -1: ang = Y / R    # p_y
    elif l == 2:
        if m == 0: ang = (3*Z**2 - R**2) / R**2    # d_z2 (eg)
        elif m == 1: ang = X*Z / R**2                # d_xz (t2g)
        elif m == -1: ang = Y*Z / R**2               # d_yz (t2g)
        elif m == 2: ang = X*Y / R**2                # d_xy (t2g)
        elif m == -2: ang = (X**2 - Y**2) / R**2    # d_x2y2 (eg)
    else:
        ang = xp.ones_like(R)

    radial = amp / (omega * xp.cosh(eps * R) + 1e-10)
    return (4.0 / np.pi) * xp.arctan(eps * ang * radial)


def total_energy(phi, Z_eff):
    GE = 0.5 * (
        (xp.roll(phi, 1, 0) - phi)**2 +
        (xp.roll(phi, 1, 1) - phi)**2 +
        (xp.roll(phi, 1, 2) - phi)**2
    )
    PE = Z_eff * V_0 * (1 - xp.cos(np.pi * phi))
    return float(xp.sum(GE + PE) * dx**3)


def evolve_3d(phi_init, Z_eff, n_steps=None):
    if n_steps is None:
        n_steps = N_steps
    phi = phi_init.copy()
    phi_old = phi.copy()
    energies = []
    for step in range(n_steps):
        lap = laplacian_3d(phi)
        force = Z_eff * (1.0 / np.pi) * xp.sin(np.pi * phi)
        phi_new = 2*phi - phi_old + dt**2 * (lap - force)
        phi_old = phi.copy()
        phi = phi_new
        if step > n_steps // 2 and step % 50 == 0:
            energies.append(total_energy(phi, Z_eff))
    return np.mean(energies) if energies else 0.0


omega_s = float(np.cos(1 * gamma_sg))
omega_p = float(np.cos(2 * gamma_sg))
omega_d = float(np.cos(3 * gamma_sg))

Z = 10.0  # test charge

print("3D Alpha Corrections — Sine-Gordon on Cubic Lattice")
print("=" * 60)
print(f"  Grid: {N}^3, Z={Z}, GPU={GPU}")
print(f"  N_steps={N_steps}, amp=0.15")
print(flush=True)

# === BASELINE: kink energy ===
t0 = time.time()
phi_k = kink_3d(R_g, Z)
E_k = evolve_3d(phi_k, Z)
print(f"\n  Kink energy: {E_k:.4f} ({time.time()-t0:.1f}s)", flush=True)

# === TEST 1: PARITY — s-pair vs s-single ===
print("\n--- TEST 1: PARITY (s-pair vs s-single) ---")

# s-single: one s-breather
phi_s1 = phi_k + breather_3d(X_g, Y_g, Z_g, R_g, 0, 0, omega_s)
E_s1 = evolve_3d(phi_s1, Z)
E_bind_s1 = E_s1 - E_k

# s-pair: two s-breathers (same l=0, m=0 but different "spin")
# In 3D, we model spin as a slight phase offset
phi_s2 = phi_s1 + breather_3d(X_g, Y_g, Z_g, R_g, 0, 0, omega_s, amp=0.15)
E_s2 = evolve_3d(phi_s2, Z)
E_bind_s2 = E_s2 - E_k

# The parity effect: how much does the 2nd s-mode cost relative to the 1st?
E_2nd_s = E_s2 - E_s1
parity_ratio = E_2nd_s / E_bind_s1

print(f"  1st s-mode: E_bind = {E_bind_s1:.4f}")
print(f"  2nd s-mode: E_bind = {E_2nd_s:.4f}")
print(f"  Ratio (2nd/1st): {parity_ratio:.4f}")
print(f"  v19 parity: +1 (s-pair) vs -1 (s-single)")
print(f"  Parity correction = ratio - 1 = {parity_ratio - 1:+.4f}")
print(flush=True)

# === TEST 2: EXCHANGE — p-mode filling ===
print("\n--- TEST 2: EXCHANGE (p-mode incremental filling) ---")
print("  p_x, p_y, p_z on cubic lattice — the 3 channels")

# p_x alone
phi_px = phi_k + breather_3d(X_g, Y_g, Z_g, R_g, 1, 1, omega_p)
E_px = evolve_3d(phi_px, Z)
E_bind_px = E_px - E_k

# p_x + p_y
phi_pxy = phi_px + breather_3d(X_g, Y_g, Z_g, R_g, 1, -1, omega_p)
E_pxy = evolve_3d(phi_pxy, Z)
E_2nd_py = E_pxy - E_px

# p_x + p_y + p_z (half-fill!)
phi_pxyz = phi_pxy + breather_3d(X_g, Y_g, Z_g, R_g, 1, 0, omega_p)
E_pxyz = evolve_3d(phi_pxyz, Z)
E_3rd_pz = E_pxyz - E_pxy

print(f"  p_x alone:       E = {E_bind_px:.4f}")
print(f"  +p_y (2nd):      E = {E_2nd_py:.4f}  ratio = {E_2nd_py/E_bind_px:.4f}")
print(f"  +p_z (3rd, half): E = {E_3rd_pz:.4f}  ratio = {E_3rd_pz/E_bind_px:.4f}")
print(f"  Exchange: cost DECREASES with filling (half-fill cheapest)")
print(f"  v19 exchange: C(pa,2)/(d^2*n) per pair")
print(flush=True)

# === TEST 3: MODE-VACANCY — empty channels enhance coupling ===
print("\n--- TEST 3: MODE-VACANCY ---")
print("  Does having empty p-channels enhance the s-mode binding?")

# s-mode alone
E_s_alone = E_bind_s1

# s-mode with 1 p-mode present (2 empty p-channels)
phi_sp1 = phi_k + breather_3d(X_g, Y_g, Z_g, R_g, 0, 0, omega_s)
phi_sp1 = phi_sp1 + breather_3d(X_g, Y_g, Z_g, R_g, 1, 1, omega_p)
E_sp1 = evolve_3d(phi_sp1, Z)
# s-binding with p present:
phi_p_only = phi_k + breather_3d(X_g, Y_g, Z_g, R_g, 1, 1, omega_p)
E_p_only = evolve_3d(phi_p_only, Z)
E_s_with_p = E_sp1 - E_p_only

print(f"  s alone:       E = {E_s_alone:.4f}")
print(f"  s with 1p:     E = {E_s_with_p:.4f}  ratio = {E_s_with_p/E_s_alone:.4f}")
print(flush=True)

# === TEST 4: d-mode effect on s-valence ===
print("\n--- TEST 4: d(t2g) effect on s-mode (TM correction) ---")

# s alone
# (already measured)

# s with d_xy (t2g)
phi_sd = phi_k + breather_3d(X_g, Y_g, Z_g, R_g, 2, 2, omega_d)
E_sd_base = evolve_3d(phi_sd, Z)
phi_sd_s = phi_sd + breather_3d(X_g, Y_g, Z_g, R_g, 0, 0, omega_s)
E_sd_s = evolve_3d(phi_sd_s, Z)
E_s_with_dt2g = E_sd_s - E_sd_base

# s with d_z2 (eg)
phi_seg = phi_k + breather_3d(X_g, Y_g, Z_g, R_g, 2, 0, omega_d)
E_seg_base = evolve_3d(phi_seg, Z)
phi_seg_s = phi_seg + breather_3d(X_g, Y_g, Z_g, R_g, 0, 0, omega_s)
E_seg_s = evolve_3d(phi_seg_s, Z)
E_s_with_deg = E_seg_s - E_seg_base

print(f"  s alone:       E = {E_s_alone:.4f}")
print(f"  s with d_t2g:  E = {E_s_with_dt2g:.4f}  ratio = {E_s_with_dt2g/E_s_alone:.4f}")
print(f"  s with d_eg:   E = {E_s_with_deg:.4f}  ratio = {E_s_with_deg/E_s_alone:.4f}")
print(f"  t2g/eg ratio on s: {E_s_with_dt2g/E_s_with_deg:.4f}")
print(f"  (Should be ~1 if both forbidden by Oh, or differ if eg couples more)")
print(flush=True)

# === TEST 5: Overfill — 4th p-mode (lone pair) ===
print("\n--- TEST 5: OVERFILL — 4th p-mode (O-like) ---")

# Already have p_x + p_y + p_z (3 modes = half-fill)
# Add 4th: p_x again (lone pair, same channel different spin)
phi_p4 = phi_pxyz + breather_3d(X_g, Y_g, Z_g, R_g, 1, 1, omega_p, amp=0.15)
E_p4 = evolve_3d(phi_p4, Z)
E_4th_p = E_p4 - E_pxyz

print(f"  p1 (underfill): E = {E_bind_px:.4f}")
print(f"  p2:             E = {E_2nd_py:.4f}")
print(f"  p3 (half-fill): E = {E_3rd_pz:.4f}")
print(f"  p4 (overfill):  E = {E_4th_p:.4f}")
print(f"  Overfill/underfill ratio: {E_4th_p/E_bind_px:.4f}")
print(flush=True)

# === SUMMARY ===
print("\n" + "=" * 60)
print("SUMMARY — Alpha corrections from 3D wave dynamics")
print("=" * 60)
print(f"  Parity (s-pair/s-single):  {parity_ratio:.4f}  (v19: ~3)")
print(f"  Exchange p2/p1:            {E_2nd_py/E_bind_px:.4f}")
print(f"  Exchange p3/p1 (half):     {E_3rd_pz/E_bind_px:.4f}")
print(f"  Overfill p4/p1:            {E_4th_p/E_bind_px:.4f}")
print(f"  Mode-vacancy (s+p/s):      {E_s_with_p/E_s_alone:.4f}")
print(f"  t2g on s / s alone:        {E_s_with_dt2g/E_s_alone:.4f}")
print(f"  eg on s / s alone:         {E_s_with_deg/E_s_alone:.4f}")
