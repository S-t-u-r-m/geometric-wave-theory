#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
2D Sine-Gordon with Angular Momentum — GPU Simulation
========================================================

Test whether delta-channel (l=2) modes anti-screen sigma-channel (l=0) modes,
as predicted by GWT: w_delta = cos(2*pi/d) = -0.5.

In 2D, angular momentum is real:
  - l=0 (s-mode/sigma): radially symmetric, can penetrate to r=0
  - l=1 (p-mode/pi): cos(theta) dependence, node at r=0
  - l=2 (d-mode/delta): cos(2*theta), pushed further out by centrifugal barrier

The "atom" = radial kink (nucleus) + angular breather modes (electrons).

EXPERIMENTS:
A) Baseline: binding energy of each l-mode to the radial kink
B) Cross-channel screening: how l=0 screens l=1, l=2, and vice versa
C) THE KEY TEST: does l=2 anti-screen l=0? (delta-sigma interaction)
D) Screening weight extraction: compare to cos(l*pi/d)

Uses CuPy for GPU acceleration on RTX 4070 Ti.
"""

import numpy as np
import time as timer

try:
    import cupy as cp
    GPU = True
    print("  GPU: CuPy + CUDA detected")
except ImportError:
    import numpy as cp
    GPU = False
    print("  GPU: Not available, using NumPy (slow)")

pi = np.pi
d_dim = 3  # spatial dimension for GWT constants

# Physics
V_0 = 1.0 / pi**2

# 2D Grid
N = 512  # 512x512 grid
L = 80.0  # physical size
dx = L / N
dt = 0.3 * dx  # CFL condition for 2D: dt < dx/sqrt(2)
x = cp.linspace(-L/2, L/2, N)
y = cp.linspace(-L/2, L/2, N)
X, Y = cp.meshgrid(x, y, indexing='ij')
R = cp.sqrt(X**2 + Y**2) + 1e-10  # radial distance
THETA = cp.arctan2(Y, X)           # polar angle

# Damping sponge at boundaries (absorb outgoing radiation)
sponge_width = L * 0.15
dist_from_edge = cp.minimum(
    cp.minimum(X - (-L/2), L/2 - X),
    cp.minimum(Y - (-L/2), L/2 - Y)
)
damping = cp.where(dist_from_edge < sponge_width,
                   0.1 * (1 - dist_from_edge / sponge_width)**2,
                   0.0)


def radial_kink(R_arr, r0=0.0, width=1.0):
    """Radial kink centered at origin. phi goes from ~1 at r=0 to 0 at r->inf.
    This is the 'nucleus' — a topological defect."""
    arg = (R_arr - r0) / width
    return (4.0 / pi) * cp.arctan(cp.exp(-cp.clip(arg, -30, 30)))


def angular_mode(R_arr, theta_arr, l, amplitude=0.3, r_peak=5.0, r_width=3.0):
    """Create an angular mode with quantum number l.
    Radial profile peaked at r_peak, angular dependence cos(l*theta).
    l=0: s-mode (sigma), l=1: p-mode (pi), l=2: d-mode (delta)."""
    # Radial envelope: Gaussian peaked at r_peak
    # For l>0, include r^l factor (centrifugal barrier pushes mode outward)
    if l == 0:
        radial = cp.exp(-(R_arr - r_peak)**2 / (2 * r_width**2))
    else:
        # r^l * exp(-(r-r_peak)^2 / ...) with node at r=0
        r_safe = cp.maximum(R_arr, 0.01)
        radial = (r_safe / r_peak)**l * cp.exp(-(R_arr - r_peak)**2 / (2 * r_width**2))

    # Angular part
    if l == 0:
        angular = cp.ones_like(theta_arr)
    else:
        angular = cp.cos(l * theta_arr)

    return amplitude * radial * angular


def total_energy(phi, dphi):
    """Total energy of the field configuration."""
    E_kin = 0.5 * cp.sum(dphi**2) * dx**2

    # Gradient energy (2D finite differences)
    grad_x = cp.zeros_like(phi)
    grad_y = cp.zeros_like(phi)
    grad_x[1:-1, :] = (phi[2:, :] - phi[:-2, :]) / (2*dx)
    grad_y[:, 1:-1] = (phi[:, 2:] - phi[:, :-2]) / (2*dx)
    E_grad = 0.5 * cp.sum(grad_x**2 + grad_y**2) * dx**2

    E_pot = V_0 * cp.sum(1.0 - cp.cos(pi * phi)) * dx**2
    return float(E_kin + E_grad + E_pot)


def energy_in_ring(phi, dphi, r_min, r_max):
    """Energy contained in an annular ring [r_min, r_max]."""
    mask = ((R >= r_min) & (R < r_max)).astype(phi.dtype)

    E_kin = 0.5 * cp.sum(dphi**2 * mask) * dx**2

    grad_x = cp.zeros_like(phi)
    grad_y = cp.zeros_like(phi)
    grad_x[1:-1, :] = (phi[2:, :] - phi[:-2, :]) / (2*dx)
    grad_y[:, 1:-1] = (phi[:, 2:] - phi[:, :-2]) / (2*dx)
    E_grad = 0.5 * cp.sum((grad_x**2 + grad_y**2) * mask) * dx**2

    E_pot = V_0 * cp.sum((1.0 - cp.cos(pi * phi)) * mask) * dx**2
    return float(E_kin + E_grad + E_pot)


def angular_projection(phi, l):
    """Project field onto angular mode l: integral of phi * cos(l*theta) over theta."""
    if l == 0:
        weight = cp.ones_like(THETA)
    else:
        weight = cp.cos(l * THETA) * 2  # factor 2 for normalization

    # Radial profile of the l-component
    proj = cp.sum(phi * weight, axis=1) * dx  # sum over y (theta integration)
    return cp.asnumpy(proj) if GPU else proj


def laplacian_2d(phi):
    """2D Laplacian with second-order finite differences."""
    lap = cp.zeros_like(phi)
    lap[1:-1, 1:-1] = (
        phi[2:, 1:-1] + phi[:-2, 1:-1] +
        phi[1:-1, 2:] + phi[1:-1, :-2] -
        4 * phi[1:-1, 1:-1]
    ) / dx**2
    return lap


def evolve(phi0, dphi0, n_steps, measure_every=200):
    """Leapfrog evolution with sponge layer damping."""
    phi = phi0.copy()
    dphi = dphi0.copy()
    energies = []

    for step in range(n_steps):
        # Leapfrog
        acc = laplacian_2d(phi) - V_0 * pi * cp.sin(pi * phi) - damping * dphi
        dphi_half = dphi + 0.5 * dt * acc
        phi = phi + dt * dphi_half
        acc_new = laplacian_2d(phi) - V_0 * pi * cp.sin(pi * phi) - damping * dphi_half
        dphi = dphi_half + 0.5 * dt * acc_new

        if step % measure_every == 0 and step > n_steps // 4:
            energies.append(total_energy(phi, dphi))

    mean_E = np.mean(energies) if energies else total_energy(phi, dphi)
    return mean_E, phi, dphi


def evolve_quick(phi0, dphi0, n_steps=6000):
    """Quick evolution returning just mean energy."""
    E, _, _ = evolve(phi0, dphi0, n_steps)
    return E


# =============================================================================
# PARAMETERS
# =============================================================================
# Kink parameters
kink_width = 2.0  # width of radial kink

# Mode parameters — where each l-mode peaks
# Higher l -> further out (centrifugal barrier)
r_peak = {0: 4.0, 1: 6.0, 2: 8.0}
r_width = {0: 2.5, 1: 2.5, 2: 3.0}
amp = 0.25  # breather amplitude (keep small for perturbative regime)

n_steps = 8000  # evolution steps

print("=" * 78)
print("  2D SINE-GORDON — Angular Momentum Mode Mixing (GPU)")
print("=" * 78)
print(f"  Grid: {N}x{N}, L={L}, dx={dx:.4f}, dt={dt:.4f}")
print(f"  Kink width: {kink_width}")
print(f"  Mode peaks: l=0 at r={r_peak[0]}, l=1 at r={r_peak[1]}, l=2 at r={r_peak[2]}")
print(f"  Amplitude: {amp}")
print(f"  Steps: {n_steps}")

t_start = timer.time()

# =============================================================================
# BASELINES
# =============================================================================
print(f"\n{'='*78}")
print("  BASELINES")
print(f"{'='*78}")

# Empty field
E_vacuum = evolve_quick(cp.zeros((N, N)), cp.zeros((N, N)), n_steps=2000)
print(f"  E_vacuum = {E_vacuum:.6f}")

# Kink alone
phi_kink = radial_kink(R, width=kink_width)
E_kink = evolve_quick(phi_kink, cp.zeros((N, N)))
print(f"  E_kink = {E_kink:.6f}")

# Free modes (no kink) — each l
free_E = {}
for l in [0, 1, 2]:
    phi_mode = angular_mode(R, THETA, l, amplitude=amp,
                            r_peak=r_peak[l], r_width=r_width[l])
    E_free = evolve_quick(phi_mode, cp.zeros((N, N)))
    free_E[l] = E_free - E_vacuum
    print(f"  E_free(l={l}) = {free_E[l]:.6f}  (peak at r={r_peak[l]})")

t_base = timer.time() - t_start
print(f"  Baselines done in {t_base:.1f}s")

# =============================================================================
# A) BINDING ENERGY: each l-mode to the radial kink
# =============================================================================
print(f"\n{'='*78}")
print("  A) BINDING ENERGY OF ANGULAR MODES TO KINK")
print("  E_bind = E(kink+mode) - E(kink) - E(mode_free)")
print(f"{'='*78}")

bind_E = {}
kink_mode_E = {}
for l in [0, 1, 2]:
    phi0 = radial_kink(R, width=kink_width) + \
           angular_mode(R, THETA, l, amplitude=amp,
                        r_peak=r_peak[l], r_width=r_width[l])
    E_km = evolve_quick(phi0, cp.zeros((N, N)))
    kink_mode_E[l] = E_km
    E_bind = E_km - E_kink - free_E[l]
    bind_E[l] = E_bind
    w_channel = np.cos(l * pi / d_dim)
    print(f"  l={l}: E_bind = {E_bind:+.6f}  (channel weight cos({l}*pi/d) = {w_channel:.4f})")

print(f"\n  Binding ratios:")
if abs(bind_E[0]) > 1e-10:
    for l in [1, 2]:
        ratio = bind_E[l] / bind_E[0]
        print(f"    E_bind(l={l}) / E_bind(l=0) = {ratio:.4f}")

t_A = timer.time() - t_start
print(f"  Part A done in {t_A:.1f}s")


# =============================================================================
# B) CROSS-CHANNEL SCREENING
# =============================================================================
print(f"\n{'='*78}")
print("  B) CROSS-CHANNEL SCREENING")
print("  How does one angular mode change another's binding to the kink?")
print("  S = E_bind(out | in present) / E_bind(out | alone)")
print("  GWT prediction: w(l_in, l_out) = cos(l_in * pi/d) * cos(l_out * pi/d)")
print(f"{'='*78}")

print(f"\n  {'l_in':>4} {'l_out':>5} {'dE_bare':>10} {'dE_screen':>10} "
      f"{'S_frac':>8} {'1-S':>8} {'GWT_pred':>10}")
print(f"  " + "-" * 65)

for l_in in [0, 1, 2]:
    for l_out in [0, 1, 2]:
        if l_in == l_out:
            continue

        # Kink + inner mode (already computed)
        # Kink + inner + outer
        phi0 = radial_kink(R, width=kink_width) + \
               angular_mode(R, THETA, l_in, amplitude=amp,
                            r_peak=r_peak[l_in], r_width=r_width[l_in]) + \
               angular_mode(R, THETA, l_out, amplitude=amp,
                            r_peak=r_peak[l_out], r_width=r_width[l_out])
        E_both = evolve_quick(phi0, cp.zeros((N, N)))

        # Binding of outer with inner present
        dE_screen = E_both - kink_mode_E[l_in] - free_E[l_out]
        dE_bare = bind_E[l_out]
        S_frac = dE_screen / dE_bare if abs(dE_bare) > 1e-10 else float('nan')
        screen = 1 - S_frac

        # GWT prediction: product of channel weights
        w_in = np.cos(l_in * pi / d_dim)
        w_out = np.cos(l_out * pi / d_dim)
        gwt_pred = w_in * w_out

        marker = ""
        if screen < -0.1:
            marker = " <-- ANTI-SCREEN!"
        elif abs(screen - 0.5) < 0.15:
            marker = " ~ w_pi"

        print(f"  l={l_in:1d}   l={l_out:1d}   {dE_bare:+10.6f} {dE_screen:+10.6f} "
              f"{S_frac:8.4f} {screen:+8.4f} {gwt_pred:+10.4f}{marker}")

t_B = timer.time() - t_start
print(f"  Part B done in {t_B:.1f}s")


# =============================================================================
# C) THE KEY TEST: l=2 screening of l=0 with varying amplitude
# =============================================================================
print(f"\n{'='*78}")
print("  C) KEY TEST: l=2 (delta) screening of l=0 (sigma)")
print("  Sweep l=2 amplitude to see if screening changes sign")
print("  GWT predicts: w_delta * w_sigma = -0.5 * 1.0 = -0.5 (anti-screening!)")
print(f"{'='*78}")

print(f"\n  {'amp_d':>6} {'dE_bare(s)':>12} {'dE_screen(s)':>13} "
      f"{'S_frac':>8} {'1-S':>8}")
print(f"  " + "-" * 55)

for amp_d in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]:
    # Kink + d-mode at this amplitude
    phi_kd = radial_kink(R, width=kink_width) + \
             angular_mode(R, THETA, 2, amplitude=amp_d,
                          r_peak=r_peak[2], r_width=r_width[2])
    E_kd = evolve_quick(phi_kd, cp.zeros((N, N)))

    # Kink + d-mode + s-mode
    phi_kds = phi_kd + angular_mode(R, THETA, 0, amplitude=amp,
                                     r_peak=r_peak[0], r_width=r_width[0])
    E_kds = evolve_quick(phi_kds, cp.zeros((N, N)))

    dE_screen = E_kds - E_kd - free_E[0]
    dE_bare = bind_E[0]
    S_frac = dE_screen / dE_bare if abs(dE_bare) > 1e-10 else float('nan')
    screen = 1 - S_frac

    marker = ""
    if screen < -0.05:
        marker = " <-- ANTI"
    elif screen < 0:
        marker = " (slightly anti)"

    print(f"  {amp_d:6.2f} {dE_bare:+12.6f} {dE_screen:+12.6f} "
          f"{S_frac:8.4f} {screen:+8.4f}{marker}")

t_C = timer.time() - t_start
print(f"  Part C done in {t_C:.1f}s")


# =============================================================================
# D) PAIRED l=2 MODE — does pairing make it transparent?
# =============================================================================
print(f"\n{'='*78}")
print("  D) PAIRED l=2 MODE SCREENING OF l=0")
print("  A d-mode PAIR (cos(2*theta) + sin(2*theta)) should be more")
print("  symmetric. Does pairing change the screening?")
print(f"{'='*78}")

# Paired d-mode: cos(2*theta) + sin(2*theta) = orthogonal pair
phi_d_pair = radial_kink(R, width=kink_width) + \
    angular_mode(R, THETA, 2, amplitude=amp, r_peak=r_peak[2], r_width=r_width[2])
# Add the sin(2*theta) component using a phase shift
phi_d_pair = phi_d_pair + amp * (R / r_peak[2])**2 * \
    cp.exp(-(R - r_peak[2])**2 / (2 * r_width[2]**2)) * cp.sin(2 * THETA)

E_kd_pair = evolve_quick(phi_d_pair, cp.zeros((N, N)))

# Add s-mode
phi_kds_pair = phi_d_pair + angular_mode(R, THETA, 0, amplitude=amp,
                                          r_peak=r_peak[0], r_width=r_width[0])
E_kds_pair = evolve_quick(phi_kds_pair, cp.zeros((N, N)))

dE_pair = E_kds_pair - E_kd_pair - free_E[0]
S_pair = dE_pair / bind_E[0] if abs(bind_E[0]) > 1e-10 else float('nan')

print(f"\n  Single l=2 screening of l=0: 1-S = {1 - (bind_E.get('B_20_S', S_frac)):.4f}")
print(f"  Paired l=2 screening of l=0: 1-S = {1 - S_pair:+.4f}")
print(f"  (If paired is more transparent -> confirms pairing rule)")


# =============================================================================
# E) SCREENING WEIGHT EXTRACTION
# =============================================================================
print(f"\n{'='*78}")
print("  E) SCREENING WEIGHT SUMMARY")
print("  Compare measured screening to GWT channel weights")
print(f"{'='*78}")

print(f"\n  GWT predictions (d={d_dim}):")
for l in [0, 1, 2]:
    w = np.cos(l * pi / d_dim)
    print(f"    w(l={l}) = cos({l}*pi/{d_dim}) = {w:.4f}")

print(f"\n  Product rule predictions:")
for l_in in [0, 1, 2]:
    for l_out in [0, 1, 2]:
        if l_in == l_out:
            continue
        w_in = np.cos(l_in * pi / d_dim)
        w_out = np.cos(l_out * pi / d_dim)
        print(f"    w({l_in},{l_out}) = {w_in:.4f} * {w_out:.4f} = {w_in * w_out:+.4f}")


# =============================================================================
# F) RADIAL BINDING PROFILE — where does each mode sit?
# =============================================================================
print(f"\n{'='*78}")
print("  F) RADIAL BINDING PROFILES")
print("  Energy density in radial shells for kink + each mode")
print(f"{'='*78}")

for l in [0, 1, 2]:
    phi0 = radial_kink(R, width=kink_width) + \
           angular_mode(R, THETA, l, amplitude=amp,
                        r_peak=r_peak[l], r_width=r_width[l])
    # Quick evolve to get settled state
    _, phi_settled, dphi_settled = evolve(phi0, cp.zeros((N, N)), 4000)

    print(f"\n  l={l} mode (peak at r={r_peak[l]}):")
    print(f"  {'r':>6} {'E_ring':>10} {'E_kink_ring':>12} {'diff':>10}")

    _, phi_k_settled, dphi_k_settled = evolve(
        radial_kink(R, width=kink_width), cp.zeros((N, N)), 4000)

    for r_lo in np.arange(0, 20, 2):
        r_hi = r_lo + 2
        E_ring = energy_in_ring(phi_settled, dphi_settled, r_lo, r_hi)
        E_k_ring = energy_in_ring(phi_k_settled, dphi_k_settled, r_lo, r_hi)
        diff = E_ring - E_k_ring
        print(f"  {r_lo:5.1f}-{r_hi:4.1f} {E_ring:10.6f} {E_k_ring:12.6f} {diff:+10.6f}")


# =============================================================================
# SUMMARY
# =============================================================================
t_total = timer.time() - t_start
print(f"\n{'='*78}")
print("  SUMMARY")
print(f"{'='*78}")
print(f"  Runtime: {t_total:.1f}s")
print(f"\n  Binding energies:")
for l in [0, 1, 2]:
    w = np.cos(l * pi / d_dim)
    print(f"    l={l} ({['sigma','pi','delta'][l]}): E_bind = {bind_E[l]:+.6f}, "
          f"channel weight = {w:.4f}")

print(f"\n  Cross-channel screening (1-S):")
print(f"    Measured vs GWT product rule:")
print(f"    (Results from section B above)")

print(f"\n  KEY QUESTION: Does l=2 anti-screen l=0?")
print(f"    GWT predicts: cos(2*pi/d) * cos(0) = {np.cos(2*pi/d_dim):.4f} (ANTI-SCREENING)")
print(f"    Simulation result: see section B and C above")
