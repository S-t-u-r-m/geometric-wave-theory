#!/usr/bin/env python3
"""
f->d wave coupling: measure angular momentum coupling on 2D lattice.

Fix from v1: add centrifugal barrier m^2/(2r^2) to preserve angular
momentum during imaginary-time relaxation. Without it, all modes
collapse to l=0 ground state.

In 2D, angular momentum m gives centrifugal barrier m^2/(2r^2).
Each mode (l,m) has a different effective potential.

We also project onto the correct angular sector after each step
to prevent mode mixing during relaxation.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import time

try:
    import cupy as cp
    xp = cp
    GPU = True
    print("GPU: CuPy available")
except ImportError:
    xp = np
    GPU = False
    print("GPU: using CPU")

# === 2D LATTICE ===
Nx = 256
Ny = 256
dx = 0.1
dt = 0.001
n_steps = 60000

x = (xp.arange(Nx) - Nx//2) * dx
y = (xp.arange(Ny) - Ny//2) * dx
X, Y = xp.meshgrid(x, y, indexing='ij')
R = xp.sqrt(X**2 + Y**2) + 1e-10
THETA = xp.arctan2(Y, X)

d = 3
Z_eff = 3.0
r_soft = 2 * dx  # softening radius

# Base kink potential
V_kink = -Z_eff / xp.sqrt(R**2 + r_soft**2)

print(f"\n2D lattice: {Nx}x{Ny}, dx={dx}")
print(f"Z_eff={Z_eff}, dt={dt}, {n_steps} steps")

def laplacian(psi, dx):
    lap = xp.zeros_like(psi)
    lap[1:-1, :] += (psi[2:, :] - 2*psi[1:-1, :] + psi[:-2, :]) / dx**2
    lap[:, 1:-1] += (psi[:, 2:] - 2*psi[:, 1:-1] + psi[:, :-2]) / dx**2
    return lap

def normalize(psi, dx):
    return psi / (xp.sqrt(xp.sum(psi**2) * dx**2) + 1e-30)

def project_angular(psi, l, THETA):
    """Project wavefunction onto angular momentum sector l.
    Keep only the cos(l*theta) component."""
    if l == 0:
        # Average over angle (keep isotropic part)
        return psi  # for l=0, no projection needed in practice
    # Multiply by cos(l*theta), average to get coefficient, multiply back
    ang = xp.cos(l * THETA)
    # The coefficient is 2 * <psi * cos(l*theta)> (Fourier)
    coeff_field = psi * ang
    # Reconstruct
    return coeff_field * ang * 2  # approximate projection

def energy(psi, V, dx):
    lap = laplacian(psi, dx)
    KE = float(-0.5 * xp.sum(psi * lap) * dx**2)
    PE = float(xp.sum(psi * V * psi) * dx**2)
    return KE + PE

def relax_mode(l, V_extra, n_e_label="", steps=None):
    """Find ground state of angular momentum mode l in potential V_kink + V_extra.
    Uses centrifugal barrier + angular projection to maintain l."""
    if steps is None:
        steps = n_steps

    # Centrifugal barrier: m^2 / (2*r^2) in 2D
    V_cent = l**2 / (2 * xp.maximum(R**2, r_soft**2))
    V_total = V_kink + V_cent + V_extra

    # Initial guess with correct angular symmetry
    ang = xp.cos(l * THETA) if l > 0 else xp.ones_like(THETA)
    rad = R**max(l, 0.5) * xp.exp(-Z_eff * R / (2 * (l + 1.5)))
    psi = ang * rad
    psi = normalize(psi, dx)

    t0 = time.time()
    for step in range(steps):
        lap = laplacian(psi, dx)
        Hpsi = -0.5 * lap + V_total * psi
        psi = psi - dt * Hpsi

        # Project onto angular sector every 100 steps to prevent drift
        if l > 0 and step % 100 == 0:
            psi = project_angular(psi, l, THETA)

        psi = normalize(psi, dx)

    E = energy(psi, V_total, dx)
    elapsed = time.time() - t0

    # Measure radial extent
    rho = psi**2
    r_mean = float(xp.sum(rho * R * dx**2) / xp.sum(rho * dx**2))

    label = f"l={l}" + (f" + {n_e_label}" if n_e_label else "")
    print(f"  {label:20s}: E = {E:.6f}, <r> = {r_mean:.3f} ({elapsed:.1f}s)")
    return psi, E, r_mean

def screening_from_mode(psi_screen, n_electrons):
    """Create screening potential from occupied mode."""
    rho = psi_screen**2 * n_electrons
    # Solve Poisson equation for 2D: nabla^2 V = -4*pi*rho
    # Use iterative relaxation (simple Gauss-Seidel)
    V = xp.zeros_like(psi_screen)
    for _ in range(2000):
        V[1:-1, 1:-1] = 0.25 * (
            V[2:, 1:-1] + V[:-2, 1:-1] + V[1:-1, 2:] + V[1:-1, :-2]
            + dx**2 * 4 * np.pi * rho[1:-1, 1:-1]
        )
    return V

# === EXPERIMENTS ===

# 1. Bare modes at each l
print(f"\n{'='*60}")
print("Bare mode energies (no screening)")
print("=" * 60)
V_zero = xp.zeros_like(R)
psi_s, E_s, r_s = relax_mode(0, V_zero, steps=n_steps)
psi_d, E_d, r_d = relax_mode(2, V_zero, steps=n_steps)
psi_f, E_f, r_f = relax_mode(3, V_zero, steps=n_steps)
psi_p, E_p, r_p = relax_mode(1, V_zero, steps=n_steps)

print(f"\n  Energies: s={E_s:.4f}, p={E_p:.4f}, d={E_d:.4f}, f={E_f:.4f}")
print(f"  Radii:   s={r_s:.3f}, p={r_p:.3f}, d={r_d:.3f}, f={r_f:.3f}")
print(f"  Centrifugal ordering: {'s < p < d < f' if r_s < r_p < r_d < r_f else 'UNEXPECTED'}")

# 2. Create screening potentials
print(f"\n{'='*60}")
print("Creating screening potentials...")
print("=" * 60)
t0 = time.time()
V_from_s = screening_from_mode(psi_s, 2)
V_from_p = screening_from_mode(psi_p, 6)
V_from_d_inner = screening_from_mode(psi_d, 10)
V_from_f = screening_from_mode(psi_f, 2)
V_from_f7 = screening_from_mode(psi_f, 7)
V_from_f14 = screening_from_mode(psi_f, 14)
print(f"  Done ({time.time()-t0:.1f}s)")

# 3. d-mode with different screening shells
print(f"\n{'='*60}")
print("d-mode (l=2) with various screening shells")
print("=" * 60)
_, E_d_s2, _ = relax_mode(2, V_from_s, "2s")
_, E_d_p6, _ = relax_mode(2, V_from_p, "6p")
_, E_d_d10, _ = relax_mode(2, V_from_d_inner, "10d")
_, E_d_f2, _ = relax_mode(2, V_from_f, "2f")
_, E_d_f7, _ = relax_mode(2, V_from_f7, "7f")
_, E_d_f14, _ = relax_mode(2, V_from_f14, "14f")

# 4. s-mode with same screening (for comparison)
print(f"\n{'='*60}")
print("s-mode (l=0) with f-screening (comparison)")
print("=" * 60)
_, E_s_f2, _ = relax_mode(0, V_from_f, "2f")
_, E_s_f14, _ = relax_mode(0, V_from_f14, "14f")

# === ANALYSIS ===
print(f"\n{'='*60}")
print("COUPLING ANALYSIS")
print("=" * 60)

dE_d_s = E_d_s2 - E_d
dE_d_p = E_d_p6 - E_d
dE_d_d = E_d_d10 - E_d
dE_d_f2 = E_d_f2 - E_d
dE_d_f7 = E_d_f7 - E_d
dE_d_f14 = E_d_f14 - E_d

dE_s_f2 = E_s_f2 - E_s
dE_s_f14 = E_s_f14 - E_s

print(f"\n  Energy shifts on d-mode (l=2):")
print(f"    + 2s:  dE = {dE_d_s:+.6f}  per e: {dE_d_s/2:+.6f}")
print(f"    + 6p:  dE = {dE_d_p:+.6f}  per e: {dE_d_p/6:+.6f}")
print(f"    + 10d: dE = {dE_d_d:+.6f}  per e: {dE_d_d/10:+.6f}")
print(f"    + 2f:  dE = {dE_d_f2:+.6f}  per e: {dE_d_f2/2:+.6f}")
print(f"    + 7f:  dE = {dE_d_f7:+.6f}  per e: {dE_d_f7/7:+.6f}")
print(f"    + 14f: dE = {dE_d_f14:+.6f}  per e: {dE_d_f14/14:+.6f}")

print(f"\n  Energy shifts on s-mode (l=0) from f:")
print(f"    + 2f:  dE = {dE_s_f2:+.6f}  per e: {dE_s_f2/2:+.6f}")
print(f"    + 14f: dE = {dE_s_f14:+.6f}  per e: {dE_s_f14/14:+.6f}")

print(f"\n  SCREENING vs ANTI-SCREENING:")
print(f"  (positive = screening/shielding, negative = anti-screening)")

# Coupling ratios relative to s->d
w_per_e = {
    's->d': dE_d_s/2,
    'p->d': dE_d_p/6,
    'd->d': dE_d_d/10,
    'f->d (2e)': dE_d_f2/2,
    'f->d (7e)': dE_d_f7/7,
    'f->d (14e)': dE_d_f14/14,
    'f->s (2e)': dE_s_f2/2,
    'f->s (14e)': dE_s_f14/14,
}

ref = w_per_e['s->d']
print(f"\n  Per-electron coupling relative to s->d:")
for name, w in w_per_e.items():
    ratio = w/ref if abs(ref) > 1e-10 else 0
    print(f"    {name:15s}: {w:+.6f} = {ratio:+.4f} * w_s")

print(f"\n  GWT model predictions:")
print(f"    s->d: w_pi = +0.5 (screen)  -> ratio = 1.0")
print(f"    d->d: w_delta = -0.5 (anti-screen) -> ratio = -1.0")
print(f"    f->d: w_f_d = -1/6 (anti-screen) -> ratio = -0.333")
print(f"    f->s: w_f_s = -1/6 (anti-screen) -> ratio = -0.333")

# KEY: is f->d DIFFERENT from f->s?
if abs(ref) > 1e-10:
    r_fd = w_per_e['f->d (2e)'] / ref
    r_fs = w_per_e['f->s (2e)'] / ref
    print(f"\n  KEY COMPARISON:")
    print(f"    f->d ratio: {r_fd:+.4f}")
    print(f"    f->s ratio: {r_fs:+.4f}")
    if abs(r_fd - r_fs) > 0.05:
        print(f"    DIFFERENT! f couples to d and s differently!")
        print(f"    This could explain Lu's anomaly.")
    else:
        print(f"    Similar. f couples equally to d and s.")
