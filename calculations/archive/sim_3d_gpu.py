"""
3D Sine-Gordon on Cubic Lattice — GPU
======================================
Simulate the GWT Lagrangian on the actual d=3 lattice:
  L = sum_<i,j> [(1/2)(phi_i - phi_j)^2 + (1/pi^2)(1 - cos(pi*phi_i))]

Breather modes have REAL angular momentum on the cubic lattice.
The t2g/eg split, angular parity, and three-body coupling should
emerge from the wave dynamics — no rules imposed.

Uses CuPy (CUDA) on RTX 4070 Ti.
"""

import numpy as np
from math import factorial
import time

try:
    import cupy as cp
    GPU = True
    print("GPU: CuPy available, using CUDA")
except ImportError:
    cp = np
    GPU = False
    print("GPU: CuPy not available, falling back to NumPy (slow)")

d = 3
gamma_sg = np.pi / (2**(d+1) * np.pi - 2)
V_0 = 1.0 / np.pi**2
xp = cp  # use gpu arrays

# === GRID ===
N = 64  # 64^3 = 262,144 points. Fits in GPU memory.
L = 8.0  # smaller box, better resolution near center
dx = 2 * L / N
x1d = np.linspace(-L, L, N, endpoint=False)
X, Y, Z_grid = np.meshgrid(x1d, x1d, x1d, indexing='ij')
R = np.sqrt(X**2 + Y**2 + Z_grid**2) + 1e-10

# Transfer to GPU
X_g = xp.asarray(X)
Y_g = xp.asarray(Y)
Z_g = xp.asarray(Z_grid)
R_g = xp.asarray(R)

dt = 0.3 * dx  # CFL for 3D: dt < dx/sqrt(3)


def laplacian_3d(phi):
    """6-point stencil Laplacian on periodic cubic lattice."""
    return (
        xp.roll(phi, 1, 0) + xp.roll(phi, -1, 0) +
        xp.roll(phi, 1, 1) + xp.roll(phi, -1, 1) +
        xp.roll(phi, 1, 2) + xp.roll(phi, -1, 2) -
        6 * phi
    ) / dx**2


def kink_3d(R, Z_eff):
    """Spherical kink (hedgehog) in 3D.
    The kink wraps around the origin with charge Z_eff.
    """
    return (4.0 / np.pi) * xp.arctan(xp.exp(-xp.sqrt(xp.asarray(Z_eff, dtype=float)) * R))


def breather_angular(X, Y, Z, R, l, m, omega, amplitude=0.05):
    """Breather mode with angular momentum (l, m) on cubic lattice.

    Angular patterns (real spherical harmonics):
      l=0: uniform (s-wave)
      l=1, m=0: z/r   (p_z)
      l=1, m=1: x/r   (p_x)
      l=1, m=-1: y/r  (p_y)
      l=2, m=0: (3z^2-r^2)/r^2  (d_z2, eg)
      l=2, m=1: xz/r^2           (d_xz, t2g)
      l=2, m=-1: yz/r^2          (d_yz, t2g)
      l=2, m=2: xy/r^2           (d_xy, t2g)
      l=2, m=-2: (x^2-y^2)/r^2  (d_x2y2, eg)

    Returns the angular pattern * radial breather envelope.
    """
    eps = float(np.sqrt(max(1.0 - omega**2, 1e-12)))

    # Angular part
    if l == 0:
        ang = xp.ones_like(R)
    elif l == 1:
        if m == 0: ang = Z / R
        elif m == 1: ang = X / R
        elif m == -1: ang = Y / R
        else: ang = xp.ones_like(R)
    elif l == 2:
        if m == 0: ang = (3*Z**2 - R**2) / R**2       # eg: d_z2
        elif m == 1: ang = X*Z / R**2                    # t2g: d_xz
        elif m == -1: ang = Y*Z / R**2                   # t2g: d_yz
        elif m == 2: ang = X*Y / R**2                    # t2g: d_xy
        elif m == -2: ang = (X**2 - Y**2) / R**2        # eg: d_x2-y2
        else: ang = xp.ones_like(R)
    elif l == 3:
        # f-modes: just use z^3/r^3 as representative
        if m == 0: ang = Z * (5*Z**2 - 3*R**2) / R**3
        else: ang = X * Y * Z / R**3  # representative mixed f
    else:
        ang = xp.ones_like(R)

    # Radial breather envelope
    radial = amplitude / (omega * xp.cosh(eps * R) + 1e-10)

    return (4.0 / np.pi) * xp.arctan(eps * ang * radial)


def total_energy_3d(phi, Z_eff):
    """Compute total energy of field configuration."""
    # Kinetic: skip (we measure at equilibrium)
    # Gradient energy: sum of (phi_i - phi_j)^2 / 2 for all neighbors
    GE = 0.5 * (
        (xp.roll(phi, 1, 0) - phi)**2 +
        (xp.roll(phi, 1, 1) - phi)**2 +
        (xp.roll(phi, 1, 2) - phi)**2
    )
    # Potential energy: Z * V_0 * (1 - cos(pi*phi))
    PE = Z_eff * V_0 * (1 - xp.cos(np.pi * phi))
    return float(xp.sum(GE + PE) * dx**3)


def evolve_3d(phi_init, Z_eff, n_steps=5000):
    """Evolve 3D sine-Gordon and return time-averaged energy."""
    phi = phi_init.copy()
    phi_old = phi.copy()
    energies = []

    for step in range(n_steps):
        lap = laplacian_3d(phi)
        force = Z_eff * (1.0 / np.pi) * xp.sin(np.pi * phi)
        phi_new = 2*phi - phi_old + dt**2 * (lap - force)
        phi_old = phi.copy()
        phi = phi_new

        if step > n_steps // 2 and step % 100 == 0:
            E = total_energy_3d(phi, Z_eff)
            energies.append(E)

    return np.mean(energies) if energies else 0.0


# === RUN ===
print(f"\n3D Sine-Gordon Simulation on {N}^3 Cubic Lattice")
print(f"=" * 55)
print(f"  Grid: {N}^3 = {N**3} points, L={L}, dx={dx:.3f}")
print(f"  dt = {dt:.4f}, GPU = {GPU}")
print()

Z_test = 5.0
omega_s = float(np.cos(1 * gamma_sg))
omega_p = float(np.cos(2 * gamma_sg))
omega_d = float(np.cos(3 * gamma_sg))

# Step 1: Kink alone
print("Step 1: Kink energy...", flush=True)
t0 = time.time()
phi_kink = kink_3d(R_g, Z_test)
E_kink = evolve_3d(phi_kink, Z_test, n_steps=3000)
print(f"  E_kink = {E_kink:.4f}  ({time.time()-t0:.1f}s)")

# Step 2: Kink + s-mode (l=0)
print("\nStep 2: s-mode (l=0) binding...", flush=True)
t0 = time.time()
phi_ks = phi_kink + breather_angular(X_g, Y_g, Z_g, R_g, l=0, m=0, omega=omega_s, amplitude=0.1)
E_ks = evolve_3d(phi_ks, Z_test, n_steps=3000)
E_s_bind = E_ks - E_kink
print(f"  E_bind(s) = {E_s_bind:.4f}  ({time.time()-t0:.1f}s)")

# Step 3: Kink + p-mode (l=1, m=0)
print("\nStep 3: p-mode (l=1) binding...", flush=True)
t0 = time.time()
phi_kp = phi_kink + breather_angular(X_g, Y_g, Z_g, R_g, l=1, m=0, omega=omega_p, amplitude=0.1)
E_kp = evolve_3d(phi_kp, Z_test, n_steps=3000)
E_p_bind = E_kp - E_kink
print(f"  E_bind(p) = {E_p_bind:.4f}  ({time.time()-t0:.1f}s)")

# Step 4: Kink + d-mode (t2g: l=2, m=2 = d_xy)
print("\nStep 4: d-mode t2g (l=2, m=2, d_xy) binding...", flush=True)
t0 = time.time()
phi_kd_t2g = phi_kink + breather_angular(X_g, Y_g, Z_g, R_g, l=2, m=2, omega=omega_d, amplitude=0.1)
E_kd_t2g = evolve_3d(phi_kd_t2g, Z_test, n_steps=3000)
E_d_t2g_bind = E_kd_t2g - E_kink
print(f"  E_bind(d_xy/t2g) = {E_d_t2g_bind:.4f}  ({time.time()-t0:.1f}s)")

# Step 5: Kink + d-mode (eg: l=2, m=0 = d_z2)
print("\nStep 5: d-mode eg (l=2, m=0, d_z2) binding...", flush=True)
t0 = time.time()
phi_kd_eg = phi_kink + breather_angular(X_g, Y_g, Z_g, R_g, l=2, m=0, omega=omega_d, amplitude=0.1)
E_kd_eg = evolve_3d(phi_kd_eg, Z_test, n_steps=3000)
E_d_eg_bind = E_kd_eg - E_kink
print(f"  E_bind(d_z2/eg) = {E_d_eg_bind:.4f}  ({time.time()-t0:.1f}s)")

# Step 6: Coupling via simultaneous modes
# Place BOTH modes at once. The interaction energy is:
# E_int = E(kink+A+B) - E(kink+A) - E(kink+B) + E(kink)
# Positive E_int = repulsion. Negative = attraction.

print("\nStep 6: d(t2g) + s coupling (interaction energy)...", flush=True)
t0 = time.time()
phi_ds = phi_kink.copy()
phi_ds = phi_ds + breather_angular(X_g, Y_g, Z_g, R_g, l=2, m=2, omega=omega_d, amplitude=0.1)
phi_ds = phi_ds + breather_angular(X_g, Y_g, Z_g, R_g, l=0, m=0, omega=omega_s, amplitude=0.1)
E_ds = evolve_3d(phi_ds, Z_test, n_steps=3000)
E_int_ds = E_ds - E_kd_t2g - E_ks + E_kink  # interaction energy
C_d_s = E_int_ds / E_s_bind if abs(E_s_bind) > 1e-10 else 0
print(f"  E_int(d_t2g, s) = {E_int_ds:.6f}")
print(f"  C(d_t2g, s) = E_int / E_s = {C_d_s:.4f}  ({time.time()-t0:.1f}s)")

print("\nStep 7: d(t2g) + p coupling (interaction energy)...", flush=True)
t0 = time.time()
phi_dp = phi_kink.copy()
phi_dp = phi_dp + breather_angular(X_g, Y_g, Z_g, R_g, l=2, m=2, omega=omega_d, amplitude=0.1)
phi_dp = phi_dp + breather_angular(X_g, Y_g, Z_g, R_g, l=1, m=0, omega=omega_p, amplitude=0.1)
E_dp = evolve_3d(phi_dp, Z_test, n_steps=3000)
E_int_dp = E_dp - E_kd_t2g - E_kp + E_kink
C_d_p = E_int_dp / E_p_bind if abs(E_p_bind) > 1e-10 else 0
print(f"  E_int(d_t2g, p) = {E_int_dp:.6f}")
print(f"  C(d_t2g, p) = E_int / E_p = {C_d_p:.4f}  ({time.time()-t0:.1f}s)")

# Step 8: d(eg) + s and d(eg) + p for comparison
print("\nStep 8: d(eg) + s coupling...", flush=True)
t0 = time.time()
phi_eg_s = phi_kink.copy()
phi_eg_s = phi_eg_s + breather_angular(X_g, Y_g, Z_g, R_g, l=2, m=0, omega=omega_d, amplitude=0.1)
phi_eg_s = phi_eg_s + breather_angular(X_g, Y_g, Z_g, R_g, l=0, m=0, omega=omega_s, amplitude=0.1)
E_eg_s = evolve_3d(phi_eg_s, Z_test, n_steps=3000)
E_int_eg_s = E_eg_s - E_kd_eg - E_ks + E_kink
C_eg_s = E_int_eg_s / E_s_bind if abs(E_s_bind) > 1e-10 else 0
print(f"  E_int(d_eg, s) = {E_int_eg_s:.6f}")
print(f"  C(d_eg, s) = {C_eg_s:.4f}  ({time.time()-t0:.1f}s)")

print("\nStep 9: d(eg) + p coupling...", flush=True)
t0 = time.time()
phi_eg_p = phi_kink.copy()
phi_eg_p = phi_eg_p + breather_angular(X_g, Y_g, Z_g, R_g, l=2, m=0, omega=omega_d, amplitude=0.1)
phi_eg_p = phi_eg_p + breather_angular(X_g, Y_g, Z_g, R_g, l=1, m=0, omega=omega_p, amplitude=0.1)
E_eg_p = evolve_3d(phi_eg_p, Z_test, n_steps=3000)
E_int_eg_p = E_eg_p - E_kd_eg - E_kp + E_kink
C_eg_p = E_int_eg_p / E_p_bind if abs(E_p_bind) > 1e-10 else 0
print(f"  E_int(d_eg, p) = {E_int_eg_p:.6f}")
print(f"  C(d_eg, p) = {C_eg_p:.4f}  ({time.time()-t0:.1f}s)")

# Step 10: s + p coupling for reference
print("\nStep 10: s + p coupling...", flush=True)
t0 = time.time()
phi_sp = phi_kink.copy()
phi_sp = phi_sp + breather_angular(X_g, Y_g, Z_g, R_g, l=0, m=0, omega=omega_s, amplitude=0.1)
phi_sp = phi_sp + breather_angular(X_g, Y_g, Z_g, R_g, l=1, m=0, omega=omega_p, amplitude=0.1)
E_sp = evolve_3d(phi_sp, Z_test, n_steps=3000)
E_int_sp = E_sp - E_ks - E_kp + E_kink
C_s_p = E_int_sp / E_p_bind if abs(E_p_bind) > 1e-10 else 0
print(f"  E_int(s, p) = {E_int_sp:.6f}")
print(f"  C(s, p) = {C_s_p:.4f}  ({time.time()-t0:.1f}s)")

# Results
print("\n" + "=" * 55)
print("RESULTS:")
print(f"  s-mode binding:     {E_s_bind:.4f}")
print(f"  p-mode binding:     {E_p_bind:.4f}")
print(f"  d(t2g) binding:     {E_d_t2g_bind:.4f}")
print(f"  d(eg) binding:      {E_d_eg_bind:.4f}")
print(f"  t2g/eg ratio:       {E_d_t2g_bind/E_d_eg_bind:.4f}")
print(f"    (v19: t2g couples at w_delta, eg at w_delta/d)")
print(f"    (expect ratio = d = {d} if t2g is d times stronger)")
print()
print(f"\n  COUPLING MATRIX (interaction energy / binding energy):")
print(f"    t2g + s: {C_d_s:+.4f}")
print(f"    t2g + p: {C_d_p:+.4f}")
print(f"    eg  + s: {C_eg_s:+.4f}")
print(f"    eg  + p: {C_eg_p:+.4f}")
print(f"    s   + p: {C_s_p:+.4f}")
print(f"\n  t2g->s / t2g->p: {C_d_s/C_d_p:.4f}" if abs(C_d_p) > 1e-6 else "")
print(f"  eg->s / eg->p:   {C_eg_s/C_eg_p:.4f}" if abs(C_eg_p) > 1e-6 else "")
print(f"  t2g / eg (on s): {C_d_s/C_eg_s:.4f}" if abs(C_eg_s) > 1e-6 else "")
print(f"  t2g / eg (on p): {C_d_p/C_eg_p:.4f}" if abs(C_eg_p) > 1e-6 else "")
