"""
9-Channel Bond Simulator — Oh Tensor Product Decomposition
============================================================
Two breathers interact on the d=3 cubic lattice.
The interaction energy is decomposed into 9 Oh channels:

  T1u x T1u = A1g(1) + Eg(2) + T1g(3) + T2g(3)
              sigma     shape   rotation  shear

Each channel has a physical role:
  A1g: sigma bond (scalar coupling)
  Eg:  LP repulsion (shape distortion)
  T1g: bond angle rigidity (rotation)
  T2g: pi enhancement (off-diagonal shear)

The NET energy across all channels = the bond energy.
No analytical formula needed — just physics on a lattice.

Uses CuPy (CUDA) on RTX 4070 Ti.
"""

import numpy as np
from math import factorial
import time

try:
    import cupy as cp
    xp = cp
    print("GPU: CuPy active")
except ImportError:
    xp = np
    print("CPU: CuPy not available")

d = 3
gamma_sg = np.pi / (2**(d+1) * np.pi - 2)
V_0 = 1.0 / np.pi**2

# Grid
N = 80
BOX = 10.0
dx = 2 * BOX / N
dt = 0.2 * dx

x1d = xp.linspace(-BOX, BOX, N, endpoint=False, dtype=np.float64)
X, Y, Z = xp.meshgrid(x1d, x1d, x1d, indexing='ij')


def laplacian(phi):
    return (xp.roll(phi,1,0) + xp.roll(phi,-1,0) +
            xp.roll(phi,1,1) + xp.roll(phi,-1,1) +
            xp.roll(phi,1,2) + xp.roll(phi,-1,2) - 6*phi) / dx**2


def energy(phi):
    """Total energy of a field configuration."""
    GE = 0.5 * ((xp.roll(phi,1,0)-phi)**2 +
                 (xp.roll(phi,1,1)-phi)**2 +
                 (xp.roll(phi,1,2)-phi)**2)
    PE = V_0 * (1 - xp.cos(np.pi * phi))
    return float(xp.sum(GE + PE) * dx**3)


def evolve(phi0, n_steps=3000, damping=0.02):
    """Damped evolution to find equilibrium energy."""
    phi = phi0.copy()
    phi_old = phi.copy()
    energies = []
    for step in range(n_steps):
        lap = laplacian(phi)
        force = (1.0/np.pi) * xp.sin(np.pi * phi)
        phi_new = (2 - damping*dt)*phi - (1 - damping*dt)*phi_old + dt**2*(lap - force)
        phi_old = phi.copy()
        phi = phi_new
        if step > n_steps * 3 // 4 and step % 50 == 0:
            energies.append(energy(phi))
    return np.mean(energies) if energies else energy(phi)


def make_breather(center_z, angular_type, amplitude, omega=None, eps=None):
    """Create a breather mode at position center_z along the z-axis.

    angular_type: 's', 'pz', 'px', 'py', 'dz2', 'dxz', 'dyz', 'dxy', 'dx2y2'
    """
    if omega is None:
        omega = float(np.cos(2 * gamma_sg))  # p-mode default
    if eps is None:
        eps = float(np.sqrt(max(1 - omega**2, 1e-12)))

    Zs = Z - center_z
    R = xp.sqrt(X**2 + Y**2 + Zs**2) + 1e-10

    # Angular part
    ang_map = {
        's':     lambda: xp.ones_like(R),
        'pz':    lambda: Zs / R,
        'px':    lambda: X / R,
        'py':    lambda: Y / R,
        'dz2':   lambda: (3*Zs**2 - R**2) / R**2,
        'dxz':   lambda: X*Zs / R**2,
        'dyz':   lambda: Y*Zs / R**2,
        'dxy':   lambda: X*Y / R**2,
        'dx2y2': lambda: (X**2 - Y**2) / R**2,
    }
    ang = ang_map.get(angular_type, ang_map['s'])()

    radial = amplitude / (omega * xp.cosh(eps * R) + 1e-10)
    return (4.0 / np.pi) * xp.arctan(eps * ang * radial)


def interaction_energy(phi_a, phi_b):
    """Compute interaction energy: E(A+B) - E(A) - E(B)."""
    E_a = evolve(phi_a)
    E_b = evolve(phi_b)
    E_ab = evolve(phi_a + phi_b)
    return E_ab - E_a - E_b


def measure_9_channels(R_sep, amplitude=1.0):
    """Measure all 9 Oh channels for two breathers at separation R_sep.

    The 9 channels of T1u x T1u:
      A1g(1): sigma (pz x pz along bond)
      Eg(2):  shape (pz x dz2, pz x dx2y2)
      T1g(3): rotation (pz x px cross-terms that give angular momentum)
      T2g(3): shear (px x px, py x py perpendicular same-type)

    We measure each by choosing appropriate angular mode combinations.
    """
    results = {}

    print(f"\n9-CHANNEL MEASUREMENT at R = {R_sep}")
    print(f"{'Channel':>12} {'Type':>6} {'Modes':>12} {'V_int':>12}")
    print("-" * 48)

    # A1g channel: sigma bond (pz + pz along bond axis)
    phi_a = make_breather(-R_sep/2, 'pz', amplitude)
    phi_b = make_breather(+R_sep/2, 'pz', amplitude)
    V_sigma = interaction_energy(phi_a, phi_b)
    results['A1g'] = V_sigma
    print(f"{'A1g':>12} {'sigma':>6} {'pz+pz':>12} {V_sigma:+12.6f}")

    # T2g channels: pi-type (perpendicular same-mode coupling)
    # px+px (both perpendicular to bond axis)
    phi_a = make_breather(-R_sep/2, 'px', amplitude)
    phi_b = make_breather(+R_sep/2, 'px', amplitude)
    V_pi_x = interaction_energy(phi_a, phi_b)
    results['T2g_xx'] = V_pi_x
    print(f"{'T2g':>12} {'pi_x':>6} {'px+px':>12} {V_pi_x:+12.6f}")

    # py+py
    phi_a = make_breather(-R_sep/2, 'py', amplitude)
    phi_b = make_breather(+R_sep/2, 'py', amplitude)
    V_pi_y = interaction_energy(phi_a, phi_b)
    results['T2g_yy'] = V_pi_y
    print(f"{'T2g':>12} {'pi_y':>6} {'py+py':>12} {V_pi_y:+12.6f}")

    # T1g channel: cross-coupling (pz x px = angular momentum)
    phi_a = make_breather(-R_sep/2, 'pz', amplitude)
    phi_b = make_breather(+R_sep/2, 'px', amplitude)
    V_cross = interaction_energy(phi_a, phi_b)
    results['T1g_zx'] = V_cross
    print(f"{'T1g':>12} {'cross':>6} {'pz+px':>12} {V_cross:+12.6f}")

    # pz x py
    phi_a = make_breather(-R_sep/2, 'pz', amplitude)
    phi_b = make_breather(+R_sep/2, 'py', amplitude)
    V_cross2 = interaction_energy(phi_a, phi_b)
    results['T1g_zy'] = V_cross2
    print(f"{'T1g':>12} {'cross':>6} {'pz+py':>12} {V_cross2:+12.6f}")

    # px x py
    phi_a = make_breather(-R_sep/2, 'px', amplitude)
    phi_b = make_breather(+R_sep/2, 'py', amplitude)
    V_cross3 = interaction_energy(phi_a, phi_b)
    results['T1g_xy'] = V_cross3
    print(f"{'T1g':>12} {'cross':>6} {'px+py':>12} {V_cross3:+12.6f}")

    # Eg channels: shape coupling (s x s for A1g comparison)
    phi_a = make_breather(-R_sep/2, 's', amplitude)
    phi_b = make_breather(+R_sep/2, 's', amplitude)
    V_ss = interaction_energy(phi_a, phi_b)
    results['s_s'] = V_ss
    print(f"{'(ref)':>12} {'s+s':>6} {'s+s':>12} {V_ss:+12.6f}")

    # FULL ORBITAL (out-of-phase) pz — LP repulsion test
    # Approximate by doubling amplitude (paired mode)
    phi_a = make_breather(-R_sep/2, 'px', amplitude*2)
    phi_b = make_breather(+R_sep/2, 'px', amplitude*2)
    V_LP = interaction_energy(phi_a, phi_b)
    results['LP_px'] = V_LP
    print(f"{'LP':>12} {'full':>6} {'2px+2px':>12} {V_LP:+12.6f}")

    # Summary
    print()
    print("CHANNEL SUMMARY:")
    print(f"  Sigma (A1g):  {V_sigma:+.6f}  {'BONDING' if V_sigma < 0 else 'repulsive'}")
    print(f"  Pi avg (T2g): {(V_pi_x+V_pi_y)/2:+.6f}  {'BONDING' if (V_pi_x+V_pi_y)/2 < 0 else 'repulsive'}")
    print(f"  Cross (T1g):  {(V_cross+V_cross2+V_cross3)/3:+.6f}  {'rotation' if abs((V_cross+V_cross2+V_cross3)/3) < abs(V_sigma)/10 else 'significant'}")
    print(f"  LP (full):    {V_LP:+.6f}  {'REPULSIVE' if V_LP > 0 else 'attractive'}")

    if abs(V_sigma) > 1e-10:
        print(f"\n  Pi/Sigma ratio: {(V_pi_x+V_pi_y)/(2*V_sigma):.3f}  (Oh predicts: {np.cos(np.pi/d):.3f})")
        print(f"  Cross/Sigma:    {abs((V_cross+V_cross2+V_cross3)/(3*V_sigma)):.3f}  (Oh predicts: ~0)")
        print(f"  LP/Sigma:       {V_LP/V_sigma:.3f}")

    return results


if __name__ == "__main__":
    print("9-Channel Bond Simulator")
    print("=" * 50)
    print(f"Grid: {N}^3, dx={dx:.3f}, BOX={BOX}")
    print()

    # Scan separations
    for R in [2.0, 3.0, 4.0]:
        t0 = time.time()
        results = measure_9_channels(R, amplitude=1.0)
        print(f"\n  Time: {time.time()-t0:.0f}s")
        print()
