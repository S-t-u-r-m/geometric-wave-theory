"""
3D Breather Interaction on Sine-Gordon Lattice
================================================
Two spherical standing wave patterns on a 3D cubic lattice.
Measure interaction energy vs separation to see if angular
interference produces oscillatory sin(kR) behavior.

Lattice equation of motion (3D):
  d^2 phi/dt^2 = Laplacian(phi) - V_0 * pi * sin(pi * phi)

where Laplacian uses 6 nearest neighbors on cubic lattice.

Breather in 3D: spherically symmetric localized oscillation.
  phi(r,t) = f(r) * cos(omega * t)
where f(r) = (4/pi) * arctan[ eta / (w * cosh(eta * r)) ] / r_factor

For p-wave: phi ~ f(r) * cos(theta) (angular dependence)
"""

import numpy as np
import time as timer

pi = np.pi
V_0 = 1.0 / pi**2
N_br = 24
gamma_br = pi / (N_br + 1)


def run_3d_simulation(N_grid, L, n1, n2, R_sep, sign2=+1,
                       l1=0, l2=0, n_steps=2000, dt_factor=0.15):
    """Run a 3D breather interaction simulation.

    Parameters:
        N_grid: grid points per dimension
        L: box size
        n1, n2: breather modes
        R_sep: separation along x-axis
        sign2: +1 for in-phase, -1 for out-of-phase
        l1, l2: angular momentum (0=s-wave, 1=p-wave)
        n_steps: evolution steps
    """
    dx = L / N_grid
    dt = dt_factor * dx  # CFL condition for 3D: dt < dx/sqrt(3)

    # Create 3D grid
    coords = np.linspace(-L/2, L/2, N_grid)
    # Use meshgrid for vectorized computation
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing='ij')

    def breather_3d(X, Y, Z, n, center_x, amplitude, l_ang):
        """3D breather profile (spherically symmetric or p-wave)."""
        w = np.cos(n * gamma_br)
        eta = np.sin(n * gamma_br)

        rx = X - center_x
        r = np.sqrt(rx**2 + Y**2 + Z**2)
        r = np.maximum(r, dx/2)  # regularize at origin

        # Radial profile (breather envelope)
        profile = amplitude * (4.0/pi) * np.arctan(eta / (w * np.cosh(eta * r)))

        # For 3D: divide by r to get proper radial dependence
        # (3D standing wave ~ f(r)/r for s-wave, f(r)*cos(theta)/r for p-wave)
        # Actually, keep it simple: use the profile as a 3D scalar field
        # The 1/r factor is already implicit in the spherical spread

        if l_ang == 1:
            # p-wave: multiply by cos(theta) = x/r along bond axis
            cos_theta = rx / r
            profile = profile * cos_theta

        return profile

    # Single breather energies
    def compute_energy(phi):
        """Total energy of 3D field."""
        # Gradient energy (6 nearest neighbors)
        gx = np.diff(phi, axis=0)
        gy = np.diff(phi, axis=1)
        gz = np.diff(phi, axis=2)
        E_grad = 0.5 * (np.sum(gx**2) + np.sum(gy**2) + np.sum(gz**2)) / dx * dx**2
        # Potential energy
        E_pot = V_0 * np.sum(1 - np.cos(pi * phi)) * dx**3
        return E_grad + E_pot

    def compute_force(phi):
        """3D lattice force = Laplacian - V_0*pi*sin(pi*phi)."""
        force = np.zeros_like(phi)
        # Laplacian via finite differences
        force[1:-1,:,:] += phi[:-2,:,:] + phi[2:,:,:] - 2*phi[1:-1,:,:]
        force[:,1:-1,:] += phi[:,:-2,:] + phi[:,2:,:] - 2*phi[:,1:-1,:]
        force[:,:,1:-1] += phi[:,:,:-2] + phi[:,:,2:] - 2*phi[:,:,1:-1]
        force /= dx**2
        # Potential
        force -= V_0 * pi * np.sin(pi * phi)
        return force

    # Set up single breather 1
    phi1_single = breather_3d(X, Y, Z, n1, 0.0, 1.0, l1)
    E1_init = compute_energy(phi1_single)

    # Evolve single breather 1
    phi = phi1_single.copy()
    dphi = np.zeros_like(phi)
    force = compute_force(phi)
    dphi += 0.5 * dt * force
    E1_samples = []
    for step in range(n_steps):
        phi += dt * dphi
        force = compute_force(phi)
        dphi += dt * force
        if step % 50 == 0 and step > n_steps//4:
            dphi_half = dphi - 0.5*dt*force
            E_kin = 0.5 * np.sum(dphi_half**2) * dx**3
            E_pot_grad = compute_energy(phi)
            E1_samples.append(E_kin + E_pot_grad)
    E1 = np.mean(E1_samples) if E1_samples else E1_init

    # Single breather 2
    phi2_single = breather_3d(X, Y, Z, n2, 0.0, sign2, l2)
    E2_init = compute_energy(phi2_single)

    phi = phi2_single.copy()
    dphi = np.zeros_like(phi)
    force = compute_force(phi)
    dphi += 0.5 * dt * force
    E2_samples = []
    for step in range(n_steps):
        phi += dt * dphi
        force = compute_force(phi)
        dphi += dt * force
        if step % 50 == 0 and step > n_steps//4:
            dphi_half = dphi - 0.5*dt*force
            E_kin = 0.5 * np.sum(dphi_half**2) * dx**3
            E_pot_grad = compute_energy(phi)
            E2_samples.append(E_kin + E_pot_grad)
    E2 = np.mean(E2_samples) if E2_samples else E2_init

    # Two-breather system
    phi_two = (breather_3d(X, Y, Z, n1, -R_sep/2, 1.0, l1) +
               breather_3d(X, Y, Z, n2, +R_sep/2, sign2, l2))

    phi = phi_two.copy()
    dphi = np.zeros_like(phi)
    force = compute_force(phi)
    dphi += 0.5 * dt * force
    E12_samples = []
    for step in range(n_steps):
        phi += dt * dphi
        force = compute_force(phi)
        dphi += dt * force
        if step % 50 == 0 and step > n_steps//4:
            dphi_half = dphi - 0.5*dt*force
            E_kin = 0.5 * np.sum(dphi_half**2) * dx**3
            E_pot_grad = compute_energy(phi)
            E12_samples.append(E_kin + E_pot_grad)
    E12 = np.mean(E12_samples) if E12_samples else compute_energy(phi_two)

    E_int = E12 - E1 - E2
    return E1, E2, E12, E_int


# =============================================================================
# MAIN
# =============================================================================

print("=" * 80)
print("  3D BREATHER INTERACTION SIMULATION")
print("=" * 80)

# Use moderate grid — 3D is expensive!
# Memory: N^3 * 8 bytes. N=60 -> 1.7 MB. N=80 -> 4 MB. N=100 -> 8 MB.
N = 60
L = 20.0
n_steps = 3000

print(f"  Grid: {N}x{N}x{N} = {N**3} points")
print(f"  Box: L={L}, dx={L/N:.3f}")
print(f"  Steps: {n_steps}")
print(f"  V_0 = {V_0:.6f}")

# --- s-wave (l=0) breathers ---
print(f"\n{'='*80}")
print(f"  s-WAVE (l=0) BREATHER INTERACTION")
print(f"{'='*80}")

for n1, n2, label in [(1,1,"mode 1+1"), (3,3,"mode 3+3"), (1,3,"mode 1+3")]:
    print(f"\n  --- ({n1},{n2}) {label}, s-wave ---")
    print(f"  {'R':>6} {'E1':>10} {'E2':>10} {'E12':>10} {'dE(+-)':>10} {'dE(++)':>10}")

    for R in [2.0, 3.0, 4.0, 5.0, 7.0, 10.0]:
        t0 = timer.time()
        E1, E2, E12_pm, dE_pm = run_3d_simulation(N, L, n1, n2, R, sign2=-1,
                                                     l1=0, l2=0, n_steps=n_steps)
        _, _, E12_pp, dE_pp = run_3d_simulation(N, L, n1, n2, R, sign2=+1,
                                                  l1=0, l2=0, n_steps=n_steps)
        elapsed = timer.time() - t0
        print(f"  {R:6.1f} {E1:10.4f} {E2:10.4f} {E12_pm:10.4f} {dE_pm:+10.4f} {dE_pp:+10.4f}  ({elapsed:.1f}s)")


# --- p-wave (l=1) breathers ---
print(f"\n{'='*80}")
print(f"  p-WAVE (l=1) BREATHER INTERACTION (along bond axis)")
print(f"{'='*80}")

for n1, n2, label in [(3,3,"mode 3+3 p-wave"), (1,3,"mode 1s+3p")]:
    l1_val = 1 if n1 >= 3 else 0
    l2_val = 1 if n2 >= 3 else 0
    print(f"\n  --- ({n1},{n2}) {label}, l=({l1_val},{l2_val}) ---")
    print(f"  {'R':>6} {'E1':>10} {'E2':>10} {'E12':>10} {'dE(+-)':>10} {'dE(++)':>10}")

    for R in [2.0, 3.0, 4.0, 5.0, 7.0, 10.0]:
        t0 = timer.time()
        E1, E2, E12_pm, dE_pm = run_3d_simulation(N, L, n1, n2, R, sign2=-1,
                                                     l1=l1_val, l2=l2_val, n_steps=n_steps)
        _, _, E12_pp, dE_pp = run_3d_simulation(N, L, n1, n2, R, sign2=+1,
                                                  l1=l1_val, l2=l2_val, n_steps=n_steps)
        elapsed = timer.time() - t0
        print(f"  {R:6.1f} {E1:10.4f} {E2:10.4f} {E12_pm:10.4f} {dE_pm:+10.4f} {dE_pp:+10.4f}  ({elapsed:.1f}s)")


# --- Mixed: s + p interaction ---
print(f"\n{'='*80}")
print(f"  MIXED s+p INTERACTION (like H_1s + C_2p)")
print(f"{'='*80}")

print(f"\n  --- mode 1(s) + mode 3(p) ---")
print(f"  {'R':>6} {'dE(+-)':>10} {'dE(++)':>10}")

for R in [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]:
    t0 = timer.time()
    _, _, _, dE_pm = run_3d_simulation(N, L, 1, 3, R, sign2=-1, l1=0, l2=1, n_steps=n_steps)
    _, _, _, dE_pp = run_3d_simulation(N, L, 1, 3, R, sign2=+1, l1=0, l2=1, n_steps=n_steps)
    elapsed = timer.time() - t0
    print(f"  {R:6.1f} {dE_pm:+10.4f} {dE_pp:+10.4f}  ({elapsed:.1f}s)")
