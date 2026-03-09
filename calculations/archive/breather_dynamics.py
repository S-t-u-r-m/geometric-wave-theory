"""
Dynamic Two-Breather Simulation
================================
Breathers are TIME-DEPENDENT — they oscillate. The static minimum
is phi=0 everywhere (no breather). Must simulate the dynamics.

Equations of motion for sine-Gordon lattice:
  d^2 phi_i / dt^2 = (phi_{i+1} + phi_{i-1} - 2*phi_i)/dx^2
                     - V_0 * pi * sin(pi * phi_i)

Use symplectic (leapfrog) integrator for energy conservation.

Breather initial conditions at t=0:
  phi(x) = (4/pi) * arctan[ eta * sin(w*0) / (w * cosh(eta*x)) ] = 0
  dphi/dt(x) = (4/pi) * eta * cos(w*0) / (w * cosh(eta*x)) * w
             = (4/pi) * eta / cosh(eta*x)

Actually, the exact sine-Gordon breather at t=0 (choosing phase=pi/2):
  phi(x,0) = (4/pi) * arctan[ (eta/w) / cosh(eta*x) ]
  dphi/dt(x,0) = 0

This is the MAX DISPLACEMENT state. Let it evolve and measure E.
"""

import numpy as np

pi = np.pi
V_0 = 1.0 / pi**2
N_br = 24
gamma_br = pi / (N_br + 1)

# Lattice setup
N = 2000
L = 80.0
dx = L / N
x = np.linspace(-L/2, L/2, N)


def breather_ic(x_arr, n, center=0.0, amp=1.0):
    """Breather at maximum displacement (velocity = 0)."""
    w = np.cos(n * gamma_br)
    eta = np.sin(n * gamma_br)
    phi = amp * (4.0/pi) * np.arctan(eta / (w * np.cosh(eta * (x_arr - center))))
    dphi = np.zeros_like(x_arr)
    return phi, dphi


def lattice_force(phi, dx_val):
    """Compute d^2phi/dt^2 from the lattice equations of motion."""
    N_pts = len(phi)
    force = np.zeros(N_pts)

    # Nearest-neighbor coupling (Laplacian)
    force[1:-1] = (phi[:-2] + phi[2:] - 2*phi[1:-1]) / dx_val**2
    force[0] = (phi[1] - phi[0]) / dx_val**2  # boundary
    force[-1] = (phi[-2] - phi[-1]) / dx_val**2

    # Sine-Gordon potential force
    force -= V_0 * pi * np.sin(pi * phi)

    return force


def total_energy(phi, dphi, dx_val):
    """Total energy = kinetic + gradient + potential."""
    E_kin = 0.5 * np.sum(dphi**2) * dx_val
    grad = np.diff(phi) / dx_val
    E_grad = 0.5 * np.sum(grad**2) * dx_val
    E_pot = V_0 * np.sum(1 - np.cos(pi * phi)) * dx_val
    return E_kin + E_grad + E_pot


def evolve(phi, dphi, dt, n_steps, dx_val):
    """Leapfrog (Stormer-Verlet) integrator."""
    # Half step velocity
    force = lattice_force(phi, dx_val)
    dphi = dphi + 0.5 * dt * force

    energies = []
    for step in range(n_steps):
        phi = phi + dt * dphi
        force = lattice_force(phi, dx_val)
        dphi = dphi + dt * force
        if step % 100 == 0:
            # Correct velocity for energy measurement
            dphi_half = dphi - 0.5 * dt * force
            E = total_energy(phi, dphi_half, dx_val)
            energies.append(E)

    # Final half step correction
    dphi = dphi - 0.5 * dt * lattice_force(phi, dx_val)
    return phi, dphi, energies


# =============================================================================
# MAIN
# =============================================================================
print("=" * 80)
print("  DYNAMIC BREATHER SIMULATION")
print(f"  Lattice: N={N}, L={L}, dx={dx:.4f}")
print(f"  V_0 = {V_0:.6f}")
print("=" * 80)

# Time step (CFL condition: dt < dx for wave equation)
dt = 0.3 * dx
T_evolve = 200.0  # total evolution time
n_steps = int(T_evolve / dt)
print(f"  dt = {dt:.5f}, T = {T_evolve}, steps = {n_steps}")

# --- Single breather energies ---
print(f"\n--- Single breather energies (time-averaged) ---")
single_E = {}

for n in [1, 2, 3, 4, 5]:
    phi, dphi = breather_ic(x, n, center=0.0)
    E_init = total_energy(phi, dphi, dx)

    phi_f, dphi_f, E_hist = evolve(phi.copy(), dphi.copy(), dt, n_steps, dx)
    E_avg = np.mean(E_hist[len(E_hist)//4:])  # skip transient
    E_std = np.std(E_hist[len(E_hist)//4:])
    E_theory = 2 * (8/pi**2) * np.sin(n * gamma_br)

    print(f"  mode {n}: E_init = {E_init:.6f}, E_avg = {E_avg:.6f} +/- {E_std:.6f}, "
          f"E_theory = {E_theory:.6f}")
    single_E[n] = E_avg


# --- Two-breather interaction ---
print(f"\n{'='*80}")
print(f"  TWO-BREATHER INTERACTION")
print(f"{'='*80}")

mode_pairs = [
    (1, 1, "H-H"),
    (1, 3, "H-2p"),
    (3, 3, "2p-2p"),
    (1, 2, "H-Li"),
]

R_values = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]

for n1, n2, label in mode_pairs:
    E1 = single_E[n1]
    E2 = single_E[n2]
    print(f"\n  --- ({n1},{n2}): {label} ---")
    print(f"  E1={E1:.6f}, E2={E2:.6f}, E1+E2={E1+E2:.6f}")
    print(f"  {'R':>6} {'E(++)':>10} {'dE(++)':>10} {'E(+-)':>10} {'dE(+-)':>10}")

    for R in R_values:
        # In-phase
        phi1_pp, dp1_pp = breather_ic(x, n1, center=-R/2, amp=+1)
        phi2_pp, dp2_pp = breather_ic(x, n2, center=+R/2, amp=+1)
        phi_pp = phi1_pp + phi2_pp
        dphi_pp = dp1_pp + dp2_pp

        _, _, E_hist_pp = evolve(phi_pp.copy(), dphi_pp.copy(), dt, n_steps, dx)
        E_pp = np.mean(E_hist_pp[len(E_hist_pp)//4:])
        dE_pp = E_pp - E1 - E2

        # Out-of-phase
        phi1_pm, dp1_pm = breather_ic(x, n1, center=-R/2, amp=+1)
        phi2_pm, dp2_pm = breather_ic(x, n2, center=+R/2, amp=-1)
        phi_pm = phi1_pm + phi2_pm
        dphi_pm = dp1_pm + dp2_pm

        _, _, E_hist_pm = evolve(phi_pm.copy(), dphi_pm.copy(), dt, n_steps, dx)
        E_pm = np.mean(E_hist_pm[len(E_hist_pm)//4:])
        dE_pm = E_pm - E1 - E2

        bond = "*" if abs(dE_pp) > 0.001 or abs(dE_pm) > 0.001 else ""
        print(f"  {R:6.1f} {E_pp:10.6f} {dE_pp:+10.6f} {E_pm:10.6f} {dE_pm:+10.6f} {bond}")


# --- Check energy conservation ---
print(f"\n--- Energy conservation check ---")
phi, dphi = breather_ic(x, 3, center=0.0)
phi_f, dphi_f, E_hist = evolve(phi.copy(), dphi.copy(), dt, n_steps, dx)
print(f"  Mode 3: E_init={E_hist[0]:.8f}, E_final={E_hist[-1]:.8f}, "
      f"drift={(E_hist[-1]-E_hist[0])/E_hist[0]*100:.4f}%")
