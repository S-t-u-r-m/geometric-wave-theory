"""
Two-Breather Dynamical Interaction on 1D Sine-Gordon Lattice
=============================================================

Previous attempt (breather_bond_sim.py) FAILED because it used static energy
minimization — the optimizer just flattened breathers to vacuum (phi=0).

Fix: Breathers are DYNAMICAL solutions. They need both phi(x,0) AND dphi/dt(x,0).
Use symplectic (leapfrog) time integration. Total energy is CONSERVED.

Bond energy = E_total(two breathers at R) - 2 * E_total(one breather)

Sine-Gordon equation on discrete lattice:
  d²phi_i/dt² = (phi_{i+1} - 2*phi_i + phi_{i-1}) - V_0 * pi * sin(pi * phi_i)

where V_0 = 1/pi² (from GWT Hamiltonian).

Exact breather solution (continuum):
  phi(x,t) = (4/pi) * arctan[ (eta/w) * sin(w*t) / cosh(eta*x) ]
  dphi/dt  = (4/pi) * (eta/cosh(eta*x)) * w*cos(w*t) / [w² + eta²*sin²(w*t)/cosh²(eta*x)]

where w = cos(n*gamma), eta = sin(n*gamma), gamma = pi/(N_br+1).
"""
import numpy as np
import time as timer

pi = np.pi
d = 3

# =============================================================================
# LATTICE PARAMETERS
# =============================================================================
V_0 = 1.0 / pi**2
N_br = 24  # floor(8*pi - 1) = 24 breather modes
gamma = pi / (N_br + 1)  # = pi/25

# Lattice setup
N = 2000        # lattice sites
L = 80.0        # total length (in lattice units)
dx = L / N
dt = 0.1 * dx   # CFL condition: dt < dx for stability
x = np.linspace(-L/2, L/2, N)


def breather_phi(x_arr, n, center=0.0, t=0.0):
    """Exact sine-Gordon breather field phi(x,t)."""
    w = np.cos(n * gamma)
    eta = np.sin(n * gamma)
    xi = eta * (x_arr - center)
    cosh_xi = np.cosh(np.clip(xi, -50, 50))  # prevent overflow
    sin_wt = np.sin(w * t)
    arg = (eta / w) * sin_wt / cosh_xi
    return (4.0 / pi) * np.arctan(arg)


def breather_dphi(x_arr, n, center=0.0, t=0.0):
    """Exact sine-Gordon breather time derivative dphi/dt(x,t)."""
    w = np.cos(n * gamma)
    eta = np.sin(n * gamma)
    xi = eta * (x_arr - center)
    cosh_xi = np.cosh(np.clip(xi, -50, 50))
    sin_wt = np.sin(w * t)
    cos_wt = np.cos(w * t)

    # dphi/dt from the arctan formula:
    # phi = (4/pi) * arctan(f) where f = (eta/w) * sin(wt) / cosh(eta*x)
    # dphi/dt = (4/pi) * df/dt / (1 + f²)
    # df/dt = eta * cos(wt) / cosh(eta*x)
    f = (eta / w) * sin_wt / cosh_xi
    dfdt = eta * cos_wt / cosh_xi
    return (4.0 / pi) * dfdt / (1 + f**2)


def total_energy(phi, dphi_dt):
    """Total energy: kinetic + gradient + potential."""
    E_kin = 0.5 * np.sum(dphi_dt**2) * dx
    grad = np.diff(phi) / dx
    E_grad = 0.5 * np.sum(grad**2) * dx
    E_pot = V_0 * np.sum(1.0 - np.cos(pi * phi)) * dx
    return E_kin + E_grad + E_pot


def leapfrog_step(phi, dphi_dt, dt_step):
    """One leapfrog step for the sine-Gordon equation."""
    # Compute acceleration: d²phi/dt² = laplacian(phi) - V_0*pi*sin(pi*phi)
    acc = np.zeros_like(phi)
    acc[1:-1] = (phi[:-2] - 2*phi[1:-1] + phi[2:]) / dx**2
    # Fixed boundary conditions (phi = 0 at boundaries)
    acc[0] = (0 - 2*phi[0] + phi[1]) / dx**2
    acc[-1] = (phi[-2] - 2*phi[-1] + 0) / dx**2
    acc -= V_0 * pi * np.sin(pi * phi)

    # Leapfrog: kick-drift-kick
    dphi_dt_half = dphi_dt + 0.5 * dt_step * acc
    phi_new = phi + dt_step * dphi_dt_half

    acc_new = np.zeros_like(phi_new)
    acc_new[1:-1] = (phi_new[:-2] - 2*phi_new[1:-1] + phi_new[2:]) / dx**2
    acc_new[0] = (0 - 2*phi_new[0] + phi_new[1]) / dx**2
    acc_new[-1] = (phi_new[-2] - 2*phi_new[-1] + 0) / dx**2
    acc_new -= V_0 * pi * np.sin(pi * phi_new)

    dphi_dt_new = dphi_dt_half + 0.5 * dt_step * acc_new
    return phi_new, dphi_dt_new


def evolve_and_measure(phi0, dphi0, n_steps, measure_interval=100):
    """Evolve and measure energy conservation."""
    phi = phi0.copy()
    dphi = dphi0.copy()
    energies = []
    E0 = total_energy(phi, dphi)
    energies.append(E0)

    for step in range(n_steps):
        phi, dphi = leapfrog_step(phi, dphi, dt)
        if (step + 1) % measure_interval == 0:
            E = total_energy(phi, dphi)
            energies.append(E)

    return np.array(energies), phi, dphi


# =============================================================================
# PART 1: Single breather energy (verify conservation)
# =============================================================================
print("=" * 80)
print("  DYNAMICAL TWO-BREATHER INTERACTION")
print("  Sine-Gordon lattice, leapfrog integration")
print("=" * 80)
print(f"  N={N}, L={L}, dx={dx:.4f}, dt={dt:.6f}")
print(f"  V_0 = 1/pi^2 = {V_0:.6f}")
print(f"  N_breather_modes = {N_br}, gamma = pi/25 = {gamma:.6f}")

# Start breathers at t = pi/(4*w) to get both phi and dphi nonzero
# (At t=0, phi=0 everywhere; at t=pi/(2w), dphi=0; pi/(4w) is in between)
t_init_factor = 0.25  # fraction of half-period

print(f"\n  --- Single breather energies ---")
print(f"  {'mode':>4} {'width':>6} {'E_total':>10} {'E_drift':>10} {'M_theory':>10}")

single_E = {}
n_evolve = 5000  # evolution steps

for n_mode in [1, 2, 3, 4, 5]:
    w_n = np.cos(n_mode * gamma)
    eta_n = np.sin(n_mode * gamma)
    width = 1.0 / eta_n
    t0 = t_init_factor * pi / w_n

    phi0 = breather_phi(x, n_mode, center=0.0, t=t0)
    dphi0 = breather_dphi(x, n_mode, center=0.0, t=t0)

    E_init = total_energy(phi0, dphi0)
    energies, _, _ = evolve_and_measure(phi0, dphi0, n_evolve, measure_interval=500)

    E_mean = np.mean(energies)
    E_drift = (energies[-1] - energies[0]) / energies[0] * 100 if energies[0] != 0 else 0

    # Theoretical mass: M_n = 2*M_kink*sin(n*gamma) = (16/pi²)*sin(n*gamma)
    M_theory = (16.0 / pi**2) * np.sin(n_mode * gamma)

    single_E[n_mode] = E_mean
    print(f"  {n_mode:4d} {width:6.2f} {E_mean:10.6f} {E_drift:+9.4f}% {M_theory:10.6f}")


# =============================================================================
# PART 2: Two-breather interaction energy vs separation
# =============================================================================
print(f"\n{'='*80}")
print(f"  TWO-BREATHER INTERACTION vs SEPARATION R")
print(f"  E_int = E(two at R) - E(one,A) - E(one,B)")
print(f"{'='*80}")

mode_pairs = [
    (1, 1, "H-H (ss)"),
    (3, 3, "2p-2p (N2,O2,F2,B2)"),
    (1, 3, "1s-2p (sp bonds)"),
    (2, 2, "2s-2s (Li-Li)"),
    (2, 3, "2s-2p (Li-F)"),
]

R_values = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]

for n1, n2, label in mode_pairs:
    w1 = np.cos(n1 * gamma); w2 = np.cos(n2 * gamma)
    t0_1 = t_init_factor * pi / w1
    t0_2 = t_init_factor * pi / w2

    E1 = single_E[n1]
    E2 = single_E[n2]

    print(f"\n  --- ({n1},{n2}): {label} ---")
    print(f"  E_single: {E1:.6f} + {E2:.6f} = {E1+E2:.6f}")
    print(f"  {'R':>6} {'E_coupled':>12} {'E_int':>12} {'E_int/E1':>10}")

    for R in R_values:
        # Set up two breathers: same phase (in-phase bonding)
        phi0 = breather_phi(x, n1, center=-R/2, t=t0_1) + \
               breather_phi(x, n2, center=+R/2, t=t0_2)
        dphi0 = breather_dphi(x, n1, center=-R/2, t=t0_1) + \
                breather_dphi(x, n2, center=+R/2, t=t0_2)

        E_init = total_energy(phi0, dphi0)
        energies, _, _ = evolve_and_measure(phi0, dphi0, n_evolve, measure_interval=500)
        E_coupled = np.mean(energies)

        E_int = E_coupled - E1 - E2
        ratio = E_int / E1 if E1 > 1e-10 else 0

        bond = " BOND" if E_int < -0.001 * E1 else ""
        print(f"  {R:6.1f} {E_coupled:12.6f} {E_int:+12.6f} {ratio:+10.4f}{bond}")


# =============================================================================
# PART 3: Detailed scan for mode (1,1) — H2 analog
# =============================================================================
print(f"\n{'='*80}")
print(f"  DETAILED SCAN: Mode (1,1) — H2 analog")
print(f"{'='*80}")

n1, n2 = 1, 1
w1 = np.cos(n1 * gamma)
t0 = t_init_factor * pi / w1
E1 = single_E[n1]

R_fine = np.arange(0.5, 20.1, 0.5)
print(f"  {'R':>6} {'E_int':>12} {'|sin(R)|':>10} {'E_int/sin':>10}")

interaction_data = []
for R in R_fine:
    phi0 = breather_phi(x, n1, center=-R/2, t=t0) + \
           breather_phi(x, n2, center=+R/2, t=t0)
    dphi0 = breather_dphi(x, n1, center=-R/2, t=t0) + \
            breather_dphi(x, n2, center=+R/2, t=t0)

    E_init = total_energy(phi0, dphi0)
    energies, _, _ = evolve_and_measure(phi0, dphi0, n_evolve, measure_interval=500)
    E_coupled = np.mean(energies)
    E_int = E_coupled - 2 * E1

    sin_R = abs(np.sin(R))
    ratio_sin = E_int / sin_R if sin_R > 0.01 else float('nan')
    interaction_data.append((R, E_int, sin_R))

    print(f"  {R:6.1f} {E_int:+12.6f} {sin_R:10.4f} {ratio_sin:+10.4f}")


# =============================================================================
# PART 4: Anti-phase (antibonding) comparison
# =============================================================================
print(f"\n{'='*80}")
print(f"  BONDING vs ANTIBONDING: Mode (3,3)")
print(f"  In-phase (+,+) vs out-of-phase (+,-)")
print(f"{'='*80}")

n1, n2 = 3, 3
w1 = np.cos(n1 * gamma)
t0 = t_init_factor * pi / w1
E1 = single_E[n1]

R_test = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
print(f"  {'R':>6} {'E_int(++)':>12} {'E_int(+-)':>12} {'split':>10}")

for R in R_test:
    # In-phase (bonding)
    phi0_pp = breather_phi(x, n1, center=-R/2, t=t0) + \
              breather_phi(x, n2, center=+R/2, t=t0)
    dphi0_pp = breather_dphi(x, n1, center=-R/2, t=t0) + \
               breather_dphi(x, n2, center=+R/2, t=t0)

    E_pp_init = total_energy(phi0_pp, dphi0_pp)
    ens_pp, _, _ = evolve_and_measure(phi0_pp, dphi0_pp, n_evolve, measure_interval=500)
    E_pp = np.mean(ens_pp) - 2*E1

    # Out-of-phase (antibonding): second breather with opposite sign
    phi0_pm = breather_phi(x, n1, center=-R/2, t=t0) - \
              breather_phi(x, n2, center=+R/2, t=t0)
    dphi0_pm = breather_dphi(x, n1, center=-R/2, t=t0) - \
               breather_dphi(x, n2, center=+R/2, t=t0)

    E_pm_init = total_energy(phi0_pm, dphi0_pm)
    ens_pm, _, _ = evolve_and_measure(phi0_pm, dphi0_pm, n_evolve, measure_interval=500)
    E_pm = np.mean(ens_pm) - 2*E1

    split = E_pm - E_pp
    bond_label = "bond<anti" if E_pp < E_pm else "anti<bond" if E_pm < E_pp else "equal"
    print(f"  {R:6.1f} {E_pp:+12.6f} {E_pm:+12.6f} {split:+10.6f}  {bond_label}")


print(f"\n{'='*80}")
print(f"  DONE. Check if interaction energy follows |sin(R)| pattern.")
print(f"{'='*80}")
