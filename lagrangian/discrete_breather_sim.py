"""
Discrete Lattice Breather Simulation
======================================
Simulates the sine-Gordon equation on a finite discrete lattice,
seeds DHN breathers at each GWT n-value, measures actual energies
vs the continuous formula, and quantifies discreteness corrections.

The sine-Gordon equation (standard form):
  phi_tt - phi_xx + sin(phi) = 0

On a discrete lattice with spacing a:
  phi_tt = (phi_{i+1} - 2*phi_i + phi_{i-1})/a^2 - sin(phi_i)

The continuous breather solution at maximum displacement (phi_t = 0):
  phi(x) = 4 * arctan[ eta/omega / cosh(eta * x) ]
  where eta = sqrt(1 - omega^2), omega = cos(n*gamma)

Energy is measured as:
  E = sum_i [ (1/2)*phi_t_i^2 + (1/2)*((phi_{i+1}-phi_i)/a)^2 + (1-cos(phi_i)) ] * a

GWT parameters:
  gamma = pi / (16*pi - 2)
  M_kink = 8/pi^2 (Planck units) => 16 in standard sine-Gordon units
  M_n = 2 * M_kink_std * sin(n*gamma)  (continuous prediction)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# GWT CONSTANTS
# ============================================================
d = 3
gamma_gwt = np.pi / (16 * np.pi - 2)  # ~ 0.0648
N_breathers = 24

# In standard sine-Gordon (phi_tt - phi_xx + sin(phi) = 0):
# soliton mass M_s = 8 (with beta=1 normalization)
M_soliton_std = 8.0

# GWT maps: M_kink = 8/pi^2 Planck units, but in the standard SG
# the soliton mass is just 8. The ratio is the unit conversion.
# Breather mass in standard SG units:
def M_breather_continuous(n):
    """Continuous DHN breather mass for mode n."""
    return 2 * M_soliton_std * np.sin(n * gamma_gwt)

# GWT fermion assignments: (n, p, name, observed_MeV)
fermions = [
    (4,  28, "muon/strange", 98.56),    # degenerate in continuous
    (5,  30, "down",         4.67),
    (7,  26, "bottom",       4183),
    (11, 27, "charm",        1271),
    (12, 24, "top",          172760),
    (13, 31, "up",           2.16),
    (16, 32, "electron",     0.511),
    (18, 27, "tau",          1776.86),
]

# ============================================================
# SIMULATION PARAMETERS
# ============================================================
# The GWT lattice has spacing a = 1 Planck length (by definition).
# In standard SG units, this IS the natural length scale.
# We test BOTH a=1 (physical Planck lattice) and a=0.1 (fine reference).

a_planck = 1.0      # Physical Planck lattice spacing
L_half = 30.0       # half-domain size (in SG natural units)
N_sites_planck = int(2 * L_half / a_planck) + 1  # = 61 sites
dt = 0.005          # time step
N_steps = 20000     # evolution steps (total time = 100 SG units)
MEASURE_AFTER = 5000  # start measuring energy after transient

# Also run fine lattice for reference
a_fine = 0.05
N_sites_fine = int(2 * L_half / a_fine) + 1

# Default lattice (Planck spacing)
a = a_planck
N_sites = N_sites_planck
x = np.linspace(-L_half, L_half, N_sites)

print(f"Lattice: {N_sites} sites, spacing a = {a:.4f}, domain [-{L_half}, {L_half}]")
print(f"Time step: {dt}, total steps: {N_steps}, measure after: {MEASURE_AFTER}")
print(f"GWT gamma = {gamma_gwt:.6f}")
print()

# ============================================================
# BREATHER INITIAL CONDITIONS
# ============================================================
def breather_ic(x_arr, n, a_spacing):
    """
    Initialize breather at maximum displacement (phi_t = 0).
    phi(x) = 4 * arctan[ eta / (omega * cosh(eta * x)) ]
    omega = cos(n * gamma), eta = sqrt(1 - omega^2) = sin(n * gamma)
    """
    omega = np.cos(n * gamma_gwt)
    eta = np.sin(n * gamma_gwt)

    if eta < 1e-10 or omega < 1e-10:
        return np.zeros_like(x_arr), np.zeros_like(x_arr)

    phi = 4.0 * np.arctan(eta / (omega * np.cosh(eta * x_arr)))
    phi_t = np.zeros_like(x_arr)
    return phi, phi_t


# ============================================================
# ENERGY MEASUREMENT
# ============================================================
def compute_energy(phi, phi_t, a_spacing):
    """
    Total energy of the discrete lattice field.
    E = sum_i [ (1/2)*phi_t^2 + (1/2)*((phi_{i+1}-phi_i)/a)^2 + (1-cos(phi_i)) ] * a
    """
    KE = 0.5 * phi_t**2
    # Gradient energy (use forward differences, last site wraps or is boundary)
    dphi = np.diff(phi) / a_spacing
    GE = np.zeros_like(phi)
    GE[:-1] = 0.5 * dphi**2
    # Potential energy
    PE = 1.0 - np.cos(phi)
    return np.sum((KE + GE + PE) * a_spacing)


# ============================================================
# SYMPLECTIC (STORMER-VERLET) INTEGRATOR
# ============================================================
def acceleration(phi, a_spacing):
    """
    phi_tt = (phi_{i+1} - 2*phi_i + phi_{i-1})/a^2 - sin(phi_i)
    Fixed boundary: phi[0] = phi[-1] = 0
    """
    acc = np.zeros_like(phi)
    # Discrete Laplacian (fixed boundaries)
    acc[1:-1] = (phi[2:] - 2*phi[1:-1] + phi[:-2]) / a_spacing**2 - np.sin(phi[1:-1])
    # Boundary sites: neighbors are 0
    acc[0] = (phi[1] - 2*phi[0] + 0) / a_spacing**2 - np.sin(phi[0])
    acc[-1] = (0 - 2*phi[-1] + phi[-2]) / a_spacing**2 - np.sin(phi[-1])
    return acc


def evolve(phi, phi_t, dt_step, a_spacing, n_steps, measure_start=0):
    """Stormer-Verlet integration. Returns energy time series."""
    energies = []
    acc = acceleration(phi, a_spacing)

    for step in range(n_steps):
        # Velocity Verlet
        phi_t += 0.5 * dt_step * acc
        phi += dt_step * phi_t
        acc = acceleration(phi, a_spacing)
        phi_t += 0.5 * dt_step * acc

        if step >= measure_start and step % 50 == 0:
            E = compute_energy(phi, phi_t, a_spacing)
            energies.append(E)

    return phi, phi_t, np.array(energies)


# ============================================================
# RUN SIMULATION FOR EACH n-VALUE
# ============================================================
print("=" * 85)
print(f"{'n':>3s}  {'Particle':>14s}  {'M_cont':>10s}  {'E_discrete':>12s}  "
      f"{'E_std':>8s}  {'Correction':>10s}  {'E_stable':>8s}")
print("=" * 85)

results = []

for n_val in sorted(set(n for n, _, _, _ in fermions)):
    # Find particle name
    names = [name for n, p, name, obs in fermions if n == n_val]
    name = names[0] if names else f"n={n_val}"

    # Continuous prediction (standard SG units)
    M_cont = M_breather_continuous(n_val)

    # Initialize breather
    phi, phi_t_field = breather_ic(x, n_val, a)

    # Initial energy (should match continuous prediction)
    E_init = compute_energy(phi, phi_t_field, a)

    # Evolve
    phi_final, phi_t_final, energies = evolve(
        phi.copy(), phi_t_field.copy(), dt, a, N_steps, MEASURE_AFTER
    )

    # Measure stable energy
    if len(energies) > 10:
        E_mean = np.mean(energies)
        E_std = np.std(energies)
    else:
        E_mean = E_init
        E_std = 0

    # Discreteness correction
    correction_pct = (E_mean - M_cont) / M_cont * 100

    print(f"{n_val:3d}  {name:>14s}  {M_cont:10.4f}  {E_mean:12.6f}  "
          f"{E_std:8.4f}  {correction_pct:+9.3f}%  {'STABLE' if E_std/E_mean < 0.01 else 'UNSTABLE'}")

    results.append({
        'n': n_val,
        'name': name,
        'M_continuous': M_cont,
        'E_discrete': E_mean,
        'E_std': E_std,
        'correction_pct': correction_pct,
        'energies': energies,
        'E_init': E_init,
    })


# ============================================================
# ANALYSIS: MAP BACK TO PHYSICAL MASSES
# ============================================================
print("\n" + "=" * 85)
print("PHYSICAL MASS PREDICTIONS WITH DISCRETENESS CORRECTIONS")
print("=" * 85)

m_Planck_MeV = 1.2209e22
gamma_sg = np.pi / (16 * np.pi - 2)

print(f"\n{'Particle':>14s}  {'n':>3s} {'p':>3s}  {'m_cont':>12s}  {'m_disc':>12s}  "
      f"{'Observed':>12s}  {'Err_cont':>9s}  {'Err_disc':>9s}")
print("-" * 95)

for n_val, p_val, name, obs_MeV in fermions:
    # Continuous mass
    m_cont = (16.0/np.pi**2) * np.sin(n_val * gamma_sg) * np.exp(-16*p_val/np.pi**2) * m_Planck_MeV

    # Find the correction factor for this n
    for r in results:
        if r['n'] == n_val:
            corr_factor = 1.0 + r['correction_pct'] / 100.0
            break

    # Discrete-corrected mass
    m_disc = m_cont * corr_factor

    err_cont = (m_cont - obs_MeV) / obs_MeV * 100
    err_disc = (m_disc - obs_MeV) / obs_MeV * 100

    print(f"{name:>14s}  {n_val:3d} {p_val:3d}  {m_cont:12.4f}  {m_disc:12.4f}  "
          f"{obs_MeV:12.4f}  {err_cont:+8.2f}%  {err_disc:+8.2f}%")


# ============================================================
# MU-STRANGE SPLITTING ANALYSIS
# ============================================================
print("\n" + "=" * 85)
print("MU-STRANGE DEGENERACY ANALYSIS")
print("=" * 85)

# Both have n=4, p=28 in continuous theory
m_cont_4_28 = (16.0/np.pi**2) * np.sin(4 * gamma_sg) * np.exp(-16*28/np.pi**2) * m_Planck_MeV

for r in results:
    if r['n'] == 4:
        corr = r['correction_pct']
        break

print(f"Continuous m(4,28) = {m_cont_4_28:.2f} MeV")
print(f"Observed muon     = 105.66 MeV")
print(f"Observed strange   = 93.4 MeV")
print(f"Observed average   = {(105.66 + 93.4)/2:.2f} MeV")
print(f"Discrete lattice correction for n=4: {corr:+.3f}%")
print(f"\nNote: The mu-strange splitting requires DIFFERENT corrections for")
print(f"free (muon) vs confined (strange) particles — same n but different")
print(f"effective lattice support. A 1D simulation gives the SAME correction")
print(f"to both. Full 3D confinement geometry needed to split them.")


# ============================================================
# FINE LATTICE COMPARISON
# ============================================================
print("\n" + "=" * 85)
print("FINE LATTICE (a=0.05) vs PLANCK LATTICE (a=1.0)")
print("=" * 85)

x_fine = np.linspace(-L_half, L_half, N_sites_fine)
print(f"{'n':>3s}  {'M_cont':>10s}  {'E_planck':>12s}  {'E_fine':>12s}  "
      f"{'Corr_planck':>12s}  {'Corr_fine':>10s}")
print("-" * 70)

for r in results:
    n_val = r['n']
    phi_f, phit_f = breather_ic(x_fine, n_val, a_fine)
    _, _, en_f = evolve(phi_f, phit_f, dt, a_fine, N_steps, MEASURE_AFTER)
    E_fine = np.mean(en_f) if len(en_f) > 5 else compute_energy(phi_f, phit_f, a_fine)
    corr_fine = (E_fine - r['M_continuous']) / r['M_continuous'] * 100
    print(f"{n_val:3d}  {r['M_continuous']:10.4f}  {r['E_discrete']:12.6f}  {E_fine:12.6f}  "
          f"{r['correction_pct']:+11.4f}%  {corr_fine:+9.4f}%")
    r['E_fine'] = E_fine
    r['corr_fine'] = corr_fine
    # The PHYSICAL correction is Planck minus fine (isolate discreteness)
    r['physical_correction'] = r['correction_pct'] - corr_fine

print(f"\n{'n':>3s}  {'Particle':>14s}  {'PHYSICAL discreteness':>22s}")
print("-" * 45)
for r in results:
    print(f"{r['n']:3d}  {r['name']:>14s}  {r['physical_correction']:+.4f}%")


# ============================================================
# CONVERGENCE TEST: vary lattice spacing
# ============================================================
print("\n" + "=" * 85)
print("CONVERGENCE TEST: n=12 (top) vs lattice spacing")
print("=" * 85)

test_n = 12
M_cont_12 = M_breather_continuous(test_n)
spacings = [1.0, 0.5, 0.25, 0.1, 0.05]

print(f"{'a (spacing)':>12s}  {'N_sites':>8s}  {'E_discrete':>12s}  {'Correction':>10s}")
print("-" * 50)

convergence_data = []
for a_test in spacings:
    N_test = int(2 * L_half / a_test) + 1
    if N_test % 2 == 0:
        N_test += 1
    x_test = np.linspace(-L_half, L_half, N_test)

    phi_test, phit_test = breather_ic(x_test, test_n, a_test)
    _, _, energies_test = evolve(phi_test, phit_test, dt, a_test, 5000, 2000)

    if len(energies_test) > 5:
        E_test = np.mean(energies_test)
    else:
        E_test = compute_energy(phi_test, phit_test, a_test)

    corr_test = (E_test - M_cont_12) / M_cont_12 * 100
    convergence_data.append((a_test, N_test, E_test, corr_test))
    print(f"{a_test:12.3f}  {N_test:8d}  {E_test:12.6f}  {corr_test:+9.4f}%")


# ============================================================
# PLOT RESULTS
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: Continuous vs Discrete energy for each n
ax1 = axes[0, 0]
n_vals = [r['n'] for r in results]
M_conts = [r['M_continuous'] for r in results]
E_discs = [r['E_discrete'] for r in results]
ax1.bar(np.arange(len(n_vals)) - 0.15, M_conts, 0.3, label='Continuous', alpha=0.8)
ax1.bar(np.arange(len(n_vals)) + 0.15, E_discs, 0.3, label='Discrete', alpha=0.8)
ax1.set_xticks(np.arange(len(n_vals)))
ax1.set_xticklabels([f"n={n}" for n in n_vals], rotation=45)
ax1.set_ylabel('Energy (SG units)')
ax1.set_title('Breather Energy: Continuous vs Discrete Lattice')
ax1.legend()

# Panel 2: Correction percentage vs n
ax2 = axes[0, 1]
corrections = [r['correction_pct'] for r in results]
colors = ['red' if abs(c) > 1 else 'green' for c in corrections]
ax2.bar(n_vals, corrections, color=colors, alpha=0.7)
ax2.axhline(y=0, color='black', linestyle='--', alpha=0.5)
ax2.set_xlabel('Breather index n')
ax2.set_ylabel('Discreteness correction (%)')
ax2.set_title('Lattice Discreteness Corrections')

# Panel 3: Energy stability over time for n=4 (mu/strange)
ax3 = axes[1, 0]
for r in results:
    if r['n'] == 4:
        ax3.plot(r['energies'], label=f"n=4 (mu/strange)")
    if r['n'] == 12:
        ax3.plot(r['energies'], label=f"n=12 (top)")
ax3.set_xlabel('Measurement index')
ax3.set_ylabel('Energy')
ax3.set_title('Energy Stability During Evolution')
ax3.legend()

# Panel 4: Convergence with lattice spacing
ax4 = axes[1, 1]
a_vals = [c[0] for c in convergence_data]
corr_vals = [c[3] for c in convergence_data]
ax4.plot(a_vals, corr_vals, 'bo-', markersize=8)
ax4.set_xlabel('Lattice spacing a')
ax4.set_ylabel('Correction (%)')
ax4.set_title(f'Convergence: n=12 correction vs spacing')
ax4.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('lagrangian/breather_simulation_results.png', dpi=150)
print(f"\nPlot saved to lagrangian/breather_simulation_results.png")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 85)
print("SUMMARY")
print("=" * 85)
print(f"Lattice spacing used: a = {a:.4f} (in SG natural units)")
print(f"Planck spacing equivalent: a_Planck = 1 (by definition)")
print(f"\nKey findings:")
for r in results:
    print(f"  n={r['n']:2d} ({r['name']:>14s}): correction = {r['correction_pct']:+.3f}%, "
          f"energy stable = {'YES' if r['E_std']/r['E_discrete'] < 0.01 else 'NO'}")

avg_corr = np.mean([abs(r['correction_pct']) for r in results])
print(f"\nAverage |correction| = {avg_corr:.3f}%")
print(f"This represents the leading-order effect of lattice discreteness")
print(f"on breather energies in 1D. The full 3D correction and confinement")
print(f"effects require a 3D simulation.")
