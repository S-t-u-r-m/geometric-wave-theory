"""
3D Cubic Lattice Breather Simulation
======================================
Simulates the sine-Gordon equation on a 3D cubic lattice to compute
discreteness corrections and the mu-strange splitting.

Key physics:
  - FREE breather (periodic BC): models leptons (muon) — full lattice support
  - CONFINED breather (hard walls): models quarks (strange) — limited lattice sites
  - Same (n=4, p=28) quantum numbers → same continuous mass
  - Different boundary conditions → different discrete corrections → SPLITTING

The 3D sine-Gordon equation on a cubic lattice:
  phi_tt = (1/a^2) * sum_{neighbors} (phi_j - phi_i) - sin(phi_i)
  where sum is over 2d=6 nearest neighbors in 3D

The breather is seeded as a spherically symmetric mode:
  phi(r) = 4 * arctan[ eta / (omega * cosh(eta * r)) ]
  where eta = sin(n*gamma), omega = cos(n*gamma)
"""

import numpy as np
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# GWT CONSTANTS
# ============================================================
d = 3
gamma_gwt = np.pi / (16 * np.pi - 2)
M_soliton_std = 8.0

def M_breather_continuous(n):
    return 2 * M_soliton_std * np.sin(n * gamma_gwt)

# All fermion n-values to test
n_values = [4, 5, 7, 11, 12, 13, 16, 18]
particle_names = {
    4: "muon/strange", 5: "down", 7: "bottom", 11: "charm",
    12: "top", 13: "up", 16: "electron", 18: "tau"
}

# ============================================================
# 3D LATTICE SETUP
# ============================================================
N = 48          # grid size per dimension (N^3 total sites)
a = 1.0         # Planck spacing
dt = 0.002      # time step (smaller for 3D stability: CFL condition)
N_steps = 10000
MEASURE_AFTER = 3000
MEASURE_EVERY = 50

# For confinement: proton radius in lattice units
# The proton is ~1 fm = ~1.2e20 Planck lengths in reality,
# but we're testing the GEOMETRY of confinement, not the absolute size.
# Use R_conf = N/4 to have a meaningful confined vs free comparison.
R_conf_values = [6, 8, 10, 12]  # test multiple confinement radii

print(f"3D Lattice: {N}^3 = {N**3} sites, spacing a = {a}")
print(f"Time step: {dt}, CFL limit ~ {a / np.sqrt(3):.3f}")
print(f"Steps: {N_steps}, measure after: {MEASURE_AFTER}")
print(f"Confinement radii to test: {R_conf_values}")
print()


# ============================================================
# 3D BREATHER INITIAL CONDITIONS (spherically symmetric)
# ============================================================
def make_breather_3d(N_grid, n, a_spacing, center=None):
    """
    Spherically symmetric breather on a 3D cubic lattice.
    phi(r) = 4 * arctan[ eta / (omega * cosh(eta * r)) ]
    at maximum displacement (phi_t = 0).
    """
    if center is None:
        center = np.array([N_grid // 2, N_grid // 2, N_grid // 2])

    omega = np.cos(n * gamma_gwt)
    eta = np.sin(n * gamma_gwt)

    phi = np.zeros((N_grid, N_grid, N_grid))
    ix, iy, iz = np.meshgrid(
        np.arange(N_grid), np.arange(N_grid), np.arange(N_grid), indexing='ij'
    )
    r = a_spacing * np.sqrt(
        (ix - center[0])**2 + (iy - center[1])**2 + (iz - center[2])**2
    )
    r = np.maximum(r, 1e-10)  # avoid division issues at center

    phi = 4.0 * np.arctan(eta / (omega * np.cosh(np.minimum(eta * r, 50))))
    phi_t = np.zeros_like(phi)
    return phi, phi_t


# ============================================================
# 3D ENERGY MEASUREMENT
# ============================================================
def compute_energy_3d(phi, phi_t, a_spacing):
    """
    Total energy on 3D cubic lattice.
    E = sum [ (1/2)*phi_t^2 + (1/2)*sum_nn((phi_j-phi_i)/a)^2/(2d) + (1-cos(phi)) ] * a^3
    Note: gradient term uses (1/2d) normalization to avoid double-counting
    """
    KE = 0.5 * phi_t**2

    # Gradient energy: (1/2) * sum over 3 axes of (dphi/dx_i)^2
    GE = np.zeros_like(phi)
    for axis in range(3):
        dphi = np.diff(phi, axis=axis) / a_spacing
        # Add squared gradient (forward differences)
        slices_lo = [slice(None)] * 3
        slices_lo[axis] = slice(None, -1)
        GE[tuple(slices_lo)] += 0.5 * dphi**2

    PE = 1.0 - np.cos(phi)
    return np.sum((KE + GE + PE)) * a_spacing**3


# ============================================================
# 3D ACCELERATION (discrete Laplacian)
# ============================================================
def acceleration_3d_free(phi, a_spacing):
    """3D discrete Laplacian with periodic boundary conditions (FREE particle)."""
    acc = -6 * phi / a_spacing**2 - np.sin(phi)
    # Add neighbor contributions with periodic BC
    for axis in range(3):
        acc += (np.roll(phi, 1, axis=axis) + np.roll(phi, -1, axis=axis)) / a_spacing**2
    return acc


def acceleration_3d_confined(phi, a_spacing, mask):
    """3D discrete Laplacian with hard-wall confinement (CONFINED particle).
    mask = 1 inside confinement region, 0 outside.
    phi is forced to 0 outside the mask.
    """
    phi_masked = phi * mask
    acc = -6 * phi_masked / a_spacing**2 - np.sin(phi_masked)
    for axis in range(3):
        neighbor = np.roll(phi_masked, 1, axis=axis) + np.roll(phi_masked, -1, axis=axis)
        acc += neighbor / a_spacing**2
    acc *= mask  # zero acceleration outside
    return acc


def make_confinement_mask(N_grid, R_conf):
    """Spherical confinement mask: 1 inside radius R_conf, 0 outside."""
    center = N_grid // 2
    ix, iy, iz = np.meshgrid(
        np.arange(N_grid), np.arange(N_grid), np.arange(N_grid), indexing='ij'
    )
    r = np.sqrt((ix - center)**2 + (iy - center)**2 + (iz - center)**2)
    return (r <= R_conf).astype(float)


# ============================================================
# 3D EVOLUTION (Stormer-Verlet)
# ============================================================
def evolve_3d(phi, phi_t, dt_step, a_spacing, n_steps, acc_func,
              measure_start=0, measure_every=50):
    """Stormer-Verlet for 3D field. Returns energy time series."""
    energies = []
    acc = acc_func(phi, a_spacing)

    for step in range(n_steps):
        phi_t += 0.5 * dt_step * acc
        phi += dt_step * phi_t
        acc = acc_func(phi, a_spacing)
        phi_t += 0.5 * dt_step * acc

        if step >= measure_start and step % measure_every == 0:
            E = compute_energy_3d(phi, phi_t, a_spacing)
            energies.append(E)

    return phi, phi_t, np.array(energies)


# ============================================================
# RUN: FREE vs CONFINED for each n-value
# ============================================================
print("=" * 90)
print("3D SIMULATION: FREE (periodic BC) vs CONFINED (hard wall)")
print("=" * 90)

# Use a single confinement radius for the main comparison
R_conf_main = 8
mask_main = make_confinement_mask(N, R_conf_main)
n_sites_confined = int(np.sum(mask_main))
print(f"\nConfinement radius: {R_conf_main} sites ({n_sites_confined} sites inside sphere)")
print(f"Total lattice: {N}^3 = {N**3} sites\n")

print(f"{'n':>3s}  {'Particle':>14s}  {'M_cont':>10s}  {'E_free':>12s}  {'E_conf':>12s}  "
      f"{'Corr_free':>10s}  {'Corr_conf':>10s}  {'Splitting':>10s}")
print("=" * 90)

results_3d = []

for n_val in n_values:
    name = particle_names[n_val]
    M_cont = M_breather_continuous(n_val)

    t0 = time.time()

    # --- FREE breather (periodic BC) ---
    phi_free, phit_free = make_breather_3d(N, n_val, a)
    E_init_free = compute_energy_3d(phi_free, phit_free, a)

    _, _, energies_free = evolve_3d(
        phi_free, phit_free, dt, a,
        N_steps, acceleration_3d_free,
        MEASURE_AFTER, MEASURE_EVERY
    )
    E_free = np.mean(energies_free) if len(energies_free) > 5 else E_init_free
    E_free_std = np.std(energies_free) if len(energies_free) > 5 else 0

    # --- CONFINED breather (hard wall) ---
    phi_conf, phit_conf = make_breather_3d(N, n_val, a)
    phi_conf *= mask_main  # enforce confinement from start
    E_init_conf = compute_energy_3d(phi_conf, phit_conf, a)

    def acc_conf(phi, a_sp):
        return acceleration_3d_confined(phi, a_sp, mask_main)

    _, _, energies_conf = evolve_3d(
        phi_conf, phit_conf, dt, a,
        N_steps, acc_conf,
        MEASURE_AFTER, MEASURE_EVERY
    )
    E_conf = np.mean(energies_conf) if len(energies_conf) > 5 else E_init_conf
    E_conf_std = np.std(energies_conf) if len(energies_conf) > 5 else 0

    elapsed = time.time() - t0

    # Corrections relative to 3D continuous
    # Note: M_cont is the 1D continuous mass. In 3D, the spherical breather
    # has a different normalization. We compare FREE vs CONFINED instead.
    corr_free = (E_free - E_init_free) / E_init_free * 100
    corr_conf = (E_conf - E_init_conf) / E_init_conf * 100

    # The splitting: relative difference between free and confined final energies
    splitting = (E_free - E_conf) / ((E_free + E_conf) / 2) * 100

    print(f"{n_val:3d}  {name:>14s}  {M_cont:10.4f}  {E_free:12.4f}  {E_conf:12.4f}  "
          f"{corr_free:+9.3f}%  {corr_conf:+9.3f}%  {splitting:+9.3f}%  ({elapsed:.1f}s)")

    results_3d.append({
        'n': n_val,
        'name': name,
        'M_continuous': M_cont,
        'E_free': E_free,
        'E_conf': E_conf,
        'E_free_std': E_free_std,
        'E_conf_std': E_conf_std,
        'corr_free': corr_free,
        'corr_conf': corr_conf,
        'splitting': splitting,
        'energies_free': energies_free,
        'energies_conf': energies_conf,
        'E_init_free': E_init_free,
        'E_init_conf': E_init_conf,
    })


# ============================================================
# MU-STRANGE SPLITTING ANALYSIS
# ============================================================
print("\n" + "=" * 90)
print("MU-STRANGE SPLITTING FROM 3D SIMULATION")
print("=" * 90)

for r in results_3d:
    if r['n'] == 4:
        r4 = r
        break

print(f"\nn=4 breather (both muon and strange):")
print(f"  Free (muon) energy:     {r4['E_free']:.6f} (std: {r4['E_free_std']:.6f})")
print(f"  Confined (strange) energy: {r4['E_conf']:.6f} (std: {r4['E_conf_std']:.6f})")
print(f"  Splitting (free-conf):  {r4['splitting']:+.3f}%")
print(f"\nObserved splitting:")
print(f"  muon = 105.66 MeV, strange = 93.4 MeV")
print(f"  Observed splitting: {(105.66 - 93.4) / ((105.66 + 93.4)/2) * 100:+.1f}%")
print(f"\nContinuous prediction: m(4,28) = 98.56 MeV (degenerate)")
print(f"  Observed average: {(105.66+93.4)/2:.1f} MeV (1% off from continuous)")


# ============================================================
# CONFINEMENT RADIUS SCAN (for n=4 only)
# ============================================================
print("\n" + "=" * 90)
print("CONFINEMENT RADIUS SCAN: n=4 (mu/strange)")
print("=" * 90)

print(f"{'R_conf':>8s}  {'N_sites':>8s}  {'E_free':>12s}  {'E_conf':>12s}  {'Splitting':>10s}")
print("-" * 55)

for R_c in R_conf_values:
    mask_rc = make_confinement_mask(N, R_c)
    n_in = int(np.sum(mask_rc))

    phi_f, phit_f = make_breather_3d(N, 4, a)
    E_f_init = compute_energy_3d(phi_f, phit_f, a)
    _, _, en_f = evolve_3d(phi_f, phit_f, dt, a, 8000, acceleration_3d_free, 3000, 50)
    E_f = np.mean(en_f) if len(en_f) > 5 else E_f_init

    phi_c, phit_c = make_breather_3d(N, 4, a)
    phi_c *= mask_rc

    def acc_rc(phi, a_sp, m=mask_rc):
        return acceleration_3d_confined(phi, a_sp, m)

    E_c_init = compute_energy_3d(phi_c, phit_c, a)
    _, _, en_c = evolve_3d(phi_c, phit_c, dt, a, 8000, acc_rc, 3000, 50)
    E_c = np.mean(en_c) if len(en_c) > 5 else E_c_init

    split = (E_f - E_c) / ((E_f + E_c) / 2) * 100
    print(f"{R_c:8d}  {n_in:8d}  {E_f:12.4f}  {E_c:12.4f}  {split:+9.3f}%")


# ============================================================
# PHYSICAL MASS PREDICTIONS WITH 3D CORRECTIONS
# ============================================================
print("\n" + "=" * 90)
print("PHYSICAL MASS PREDICTIONS: 3D FREE vs CONFINED CORRECTIONS")
print("=" * 90)

m_Planck_MeV = 1.2209e22
gamma_sg = np.pi / (16 * np.pi - 2)

# Fermion data: (n, p, name, obs_MeV, type)
fermion_data = [
    (4,  28, "muon",     105.66,  "free"),
    (4,  28, "strange",  93.4,    "confined"),
    (5,  30, "down",     4.67,    "confined"),
    (7,  26, "bottom",   4183,    "confined"),
    (11, 27, "charm",    1271,    "confined"),
    (12, 24, "top",      172760,  "confined"),
    (13, 31, "up",       2.16,    "confined"),
    (16, 32, "electron", 0.511,   "free"),
    (18, 27, "tau",      1776.86, "free"),
]

print(f"\n{'Particle':>10s}  {'n':>3s} {'p':>3s}  {'Type':>8s}  {'m_cont':>12s}  "
      f"{'m_3D':>12s}  {'Observed':>12s}  {'Err_cont':>8s}  {'Err_3D':>8s}")
print("-" * 100)

for n_val, p_val, name, obs_MeV, ptype in fermion_data:
    m_cont = (16.0/np.pi**2) * np.sin(n_val * gamma_sg) * np.exp(-16*p_val/np.pi**2) * m_Planck_MeV

    # Find 3D correction for this n
    for r in results_3d:
        if r['n'] == n_val:
            if ptype == "free":
                corr_3d = r['corr_free']
            else:
                corr_3d = r['corr_conf']
            break

    m_3d = m_cont * (1 + corr_3d / 100)
    err_cont = (m_cont - obs_MeV) / obs_MeV * 100
    err_3d = (m_3d - obs_MeV) / obs_MeV * 100

    print(f"{name:>10s}  {n_val:3d} {p_val:3d}  {ptype:>8s}  {m_cont:12.4f}  "
          f"{m_3d:12.4f}  {obs_MeV:12.4f}  {err_cont:+7.2f}%  {err_3d:+7.2f}%")


# ============================================================
# PLOTS
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: Free vs Confined energy for each n
ax1 = axes[0, 0]
ns = [r['n'] for r in results_3d]
e_free = [r['E_free'] for r in results_3d]
e_conf = [r['E_conf'] for r in results_3d]
x_pos = np.arange(len(ns))
ax1.bar(x_pos - 0.15, e_free, 0.3, label='Free (lepton)', alpha=0.8, color='blue')
ax1.bar(x_pos + 0.15, e_conf, 0.3, label='Confined (quark)', alpha=0.8, color='red')
ax1.set_xticks(x_pos)
ax1.set_xticklabels([f"n={n}" for n in ns], rotation=45)
ax1.set_ylabel('Energy (SG units)')
ax1.set_title('3D Breather Energy: Free vs Confined')
ax1.legend()

# Panel 2: Splitting percentage
ax2 = axes[0, 1]
splits = [r['splitting'] for r in results_3d]
ax2.bar(ns, splits, color='purple', alpha=0.7)
ax2.axhline(y=0, color='black', linestyle='--', alpha=0.5)
# Mark the observed mu-strange splitting
ax2.axhline(y=12.3, color='green', linestyle=':', alpha=0.5, label='Observed mu-s (12.3%)')
ax2.set_xlabel('Breather index n')
ax2.set_ylabel('Splitting free-confined (%)')
ax2.set_title('Free vs Confined Splitting')
ax2.legend()

# Panel 3: Energy stability for n=4
ax3 = axes[1, 0]
for r in results_3d:
    if r['n'] == 4:
        ax3.plot(r['energies_free'], label='n=4 free (muon)', color='blue')
        ax3.plot(r['energies_conf'], label='n=4 confined (strange)', color='red')
ax3.set_xlabel('Measurement index')
ax3.set_ylabel('Energy')
ax3.set_title('n=4 Energy Stability: Free vs Confined')
ax3.legend()

# Panel 4: Correction % for free and confined
ax4 = axes[1, 1]
corr_f = [r['corr_free'] for r in results_3d]
corr_c = [r['corr_conf'] for r in results_3d]
ax4.plot(ns, corr_f, 'bo-', label='Free (lepton)', markersize=8)
ax4.plot(ns, corr_c, 'rs-', label='Confined (quark)', markersize=8)
ax4.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax4.set_xlabel('Breather index n')
ax4.set_ylabel('Energy drift from initial (%)')
ax4.set_title('3D Discreteness Corrections')
ax4.legend()

plt.tight_layout()
plt.savefig('lagrangian/breather_3d_results.png', dpi=150)
print(f"\nPlot saved to lagrangian/breather_3d_results.png")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 90)
print("3D SIMULATION SUMMARY")
print("=" * 90)
for r in results_3d:
    stable_f = "STABLE" if r['E_free_std'] / max(r['E_free'], 1e-10) < 0.01 else "UNSTABLE"
    stable_c = "STABLE" if r['E_conf_std'] / max(r['E_conf'], 1e-10) < 0.01 else "UNSTABLE"
    print(f"  n={r['n']:2d} ({r['name']:>14s}): free={stable_f}, conf={stable_c}, "
          f"splitting={r['splitting']:+.3f}%")
