"""
Two-Breather Bond Energy on Discrete Lattice
=============================================
Does bond energy EMERGE from two breathers interacting on the d=3 lattice?

Setup:
  - Truly discrete 1D lattice (a=1, periodic BC)
  - Two n=1 breathers separated by R lattice sites
  - Opposite phase (+,-) = bonding, same phase (+,+) = antibonding
  - Measure: E_pair(R) - 2*E_single = interaction energy V(R)
  - Scan R from 2 to 40 sites

The bond energy should emerge from the Lagrangian without any formula.
If V(R) shows a Morse-like well, GWT predicts D_e = pi/d^2 * E_harm.

EOM: phi_ddot_i = (phi_{i+1} + phi_{i-1} - 2*phi_i) - (1/pi)*sin(pi*phi_i)
"""
import sys, io, os, time
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "breather_bond_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("TWO-BREATHER BOND ENERGY ON DISCRETE LATTICE")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"gamma = {gamma:.10f}")
report("")

# ============================================================
# DISCRETE 1D LATTICE ENGINE
# ============================================================
N = 512  # lattice sites
dt = 0.05
center = N // 2

def evolve_and_measure(phi_init, phi_dot_init, N_periods=30, omega_ref=None):
    """Evolve a 1D discrete lattice and measure time-averaged energy.

    Returns: E_avg, E_std, phi_final, omega_measured
    """
    phi = phi_init.copy()
    phi_old = phi - dt * phi_dot_init

    if omega_ref is None:
        omega_ref = np.cos(gamma)  # n=1 breather frequency
    period = 2*PI / omega_ref
    N_steps = int(N_periods * period / dt)
    N_steps = min(N_steps, 1500000)

    # Let the system settle for 5 periods before measuring
    settle_steps = int(5 * period / dt)
    rec = max(1, (N_steps - settle_steps) // 5000)

    energies = []
    ts_center = []

    for step in range(N_steps):
        # Discrete 1D Laplacian with periodic BC
        lap = np.roll(phi, 1) + np.roll(phi, -1) - 2*phi
        force = (1.0/PI) * np.sin(PI * phi)
        phi_new = 2*phi - phi_old + dt**2 * (lap - force)
        phi_old = phi.copy()
        phi = phi_new

        if step >= settle_steps and step % rec == 0:
            # Total energy (conserved quantity)
            KE = 0.5 * np.sum(((phi - phi_old)/dt)**2)
            PE_site = np.sum((1.0/PI**2) * (1 - np.cos(PI * phi)))
            PE_coup = 0.5 * np.sum((np.roll(phi, -1) - phi)**2)
            E_total = KE + PE_site + PE_coup
            energies.append(E_total)
            ts_center.append(phi[center])

    energies = np.array(energies)
    E_avg = np.mean(energies)
    E_std = np.std(energies)

    # Measure frequency
    ts = np.array(ts_center)
    ts_mean = ts - np.mean(ts)
    crossings = []
    for i in range(len(ts_mean)-1):
        if ts_mean[i] <= 0 and ts_mean[i+1] > 0:
            t_cross = i + (-ts_mean[i]) / (ts_mean[i+1] - ts_mean[i] + 1e-30)
            crossings.append(t_cross)
    if len(crossings) >= 3:
        periods_list = np.diff(crossings) * dt * rec
        T_meas = np.median(periods_list)
        omega_meas = 2*PI / T_meas
    else:
        omega_meas = 0

    return E_avg, E_std, phi, omega_meas

# ============================================================
# SINGLE BREATHER REFERENCE ENERGY
# ============================================================
report("STEP 1: Single breather energy (n=1)")
report("-" * 55)

# Test multiple breather modes
for n_mode in [1, 2, 3, 4]:
    omega_n = np.cos(n_mode * gamma)
    eps_n = np.sin(n_mode * gamma)

    x = np.arange(N, dtype=np.float64) - center
    phi_init = np.zeros(N)
    phi_dot_init = (4.0/PI) * eps_n / (omega_n * np.cosh(eps_n * x) + 1e-30)

    E_single, E_std, _, omega_meas = evolve_and_measure(
        phi_init, phi_dot_init, N_periods=30, omega_ref=omega_n)

    report(f"  n={n_mode}: E_single = {E_single:.6f} +/- {E_std:.6f}, "
           f"omega = {omega_meas:.6f} (pred {omega_n:.6f})")

# Use n=1 as the reference for H2-like bonds
omega_1 = np.cos(gamma)
eps_1 = np.sin(gamma)

x = np.arange(N, dtype=np.float64) - center
phi_single = np.zeros(N)
v_single = (4.0/PI) * eps_1 / (omega_1 * np.cosh(eps_1 * x) + 1e-30)
E_single_ref, E_std_ref, _, _ = evolve_and_measure(
    phi_single, v_single, N_periods=40, omega_ref=omega_1)

report(f"\n  Reference: E_single(n=1) = {E_single_ref:.6f} +/- {E_std_ref:.6f}")
report("")

# ============================================================
# TWO-BREATHER SCAN: Opposite phase (bonding)
# ============================================================
report("STEP 2: Two-breather interaction V(R)")
report("-" * 55)
report("Opposite phase (+,-) = bonding")
report("Same phase (+,+) = antibonding")
report("")

R_values = list(range(3, 42, 2)) + [50, 60, 80]  # site separations

report(f"{'R':>4} {'E_bond':>12} {'E_anti':>12} {'V_bond':>12} {'V_anti':>12} "
       f"{'V_split':>12} {'E_std_b':>10}")
report("-" * 80)

results = []

for R_sep in R_values:
    # Place two breathers at center +/- R_sep/2
    pos_A = center - R_sep // 2
    pos_B = center + R_sep // 2

    x = np.arange(N, dtype=np.float64)

    # Breather A centered at pos_A
    x_A = x - pos_A
    v_A = (4.0/PI) * eps_1 / (omega_1 * np.cosh(eps_1 * x_A) + 1e-30)

    # Breather B centered at pos_B
    x_B = x - pos_B
    v_B = (4.0/PI) * eps_1 / (omega_1 * np.cosh(eps_1 * x_B) + 1e-30)

    # BONDING: opposite phase (A positive, B negative velocity)
    phi_bond = np.zeros(N)
    v_bond = v_A - v_B  # opposite phase

    E_bond, E_std_b, _, _ = evolve_and_measure(
        phi_bond, v_bond, N_periods=25, omega_ref=omega_1)

    # ANTIBONDING: same phase
    v_anti = v_A + v_B  # same phase

    E_anti, E_std_a, _, _ = evolve_and_measure(
        phi_bond, v_anti, N_periods=25, omega_ref=omega_1)

    # Interaction energies
    V_bond = E_bond - 2 * E_single_ref
    V_anti = E_anti - 2 * E_single_ref
    V_split = V_anti - V_bond  # bonding-antibonding splitting

    results.append((R_sep, E_bond, E_anti, V_bond, V_anti, V_split, E_std_b))

    report(f"{R_sep:4d} {E_bond:12.6f} {E_anti:12.6f} {V_bond:+12.6f} {V_anti:+12.6f} "
           f"{V_split:12.6f} {E_std_b:10.6f}")

report("")

# ============================================================
# ANALYSIS: Extract bond parameters
# ============================================================
report("STEP 3: Bond energy analysis")
report("-" * 55)

R_arr = np.array([r[0] for r in results])
V_bond_arr = np.array([r[3] for r in results])
V_anti_arr = np.array([r[4] for r in results])
V_split_arr = np.array([r[5] for r in results])

# Find minimum of bonding potential
i_min = np.argmin(V_bond_arr)
R_eq = R_arr[i_min]
D_e_lattice = -V_bond_arr[i_min]  # well depth in lattice units

report(f"Bonding potential minimum:")
report(f"  R_eq = {R_eq} lattice sites")
report(f"  D_e = {D_e_lattice:.6f} (lattice units)")
report(f"  V_bond(R_eq) = {V_bond_arr[i_min]:.6f}")
report("")

# Antibonding at R_eq
V_anti_at_eq = V_anti_arr[i_min]
report(f"At R_eq = {R_eq}:")
report(f"  V_bond = {V_bond_arr[i_min]:+.6f}")
report(f"  V_anti = {V_anti_at_eq:+.6f}")
report(f"  Splitting = {V_anti_at_eq - V_bond_arr[i_min]:.6f}")
report("")

# Check asymptotic behavior (should go to 0 at large R)
report(f"Asymptotic (R={R_arr[-1]}): V_bond = {V_bond_arr[-1]:+.6f}, V_anti = {V_anti_arr[-1]:+.6f}")
report("")

# The key ratio: bonding/antibonding splitting
# In GWT, antibonding should be exactly 2x stronger than bonding (f_anti = 6/5)
if D_e_lattice > 0 and V_anti_at_eq > 0:
    ratio = V_anti_at_eq / D_e_lattice
    report(f"Antibonding/bonding ratio: {ratio:.3f}")
    report(f"  GWT predicts f_anti = 2d/(2d-1) = {2*d/(2*d-1):.3f}")
    report("")

# Energy scale: single breather energy
report(f"Single breather energy: {E_single_ref:.6f} (lattice units)")
report(f"Bond depth / single energy: {D_e_lattice / E_single_ref:.6f}")
report("")

# The Pöschl-Teller parameter
# s = (-1 + sqrt(1 + 8/pi^2)) / 2 from the source of truth
s = (-1 + np.sqrt(1 + 8/PI**2)) / 2
report(f"Pöschl-Teller s = {s:.6f}")
report(f"s^2 = {s**2:.6f}")
report(f"D_e / (2*s^2) = {D_e_lattice / (2*s**2):.6f} (should relate to E_single)")
report("")

# GWT prediction for H2
alpha_bare = np.exp(-(2/6) * (2**7/PI**2 + np.log(6)))
E_H = alpha_bare**2 / 2 * 0.511e6  # eV
D_e_pred = PI / d**2 * E_H  # eV
report(f"GWT analytical prediction:")
report(f"  D_e(H2) = pi/d^2 * E_H = {D_e_pred:.3f} eV")
report(f"  R_eq(H2) = 1.401 Bohr")
report("")

# Can we find a conversion factor?
# If D_e_lattice (lattice units) maps to D_e_pred (eV):
if D_e_lattice > 0:
    eV_per_lattice = D_e_pred / D_e_lattice
    report(f"If this IS H2: conversion = {eV_per_lattice:.3f} eV per lattice unit")
    report(f"  Single breather = {E_single_ref * eV_per_lattice:.1f} eV = "
           f"{E_single_ref * eV_per_lattice / 1e6:.4f} MeV")
    report(f"  (electron mass = 0.511 MeV)")
report("")

report("=" * 70)
report("KEY QUESTION: Does the potential have a WELL?")
report("If V_bond goes negative and comes back to 0, a bond forms.")
report("The well depth IS the bond energy. No formula needed.")
report("=" * 70)
report("")

# Summary table of potential curve
report("POTENTIAL CURVE:")
report(f"{'R':>4} {'V_bond':>12} {'V_anti':>12} {'split':>12}")
for r in results:
    marker = " <-- min" if r[0] == R_eq else ""
    report(f"{r[0]:4d} {r[3]:+12.6f} {r[4]:+12.6f} {r[5]:12.6f}{marker}")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))

log.close()
print(f"\nResults saved to: {outfile}")
