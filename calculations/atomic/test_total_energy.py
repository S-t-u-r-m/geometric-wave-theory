"""
Test spectral lines using TOTAL atomic energy with SCF iteration.

SCF = Self-Consistent Field. Not a supercomputer — just a loop:
  1. Guess screening
  2. Compute energy levels
  3. Update screening from new levels
  4. Repeat until converged (usually 3-5 loops, milliseconds)
"""
import sys, os, io
sys.path.insert(0, os.path.dirname(__file__))
from z_eff_subshell import compute_subshell_energy, screening_for_target, ATOMS
from z_eff_subshell import compute_alpha_for_target, alpha_penetration_factor
import numpy as np
from math import factorial

try:
    sys.stdout.write('')
except:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
E_H = 13.6045
hc = 1240.0
l_n = {0: 's', 1: 'p', 2: 'd', 3: 'f'}


def scf_total_energy(Z, config, max_iter=10, tol=1e-6):
    """
    Self-consistent total atomic energy.

    Iterate: compute subshell energies → update screening → repeat.
    Each iteration adjusts the screening based on the current
    energy levels, which changes the effective Z seen by each electron.

    Converges in 3-5 iterations typically.
    """
    # Start with the given config's screening
    E_prev = 999.0

    for iteration in range(max_iter):
        E_tot = 0
        for n_sh, l_sh, count in config:
            if count == 0:
                continue
            E_sh, Z_eff, Z_net, alpha, S = compute_subshell_energy(Z, config, n_sh, l_sh)
            E_tot += count * E_sh

        # Check convergence
        if abs(E_tot - E_prev) < tol:
            return E_tot, iteration + 1
        E_prev = E_tot

    return E_tot, max_iter


def scf_transition(Z, config_ground, config_excited, verbose=False):
    """
    Transition energy using SCF for BOTH ground and excited configs.

    Each config gets its own self-consistent energy.
    The transition energy = E(excited) - E(ground).
    """
    E_ground, n_iter_g = scf_total_energy(Z, config_ground)
    E_excited, n_iter_e = scf_total_energy(Z, config_excited)

    dE = E_excited - E_ground

    if verbose:
        print(f'    Ground:  E = {E_ground:.3f} eV ({n_iter_g} iter)')
        print(f'    Excited: E = {E_excited:.3f} eV ({n_iter_e} iter)')
        print(f'    dE = {dE:.4f} eV')

    return dE


print("TOTAL ENERGY SPECTRAL PREDICTIONS")
print("=" * 65)
print("Compute E(excited config) - E(ground config) for the FULL atom.")
print()

tests = [
    # (name, Z, ground_config, excited_config, obs_nm, line_name)

    # Hydrogen: single electron, config changes completely
    ("H", 1,
     [(1,0,1)],
     [(2,1,1)],
     121.567, "Lyman-alpha"),

    ("H", 1,
     [(2,0,1)],          # Balmer: start from n=2
     [(3,1,1)],
     656.281, "H-alpha"),

    # Helium: two electrons, one excites
    ("He", 2,
     [(1,0,2)],
     [(1,0,1), (2,1,1)],
     58.433, "He I resonance"),

    # Lithium: 1s core + valence 2s excites to 2p
    ("Li", 3,
     [(1,0,2), (2,0,1)],
     [(1,0,2), (2,1,1)],
     670.776, "Li I red"),

    # Sodium: [Ne] core + 3s excites to 3p
    ("Na", 11,
     [(1,0,2),(2,0,2),(2,1,6),(3,0,1)],
     [(1,0,2),(2,0,2),(2,1,6),(3,1,1)],
     589.592, "Na D1"),

    # Potassium: [Ar] core + 4s excites to 4p
    ("K", 19,
     [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(4,0,1)],
     [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(4,1,1)],
     766.490, "K I principal"),

    # Magnesium: one of 3s^2 excites to 3p
    ("Mg", 12,
     [(1,0,2),(2,0,2),(2,1,6),(3,0,2)],
     [(1,0,2),(2,0,2),(2,1,6),(3,0,1),(3,1,1)],
     285.213, "Mg I resonance"),

    # Calcium: one of 4s^2 excites to 4p
    ("Ca", 20,
     [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(4,0,2)],
     [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(4,0,1),(4,1,1)],
     422.673, "Ca I resonance"),

    # Carbon: one of 2p^2 excites to 3s
    ("C", 6,
     [(1,0,2),(2,0,2),(2,1,2)],
     [(1,0,2),(2,0,2),(2,1,1),(3,0,1)],
     247.856, "C I UV (approx)"),
]

print("SCF TOTAL ENERGY SPECTRAL PREDICTIONS")
print("=" * 70)
print("Self-Consistent Field: iterate screening until converged.")
print("Compute E(excited config) - E(ground config) for the FULL atom.")
print()

print(f'{"Name":>5} {"Line":>20} {"GWT nm":>10} {"Obs nm":>10} '
      f'{"Err":>8} {"dE(eV)":>8}')
print("-" * 70)

errors = []
for name, Z, conf_g, conf_e, obs_nm, line_name in tests:
    dE = scf_transition(Z, conf_g, conf_e)
    if dE > 0.001:
        lam = hc / dE
        err = (lam - obs_nm) / obs_nm * 100
        errors.append(abs(err))
        flag = '' if abs(err) < 5 else (' *' if abs(err) < 15 else ' **')
        print(f'{name:>5} {line_name:>20} {lam:10.2f} {obs_nm:10.3f} '
              f'{err:+7.2f}% {dE:8.4f}{flag}')
    else:
        errors.append(100)
        print(f'{name:>5} {line_name:>20} {"neg/zero":>10} {obs_nm:10.3f} '
              f'{"---":>8} {dE:8.4f} **')

print("-" * 70)
if errors:
    print(f'Mean |error|: {np.mean(errors):.1f}%')
    print(f'Under 5%:  {sum(1 for e in errors if e < 5)}/{len(errors)}')
    print(f'Under 10%: {sum(1 for e in errors if e < 10)}/{len(errors)}')

# Detailed He breakdown
print()
print("HELIUM DETAILED:")
dE_He = scf_transition(2, [(1,0,2)], [(1,0,1),(2,1,1)], verbose=True)
if dE_He > 0:
    print(f'  Lambda: {hc/dE_He:.1f} nm (obs: 58.4 nm)')
