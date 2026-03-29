"""Test spectral lines using TOTAL atomic energy difference."""
import sys, os, io
sys.path.insert(0, os.path.dirname(__file__))
from z_eff_subshell import compute_subshell_energy, ATOMS
import numpy as np

try:
    sys.stdout.write('')
except:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

E_H = 13.6045
hc = 1240.0
l_n = {0:'s', 1:'p', 2:'d', 3:'f'}


def total_energy(Z, config):
    """Total atomic energy = sum of all subshell energies."""
    E_tot = 0
    for n, l, count in config:
        E_sh, *_ = compute_subshell_energy(Z, config, n, l)
        E_tot += count * E_sh
    return E_tot


def transition_energy(Z, config_ground, config_excited):
    """Energy difference between excited and ground configs."""
    return total_energy(Z, config_excited) - total_energy(Z, config_ground)


print("TOTAL ENERGY SPECTRAL PREDICTIONS")
print("=" * 65)
print("Compute E(excited config) - E(ground config) for the FULL atom.")
print()

tests = [
    # (name, Z, ground_config, excited_config, obs_nm, line_name)
    ("H", 1,
     [(1,0,1)],
     [(2,1,1)],
     121.567, "Lyman-alpha"),

    ("H", 1,
     [(1,0,1)],
     [(3,1,1)],
     656.281, "H-alpha (from n=1)"),

    ("He", 2,
     [(1,0,2)],
     [(1,0,1), (2,1,1)],
     58.433, "He I resonance"),

    ("Li", 3,
     [(1,0,2), (2,0,1)],
     [(1,0,2), (2,1,1)],
     670.776, "Li I red"),

    ("Na", 11,
     [(1,0,2),(2,0,2),(2,1,6),(3,0,1)],
     [(1,0,2),(2,0,2),(2,1,6),(3,1,1)],
     589.592, "Na D1"),

    ("K", 19,
     [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(4,0,1)],
     [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(4,1,1)],
     766.490, "K I principal"),

    ("Mg", 12,
     [(1,0,2),(2,0,2),(2,1,6),(3,0,2)],
     [(1,0,2),(2,0,2),(2,1,6),(3,0,1),(3,1,1)],
     285.213, "Mg I resonance"),

    ("Ca", 20,
     [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(4,0,2)],
     [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(4,0,1),(4,1,1)],
     422.673, "Ca I resonance"),
]

print(f'{"Name":>5} {"Line":>20} {"GWT nm":>10} {"Obs nm":>10} {"Err":>8} {"dE(eV)":>8}')
print("-" * 65)

for name, Z, conf_g, conf_e, obs_nm, line_name in tests:
    dE = transition_energy(Z, conf_g, conf_e)
    if dE > 0.001:
        lam = hc / dE
        err = (lam - obs_nm) / obs_nm * 100
        flag = '' if abs(err) < 5 else (' *' if abs(err) < 15 else ' **')
        print(f'{name:>5} {line_name:>20} {lam:10.2f} {obs_nm:10.3f} '
              f'{err:+7.2f}% {dE:8.4f}{flag}')
    else:
        print(f'{name:>5} {line_name:>20} {"neg/zero":>10} {obs_nm:10.3f} '
              f'{"---":>8} {dE:8.4f} **')

print("-" * 65)
