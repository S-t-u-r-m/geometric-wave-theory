"""Test GWT spectral predictions against precisely known lines."""
import sys, os, io
sys.path.insert(0, os.path.dirname(__file__))

# Import BEFORE wrapping stdout (atomic_data may trigger a wrapper)
import numpy as np
from z_eff_subshell import compute_all_levels, ATOMS

# Ensure stdout works
try:
    sys.stdout.write('')
except:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

hc = 1240.0
l_names = {0: 's', 1: 'p', 2: 'd', 3: 'f'}

# Well-known spectral lines (strongest transitions)
test_lines = [
    ('H',  2,1, 1,0,  121.567, 'Lyman-alpha'),
    ('H',  3,1, 1,0,  102.573, 'Lyman-beta'),
    ('H',  3,1, 2,0,  656.281, 'H-alpha'),
    ('H',  4,1, 2,0,  486.135, 'H-beta'),
    ('H',  5,1, 2,0,  434.047, 'H-gamma'),
    ('H',  6,1, 2,0,  410.174, 'H-delta'),
    ('He', 2,1, 1,0,   58.433, 'He I resonance'),
    ('Li', 2,1, 2,0,  670.776, 'Li I red doublet'),
    ('Na', 3,1, 3,0,  589.592, 'Na D1'),
    ('K',  4,1, 4,0,  766.490, 'K I principal'),
    ('Ca', 4,1, 4,0,  422.673, 'Ca I resonance'),
    ('Mg', 3,1, 3,0,  285.213, 'Mg I resonance'),
]

print('GWT SPECTRAL LINE PREDICTIONS vs OBSERVED')
print('=' * 78)
print(f'{"Elem":>5} {"Transition":>10} {"GWT nm":>10} {"Obs nm":>10} '
      f'{"Err":>8} {"Z":>3} {"Name"}')
print('-' * 78)

errors = []
for element, n_hi, l_hi, n_lo, l_lo, obs_nm, name in test_lines:
    Z = ATOMS[element]['Z']
    levels = compute_all_levels(element, max_n=max(n_hi, n_lo) + 2)

    key_hi = (n_hi, l_hi)
    key_lo = (n_lo, l_lo)

    if key_hi not in levels or key_lo not in levels:
        print(f'{element:>5} {n_hi}{l_names[l_hi]}->{n_lo}{l_names[l_lo]:>6} '
              f'{"SKIP":>10} {obs_nm:10.3f} {"---":>8} {Z:3d} {name}')
        continue

    dE = levels[key_hi]['E'] - levels[key_lo]['E']
    if dE <= 0.001:
        print(f'{element:>5} {n_hi}{l_names[l_hi]}->{n_lo}{l_names[l_lo]:>6} '
              f'{"dE<=0":>10} {obs_nm:10.3f} {"---":>8} {Z:3d} {name}')
        continue

    lam = hc / dE
    err = (lam - obs_nm) / obs_nm * 100
    errors.append((Z, element, abs(err)))
    flag = '' if abs(err) < 2 else (' *' if abs(err) < 10 else ' **')
    print(f'{element:>5} {n_hi}{l_names[l_hi]}->{n_lo}{l_names[l_lo]:>6} '
          f'{lam:10.2f} {obs_nm:10.3f} {err:+7.2f}% {Z:3d} {name}{flag}')

print('-' * 78)

# Summary by Z
print(f'\nERROR TREND vs ATOMIC NUMBER:')
for Z, elem, err in sorted(errors):
    bar = '*' * min(50, int(err * 2))
    print(f'  Z={Z:2d} ({elem:>2s}): {err:6.2f}% {bar}')

mean_err = np.mean([e[2] for e in errors])
print(f'\nMean |error|: {mean_err:.2f}%')
print(f'Under 1%:  {sum(1 for _,_,e in errors if e < 1)}/{len(errors)}')
print(f'Under 5%:  {sum(1 for _,_,e in errors if e < 5)}/{len(errors)}')
print(f'Under 10%: {sum(1 for _,_,e in errors if e < 10)}/{len(errors)}')
