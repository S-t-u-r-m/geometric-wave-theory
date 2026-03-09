"""
v23 Corrected -- Systematic parameter search with FIXED data
=============================================================
Re-runs the bond_3d_harmonic_v23.py systematic search with:
  1. Corrected Cl Z_eff = 6.1161 (was 4.8864)
  2. Corrected phase: b = 1 + beta*h (no (1-2l) factor)
  3. |sin(phase)| instead of sin(phase)
  4. All 24 molecules (was 12)
  5. CN BO=2.5, NO pi_anti=1

Goal: see if corrected data reveals different optimal parameters
than what we found before with bad data.
"""

import numpy as np

pi = np.pi
E_H = 13.6057
d = 3

# Corrected Z_eff
Z_eff = {
    'H_1s': 1.0000, 'Li_2s': 1.2792, 'B_2p': 2.4214,
    'C_2p': 3.1358, 'N_2p': 3.8340, 'O_2p': 4.4532,
    'F_2p': 5.0998, 'Na_3s': 2.5074, 'Cl_3p': 6.1161,  # CORRECTED
}

def orbital_energy(orb):
    n = int(orb.split('_')[1][0])
    return E_H * (Z_eff[orb] / n)**2

def get_n(orb): return int(orb.split('_')[1][0])
def get_l(orb): return {'s': 0, 'p': 1}[orb.split('_')[1][1]]

# All 24 molecules with corrected bond orders
molecules = [
    ('H2',   1.401,  4.745, [('ss',1)],                                         'H_1s',  'H_1s'),
    ('Li2',  5.051,  1.056, [('ss',1)],                                         'Li_2s', 'Li_2s'),
    ('B2',   3.005,  3.02,  [('pi',2)],                                         'B_2p',  'B_2p'),
    ('C2',   2.348,  6.32,  [('pi',2),('sp_sigma',1),('sp_sigma_anti',1)],      'C_2p',  'C_2p'),
    ('N2',   2.074,  9.759, [('pp_sigma',1),('pi',2)],                          'N_2p',  'N_2p'),
    ('O2',   2.282,  5.213, [('pp_sigma',1),('pi',2),('pi_anti',1)],            'O_2p',  'O_2p'),
    ('F2',   2.668,  1.660, [('pp_sigma',1),('pi',2),('pi_anti',2)],            'F_2p',  'F_2p'),
    ('Na2',  5.818,  0.746, [('ss',1)],                                         'Na_3s', 'Na_3s'),
    ('Cl2',  3.757,  2.514, [('pp_sigma',1),('pi',2),('pi_anti',2)],            'Cl_3p', 'Cl_3p'),
    ('HF',   1.733,  5.869, [('sp',1)],                                         'H_1s',  'F_2p'),
    ('CO',   2.132, 11.225, [('pp_sigma',1),('pi',2)],                          'C_2p',  'O_2p'),
    ('NO',   2.175,  6.497, [('pp_sigma',1),('pi',2),('pi_anti',1)],            'N_2p',  'O_2p'),
    # Blind test 12
    ('OH',   1.834,  4.392, [('sp',1)],                                         'O_2p',  'H_1s'),
    ('HCl',  2.409,  4.434, [('sp',1)],                                         'H_1s',  'Cl_3p'),
    ('LiH',  3.015,  2.515, [('ss',1)],                                         'Li_2s', 'H_1s'),
    ('LiF',  2.955,  5.939, [('sp',1)],                                         'Li_2s', 'F_2p'),
    ('BH',   2.329,  3.42,  [('sp',1)],                                         'B_2p',  'H_1s'),
    ('CH',   2.116,  3.47,  [('sp',1)],                                         'C_2p',  'H_1s'),
    ('NH',   1.958,  3.57,  [('sp',1)],                                         'N_2p',  'H_1s'),
    ('BF',   2.386,  7.81,  [('pp_sigma',1),('pi',2)],                          'B_2p',  'F_2p'),
    ('CN',   2.214,  7.72,  [('pp_sigma',0.5),('pi',2)],                        'C_2p',  'N_2p'),
    ('NaH',  3.566,  1.97,  [('ss',1)],                                         'Na_3s', 'H_1s'),
    ('NaCl', 4.461,  4.23,  [('sp',1)],                                         'Na_3s', 'Cl_3p'),
    ('H2O',  1.809,  5.117, [('sp',1)],                                         'O_2p',  'H_1s'),
]

f_pi = d**2 / (d**2 + 1)  # 9/10 -- derived from d=3


def compute_full(mol_data, alpha, beta, c_ionic, f_anti_partial, f_anti_full=1.0):
    name, R, De, bonds, orb1, orb2 = mol_data
    n1, n2 = get_n(orb1), get_n(orb2)
    l1, l2 = get_l(orb1), get_l(orb2)
    h1 = min(n1-l1-1, 1)
    h2 = min(n2-l2-1, 1)

    # CORRECTED: b = 1 + beta*h (no (1-2l) factor)
    a1 = 2 + (1-2*l1)*alpha*h1
    a2 = 2 + (1-2*l2)*alpha*h2
    b1 = 1 + beta*h1
    b2 = 1 + beta*h2

    E1 = E_H / n1**a1
    E2 = E_H / n2**a2
    E_scale = np.sqrt(E1 * E2)
    sigma_phase = R/n1**b1 + R/n2**b2

    n_pi_bond = sum(c for b,c in bonds if 'pi' in b and 'anti' not in b)
    n_pi_anti = sum(c for b,c in bonds if 'pi' in b and 'anti' in b)
    is_full = (n_pi_anti >= n_pi_bond) if n_pi_bond > 0 else True

    Dc = 0
    for btype, count in bonds:
        if 'sigma' in btype or btype in ('ss', 'sp'):
            phase = sigma_phase
        else:
            phase = sigma_phase * f_pi

        # CORRECTED: |sin(phase)| instead of sin(phase)
        contrib = (pi/3) * E_scale * abs(np.sin(phase))

        if 'anti' in btype:
            if 'sigma' in btype or is_full:
                f_a = f_anti_full
            else:
                f_a = f_anti_partial
            Dc -= count * f_a * contrib
        else:
            Dc += count * contrib

    eps1 = orbital_energy(orb1)
    eps2 = orbital_energy(orb2)
    delta = abs(eps1 - eps2)
    V = max(abs(Dc), 0.01)
    q = delta / np.sqrt(delta**2 + (2*V)**2) if delta > 0 else 0
    D_ionic = c_ionic * q**2 * 2 * E_H / R

    return Dc + D_ionic


def evaluate(alpha, beta, c_ionic, f_anti_p, verbose=False):
    results = []
    for mol in molecules:
        D = compute_full(mol, alpha, beta, c_ionic, f_anti_p)
        err = (D - mol[2]) / mol[2]
        results.append((mol[0], D, mol[2], err))
    errors = [abs(r[3]) for r in results]
    n5 = sum(1 for e in errors if e < 0.05)
    n10 = sum(1 for e in errors if e < 0.10)
    avg = np.mean(errors) * 100
    if verbose:
        for name, Dc, Dobs, err in results:
            flag = '***' if abs(err) < 0.02 else ' **' if abs(err) < 0.05 else '  *' if abs(err) < 0.10 else ''
            print(f"  {flag} {name:5s} {Dc:7.3f} {Dobs:7.3f} {err:+8.1%}")
    return n5, n10, avg, results


print("=" * 75)
print("  v23 CORRECTED: Systematic search with fixed Cl, phase, |sin|")
print("=" * 75)

# First: show current d=3 formula
print("\n--- Current d=3 formula (alpha=7/10, beta=19/20) ---")
alpha_cur = 1 - f_pi/d  # 7/10
beta_cur = (1 + f_pi)/2  # 19/20
c_ionic_cur = 1/(2*d+1)  # 1/7
f_anti_cur = 2*d/(2*d-1)  # 6/5

n5, n10, avg, _ = evaluate(alpha_cur, beta_cur, c_ionic_cur, f_anti_cur, verbose=True)
print(f"  within 5%: {n5}/24,  within 10%: {n10}/24,  avg = {avg:.1f}%")


# =====================================================================
# SYSTEMATIC SEARCH with expanded candidates
# =====================================================================
print("\n" + "=" * 75)
print("  SYSTEMATIC SEARCH: All geometric parameter combinations")
print("=" * 75)

alpha_candidates = {
    '2/d=2/3':    2/d,
    '(d-1)/d=2/3': (d-1)/d,  # same as 2/d for d=3
    '1-f_pi/d=7/10': 1-f_pi/d,
    '3/4':        3/4,
    'ln2':        np.log(2),
    '1/sqrt2':    1/np.sqrt(2),
    '2/pi':       2/pi,
    '1-1/pi':     1-1/pi,
    'pi/4-1/4':   (pi-1)/4,
    '(d²-1)/d²=8/9': (d**2-1)/d**2,
}

beta_candidates = {
    '1':          1.0,
    '(1+f_pi)/2=19/20': (1+f_pi)/2,
    '(d-1)/d=2/3': (d-1)/d,
    'pi/d=pi/3': pi/d,
    'd/(d+1)=3/4': d/(d+1),
    'f_pi=9/10': f_pi,
    '(2d-1)/(2d)=5/6': (2*d-1)/(2*d),
}

c_ionic_candidates = {
    '0':          0,
    '1/(2d+1)=1/7': 1/(2*d+1),
    '1/(2d)=1/6': 1/(2*d),
    '1/d²=1/9':  1/(d**2),
    '1/(2pi)':   1/(2*pi),
    '1/(d+4)=1/7': 1/(d+4),  # same as 1/7
    '1/(4pi)':   1/(4*pi),
    '2/d!=1/3':  2/6,
}

f_anti_candidates = {
    '1':          1.0,
    '2d/(2d-1)=6/5': 2*d/(2*d-1),
    '(d+1)/d=4/3': (d+1)/d,
    'd/(d-1)=3/2': d/(d-1),
    '1+1/d²=10/9': 1+1/d**2,
    'pi/e':       pi/np.e,
}

results_list = []

for a_name, a_val in alpha_candidates.items():
    for b_name, b_val in beta_candidates.items():
        for c_name, c_val in c_ionic_candidates.items():
            for f_name, f_val in f_anti_candidates.items():
                n5, n10, avg, _ = evaluate(a_val, b_val, c_val, f_val)
                results_list.append((n5, n10, avg, a_name, b_name, c_name, f_name,
                                      a_val, b_val, c_val, f_val))

results_list.sort(key=lambda x: (-x[0], -x[1], x[2]))

print(f"\n  Top 30 combinations (sorted by w5, then w10, then avg):")
print(f"  {'w5':>3s} {'w10':>3s} {'avg':>6s}  {'alpha':>20s} {'beta':>22s} {'c_ionic':>15s} {'f_anti':>18s}")
print("  " + "-" * 95)

seen = set()
count = 0
for n5, n10, avg, an, bn, cn, fn, av, bv, cv, fv in results_list:
    key = (round(av,6), round(bv,6), round(cv,6), round(fv,6))
    if key in seen:
        continue
    seen.add(key)
    cur = " <-- CURRENT" if abs(av-alpha_cur)<0.001 and abs(bv-beta_cur)<0.001 else ""
    print(f"  {n5:3d} {n10:3d} {avg:5.1f}%  {an:>20s} {bn:>22s} {cn:>15s} {fn:>18s}{cur}")
    count += 1
    if count >= 30:
        break


# =====================================================================
# FINE-GRAINED SCAN around current optimum
# =====================================================================
print("\n" + "=" * 75)
print("  FINE SCAN: alpha and beta near current optimum")
print("=" * 75)

best_n5 = 0
best_result = None

for a_val in np.arange(0.60, 0.80, 0.01):
    for b_val in np.arange(0.85, 1.05, 0.01):
        n5, n10, avg, _ = evaluate(a_val, b_val, c_ionic_cur, f_anti_cur)
        if n5 > best_n5 or (n5 == best_n5 and avg < best_result[2]):
            best_n5 = n5
            best_result = (n5, n10, avg, a_val, b_val)

n5, n10, avg, a_best, b_best = best_result
print(f"  Best: alpha={a_best:.2f}, beta={b_best:.2f} -> w5={n5}, w10={n10}, avg={avg:.1f}%")
print(f"\n  Detail:")
evaluate(a_best, b_best, c_ionic_cur, f_anti_cur, verbose=True)
print(f"  w5={n5}, w10={n10}, avg={avg:.1f}%")


# =====================================================================
# COMPARE: old (wrong) vs new (corrected) for the 12 original molecules
# =====================================================================
print("\n" + "=" * 75)
print("  COMPARE: Impact of corrections on original 12 molecules")
print("=" * 75)

# Old formula (wrong Cl, (1-2l) in beta, sin not |sin|)
Z_eff_old = dict(Z_eff)
Z_eff_old['Cl_3p'] = 4.8864  # wrong value

def compute_old(mol_data, alpha, beta, c_ionic, f_anti):
    """Old formula: (1-2l) in beta, sin() not |sin|, old Cl"""
    name, R, De, bonds, orb1, orb2 = mol_data
    n1, n2 = get_n(orb1), get_n(orb2)
    l1, l2 = get_l(orb1), get_l(orb2)
    h1, h2 = min(n1-l1-1,1), min(n2-l2-1,1)

    a1 = 2 + (1-2*l1)*alpha*h1
    a2 = 2 + (1-2*l2)*alpha*h2
    b1 = 1 + (1-2*l1)*beta*h1  # OLD: with (1-2l)
    b2 = 1 + (1-2*l2)*beta*h2

    z1_old = Z_eff_old[orb1]
    z2_old = Z_eff_old[orb2]

    E1 = E_H * (z1_old/n1)**2 / (z1_old/n1)**(a1-2) if a1 != 2 else E_H*(z1_old/n1)**2
    # Actually the old formula uses E_H/n^a, not Z_eff in E_scale
    E1 = E_H / n1**a1
    E2 = E_H / n2**a2
    E_scale = np.sqrt(E1 * E2)
    sp = R/n1**b1 + R/n2**b2

    n_pb = sum(c for b,c in bonds if 'pi' in b and 'anti' not in b)
    n_pa = sum(c for b,c in bonds if 'pi' in b and 'anti' in b)
    is_full = (n_pa >= n_pb) if n_pb > 0 else True

    Dc = 0
    for btype, count in bonds:
        phase = sp if ('sigma' in btype or btype in ('ss','sp')) else sp * f_pi
        contrib = (pi/3) * E_scale * np.sin(phase)  # OLD: sin, not |sin|
        if 'anti' in btype:
            fa = 1.0 if ('sigma' in btype or is_full) else f_anti
            Dc -= count * fa * contrib
        else:
            Dc += count * contrib

    # Old uses Z_eff_old for orbital energies
    eps1 = E_H * (z1_old/n1)**2
    eps2 = E_H * (z2_old/n2)**2
    delta = abs(eps1 - eps2)
    V = max(abs(Dc), 0.01)
    q = delta / np.sqrt(delta**2 + (2*V)**2) if delta > 0 else 0
    D_ion = c_ionic * q**2 * 2 * E_H / R

    return Dc + D_ion

print(f"\n  {'Mol':<7} {'OLD':>10} {'NEW':>10} {'exp':>8}  {'old_err':>8} {'new_err':>8} {'change':>8}")
print("  " + "-" * 70)

for mol in molecules[:12]:  # original 12 only
    name = mol[0]
    D_old = compute_old(mol, alpha_cur, beta_cur, c_ionic_cur, f_anti_cur)
    D_new = compute_full(mol, alpha_cur, beta_cur, c_ionic_cur, f_anti_cur)
    De = mol[2]
    e_old = (D_old - De)/De * 100
    e_new = (D_new - De)/De * 100
    better = "BETTER" if abs(e_new) < abs(e_old) else "worse" if abs(e_new) > abs(e_old) else "same"
    print(f"  {name:<7} {D_old:10.3f} {D_new:10.3f} {De:8.3f}  {e_old:+7.1f}% {e_new:+7.1f}% {better:>8}")


# =====================================================================
# What if we use the OLD (1-2l) factor in phase with corrected Cl?
# =====================================================================
print("\n" + "=" * 75)
print("  TEST: Old (1-2l) in phase but with corrected Cl and |sin|")
print("=" * 75)

def compute_old_phase(mol_data, alpha, beta, c_ionic, f_anti):
    """Old phase (1-2l)*beta*h, but corrected Cl and |sin|"""
    name, R, De, bonds, orb1, orb2 = mol_data
    n1, n2 = get_n(orb1), get_n(orb2)
    l1, l2 = get_l(orb1), get_l(orb2)
    h1, h2 = min(n1-l1-1,1), min(n2-l2-1,1)

    a1 = 2 + (1-2*l1)*alpha*h1
    a2 = 2 + (1-2*l2)*alpha*h2
    b1 = 1 + (1-2*l1)*beta*h1  # OLD: with (1-2l)
    b2 = 1 + (1-2*l2)*beta*h2

    E1 = E_H / n1**a1
    E2 = E_H / n2**a2
    E_scale = np.sqrt(E1 * E2)
    sp = R/n1**b1 + R/n2**b2

    n_pb = sum(c for b,c in bonds if 'pi' in b and 'anti' not in b)
    n_pa = sum(c for b,c in bonds if 'pi' in b and 'anti' in b)
    is_full = (n_pa >= n_pb) if n_pb > 0 else True

    Dc = 0
    for btype, count in bonds:
        phase = sp if ('sigma' in btype or btype in ('ss','sp')) else sp * f_pi
        contrib = (pi/3) * E_scale * abs(np.sin(phase))  # CORRECTED: |sin|
        if 'anti' in btype:
            fa = 1.0 if ('sigma' in btype or is_full) else f_anti
            Dc -= count * fa * contrib
        else:
            Dc += count * contrib

    eps1 = orbital_energy(orb1)
    eps2 = orbital_energy(orb2)
    delta = abs(eps1 - eps2)
    V = max(abs(Dc), 0.01)
    q = delta / np.sqrt(delta**2 + (2*V)**2) if delta > 0 else 0
    D_ion = c_ionic * q**2 * 2 * E_H / R

    return Dc + D_ion

errs_old_phase = []
errs_new_phase = []
print(f"\n  {'Mol':<7} {'old_ph':>10} {'new_ph':>10} {'exp':>8}  {'old_err':>8} {'new_err':>8}")
print("  " + "-" * 60)

for mol in molecules:
    name = mol[0]
    D_op = compute_old_phase(mol, alpha_cur, beta_cur, c_ionic_cur, f_anti_cur)
    D_np = compute_full(mol, alpha_cur, beta_cur, c_ionic_cur, f_anti_cur)
    De = mol[2]
    e_op = (D_op - De)/De * 100
    e_np = (D_np - De)/De * 100
    errs_old_phase.append(abs(e_op))
    errs_new_phase.append(abs(e_np))
    better = "<-" if abs(e_op) < abs(e_np) - 1 else "->" if abs(e_np) < abs(e_op) - 1 else ""
    print(f"  {name:<7} {D_op:10.3f} {D_np:10.3f} {De:8.3f}  {e_op:+7.1f}% {e_np:+7.1f}% {better}")

print(f"\n  Old phase: avg={np.mean(errs_old_phase):.1f}%, w5={sum(1 for e in errs_old_phase if e<5)}, w10={sum(1 for e in errs_old_phase if e<10)}")
print(f"  New phase: avg={np.mean(errs_new_phase):.1f}%, w5={sum(1 for e in errs_new_phase if e<5)}, w10={sum(1 for e in errs_new_phase if e<10)}")
