"""
HONEST MODEL: What does GWT actually predict from d=3 alone?

No fitting, no scanning, no numerology.
Just the derived constants and simple physical arguments.

DERIVED from d=3:
  C_bond = pi/3, f_pi = 9/10, alpha = 7/10, beta = 19/20
  f_anti = 6/5, c_ionic = 1/7

PHYSICAL ARGUMENT for |sin|:
  sin(phase) is the overlap integral of two standing waves.
  When phase > pi, there is a node between the atoms.
  But the BONDING energy depends on the overlap AMPLITUDE,
  not its sign. The sign determines bonding vs antibonding
  for a GIVEN orbital, but we already specify that via bond_list.
  Therefore: use |sin(phase)| for the overlap strength.
  No damping factor needed if the physical argument is correct.

This gives ONE formula with ZERO free parameters:
  D = (pi/3) * sqrt(E1*E2) * sum[|sin(phase)|] + (1/7)*q^2*2E_H/R
"""

import numpy as np

pi = np.pi
E_H = 13.6057
d = 3
C_bond = pi / d
f_pi = d**2 / (d**2 + 1)
alpha = 1 - f_pi / d
beta = (1 + f_pi) / 2
f_anti = 2*d / (2*d - 1)
c_ionic = 1.0 / (2*d + 1)

Z_eff = {
    'H_1s': 1.0, 'Li_2s': 1.2792, 'B_2p': 2.4214,
    'C_2p': 3.1358, 'N_2p': 3.8340, 'O_2p': 4.4532,
    'F_2p': 5.0998, 'Na_3s': 2.5074, 'Cl_3p': 6.1161,
}

def get_n(o): return int(o.split('_')[1][0])
def get_l(o): return {'s': 0, 'p': 1}[o.split('_')[1][1]]
def has_nodes(o): return min(get_n(o) - get_l(o) - 1, 1)
def orb_e(o): return E_H * (Z_eff[o] / get_n(o))**2
def params(o):
    n, l = get_n(o), get_l(o)
    h = has_nodes(o)
    return n, l, h, 2 + (1 - 2*l)*alpha*h, 1 + (1 - 2*l)*beta*h

molecules = [
    ('H2',   1.401, 4.745, [('ss', 1)],                                        'H_1s',  'H_1s'),
    ('Li2',  5.051, 1.056, [('ss', 1)],                                        'Li_2s', 'Li_2s'),
    ('B2',   3.005, 3.02,  [('pi', 2)],                                        'B_2p',  'B_2p'),
    ('C2',   2.348, 6.32,  [('pi', 2), ('sp_sigma', 1), ('sp_sigma_anti', 1)], 'C_2p',  'C_2p'),
    ('N2',   2.074, 9.759, [('pp_sigma', 1), ('pi', 2)],                       'N_2p',  'N_2p'),
    ('O2',   2.282, 5.213, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)],       'O_2p',  'O_2p'),
    ('F2',   2.668, 1.660, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)],       'F_2p',  'F_2p'),
    ('Na2',  5.818, 0.746, [('ss', 1)],                                        'Na_3s', 'Na_3s'),
    ('Cl2',  3.757, 2.514, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)],       'Cl_3p', 'Cl_3p'),
    ('HF',   1.733, 5.869, [('sp', 1)],                                        'H_1s',  'F_2p'),
    ('CO',   2.132, 11.225,[('pp_sigma', 1), ('pi', 2)],                       'C_2p',  'O_2p'),
    ('NO',   2.175, 6.497, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)],       'N_2p',  'O_2p'),
    ('OH',   1.834, 4.392, [('sp', 1)],                                        'O_2p',  'H_1s'),
    ('HCl',  2.409, 4.434, [('sp', 1)],                                        'H_1s',  'Cl_3p'),
    ('LiH',  3.015, 2.515, [('ss', 1)],                                        'Li_2s', 'H_1s'),
    ('LiF',  2.955, 5.939, [('sp', 1)],                                        'Li_2s', 'F_2p'),
    ('BH',   2.329, 3.42,  [('sp', 1)],                                        'B_2p',  'H_1s'),
    ('CH',   2.116, 3.47,  [('sp', 1)],                                        'C_2p',  'H_1s'),
    ('NH',   1.958, 3.57,  [('sp', 1)],                                        'N_2p',  'H_1s'),
    ('BF',   2.386, 7.81,  [('pp_sigma', 1), ('pi', 2)],                       'B_2p',  'F_2p'),
    ('CN',   2.214, 7.72,  [('pp_sigma', 1), ('pi', 2)],                       'C_2p',  'N_2p'),
    ('NaH',  3.566, 1.97,  [('ss', 1)],                                        'Na_3s', 'H_1s'),
    ('NaCl', 4.461, 4.23,  [('sp', 1)],                                        'Na_3s', 'Cl_3p'),
    ('H2O',  1.809, 5.117, [('sp', 1)],                                        'O_2p',  'H_1s'),
]


def compute_energy(mol, use_abs_sin=False):
    nm, R, De, bonds, o1, o2 = mol
    n1, l1, h1, a1, b1 = params(o1)
    n2, l2, h2, a2, b2 = params(o2)
    Es = np.sqrt(E_H / n1**a1 * E_H / n2**a2)
    sp = R / n1**b1 + R / n2**b2

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    Dc = 0
    for bt, cnt in bonds:
        ph = sp if ('sigma' in bt or bt in ('ss', 'sp')) else sp * f_pi
        if use_abs_sin:
            sv = abs(np.sin(ph))
        else:
            sv = np.sin(ph)
        cont = C_bond * Es * sv
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            Dc -= cnt * fa * cont
        else:
            Dc += cnt * cont

    e1, e2 = orb_e(o1), orb_e(o2)
    de = abs(e1 - e2)
    Vc = max(abs(Dc), 0.01)
    q = de / np.sqrt(de**2 + (2*Vc)**2) if de > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / R
    return Dc + Di, Dc, Di, q, sp


# =============================================================================
# MODEL A: Original sin(phase) -- current formula
# =============================================================================
print("=" * 90)
print("  MODEL A: sin(phase)  [current formula]")
print("  D = (pi/3) * sqrt(E1*E2) * sin(phase) + (1/7)*q^2*2E_H/R")
print("=" * 90)
print()

print(f"{'Mol':<7} {'De_exp':>7} {'D_pred':>7} {'err%':>7} {'D_cov':>7} {'D_ion':>7} {'ph/pi':>6}")
print("-" * 55)

errs_A = []
for mol in molecules:
    nm, De = mol[0], mol[2]
    Dt, Dc, Di, q, sp = compute_energy(mol, use_abs_sin=False)
    err = (Dt - De) / De * 100
    errs_A.append(abs(err))
    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
    if Dt < 0 and De > 0:
        flag = ' NEG'
    print(f"{nm:<7} {De:7.3f} {Dt:+7.3f} {err:+6.1f}% {Dc:+7.3f} {Di:7.3f} {sp/pi:6.3f} {flag}")

neg_A = sum(1 for mol in molecules if compute_energy(mol, False)[0] < 0 and mol[2] > 0)
w5_A = sum(1 for e in errs_A if e < 5)
w10_A = sum(1 for e in errs_A if e < 10)
print(f"\nAvg = {np.mean(errs_A):.1f}%, <5%: {w5_A}/24, <10%: {w10_A}/24, D<0: {neg_A}/24")


# =============================================================================
# MODEL B: |sin(phase)| -- absolute value everywhere
# =============================================================================
print()
print("=" * 90)
print("  MODEL B: |sin(phase)|  [overlap amplitude]")
print("  D = (pi/3) * sqrt(E1*E2) * |sin(phase)| + (1/7)*q^2*2E_H/R")
print("  Physical argument: overlap energy = amplitude, not signed integral")
print("=" * 90)
print()

print(f"{'Mol':<7} {'De_exp':>7} {'D_pred':>7} {'err%':>7} {'D_cov':>7} {'D_ion':>7} {'ph/pi':>6}")
print("-" * 55)

errs_B = []
for mol in molecules:
    nm, De = mol[0], mol[2]
    Dt, Dc, Di, q, sp = compute_energy(mol, use_abs_sin=True)
    err = (Dt - De) / De * 100
    errs_B.append(abs(err))
    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
    print(f"{nm:<7} {De:7.3f} {Dt:+7.3f} {err:+6.1f}% {Dc:+7.3f} {Di:7.3f} {sp/pi:6.3f} {flag}")

neg_B = sum(1 for mol in molecules if compute_energy(mol, True)[0] < 0 and mol[2] > 0)
w5_B = sum(1 for e in errs_B if e < 5)
w10_B = sum(1 for e in errs_B if e < 10)
w15_B = sum(1 for e in errs_B if e < 15)
print(f"\nAvg = {np.mean(errs_B):.1f}%, <5%: {w5_B}/24, <10%: {w10_B}/24, <15%: {w15_B}/24, D<0: {neg_B}/24")


# =============================================================================
# HEAD TO HEAD
# =============================================================================
print()
print("=" * 90)
print("  HEAD TO HEAD: sin vs |sin|")
print("=" * 90)
print()

print(f"{'Mol':<7} {'err_sin':>8} {'err_|sin|':>9} {'better':>7} {'ph/pi':>6}")
print("-" * 45)

wins_A = 0; wins_B = 0; ties = 0
for i, mol in enumerate(molecules):
    nm = mol[0]
    _, _, _, _, sp = compute_energy(mol, False)
    eA = errs_A[i]
    eB = errs_B[i]
    if eA < eB - 1:
        better = 'sin'
        wins_A += 1
    elif eB < eA - 1:
        better = '|sin|'
        wins_B += 1
    else:
        better = '~'
        ties += 1
    print(f"{nm:<7} {eA:7.1f}% {eB:8.1f}% {better:>7} {sp/pi:6.3f}")

print(f"\nsin wins: {wins_A}, |sin| wins: {wins_B}, ties: {ties}")
print()

# =============================================================================
# WHAT MODEL B GETS WRONG AND WHY
# =============================================================================
print("=" * 90)
print("  ANALYSIS: Where does |sin| fail?")
print("=" * 90)
print()
print("Molecules where |sin| error > 20%:")
print()
for i, mol in enumerate(molecules):
    if errs_B[i] > 20:
        nm, R, De = mol[0], mol[1], mol[2]
        Dt, Dc, Di, q, sp = compute_energy(mol, use_abs_sin=True)
        err = (Dt - De) / De * 100
        e1, e2 = orb_e(mol[4]), orb_e(mol[5])
        E_ratio = max(e1, e2) / min(e1, e2)
        cov_pct = Dc / Dt * 100 if Dt > 0.01 else 0
        ion_pct = Di / Dt * 100 if Dt > 0.01 else 0
        print(f"  {nm}: err={err:+.1f}%, E_ratio={E_ratio:.1f}, "
              f"cov={cov_pct:.0f}%, ion={ion_pct:.0f}%, ph/pi={sp/pi:.3f}")
