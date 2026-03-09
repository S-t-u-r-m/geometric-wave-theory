"""
PHASE FIX TEST: Decouple energy and phase sign conventions
===========================================================
Hypothesis: The (1-2l) sign flip is correct for ENERGY but wrong for PHASE.

Energy: p-orbitals with nodes have MORE KE (centrifugal) -> decrease a -> correct
Phase:  nodes always EXTEND the orbital spatially -> increase b -> sign should NOT flip

Test: keep a = 2 + (1-2l)*alpha*h for energy
      use  b = 1 + beta*h for phase (no sign dependence on l)
      or   b = 1 + |1-2l|*beta*h (same result since |1-2l|=1 always)
"""
import numpy as np

pi = np.pi
E_H = 13.6057
d = 3
C_bond = pi / d
f_pi = d**2 / (d**2 + 1)
alpha_n = 1 - f_pi / d       # 7/10
beta_n = (1 + f_pi) / 2      # 19/20
f_anti = 2*d / (2*d - 1)     # 6/5
c_ionic = 1.0 / (2*d + 1)    # 1/7

Z_eff = {
    'H_1s': 1.0, 'Li_2s': 1.2792, 'B_2p': 2.4214,
    'C_2p': 3.1358, 'N_2p': 3.8340, 'O_2p': 4.4532,
    'F_2p': 5.0998, 'Na_3s': 2.5074, 'Cl_3p': 6.1161,
}

molecules = [
    ('H2',   1.401,  4.745, [('ss', 1)], 'H_1s',  'H_1s'),
    ('Li2',  5.051,  1.056, [('ss', 1)], 'Li_2s', 'Li_2s'),
    ('B2',   3.005,  3.02,  [('pi', 2)], 'B_2p',  'B_2p'),
    ('C2',   2.348,  6.32,  [('pi', 2), ('sp_sigma', 1), ('sp_sigma_anti', 1)], 'C_2p',  'C_2p'),
    ('N2',   2.074,  9.759, [('pp_sigma', 1), ('pi', 2)], 'N_2p',  'N_2p'),
    ('O2',   2.282,  5.213, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)], 'O_2p',  'O_2p'),
    ('F2',   2.668,  1.660, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)], 'F_2p',  'F_2p'),
    ('Na2',  5.818,  0.746, [('ss', 1)], 'Na_3s', 'Na_3s'),
    ('Cl2',  3.757,  2.514, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)], 'Cl_3p', 'Cl_3p'),
    ('HF',   1.733,  5.869, [('sp', 1)], 'H_1s',  'F_2p'),
    ('CO',   2.132, 11.225, [('pp_sigma', 1), ('pi', 2)], 'C_2p',  'O_2p'),
    ('NO',   2.175,  6.497, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)], 'N_2p',  'O_2p'),
    ('OH',   1.834,  4.392, [('sp', 1)], 'O_2p',  'H_1s'),
    ('HCl',  2.409,  4.434, [('sp', 1)], 'H_1s',  'Cl_3p'),
    ('LiH',  3.015,  2.515, [('ss', 1)], 'Li_2s', 'H_1s'),
    ('LiF',  2.955,  5.939, [('sp', 1)], 'Li_2s', 'F_2p'),
    ('BH',   2.329,  3.42,  [('sp', 1)], 'B_2p',  'H_1s'),
    ('CH',   2.116,  3.47,  [('sp', 1)], 'C_2p',  'H_1s'),
    ('NH',   1.958,  3.57,  [('sp', 1)], 'N_2p',  'H_1s'),
    ('BF',   2.386,  7.81,  [('pp_sigma', 1), ('pi', 2)], 'B_2p',  'F_2p'),
    ('CN',   2.214,  7.72,  [('pp_sigma', 1), ('pi', 2)], 'C_2p',  'N_2p'),
    ('NaH',  3.566,  1.97,  [('ss', 1)], 'Na_3s', 'H_1s'),
    ('NaCl', 4.461,  4.23,  [('sp', 1)], 'Na_3s', 'Cl_3p'),
    ('H2O',  1.809,  5.117, [('sp', 1)], 'O_2p',  'H_1s'),
]


def get_n(orb): return int(orb.split('_')[1][0])
def get_l(orb): return {'s': 0, 'p': 1}[orb.split('_')[1][1]]
def has_nodes(orb): return min(get_n(orb) - get_l(orb) - 1, 1)
def orb_energy(orb): return E_H * (Z_eff[orb] / get_n(orb))**2


def compute(mol, phase_model='old'):
    """
    phase_model:
      'old': b = 1 + (1-2l)*beta*h  (current formula)
      'abs': b = 1 + |1-2l|*beta*h = 1 + beta*h  (always positive node correction)
      'no_node': b = 1  (no node correction at all)
    """
    name, R, De_exp, bonds, orb1, orb2 = mol

    for orb in [orb1, orb2]:
        n = get_n(orb)
        l = get_l(orb)
        h = has_nodes(orb)

    # Energy exponents: keep (1-2l) sign flip (physically correct for KE)
    n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
    n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)
    a1 = 2 + (1 - 2*l1) * alpha_n * h1
    a2 = 2 + (1 - 2*l2) * alpha_n * h2

    # Phase exponents: depends on model
    if phase_model == 'old':
        b1 = 1 + (1 - 2*l1) * beta_n * h1
        b2 = 1 + (1 - 2*l2) * beta_n * h2
    elif phase_model == 'abs':
        b1 = 1 + beta_n * h1  # always positive (|1-2l| = 1)
        b2 = 1 + beta_n * h2
    elif phase_model == 'no_node':
        b1 = 1.0
        b2 = 1.0
    else:
        raise ValueError(f"Unknown phase_model: {phase_model}")

    E_scale = np.sqrt(E_H / n1**a1 * E_H / n2**a2)
    sigma_phase = R * (1.0 / n1**b1 + 1.0 / n2**b2)

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss', 'sp')) else sigma_phase * f_pi
        cont = C_bond * E_scale * abs(np.sin(ph))
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    eps1 = E_H * (Z_eff[orb1] / n1)**2
    eps2 = E_H * (Z_eff[orb2] / n2)**2
    de = abs(eps1 - eps2)
    Vc = max(abs(D_cov), 0.01)
    q = de / np.sqrt(de**2 + (2*Vc)**2) if de > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / R

    return D_cov + Di, D_cov, Di, sigma_phase, q


# =============================================================================
# Show orbital parameters for each model
# =============================================================================
print("=" * 85)
print("  ORBITAL PARAMETERS: old vs fixed phase exponent")
print("=" * 85)
print()
print(f"{'orbital':>8} {'n':>2} {'l':>2} {'h':>2} | {'a':>5} | {'b_old':>6} {'k_old':>7} | {'b_fix':>6} {'k_fix':>7} | {'rB':>6}")
print("-" * 85)

for orb in sorted(Z_eff.keys()):
    n = get_n(orb)
    l = get_l(orb)
    h = has_nodes(orb)
    a = 2 + (1 - 2*l) * alpha_n * h
    b_old = 1 + (1 - 2*l) * beta_n * h
    b_fix = 1 + beta_n * h
    k_old = 1.0 / n**b_old
    k_fix = 1.0 / n**b_fix
    rB = n**2 / Z_eff[orb]
    print(f"{orb:>8} {n:2d} {l:2d} {h:2d} | {a:5.2f} | {b_old:6.3f} {k_old:7.4f} | {b_fix:6.3f} {k_fix:7.4f} | {rB:6.3f}")

print()
print("Key change: Cl_3p goes from k=0.947 to k=0.117 (now matches Na_3s scale)")
print("            This means Cl_3p is treated as a LARGE orbital (correct!)")


# =============================================================================
# Compare models
# =============================================================================
models = ['old', 'abs', 'no_node']
model_labels = {
    'old': 'OLD: b=1+(1-2l)*beta*h',
    'abs': 'FIX: b=1+beta*h (no sign flip)',
    'no_node': 'NONE: b=1 (no node correction)',
}

for model in models:
    print()
    print("=" * 95)
    print(f"  {model_labels[model]}")
    print("=" * 95)
    print()
    print(f"{'Mol':<6} {'De_exp':>7} {'D_tot':>7} {'err%':>7} {'D_cov':>7} {'D_ion':>7} {'ph/pi':>6} {'q':>5}")
    print("-" * 60)

    errs = []
    for mol in molecules:
        name = mol[0]
        De_exp = mol[2]
        D_tot, D_cov, D_ion, phase, q = compute(mol, phase_model=model)
        err = (D_tot - De_exp) / De_exp * 100
        errs.append(abs(err))
        flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
        print(f"{name:<6} {De_exp:7.3f} {D_tot:7.3f} {err:+6.1f}% {D_cov:7.3f} {D_ion:7.3f} "
              f"{phase/pi:6.3f} {q:5.3f} {flag}")

    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"\nAvg={np.mean(errs):.1f}%, Med={np.median(errs):.1f}%, "
          f"<5%:{w5}/24, <10%:{w10}/24, <20%:{w20}/24")


# =============================================================================
# Head-to-head comparison: old vs fix
# =============================================================================
print()
print("=" * 95)
print("  HEAD TO HEAD: old vs fix")
print("=" * 95)
print()
print(f"{'Mol':<6} {'De_exp':>7} {'e_old':>7} {'e_fix':>7} {'verdict':>8} {'ph_old':>7} {'ph_fix':>7}")
print("-" * 60)

n_better = 0
n_worse = 0
n_same = 0
for mol in molecules:
    name = mol[0]
    De_exp = mol[2]
    D_old = compute(mol, 'old')
    D_fix = compute(mol, 'abs')
    e_old = abs((D_old[0] - De_exp) / De_exp * 100)
    e_fix = abs((D_fix[0] - De_exp) / De_exp * 100)

    if e_fix < e_old - 1:
        v = 'BETTER'
        n_better += 1
    elif e_fix > e_old + 1:
        v = 'WORSE'
        n_worse += 1
    else:
        v = '~'
        n_same += 1

    print(f"{name:<6} {De_exp:7.3f} {e_old:6.1f}% {e_fix:6.1f}% {v:>8} "
          f"{D_old[3]/pi:7.3f} {D_fix[3]/pi:7.3f}")

print(f"\nBetter: {n_better}, Worse: {n_worse}, Same: {n_same}")


# =============================================================================
# What if we also fix the energy exponent?
# =============================================================================
print()
print("=" * 95)
print("  VARIANT: What if energy exponent ALSO uses |1-2l|?")
print("  a = 2 + alpha*h, b = 1 + beta*h (no sign flip for either)")
print("=" * 95)
print()

print(f"{'Mol':<6} {'De_exp':>7} {'D_tot':>7} {'err%':>7} {'ph/pi':>6}")
print("-" * 40)

errs_both = []
for mol in molecules:
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
    n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)

    # BOTH exponents use |1-2l| = 1 (no sign flip)
    a1 = 2 + alpha_n * h1
    a2 = 2 + alpha_n * h2
    b1 = 1 + beta_n * h1
    b2 = 1 + beta_n * h2

    E_scale = np.sqrt(E_H / n1**a1 * E_H / n2**a2)
    sigma_phase = R * (1.0 / n1**b1 + 1.0 / n2**b2)

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss', 'sp')) else sigma_phase * f_pi
        cont = C_bond * E_scale * abs(np.sin(ph))
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    eps1 = E_H * (Z_eff[orb1] / n1)**2
    eps2 = E_H * (Z_eff[orb2] / n2)**2
    de = abs(eps1 - eps2)
    Vc = max(abs(D_cov), 0.01)
    q = de / np.sqrt(de**2 + (2*Vc)**2) if de > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / R
    D_tot = D_cov + Di
    err = (D_tot - De_exp) / De_exp * 100
    errs_both.append(abs(err))
    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
    print(f"{name:<6} {De_exp:7.3f} {D_tot:7.3f} {err:+6.1f}% {sigma_phase/pi:6.3f} {flag}")

w5 = sum(1 for e in errs_both if e < 5)
w10 = sum(1 for e in errs_both if e < 10)
print(f"\nAvg={np.mean(errs_both):.1f}%, Med={np.median(errs_both):.1f}%, "
      f"<5%:{w5}/24, <10%:{w10}/24")
