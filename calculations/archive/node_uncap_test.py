"""
Test: Uncap has_nodes so Na_3s gets h=2 (two radial nodes).
Combined with radical pi weakening and asymmetric node exponent x=1/d.

Currently has_nodes = min(n-l-1, 1) -- caps at 1.
New: has_nodes = n-l-1 (actual radial node count).

This affects:
  Li_2s: n=2, l=0 -> h=1 (unchanged)
  Na_3s: n=3, l=0 -> h=2 (was 1, now 2)
  All p-orbitals: n=2, l=1 -> h=0 (unchanged)
  Cl_3p: n=3, l=1 -> h=1 (unchanged)
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
overlap_floor = 1.0 / (d+1)
ionic_threshold = 1.0 / d**3
c_ionic_enhanced = d / (2*d + 1)
phase_ext_power = d - 1

Z_eff = {
    'H_1s': 1.0000, 'Li_2s': 1.2792, 'B_2p': 2.4214,
    'C_2p': 3.1358, 'N_2p': 3.8340, 'O_2p': 4.4532,
    'F_2p': 5.0998, 'Na_3s': 2.5074, 'Cl_3p': 6.1161,
}

def get_n(orb): return int(orb.split('_')[1][0])
def get_l(orb): return {'s': 0, 'p': 1}[orb.split('_')[1][1]]

def has_nodes_capped(orb):
    n, l = get_n(orb), get_l(orb)
    return min(n - l - 1, 1)

def has_nodes_uncapped(orb):
    n, l = get_n(orb), get_l(orb)
    return n - l - 1

def orbital_energy(orb):
    n, z = get_n(orb), Z_eff[orb]
    return E_H * (z / n)**2

molecules = [
    ('H2',   1.401,  4.745, [('ss', 1)],                                           'H_1s',  'H_1s',   2),
    ('Li2',  5.051,  1.056, [('ss', 1)],                                           'Li_2s', 'Li_2s',   6),
    ('B2',   3.005,  3.02,  [('pi', 2)],                                           'B_2p',  'B_2p',  10),
    ('C2',   2.348,  6.32,  [('pi', 2), ('sp_sigma', 1), ('sp_sigma_anti', 1)],    'C_2p',  'C_2p',  12),
    ('N2',   2.074,  9.759, [('pp_sigma', 1), ('pi', 2)],                          'N_2p',  'N_2p',  14),
    ('O2',   2.282,  5.213, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)],          'O_2p',  'O_2p',  16),
    ('F2',   2.668,  1.660, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)],          'F_2p',  'F_2p',  18),
    ('Na2',  5.818,  0.746, [('ss', 1)],                                           'Na_3s', 'Na_3s', 22),
    ('Cl2',  3.757,  2.514, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)],          'Cl_3p', 'Cl_3p', 34),
    ('HF',   1.733,  5.869, [('sp', 1)],                                           'H_1s',  'F_2p',  10),
    ('CO',   2.132, 11.225, [('pp_sigma', 1), ('pi', 2)],                          'C_2p',  'O_2p',  14),
    ('NO',   2.175,  6.497, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)],          'N_2p',  'O_2p',  15),
    ('OH',   1.834,  4.392, [('sp', 1)],                                           'O_2p',  'H_1s',   9),
    ('HCl',  2.409,  4.434, [('sp', 1)],                                           'H_1s',  'Cl_3p', 18),
    ('LiH',  3.015,  2.515, [('ss', 1)],                                           'Li_2s', 'H_1s',   4),
    ('LiF',  2.955,  5.939, [('sp', 1)],                                           'Li_2s', 'F_2p',  12),
    ('BH',   2.329,  3.42,  [('sp', 1)],                                           'B_2p',  'H_1s',   6),
    ('CH',   2.116,  3.47,  [('sp', 1)],                                           'C_2p',  'H_1s',   7),
    ('NH',   1.958,  3.57,  [('sp', 1)],                                           'N_2p',  'H_1s',   8),
    ('BF',   2.386,  7.81,  [('pp_sigma', 1), ('pi', 2)],                          'B_2p',  'F_2p',  14),
    ('CN',   2.214,  7.738, [('pp_sigma', 1), ('pi', 2)],                          'C_2p',  'N_2p',  13),
    ('NaH',  3.566,  1.97,  [('ss', 1)],                                           'Na_3s', 'H_1s',  12),
    ('NaCl', 4.461,  4.23,  [('sp', 1)],                                           'Na_3s', 'Cl_3p', 28),
    ('H2O',  1.809,  5.117, [('sp', 1)],                                           'O_2p',  'H_1s',  10),
]

def sigma_half_filled(bonds, ne, orb1, orb2):
    has_pp = any(bt == 'pp_sigma' for bt, c in bonds)
    if not has_pp: return False
    l1, l2 = get_l(orb1), get_l(orb2)
    if l1 != 1 or l2 != 1: return False
    n1, n2 = get_n(orb1), get_n(orb2)
    if n1 == 2 and n2 == 2: core = 8
    elif n1 == 3 and n2 == 3: core = 24
    else: return False
    ne_pp = ne - core
    return (ne_pp % 2 == 1) and (ne_pp <= 6)


def compute(mol, uncap_nodes=False, asymm_exp=0.0, radical_pi=False):
    name, R, De_exp, bonds, orb1, orb2, ne = mol
    n1, l1 = get_n(orb1), get_l(orb1)
    n2_, l2 = get_n(orb2), get_l(orb2)

    if uncap_nodes:
        h1, h2 = has_nodes_uncapped(orb1), has_nodes_uncapped(orb2)
    else:
        h1, h2 = has_nodes_capped(orb1), has_nodes_capped(orb2)

    a1 = 2 + (1 - 2*l1) * alpha * min(h1, 1)  # energy correction still caps at 1
    a2 = 2 + (1 - 2*l2) * alpha * min(h2, 1)
    b1 = 1 + beta * min(h1, 1)  # phase correction still caps at 1
    b2 = 1 + beta * min(h2, 1)

    E1 = E_H / n1**a1
    E2 = E_H / n2_**a2
    E_scale = np.sqrt(E1 * E2)
    sigma_phase = R / n1**b1 + R / n2_**b2

    z1, z2 = Z_eff[orb1], Z_eff[orb2]
    is_het = (orb1 != orb2)

    if l1 == 1 and l2 == 1 and n1 == n2_ and is_het:
        base = 2*np.sqrt(z1*z2)/(z1+z2)
        sigma_phase = sigma_phase / base**phase_ext_power

    half_sigma = sigma_half_filled(bonds, ne, orb1, orb2)

    n_pi_bond = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    n_pi_anti = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    is_full_anti = (n_pi_anti >= n_pi_bond) if n_pi_bond > 0 else True

    D_cov = 0
    for btype, count in bonds:
        if 'sigma' in btype or btype in ('ss', 'sp'):
            phase = sigma_phase
        else:
            phase = sigma_phase * f_pi

        S = abs(np.sin(phase))

        # Node correction with uncapped h values for the exponent
        if (h1 + h2 > 0) and phase > pi:
            n_lobes = int(np.ceil(phase / pi))
            # Asymmetric node penalty: extra exponent when h1 != h2
            if asymm_exp > 0 and h1 != h2:
                node_diff = abs(h1 - h2)  # how asymmetric
                exp = 1.0 + node_diff * asymm_exp
                S = S / (n_lobes ** exp)
            else:
                S = S / n_lobes

        if S < overlap_floor:
            S = overlap_floor

        effective_count = count
        if half_sigma and btype == 'pp_sigma' and 'anti' not in btype:
            effective_count = count * 0.5

        # Radical pi weakening
        if radical_pi and half_sigma and 'pi' in btype and 'anti' not in btype:
            ne_pp = ne - 8
            effective_count = count * (ne_pp - 1) / ne_pp

        contribution = C_bond * E_scale * S

        if 'anti' in btype:
            if 'sigma' in btype or is_full_anti:
                f_a = 1.0
            else:
                f_a = f_anti
            D_cov -= effective_count * f_a * contribution
        else:
            D_cov += effective_count * contribution

    eps1, eps2 = orbital_energy(orb1), orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    ratio_ic = abs(D_cov) / delta_eps if delta_eps > 0.01 else float('inf')
    c_ion = c_ionic_enhanced if ratio_ic < ionic_threshold else c_ionic
    D_ion = c_ion * q**2 * 2 * E_H / R

    return D_cov + D_ion


# =====================================================================
# Show what uncapping changes
# =====================================================================
print("=" * 80)
print("  NODE COUNT COMPARISON: capped vs uncapped")
print("=" * 80)

for orb in ['H_1s', 'Li_2s', 'B_2p', 'C_2p', 'N_2p', 'O_2p', 'F_2p', 'Na_3s', 'Cl_3p']:
    hc = has_nodes_capped(orb)
    hu = has_nodes_uncapped(orb)
    flag = " <-- CHANGED" if hc != hu else ""
    print(f"  {orb:<8}: capped={hc}, uncapped={hu}{flag}")


# =====================================================================
# Which molecules are affected by uncapping?
# =====================================================================
print(f"\n{'='*80}")
print("  AFFECTED MOLECULES (those using Na_3s)")
print("=" * 80)

for mol in molecules:
    name, R, De_exp, bonds, orb1, orb2, ne = mol
    h1c, h2c = has_nodes_capped(orb1), has_nodes_capped(orb2)
    h1u, h2u = has_nodes_uncapped(orb1), has_nodes_uncapped(orb2)
    if h1c != h1u or h2c != h2u:
        asym_c = abs(h1c - h2c)
        asym_u = abs(h1u - h2u)
        print(f"  {name:<7}: ({orb1}, {orb2})")
        print(f"    capped:   h=({h1c},{h2c}), |diff|={asym_c}")
        print(f"    uncapped: h=({h1u},{h2u}), |diff|={asym_u}")


# =====================================================================
# MAIN TEST: scan asymm_exp with x = 1/d
# =====================================================================
print(f"\n{'='*80}")
print(f"  SCAN: uncapped nodes + radical_pi + asymm_exp = 1/d = {1/d:.4f}")
print(f"{'='*80}")

x = 1.0/d  # = 1/3

configs = [
    ("V5 current (capped, no radical_pi)", False, 0.0, False),
    ("+ radical_pi only", False, 0.0, True),
    ("+ asymm x=1/d (capped)", False, x, True),
    ("+ asymm x=1/d (UNCAPPED)", True, x, True),
]

for label, uncap, asymm, rad_pi in configs:
    errs = []
    details = {}
    for mol in molecules:
        De_pred = compute(mol, uncap_nodes=uncap, asymm_exp=asymm, radical_pi=rad_pi)
        err = (De_pred - mol[2]) / mol[2] * 100
        errs.append(abs(err))
        details[mol[0]] = err

    avg = np.mean(errs)
    w2 = sum(1 for e in errs if e < 2)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    print(f"\n  {label}:")
    print(f"    avg={avg:.2f}% w2={w2}/24 w5={w5}/24 w10={w10}/24")
    print(f"    CN={details['CN']:+.1f}% LiH={details['LiH']:+.1f}% "
          f"NaH={details['NaH']:+.1f}% Na2={details['Na2']:+.1f}% "
          f"Li2={details['Li2']:+.1f}% NaCl={details['NaCl']:+.1f}% "
          f"LiF={details['LiF']:+.1f}%")


# =====================================================================
# DETAILED: uncapped + radical_pi + x=1/d — full molecule list
# =====================================================================
print(f"\n{'='*80}")
print(f"  FULL RESULTS: uncapped nodes + radical_pi + asymm x=1/d")
print(f"{'='*80}")

print(f"\n  {'Mol':<7} {'De_exp':>7} {'De_pred':>8} {'err':>7} {'h1':>3} {'h2':>3} {'|dh|':>4}")
print(f"  {'-'*45}")

all_errs = []
for mol in molecules:
    name, R, De_exp, bonds, orb1, orb2, ne = mol
    De_pred = compute(mol, uncap_nodes=True, asymm_exp=1.0/d, radical_pi=True)
    err = (De_pred - De_exp) / De_exp * 100
    all_errs.append(abs(err))
    h1 = has_nodes_uncapped(orb1)
    h2 = has_nodes_uncapped(orb2)
    dh = abs(h1 - h2)
    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
    print(f"  {name:<7} {De_exp:7.3f} {De_pred:8.3f} {err:+6.1f}% {h1:3d} {h2:3d} {dh:4d} {flag}")

print(f"  {'-'*45}")
avg = np.mean(all_errs)
w2 = sum(1 for e in all_errs if e < 2)
w5 = sum(1 for e in all_errs if e < 5)
w10 = sum(1 for e in all_errs if e < 10)
w20 = sum(1 for e in all_errs if e < 20)
print(f"  avg={avg:.2f}% med={np.median(all_errs):.2f}% "
      f"w2={w2}/24 w5={w5}/24 w10={w10}/24 w20={w20}/24")


# =====================================================================
# ALSO TEST: what if the phase/energy corrections (alpha, beta)
# also use uncapped h? (more aggressive change)
# =====================================================================
print(f"\n{'='*80}")
print(f"  BONUS: What if alpha/beta also use uncapped h?")
print(f"  (uses h directly in energy and phase exponents)")
print(f"{'='*80}")

def compute_full_uncap(mol, asymm_exp=0.0, radical_pi=False):
    """Fully uncapped: alpha, beta, AND node correction all use true h."""
    name, R, De_exp, bonds, orb1, orb2, ne = mol
    n1, l1 = get_n(orb1), get_l(orb1)
    n2_, l2 = get_n(orb2), get_l(orb2)
    h1 = has_nodes_uncapped(orb1)
    h2 = has_nodes_uncapped(orb2)

    # Use actual h values (not capped) in alpha/beta
    a1 = 2 + (1 - 2*l1) * alpha * h1
    a2 = 2 + (1 - 2*l2) * alpha * h2
    b1 = 1 + beta * h1
    b2 = 1 + beta * h2

    E1 = E_H / n1**a1
    E2 = E_H / n2_**a2
    E_scale = np.sqrt(E1 * E2)
    sigma_phase = R / n1**b1 + R / n2_**b2

    z1, z2 = Z_eff[orb1], Z_eff[orb2]
    is_het = (orb1 != orb2)

    if l1 == 1 and l2 == 1 and n1 == n2_ and is_het:
        base = 2*np.sqrt(z1*z2)/(z1+z2)
        sigma_phase = sigma_phase / base**phase_ext_power

    half_sigma = sigma_half_filled(bonds, ne, orb1, orb2)

    n_pi_bond = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    n_pi_anti = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    is_full_anti = (n_pi_anti >= n_pi_bond) if n_pi_bond > 0 else True

    D_cov = 0
    for btype, count in bonds:
        if 'sigma' in btype or btype in ('ss', 'sp'):
            phase = sigma_phase
        else:
            phase = sigma_phase * f_pi

        S = abs(np.sin(phase))

        if (h1 + h2 > 0) and phase > pi:
            n_lobes = int(np.ceil(phase / pi))
            if asymm_exp > 0 and h1 != h2:
                node_diff = abs(h1 - h2)
                exp = 1.0 + node_diff * asymm_exp
                S = S / (n_lobes ** exp)
            else:
                S = S / n_lobes

        if S < overlap_floor:
            S = overlap_floor

        effective_count = count
        if half_sigma and btype == 'pp_sigma' and 'anti' not in btype:
            effective_count = count * 0.5

        if radical_pi and half_sigma and 'pi' in btype and 'anti' not in btype:
            ne_pp = ne - 8
            effective_count = count * (ne_pp - 1) / ne_pp

        contribution = C_bond * E_scale * S

        if 'anti' in btype:
            if 'sigma' in btype or is_full_anti:
                f_a = 1.0
            else:
                f_a = f_anti
            D_cov -= effective_count * f_a * contribution
        else:
            D_cov += effective_count * contribution

    eps1, eps2 = orbital_energy(orb1), orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    ratio_ic = abs(D_cov) / delta_eps if delta_eps > 0.01 else float('inf')
    c_ion = c_ionic_enhanced if ratio_ic < ionic_threshold else c_ionic
    D_ion = c_ion * q**2 * 2 * E_H / R

    return D_cov + D_ion

print(f"\n  {'Mol':<7} {'De_exp':>7} {'De_pred':>8} {'err':>7}")
print(f"  {'-'*30}")

all_errs2 = []
for mol in molecules:
    De_pred = compute_full_uncap(mol, asymm_exp=1.0/d, radical_pi=True)
    err = (De_pred - mol[2]) / mol[2] * 100
    all_errs2.append(abs(err))
    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
    print(f"  {mol[0]:<7} {mol[2]:7.3f} {De_pred:8.3f} {err:+6.1f}% {flag}")

print(f"  {'-'*30}")
avg2 = np.mean(all_errs2)
w2_2 = sum(1 for e in all_errs2 if e < 2)
w5_2 = sum(1 for e in all_errs2 if e < 5)
w10_2 = sum(1 for e in all_errs2 if e < 10)
print(f"  avg={avg2:.2f}% w2={w2_2}/24 w5={w5_2}/24 w10={w10_2}/24")
