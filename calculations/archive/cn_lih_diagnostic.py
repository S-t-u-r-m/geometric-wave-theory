"""
Diagnostic: decompose CN and LiH predictions to find the overshoot source.
Also test targeted fixes for both outliers.
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
def has_nodes(orb):
    n, l = get_n(orb), get_l(orb)
    return min(n - l - 1, 1)
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
    if not has_pp:
        return False
    l1, l2 = get_l(orb1), get_l(orb2)
    if l1 != 1 or l2 != 1:
        return False
    n1, n2 = get_n(orb1), get_n(orb2)
    if n1 == 2 and n2 == 2: core = 8
    elif n1 == 3 and n2 == 3: core = 24
    else: return False
    ne_pp = ne - core
    return (ne_pp % 2 == 1) and (ne_pp <= 6)


# =====================================================================
# SECTION 1: DECOMPOSE TRIPLE BONDS
# =====================================================================
print("=" * 90)
print("  TRIPLE BOND DECOMPOSITION (V5)")
print("=" * 90)

triple_bonds = ['N2', 'CO', 'BF', 'CN', 'NO']

for mol in molecules:
    name = mol[0]
    if name not in triple_bonds:
        continue

    name, R, De_exp, bonds, orb1, orb2, ne = mol
    n1, l1 = get_n(orb1), get_l(orb1)
    n2_, l2 = get_n(orb2), get_l(orb2)
    h1, h2 = has_nodes(orb1), has_nodes(orb2)

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

    phase_ext = 1.0
    if l1 == 1 and l2 == 1 and n1 == n2_ and is_het:
        base = 2*np.sqrt(z1*z2)/(z1+z2)
        phase_ext = 1/base**phase_ext_power
        sigma_phase *= phase_ext

    half_sigma = sigma_half_filled(bonds, ne, orb1, orb2)

    pi_phase = sigma_phase * f_pi
    S_sigma = abs(np.sin(sigma_phase))
    S_pi = abs(np.sin(pi_phase))

    # Decompose
    D_sigma = 0
    D_pi = 0
    D_pi_anti = 0
    for btype, count in bonds:
        if 'sigma' in btype or btype in ('ss', 'sp'):
            S = S_sigma
            eff_count = count * 0.5 if (half_sigma and btype == 'pp_sigma') else count
            contrib = eff_count * C_bond * E_scale * S
            if 'anti' in btype:
                D_sigma -= contrib
            else:
                D_sigma += contrib
        elif 'pi' in btype:
            S = S_pi
            contrib = count * C_bond * E_scale * S
            if 'anti' in btype:
                n_pi_bond = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
                n_pi_anti = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
                is_full_anti = (n_pi_anti >= n_pi_bond) if n_pi_bond > 0 else True
                f_a = 1.0 if is_full_anti else f_anti
                D_pi_anti -= count * f_a * C_bond * E_scale * S
            else:
                D_pi += contrib

    D_cov = D_sigma + D_pi + D_pi_anti

    eps1, eps2 = orbital_energy(orb1), orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    ratio_ic = abs(D_cov) / delta_eps if delta_eps > 0.01 else float('inf')
    c_ion = c_ionic_enhanced if ratio_ic < ionic_threshold else c_ionic
    D_ion = c_ion * q**2 * 2 * E_H / R

    De_pred = D_cov + D_ion
    err = (De_pred - De_exp) / De_exp * 100

    print(f"\n  {name} (R={R}, ne={ne}, {'radical' if half_sigma else 'closed-shell'}):")
    print(f"    E_scale = {E_scale:.4f} eV")
    print(f"    sigma_phase = {sigma_phase:.4f} (phase_ext={phase_ext:.4f})")
    print(f"    pi_phase    = {pi_phase:.4f}")
    print(f"    |sin(sigma)| = {S_sigma:.4f}")
    print(f"    |sin(pi)|    = {S_pi:.4f}")
    print(f"    D_sigma = {D_sigma:+.4f} eV  (count={0.5 if half_sigma else 1})")
    print(f"    D_pi    = {D_pi:+.4f} eV  (count=2)")
    print(f"    D_anti  = {D_pi_anti:+.4f} eV")
    print(f"    D_cov   = {D_cov:+.4f} eV")
    print(f"    D_ion   = {D_ion:+.4f} eV  (q={q:.4f})")
    print(f"    De_pred = {De_pred:.4f} eV  vs  De_exp = {De_exp:.4f} eV")
    print(f"    Error   = {err:+.1f}%")
    print(f"    Overshoot = {De_pred - De_exp:+.4f} eV")


# =====================================================================
# SECTION 2: DECOMPOSE LiH
# =====================================================================
print(f"\n{'='*90}")
print("  LiH DECOMPOSITION")
print("=" * 90)

for mol in molecules:
    if mol[0] != 'LiH':
        continue

    name, R, De_exp, bonds, orb1, orb2, ne = mol
    n1, l1 = get_n(orb1), get_l(orb1)
    n2_, l2 = get_n(orb2), get_l(orb2)
    h1, h2 = has_nodes(orb1), has_nodes(orb2)

    a1 = 2 + (1 - 2*l1) * alpha * h1
    a2 = 2 + (1 - 2*l2) * alpha * h2
    b1 = 1 + beta * h1
    b2 = 1 + beta * h2

    E1 = E_H / n1**a1
    E2 = E_H / n2_**a2
    E_scale = np.sqrt(E1 * E2)
    sigma_phase = R / n1**b1 + R / n2_**b2

    n_lobes = int(np.ceil(sigma_phase / pi)) if sigma_phase > pi else 1
    S = abs(np.sin(sigma_phase))
    S_raw = S
    if h1 + h2 > 0 and sigma_phase > pi:
        S = S / n_lobes
    S_floor = max(S, overlap_floor)

    D_cov = C_bond * E_scale * S_floor

    eps1, eps2 = orbital_energy(orb1), orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    ratio_ic = abs(D_cov) / delta_eps if delta_eps > 0.01 else float('inf')
    c_ion = c_ionic_enhanced if ratio_ic < ionic_threshold else c_ionic
    D_ion = c_ion * q**2 * 2 * E_H / R

    De_pred = D_cov + D_ion
    err = (De_pred - De_exp) / De_exp * 100

    print(f"\n  LiH (R={R}, {orb1} + {orb2}):")
    print(f"    n1={n1}, n2={n2_}, h1={h1}, h2={h2}")
    print(f"    b1={b1:.4f}, b2={b2:.4f}")
    print(f"    E1={E1:.4f}, E2={E2:.4f}, E_scale={E_scale:.4f}")
    print(f"    sigma_phase = R/n1^b1 + R/n2^b2 = {R}/{n1}^{b1:.3f} + {R}/{n2_}^{b2:.3f}")
    print(f"                = {R/n1**b1:.4f} + {R/n2_**b2:.4f} = {sigma_phase:.4f}")
    print(f"    phase/pi = {sigma_phase/pi:.4f}")
    print(f"    |sin(phase)| = {S_raw:.4f}")
    print(f"    n_lobes = {n_lobes}")
    print(f"    S after node correction = {S:.4f}")
    print(f"    S after floor = {S_floor:.4f}")
    print(f"    D_cov = {D_cov:.4f} eV")
    print(f"    D_ion = {D_ion:.4f} eV  (q={q:.4f}, c_ion={c_ion:.4f})")
    print(f"    De_pred = {De_pred:.4f} eV  vs  De_exp = {De_exp:.4f}")
    print(f"    Error = {err:+.1f}%")
    print(f"    Overshoot = {De_pred - De_exp:+.4f} eV")


# =====================================================================
# SECTION 3: TEST FIX IDEAS
# =====================================================================
print(f"\n{'='*90}")
print("  TEST FIX IDEAS")
print("=" * 90)

# --- FIX A: For radicals, scale ALL pp bonds by ne_pp/(ne_pp+1) ---
# Only CN triggers (ne_pp=5, odd, <=6)
# Factor = 5/6 = 0.833

print(f"\n  FIX A: Radical weakening of ALL bonds: D_cov *= ne_pp/(ne_pp+1)")
print(f"  Only triggers when ne_pp odd AND <=6 (same as half-sigma)")

# --- FIX B: For radicals, pi count *= (ne_pp-1)/ne_pp ---
# CN: ne_pp=5, pi_scale = 4/5 = 0.8

print(f"\n  FIX B: Radical pi weakening: pi count *= (ne_pp-1)/ne_pp")

# --- FIX C: Continuous node correction: S /= (phase/pi) instead of ceil ---
# Affects LiH (and NaH, Li2, Na2)

print(f"\n  FIX C: Continuous node correction: S /= max(1, phase/pi)")

# --- FIX D: Asymmetric node penalty: extra /sqrt(2) when h1 != h2 ---

print(f"\n  FIX D: Asymmetric node penalty: S /= sqrt(2) extra when h1 != h2")

# --- FIX E: For LiH, the node correction exponent: S /= n_lobes^(1 + |h1-h2|*x) ---

print(f"\n  FIX E: Node exponent: S /= n_lobes^(1+x) when h1 != h2")


def compute_with_fix(mol, fix='none', fix_param=1.0):
    """Compute De with a specific fix applied."""
    name, R, De_exp, bonds, orb1, orb2, ne = mol
    n1, l1 = get_n(orb1), get_l(orb1)
    n2_, l2 = get_n(orb2), get_l(orb2)
    h1, h2 = has_nodes(orb1), has_nodes(orb2)

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

    # Radical weakening factor for pi bonds
    radical_pi_scale = 1.0
    if fix == 'radical_all' and half_sigma:
        ne_pp = ne - 8
        radical_pi_scale = ne_pp / (ne_pp + 1)
    elif fix == 'radical_pi' and half_sigma:
        ne_pp = ne - 8
        radical_pi_scale = (ne_pp - 1) / ne_pp

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

        # Node correction
        if fix == 'continuous_node' and (h1+h2 > 0) and phase > pi:
            S = S / max(1.0, phase / pi)
        elif fix == 'asymmetric_node' and (h1+h2 > 0) and phase > pi:
            n_lobes = int(np.ceil(phase / pi))
            if h1 != h2:
                S = S / (n_lobes ** fix_param)
            else:
                S = S / n_lobes
        elif fix == 'node_exp' and (h1+h2 > 0) and phase > pi:
            n_lobes = int(np.ceil(phase / pi))
            exp = 1.0 + abs(h1-h2) * fix_param
            S = S / (n_lobes ** exp)
        elif (h1+h2 > 0) and phase > pi:
            n_lobes = int(np.ceil(phase / pi))
            S = S / n_lobes

        if S < overlap_floor:
            S = overlap_floor

        effective_count = count
        if half_sigma and btype == 'pp_sigma' and 'anti' not in btype:
            effective_count = count * 0.5

        # Apply radical pi scaling
        if 'pi' in btype and 'anti' not in btype:
            effective_count = count * radical_pi_scale

        contribution = C_bond * E_scale * S

        if 'anti' in btype:
            if 'sigma' in btype or is_full_anti:
                f_a = 1.0
            else:
                f_a = f_anti
            D_cov -= effective_count * f_a * contribution
        else:
            D_cov += effective_count * contribution

    # Apply fix A: scale entire D_cov
    if fix == 'radical_all' and half_sigma:
        ne_pp = ne - 8
        D_cov *= ne_pp / (ne_pp + 1)

    # Ionic
    eps1, eps2 = orbital_energy(orb1), orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    ratio_ic = abs(D_cov) / delta_eps if delta_eps > 0.01 else float('inf')
    c_ion = c_ionic_enhanced if ratio_ic < ionic_threshold else c_ionic
    D_ion = c_ion * q**2 * 2 * E_H / R

    return D_cov + D_ion


# Run all fixes and compare
fixes = [
    ('V5 (current)', 'none', 1.0),
    ('Fix A: radical D_cov*5/6', 'radical_all', 1.0),
    ('Fix B: radical pi*4/5', 'radical_pi', 1.0),
    ('Fix C: continuous node', 'continuous_node', 1.0),
    ('Fix D: asymm node ^1.5', 'asymmetric_node', 1.5),
    ('Fix D: asymm node ^2.0', 'asymmetric_node', 2.0),
    ('Fix E: node_exp x=0.3', 'node_exp', 0.3),
    ('Fix E: node_exp x=0.5', 'node_exp', 0.5),
    ('Fix E: node_exp x=0.7', 'node_exp', 0.7),
    ('Fix E: node_exp x=1.0', 'node_exp', 1.0),
]

print(f"\n{'='*90}")
print(f"  IMPACT OF EACH FIX ON ALL 24 MOLECULES")
print(f"{'='*90}")

for fix_label, fix_type, fix_param in fixes:
    errs = []
    cn_err = None
    lih_err = None
    for mol in molecules:
        De_pred = compute_with_fix(mol, fix_type, fix_param)
        err = abs((De_pred - mol[2]) / mol[2] * 100)
        errs.append(err)
        if mol[0] == 'CN':
            cn_err = (De_pred - mol[2]) / mol[2] * 100
        if mol[0] == 'LiH':
            lih_err = (De_pred - mol[2]) / mol[2] * 100

    avg = np.mean(errs)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)

    print(f"  {fix_label:<30}: avg={avg:.1f}% w5={w5}/24 w10={w10}/24  "
          f"CN={cn_err:+.1f}% LiH={lih_err:+.1f}%")


# =====================================================================
# SECTION 4: Best combined fix
# =====================================================================
print(f"\n{'='*90}")
print(f"  COMBINED FIXES: radical pi + node exponent scan")
print(f"{'='*90}")

# Try combining radical_pi with node_exp
print(f"\n  Combining: pi_scale for radicals + node_exp for asymmetric nodes")
print(f"  Scanning node_exp parameter x")

best_avg = 100
best_x = 0
best_combo = None

for x in np.arange(0.0, 2.01, 0.1):
    errs = []
    details = {}
    for mol in molecules:
        name, R, De_exp, bonds, orb1, orb2, ne = mol
        n1, l1 = get_n(orb1), get_l(orb1)
        n2_, l2 = get_n(orb2), get_l(orb2)
        h1, h2 = has_nodes(orb1), has_nodes(orb2)

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

            # Node correction with asymmetric exponent
            if (h1+h2 > 0) and phase > pi:
                n_lobes = int(np.ceil(phase / pi))
                exp = 1.0 + abs(h1-h2) * x
                S = S / (n_lobes ** exp)

            if S < overlap_floor:
                S = overlap_floor

            effective_count = count
            if half_sigma and btype == 'pp_sigma' and 'anti' not in btype:
                effective_count = count * 0.5

            # Radical pi weakening
            if half_sigma and 'pi' in btype and 'anti' not in btype:
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

        De_pred = D_cov + D_ion
        err_pct = (De_pred - De_exp) / De_exp * 100
        errs.append(abs(err_pct))
        details[name] = err_pct

    avg = np.mean(errs)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w2 = sum(1 for e in errs if e < 2)

    if avg < best_avg:
        best_avg = avg
        best_x = x
        best_combo = details.copy()

    if x in [0.0, 0.5, 1.0, 1.5, 2.0] or abs(avg - best_avg) < 0.05:
        print(f"  x={x:.1f}: avg={avg:.2f}% w2={w2} w5={w5} w10={w10}  "
              f"CN={details['CN']:+.1f}% LiH={details['LiH']:+.1f}% "
              f"Li2={details['Li2']:+.1f}% NaH={details['NaH']:+.1f}% Na2={details['Na2']:+.1f}%")


print(f"\n  Best: x={best_x:.1f}, avg={best_avg:.2f}%")

# Show full results for best combo
print(f"\n  Full results for best combo (x={best_x:.1f}):")
print(f"  {'Mol':<7} {'err':>7}")
for mol in molecules:
    name = mol[0]
    err = best_combo[name]
    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else '   '
    print(f"  {name:<7} {err:+7.1f}% {flag}")
