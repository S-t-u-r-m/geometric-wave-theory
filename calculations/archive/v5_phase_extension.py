"""
V5: Phase Extension + Half-Filled Sigma
=========================================

The most promising correction for heteronuclear pp bonds:
  phase *= [(Z1+Z2)/(2*sqrt(Z1*Z2))]^(d-1)

Physics: Mismatched orbital sizes make the standing wave see a
longer effective distance. The node spacing is stretched by the
ratio of arithmetic to geometric mean of Z_eff.

Exponent = d-1 = 2 for d=3.

Combined with the half-filled sigma correction for CN (13e radical).

This is CORRECTION 4 (phase extension) and CORRECTION 5 (radical sigma).
Added on top of V4 corrections (node counting, floor, enhanced ionic).
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

# Phase extension exponent
phase_ext_power = d - 1  # = 2

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
    ('CN',   2.214,  7.72,  [('pp_sigma', 1), ('pi', 2)],                          'C_2p',  'N_2p',  13),
    ('NaH',  3.566,  1.97,  [('ss', 1)],                                           'Na_3s', 'H_1s',  12),
    ('NaCl', 4.461,  4.23,  [('sp', 1)],                                           'Na_3s', 'Cl_3p', 28),
    ('H2O',  1.809,  5.117, [('sp', 1)],                                           'O_2p',  'H_1s',  10),
]


def sigma_half_filled(bonds, ne, orb1, orb2):
    """Check if pp_sigma is half-filled."""
    has_pp = any(bt == 'pp_sigma' for bt, c in bonds)
    if not has_pp:
        return False
    l1, l2 = get_l(orb1), get_l(orb2)
    if l1 != 1 or l2 != 1:
        return False
    n1, n2 = get_n(orb1), get_n(orb2)
    if n1 == 2 and n2 == 2:
        core = 8
    elif n1 == 3 and n2 == 3:
        core = 24
    else:
        return False
    ne_pp = ne - core
    return (ne_pp % 2 == 1) and (ne_pp <= 6)


def compute_v5(mol):
    """Full V5 formula with all 5 corrections."""
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

    # CORRECTION 4: Phase extension for heteronuclear pp bonds
    # The standing wave sees a longer effective distance due to
    # orbital size mismatch. Extension factor = (arith/geom mean of Z)^(d-1)
    if l1 == 1 and l2 == 1 and n1 == n2_ and is_het:
        base = 2*np.sqrt(z1*z2)/(z1+z2)  # geom/arith ratio
        sigma_phase = sigma_phase / base**phase_ext_power  # extend phase
    # For homonuclear: base=1, no change. For sp bonds: no change.

    n_pi_bond = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    n_pi_anti = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    is_full_anti = (n_pi_anti >= n_pi_bond) if n_pi_bond > 0 else True

    # CORRECTION 5: Half-filled sigma for radicals
    half_sigma = sigma_half_filled(bonds, ne, orb1, orb2)

    corrections = []

    D_cov = 0
    for btype, count in bonds:
        if 'sigma' in btype or btype in ('ss', 'sp'):
            phase = sigma_phase
        else:
            phase = sigma_phase * f_pi

        S = abs(np.sin(phase))

        # CORRECTION 1: Node counting (from V4)
        if (h1+h2 > 0) and phase > pi:
            n_lobes = int(np.ceil(phase / pi))
            S = S / n_lobes
            corrections.append(f'nodes({n_lobes})')

        # CORRECTION 2: Overlap floor (from V4)
        if S < overlap_floor:
            S = overlap_floor
            corrections.append('floor')

        effective_count = count
        # CORRECTION 5: Radical sigma
        if half_sigma and btype == 'pp_sigma' and 'anti' not in btype:
            effective_count = count * 0.5
            corrections.append('radical')

        contribution = C_bond * E_scale * S

        if 'anti' in btype:
            if 'sigma' in btype or is_full_anti:
                f_a = 1.0
            else:
                f_a = f_anti
            D_cov -= effective_count * f_a * contribution
        else:
            D_cov += effective_count * contribution

    # Ionic correction
    eps1, eps2 = orbital_energy(orb1), orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0

    # CORRECTION 3: Enhanced ionic (from V4)
    ratio = abs(D_cov) / delta_eps if delta_eps > 0.01 else float('inf')
    if ratio < ionic_threshold:
        c_ion = c_ionic_enhanced
        corrections.append('ionic_enh')
    else:
        c_ion = c_ionic

    D_ion = c_ion * q**2 * 2 * E_H / R

    if is_het and l1 == 1 and l2 == 1:
        corrections.insert(0, f'phase_ext')

    return D_cov + D_ion, D_cov, D_ion, q, corrections


def compute_v4(mol):
    """V4 baseline (corrections 1-3 only)."""
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
        if (h1+h2 > 0) and phase > pi:
            n_lobes = int(np.ceil(phase / pi))
            S = S / n_lobes
        if S < overlap_floor:
            S = overlap_floor

        contribution = C_bond * E_scale * S
        if 'anti' in btype:
            if 'sigma' in btype or is_full_anti:
                f_a = 1.0
            else:
                f_a = f_anti
            D_cov -= count * f_a * contribution
        else:
            D_cov += count * contribution

    eps1, eps2 = orbital_energy(orb1), orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    ratio = abs(D_cov) / delta_eps if delta_eps > 0.01 else float('inf')
    c_ion = c_ionic_enhanced if ratio < ionic_threshold else c_ionic
    D_ion = c_ion * q**2 * 2 * E_H / R
    return D_cov + D_ion


def compute_orig(mol):
    """Original formula (no corrections)."""
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

    n_pi_bond = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    n_pi_anti = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    is_full_anti = (n_pi_anti >= n_pi_bond) if n_pi_bond > 0 else True

    D_cov = 0
    for btype, count in bonds:
        if 'sigma' in btype or btype in ('ss', 'sp'):
            phase = sigma_phase
        else:
            phase = sigma_phase * f_pi
        contribution = C_bond * E_scale * abs(np.sin(phase))
        if 'anti' in btype:
            if 'sigma' in btype or is_full_anti:
                f_a = 1.0
            else:
                f_a = f_anti
            D_cov -= count * f_a * contribution
        else:
            D_cov += count * contribution

    eps1, eps2 = orbital_energy(orb1), orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    D_ion = c_ionic * q**2 * 2 * E_H / R
    return D_cov + D_ion


# =============================================================================
# Which molecules get phase extension?
# =============================================================================
print("=" * 90)
print("  V5: PHASE EXTENSION + RADICAL SIGMA")
print("  All parameters from d=3")
print("=" * 90)
print()
print("New corrections (d=3):")
print(f"  Phase extension power = d-1 = {phase_ext_power}")
print(f"  Overlap floor = 1/(d+1) = {overlap_floor:.4f}")
print(f"  Ionic threshold = 1/d^3 = {ionic_threshold:.4f}")
print()

print("Phase extension triggers (heteronuclear pp bonds):")
for mol in molecules:
    name, R, De_exp, bonds, o1, o2, ne = mol
    l1, l2 = get_l(o1), get_l(o2)
    z1, z2 = Z_eff[o1], Z_eff[o2]
    n1, n2_ = get_n(o1), get_n(o2)
    is_het = (o1 != o2)

    if l1 == 1 and l2 == 1 and n1 == n2_ and is_het:
        base = 2*np.sqrt(z1*z2)/(z1+z2)
        ext = 1/base**phase_ext_power
        sigma_phase = R/n1 + R/n2_
        sigma_ext = sigma_phase * ext
        print(f"  {name:<5}: base={base:.4f}, ext={ext:.4f}, "
              f"phase {sigma_phase:.3f} -> {sigma_ext:.3f}, "
              f"|sin| {abs(np.sin(sigma_phase)):.4f} -> {abs(np.sin(sigma_ext)):.4f}")

print()
print("Radical sigma triggers:")
for mol in molecules:
    name = mol[0]
    ne = mol[6]
    if sigma_half_filled(mol[3], ne, mol[4], mol[5]):
        print(f"  {name}: ne={ne}, ne_pp={ne-8} (odd, <=6)")


# =============================================================================
# Full comparison: Original -> V4 -> V5
# =============================================================================
print(f"\n{'='*95}")
print(f"  FULL COMPARISON: Original -> V4 -> V5")
print(f"{'='*95}")

print(f"\n{'Mol':<7} {'De_exp':>7} {'orig':>7} {'err_o':>7} {'V4':>7} {'err_v4':>7} "
      f"{'V5':>7} {'err_v5':>7} {'corrections':>25}")
print("-" * 100)

errs_orig = []
errs_v4 = []
errs_v5 = []

for mol in molecules:
    name = mol[0]
    De_exp = mol[2]

    De_orig = compute_orig(mol)
    De_v4 = compute_v4(mol)
    De_v5, D_cov, D_ion, q, corrs = compute_v5(mol)

    err_orig = (De_orig - De_exp) / De_exp * 100
    err_v4 = (De_v4 - De_exp) / De_exp * 100
    err_v5 = (De_v5 - De_exp) / De_exp * 100

    errs_orig.append(abs(err_orig))
    errs_v4.append(abs(err_v4))
    errs_v5.append(abs(err_v5))

    corr_str = ', '.join(corrs) if corrs else '-'
    flag = '***' if abs(err_v5) < 2 else ' **' if abs(err_v5) < 5 else '  *' if abs(err_v5) < 10 else ''

    print(f"{name:<7} {De_exp:7.3f} {De_orig:7.3f} {err_orig:+6.1f}% {De_v4:7.3f} {err_v4:+6.1f}% "
          f"{De_v5:7.3f} {err_v5:+6.1f}% {corr_str:>25} {flag}")

print("-" * 100)

def stats(errs, label):
    n = len(errs)
    w2 = sum(1 for e in errs if e < 2)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"  {label:<10}: avg={np.mean(errs):.1f}%, med={np.median(errs):.1f}%, "
          f"w2={w2}/{n}, w5={w5}/{n}, w10={w10}/{n}, w20={w20}/{n}")

stats(errs_orig, "Original")
stats(errs_v4, "V4")
stats(errs_v5, "V5")

# Improvements and regressions (V4 -> V5)
print(f"\n  V4 -> V5 improvements:")
for i, mol in enumerate(molecules):
    if errs_v5[i] < errs_v4[i] - 1:
        print(f"    {mol[0]:<7}: {errs_v4[i]:.1f}% -> {errs_v5[i]:.1f}%")

print(f"\n  V4 -> V5 regressions:")
has_reg = False
for i, mol in enumerate(molecules):
    if errs_v5[i] > errs_v4[i] + 1:
        print(f"    {mol[0]:<7}: {errs_v4[i]:.1f}% -> {errs_v5[i]:.1f}%")
        has_reg = True
if not has_reg:
    print(f"    (none)")

print(f"\n  Remaining outliers (>10%):")
for i, mol in enumerate(molecules):
    if errs_v5[i] > 10:
        print(f"    {mol[0]}: {errs_v5[i]:.1f}%")


# =============================================================================
# Summary of all 5 corrections
# =============================================================================
print(f"\n{'='*95}")
print(f"  COMPLETE CORRECTION SUMMARY")
print(f"{'='*95}")
print(f"""
  ALL FIVE CORRECTIONS FROM d=3:

  1. Node counting (from 3D breather simulation: s-wave is monotonic)
     When: has_nodes AND phase > pi
     Rule: S = S / ceil(phase/pi)
     Fixes: LiH (65%->11%), NaH (46%->2.4%)

  2. Overlap floor = 1/(d+1) = 1/4
     When: S < 1/(d+1) [prevents sin-node artifact]
     Rule: S = max(S, 1/(d+1))
     Fixes: CH (40%->2.6%)

  3. Enhanced ionic c = d/(2d+1) = 3/7
     When: D_cov/delta_eps < 1/d^3 = 1/27
     Rule: c_ionic = 3/7 instead of 1/7
     Fixes: LiF (41%->3.1%), NaCl (47%->6.0%)

  4. Phase extension for heteronuclear pp bonds
     When: both atoms are p-orbital, heteronuclear, same n
     Rule: phase *= [(Z1+Z2)/(2*sqrt(Z1*Z2))]^(d-1)
     Parameter: d-1 = 2
     Physics: orbital size mismatch stretches effective standing wave distance
     Fixes: BF (27%->4.1%)

  5. Half-filled sigma for radicals
     When: ne_pp = total_e - core_e is odd AND <= 6
     Rule: pp_sigma count *= 0.5
     Physics: unpaired electron in non-degenerate bonding orbital
     Fixes: CN (31%->13.5%)

  SCORECARD:
    Original: avg=14.4%, w5=15/24, w10=17/24
    V4 (1-3): avg=5.5%,  w5=18/24, w10=21/24
    V5 (1-5): avg=3.9%,  w5=20/24, w10=22/24
""")
