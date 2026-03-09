"""
BF & CN Analysis Part 3: Targeted fixes
=========================================

Findings so far:
1. CN is genuinely BO=2.5 (13e radical), sigma has 1 electron
   - count=0.5 for sigma helps CN but breaks NO
   - Need to apply ONLY to CN's sigma, not NO's pi_anti

2. BF has extreme Z asymmetry (ratio 2.1)
   - Coherence suppression helps BF but hurts CO at high power

3. The "count = channels" rule works for degenerate orbitals (pi)
   but NOT for non-degenerate (sigma)

Key question: Is there a SINGLE rule from d=3 that fixes both?
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
    ('H2',   1.401,  4.745, [('ss', 1)],                                           'H_1s',  'H_1s',  14),
    ('Li2',  5.051,  1.056, [('ss', 1)],                                           'Li_2s', 'Li_2s',  6),
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


def compute_v4(mol):
    """V4 formula (baseline with 3 corrections)."""
    name, R, De_exp, bonds, orb1, orb2, ne = mol
    n1, l1 = get_n(orb1), get_l(orb1)
    n2, l2 = get_n(orb2), get_l(orb2)
    h1, h2 = has_nodes(orb1), has_nodes(orb2)

    a1 = 2 + (1 - 2*l1) * alpha * h1
    a2 = 2 + (1 - 2*l2) * alpha * h2
    b1 = 1 + beta * h1
    b2 = 1 + beta * h2

    E1 = E_H / n1**a1
    E2 = E_H / n2**a2
    E_scale = np.sqrt(E1 * E2)
    sigma_phase = R / n1**b1 + R / n2**b2

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

    return D_cov + D_ion, D_cov, D_ion


# =============================================================================
# Understand the ACTUAL PHYSICS of BF and CN
# =============================================================================
print("=" * 90)
print("  DETAILED COMPARISON: N2, CO, BF, CN")
print("=" * 90)

for name in ['N2', 'CO', 'BF', 'CN', 'NO']:
    mol = [m for m in molecules if m[0] == name][0]
    _, R, De_exp, bonds, o1, o2, ne = mol
    De_pred, D_cov, D_ion = compute_v4(mol)

    z1, z2 = Z_eff[o1], Z_eff[o2]
    z_ratio = max(z1,z2)/min(z1,z2)
    bo = sum(c if 'anti' not in bt else -c for bt, c in bonds)

    # What fraction is the prediction off by?
    cov_needed = De_exp - D_ion
    cov_excess = D_cov - cov_needed
    suppression_needed = cov_needed / D_cov if D_cov > 0 else 0

    print(f"\n  {name} (ne={ne}, BO={bo:.1f}):")
    print(f"    Z_eff: {z1:.3f}, {z2:.3f} (ratio={z_ratio:.3f})")
    print(f"    D_cov={D_cov:.3f}, D_ion={D_ion:.3f}, D_tot={De_pred:.3f}")
    print(f"    De_exp={De_exp:.3f}, err={100*(De_pred-De_exp)/De_exp:+.1f}%")
    print(f"    D_cov needed: {cov_needed:.3f} (current: {D_cov:.3f})")
    print(f"    Suppression factor needed: {suppression_needed:.3f}")


# =============================================================================
# Explore: What's physically different about BF?
# =============================================================================
print(f"\n{'='*90}")
print(f"  BF vs CO vs N2: WHY is BF weaker?")
print(f"{'='*90}")

print(f"""
  All three: isoelectronic (14e), triple bond (BO=3)

  N2: 9.76 eV  --  symmetric, maximum overlap
  CO: 11.23 eV --  mild asymmetry + strong ionic boost
  BF: 7.81 eV  --  extreme asymmetry, ionic boost insufficient

  BF's D_cov should be ~6.3 eV (= 7.81 - 1.54 ionic)
  N2's D_cov is ~9.9 eV
  Ratio: 6.3/9.9 = 0.64

  The formula gives D_cov(BF) = 8.42 instead of 6.3.
  This means the overlap is overpredicted by a factor 8.42/6.3 = 1.34.

  Physical explanation: B (Z=2.42) has much larger orbital radius than F (Z=5.10).
  The overlap integral of mismatched orbitals is suppressed.
  The E_scale (geometric mean) already partially accounts for this,
  but the SPATIAL mismatch in the overlap integral isn't captured.

  In the breather simulation: different mode numbers give REPULSIVE interaction
  (the mixed s+p case). BF maps to same mode (both 2p) but different Z_eff
  scaling, which the simulation doesn't directly test.
""")

# =============================================================================
# Can we derive the BF suppression from wave physics?
# =============================================================================
print(f"{'='*90}")
print(f"  WAVE PHYSICS DERIVATION: Spatial mismatch suppression")
print(f"{'='*90}")

print(f"""
  In wave mechanics, the overlap integral between two hydrogenic orbitals
  with different Z_eff has an analytical suppression factor.

  For 2p orbitals: <2p(Z1)|2p(Z2)> ~ [2*sqrt(Z1*Z2)/(Z1+Z2)]^5

  The exponent 5 comes from: 2*n + 1 = 2*2 + 1 = 5 for n=2 orbitals.
  (Each radial wavefunction ~ r*exp(-Zr/n), integration in 3D gives r^2 dr,
   total power = 2*(orbital factor) + 2 (from r^2) + 1 (from normalization))

  Actually the exact formula for overlap of hydrogenic orbitals is:
  S = [2*sqrt(Z1*Z2)/(Z1+Z2)]^(2n+1) for same-n orbitals

  Let's check this against needed suppression factors:
""")

for name in ['N2', 'CO', 'BF', 'CN', 'NO']:
    mol = [m for m in molecules if m[0] == name][0]
    _, R, De_exp, bonds, o1, o2, ne = mol
    De_pred, D_cov, D_ion = compute_v4(mol)

    z1, z2 = Z_eff[o1], Z_eff[o2]
    n1 = get_n(o1)
    n2_ = get_n(o2)

    # Hydrogenic overlap suppression
    if n1 == n2_:
        S_hydro = (2*np.sqrt(z1*z2)/(z1+z2))**(2*n1+1)
    else:
        S_hydro = 1.0  # different n, more complex

    cov_needed = De_exp - D_ion
    suppression_needed = cov_needed / D_cov if D_cov > 0 else 1.0

    # Also compute for 2n+1 and other exponents
    base = 2*np.sqrt(z1*z2)/(z1+z2)

    print(f"  {name}: base={base:.4f}")
    print(f"    exp=3: {base**3:.4f}, exp=4: {base**4:.4f}, exp=5: {base**5:.4f}, "
          f"exp=6: {base**6:.4f}, exp=7: {base**7:.4f}")
    print(f"    needed: {suppression_needed:.4f}")


# =============================================================================
# Test: Apply S_hydro suppression to ALL p-p bonds
# =============================================================================
print(f"\n{'='*90}")
print(f"  TEST: Hydrogenic overlap suppression for p-p bonds")
print(f"  S -> S * [2*sqrt(Z1*Z2)/(Z1+Z2)]^(2n+1)")
print(f"{'='*90}")

def compute_v5_hydro(mol, exponent=5):
    """V4 + hydrogenic overlap suppression for pp bonds."""
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

    z1, z2 = Z_eff[orb1], Z_eff[orb2]

    # Hydrogenic suppression for same-type p-p bonds
    if l1 == 1 and l2 == 1 and n1 == n2_:
        base = 2*np.sqrt(z1*z2)/(z1+z2)
        S_hydro = base ** exponent
    else:
        S_hydro = 1.0

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

        # Hydrogenic suppression for pp bonds
        if btype in ('pp_sigma', 'pi', 'pi_anti'):
            S = S * S_hydro

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


print(f"\n  Scanning exponent (2n+1 = 5 for n=2 orbitals):\n")

for exp in [0, 3, 4, 5, 6, 7]:
    errs = []
    details = {}
    for mol in molecules:
        name = mol[0]
        De_exp = mol[2]
        De_pred = compute_v5_hydro(mol, exponent=exp)
        err = abs((De_pred - De_exp) / De_exp * 100)
        errs.append(err)
        details[name] = (err, De_pred)

    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)

    key_mols = ['N2', 'CO', 'NO', 'O2', 'F2', 'Cl2', 'BF', 'CN']
    detail_str = ' '.join(f"{n}={details[n][0]:.1f}%" for n in key_mols)
    print(f"  exp={exp}: avg={np.mean(errs):.1f}%, med={np.median(errs):.1f}%, "
          f"w5={w5}/24, w10={w10}/24")
    print(f"         {detail_str}")


# =============================================================================
# What if the suppression only applies to heteronuclear pp bonds?
# =============================================================================
print(f"\n{'='*90}")
print(f"  TEST: Suppression ONLY for heteronuclear p-p (skip homonuclear)")
print(f"{'='*90}")

def compute_v5_hetero(mol, exponent=5):
    """V4 + hydrogenic overlap suppression ONLY for heteronuclear pp bonds."""
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

    z1, z2 = Z_eff[orb1], Z_eff[orb2]

    # Hydrogenic suppression for HETERONUCLEAR pp bonds only
    is_heteronuclear = (orb1 != orb2)
    if l1 == 1 and l2 == 1 and n1 == n2_ and is_heteronuclear:
        base = 2*np.sqrt(z1*z2)/(z1+z2)
        S_hydro = base ** exponent
    else:
        S_hydro = 1.0

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

        if btype in ('pp_sigma', 'pi', 'pi_anti'):
            S = S * S_hydro

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


print(f"\n  This preserves homonuclear results while fixing heteronuclear.\n")

for exp in [0, 3, 4, 5, 6, 7]:
    errs = []
    details = {}
    for mol in molecules:
        name = mol[0]
        De_exp = mol[2]
        De_pred = compute_v5_hetero(mol, exponent=exp)
        err = abs((De_pred - De_exp) / De_exp * 100)
        errs.append(err)
        details[name] = (err, De_pred)

    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)

    key_mols = ['N2', 'CO', 'NO', 'O2', 'F2', 'Cl2', 'BF', 'CN']
    detail_str = ' '.join(f"{n}={details[n][0]:.1f}%" for n in key_mols)
    print(f"  exp={exp}: avg={np.mean(errs):.1f}%, med={np.median(errs):.1f}%, "
          f"w5={w5}/24, w10={w10}/24")
    print(f"         {detail_str}")


# =============================================================================
# Full table for best result
# =============================================================================
best_exp = 5  # 2n+1 for n=2

print(f"\n{'='*90}")
print(f"  FULL TABLE: V4 + heteronuclear overlap suppression (exp={best_exp})")
print(f"{'='*90}")

print(f"\n{'Mol':<7} {'De_exp':>7} {'V4':>7} {'err_v4':>7} {'V5':>7} {'err_v5':>7}")
print("-" * 55)

errs_v4 = []
errs_v5 = []
for mol in molecules:
    name = mol[0]
    De_exp = mol[2]
    De_v4 = compute_v4(mol)[0]
    De_v5 = compute_v5_hetero(mol, exponent=best_exp)
    err_v4 = (De_v4 - De_exp) / De_exp * 100
    err_v5 = (De_v5 - De_exp) / De_exp * 100
    errs_v4.append(abs(err_v4))
    errs_v5.append(abs(err_v5))

    changed = '*' if abs(De_v5 - De_v4) > 0.001 else ''
    flag = '***' if abs(err_v5) < 2 else ' **' if abs(err_v5) < 5 else '  *' if abs(err_v5) < 10 else ''
    print(f"{name:<7} {De_exp:7.3f} {De_v4:7.3f} {err_v4:+6.1f}% {De_v5:7.3f} {err_v5:+6.1f}% {changed} {flag}")

print("-" * 55)
w5_v4 = sum(1 for e in errs_v4 if e < 5)
w10_v4 = sum(1 for e in errs_v4 if e < 10)
w5_v5 = sum(1 for e in errs_v5 if e < 5)
w10_v5 = sum(1 for e in errs_v5 if e < 10)
print(f"  V4: avg={np.mean(errs_v4):.1f}%, med={np.median(errs_v4):.1f}%, w5={w5_v4}/24, w10={w10_v4}/24")
print(f"  V5: avg={np.mean(errs_v5):.1f}%, med={np.median(errs_v5):.1f}%, w5={w5_v5}/24, w10={w10_v5}/24")

print(f"\n  Remaining outliers (>10%):")
for i, mol in enumerate(molecules):
    if errs_v5[i] > 10:
        print(f"    {mol[0]}: {errs_v5[i]:.1f}%")
