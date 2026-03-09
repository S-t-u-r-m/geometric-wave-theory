"""
BF & CN Analysis Part 4: Combined corrections
===============================================

Two independent failure modes, two independent corrections:

1. BF: Hydrogenic overlap suppression [2*sqrt(Z1*Z2)/(Z1+Z2)]^(2n+1)
   for heteronuclear p-p bonds. Exponent = 2n+1 = 5 for n=2 orbitals.
   This is the SPATIAL MISMATCH correction.

2. CN: Half-filled sigma orbital. The question is how to implement
   this without breaking NO.

The problem with the naive approach (count=0.5 for CN sigma):
- CN: sigma has 1 electron -> count=0.5 -> helps (31% -> 14.5%)
- NO: pi_anti has 1 electron -> count=0.5 -> HURTS (5.3% -> 31.9%)

But NO and CN are NOT the same situation:
- CN: bonding orbital half-filled -> less bonding
- NO: antibonding orbital half-filled -> less antibonding -> MORE bonding

In the current formula, NO uses ('pi_anti', 1) with f_anti = 6/5.
The f_anti enhancement (from partial cancellation) already captures the
partial filling. So NO doesn't need a count correction.

For CN: the sigma is a NON-DEGENERATE bonding orbital with 1 electron.
Unlike B2's pi (which has 2 degenerate channels), sigma is a single channel.
One electron in a single bonding channel = half the pair energy.

Rule: for non-degenerate bonding orbitals, count *= (electrons/2).
This affects only CN's pp_sigma in our set.

NO's pi_anti is DEGENERATE (two pi* orbitals). One electron fills one of
two channels = the count is already right at 1 (= channels_filled, not pairs).

Let's test this interpretation.
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

# Standard molecule data with total electron counts
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


def compute_v5(mol, hydro_exp=5, radical_sigma=True):
    """
    V4 + two new corrections:
    1. Hydrogenic overlap suppression for heteronuclear pp bonds
    2. Half-count for half-filled non-degenerate bonding orbitals (CN only)
    """
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
    is_het = (orb1 != orb2)

    # Correction 4: Hydrogenic overlap suppression
    if l1 == 1 and l2 == 1 and n1 == n2_ and is_het and hydro_exp > 0:
        base = 2*np.sqrt(z1*z2)/(z1+z2)
        S_hydro = base ** hydro_exp
    else:
        S_hydro = 1.0

    # Correction 5: Odd-electron sigma (half-filled non-degenerate bonding)
    # CN has 13 electrons -> odd -> sigma_2p has 1 electron
    # Rule: if total electron count is odd AND molecule has pp_sigma bonding,
    # the HOMO sigma is half-filled -> count *= 0.5
    is_odd_electron = (ne % 2 == 1)
    has_pp_sigma_bonding = any(bt == 'pp_sigma' for bt, c in bonds)

    D_cov = 0
    for btype, count in bonds:
        if 'sigma' in btype or btype in ('ss', 'sp'):
            phase = sigma_phase
        else:
            phase = sigma_phase * f_pi

        S = abs(np.sin(phase))

        # V4 correction 1: Node counting
        if (h1+h2 > 0) and phase > pi:
            n_lobes = int(np.ceil(phase / pi))
            S = S / n_lobes

        # V4 correction 2: Overlap floor
        if S < overlap_floor:
            S = overlap_floor

        # New correction 4: Hydrogenic suppression for pp bonds
        if btype in ('pp_sigma', 'pi', 'pi_anti'):
            S = S * S_hydro

        effective_count = count
        # New correction 5: Radical sigma
        if radical_sigma and is_odd_electron and btype == 'pp_sigma' and 'anti' not in btype:
            effective_count = count * 0.5

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
    ratio = abs(D_cov) / delta_eps if delta_eps > 0.01 else float('inf')
    c_ion = c_ionic_enhanced if ratio < ionic_threshold else c_ionic
    D_ion = c_ion * q**2 * 2 * E_H / R

    return D_cov + D_ion, D_cov, D_ion, q


# =============================================================================
# Check which molecules are odd-electron with pp_sigma
# =============================================================================
print("=" * 90)
print("  ODD-ELECTRON MOLECULES WITH pp_sigma")
print("=" * 90)
print()

for mol in molecules:
    name, R, De_exp, bonds, o1, o2, ne = mol
    is_odd = ne % 2 == 1
    has_pp = any(bt == 'pp_sigma' for bt, c in bonds)
    if is_odd or has_pp:
        print(f"  {name:<7} ne={ne:2d}  odd={is_odd}  pp_sigma={has_pp}  "
              f"{'-> RADICAL CORRECTION' if is_odd and has_pp else ''}")


# =============================================================================
# Test V5 with both corrections
# =============================================================================
print(f"\n{'='*90}")
print(f"  V5: V4 + hydrogenic(exp=5) + radical_sigma")
print(f"{'='*90}")

print(f"\n{'Mol':<7} {'De_exp':>7} {'V4':>7} {'err_v4':>7} {'V5':>7} {'err_v5':>7} {'corrections':>20}")
print("-" * 75)

errs_v4 = []
errs_v5 = []
for mol in molecules:
    name = mol[0]
    De_exp = mol[2]
    ne = mol[6]

    De_v4 = compute_v5(mol, hydro_exp=0, radical_sigma=False)[0]
    De_v5, D_cov, D_ion, q = compute_v5(mol, hydro_exp=5, radical_sigma=True)
    err_v4 = (De_v4 - De_exp) / De_exp * 100
    err_v5 = (De_v5 - De_exp) / De_exp * 100
    errs_v4.append(abs(err_v4))
    errs_v5.append(abs(err_v5))

    corrs = []
    z1, z2 = Z_eff[mol[4]], Z_eff[mol[5]]
    l1, l2 = get_l(mol[4]), get_l(mol[5])
    if l1 == 1 and l2 == 1 and mol[4] != mol[5]:
        base = 2*np.sqrt(z1*z2)/(z1+z2)
        corrs.append(f'hydro({base**5:.3f})')
    if ne % 2 == 1 and any(bt == 'pp_sigma' for bt, c in mol[3]):
        corrs.append('radical')
    corr_str = ', '.join(corrs) if corrs else '-'

    flag = '***' if abs(err_v5) < 2 else ' **' if abs(err_v5) < 5 else '  *' if abs(err_v5) < 10 else ''
    print(f"{name:<7} {De_exp:7.3f} {De_v4:7.3f} {err_v4:+6.1f}% {De_v5:7.3f} {err_v5:+6.1f}% {corr_str:>20} {flag}")

print("-" * 75)
w5_v4 = sum(1 for e in errs_v4 if e < 5)
w10_v4 = sum(1 for e in errs_v4 if e < 10)
w5_v5 = sum(1 for e in errs_v5 if e < 5)
w10_v5 = sum(1 for e in errs_v5 if e < 10)
w20_v5 = sum(1 for e in errs_v5 if e < 20)
print(f"  V4: avg={np.mean(errs_v4):.1f}%, med={np.median(errs_v4):.1f}%, w5={w5_v4}/24, w10={w10_v4}/24")
print(f"  V5: avg={np.mean(errs_v5):.1f}%, med={np.median(errs_v5):.1f}%, "
      f"w5={w5_v5}/24, w10={w10_v5}/24, w20={w20_v5}/24")

print(f"\n  Improvements:")
for i, mol in enumerate(molecules):
    if errs_v5[i] < errs_v4[i] - 1:
        print(f"    {mol[0]:<7}: {errs_v4[i]:.1f}% -> {errs_v5[i]:.1f}%")

print(f"\n  Regressions:")
for i, mol in enumerate(molecules):
    if errs_v5[i] > errs_v4[i] + 1:
        print(f"    {mol[0]:<7}: {errs_v4[i]:.1f}% -> {errs_v5[i]:.1f}%")

print(f"\n  Remaining outliers (>10%):")
for i, mol in enumerate(molecules):
    if errs_v5[i] > 10:
        print(f"    {mol[0]}: {errs_v5[i]:.1f}%")


# =============================================================================
# But wait: does the hydrogenic suppression also affect NO?
# Let's separate the effects
# =============================================================================
print(f"\n{'='*90}")
print(f"  ABLATION STUDY: Each correction independently")
print(f"{'='*90}")

configs = [
    ("V4 (baseline)",       0, False),
    ("+ hydro only (exp=5)", 5, False),
    ("+ radical only",       0, True),
    ("+ both",               5, True),
]

for label, h_exp, rad in configs:
    errs = []
    details = {}
    for mol in molecules:
        De_pred = compute_v5(mol, hydro_exp=h_exp, radical_sigma=rad)[0]
        err = abs((De_pred - mol[2]) / mol[2] * 100)
        errs.append(err)
        details[mol[0]] = err

    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    print(f"  {label:<25}: avg={np.mean(errs):.1f}%, med={np.median(errs):.1f}%, "
          f"w5={w5}/24, w10={w10}/24 | "
          f"BF={details['BF']:.1f}% CN={details['CN']:.1f}% CO={details['CO']:.1f}% NO={details['NO']:.1f}%")


# =============================================================================
# The CO regression: hydro suppresses CO by 7.2% when it was at -1.2%
# Can we reduce the exponent to balance?
# =============================================================================
print(f"\n{'='*90}")
print(f"  EXPONENT SCAN with radical correction ON")
print(f"{'='*90}")

for exp in range(0, 8):
    errs = []
    details = {}
    for mol in molecules:
        De_pred = compute_v5(mol, hydro_exp=exp, radical_sigma=True)[0]
        err = abs((De_pred - mol[2]) / mol[2] * 100)
        errs.append(err)
        details[mol[0]] = err

    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    print(f"  exp={exp}: avg={np.mean(errs):.1f}%, med={np.median(errs):.1f}%, "
          f"w5={w5}/24, w10={w10}/24 | "
          f"BF={details['BF']:.1f}% CN={details['CN']:.1f}% CO={details['CO']:.1f}% "
          f"NO={details['NO']:.1f}% LiH={details['LiH']:.1f}%")


# =============================================================================
# Summary of the CN problem
# =============================================================================
print(f"\n{'='*90}")
print(f"  SUMMARY: CN RESIDUAL ERROR")
print(f"{'='*90}")

# With radical sigma + hydro(5):
De_cn_v5, D_cov_cn, D_ion_cn, q_cn = compute_v5(
    ('CN', 2.214, 7.72, [('pp_sigma', 1), ('pi', 2)], 'C_2p', 'N_2p', 13),
    hydro_exp=5, radical_sigma=True)

De_cn_exp = 7.72

print(f"""
  CN with all corrections:
    D_cov  = {D_cov_cn:.3f}  (sigma at half-count + hydro suppression)
    D_ion  = {D_ion_cn:.3f}
    D_tot  = {De_cn_v5:.3f}  vs exp {De_cn_exp:.3f}
    error  = {100*(De_cn_v5-De_cn_exp)/De_cn_exp:+.1f}%

  The radical correction reduces sigma by half: 2.85 -> 1.43
  The hydro suppression is negligible: base^5 = 0.975

  Remaining overprediction of ~{100*(De_cn_v5-De_cn_exp)/De_cn_exp:.0f}% comes from
  the pi bonds being too strong.

  In CN, the pi bonds also have asymmetric electron populations:
  CN has 4 electrons in pi (2 in each pi channel) - these are FULLY occupied.
  So the pi bonds are genuinely at full strength. The overprediction
  can't be fixed by electron counting alone.

  Possible explanations for residual CN error:
  1. The E_scale for CN (sqrt(E_C * E_N) = {np.sqrt(orbital_energy('C_2p')*orbital_energy('N_2p')):.2f})
     overpredicts the actual coupling strength
  2. The phase is slightly wrong for an asymmetric system
  3. Some correlation effect reduces pi bonding in 13e systems
""")
