"""
Per-Orbital Bond Coupling Test
==============================
Instead of one universal E_scale for all bond types in a molecule,
compute each bond's contribution using the SPECIFIC orbital pair involved.

Key difference from bond_predictions.py:
  - Triple bonds: sigma may use s+p or p+p; pi always uses p+p
  - Each bond type gets its own E_scale AND phase from its orbital pair
  - This tests whether "separate wave couplings" give better results

Uses corrected Cl Z_eff = 6.1161 and corrected phase formula b = 1 + beta*h.
"""

import numpy as np

pi = np.pi
E_H = 13.6057  # eV

# Parameters from d=3
d = 3
C_bond = pi / d
f_pi   = d**2 / (d**2 + 1)    # 9/10
alpha  = 1 - f_pi / d          # 7/10
beta   = (1 + f_pi) / 2        # 19/20
f_anti = 2*d / (2*d - 1)       # 6/5
c_ionic = 1.0 / (2*d + 1)      # 1/7

# Clementi-Raimondi Z_eff -- FULL table including 2s orbitals
Z_eff = {
    'H_1s':  1.0000,
    'Li_2s': 1.2792,
    'B_2s':  2.5762,  'B_2p':  2.4214,
    'C_2s':  3.2166,  'C_2p':  3.1358,
    'N_2s':  3.8474,  'N_2p':  3.8340,
    'O_2s':  4.4916,  'O_2p':  4.4532,
    'F_2s':  5.1276,  'F_2p':  5.0998,
    'Na_3s': 2.5074,
    'Cl_3s': 6.3682,  'Cl_3p': 6.1161,
}


def get_n(orb):
    return int(orb.split('_')[1][0])

def get_l(orb):
    return {'s': 0, 'p': 1}[orb.split('_')[1][1]]

def has_nodes(orb):
    n, l = get_n(orb), get_l(orb)
    return min(n - l - 1, 1)

def orbital_energy(orb):
    n, z = get_n(orb), Z_eff[orb]
    return E_H * (z / n)**2

def bond_order(bonds):
    bo = 0
    for bt, cnt in bonds:
        bo += cnt if 'anti' not in bt else -cnt
    return bo


# =============================================================================
# MOLECULE DATA -- per-orbital version
# =============================================================================
# Each bond now specifies its OWN orbital pair: (bond_type, count, orb_a, orb_b)
# This is the key difference from bond_predictions.py

molecules = [
    # --- Original 12 ---
    ('H2',   1.401,  4.745, [('ss', 1, 'H_1s', 'H_1s')]),
    ('Li2',  5.051,  1.056, [('ss', 1, 'Li_2s', 'Li_2s')]),
    ('B2',   3.005,  3.02,  [('pi', 2, 'B_2p', 'B_2p')]),
    ('C2',   2.348,  6.32,  [('pi', 2, 'C_2p', 'C_2p'),
                              ('sp_sigma', 1, 'C_2s', 'C_2p'),
                              ('sp_sigma_anti', 1, 'C_2s', 'C_2p')]),
    ('N2',   2.074,  9.759, [('pp_sigma', 1, 'N_2p', 'N_2p'),
                              ('pi', 2, 'N_2p', 'N_2p')]),
    ('O2',   2.282,  5.213, [('pp_sigma', 1, 'O_2p', 'O_2p'),
                              ('pi', 2, 'O_2p', 'O_2p'),
                              ('pi_anti', 1, 'O_2p', 'O_2p')]),
    ('F2',   2.668,  1.660, [('pp_sigma', 1, 'F_2p', 'F_2p'),
                              ('pi', 2, 'F_2p', 'F_2p'),
                              ('pi_anti', 2, 'F_2p', 'F_2p')]),
    ('Na2',  5.818,  0.746, [('ss', 1, 'Na_3s', 'Na_3s')]),
    ('Cl2',  3.757,  2.514, [('pp_sigma', 1, 'Cl_3p', 'Cl_3p'),
                              ('pi', 2, 'Cl_3p', 'Cl_3p'),
                              ('pi_anti', 2, 'Cl_3p', 'Cl_3p')]),
    ('HF',   1.733,  5.869, [('sp', 1, 'H_1s', 'F_2p')]),
    ('CO',   2.132, 11.225, [('pp_sigma', 1, 'C_2p', 'O_2p'),
                              ('pi', 2, 'C_2p', 'O_2p')]),
    ('NO',   2.175,  6.497, [('pp_sigma', 1, 'N_2p', 'O_2p'),
                              ('pi', 2, 'N_2p', 'O_2p'),
                              ('pi_anti', 1, 'N_2p', 'O_2p')]),
    # --- Blind test 12 ---
    ('OH',   1.834,  4.392, [('sp', 1, 'O_2p', 'H_1s')]),
    ('HCl',  2.409,  4.434, [('sp', 1, 'H_1s', 'Cl_3p')]),
    ('LiH',  3.015,  2.515, [('ss', 1, 'Li_2s', 'H_1s')]),
    ('LiF',  2.955,  5.939, [('sp', 1, 'Li_2s', 'F_2p')]),
    ('BH',   2.329,  3.42,  [('sp', 1, 'B_2p', 'H_1s')]),
    ('CH',   2.116,  3.47,  [('sp', 1, 'C_2p', 'H_1s')]),
    ('NH',   1.958,  3.57,  [('sp', 1, 'N_2p', 'H_1s')]),
    ('BF',   2.386,  7.81,  [('pp_sigma', 1, 'B_2p', 'F_2p'),
                              ('pi', 2, 'B_2p', 'F_2p')]),
    ('CN',   2.214,  7.72,  [('pp_sigma', 0.5, 'C_2p', 'N_2p'),
                              ('pi', 2, 'C_2p', 'N_2p')]),
    ('NaH',  3.566,  1.97,  [('ss', 1, 'Na_3s', 'H_1s')]),
    ('NaCl', 4.461,  4.23,  [('sp', 1, 'Na_3s', 'Cl_3p')]),
    ('H2O',  1.809,  5.117, [('sp', 1, 'O_2p', 'H_1s')]),
]


def compute_energy_v1(name, R, De_exp, bonds):
    """Version 1: Per-orbital E_scale and phase, same formula otherwise."""

    # Antibonding classification (need total pi counts)
    n_pi_bond = sum(c for bt, c, *_ in bonds if 'pi' in bt and 'anti' not in bt)
    n_pi_anti = sum(c for bt, c, *_ in bonds if 'pi' in bt and 'anti' in bt)
    is_full_anti = (n_pi_anti >= n_pi_bond) if n_pi_bond > 0 else True

    D_cov = 0
    # For ionic: use the "primary" orbital pair (first bond's orbitals)
    orb1_primary = bonds[0][2]
    orb2_primary = bonds[0][3]

    for entry in bonds:
        btype, count = entry[0], entry[1]
        orb_a, orb_b = entry[2], entry[3]

        n_a, l_a = get_n(orb_a), get_l(orb_a)
        n_b, l_b = get_n(orb_b), get_l(orb_b)
        h_a, h_b = has_nodes(orb_a), has_nodes(orb_b)

        # Per-orbital exponents
        a_a = 2 + (1 - 2*l_a) * alpha * h_a
        a_b = 2 + (1 - 2*l_b) * alpha * h_b
        b_a = 1 + beta * h_a
        b_b = 1 + beta * h_b

        # Per-orbital E_scale
        E_a = E_H / n_a**a_a
        E_b = E_H / n_b**a_b
        E_scale = np.sqrt(E_a * E_b)

        # Per-orbital phase
        sigma_phase = R / n_a**b_a + R / n_b**b_b

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

    # Ionic from primary orbitals
    eps1 = orbital_energy(orb1_primary)
    eps2 = orbital_energy(orb2_primary)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    D_ionic = c_ionic * q**2 * 2 * E_H / R

    return D_cov + D_ionic, D_cov, D_ionic, q


# =============================================================================
# VERSION 2: Sigma bonds use s orbitals where available
# For triple bonds, the sigma is often an s-s or s-p bond
# =============================================================================

# Alternative orbital assignments: sigma uses s orbitals
molecules_v2 = [
    ('H2',   1.401,  4.745, [('ss', 1, 'H_1s', 'H_1s')]),
    ('Li2',  5.051,  1.056, [('ss', 1, 'Li_2s', 'Li_2s')]),
    ('B2',   3.005,  3.02,  [('pi', 2, 'B_2p', 'B_2p')]),
    # C2: sigma from 2s overlap, pi from 2p
    ('C2',   2.348,  6.32,  [('pi', 2, 'C_2p', 'C_2p'),
                              ('ss_sigma', 1, 'C_2s', 'C_2s'),
                              ('ss_sigma_anti', 1, 'C_2s', 'C_2s')]),
    # N2: sigma from 2s mixing + pp_sigma; but MO theory says
    # the 3sigma_g is sp hybrid. Try: sigma from 2s, pi from 2p
    ('N2',   2.074,  9.759, [('ss_sigma', 1, 'N_2s', 'N_2s'),
                              ('pi', 2, 'N_2p', 'N_2p')]),
    ('O2',   2.282,  5.213, [('pp_sigma', 1, 'O_2p', 'O_2p'),
                              ('pi', 2, 'O_2p', 'O_2p'),
                              ('pi_anti', 1, 'O_2p', 'O_2p')]),
    ('F2',   2.668,  1.660, [('pp_sigma', 1, 'F_2p', 'F_2p'),
                              ('pi', 2, 'F_2p', 'F_2p'),
                              ('pi_anti', 2, 'F_2p', 'F_2p')]),
    ('Na2',  5.818,  0.746, [('ss', 1, 'Na_3s', 'Na_3s')]),
    ('Cl2',  3.757,  2.514, [('pp_sigma', 1, 'Cl_3p', 'Cl_3p'),
                              ('pi', 2, 'Cl_3p', 'Cl_3p'),
                              ('pi_anti', 2, 'Cl_3p', 'Cl_3p')]),
    ('HF',   1.733,  5.869, [('sp', 1, 'H_1s', 'F_2p')]),
    # CO: sigma from C_2s + O_2p (sp hybrid), pi from 2p+2p
    ('CO',   2.132, 11.225, [('sp_sigma', 1, 'C_2s', 'O_2p'),
                              ('pi', 2, 'C_2p', 'O_2p')]),
    ('NO',   2.175,  6.497, [('sp_sigma', 1, 'N_2s', 'O_2p'),
                              ('pi', 2, 'N_2p', 'O_2p'),
                              ('pi_anti', 1, 'N_2p', 'O_2p')]),
    # --- Blind test 12 ---
    ('OH',   1.834,  4.392, [('sp', 1, 'O_2p', 'H_1s')]),
    ('HCl',  2.409,  4.434, [('sp', 1, 'H_1s', 'Cl_3p')]),
    ('LiH',  3.015,  2.515, [('ss', 1, 'Li_2s', 'H_1s')]),
    ('LiF',  2.955,  5.939, [('sp', 1, 'Li_2s', 'F_2p')]),
    ('BH',   2.329,  3.42,  [('sp', 1, 'B_2p', 'H_1s')]),
    ('CH',   2.116,  3.47,  [('sp', 1, 'C_2p', 'H_1s')]),
    ('NH',   1.958,  3.57,  [('sp', 1, 'N_2p', 'H_1s')]),
    # BF: sigma from B_2s + F_2p, pi from 2p+2p
    ('BF',   2.386,  7.81,  [('sp_sigma', 1, 'B_2s', 'F_2p'),
                              ('pi', 2, 'B_2p', 'F_2p')]),
    # CN: sigma from C_2s + N_2p (half-filled), pi from 2p+2p
    ('CN',   2.214,  7.72,  [('sp_sigma', 0.5, 'C_2s', 'N_2p'),
                              ('pi', 2, 'C_2p', 'N_2p')]),
    ('NaH',  3.566,  1.97,  [('ss', 1, 'Na_3s', 'H_1s')]),
    ('NaCl', 4.461,  4.23,  [('sp', 1, 'Na_3s', 'Cl_3p')]),
    ('H2O',  1.809,  5.117, [('sp', 1, 'O_2p', 'H_1s')]),
]


# =============================================================================
# VERSION 3: Full MO-inspired orbital assignments
# Each MO maps to a specific orbital pair
# =============================================================================

molecules_v3 = [
    ('H2',   1.401,  4.745, [('ss', 1, 'H_1s', 'H_1s')]),
    ('Li2',  5.051,  1.056, [('ss', 1, 'Li_2s', 'Li_2s')]),
    ('B2',   3.005,  3.02,  [('pi', 2, 'B_2p', 'B_2p')]),
    # C2: 2sigma_g(2s-2s), 2sigma_u*(2s-2s), 1pi_u(2p-2p)x2
    ('C2',   2.348,  6.32,  [('pi', 2, 'C_2p', 'C_2p'),
                              ('ss_sigma', 1, 'C_2s', 'C_2s'),
                              ('ss_sigma_anti', 1, 'C_2s', 'C_2s')]),
    # N2: 2sigma_g(2s), 2sigma_u*(2s), 1pi_u(2p)x2, 3sigma_g(2p)
    # Net: 1 sigma(2s) + 2 pi(2p) + 1 sigma(2p) - 1 anti(2s) = BO 3
    # But 2s sigma + anti cancel. Net bonding = pp_sigma + 2*pi
    ('N2',   2.074,  9.759, [('pp_sigma', 1, 'N_2p', 'N_2p'),
                              ('pi', 2, 'N_2p', 'N_2p')]),
    ('O2',   2.282,  5.213, [('pp_sigma', 1, 'O_2p', 'O_2p'),
                              ('pi', 2, 'O_2p', 'O_2p'),
                              ('pi_anti', 1, 'O_2p', 'O_2p')]),
    ('F2',   2.668,  1.660, [('pp_sigma', 1, 'F_2p', 'F_2p'),
                              ('pi', 2, 'F_2p', 'F_2p'),
                              ('pi_anti', 2, 'F_2p', 'F_2p')]),
    ('Na2',  5.818,  0.746, [('ss', 1, 'Na_3s', 'Na_3s')]),
    ('Cl2',  3.757,  2.514, [('pp_sigma', 1, 'Cl_3p', 'Cl_3p'),
                              ('pi', 2, 'Cl_3p', 'Cl_3p'),
                              ('pi_anti', 2, 'Cl_3p', 'Cl_3p')]),
    ('HF',   1.733,  5.869, [('sp', 1, 'H_1s', 'F_2p')]),
    # CO: 3sigma has sp mixing. Try: sigma=C_2p+O_2p, pi=C_2p+O_2p
    # (In CO, the 3sigma is actually mostly C_2p character)
    # Plus 2sigma/2sigma* from 2s cancel
    ('CO',   2.132, 11.225, [('pp_sigma', 1, 'C_2p', 'O_2p'),
                              ('pi', 2, 'C_2p', 'O_2p')]),
    ('NO',   2.175,  6.497, [('pp_sigma', 1, 'N_2p', 'O_2p'),
                              ('pi', 2, 'N_2p', 'O_2p'),
                              ('pi_anti', 1, 'N_2p', 'O_2p')]),
    ('OH',   1.834,  4.392, [('sp', 1, 'O_2p', 'H_1s')]),
    ('HCl',  2.409,  4.434, [('sp', 1, 'H_1s', 'Cl_3p')]),
    ('LiH',  3.015,  2.515, [('ss', 1, 'Li_2s', 'H_1s')]),
    ('LiF',  2.955,  5.939, [('sp', 1, 'Li_2s', 'F_2p')]),
    ('BH',   2.329,  3.42,  [('sp', 1, 'B_2p', 'H_1s')]),
    ('CH',   2.116,  3.47,  [('sp', 1, 'C_2p', 'H_1s')]),
    ('NH',   1.958,  3.57,  [('sp', 1, 'N_2p', 'H_1s')]),
    # BF: sigma mostly B_2p character (like CO)
    ('BF',   2.386,  7.81,  [('pp_sigma', 1, 'B_2p', 'F_2p'),
                              ('pi', 2, 'B_2p', 'F_2p')]),
    ('CN',   2.214,  7.72,  [('pp_sigma', 0.5, 'C_2p', 'N_2p'),
                              ('pi', 2, 'C_2p', 'N_2p')]),
    ('NaH',  3.566,  1.97,  [('ss', 1, 'Na_3s', 'H_1s')]),
    ('NaCl', 4.461,  4.23,  [('sp', 1, 'Na_3s', 'Cl_3p')]),
    ('H2O',  1.809,  5.117, [('sp', 1, 'O_2p', 'H_1s')]),
]


# =============================================================================
# RUN ALL VERSIONS
# =============================================================================

def run_test(label, mol_list):
    print(f"\n{'='*85}")
    print(f"  {label}")
    print(f"{'='*85}")

    header = f"{'Mol':<7} {'De_exp':>7} {'De_pred':>7} {'err%':>7}  {'D_cov':>7} {'D_ion':>6} {'q':>5} {'BO':>4}"
    print(header)
    print("-" * 70)

    errs = []
    for mol in mol_list:
        name, R, De_exp, bonds = mol
        De_pred, D_cov, D_ion, q = compute_energy_v1(name, R, De_exp, bonds)
        bo = bond_order([(b[0], b[1]) for b in bonds])
        err = (De_pred - De_exp) / De_exp * 100
        errs.append(abs(err))

        flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
        print(f"{name:<7} {De_exp:7.3f} {De_pred:7.3f} {err:+6.1f}%  "
              f"{D_cov:7.3f} {D_ion:6.3f} {q:5.3f} {bo:4.1f} {flag}")

    print("-" * 70)
    n = len(errs)
    w2 = sum(1 for e in errs if e < 2)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    print(f"  avg = {np.mean(errs):.1f}%, median = {np.median(errs):.1f}%")
    print(f"  within 2%: {w2}/{n},  within 5%: {w5}/{n},  within 10%: {w10}/{n}")

    return errs


# V1: Same orbital assignments as bond_predictions.py (baseline with CN BO=2.5)
print("Per-orbital coupling test: does using different orbitals per bond type help?")
print(f"Parameters: C={C_bond:.4f}, f_pi={f_pi}, alpha={alpha}, beta={beta}")
print(f"            f_anti={f_anti}, c_ionic={c_ionic:.4f}")

errs_v1 = run_test("V1: BASELINE (same orbitals as bond_predictions.py, CN BO=2.5)", molecules)
errs_v2 = run_test("V2: SIGMA USES s ORBITALS (N2/CO/BF/CN sigma from 2s)", molecules_v2)
errs_v3 = run_test("V3: FULL MO (sigma=pp for all, 2s bonds cancel)", molecules_v3)


# =============================================================================
# COMPARISON TABLE
# =============================================================================
print(f"\n{'='*85}")
print(f"  COMPARISON: Which orbital assignment works best?")
print(f"{'='*85}")
print(f"\n{'Mol':<7} {'V1(base)':>10} {'V2(s-sig)':>10} {'V3(MO)':>10}  {'best':>6}")
print("-" * 50)

names = [m[0] for m in molecules]
for i, name in enumerate(names):
    e1, e2, e3 = errs_v1[i], errs_v2[i], errs_v3[i]
    best = 'V1' if e1 <= e2 and e1 <= e3 else 'V2' if e2 <= e3 else 'V3'
    # Mark if they differ
    diff = '' if abs(e1 - e2) < 0.01 and abs(e1 - e3) < 0.01 else ' <--'
    print(f"{name:<7} {e1:9.1f}% {e2:9.1f}% {e3:9.1f}%  {best:>6}{diff}")

print(f"\n{'':7} {'V1':>10} {'V2':>10} {'V3':>10}")
print(f"{'avg':7} {np.mean(errs_v1):9.1f}% {np.mean(errs_v2):9.1f}% {np.mean(errs_v3):9.1f}%")
print(f"{'w5':7} {sum(1 for e in errs_v1 if e<5):>9} {sum(1 for e in errs_v2 if e<5):>9} {sum(1 for e in errs_v3 if e<5):>9}")
print(f"{'w10':7} {sum(1 for e in errs_v1 if e<10):>9} {sum(1 for e in errs_v2 if e<10):>9} {sum(1 for e in errs_v3 if e<10):>9}")


# =============================================================================
# DETAILED: Show E_scale and phase for each bond in changed molecules
# =============================================================================
print(f"\n{'='*85}")
print(f"  DETAIL: Per-bond E_scale and phase for changed molecules")
print(f"{'='*85}")

changed_mols = ['C2', 'N2', 'CO', 'NO', 'BF', 'CN']
for version_label, mol_list in [("V1 (baseline)", molecules), ("V2 (s-sigma)", molecules_v2)]:
    print(f"\n--- {version_label} ---")
    for mol in mol_list:
        name = mol[0]
        if name not in changed_mols:
            continue
        R = mol[1]
        print(f"\n  {name} (R = {R:.3f} bohr):")
        for entry in mol[3]:
            btype, count = entry[0], entry[1]
            orb_a, orb_b = entry[2], entry[3]

            n_a, l_a = get_n(orb_a), get_l(orb_a)
            n_b, l_b = get_n(orb_b), get_l(orb_b)
            h_a, h_b = has_nodes(orb_a), has_nodes(orb_b)

            a_a = 2 + (1 - 2*l_a) * alpha * h_a
            a_b = 2 + (1 - 2*l_b) * alpha * h_b
            b_a = 1 + beta * h_a
            b_b = 1 + beta * h_b

            E_a = E_H / n_a**a_a
            E_b = E_H / n_b**a_b
            E_scale = np.sqrt(E_a * E_b)

            sigma_phase = R / n_a**b_a + R / n_b**b_b
            if 'sigma' in btype or btype in ('ss', 'sp'):
                phase = sigma_phase
            else:
                phase = sigma_phase * f_pi

            print(f"    {btype:>15} x{count}: {orb_a}+{orb_b}  "
                  f"E_s={E_scale:.3f}  ph={phase:.3f} ({phase/pi:.3f}pi)  "
                  f"|sin|={abs(np.sin(phase)):.4f}  "
                  f"contrib={C_bond * E_scale * abs(np.sin(phase)) * count:.3f}")
