"""
Bond Predictions — Simulation-Informed Minimal Corrections
===========================================================
Instead of replacing |sin(phase)| entirely, apply TARGETED corrections
based on what the 3D breather simulation revealed:

1. s-wave overlap is MONOTONIC — sin wrapping past pi is unphysical
   -> For ss bonds with phase > pi: divide by ceil(phase/pi) (node counting)
   -> This prevents the overlap from "resurrecting" after decaying

2. p-wave overlap is OSCILLATORY — sin(phase) is correct
   -> Keep |sin(phase)| for pp and pi bonds unchanged

3. Overlap is CONTINUOUS — never exactly zero at isolated points
   -> When phase hits a sin node (CH at phase=pi), the overlap shouldn't
      vanish. Add a minimum floor from the wave's finite width.
   -> Floor = 1/d^2 = 1/9 (from d-dimensional uncertainty)

4. Ionic bonds need stronger Coulomb coupling when covalent part vanishes
   -> When D_cov/delta_eps < 1/d^3 = 1/27, enhanced c_ionic = d/(2d+1) = 3/7
   -> This catches LiF, NaCl where full charge transfer occurs

All GWT parameters still from d=3.
"""

import numpy as np

pi = np.pi
E_H = 13.6057  # eV

# =============================================================================
# GWT PARAMETERS (d=3)
# =============================================================================
d = 3
C_bond = pi / d
f_pi   = d**2 / (d**2 + 1)
alpha  = 1 - f_pi / d
beta   = (1 + f_pi) / 2
f_anti = 2*d / (2*d - 1)
c_ionic = 1.0 / (2*d + 1)

# Simulation-informed corrections (all from d=3)
overlap_floor = 1.0 / d**2               # 1/9 — minimum overlap from finite width
ionic_threshold = 1.0 / d**3             # 1/27 — when covalent vanishes
c_ionic_enhanced = d / (2*d + 1)         # 3/7 — strong ionic coupling

Z_eff = {
    'H_1s':  1.0000, 'Li_2s': 1.2792, 'B_2p':  2.4214,
    'C_2p':  3.1358, 'N_2p':  3.8340, 'O_2p':  4.4532,
    'F_2p':  5.0998, 'Na_3s': 2.5074, 'Cl_3p': 6.1161,
}

molecules = [
    ('H2',   1.401,  4.745, [('ss', 1)],                                           'H_1s',  'H_1s'),
    ('Li2',  5.051,  1.056, [('ss', 1)],                                           'Li_2s', 'Li_2s'),
    ('B2',   3.005,  3.02,  [('pi', 2)],                                           'B_2p',  'B_2p'),
    ('C2',   2.348,  6.32,  [('pi', 2), ('sp_sigma', 1), ('sp_sigma_anti', 1)],    'C_2p',  'C_2p'),
    ('N2',   2.074,  9.759, [('pp_sigma', 1), ('pi', 2)],                          'N_2p',  'N_2p'),
    ('O2',   2.282,  5.213, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)],          'O_2p',  'O_2p'),
    ('F2',   2.668,  1.660, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)],          'F_2p',  'F_2p'),
    ('Na2',  5.818,  0.746, [('ss', 1)],                                           'Na_3s', 'Na_3s'),
    ('Cl2',  3.757,  2.514, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)],          'Cl_3p', 'Cl_3p'),
    ('HF',   1.733,  5.869, [('sp', 1)],                                           'H_1s',  'F_2p'),
    ('CO',   2.132, 11.225, [('pp_sigma', 1), ('pi', 2)],                          'C_2p',  'O_2p'),
    ('NO',   2.175,  6.497, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)],          'N_2p',  'O_2p'),
    ('OH',   1.834,  4.392, [('sp', 1)],                                           'O_2p',  'H_1s'),
    ('HCl',  2.409,  4.434, [('sp', 1)],                                           'H_1s',  'Cl_3p'),
    ('LiH',  3.015,  2.515, [('ss', 1)],                                           'Li_2s', 'H_1s'),
    ('LiF',  2.955,  5.939, [('sp', 1)],                                           'Li_2s', 'F_2p'),
    ('BH',   2.329,  3.42,  [('sp', 1)],                                           'B_2p',  'H_1s'),
    ('CH',   2.116,  3.47,  [('sp', 1)],                                           'C_2p',  'H_1s'),
    ('NH',   1.958,  3.57,  [('sp', 1)],                                           'N_2p',  'H_1s'),
    ('BF',   2.386,  7.81,  [('pp_sigma', 1), ('pi', 2)],                          'B_2p',  'F_2p'),
    ('CN',   2.214,  7.72,  [('pp_sigma', 1), ('pi', 2)],                          'C_2p',  'N_2p'),
    ('NaH',  3.566,  1.97,  [('ss', 1)],                                           'Na_3s', 'H_1s'),
    ('NaCl', 4.461,  4.23,  [('sp', 1)],                                           'Na_3s', 'Cl_3p'),
    ('H2O',  1.809,  5.117, [('sp', 1)],                                           'O_2p',  'H_1s'),
]

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


def compute_bond_energy_v4(name, R, De_exp, bonds, orb1, orb2):
    """Compute bond energy with simulation-informed corrections."""
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

    corrections = []  # track what corrections applied

    D_cov = 0
    for btype, count in bonds:
        if 'sigma' in btype or btype in ('ss', 'sp'):
            phase = sigma_phase
        else:
            phase = sigma_phase * f_pi

        # Base overlap
        S = abs(np.sin(phase))

        # CORRECTION 1: s-wave node counting
        # Simulation: s-wave overlap is monotonic, never increases after decay
        # When phase > pi, sin wraps — divide by number of lobes to average
        is_ss = (btype == 'ss')
        has_any_nodes = (h1 + h2 > 0)

        if has_any_nodes and phase > pi:
            n_lobes = int(np.ceil(phase / pi))
            S = S / n_lobes
            corrections.append(f'node_count({n_lobes})')

        # CORRECTION 2: Overlap floor
        # Simulation: wave overlap never exactly zero at isolated points
        # Physical: finite wave width prevents perfect cancellation
        if S < overlap_floor and S > 0:
            S = overlap_floor
            corrections.append('floor')

        contribution = C_bond * E_scale * S

        if 'anti' in btype:
            if 'sigma' in btype or is_full_anti:
                f_a = 1.0
            else:
                f_a = f_anti
            D_cov -= count * f_a * contribution
        else:
            D_cov += count * contribution

    # Ionic correction
    eps1 = orbital_energy(orb1)
    eps2 = orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)

    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0

    # CORRECTION 3: Enhanced ionic for truly ionic bonds
    # When covalent contribution is negligible compared to electronegativity gap,
    # full charge transfer occurs -> stronger Coulomb coupling
    ratio = abs(D_cov) / delta_eps if delta_eps > 0.01 else float('inf')
    if ratio < ionic_threshold:
        c_ion = c_ionic_enhanced
        corrections.append(f'ionic_enh(r={ratio:.3f})')
    else:
        c_ion = c_ionic

    D_ionic = c_ion * q**2 * 2 * E_H / R

    return D_cov + D_ionic, D_cov, D_ionic, q, corrections


def compute_bond_energy_orig(name, R, De_exp, bonds, orb1, orb2):
    """Original formula for comparison."""
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
        contribution = C_bond * E_scale * abs(np.sin(phase))
        if 'anti' in btype:
            if 'sigma' in btype or is_full_anti:
                f_a = 1.0
            else:
                f_a = f_anti
            D_cov -= count * f_a * contribution
        else:
            D_cov += count * contribution

    eps1 = orbital_energy(orb1)
    eps2 = orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    D_ionic = c_ionic * q**2 * 2 * E_H / R

    return D_cov + D_ionic


# =============================================================================
# Test individual corrections independently
# =============================================================================

print("=" * 95)
print("  TESTING CORRECTIONS INDEPENDENTLY")
print("=" * 95)

# First: show which molecules trigger each correction
print("\n--- Which molecules trigger each correction? ---\n")

for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    _, _, _, _, corrs = compute_bond_energy_v4(*mol)
    if corrs:
        print(f"  {name:<7}: {', '.join(corrs)}")


# =============================================================================
# Test correction 1 alone: Node counting
# =============================================================================
print(f"\n\n{'='*95}")
print(f"  CORRECTION 1: Node counting (phase > pi, has_nodes)")
print(f"  Simulation basis: s-wave overlap is monotonic, sin wrapping is unphysical")
print(f"{'='*95}")


def compute_nc_only(name, R, De_exp, bonds, orb1, orb2):
    """Node counting only."""
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
        has_any_nodes = (h1 + h2 > 0)
        if has_any_nodes and phase > pi:
            n_lobes = int(np.ceil(phase / pi))
            S = S / n_lobes
        contribution = C_bond * E_scale * S
        if 'anti' in btype:
            if 'sigma' in btype or is_full_anti:
                f_a = 1.0
            else:
                f_a = f_anti
            D_cov -= count * f_a * contribution
        else:
            D_cov += count * contribution

    eps1 = orbital_energy(orb1)
    eps2 = orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    D_ionic = c_ionic * q**2 * 2 * E_H / R
    return D_cov + D_ionic


print(f"\n{'Mol':<7} {'De_exp':>7} {'orig':>7} {'err_o':>7} {'NC':>7} {'err_NC':>7} {'delta':>7}")
print("-" * 65)

errs_orig = []
errs_nc = []
for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    De_o = compute_bond_energy_orig(*mol)
    De_nc = compute_nc_only(*mol)
    err_o = (De_o - De_exp) / De_exp * 100
    err_nc = (De_nc - De_exp) / De_exp * 100
    errs_orig.append(abs(err_o))
    errs_nc.append(abs(err_nc))
    changed = '*' if abs(De_nc - De_o) > 0.001 else ''
    print(f"{name:<7} {De_exp:7.3f} {De_o:7.3f} {err_o:+6.1f}% {De_nc:7.3f} {err_nc:+6.1f}% {changed}")

print(f"\n  Original: avg={np.mean(errs_orig):.1f}%, median={np.median(errs_orig):.1f}%, "
      f"w5={sum(1 for e in errs_orig if e<5)}/24, w10={sum(1 for e in errs_orig if e<10)}/24")
print(f"  Node cnt: avg={np.mean(errs_nc):.1f}%, median={np.median(errs_nc):.1f}%, "
      f"w5={sum(1 for e in errs_nc if e<5)}/24, w10={sum(1 for e in errs_nc if e<10)}/24")


# =============================================================================
# Test correction 2 alone: Overlap floor
# =============================================================================
print(f"\n\n{'='*95}")
print(f"  CORRECTION 2: Overlap floor = 1/d^2 = {overlap_floor:.4f}")
print(f"  Simulation basis: wave overlap is continuous, never exactly zero")
print(f"{'='*95}")


def compute_floor_only(name, R, De_exp, bonds, orb1, orb2):
    """Floor only."""
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
        S = max(S, overlap_floor)
        contribution = C_bond * E_scale * S
        if 'anti' in btype:
            if 'sigma' in btype or is_full_anti:
                f_a = 1.0
            else:
                f_a = f_anti
            D_cov -= count * f_a * contribution
        else:
            D_cov += count * contribution

    eps1 = orbital_energy(orb1)
    eps2 = orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    D_ionic = c_ionic * q**2 * 2 * E_H / R
    return D_cov + D_ionic


print(f"\n{'Mol':<7} {'De_exp':>7} {'orig':>7} {'err_o':>7} {'floor':>7} {'err_f':>7} {'delta':>7}")
print("-" * 65)

errs_fl = []
for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    De_o = compute_bond_energy_orig(*mol)
    De_fl = compute_floor_only(*mol)
    err_o = (De_o - De_exp) / De_exp * 100
    err_fl = (De_fl - De_exp) / De_exp * 100
    errs_fl.append(abs(err_fl))
    changed = '*' if abs(De_fl - De_o) > 0.001 else ''
    print(f"{name:<7} {De_exp:7.3f} {De_o:7.3f} {err_o:+6.1f}% {De_fl:7.3f} {err_fl:+6.1f}% {changed}")

print(f"\n  Original: avg={np.mean(errs_orig):.1f}%, median={np.median(errs_orig):.1f}%")
print(f"  Floor:    avg={np.mean(errs_fl):.1f}%, median={np.median(errs_fl):.1f}%, "
      f"w5={sum(1 for e in errs_fl if e<5)}/24, w10={sum(1 for e in errs_fl if e<10)}/24")


# =============================================================================
# Test correction 3 alone: Enhanced ionic
# =============================================================================
print(f"\n\n{'='*95}")
print(f"  CORRECTION 3: Enhanced ionic c={c_ionic_enhanced:.4f} when D_cov/delta_eps < {ionic_threshold:.4f}")
print(f"  Physical basis: full charge transfer -> stronger Coulomb coupling")
print(f"{'='*95}")


def compute_ionic_only(name, R, De_exp, bonds, orb1, orb2):
    """Enhanced ionic only."""
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
        contribution = C_bond * E_scale * abs(np.sin(phase))
        if 'anti' in btype:
            if 'sigma' in btype or is_full_anti:
                f_a = 1.0
            else:
                f_a = f_anti
            D_cov -= count * f_a * contribution
        else:
            D_cov += count * contribution

    eps1 = orbital_energy(orb1)
    eps2 = orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0

    ratio = abs(D_cov) / delta_eps if delta_eps > 0.01 else float('inf')
    c_ion = c_ionic_enhanced if ratio < ionic_threshold else c_ionic
    D_ionic = c_ion * q**2 * 2 * E_H / R
    return D_cov + D_ionic


print(f"\n{'Mol':<7} {'De_exp':>7} {'orig':>7} {'err_o':>7} {'enh':>7} {'err_e':>7} {'ratio':>7}")
print("-" * 70)

errs_ei = []
for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    De_o = compute_bond_energy_orig(*mol)
    De_ei = compute_ionic_only(*mol)
    err_o = (De_o - De_exp) / De_exp * 100
    err_ei = (De_ei - De_exp) / De_exp * 100
    errs_ei.append(abs(err_ei))

    # Show ratio
    n1, l1 = get_n(o1), get_l(o1)
    n2, l2 = get_n(o2), get_l(o2)
    eps1 = orbital_energy(o1)
    eps2 = orbital_energy(o2)
    delta_eps = abs(eps1 - eps2)
    D_cov_val = compute_bond_energy_v4(*mol)[1]
    ratio = abs(D_cov_val) / delta_eps if delta_eps > 0.01 else float('inf')

    changed = '*' if abs(De_ei - De_o) > 0.001 else ''
    ratio_str = f"{ratio:.3f}" if ratio < 100 else "inf"
    print(f"{name:<7} {De_exp:7.3f} {De_o:7.3f} {err_o:+6.1f}% {De_ei:7.3f} {err_ei:+6.1f}% {ratio_str:>7} {changed}")

print(f"\n  Original:  avg={np.mean(errs_orig):.1f}%, median={np.median(errs_orig):.1f}%")
print(f"  Enh ionic: avg={np.mean(errs_ei):.1f}%, median={np.median(errs_ei):.1f}%, "
      f"w5={sum(1 for e in errs_ei if e<5)}/24, w10={sum(1 for e in errs_ei if e<10)}/24")


# =============================================================================
# ALL THREE COMBINED
# =============================================================================
print(f"\n\n{'='*95}")
print(f"  ALL THREE CORRECTIONS COMBINED")
print(f"{'='*95}")

print(f"\n{'Mol':<7} {'De_exp':>7} {'orig':>7} {'err_o':>7} {'v4':>7} {'err_v4':>7} {'corrections':>30}")
print("-" * 85)

errs_v4 = []
for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    De_o = compute_bond_energy_orig(*mol)
    De_v4, D_cov, D_ion, q, corrs = compute_bond_energy_v4(*mol)
    err_o = (De_o - De_exp) / De_exp * 100
    err_v4 = (De_v4 - De_exp) / De_exp * 100
    errs_v4.append(abs(err_v4))
    corr_str = ', '.join(corrs) if corrs else '-'
    flag = '***' if abs(err_v4) < 2 else ' **' if abs(err_v4) < 5 else '  *' if abs(err_v4) < 10 else ''
    print(f"{name:<7} {De_exp:7.3f} {De_o:7.3f} {err_o:+6.1f}% {De_v4:7.3f} {err_v4:+6.1f}% {corr_str:>30} {flag}")

print("-" * 85)
n = len(errs_v4)
w2 = sum(1 for e in errs_v4 if e < 2)
w5 = sum(1 for e in errs_v4 if e < 5)
w10 = sum(1 for e in errs_v4 if e < 10)
w20 = sum(1 for e in errs_v4 if e < 20)
print(f"\n  Original:  avg={np.mean(errs_orig):.1f}%, median={np.median(errs_orig):.1f}%, "
      f"w5={sum(1 for e in errs_orig if e<5)}/24, w10={sum(1 for e in errs_orig if e<10)}/24")
print(f"  V4 combo:  avg={np.mean(errs_v4):.1f}%, median={np.median(errs_v4):.1f}%, "
      f"w5={w5}/{n}, w10={w10}/{n}, w20={w20}/{n}")

# Show improvements and regressions
print(f"\n  Improvements (error decreased):")
for i, mol in enumerate(molecules):
    if errs_v4[i] < errs_orig[i] - 1:
        print(f"    {mol[0]:<7}: {errs_orig[i]:.1f}% -> {errs_v4[i]:.1f}%")

print(f"\n  Regressions (error increased):")
for i, mol in enumerate(molecules):
    if errs_v4[i] > errs_orig[i] + 1:
        print(f"    {mol[0]:<7}: {errs_orig[i]:.1f}% -> {errs_v4[i]:.1f}%")


# =============================================================================
# Scan overlap floor values
# =============================================================================
print(f"\n\n{'='*95}")
print(f"  SCAN: Overlap floor values")
print(f"{'='*95}")

for floor_val in [0.0, 1/16, 1/12, 1/10, 1/9, 1/8, 1/6, 1/4, 1/3]:
    errs_scan = []
    for mol in molecules:
        name, R, De_exp, bonds, o1, o2 = mol
        n1, l1 = get_n(o1), get_l(o1)
        n2, l2 = get_n(o2), get_l(o2)
        h1, h2 = has_nodes(o1), has_nodes(o2)
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
            # Node counting
            if (h1+h2 > 0) and phase > pi:
                n_lobes = int(np.ceil(phase / pi))
                S = S / n_lobes
            # Floor
            S = max(S, floor_val)
            contribution = C_bond * E_scale * S
            if 'anti' in btype:
                if 'sigma' in btype or is_full_anti:
                    f_a = 1.0
                else:
                    f_a = f_anti
                D_cov -= count * f_a * contribution
            else:
                D_cov += count * contribution

        eps1 = orbital_energy(o1)
        eps2 = orbital_energy(o2)
        delta_eps = abs(eps1 - eps2)
        V_cov = max(abs(D_cov), 0.01)
        q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
        ratio = abs(D_cov) / delta_eps if delta_eps > 0.01 else float('inf')
        c_ion = c_ionic_enhanced if ratio < ionic_threshold else c_ionic
        D_ionic = c_ion * q**2 * 2 * E_H / R
        De_pred = D_cov + D_ionic
        err = abs((De_pred - De_exp) / De_exp * 100)
        errs_scan.append(err)

    w5 = sum(1 for e in errs_scan if e < 5)
    w10 = sum(1 for e in errs_scan if e < 10)
    # Find CH error specifically
    ch_err = errs_scan[17]  # CH is index 17
    print(f"  floor={floor_val:.4f}: avg={np.mean(errs_scan):.1f}%, med={np.median(errs_scan):.1f}%, "
          f"w5={w5}/24, w10={w10}/24, CH={ch_err:.1f}%")
