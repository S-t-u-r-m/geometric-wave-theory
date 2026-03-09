"""
CO-BF Unified Analysis
=======================

The hydrogenic suppression [2*sqrt(Z1Z2)/(Z1+Z2)]^5 fixes BF but hurts CO.
There MUST be a rational explanation. Let's find it.

Key question: Why does the base formula work for CO (1.2% error)?
CO has Z_ratio=1.42, so the spatial mismatch is real — but E_scale already
captures some of it via the geometric mean sqrt(E1*E2).

Hypothesis: The base formula's accuracy for CO comes from E_scale already
including PART of the spatial correction. Adding the full hydrogenic
factor double-counts.

Let's decompose exactly what E_scale does vs what spatial overlap does.
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


# =============================================================================
# DECOMPOSITION: What does E_scale already capture?
# =============================================================================
print("=" * 90)
print("  E_SCALE DECOMPOSITION: geometric mean vs arithmetic mean")
print("=" * 90)
print(f"""
  E_scale = sqrt(E1 * E2) = geometric mean
  E_arith = (E1 + E2) / 2 = arithmetic mean

  For same-n p-orbitals: E = E_H * (Z/n)^2
  E_scale = E_H * Z1*Z2/n^2
  E_arith = E_H * (Z1^2 + Z2^2) / (2*n^2)

  Ratio = E_scale/E_arith = 2*Z1*Z2/(Z1^2+Z2^2)

  This is NOT the same as [2*sqrt(Z1Z2)/(Z1+Z2)]^2 = 4*Z1*Z2/(Z1+Z2)^2

  E_scale/E_arith = 2*Z1*Z2/(Z1^2+Z2^2) = 2/(Z1/Z2 + Z2/Z1)
  Hydrogenic_base^2 = 4*Z1*Z2/(Z1+Z2)^2
""")

print(f"{'Mol':<5} {'Z1':>6} {'Z2':>6} {'ratio':>6} {'E_geom/E_arith':>15} {'hydro_base^2':>13} {'geom covers':>12}")
print("-" * 75)

for name in ['N2', 'CO', 'NO', 'O2', 'BF', 'CN', 'B2', 'F2', 'Cl2']:
    mol = [m for m in molecules if m[0] == name][0]
    z1, z2 = Z_eff[mol[4]], Z_eff[mol[5]]
    ratio = max(z1,z2)/min(z1,z2)
    e_ratio = 2*z1*z2/(z1**2 + z2**2)  # geometric/arithmetic for energies
    h_base2 = 4*z1*z2/(z1+z2)**2  # hydrogenic coherence base^2

    # How much of the hydrogenic suppression does E_scale already capture?
    # If total suppression needed is base^5, and E_scale already gives e_ratio
    # vs geometric, then the "extra" suppression is base^5 / e_ratio
    h_base5 = h_base2 ** 2.5
    coverage = (1 - e_ratio) / (1 - h_base5) if h_base5 < 1 else 0

    print(f"{name:<5} {z1:6.3f} {z2:6.3f} {ratio:6.3f} {e_ratio:15.6f} {h_base2:13.6f} {coverage:12.3f}")


# =============================================================================
# HYPOTHESIS 1: Apply suppression only to PI, not SIGMA
# =============================================================================
print(f"\n{'='*90}")
print(f"  HYPOTHESIS 1: Suppression on PI bonds only")
print(f"  Physical argument: sigma is head-on (less affected by size mismatch)")
print(f"  pi is side-by-side (more sensitive to spatial extent)")
print(f"{'='*90}")

def compute_test(mol, hydro_exp=0, pi_only=False, use_radical=True):
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

    if l1 == 1 and l2 == 1 and n1 == n2_ and is_het and hydro_exp > 0:
        base = 2*np.sqrt(z1*z2)/(z1+z2)
        S_hydro = base ** hydro_exp
    else:
        S_hydro = 1.0

    half_sigma = use_radical and sigma_half_filled(bonds, ne, orb1, orb2)

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

        # Apply hydrogenic suppression
        if S_hydro < 1.0:
            if pi_only:
                # Only suppress pi bonds, not sigma
                if btype in ('pi', 'pi_anti'):
                    S = S * S_hydro
            else:
                if btype in ('pp_sigma', 'pi', 'pi_anti'):
                    S = S * S_hydro

        effective_count = count
        if half_sigma and btype == 'pp_sigma' and 'anti' not in btype:
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

    return D_cov + D_ion

print(f"\n  All bonds suppressed:")
for exp in [0, 3, 4, 5, 6]:
    errs = {m[0]: abs((compute_test(m, hydro_exp=exp, pi_only=False)-m[2])/m[2]*100) for m in molecules}
    ea = list(errs.values())
    w5 = sum(1 for e in ea if e < 5)
    w10 = sum(1 for e in ea if e < 10)
    print(f"  exp={exp}: avg={np.mean(ea):.1f}% w5={w5} w10={w10} | "
          f"BF={errs['BF']:.1f}% CO={errs['CO']:.1f}% NO={errs['NO']:.1f}% CN={errs['CN']:.1f}%")

print(f"\n  Pi-only suppressed:")
for exp in [0, 3, 4, 5, 6, 7, 8]:
    errs = {m[0]: abs((compute_test(m, hydro_exp=exp, pi_only=True)-m[2])/m[2]*100) for m in molecules}
    ea = list(errs.values())
    w5 = sum(1 for e in ea if e < 5)
    w10 = sum(1 for e in ea if e < 10)
    print(f"  exp={exp}: avg={np.mean(ea):.1f}% w5={w5} w10={w10} | "
          f"BF={errs['BF']:.1f}% CO={errs['CO']:.1f}% NO={errs['NO']:.1f}% CN={errs['CN']:.1f}%")


# =============================================================================
# HYPOTHESIS 2: What if the phase should depend on Z_eff?
# =============================================================================
print(f"\n{'='*90}")
print(f"  HYPOTHESIS 2: Z_eff-dependent phase")
print(f"  Current: phase = R/n^b + R/n^b (same for all Z)")
print(f"  Idea: the effective wavelength depends on Z through the Bohr radius")
print(f"{'='*90}")

# The Bohr radius for orbital (n,l,Z) is r ~ n^2/Z
# The effective wavenumber is k ~ Z/n^2 (reciprocal of radius)
# If phase ~ sum of k_i * R: phase = R * (Z1/n1^2 + Z2/n2^2)
# But we need to match the homonuclear result: phase = R/n + R/n = 2R/n
# For homonuclear Z1=Z2=Z: phase = R*(2Z/n^2) vs current 2R/n
# These differ by a factor Z/n. For n=2, Z~3-5, this is ~1.5-2.5
# Too large — the overall scale is wrong.

# Better: use Z_eff to MODIFY the phase, not replace it.
# phase = R/n^b * Z_eff^gamma + R/n^b * Z_eff^gamma
# For homonuclear: phase = 2*R/n^b * Z^gamma
# This changes homonuclear results unless gamma=0!

# What if we use the RATIO of Z values to modify the phase?
# For heteronuclear: phase *= [2*sqrt(Z1*Z2)/(Z1+Z2)]^delta
# This equals 1 for homonuclear, <1 for heteronuclear
# i.e., heteronuclear phase is SHORTER, giving larger |sin|

print(f"\n  Test: phase_hetero = phase * [2*sqrt(Z1Z2)/(Z1+Z2)]^delta")
print(f"  This SHORTENS the phase for heteronuclear -> |sin| changes")

def compute_phase_mod(mol, delta=0, hydro_exp=0, use_radical=True):
    """Test Z-dependent phase modification."""
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

    # Phase modification for heteronuclear
    z1, z2 = Z_eff[orb1], Z_eff[orb2]
    is_het = (orb1 != orb2)
    if l1 == 1 and l2 == 1 and n1 == n2_ and is_het and delta != 0:
        base = 2*np.sqrt(z1*z2)/(z1+z2)
        sigma_phase *= base ** delta

    n_pi_bond = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    n_pi_anti = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    is_full_anti = (n_pi_anti >= n_pi_bond) if n_pi_bond > 0 else True

    half_sigma = use_radical and sigma_half_filled(bonds, ne, orb1, orb2)

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

        effective_count = count
        if half_sigma and btype == 'pp_sigma' and 'anti' not in btype:
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

    return D_cov + D_ion

print()
for delta_val in [0, 1, 2, 3, 4, 5]:
    errs = {m[0]: abs((compute_phase_mod(m, delta=delta_val)-m[2])/m[2]*100) for m in molecules}
    ea = list(errs.values())
    w5 = sum(1 for e in ea if e < 5)
    w10 = sum(1 for e in ea if e < 10)
    print(f"  delta={delta_val}: avg={np.mean(ea):.1f}% w5={w5} w10={w10} | "
          f"BF={errs['BF']:.1f}% CO={errs['CO']:.1f}% NO={errs['NO']:.1f}% CN={errs['CN']:.1f}%")


# =============================================================================
# HYPOTHESIS 3: Use HARMONIC mean instead of geometric for E_scale
# =============================================================================
print(f"\n{'='*90}")
print(f"  HYPOTHESIS 3: Harmonic mean E_scale (includes natural suppression)")
print(f"  E_harm = 2*E1*E2/(E1+E2) = E_geom * [2*sqrt(E1E2)/(E1+E2)]")
print(f"  For same-n orbs: ratio = [2*sqrt(Z1Z2)/(Z1+Z2)]^2")
print(f"{'='*90}")

def compute_harmonic(mol, use_radical=True):
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
    # HARMONIC mean instead of geometric
    E_scale = 2*E1*E2/(E1+E2) if (E1+E2) > 0 else 0
    sigma_phase = R / n1**b1 + R / n2_**b2

    n_pi_bond = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    n_pi_anti = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    is_full_anti = (n_pi_anti >= n_pi_bond) if n_pi_bond > 0 else True

    half_sigma = use_radical and sigma_half_filled(bonds, ne, orb1, orb2)

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

        effective_count = count
        if half_sigma and btype == 'pp_sigma' and 'anti' not in btype:
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

    return D_cov + D_ion

print(f"\n  Harmonic E_scale changes for pp bonds:")
for name in ['N2', 'CO', 'NO', 'O2', 'BF', 'CN', 'F2', 'Cl2']:
    mol = [m for m in molecules if m[0] == name][0]
    z1, z2 = Z_eff[mol[4]], Z_eff[mol[5]]
    E1 = orbital_energy(mol[4])
    E2 = orbital_energy(mol[5])
    e_geom = np.sqrt(E1*E2)
    e_harm = 2*E1*E2/(E1+E2) if (E1+E2) > 0 else 0
    ratio = e_harm/e_geom
    print(f"  {name}: E_geom={e_geom:.2f}, E_harm={e_harm:.2f}, ratio={ratio:.4f}")

print(f"\n  Full results with harmonic E_scale:")
errs_harm = {m[0]: abs((compute_harmonic(m)-m[2])/m[2]*100) for m in molecules}
ea = list(errs_harm.values())
w5 = sum(1 for e in ea if e < 5)
w10 = sum(1 for e in ea if e < 10)
print(f"  avg={np.mean(ea):.1f}%, med={np.median(ea):.1f}%, w5={w5}/24, w10={w10}/24")
for name in ['N2', 'CO', 'NO', 'BF', 'CN', 'HF', 'OH', 'LiF', 'NaCl']:
    print(f"    {name}: {errs_harm[name]:.1f}%")


# =============================================================================
# HYPOTHESIS 4: E_scale interpolates between geometric and harmonic
# =============================================================================
print(f"\n{'='*90}")
print(f"  HYPOTHESIS 4: E_scale = geom * base^p (generalized mean)")
print(f"  p=0: geometric (current)")
print(f"  p=2: harmonic")
print(f"  p=5: full hydrogenic")
print(f"  Key: applied to E_scale, not as separate factor")
print(f"{'='*90}")

def compute_gen_mean(mol, p_exp=0, use_radical=True, pp_only=True):
    """E_scale with generalized Z-suppression built in."""
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
    E_scale_base = np.sqrt(E1 * E2)

    z1, z2 = Z_eff[orb1], Z_eff[orb2]
    is_het = (orb1 != orb2)

    # Z-suppression built into E_scale
    if l1 == 1 and l2 == 1 and n1 == n2_ and is_het and p_exp > 0:
        base = 2*np.sqrt(z1*z2)/(z1+z2)
        E_scale_pp = E_scale_base * base**p_exp
    else:
        E_scale_pp = E_scale_base

    sigma_phase = R / n1**b1 + R / n2_**b2

    n_pi_bond = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    n_pi_anti = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    is_full_anti = (n_pi_anti >= n_pi_bond) if n_pi_bond > 0 else True

    half_sigma = use_radical and sigma_half_filled(bonds, ne, orb1, orb2)

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

        effective_count = count
        if half_sigma and btype == 'pp_sigma' and 'anti' not in btype:
            effective_count = count * 0.5

        # Use pp E_scale for pp bonds, base E_scale for others
        if pp_only and btype in ('pp_sigma', 'pi', 'pi_anti'):
            e_s = E_scale_pp
        else:
            e_s = E_scale_base

        contribution = C_bond * e_s * S
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

    return D_cov + D_ion, D_cov, D_ion

# Note: E_scale * base^p applied to pp bonds is mathematically identical
# to multiplying S by base^p. So this is the same as the hydrogenic
# suppression. The difference is only conceptual (where we attribute
# the correction).

# But what if we apply it to E_scale AND have it affect the ionic
# calculation differently?

# Actually, the key insight might be: the hydrogenic suppression
# should affect D_cov (through E_scale), which then changes q,
# which changes D_ion. If we apply it BEFORE computing q, the
# ionic term compensates more. Let's check.

print(f"\n  Note: E_scale*base^p on pp bonds is same as S*base^p")
print(f"  But when applied to E_scale, it affects D_cov BEFORE the q calculation")
print(f"  giving more ionic compensation for strongly suppressed bonds.\n")

for p in [0, 2, 3, 4, 5]:
    errs = {}
    for m in molecules:
        De_pred, D_cov, D_ion = compute_gen_mean(m, p_exp=p)
        err = abs((De_pred - m[2]) / m[2] * 100)
        errs[m[0]] = (err, D_cov, D_ion)

    ea = [e[0] for e in errs.values()]
    w5 = sum(1 for e in ea if e < 5)
    w10 = sum(1 for e in ea if e < 10)
    print(f"  p={p}: avg={np.mean(ea):.1f}%, w5={w5}, w10={w10} | "
          f"BF={errs['BF'][0]:.1f}%(Dcov={errs['BF'][1]:.2f},Dion={errs['BF'][2]:.2f}) "
          f"CO={errs['CO'][0]:.1f}%(Dcov={errs['CO'][1]:.2f},Dion={errs['CO'][2]:.2f}) "
          f"CN={errs['CN'][0]:.1f}%")


# =============================================================================
# HYPOTHESIS 5: Back to basics — what if the formula IS right
# and BF/CN bond descriptions are wrong?
# =============================================================================
print(f"\n{'='*90}")
print(f"  HYPOTHESIS 5: Are the BF/CN bond lists correct?")
print(f"{'='*90}")

print(f"""
  BF: We use ('pp_sigma', 1), ('pi', 2) = triple bond (BO=3)
  CN: We use ('pp_sigma', 1), ('pi', 2) = triple bond (BO=3)

  But BF has unusual bonding:
  - B has only 3 valence electrons, F has 7
  - The triple bond requires "dative" back-donation from F to B
  - The sigma bond is a dative bond (both electrons from N/F lone pair)
  - This is physically different from a shared-pair sigma

  CN has 13 electrons:
  - sigma_2p is singly occupied (BO=2.5, we already handle this)
  - BUT even the pi bonds may be weakened by the unpaired sigma electron

  What if BF should be described as:
  - ('pp_sigma', 1), ('pi', 1) = BO=2 (only one strong pi, one weak dative pi)?
  - Or ('sp', 1), ('pi', 2) = different sigma character?

  Let's test: what bond order gives the right answer for BF?
""")

for bo_sigma in [0.0, 0.25, 0.5, 0.75, 1.0]:
    for bo_pi in [0.5, 1.0, 1.5, 2.0]:
        bonds_test = []
        if bo_sigma > 0:
            bonds_test.append(('pp_sigma', bo_sigma))
        if bo_pi > 0:
            bonds_test.append(('pi', bo_pi))
        mol_test = ('BF', 2.386, 7.81, bonds_test, 'B_2p', 'F_2p', 14)
        De = compute_gen_mean(mol_test, p_exp=0, use_radical=False)[0]
        err = (De - 7.81) / 7.81 * 100
        if abs(err) < 5:
            print(f"  BF with sigma={bo_sigma}, pi={bo_pi}: De={De:.3f}, err={err:+.1f}%  BO={bo_sigma+bo_pi:.1f}")


# =============================================================================
# HYPOTHESIS 6: The suppression should use the n-dependent Bohr radius
# =============================================================================
print(f"\n{'='*90}")
print(f"  HYPOTHESIS 6: Phase should use Z_eff-corrected wavelength")
print(f"  Instead of phase = R/n + R/n, use phase = R*Z1/n^2 + R*Z2/n^2")
print(f"  scaled to match homonuclear")
print(f"{'='*90}")

# For homonuclear with Z, phase = R/n (current) = 2R/n (sum of both atoms)
# New phase = 2*R*Z/n^2
# To match: 2R/n = 2R*Z/n^2 -> Z = n
# For n=2: matches when Z=2. For Z > 2, new phase is larger.

# What if we define:
# phase = R * [sqrt(Z1*Z2)/n^2] * C_phase
# where C_phase is chosen to match homonuclear at the average Z?
# For homonuclear: phase = R * Z/n^2 * C_phase = R/n -> C_phase = n/Z
# This is circular...

# Better: the effective wavenumber in the bond is the harmonic mean:
# k_eff = 2*k1*k2/(k1+k2) where k_i = Z_i/n_i^2
# For homonuclear: k_eff = k = Z/n^2
# For heteronuclear: k_eff = 2*Z1*Z2/(n^2*(Z1+Z2))
# Ratio to homonuclear: k_eff/k_homo = 2*Z1*Z2/((Z1+Z2)*Z_avg)
# This is complicated. Let me just test numerically.

# Simple approach: phase *= Z_geometric/Z_arithmetic
# = sqrt(Z1*Z2) / ((Z1+Z2)/2) = 2*sqrt(Z1*Z2)/(Z1+Z2)
# This is the "base" factor raised to power 1.

# For BF: base = 0.934, so phase is shortened by 6.6%
# Since phase for BF is ~2.39 (close to 3pi/4), shortening it
# moves sin(phase) to a DIFFERENT value.

# Actually, the current BF phase = 2.386 (= R since n=2, b=1 for both)
# sin(2.386) = 0.686
# sin(2.386 * 0.934) = sin(2.228) = sin(2.228) = 0.800 -> LARGER
# This would INCREASE D_cov, making it worse!

# What about LENGTHENING the phase? Use arithmetic/geometric:
# phase *= (Z1+Z2)/(2*sqrt(Z1*Z2)) = 1/base

# For BF: phase *= 1/0.934 = 1.071
# sin(2.386 * 1.071) = sin(2.555) = 0.548 -> SMALLER -> helps!

print(f"\n  Test: phase *= (Z1+Z2)/(2*sqrt(Z1*Z2)) for pp bonds")
print(f"  This LENGTHENS phase for heteronuclear -> changes |sin|")

for mol in molecules:
    name, R, De_exp, bonds, o1, o2, ne = mol
    if name not in ['N2', 'CO', 'BF', 'CN', 'NO']:
        continue
    z1, z2 = Z_eff[o1], Z_eff[o2]
    l1, l2 = get_l(o1), get_l(o2)
    if l1 != 1 or l2 != 1:
        continue

    base = 2*np.sqrt(z1*z2)/(z1+z2)
    n1 = get_n(o1)
    n2_ = get_n(o2)
    sigma_phase = R/n1 + R/n2_

    # Extended phase
    sigma_ext = sigma_phase / base
    pi_ext = sigma_ext * f_pi

    print(f"\n  {name}: Z=({z1:.3f},{z2:.3f}), base={base:.4f}")
    print(f"    sigma: {sigma_phase:.4f} -> {sigma_ext:.4f}")
    print(f"    |sin_sigma|: {abs(np.sin(sigma_phase)):.4f} -> {abs(np.sin(sigma_ext)):.4f}")
    print(f"    |sin_pi|:    {abs(np.sin(sigma_phase*f_pi)):.4f} -> {abs(np.sin(pi_ext)):.4f}")


# Now test with phase extension
print(f"\n  Scan: phase *= 1/base^delta for pp heteronuclear")
for delta_val in [0, 0.5, 1, 1.5, 2, 3]:
    errs = {}
    for m in molecules:
        name, R, De_exp, bonds, orb1, orb2, ne = m
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

        # Phase extension for heteronuclear pp
        if l1 == 1 and l2 == 1 and n1 == n2_ and is_het and delta_val > 0:
            base = 2*np.sqrt(z1*z2)/(z1+z2)
            sigma_phase_eff = sigma_phase / base**delta_val
        else:
            sigma_phase_eff = sigma_phase

        n_pi_bond = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
        n_pi_anti = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
        is_full_anti = (n_pi_anti >= n_pi_bond) if n_pi_bond > 0 else True

        half_sigma = sigma_half_filled(bonds, ne, orb1, orb2)

        D_cov = 0
        for btype, count in bonds:
            if 'sigma' in btype or btype in ('ss', 'sp'):
                phase = sigma_phase_eff
            else:
                phase = sigma_phase_eff * f_pi

            S = abs(np.sin(phase))
            if (h1+h2 > 0) and phase > pi:
                n_lobes = int(np.ceil(phase / pi))
                S = S / n_lobes
            if S < overlap_floor:
                S = overlap_floor

            effective_count = count
            if half_sigma and btype == 'pp_sigma' and 'anti' not in btype:
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
        De_pred = D_cov + D_ion
        err = abs((De_pred - De_exp) / De_exp * 100)
        errs[name] = err

    ea = list(errs.values())
    w5 = sum(1 for e in ea if e < 5)
    w10 = sum(1 for e in ea if e < 10)
    print(f"  delta={delta_val:.1f}: avg={np.mean(ea):.1f}% w5={w5} w10={w10} | "
          f"BF={errs['BF']:.1f}% CO={errs['CO']:.1f}% NO={errs['NO']:.1f}% CN={errs['CN']:.1f}%")
