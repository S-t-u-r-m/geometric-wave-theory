"""
BF and CN Failure Analysis
===========================
Both are heteronuclear triple bonds (pp_sigma + 2*pi) that the V4 formula
OVERPREDICTS by ~27-31%. Let's understand exactly why.

Questions:
1. Which bond component is too large? sigma or pi?
2. How do they compare to homonuclear triple bonds N2 and CO?
3. Is the E_scale (geometric mean) the issue for asymmetric Z_eff?
4. Does the phase structure differ from N2/CO?
5. Can simulation results tell us anything?
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

# =============================================================================
# The 5 triple bonds in our set
# =============================================================================
triple_bonds = [
    ('N2',  2.074, 9.759, 'N_2p',  'N_2p'),
    ('CO',  2.132, 11.225, 'C_2p',  'O_2p'),
    ('BF',  2.386, 7.81,  'B_2p',  'F_2p'),
    ('CN',  2.214, 7.72,  'C_2p',  'N_2p'),
    ('NO',  2.175, 6.497, 'N_2p',  'O_2p'),  # BO=2, but triple-like structure
]

bonds_triple = [('pp_sigma', 1), ('pi', 2)]
bonds_no = [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)]

print("=" * 90)
print("  DETAILED BREAKDOWN OF TRIPLE BONDS")
print("=" * 90)

print(f"\n{'Mol':<5} {'R':>5} {'Z1':>6} {'Z2':>6} {'Z2/Z1':>6} {'E1':>7} {'E2':>7} {'Escale':>7} {'E_ratio':>7}")
print("-" * 70)

for name, R, De_exp, o1, o2 in triple_bonds:
    z1, z2 = Z_eff[o1], Z_eff[o2]
    e1 = orbital_energy(o1)
    e2 = orbital_energy(o2)
    e_scale = np.sqrt(e1 * e2)
    e_arith = (e1 + e2) / 2
    print(f"{name:<5} {R:5.3f} {z1:6.3f} {z2:6.3f} {max(z1,z2)/min(z1,z2):6.3f} "
          f"{e1:7.3f} {e2:7.3f} {e_scale:7.3f} {e_scale/e_arith:7.3f}")

print("\nNote: E_ratio = geometric_mean / arithmetic_mean (always <= 1)")
print("More asymmetric Z -> ratio further below 1 -> geom mean overestimates less")

# =============================================================================
# Phase and overlap breakdown
# =============================================================================
print(f"\n{'='*90}")
print(f"  PHASE AND OVERLAP BREAKDOWN")
print(f"{'='*90}")

print(f"\n{'Mol':<5} {'sig_ph':>7} {'pi_ph':>7} {'|sin_s|':>8} {'|sin_p|':>8} "
      f"{'D_sig':>7} {'D_2pi':>7} {'D_cov':>7} {'D_ion':>6} {'D_tot':>7} {'De_exp':>7} {'err%':>6}")
print("-" * 95)

for name, R, De_exp, o1, o2 in triple_bonds:
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
    pi_phase = sigma_phase * f_pi

    sin_s = abs(np.sin(sigma_phase))
    sin_p = abs(np.sin(pi_phase))

    D_sig = C_bond * E_scale * sin_s
    D_2pi = 2 * C_bond * E_scale * sin_p

    bonds = bonds_no if name == 'NO' else bonds_triple

    if name == 'NO':
        # NO has 1 pi_anti
        D_anti = f_anti * C_bond * E_scale * sin_p
        D_cov = D_sig + D_2pi - D_anti
    else:
        D_cov = D_sig + D_2pi

    eps1, eps2 = orbital_energy(o1), orbital_energy(o2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    D_ion = c_ionic * q**2 * 2 * E_H / R

    D_tot = D_cov + D_ion
    err = (D_tot - De_exp) / De_exp * 100

    print(f"{name:<5} {sigma_phase:7.3f} {pi_phase:7.3f} {sin_s:8.4f} {sin_p:8.4f} "
          f"{D_sig:7.3f} {D_2pi:7.3f} {D_cov:7.3f} {D_ion:6.3f} {D_tot:7.3f} {De_exp:7.3f} {err:+5.1f}%")

# =============================================================================
# Key question: What SHOULD the energy be?
# =============================================================================
print(f"\n{'='*90}")
print(f"  WHAT CORRECTION WOULD FIX BF AND CN?")
print(f"{'='*90}")

print(f"\n  For each molecule, what effective overlap S_eff would match experiment?\n")

print(f"{'Mol':<5} {'De_exp':>7} {'D_ion':>6} {'D_cov_need':>10} {'D_cov_got':>10} "
      f"{'ratio':>6} {'S_eff':>6} {'S_orig':>6}")
print("-" * 70)

for name, R, De_exp, o1, o2 in triple_bonds:
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
    pi_phase = sigma_phase * f_pi

    sin_s = abs(np.sin(sigma_phase))
    sin_p = abs(np.sin(pi_phase))

    bonds = bonds_no if name == 'NO' else bonds_triple

    if name == 'NO':
        D_anti = f_anti * C_bond * E_scale * sin_p
        D_cov = C_bond * E_scale * sin_s + 2 * C_bond * E_scale * sin_p - D_anti
    else:
        D_cov = C_bond * E_scale * (sin_s + 2 * sin_p)

    eps1, eps2 = orbital_energy(o1), orbital_energy(o2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    D_ion = c_ionic * q**2 * 2 * E_H / R

    D_cov_need = De_exp - D_ion
    ratio = D_cov_need / D_cov if D_cov > 0 else 0

    # What would the effective S be?
    # D_cov = C * E_scale * (sin_s + 2*sin_p) for triple bonds
    # Need D_cov_need = C * E_scale * S_eff * 3
    S_eff = D_cov_need / (3 * C_bond * E_scale) if E_scale > 0 else 0
    S_orig = (sin_s + 2*sin_p) / 3

    print(f"{name:<5} {De_exp:7.3f} {D_ion:6.3f} {D_cov_need:10.3f} {D_cov:10.3f} "
          f"{ratio:6.3f} {S_eff:6.3f} {S_orig:6.3f}")

# =============================================================================
# The Z asymmetry angle
# =============================================================================
print(f"\n{'='*90}")
print(f"  Z_eff ASYMMETRY ANALYSIS")
print(f"{'='*90}")
print(f"\n  Hypothesis: heteronuclear overlap is suppressed by Z mismatch")
print(f"  In wave mechanics: waves of different wavelength interfere less efficiently")
print(f"  lambda ~ 1/Z_eff, so mismatch -> reduced coherent overlap\n")

for name, R, De_exp, o1, o2 in triple_bonds:
    z1, z2 = Z_eff[o1], Z_eff[o2]
    z_ratio = max(z1, z2) / min(z1, z2)

    # Various suppression factors
    harm_mean_ratio = 2*z1*z2 / (z1+z2)**2  # = (harmonic/arithmetic)^?
    # Actually: 2*z1*z2/(z1+z2)^2 = 1 - (z1-z2)^2/(z1+z2)^2
    # This equals (harmonic_mean / geometric_mean)^2 ... no
    # 4*z1*z2/(z1+z2)^2 is better known
    coherence = 4*z1*z2 / (z1+z2)**2  # = 1 for equal, <1 for unequal

    print(f"  {name}: Z=({z1:.3f}, {z2:.3f}), ratio={z_ratio:.3f}, "
          f"coherence=4Z1Z2/(Z1+Z2)^2 = {coherence:.4f}")

# =============================================================================
# Test: coherence suppression factor
# =============================================================================
print(f"\n{'='*90}")
print(f"  TEST: Coherence suppression for heteronuclear bonds")
print(f"  S -> S * [4*Z1*Z2 / (Z1+Z2)^2]  for p-p bonds")
print(f"{'='*90}")

all_mols = [
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

overlap_floor = 1.0 / (d+1)  # 1/4
ionic_threshold = 1.0 / d**3
c_ionic_enhanced = d / (2*d + 1)

def compute_with_coherence(mol, use_coherence=True):
    name, R, De_exp, bonds, orb1, orb2 = mol
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

    # Coherence factor for p-p bonds
    z1, z2 = Z_eff[orb1], Z_eff[orb2]
    coherence = 4*z1*z2 / (z1+z2)**2  # 1 for homonuclear

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
        if S < overlap_floor:
            S = overlap_floor

        # Coherence suppression for pp bonds
        if use_coherence and l1 == 1 and l2 == 1 and btype in ('pp_sigma', 'pi', 'pi_anti'):
            S = S * coherence

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

print(f"\n{'Mol':<7} {'De_exp':>7} {'V4':>7} {'err_v4':>7} {'V4+coh':>7} {'err_coh':>7} {'coh':>6}")
print("-" * 60)

errs_v4 = []
errs_coh = []
for mol in all_mols:
    name = mol[0]
    De_exp = mol[2]
    De_v4 = compute_with_coherence(mol, use_coherence=False)
    De_coh = compute_with_coherence(mol, use_coherence=True)
    err_v4 = (De_v4 - De_exp) / De_exp * 100
    err_coh = (De_coh - De_exp) / De_exp * 100
    errs_v4.append(abs(err_v4))
    errs_coh.append(abs(err_coh))

    z1, z2 = Z_eff[mol[4]], Z_eff[mol[5]]
    coh = 4*z1*z2 / (z1+z2)**2

    changed = '*' if abs(De_coh - De_v4) > 0.001 else ''
    print(f"{name:<7} {De_exp:7.3f} {De_v4:7.3f} {err_v4:+6.1f}% {De_coh:7.3f} {err_coh:+6.1f}% {coh:6.4f} {changed}")

print("-" * 60)
print(f"  V4:      avg={np.mean(errs_v4):.1f}%, med={np.median(errs_v4):.1f}%, "
      f"w5={sum(1 for e in errs_v4 if e<5)}/24, w10={sum(1 for e in errs_v4 if e<10)}/24")
print(f"  V4+coh:  avg={np.mean(errs_coh):.1f}%, med={np.median(errs_coh):.1f}%, "
      f"w5={sum(1 for e in errs_coh if e<5)}/24, w10={sum(1 for e in errs_coh if e<10)}/24")

# =============================================================================
# Alternative: Does the SIMULATION tell us about heteronuclear coupling?
# =============================================================================
print(f"\n{'='*90}")
print(f"  SIMULATION COMPARISON: mode (1,3) vs mode (3,3)")
print(f"{'='*90}")
print(f"""
  From breather_3d_gpu_output.log:
  - Mode (3,3) p-wave: oscillatory, peak at R~7
  - Mode (1,3) mixed s+p: ALWAYS REPULSIVE
  - Mode (1,3) s-wave: not directly tested for p-wave

  BF maps to: B(n=2,Z=2.42) + F(n=2,Z=5.10) — SAME n, DIFFERENT Z
  CN maps to: C(n=2,Z=3.14) + N(n=2,Z=3.83) — SAME n, DIFFERENT Z

  The simulation tested different MODE NUMBERS (n=1 vs n=3), which is
  a proxy for different wavenumbers. But BF/CN have same n, different Z.

  What matters for coherent overlap: wavelength match.
  In breather picture: wavelength ~ 1/eta = 1/sin(n*gamma)
  For SAME n but different Z_eff: the spatial scale differs because
  the radial wavefunction scales as r -> Z*r.

  So BF with Z_ratio = 5.10/2.42 = 2.1 has VERY different spatial scales.
  CN with Z_ratio = 3.83/3.14 = 1.22 has mild mismatch.

  This supports the coherence suppression hypothesis, but the 4Z1Z2/(Z1+Z2)^2
  factor may not be strong enough for BF.
""")

# =============================================================================
# Try stronger suppression: (2*sqrt(Z1*Z2)/(Z1+Z2))^p for various p
# =============================================================================
print(f"{'='*90}")
print(f"  SCAN: Coherence power factor")
print(f"  S -> S * [2*sqrt(Z1*Z2)/(Z1+Z2)]^p")
print(f"{'='*90}")

# Note: 2*sqrt(Z1*Z2)/(Z1+Z2) = sqrt(4*Z1*Z2/(Z1+Z2)^2) = sqrt(coherence)
# So coherence^(p/2) is the same as [2*sqrt(Z1*Z2)/(Z1+Z2)]^p

for p in [0, 1, 2, 3, 4, 5, 6]:
    errs_p = []
    bf_err = cn_err = co_err = n2_err = no_err = 0
    for mol in all_mols:
        name, R, De_exp, bonds, orb1, orb2 = mol
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
        coh = (4*z1*z2 / (z1+z2)**2) ** (p/2)

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

            # Coherence for pp bonds
            if p > 0 and l1 == 1 and l2 == 1 and btype in ('pp_sigma', 'pi', 'pi_anti'):
                S = S * coh

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
        De_pred = D_cov + D_ion
        err = abs((De_pred - De_exp) / De_exp * 100)
        errs_p.append(err)

        if name == 'BF': bf_err = err
        elif name == 'CN': cn_err = err
        elif name == 'CO': co_err = err
        elif name == 'N2': n2_err = err
        elif name == 'NO': no_err = err

    w5 = sum(1 for e in errs_p if e < 5)
    w10 = sum(1 for e in errs_p if e < 10)
    print(f"  p={p}: avg={np.mean(errs_p):.1f}%, med={np.median(errs_p):.1f}%, "
          f"w5={w5}/24, w10={w10}/24 | "
          f"N2={n2_err:.1f}% CO={co_err:.1f}% NO={no_err:.1f}% BF={bf_err:.1f}% CN={cn_err:.1f}%")
