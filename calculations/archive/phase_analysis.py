"""
PHASE ANALYSIS: What phase does each molecule NEED?
====================================================
The ionic coupling problem is really a PHASE problem.
Molecules in the ION regime have b ≈ 0.05 for p-orbitals
with nodes, making the effective wavevector unphysically large.

Strategy:
1. For each molecule, find the phase that matches experiment
   (assuming the rest of the formula is correct)
2. Compare needed vs computed phase
3. Look for a corrected b formula that gives the right phases

The key physics: b controls the effective de Broglie wavevector
  k_eff = 1/n^b
For s-orbitals: b > 1 (nodes extend wavelength) — WORKS
For p-orbitals: b = 1 - 0.95 = 0.05 (angular momentum shrinks wavelength)
  But this makes k(3p) ≈ k(1s), which is unphysical.
  3p orbital is MUCH larger than 1s.
"""
import numpy as np
from scipy.optimize import brentq

pi = np.pi
E_H = 13.6057
d = 3
C_bond = pi / d
f_pi = d**2 / (d**2 + 1)
alpha_n = 1 - f_pi / d
beta_n = (1 + f_pi) / 2
f_anti = 2*d / (2*d - 1)
c_ionic = 1.0 / (2*d + 1)

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

def get_params(orb):
    n, l, h = get_n(orb), get_l(orb), has_nodes(orb)
    return n, l, h, 2 + (1-2*l)*alpha_n*h, 1 + (1-2*l)*beta_n*h


def compute_De(mol, sigma_phase_override=None):
    """Compute De with optional phase override."""
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1, a1, b1 = get_params(orb1)
    n2, l2, h2, a2, b2 = get_params(orb2)
    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)

    if sigma_phase_override is not None:
        sigma_phase = sigma_phase_override
    else:
        sigma_phase = R * (1.0/n1**b1 + 1.0/n2**b2)

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

    eps1 = E_H * (Z_eff[orb1]/n1)**2
    eps2 = E_H * (Z_eff[orb2]/n2)**2
    de = abs(eps1 - eps2)
    Vc = max(abs(D_cov), 0.01)
    q = de / np.sqrt(de**2 + (2*Vc)**2) if de > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / R

    return D_cov + Di, D_cov, Di, sigma_phase


# =============================================================================
# PART 1: What phase does each molecule need?
# =============================================================================
print("=" * 105)
print("  PHASE DIAGNOSIS: Current vs Needed phase for each molecule")
print("=" * 105)
print()
print(f"{'Mol':<6} {'orbs':>10} {'De_exp':>7} {'ph_now':>7} {'ph/pi':>6} {'k1':>6} {'k2':>6} "
      f"{'b1':>5} {'b2':>5} {'ph_need':>7} {'pn/pi':>6} {'k1_n':>6} {'k2_n':>6}")
print("-" * 105)

for mol in molecules:
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1, a1, b1 = get_params(orb1)
    n2, l2, h2, a2, b2 = get_params(orb2)

    k1 = 1.0 / n1**b1
    k2 = 1.0 / n2**b2
    sigma_phase = R * (k1 + k2)

    # Find phase that matches experiment
    # De(phase) is oscillatory, so search near current phase
    De_curr = compute_De(mol)[0]
    err_curr = (De_curr - De_exp) / De_exp * 100

    # Search for the nearest phase that gives the right De
    best_phase = sigma_phase
    best_diff = abs(De_curr - De_exp)

    for test_phase in np.arange(0.1, 4*pi, 0.001):
        De_test = compute_De(mol, sigma_phase_override=test_phase)[0]
        diff = abs(De_test - De_exp)
        # Prefer phases close to current
        if diff < best_diff * 0.99:
            # But also penalize very different phases
            if diff < 0.01 or diff < best_diff * 0.5:
                best_phase = test_phase
                best_diff = diff

    # What k values would give this phase?
    # phase = R * (k1 + k2), so k1+k2 = phase/R
    # If we keep the ratio k1/k2 the same:
    k_sum_needed = best_phase / R
    k_ratio = k1 / k2 if k2 > 0.001 else 1
    k2_needed = k_sum_needed / (1 + k_ratio)
    k1_needed = k2_needed * k_ratio

    orb_str = f"{orb1[-2:]}/{orb2[-2:]}"
    print(f"{name:<6} {orb_str:>10} {De_exp:7.3f} {sigma_phase:7.3f} {sigma_phase/pi:6.3f} "
          f"{k1:6.3f} {k2:6.3f} {b1:5.2f} {b2:5.2f} "
          f"{best_phase:7.3f} {best_phase/pi:6.3f} {k1_needed:6.3f} {k2_needed:6.3f}")


# =============================================================================
# PART 2: Examine the b exponent in detail
# =============================================================================
print()
print("=" * 105)
print("  B EXPONENT ANALYSIS")
print("=" * 105)
print()

print("Current formula: b = 1 + (1-2l) * beta * has_nodes")
print(f"  beta = {beta_n:.4f}")
print()
print(f"{'orbital':>8} {'n':>2} {'l':>2} {'h':>2} {'b':>6} {'k=1/n^b':>8} {'rB=n^2/Z':>8} {'1/rB':>8}")
print("-" * 55)

for orb in sorted(Z_eff.keys()):
    n, l, h, a, b = get_params(orb)
    k = 1.0 / n**b
    rB = n**2 / Z_eff[orb]
    print(f"{orb:>8} {n:2d} {l:2d} {h:2d} {b:6.3f} {k:8.4f} {rB:8.3f} {1/rB:8.4f}")

print()
print("Notice: Cl_3p has b=0.050, k=0.9447 ~ k(H_1s)=1.0")
print("        But Cl_3p Bohr radius = 1.842 >> H_1s radius = 1.000")
print("        The wavevector should be INVERSELY related to orbital size!")
print()

# What if k = Z_eff / n^2 = 1/r_B ?  (wavevector ~ inverse Bohr radius)
print("=" * 105)
print("  ALTERNATIVE: k = c / r_B = c * Z_eff / n^2  (wavevector ~ 1/orbital_size)")
print("=" * 105)
print()

# Try different k definitions
k_models = {
    'current (1/n^b)': lambda orb: 1.0 / get_n(orb)**get_params(orb)[4],
    'Z/n^2 (1/r_B)': lambda orb: Z_eff[orb] / get_n(orb)**2,
    'Z/n (1/r_B^0.5)': lambda orb: Z_eff[orb] / get_n(orb),
    '1/n': lambda orb: 1.0 / get_n(orb),
    'sqrt(Z)/n': lambda orb: np.sqrt(Z_eff[orb]) / get_n(orb),
    'Z^0.7/n^1.4': lambda orb: Z_eff[orb]**0.7 / get_n(orb)**1.4,
}

for k_name, k_func in k_models.items():
    print(f"\n--- Model: k = {k_name} ---")
    errs = []
    for mol in molecules:
        name, R, De_exp, bonds, orb1, orb2 = mol
        n1, l1, h1, a1, _ = get_params(orb1)
        n2, l2, h2, a2, _ = get_params(orb2)
        E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)

        k1 = k_func(orb1)
        k2 = k_func(orb2)
        sigma_phase = R * (k1 + k2)

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

        eps1 = orb_energy(orb1)
        eps2 = orb_energy(orb2)
        de = abs(eps1 - eps2)
        Vc = max(abs(D_cov), 0.01)
        q = de / np.sqrt(de**2 + (2*Vc)**2) if de > 0 else 0
        Di = c_ionic * q**2 * 2 * E_H / R
        D_tot = D_cov + Di
        err = (D_tot - De_exp) / De_exp * 100
        errs.append(abs(err))

    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"  avg={np.mean(errs):.1f}%, med={np.median(errs):.1f}%, "
          f"<5%:{w5}/24, <10%:{w10}/24, <20%:{w20}/24")


# =============================================================================
# PART 3: What if we use k = Z_eff^alpha / n^beta with d-derived exponents?
# =============================================================================
print()
print("=" * 105)
print("  SCAN: k = Z^a / n^b with (a, b) derived from d=3")
print("=" * 105)
print()

# Test various d-fraction exponent pairs
d_fracs = {
    '0': 0, '1/(d+1)': 1/(d+1), '1/d': 1/d, '2/d': 2/d,
    '(d-1)/d': (d-1)/d, '1': 1, 'd/(d+1)': d/(d+1),
    '(d+1)/d': (d+1)/d, 'd/(d-1)': d/(d-1), '2': 2,
    '1/2': 0.5, '3/2': 1.5
}

results = []
for a_name, a_val in d_fracs.items():
    for b_name, b_val in d_fracs.items():
        if b_val == 0:
            continue  # k = Z^a / n^0 = Z^a, no n dependence

        errs = []
        for mol in molecules:
            name, R, De_exp, bonds, orb1, orb2 = mol
            n1, l1, h1, a1, _ = get_params(orb1)
            n2, l2, h2, a2, _ = get_params(orb2)
            E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)

            k1 = Z_eff[orb1]**a_val / get_n(orb1)**b_val
            k2 = Z_eff[orb2]**a_val / get_n(orb2)**b_val
            sigma_phase = R * (k1 + k2)

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

            eps1 = orb_energy(orb1)
            eps2 = orb_energy(orb2)
            de = abs(eps1 - eps2)
            Vc = max(abs(D_cov), 0.01)
            q = de / np.sqrt(de**2 + (2*Vc)**2) if de > 0 else 0
            Di = c_ionic * q**2 * 2 * E_H / R
            D_tot = D_cov + Di
            err = (D_tot - De_exp) / De_exp * 100
            errs.append(abs(err))

        avg = np.mean(errs)
        w5 = sum(1 for e in errs if e < 5)
        w10 = sum(1 for e in errs if e < 10)
        results.append((avg, a_name, a_val, b_name, b_val, w5, w10))

results.sort()
print(f"{'a (Z^a)':>12} {'b (n^b)':>12} {'avg%':>7} {'<5%':>5} {'<10%':>5}")
print("-" * 50)
for avg, a_name, a_val, b_name, b_val, w5, w10 in results[:15]:
    print(f"{a_name:>12} {b_name:>12} {avg:7.1f} {w5:5d} {w10:5d}")

# Show the best model in detail
print()
avg, a_name, a_val, b_name, b_val, _, _ = results[0]
print(f"Best: k = Z^{a_name} / n^{b_name} = Z^{a_val:.4f} / n^{b_val:.4f}")
print()
print(f"{'Mol':<6} {'De_exp':>7} {'De_pred':>7} {'err%':>7} {'ph/pi':>6} {'k1':>6} {'k2':>6}")
print("-" * 55)

for mol in molecules:
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1, a1, _ = get_params(orb1)
    n2, l2, h2, a2, _ = get_params(orb2)
    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)

    k1 = Z_eff[orb1]**a_val / get_n(orb1)**b_val
    k2 = Z_eff[orb2]**a_val / get_n(orb2)**b_val
    sigma_phase = R * (k1 + k2)

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

    eps1 = orb_energy(orb1)
    eps2 = orb_energy(orb2)
    de = abs(eps1 - eps2)
    Vc = max(abs(D_cov), 0.01)
    q = de / np.sqrt(de**2 + (2*Vc)**2) if de > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / R
    D_tot = D_cov + Di
    err = (D_tot - De_exp) / De_exp * 100
    print(f"{name:<6} {De_exp:7.3f} {D_tot:7.3f} {err:+6.1f}% {sigma_phase/pi:6.3f} {k1:6.3f} {k2:6.3f}")
