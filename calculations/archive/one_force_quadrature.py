"""
ONE-FORCE QUADRATURE MODEL
============================
D = sqrt(D_cov^2 + (dE/d)^2)

Physics: Two coupled waves with frequency mismatch.
  D_cov = resonant coupling = C * E_scale * |sin(phase)|
  dE/d  = off-resonance correction = energy mismatch / d

The total is the MAGNITUDE of the coupling vector:
  sqrt(resonant^2 + detuning^2)

This is ONE force (wave resonance) with two components.
The 1/d comes from the detuning spreading over d spatial dimensions.

dE_formula = |E_H/n1^a - E_H/n2^a|  (no Z_eff — fully self-consistent!)

For homonuclear: dE = 0, so D = |D_cov|. Same as current.
For heteronuclear: dE > 0 adds a guaranteed positive contribution.
"""
import numpy as np

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


def compute_cov(mol):
    """Compute covalent-only energy with phase fix."""
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
    n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)

    a1 = 2 + (1 - 2*l1) * alpha_n * h1
    a2 = 2 + (1 - 2*l2) * alpha_n * h2
    b1 = 1 + beta_n * h1
    b2 = 1 + beta_n * h2

    E_scale = np.sqrt(E_H / n1**a1 * E_H / n2**a2)
    k1 = 1.0 / n1**b1
    k2 = 1.0 / n2**b2
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

    # Formula energy difference (no Z_eff!)
    dE_formula = abs(E_H / n1**a1 - E_H / n2**a2)

    # Orbital energy difference (uses Z_eff)
    dE_orbital = abs(orb_energy(orb1) - orb_energy(orb2))

    return D_cov, dE_formula, dE_orbital, sigma_phase


# =============================================================================
# SCAN c values for D = sqrt(D_cov^2 + c * dE^2)
# =============================================================================
print("=" * 95)
print("  SCANNING: D = sqrt(D_cov^2 + c * dE_formula^2)")
print("  dE_formula = |E_H/n1^a - E_H/n2^a|  (self-consistent, no Z_eff)")
print("=" * 95)
print()

d_fracs = [
    ('0 (cov only)', 0),
    ('1/(2d+1)^2 = 1/49', 1/(2*d+1)**2),
    ('1/(d^2+1) = 1/10', 1/(d**2+1)),
    ('1/d^2 = 1/9', 1/d**2),
    ('1/(2d) = 1/6', 1/(2*d)),
    ('1/(d+1) = 1/4', 1/(d+1)),
    ('1/d = 1/3', 1/d),
    ('1/(d*(d+1)) = 1/12', 1/(d*(d+1))),
    ('f_pi/(d^2+1) = 9/100', f_pi/(d**2+1)),
    ('1/(2*(d^2+1)) = 1/20', 1/(2*(d**2+1))),
    ('1/(4*d) = 1/12', 1/(4*d)),
]

print(f"{'c (label)':>30} {'avg%':>7} {'med%':>7} {'<5%':>5} {'<10%':>5} {'<20%':>5}")
print("-" * 65)

best = (999, '', 0)
for label, c_val in d_fracs:
    errs = []
    for mol in molecules:
        D_cov, dE, _, phase = compute_cov(mol)
        D_tot = np.sqrt(D_cov**2 + c_val * dE**2)
        err = abs((D_tot - mol[2]) / mol[2] * 100)
        errs.append(err)

    avg = np.mean(errs)
    med = np.median(errs)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"{label:>30} {avg:7.1f} {med:7.1f} {w5:5d} {w10:5d} {w20:5d}")
    if avg < best[0]:
        best = (avg, label, c_val)


# =============================================================================
# DETAILED VIEW: D = sqrt(D_cov^2 + dE^2/d^2)  [c = 1/9]
# =============================================================================
c_best = 1/d**2
print()
print("=" * 95)
print(f"  DETAIL: D = sqrt(D_cov^2 + dE^2/d^2)  with c = 1/d^2 = 1/{d**2}")
print("=" * 95)
print()

print(f"{'Mol':<6} {'De_exp':>7} {'D_cov':>7} {'dE':>7} {'dE/d':>7} {'D_pred':>7} "
      f"{'err%':>7}  {'err_old%':>8}")
print("-" * 75)

errs_new = []
errs_old = []
for mol in molecules:
    name, R, De_exp = mol[0], mol[1], mol[2]
    orb1, orb2 = mol[4], mol[5]
    D_cov, dE, dE_orb, phase = compute_cov(mol)

    # New model: quadrature
    D_new = np.sqrt(D_cov**2 + c_best * dE**2)
    err_new = (D_new - De_exp) / De_exp * 100
    errs_new.append(abs(err_new))

    # Old model: cov + ionic
    V = max(abs(D_cov), 0.01)
    q = dE_orb / np.sqrt(dE_orb**2 + (2*V)**2) if dE_orb > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / R
    D_old = D_cov + Di
    err_old = (D_old - De_exp) / De_exp * 100
    errs_old.append(abs(err_old))

    flag = '***' if abs(err_new) < 2 else ' **' if abs(err_new) < 5 else '  *' if abs(err_new) < 10 else ''
    v = 'BETTER' if abs(err_new) < abs(err_old) - 1 else 'WORSE' if abs(err_new) > abs(err_old) + 1 else '~'

    print(f"{name:<6} {De_exp:7.3f} {D_cov:7.3f} {dE:7.3f} {dE/d:7.3f} "
          f"{D_new:7.3f} {err_new:+6.1f}% {err_old:+7.1f}% {flag} {v}")

w5_n = sum(1 for e in errs_new if e < 5)
w10_n = sum(1 for e in errs_new if e < 10)
w20_n = sum(1 for e in errs_new if e < 20)
w5_o = sum(1 for e in errs_old if e < 5)
w10_o = sum(1 for e in errs_old if e < 10)
w20_o = sum(1 for e in errs_old if e < 20)

print()
print(f"  NEW (quadrature):  avg={np.mean(errs_new):.1f}%, med={np.median(errs_new):.1f}%, "
      f"<5%:{w5_n}/24, <10%:{w10_n}/24, <20%:{w20_n}/24")
print(f"  OLD (cov+ionic):   avg={np.mean(errs_old):.1f}%, med={np.median(errs_old):.1f}%, "
      f"<5%:{w5_o}/24, <10%:{w10_o}/24, <20%:{w20_o}/24")

nb = sum(1 for i in range(24) if errs_new[i] < errs_old[i] - 1)
nw = sum(1 for i in range(24) if errs_new[i] > errs_old[i] + 1)
print(f"  Better: {nb}, Worse: {nw}, Same: {24-nb-nw}")


# =============================================================================
# INTERPRETATION
# =============================================================================
print()
print("=" * 95)
print("  PHYSICAL INTERPRETATION")
print("=" * 95)
print()
print("D = sqrt(D_cov^2 + (dE/d)^2)")
print()
print("This is the MAGNITUDE of a 2D coupling vector:")
print("  Component 1 (resonant): D_cov = C * E_scale * |sin(phase)|")
print("  Component 2 (detuning): dE/d  = energy mismatch / d")
print()
print("The two components are ORTHOGONAL (independent):")
print("  - Resonant: wave overlap along the bond axis")
print("  - Detuning: frequency mismatch, spread over d dimensions")
print()
print("This is ONE FORCE (wave resonance) in a 2D coupling space.")
print("The sqrt adds them in quadrature, like Pythagoras.")
print()
print(f"All constants from d={d}:")
print(f"  C     = pi/d       = {pi/d:.6f}")
print(f"  f_pi  = d^2/(d^2+1)= {f_pi:.6f}")
print(f"  alpha = 7/10       = {alpha_n:.6f}")
print(f"  beta  = 19/20      = {beta_n:.6f}")
print(f"  f_anti= 6/5        = {f_anti:.6f}")
print(f"  c_det = 1/d^2      = {1/d**2:.6f} (detuning coefficient)")
print()
print("NO Z_eff in the formula! Fully self-consistent from GWT.")
print("dE only depends on quantum numbers (n, l) through the a exponent.")
print()

# Check which atoms have which dE
print("Formula energies E_H/n^a:")
for orb in sorted(Z_eff.keys()):
    n = get_n(orb); l = get_l(orb); h = has_nodes(orb)
    a = 2 + (1 - 2*l) * alpha_n * h
    E = E_H / n**a
    print(f"  {orb:>8}: n={n}, a={a:.2f}, E={E:.4f} eV")

print()
print("dE values between orbital types:")
orbs = sorted(Z_eff.keys())
unique_types = {}
for orb in orbs:
    n = get_n(orb); l = get_l(orb); h = has_nodes(orb)
    a = 2 + (1 - 2*l) * alpha_n * h
    E = E_H / n**a
    key = f"n={n},a={a:.2f}"
    if key not in unique_types:
        unique_types[key] = (orb, E)

types = list(unique_types.values())
print(f"{'Type A':>12} {'Type B':>12} {'dE':>8} {'dE/d':>8}")
print("-" * 45)
for i in range(len(types)):
    for j in range(i+1, len(types)):
        orb_a, E_a = types[i]
        orb_b, E_b = types[j]
        dE = abs(E_a - E_b)
        print(f"{orb_a:>12} {orb_b:>12} {dE:8.4f} {dE/d:8.4f}")
