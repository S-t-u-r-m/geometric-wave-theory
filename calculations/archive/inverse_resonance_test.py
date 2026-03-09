"""
INVERSE RESONANCE MODEL TEST
==============================
Key insight: the detuning correction should be LARGE when D_cov is small
(weak resonance) and SMALL when D_cov is large (strong resonance).

Models tested:
  A) D = (D_cov^2 + (dE/d)^2) / D_cov           -- inverse resonance
  B) D = D_cov + (dE/d)^2 / (D_cov + epsilon)    -- regularized inverse
  C) D = D_cov * (1 + (dE/(d*D_cov))^2)^(1/2)    -- relative quadrature (= A when expanded)
  D) D = D_cov + c * dE^2 / E_scale               -- perturbation theory (dE^2/V)
  E) D = sqrt(D_cov^2 + c * dE^2 * sin_deficit^2) -- sin-weighted detuning
  F) D = D_cov + c * dE * (1 - |sin(ph)|)         -- linear detuning * sin deficit
  G) D = sqrt(D_cov^2 + c * dE^2) with c adaptive -- c scales with sin deficit
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
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
    n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)

    a1 = 2 + (1 - 2*l1) * alpha_n * h1
    a2 = 2 + (1 - 2*l2) * alpha_n * h2
    b1 = 1 + beta_n * h1
    b2 = 1 + beta_n * h2

    E1 = E_H / n1**a1
    E2 = E_H / n2**a2
    E_scale = np.sqrt(E1 * E2)
    k1 = 1.0 / n1**b1
    k2 = 1.0 / n2**b2
    sigma_phase = R * (k1 + k2)

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    sin_sigma = abs(np.sin(sigma_phase))
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss', 'sp')) else sigma_phase * f_pi
        cont = C_bond * E_scale * abs(np.sin(ph))
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    dE = abs(E1 - E2)
    dE_orb = abs(orb_energy(orb1) - orb_energy(orb2))
    return D_cov, dE, dE_orb, E_scale, sigma_phase, sin_sigma


def old_model(mol):
    D_cov, dE, dE_orb, E_scale, phase, sin_s = compute_cov(mol)
    V = max(abs(D_cov), 0.01)
    q = dE_orb / np.sqrt(dE_orb**2 + (2*V)**2) if dE_orb > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / mol[1]
    return D_cov + Di


# =============================================================================
# MODEL DEFINITIONS
# =============================================================================
def model_A(mol, c=1.0/d**2):
    """Inverse resonance: D = (D_cov^2 + c*dE^2) / |D_cov|"""
    D_cov, dE, _, E_scale, phase, sin_s = compute_cov(mol)
    Dc = abs(D_cov)
    if Dc < 0.01:
        # Fallback for near-zero cov
        return np.sqrt(D_cov**2 + c * dE**2)
    return (D_cov**2 + c * dE**2) / Dc

def model_B(mol, c=1.0/d**2):
    """Regularized inverse: D = |D_cov| + c*dE^2 / (|D_cov| + eps)"""
    D_cov, dE, _, E_scale, phase, sin_s = compute_cov(mol)
    Dc = abs(D_cov)
    eps_vals = [0.5, 1.0, 2.0, E_scale * 0.1]
    # Use eps = 1.0 as default
    eps = 1.0
    return Dc + c * dE**2 / (Dc + eps)

def model_D(mol, c=1.0/d**2):
    """Perturbation: D = D_cov + c * dE^2 / E_scale"""
    D_cov, dE, _, E_scale, phase, sin_s = compute_cov(mol)
    return abs(D_cov) + c * dE**2 / E_scale

def model_E(mol, c=1.0/d**2):
    """Sin-weighted detuning: sqrt(D_cov^2 + c*dE^2*(1-|sin|)^2)"""
    D_cov, dE, _, E_scale, phase, sin_s = compute_cov(mol)
    sin_deficit = (1.0 - sin_s)**2
    return np.sqrt(D_cov**2 + c * dE**2 * sin_deficit)

def model_F(mol, c=1.0/d):
    """Linear detuning * sin deficit: D_cov + c*dE*(1-|sin|)"""
    D_cov, dE, _, E_scale, phase, sin_s = compute_cov(mol)
    return abs(D_cov) + c * dE * (1.0 - sin_s)

def model_G(mol, c_base=1.0/d**2):
    """Adaptive c: c = c_base * (1-|sin|), larger when sin is small"""
    D_cov, dE, _, E_scale, phase, sin_s = compute_cov(mol)
    c_eff = c_base * (1.0 - sin_s)
    return np.sqrt(D_cov**2 + c_eff * dE**2)


# =============================================================================
# SCAN EPSILON VALUES FOR MODEL B
# =============================================================================
print("=" * 95)
print("  SCANNING EPSILON for Model B: D = |D_cov| + dE^2/(d^2*(|D_cov| + eps))")
print("=" * 95)

eps_list = [0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
print(f"\n{'eps':>8} {'avg%':>7} {'med%':>7} {'<5%':>5} {'<10%':>5} {'<20%':>5}")
print("-" * 40)
for eps in eps_list:
    errs = []
    for mol in molecules:
        D_cov, dE, _, E_scale, _, _ = compute_cov(mol)
        Dc = abs(D_cov)
        D_pred = Dc + dE**2 / (d**2 * (Dc + eps))
        errs.append(abs((D_pred - mol[2]) / mol[2] * 100))
    print(f"{eps:8.2f} {np.mean(errs):7.1f} {np.median(errs):7.1f} "
          f"{sum(1 for e in errs if e<5):5d} {sum(1 for e in errs if e<10):5d} "
          f"{sum(1 for e in errs if e<20):5d}")


# =============================================================================
# HEAD-TO-HEAD: ALL MODELS
# =============================================================================
models = {
    'A: inv_res':    model_A,
    'D: dE^2/Escl':  model_D,
    'E: sin_wt':     model_E,
    'F: lin*deficit': model_F,
    'G: adaptive_c':  model_G,
}

print()
print("=" * 95)
print("  HEAD-TO-HEAD COMPARISON: All models vs old (cov+ionic)")
print("=" * 95)

print(f"\n{'Model':>20} {'avg%':>7} {'med%':>7} {'<5%':>5} {'<10%':>5} {'<20%':>5}")
print("-" * 55)

# Old model baseline
errs_old = [abs((old_model(mol) - mol[2]) / mol[2] * 100) for mol in molecules]
print(f"{'OLD (cov+ionic)':>20} {np.mean(errs_old):7.1f} {np.median(errs_old):7.1f} "
      f"{sum(1 for e in errs_old if e<5):5d} {sum(1 for e in errs_old if e<10):5d} "
      f"{sum(1 for e in errs_old if e<20):5d}")

for name, func in models.items():
    errs = [abs((func(mol) - mol[2]) / mol[2] * 100) for mol in molecules]
    print(f"{name:>20} {np.mean(errs):7.1f} {np.median(errs):7.1f} "
          f"{sum(1 for e in errs if e<5):5d} {sum(1 for e in errs if e<10):5d} "
          f"{sum(1 for e in errs if e<20):5d}")

# Best epsilon model B
best_eps = None
best_avg = 999
for eps in [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]:
    errs = []
    for mol in molecules:
        D_cov, dE, _, E_scale, _, _ = compute_cov(mol)
        Dc = abs(D_cov)
        D_pred = Dc + dE**2 / (d**2 * (Dc + eps))
        errs.append(abs((D_pred - mol[2]) / mol[2] * 100))
    if np.mean(errs) < best_avg:
        best_avg = np.mean(errs)
        best_eps = eps

errs_B = []
for mol in molecules:
    D_cov, dE, _, E_scale, _, _ = compute_cov(mol)
    Dc = abs(D_cov)
    D_pred = Dc + dE**2 / (d**2 * (Dc + best_eps))
    errs_B.append(abs((D_pred - mol[2]) / mol[2] * 100))
print(f"{'B: inv(eps='+str(best_eps)+')':>20} {np.mean(errs_B):7.1f} {np.median(errs_B):7.1f} "
      f"{sum(1 for e in errs_B if e<5):5d} {sum(1 for e in errs_B if e<10):5d} "
      f"{sum(1 for e in errs_B if e<20):5d}")


# =============================================================================
# DETAILED VIEW OF BEST NON-OLD MODEL
# =============================================================================
# Find best
all_results = {}
for name, func in models.items():
    errs = [abs((func(mol) - mol[2]) / mol[2] * 100) for mol in molecules]
    all_results[name] = np.mean(errs)

# Add model B
all_results[f'B: inv(eps={best_eps})'] = np.mean(errs_B)

best_name = min(all_results, key=all_results.get)

print()
print("=" * 95)
print(f"  DETAIL: Best model = {best_name}")
print("=" * 95)
print()

print(f"{'Mol':<6} {'De_exp':>7} {'D_cov':>7} {'dE':>7} {'D_new':>7} {'D_old':>7} "
      f"{'err_new%':>8} {'err_old%':>8} {'verdict':>8}")
print("-" * 80)

for mol in molecules:
    name = mol[0]
    De_exp = mol[2]
    D_cov, dE, dE_orb, E_scale, phase, sin_s = compute_cov(mol)

    if 'B:' in best_name:
        Dc = abs(D_cov)
        D_new = Dc + dE**2 / (d**2 * (Dc + best_eps))
    else:
        func = models[best_name]
        D_new = func(mol)

    D_old_val = old_model(mol)

    err_new = (D_new - De_exp) / De_exp * 100
    err_old = (D_old_val - De_exp) / De_exp * 100

    v = 'BETTER' if abs(err_new) < abs(err_old) - 1 else 'WORSE' if abs(err_new) > abs(err_old) + 1 else '~'
    flag = '***' if abs(err_new) < 2 else ' **' if abs(err_new) < 5 else '  *' if abs(err_new) < 10 else ''

    print(f"{name:<6} {De_exp:7.3f} {D_cov:7.3f} {dE:7.3f} {D_new:7.3f} {D_old_val:7.3f} "
          f"{err_new:+7.1f}% {err_old:+7.1f}% {flag} {v}")


# =============================================================================
# SCAN: c * dE * (1-|sin|) for various c
# =============================================================================
print()
print("=" * 95)
print("  SCANNING: D = |D_cov| + c * dE * (1 - |sin(phase)|)")
print("  (linear detuning weighted by sin deficit)")
print("=" * 95)

c_vals = [0.05, 0.1, 0.15, 0.2, 0.25, 1.0/(2*d), 1.0/d, 0.5, 1.0/pi]
print(f"\n{'c':>8} {'avg%':>7} {'med%':>7} {'<5%':>5} {'<10%':>5} {'<20%':>5}")
print("-" * 40)
for c_val in sorted(c_vals):
    errs = []
    for mol in molecules:
        D_cov, dE, _, E_scale, phase, sin_s = compute_cov(mol)
        D_pred = abs(D_cov) + c_val * dE * (1.0 - sin_s)
        errs.append(abs((D_pred - mol[2]) / mol[2] * 100))
    print(f"{c_val:8.4f} {np.mean(errs):7.1f} {np.median(errs):7.1f} "
          f"{sum(1 for e in errs if e<5):5d} {sum(1 for e in errs if e<10):5d} "
          f"{sum(1 for e in errs if e<20):5d}")


# =============================================================================
# DETAIL: Model F with best c
# =============================================================================
print()
print("=" * 95)
print("  DETAIL: D = |D_cov| + (1/pi) * dE * (1 - |sin(phase)|)")
print("=" * 95)
print()

c_F = 1.0/pi
print(f"{'Mol':<6} {'De_exp':>7} {'D_cov':>7} {'dE':>7} {'|sin|':>7} {'1-|sin|':>7} "
      f"{'corr':>7} {'D_pred':>7} {'err%':>7} {'err_old%':>8}")
print("-" * 85)

for mol in molecules:
    name = mol[0]
    De_exp = mol[2]
    D_cov, dE, dE_orb, E_scale, phase, sin_s = compute_cov(mol)
    correction = c_F * dE * (1.0 - sin_s)
    D_pred = abs(D_cov) + correction
    D_old_val = old_model(mol)
    err_new = (D_pred - De_exp) / De_exp * 100
    err_old = (D_old_val - De_exp) / De_exp * 100
    v = 'BETTER' if abs(err_new) < abs(err_old) - 1 else 'WORSE' if abs(err_new) > abs(err_old) + 1 else '~'
    flag = '***' if abs(err_new) < 2 else ' **' if abs(err_new) < 5 else '  *' if abs(err_new) < 10 else ''
    print(f"{name:<6} {De_exp:7.3f} {D_cov:7.3f} {dE:7.3f} {sin_s:7.4f} {1-sin_s:7.4f} "
          f"{correction:7.3f} {D_pred:7.3f} {err_new:+6.1f}% {err_old:+7.1f}% {flag} {v}")
