"""
TWO-MODE BONDING TEST
=====================
Testing the insight that bonds have two wave coupling modes:
  1. OVERLAPPING harmonics (covalent) - waves share space, sin(phase) > 0
  2. INLINE harmonics (ionic) - waves side-by-side, charge asymmetry -> Coulomb

The regime is determined by the standing wave phase:
  phase < pi: overlap mode (covalent resonance intact)
  phase > pi: node between atoms, waves separated (ionic regime)

This is the GWT interpretation: ONE wave mechanism, TWO modes.
"""
import numpy as np
import math

d = 3
pi = np.pi
alpha_gwt = d**2 / (2**(d+1) * math.factorial(d+2)**(1.0/(d+1)) * pi**((d**2+d-1)/(d+1)))
F_mode = 2*d * pi**(2*d - 1)
m_e_gwt = F_mode * alpha_gwt**12 * 1.2209e22
E_H = alpha_gwt**2 * m_e_gwt * 1e6 / 2  # GWT-derived

f_pi   = d**2 / (d**2 + 1)    # 9/10
alpha_n = 1 - f_pi / d         # 7/10
beta_n  = (1 + f_pi) / 2       # 19/20
C_bond = pi / d                 # pi/3
f_anti = 2*d / (2*d - 1)       # 6/5

Z_eff = {
    'H_1s': 1.0000, 'Li_2s': 1.2792, 'B_2p': 2.4214,
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

def get_params(orb):
    n, l, h = get_n(orb), get_l(orb), has_nodes(orb)
    a = 2 + (1 - 2*l) * alpha_n * h
    b = 1 + (1 - 2*l) * beta_n * h
    return n, l, h, a, b


def compute_old(mol):
    """Current formula: |sin(phase)| + weak ionic"""
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1, a1, b1 = get_params(orb1)
    n2, l2, h2, a2, b2 = get_params(orb2)
    E_scale = np.sqrt(E_H / n1**a1 * E_H / n2**a2)
    sigma_phase = R * (1.0 / n1**b1 + 1.0 / n2**b2)
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
    eps1 = E_H * (Z_eff[orb1] / n1)**2
    eps2 = E_H * (Z_eff[orb2] / n2)**2
    de = abs(eps1 - eps2)
    Vc = max(abs(D_cov), 0.01)
    q = de / np.sqrt(de**2 + (2 * Vc)**2) if de > 0 else 0
    Di = (1.0 / (2 * d + 1)) * q**2 * 2 * E_H / R
    return D_cov + Di, sigma_phase


def compute_smooth(mol):
    """
    Smooth two-mode model:
    1. Covalent: |sin(phase)| overlap, damped past node (phase > pi)
    2. Ionic: q-dependent coefficient interpolating from screened to unscreened
       c(q) = (d + (d+1)*q) / (d*(2d+1)) = (3+4q)/21
       At q=0: 1/7 (fully screened by covalent overlap)
       At q=1: 1/3 (no covalent screening, waves are inline)
    3. Phase damping: when phase > pi, a node exists between atoms.
       The covalent overlap weakens. Damping = exp(-(phase-pi)/pi).
    """
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1, a1, b1 = get_params(orb1)
    n2, l2, h2, a2, b2 = get_params(orb2)
    E_scale = np.sqrt(E_H / n1**a1 * E_H / n2**a2)
    sigma_phase = R * (1.0 / n1**b1 + 1.0 / n2**b2)
    eps1 = E_H * (Z_eff[orb1] / n1)**2
    eps2 = E_H * (Z_eff[orb2] / n2)**2
    de = abs(eps1 - eps2)
    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    # COVALENT: overlapping harmonics with |sin| (keeps current formula)
    D_cov = 0
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss', 'sp')) else sigma_phase * f_pi
        cont = C_bond * E_scale * abs(np.sin(ph))
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    # Phase damping: node between atoms weakens covalent overlap
    if sigma_phase > pi:
        damping = np.exp(-(sigma_phase - pi) / pi)
        D_cov = D_cov * damping

    # IONIC: charge transfer (q from 2-level model, same as current)
    Vc = max(abs(D_cov), 0.01)
    q = de / np.sqrt(de**2 + (2 * Vc)**2) if de > 0 else 0

    # q-dependent ionic coefficient: interpolate screened -> unscreened
    # c(q) = (d + (d+1)*q) / (d*(2d+1))
    # Physical: more charge transfer = less covalent screening of Coulomb
    c_ionic_q = (d + (d + 1) * q) / (d * (2 * d + 1))

    D_ionic = c_ionic_q * q**2 * 2 * E_H / R
    regime = 'COV' if sigma_phase < pi else 'ION'
    return D_cov + D_ionic, regime, sigma_phase, q, D_cov, D_ionic, c_ionic_q


# ================================================================
# COMPARE: old vs smooth two-mode
# ================================================================
print("=" * 90)
print("  SMOOTH TWO-MODE: |sin| + phase damping + q-scaled ionic coefficient")
print("  c(q) = (d + (d+1)*q) / (d*(2d+1)) = (3+4q)/21")
print("  c(q=0) = 1/7 (screened)  c(q=1) = 1/3 (unscreened)")
print("  Phase damping: D_cov *= exp(-(phase-pi)/pi) when phase > pi")
print("=" * 90)
print()

print(f"{'Mol':<6} {'De_obs':>6} {'old':>6} {'e_old':>7} {'new':>6} {'e_new':>7}"
      f" {'reg':>4} {'ph/pi':>6} {'q':>5} {'c(q)':>6} {'D_cov':>6} {'D_ion':>6}  verdict")
print('-' * 100)

old_errs = []
new_errs = []
for mol in molecules:
    De_old, phase_old = compute_old(mol)
    De_new, regime, phase_new, q, Dc, Di, cq = compute_smooth(mol)
    eo = (De_old - mol[2]) / mol[2] * 100
    en = (De_new - mol[2]) / mol[2] * 100
    old_errs.append(abs(eo))
    new_errs.append(abs(en))
    v = 'BETTER' if abs(en) < abs(eo) - 1 else ('WORSE' if abs(en) > abs(eo) + 1 else '~')
    print(f"{mol[0]:<6} {mol[2]:6.3f} {De_old:6.3f} {eo:+6.1f}% {De_new:6.3f} {en:+6.1f}%"
          f" {regime:>4} {phase_old/pi:6.3f} {q:5.3f} {cq:6.4f} {Dc:6.3f} {Di:6.3f}  {v}")

print()
print(f"OLD:  avg={np.mean(old_errs):.1f}%, med={np.median(old_errs):.1f}%, "
      f"<5%: {sum(1 for e in old_errs if e<5)}/24, "
      f"<10%: {sum(1 for e in old_errs if e<10)}/24, "
      f"<20%: {sum(1 for e in old_errs if e<20)}/24")
print(f"NEW:  avg={np.mean(new_errs):.1f}%, med={np.median(new_errs):.1f}%, "
      f"<5%: {sum(1 for e in new_errs if e<5)}/24, "
      f"<10%: {sum(1 for e in new_errs if e<10)}/24, "
      f"<20%: {sum(1 for e in new_errs if e<20)}/24")
nb = sum(1 for i in range(24) if new_errs[i] < old_errs[i] - 1)
nw = sum(1 for i in range(24) if new_errs[i] > old_errs[i] + 1)
print(f"Better: {nb}, Worse: {nw}, ~Same: {24 - nb - nw}")
print()

# ================================================================
# Also test different damping functions
# ================================================================
print("=" * 90)
print("  DAMPING FUNCTION COMPARISON")
print("=" * 90)
print()

def compute_with_damping(mol, damp_func):
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1, a1, b1 = get_params(orb1)
    n2, l2, h2, a2, b2 = get_params(orb2)
    E_scale = np.sqrt(E_H / n1**a1 * E_H / n2**a2)
    sigma_phase = R * (1.0 / n1**b1 + 1.0 / n2**b2)
    eps1 = E_H * (Z_eff[orb1] / n1)**2
    eps2 = E_H * (Z_eff[orb2] / n2)**2
    de = abs(eps1 - eps2)
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
    # Apply damping
    D_cov = D_cov * damp_func(sigma_phase)
    # q-scaled ionic
    Vc = max(abs(D_cov), 0.01)
    q = de / np.sqrt(de**2 + (2 * Vc)**2) if de > 0 else 0
    c_ionic_q = (d + (d + 1) * q) / (d * (2 * d + 1))
    D_ionic = c_ionic_q * q**2 * 2 * E_H / R
    return D_cov + D_ionic

damping_funcs = {
    'none (old formula + q-scaled c)':
        lambda ph: 1.0,
    'exp(-(ph-pi)/pi) for ph>pi':
        lambda ph: np.exp(-(ph - pi) / pi) if ph > pi else 1.0,
    'exp(-2(ph-pi)/pi) for ph>pi':
        lambda ph: np.exp(-2 * (ph - pi) / pi) if ph > pi else 1.0,
    '1/d^floor(ph/pi) [mode suppression]':
        lambda ph: (1.0 / d)**max(int(ph / pi) - 0, 0) if ph > pi else 1.0,
    'sin(pi^2/ph)/sin(pi) for ph>pi [phase folding]':
        lambda ph: abs(np.sin(pi**2 / ph)) if ph > pi else 1.0,
    'pi/ph for ph>pi [1/phase decay]':
        lambda ph: pi / ph if ph > pi else 1.0,
    '(pi/ph)^2 for ph>pi':
        lambda ph: (pi / ph)**2 if ph > pi else 1.0,
}

print(f"{'Damping':<45s} {'avg':>6} {'med':>6} {'<5%':>5} {'<10%':>5} {'<20%':>5}")
print('-' * 80)
for label, dfunc in damping_funcs.items():
    errs = []
    for mol in molecules:
        De_pred = compute_with_damping(mol, dfunc)
        errs.append(abs((De_pred - mol[2]) / mol[2] * 100))
    avg = np.mean(errs)
    med = np.median(errs)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"{label:<45s} {avg:5.1f}% {med:5.1f}% {w5:3d}/24 {w10:3d}/24 {w20:3d}/24")
