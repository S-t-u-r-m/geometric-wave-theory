"""
ONE FORCE HYPOTHESIS
====================
What if there's only ONE force (standing wave resonance)?
The "ionic" term might be compensating for errors in the covalent formula.

Key question: What does each molecule NEED from the covalent formula alone?
  D_need = De_exp  (no ionic correction)

For homonuclear molecules, the current formula already works perfectly
(ionic = 0 anyway). The "ionic" term only matters for heteronuclear.

So the question is: can we modify the covalent formula to handle
heteronuclear molecules WITHOUT a separate ionic term?

Physics to consider:
- Two standing waves with different wavelengths overlapping
- The overlap integral depends on BOTH wavevectors, not just their sum
- Electronegativity difference = different wave amplitudes at the midpoint
- This amplitude difference IS the "charge transfer" - but it's still
  wave resonance, not a separate Coulomb force
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
def bohr_radius(orb): return get_n(orb)**2 / Z_eff[orb]


# =============================================================================
# PART 1: What multiplier does each heteronuclear molecule need?
# =============================================================================
print("=" * 100)
print("  ONE FORCE: What enhancement factor does each molecule need?")
print("  If D_cov_current is the current covalent prediction (no ionic),")
print("  what factor f gives f * D_cov = De_exp?")
print("=" * 100)
print()

def compute_cov_only(mol, use_phase_fix=True):
    """Compute covalent energy only, no ionic term."""
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
    n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)

    a1 = 2 + (1 - 2*l1) * alpha_n * h1
    a2 = 2 + (1 - 2*l2) * alpha_n * h2

    if use_phase_fix:
        # Fix: no sign flip for phase
        b1 = 1 + beta_n * h1
        b2 = 1 + beta_n * h2
    else:
        b1 = 1 + (1 - 2*l1) * beta_n * h1
        b2 = 1 + (1 - 2*l2) * beta_n * h2

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

    return D_cov, sigma_phase, k1, k2, E_scale


print(f"{'Mol':<6} {'De_exp':>7} {'D_cov':>7} {'f_need':>7} {'k1':>6} {'k2':>6} "
      f"{'k2/k1':>6} {'E1':>6} {'E2':>6} {'E2/E1':>6} {'rB1':>5} {'rB2':>5}")
print("-" * 100)

enhancement_data = []
for mol in molecules:
    name, R, De_exp, bonds, orb1, orb2 = mol
    D_cov, phase, k1, k2, E_scale = compute_cov_only(mol, use_phase_fix=True)

    f_need = De_exp / D_cov if abs(D_cov) > 0.01 else float('inf')
    e1 = orb_energy(orb1)
    e2 = orb_energy(orb2)
    rB1 = bohr_radius(orb1)
    rB2 = bohr_radius(orb2)

    # Ensure e1 <= e2 for ratio
    e_ratio = max(e1, e2) / min(e1, e2) if min(e1, e2) > 0.01 else 0
    k_ratio = max(k1, k2) / min(k1, k2) if min(k1, k2) > 0.01 else 0

    is_homo = (orb1 == orb2)
    flag = '' if is_homo else ' <-- hetero'

    print(f"{name:<6} {De_exp:7.3f} {D_cov:7.3f} {f_need:7.3f} {k1:6.3f} {k2:6.3f} "
          f"{k_ratio:6.2f} {e1:6.2f} {e2:6.2f} {e_ratio:6.2f} {rB1:5.2f} {rB2:5.2f}{flag}")

    if not is_homo:
        enhancement_data.append({
            'name': name, 'f_need': f_need, 'De_exp': De_exp, 'D_cov': D_cov,
            'k1': k1, 'k2': k2, 'k_ratio': k_ratio,
            'e1': e1, 'e2': e2, 'e_ratio': e_ratio,
            'rB1': rB1, 'rB2': rB2,
            'R': R, 'phase': phase, 'E_scale': E_scale,
            'orb1': orb1, 'orb2': orb2,
        })


# =============================================================================
# PART 2: What correlates with the enhancement factor?
# =============================================================================
print()
print("=" * 100)
print("  CORRELATION: What determines the enhancement factor for heteronuclear?")
print("=" * 100)
print()

f_vals = np.array([d['f_need'] for d in enhancement_data])
names = [d['name'] for d in enhancement_data]

# Test quantities
quantities = {
    'k_ratio (max/min)': [d['k_ratio'] for d in enhancement_data],
    'E_ratio': [d['e_ratio'] for d in enhancement_data],
    'sqrt(E_ratio)': [np.sqrt(d['e_ratio']) for d in enhancement_data],
    'rB_ratio': [max(d['rB1'],d['rB2'])/min(d['rB1'],d['rB2']) for d in enhancement_data],
    'de/E_gm': [abs(d['e1']-d['e2'])/np.sqrt(d['e1']*d['e2']) for d in enhancement_data],
    'phase/pi': [d['phase']/pi for d in enhancement_data],
    '|sin(phase)|': [abs(np.sin(d['phase'])) for d in enhancement_data],
    'R': [d['R'] for d in enhancement_data],
    'R/rB_gm': [d['R']/np.sqrt(d['rB1']*d['rB2']) for d in enhancement_data],
    'E_scale': [d['E_scale'] for d in enhancement_data],
    '1/|sin(ph)|': [1/max(abs(np.sin(d['phase'])), 0.01) for d in enhancement_data],
    'k1+k2': [d['k1']+d['k2'] for d in enhancement_data],
    '2*k_gm/(k1+k2)': [2*np.sqrt(d['k1']*d['k2'])/(d['k1']+d['k2']) for d in enhancement_data],
}

print(f"{'Quantity':<25} {'corr':>7}")
print("-" * 40)

correlations = []
for qname, qvals in quantities.items():
    qv = np.array(qvals)
    if np.std(qv) > 1e-10 and np.std(f_vals) > 1e-10:
        corr = np.corrcoef(f_vals, qv)[0, 1]
    else:
        corr = 0
    correlations.append((abs(corr), qname, corr, qv))

correlations.sort(reverse=True)
for _, qname, corr, _ in correlations:
    print(f"{qname:<25} {corr:+7.3f}")

# Show the top correlator in detail
print()
print("Enhancement factors vs top correlators:")
print(f"{'Mol':<6} {'f_need':>7} {'E_ratio':>8} {'de/E_gm':>8} {'k_ratio':>8} {'ph/pi':>6}")
print("-" * 50)
for d in enhancement_data:
    de_Egm = abs(d['e1']-d['e2'])/np.sqrt(d['e1']*d['e2'])
    print(f"{d['name']:<6} {d['f_need']:7.3f} {d['e_ratio']:8.2f} {de_Egm:8.3f} "
          f"{d['k_ratio']:8.2f} {d['phase']/pi:6.3f}")


# =============================================================================
# PART 3: Physical model - overlap of two different standing waves
# =============================================================================
print()
print("=" * 100)
print("  PHYSICAL MODEL: Two standing waves with different k")
print("=" * 100)
print()
print("Current: sin(k1*R + k2*R) = sin((k1+k2)*R)")
print("  This is the phase accumulated by traveling across BOTH wavefunctions")
print()
print("Alternative: sin(k1*R) * sin(k2*R) / sin(k_ref*R)")
print("  Each wave contributes its own phase at the bond distance")
print("  For homonuclear (k1=k2=k): sin^2(kR)/sin(kR) = sin(kR) = same as current")
print()
print("Alternative 2: sin(k1*R/2) * sin(k2*R/2)")
print("  Overlap at midpoint of each wave")
print("  For homonuclear: sin^2(kR/2) != sin(2kR). Different.")
print()
print("Alternative 3: sin(sqrt(k1*k2) * R) * something")
print("  Geometric mean wavevector")
print()

# Test: what if phase uses geometric mean of k instead of sum?
# For homonuclear: 2*k*R (current) vs 2*sqrt(k*k)*R = 2*k*R (same!)
# For heteronuclear: (k1+k2)*R vs 2*sqrt(k1*k2)*R
# Difference: arithmetic mean vs geometric mean, scaled by 2R

print("Test A: phase = 2*R*sqrt(k1*k2)  (geometric mean wavevector)")
print("  For homonuclear: 2*R*k = same as current")
print("  For heteronuclear: uses geometric mean instead of arithmetic")
print()

def compute_model(mol, phase_func, label=''):
    """Generic computation with custom phase function."""
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
    sigma_phase = phase_func(R, k1, k2)

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

    return D_cov, sigma_phase


# Phase models to test (all reduce to 2kR for homonuclear)
phase_models = {
    'current: R*(k1+k2)': lambda R, k1, k2: R * (k1 + k2),
    'geom: 2R*sqrt(k1*k2)': lambda R, k1, k2: 2 * R * np.sqrt(k1 * k2),
    'harm: 4R*k1*k2/(k1+k2)': lambda R, k1, k2: 4 * R * k1 * k2 / (k1 + k2),
    'min: 2R*min(k1,k2)': lambda R, k1, k2: 2 * R * min(k1, k2),
    'max: 2R*max(k1,k2)': lambda R, k1, k2: 2 * R * max(k1, k2),
}

for pname, pfunc in phase_models.items():
    print(f"--- {pname} ---")
    errs = []
    for mol in molecules:
        D_cov, phase = compute_model(mol, pfunc)
        err = abs((D_cov - mol[2]) / mol[2] * 100)
        errs.append(err)

    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"  avg={np.mean(errs):.1f}%, med={np.median(errs):.1f}%, "
          f"<5%:{w5}/24, <10%:{w10}/24, <20%:{w20}/24")
    print()


# =============================================================================
# PART 4: Product overlap model
# =============================================================================
print("=" * 100)
print("  PRODUCT OVERLAP: D = C * E_scale * f(k1,k2,R)")
print("  Instead of sin((k1+k2)*R), what about sin(k1*R)*sin(k2*R)?")
print("=" * 100)
print()

def compute_product(mol, norm_func=None):
    """Product overlap model."""
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

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        if 'sigma' in bt or bt in ('ss', 'sp'):
            ph1 = k1 * R
            ph2 = k2 * R
        else:
            ph1 = k1 * R * f_pi
            ph2 = k2 * R * f_pi

        # Product overlap
        s1 = np.sin(ph1)
        s2 = np.sin(ph2)
        overlap = abs(s1 * s2)

        # Normalization to match homonuclear case
        if norm_func:
            overlap = norm_func(s1, s2, ph1, ph2)

        cont = C_bond * E_scale * overlap
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    return D_cov


# For homonuclear with current formula: sin(2kR)
# For product: sin(kR) * sin(kR) = sin^2(kR)
# sin(2kR) = 2*sin(kR)*cos(kR)
# So product = sin^2(kR) while current = 2*sin(kR)*cos(kR)
# These are different! Need normalization.
#
# Option: use 2*|sin(k1*R)*sin(k2*R)| / max(|sin(k1*R)|, |sin(k2*R)|)
# For homonuclear: 2*sin^2(kR)/sin(kR) = 2*sin(kR). Not right either.
#
# Better: sin(k1*R + k2*R) = sin(k1*R)*cos(k2*R) + cos(k1*R)*sin(k2*R)
# This is the angle addition formula! So the current formula IS the sum.
# The product sin(k1*R)*sin(k2*R) = [cos((k1-k2)*R) - cos((k1+k2)*R)] / 2
#
# For homonuclear: [cos(0) - cos(2kR)] / 2 = [1 - cos(2kR)] / 2 = sin^2(kR)
# Current homonuclear: sin(2kR)
#
# So sin^2(kR) vs sin(2kR) = 2*sin(kR)*cos(kR). Very different functions!
# At kR = pi/4 (H2's sweet spot ~0.89*pi/2): sin^2 = 0.5, sin(2kR) = 1.0
# The product model would need 2x the prefactor.

# Let's just try various normalizations
norms = {
    'raw |s1*s2|': lambda s1, s2, p1, p2: abs(s1 * s2),
    '2*|s1*s2|': lambda s1, s2, p1, p2: 2 * abs(s1 * s2),
    '|s1*s2|^0.5': lambda s1, s2, p1, p2: np.sqrt(abs(s1 * s2)),
    'sqrt(s1^2+s2^2)/sqrt2': lambda s1, s2, p1, p2: np.sqrt((s1**2 + s2**2) / 2),
    '|sin(2*sqrt(p1*p2))|': lambda s1, s2, p1, p2: abs(np.sin(2 * np.sqrt(p1 * p2))),
}

for nname, nfunc in norms.items():
    errs = []
    for mol in molecules:
        D_cov = compute_product(mol, norm_func=nfunc)
        err = abs((D_cov - mol[2]) / mol[2] * 100)
        errs.append(err)

    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"  {nname:<30}: avg={np.mean(errs):.1f}%, med={np.median(errs):.1f}%, "
          f"<5%:{w5}/24, <10%:{w10}/24, <20%:{w20}/24")


# =============================================================================
# PART 5: What if E_scale uses Z_eff directly?
# =============================================================================
print()
print("=" * 100)
print("  ENERGY SCALE: Current E_H/n^a vs orbital energy E_H*(Z/n)^2")
print("  What if the coupling uses the ACTUAL orbital energy?")
print("=" * 100)
print()

# Compare E_H/n^a with E_H*(Z/n)^2
print(f"{'orbital':>8} {'E_H/n^a':>8} {'E_orb':>8} {'ratio':>8}")
print("-" * 40)
for orb in sorted(Z_eff.keys()):
    n = get_n(orb)
    l = get_l(orb)
    h = has_nodes(orb)
    a = 2 + (1 - 2*l) * alpha_n * h
    E_formula = E_H / n**a
    E_orb = orb_energy(orb)
    ratio = E_orb / E_formula
    print(f"{orb:>8} {E_formula:8.3f} {E_orb:8.3f} {ratio:8.3f}")

print()
print("Testing: D = C * sqrt(E_orb1 * E_orb2) * |sin(phase)| (uses Z_eff)")
print()

def compute_with_orbital_E(mol, use_phase_fix=True):
    """Use actual orbital energies instead of E_H/n^a."""
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
    n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)

    # Use actual orbital energies
    E1 = orb_energy(orb1)
    E2 = orb_energy(orb2)
    E_scale = np.sqrt(E1 * E2)

    if use_phase_fix:
        b1 = 1 + beta_n * h1
        b2 = 1 + beta_n * h2
    else:
        b1 = 1 + (1 - 2*l1) * beta_n * h1
        b2 = 1 + (1 - 2*l2) * beta_n * h2

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

    return D_cov, sigma_phase


# Need to find the right prefactor since E_scale is now different
# For H2: E_orb = 13.6 = E_H/1^2, same as before. C = pi/3 still works.
# But for others, E_orb includes Z_eff^2 factor.
# Need C' such that C' * sqrt(E_orb1 * E_orb2) matches experiment.

# For homonuclear: C' * E_orb * |sin(phase)| = De_exp
# H2: C' * 13.6 * 0.331 = 4.745 -> C' = 1.054 = pi/3 (same!)

# For HF: E_orb(H) = 13.6, E_orb(F) = 13.6*(5.1/2)^2 = 88.56
# sqrt(13.6 * 88.56) = 34.71
# Current D_cov(HF) = 3.675 with E_scale = sqrt(13.6*13.6) = 13.6
# With E_orb scale: D_cov = 3.675 * 34.71/13.6 = 9.37! Way too much.
# So C needs to decrease, or the formula is different.

# Actually the coupling should use RELATIVE energies, not absolute.
# The resonance strength depends on how similar the two waves are.
# Use the HARMONIC mean: 2*E1*E2/(E1+E2) instead of sqrt(E1*E2)?
# Or: E_scale = E_H * sqrt(Z1*Z2) / sqrt(n1*n2)?

# Let me just scan C prefactors with orbital energies
print("Scanning C prefactor with E_orb scale:")
dd = 3  # use dd to avoid shadowing
for C_test_frac, C_label in [
    (1.0/(dd**2), '1/d^2'),
    (1.0/dd, '1/d'),
    (pi/dd**2, 'pi/d^2'),
    (1.0/(dd*(dd+1)), '1/(d*(d+1))'),
    (pi/(dd*(dd+1)), 'pi/(d*(d+1))'),
    (1.0/(dd**2+1), '1/(d^2+1)'),
    (pi/(dd**2+1), 'pi/(d^2+1)'),
]:
    errs = []
    for mol in molecules:
        name, R, De_exp, bonds, orb1, orb2 = mol
        n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
        n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)

        E1 = orb_energy(orb1)
        E2 = orb_energy(orb2)
        E_scale = np.sqrt(E1 * E2)

        b1 = 1 + beta_n * h1
        b2 = 1 + beta_n * h2
        k1 = 1.0 / n1**b1
        k2 = 1.0 / n2**b2
        sigma_phase = R * (k1 + k2)

        npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
        npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
        ifa = (npa >= npb) if npb > 0 else True

        D_cov = 0
        for bt, cnt in bonds:
            ph = sigma_phase if ('sigma' in bt or bt in ('ss', 'sp')) else sigma_phase * f_pi
            cont = C_test_frac * E_scale * abs(np.sin(ph))
            if 'anti' in bt:
                fa = 1.0 if ('sigma' in bt or ifa) else f_anti
                D_cov -= cnt * fa * cont
            else:
                D_cov += cnt * cont
        err = abs((D_cov - De_exp) / De_exp * 100)
        errs.append(err)

    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    print(f"  C={C_test_frac:.6f} ({C_label:>15s}): avg={np.mean(errs):.1f}%, "
          f"med={np.median(errs):.1f}%, <5%:{w5}/24, <10%:{w10}/24")


# =============================================================================
# PART 6: The fundamental question - what IS the overlap integral?
# =============================================================================
print()
print("=" * 100)
print("  RETHINKING: What does GWT predict for the overlap integral?")
print("=" * 100)
print()
print("In GWT, a bond is the resonance between two standing waves.")
print("The resonance energy depends on the OVERLAP of the two wave patterns.")
print()
print("For two identical atoms (k1=k2=k):")
print("  psi_1 = sin(k*r)  centered at atom 1")
print("  psi_2 = sin(k*r)  centered at atom 2")
print("  Overlap ~ sin(k*R) where R is the separation")
print("  Energy ~ C * E_scale * sin(k*R)")
print()
print("For two DIFFERENT atoms (k1 != k2):")
print("  psi_1 = sin(k1*r)  centered at atom 1")
print("  psi_2 = sin(k2*r)  centered at atom 2")
print("  Overlap at atom 2's location: sin(k1*R)")
print("  Overlap at atom 1's location: sin(k2*R)")
print("  Total overlap ~ sin(k1*R) + sin(k2*R) = 2*sin((k1+k2)*R/2)*cos((k1-k2)*R/2)")
print()
print("This is the SUM formula! Let me test it:")
print("  D = C * E_scale * |sin(k1*R) + sin(k2*R)| / 2")
print("  = C * E_scale * |sin((k1+k2)*R/2) * cos((k1-k2)*R/2)|")
print("  For homonuclear: |sin(kR) * cos(0)| = |sin(kR)|")
print("  Current formula: |sin(2kR)| = |2*sin(kR)*cos(kR)|")
print("  Different by factor of 2*cos(kR)!")
print()

# The sum formula: sin(k1*R) + sin(k2*R) = 2*sin((k1+k2)*R/2)*cos((k1-k2)*R/2)
# For homonuclear: 2*sin(kR)*1 = 2*sin(kR)
# Current: sin(2kR) = 2*sin(kR)*cos(kR)
# Ratio: 1/cos(kR) ... hmm

# Actually wait - for homonuclear with current formula:
# phase = R*(k+k) = 2kR, so sin(2kR)
# For the sum model: sin(kR) + sin(kR) = 2*sin(kR)
# These are different: sin(2kR) = 2*sin(kR)*cos(kR) vs 2*sin(kR)
# For H2: kR = 1.401, sin(2*1.401) = sin(2.802) = 0.331
#   vs 2*sin(1.401) = 2*0.9855 = 1.971
# Very different!

# So the "sum of individual overlaps" model gives a completely different
# functional form. Let me test it with the right prefactor.

# For H2: D = C' * 13.6 * 2*sin(1.401) = C' * 13.6 * 1.971 = C' * 26.81
# Need D = 4.745 -> C' = 0.177 = close to 1/(2*pi) = 0.159 or 1/d! = 1/6 = 0.167

print("Testing sum-of-overlaps: D = C * E_scale * |sin(k1*R) + sin(k2*R)|")
print()

for C_test, C_label in [
    (1/(2*pi), '1/(2*pi)'),
    (1.0/6, '1/d!'),
    (1.0/(3+1), '1/(d+1)'),
    (1.0/(2*3), '1/(2d)'),
    (1.0/(3**2-3), '1/(d^2-d)'),
    (pi/(3**2*(3+1)), 'pi/(d^2*(d+1))'),
]:
    errs = []
    for mol in molecules:
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

        npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
        npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
        ifa = (npa >= npb) if npb > 0 else True

        D_cov = 0
        for bt, cnt in bonds:
            if 'sigma' in bt or bt in ('ss', 'sp'):
                ph1 = k1 * R
                ph2 = k2 * R
            else:
                ph1 = k1 * R * f_pi
                ph2 = k2 * R * f_pi

            overlap = abs(np.sin(ph1) + np.sin(ph2))
            cont = C_test * E_scale * overlap
            if 'anti' in bt:
                fa = 1.0 if ('sigma' in bt or ifa) else f_anti
                D_cov -= cnt * fa * cont
            else:
                D_cov += cnt * cont

        err = abs((D_cov - De_exp) / De_exp * 100)
        errs.append(err)

    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"  C={C_test:.6f} ({C_label:>18s}): avg={np.mean(errs):.1f}%, "
          f"med={np.median(errs):.1f}%, <5%:{w5}/24, <10%:{w10}/24, <20%:{w20}/24")


# Best model detail
print()
print("Detail for best sum-of-overlaps model:")
C_best = 1.0 / 6  # try 1/(2d) = 1/6

print(f"{'Mol':<6} {'De_exp':>7} {'D_pred':>7} {'err%':>7} {'sin1':>6} {'sin2':>6} {'sum':>6}")
print("-" * 55)

for mol in molecules:
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

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    ph1_main = k1 * R
    ph2_main = k2 * R
    for bt, cnt in bonds:
        if 'sigma' in bt or bt in ('ss', 'sp'):
            ph1 = k1 * R
            ph2 = k2 * R
        else:
            ph1 = k1 * R * f_pi
            ph2 = k2 * R * f_pi

        overlap = abs(np.sin(ph1) + np.sin(ph2))
        cont = C_best * E_scale * overlap
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    err = (D_cov - De_exp) / De_exp * 100
    s1 = np.sin(ph1_main)
    s2 = np.sin(ph2_main)
    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
    print(f"{name:<6} {De_exp:7.3f} {D_cov:7.3f} {err:+6.1f}% {s1:6.3f} {s2:6.3f} "
          f"{abs(s1+s2):6.3f} {flag}")
