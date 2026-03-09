"""
IONIC COUPLING ANALYSIS
========================
Why does c_ion vary across ION-regime molecules?

For each molecule in the ION regime (phase > pi, heteronuclear),
compute the c_ion that gives exact agreement with experiment,
then look for correlations with orbital properties.

Physical hypothesis: The 1/(2d+1) = 1/7 coupling assumes point charges.
Real orbitals have spatial extent. The effective Coulomb coupling should
depend on the orbital geometry relative to the bond axis.
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
c_ionic_cov = 1.0 / (2*d + 1)

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
def bohr_radius(orb): return get_n(orb)**2 / Z_eff[orb]
def orb_energy(orb): return E_H * (Z_eff[orb] / get_n(orb))**2

def get_params(orb):
    n, l, h = get_n(orb), get_l(orb), has_nodes(orb)
    return n, 2 + (1-2*l)*alpha_n*h, 1 + (1-2*l)*beta_n*h


def compute_base(mol):
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, a1, b1 = get_params(orb1)
    n2, a2, b2 = get_params(orb2)
    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)
    sigma_phase = R * (1.0/n1**b1 + 1.0/n2**b2)
    eps1 = E_H * (Z_eff[orb1]/n1)**2
    eps2 = E_H * (Z_eff[orb2]/n2)**2
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

    return R, D_cov, sigma_phase, de, E_scale


# For each molecule, compute what c_ion would give exact agreement
print("=" * 100)
print("  IONIC COUPLING ANALYSIS: What c_ion does each molecule need?")
print("=" * 100)
print()

# First: which molecules are in COV vs ION regime?
print(f"{'Mol':<6} {'regime':>5} {'ph/pi':>6} {'De_exp':>7} {'D_cov':>7} {'D_need':>7} "
      f"{'q':>5} {'c_need':>7} {'r_B1':>5} {'r_B2':>5} {'R/r_gm':>6} {'E_rat':>6}")
print("-" * 100)

ion_data = []
for mol in molecules:
    name, R, De_exp, bonds, orb1, orb2 = mol
    R_val, D_cov, sigma_phase, de, E_scale = compute_base(mol)

    Vc = max(abs(D_cov), 0.01)
    q = de / np.sqrt(de**2 + (2*Vc)**2) if de > 0 else 0

    regime = 'COV' if sigma_phase < pi else 'ION'
    is_homo = (de < 0.01)
    if is_homo and regime == 'ION':
        regime = 'HOM'  # homonuclear in ION phase range

    # What D_ionic is needed?
    D_need = De_exp - D_cov

    # What c_ion gives this?
    coulomb_base = q**2 * 2 * E_H / R_val
    c_need = D_need / coulomb_base if coulomb_base > 0.001 else float('inf')

    # Orbital properties
    r_B1 = bohr_radius(orb1)
    r_B2 = bohr_radius(orb2)
    r_gm = np.sqrt(r_B1 * r_B2)  # geometric mean Bohr radius
    e1 = orb_energy(orb1)
    e2 = orb_energy(orb2)
    E_ratio = max(e1, e2) / min(e1, e2) if min(e1, e2) > 0 else 0

    flag = ''
    if regime == 'ION' and not is_homo:
        flag = ' <--'
        ion_data.append({
            'name': name, 'R': R_val, 'De_exp': De_exp, 'D_cov': D_cov,
            'D_need': D_need, 'q': q, 'c_need': c_need,
            'r_B1': r_B1, 'r_B2': r_B2, 'r_gm': r_gm,
            'E_ratio': E_ratio, 'phase': sigma_phase,
            'de': de, 'orb1': orb1, 'orb2': orb2,
            'n1': get_n(orb1), 'n2': get_n(orb2),
            'l1': get_l(orb1), 'l2': get_l(orb2),
        })

    c_str = f"{c_need:7.3f}" if abs(c_need) < 100 else "   inf"
    print(f"{name:<6} {regime:>5} {sigma_phase/pi:6.3f} {De_exp:7.3f} {D_cov:7.3f} {D_need:7.3f} "
          f"{q:5.3f} {c_str} {r_B1:5.3f} {r_B2:5.3f} {R_val/r_gm:6.3f} {E_ratio:6.2f}{flag}")


# Now focus on ION molecules and look for patterns
print()
print("=" * 100)
print("  ION-REGIME MOLECULES: Searching for c_ion pattern")
print("=" * 100)
print()

if ion_data:
    print("Correlation analysis:")
    print()

    c_vals = np.array([d['c_need'] for d in ion_data])
    names = [d['name'] for d in ion_data]

    # Test various quantities
    quantities = {
        'R/r_gm': [d['R'] / d['r_gm'] for d in ion_data],
        'r_gm/R': [d['r_gm'] / d['R'] for d in ion_data],
        'r_B1/R': [d['r_B1'] / d['R'] for d in ion_data],
        'r_B2/R': [d['r_B2'] / d['R'] for d in ion_data],
        'max(r_B)/R': [max(d['r_B1'], d['r_B2']) / d['R'] for d in ion_data],
        'min(r_B)/R': [min(d['r_B1'], d['r_B2']) / d['R'] for d in ion_data],
        'r_B1*r_B2/R^2': [d['r_B1'] * d['r_B2'] / d['R']**2 for d in ion_data],
        'E_ratio': [d['E_ratio'] for d in ion_data],
        '1/E_ratio': [1/d['E_ratio'] for d in ion_data],
        'q': [d['q'] for d in ion_data],
        'q^2': [d['q']**2 for d in ion_data],
        'de/E_H': [d['de'] / E_H for d in ion_data],
        'phase/pi - 1': [d['phase']/pi - 1 for d in ion_data],
        'n1+n2': [d['n1'] + d['n2'] for d in ion_data],
        'n1*n2': [d['n1'] * d['n2'] for d in ion_data],
        'l1+l2': [d['l1'] + d['l2'] for d in ion_data],
        '|l1-l2|': [abs(d['l1'] - d['l2']) for d in ion_data],
        'D_cov/De_exp': [d['D_cov'] / d['De_exp'] for d in ion_data],
        'D_need/D_cov': [d['D_need'] / d['D_cov'] if abs(d['D_cov']) > 0.01 else 0 for d in ion_data],
    }

    print(f"{'Quantity':<20} {'corr':>7}  values per molecule")
    print("-" * 80)

    correlations = []
    for qname, qvals in quantities.items():
        qv = np.array(qvals)
        if np.std(qv) > 1e-10 and np.std(c_vals) > 1e-10:
            corr = np.corrcoef(c_vals, qv)[0, 1]
        else:
            corr = 0
        correlations.append((abs(corr), qname, corr, qv))

    correlations.sort(reverse=True)
    for _, qname, corr, qv in correlations:
        vals_str = ", ".join(f"{v:.3f}" for v in qv)
        print(f"{qname:<20} {corr:+7.3f}  [{vals_str}]")

    print()
    print("c_need values:      [" + ", ".join(f"{c:.3f}" for c in c_vals) + "]")
    print("molecules:          [" + ", ".join(names) + "]")

    # Top 3 correlations — try fitting
    print()
    print("=" * 100)
    print("  LINEAR FITS: c_ion = a * X + b")
    print("=" * 100)
    print()

    for rank, (acorr, qname, corr, qv) in enumerate(correlations[:5]):
        if len(qv) < 2:
            continue
        # Linear fit
        coeffs = np.polyfit(qv, c_vals, 1)
        fitted = np.polyval(coeffs, qv)
        residuals = c_vals - fitted
        rmse = np.sqrt(np.mean(residuals**2))

        print(f"Rank {rank+1}: c_ion = {coeffs[0]:+.4f} * {qname} + {coeffs[1]:+.4f}  "
              f"(corr={corr:+.3f}, RMSE={rmse:.3f})")
        for i, d in enumerate(ion_data):
            print(f"  {d['name']:<6}: c_need={c_vals[i]:+.3f}, c_fit={fitted[i]:+.3f}, "
                  f"X={qv[i]:.3f}, resid={residuals[i]:+.3f}")
        print()

    # Physical model: Coulomb screening by orbital size
    # Point charge approximation breaks when orbital extends beyond R
    # Screening factor ~ 1 - (r_orbital / R) for r < R
    print("=" * 100)
    print("  PHYSICAL MODEL: Orbital size screening")
    print("=" * 100)
    print()
    print("Hypothesis: c_ion = c_base * f(r_B/R)")
    print("  When orbitals are small vs R: point-charge limit, c large")
    print("  When orbitals are large vs R: charge spread out, c small")
    print()

    for d_item in ion_data:
        r1, r2, R = d_item['r_B1'], d_item['r_B2'], d_item['R']
        # Various screening models
        s1 = 1 - r1/R  # simple linear screening
        s2 = 1 - r2/R
        s_avg = 1 - (r1+r2)/(2*R)
        s_prod = (1 - r1/R) * (1 - r2/R)
        s_gm = 1 - np.sqrt(r1*r2)/R

        print(f"  {d_item['name']:<6}: r1={r1:.3f}, r2={r2:.3f}, R={R:.3f}  "
              f"1-r1/R={s1:+.3f}, 1-r2/R={s2:+.3f}, "
              f"1-r_gm/R={s_gm:+.3f}, c_need={d_item['c_need']:+.3f}")


# Also check COV molecules: do they all work with c=1/7?
print()
print("=" * 100)
print("  COV-REGIME CHECK: How well does c = 1/7 work?")
print("=" * 100)
print()

print(f"{'Mol':<6} {'De_exp':>7} {'D_cov':>7} {'D_ion':>7} {'D_tot':>7} {'err%':>7}")
print("-" * 50)

for mol in molecules:
    name, R, De_exp, bonds, orb1, orb2 = mol
    R_val, D_cov, sigma_phase, de, E_scale = compute_base(mol)

    if sigma_phase >= pi:
        continue

    Vc = max(abs(D_cov), 0.01)
    q = de / np.sqrt(de**2 + (2*Vc)**2) if de > 0 else 0
    D_ion = c_ionic_cov * q**2 * 2 * E_H / R_val
    D_tot = D_cov + D_ion
    err = (D_tot - De_exp) / De_exp * 100
    print(f"{name:<6} {De_exp:7.3f} {D_cov:7.3f} {D_ion:7.3f} {D_tot:7.3f} {err:+6.1f}%")
