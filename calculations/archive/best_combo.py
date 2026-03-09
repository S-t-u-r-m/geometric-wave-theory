"""
Best Combination Test: Node counting + Tuned enhanced ionic
============================================================
Node counting: |sin(phase)| / ceil(phase/pi) when has_nodes AND phase > pi
Enhanced ionic: c = 3/7 when D_cov/delta_eps < threshold (truly ionic limit)

Key finding from separate_functions.py:
  thr=0.03 gives avg=7.4%, w5=17, w10=20 (vs baseline 13.7%, w5=15, w10=17)

This script finds the exact D_cov/delta_eps ratios for all molecules and
identifies what d=3 threshold separates ionic from covalent bonds.
"""

import numpy as np

pi = np.pi
E_H = 13.6057
d = 3

C_bond  = pi / d
f_pi    = d**2 / (d**2 + 1)
alpha   = 1 - f_pi / d
beta    = (1 + f_pi) / 2
f_anti  = 2*d / (2*d - 1)
c_ionic_base = 1.0 / (2*d + 1)   # 1/7
c_ionic_full = float(d) / (2*d + 1)  # 3/7

Z_eff = {
    'H_1s': 1.0000, 'Li_2s': 1.2792, 'B_2p': 2.4214,
    'C_2p': 3.1358, 'N_2p': 3.8340, 'O_2p': 4.4532,
    'F_2p': 5.0998, 'Na_3s': 2.5074, 'Cl_3p': 6.1161,
}

def get_n(orb): return int(orb.split('_')[1][0])
def get_l(orb): return {'s': 0, 'p': 1}[orb.split('_')[1][1]]
def has_nodes(orb): return min(get_n(orb) - get_l(orb) - 1, 1)
def orbital_energy(orb): return E_H * (Z_eff[orb] / get_n(orb))**2

molecules = [
    ('H2',   1.401,  4.745, [('ss',1)],                                         'H_1s',  'H_1s'),
    ('Li2',  5.051,  1.056, [('ss',1)],                                         'Li_2s', 'Li_2s'),
    ('B2',   3.005,  3.02,  [('pi',2)],                                         'B_2p',  'B_2p'),
    ('C2',   2.348,  6.32,  [('pi',2),('sp_sigma',1),('sp_sigma_anti',1)],      'C_2p',  'C_2p'),
    ('N2',   2.074,  9.759, [('pp_sigma',1),('pi',2)],                          'N_2p',  'N_2p'),
    ('O2',   2.282,  5.213, [('pp_sigma',1),('pi',2),('pi_anti',1)],            'O_2p',  'O_2p'),
    ('F2',   2.668,  1.660, [('pp_sigma',1),('pi',2),('pi_anti',2)],            'F_2p',  'F_2p'),
    ('Na2',  5.818,  0.746, [('ss',1)],                                         'Na_3s', 'Na_3s'),
    ('Cl2',  3.757,  2.514, [('pp_sigma',1),('pi',2),('pi_anti',2)],            'Cl_3p', 'Cl_3p'),
    ('HF',   1.733,  5.869, [('sp',1)],                                         'H_1s',  'F_2p'),
    ('CO',   2.132, 11.225, [('pp_sigma',1),('pi',2)],                          'C_2p',  'O_2p'),
    ('NO',   2.175,  6.497, [('pp_sigma',1),('pi',2),('pi_anti',1)],            'N_2p',  'O_2p'),
    ('OH',   1.834,  4.392, [('sp',1)],                                         'O_2p',  'H_1s'),
    ('HCl',  2.409,  4.434, [('sp',1)],                                         'H_1s',  'Cl_3p'),
    ('LiH',  3.015,  2.515, [('ss',1)],                                         'Li_2s', 'H_1s'),
    ('LiF',  2.955,  5.939, [('sp',1)],                                         'Li_2s', 'F_2p'),
    ('BH',   2.329,  3.42,  [('sp',1)],                                         'B_2p',  'H_1s'),
    ('CH',   2.116,  3.47,  [('sp',1)],                                         'C_2p',  'H_1s'),
    ('NH',   1.958,  3.57,  [('sp',1)],                                         'N_2p',  'H_1s'),
    ('BF',   2.386,  7.81,  [('pp_sigma',1),('pi',2)],                          'B_2p',  'F_2p'),
    ('CN',   2.214,  7.72,  [('pp_sigma',0.5),('pi',2)],                        'C_2p',  'N_2p'),
    ('NaH',  3.566,  1.97,  [('ss',1)],                                         'Na_3s', 'H_1s'),
    ('NaCl', 4.461,  4.23,  [('sp',1)],                                         'Na_3s', 'Cl_3p'),
    ('H2O',  1.809,  5.117, [('sp',1)],                                         'O_2p',  'H_1s'),
]


def compute(mol, threshold=None):
    """Compute with node counting + optional enhanced ionic."""
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
    n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)

    a1 = 2 + (1-2*l1)*alpha*h1; a2 = 2 + (1-2*l2)*alpha*h2
    b1 = 1 + beta*h1; b2 = 1 + beta*h2

    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)
    sigma_phase = R/n1**b1 + R/n2**b2
    either_nodes = (h1 > 0 or h2 > 0)

    n_pb = sum(c for bt,c in bonds if 'pi' in bt and 'anti' not in bt)
    n_pa = sum(c for bt,c in bonds if 'pi' in bt and 'anti' in bt)
    is_full = (n_pa >= n_pb) if n_pb > 0 else True

    Dc = 0
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss','sp')) else sigma_phase * f_pi
        ov = abs(np.sin(ph))
        # Node counting
        if either_nodes and ph > pi:
            ov /= int(np.ceil(ph / pi))
        ct = C_bond * E_scale * ov
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or is_full) else f_anti
            Dc -= cnt * fa * ct
        else:
            Dc += cnt * ct

    eps1, eps2 = orbital_energy(orb1), orbital_energy(orb2)
    delta = abs(eps1 - eps2)
    V = max(abs(Dc), 0.01)
    q = delta / np.sqrt(delta**2 + (2*V)**2) if delta > 0 else 0

    ratio = abs(Dc) / delta if delta > 0 else float('inf')

    # Enhanced ionic
    enhanced = False
    if threshold is not None and q > 0.99 and ratio < threshold:
        c_ion = c_ionic_full
        enhanced = True
    else:
        c_ion = c_ionic_base

    D_ion = c_ion * q**2 * 2 * E_H / R
    return Dc + D_ion, Dc, D_ion, q, ratio, enhanced, sigma_phase


# =============================================================================
# SHOW ALL D_cov/delta_eps RATIOS (sorted) to find natural threshold
# =============================================================================
print("=" * 80)
print("  D_cov / delta_eps RATIO FOR ALL HETERONUCLEAR MOLECULES")
print("  (sorted by ratio -- this is the ionicity measure)")
print("=" * 80)

het_data = []
for mol in molecules:
    name, R, De_exp, bonds, orb1, orb2 = mol
    if orb1 == orb2:
        continue
    D_tot, Dc, D_ion, q, ratio, _, phase = compute(mol)
    het_data.append((ratio, name, Dc, orbital_energy(orb1), orbital_energy(orb2), q, phase, De_exp, D_tot))

het_data.sort()

print(f"\n  {'Mol':<7} {'D_cov':>7} {'delta':>7} {'ratio':>8} {'q':>6} {'ph/pi':>6} {'type':>12}")
print("  " + "-" * 60)
for ratio, name, Dc, e1, e2, q, ph, De_exp, D_tot in het_data:
    delta = abs(e1 - e2)
    typ = "IONIC" if ratio < 0.04 else "covalent"
    print(f"  {name:<7} {Dc:7.3f} {delta:7.1f} {ratio:8.4f} {q:6.3f} {ph/pi:6.3f} {typ:>12}")

print(f"\n  Gap analysis:")
print(f"    NaCl ratio = {het_data[2][0]:.4f}")
print(f"    HF   ratio = {het_data[3][0]:.4f}  <-- threshold must be between these")


# =============================================================================
# d=3 THRESHOLD CANDIDATES
# =============================================================================
print(f"\n{'='*80}")
print(f"  d=3 THRESHOLD CANDIDATES")
print(f"{'='*80}")

candidates = {
    '1/(2d(d+1))=1/24':    1/(2*d*(d+1)),
    '1/(d^2+d)=1/12':      1/(d**2 + d),
    '1/(d(2d+1))=1/21':    1/(d*(2*d+1)),
    '1/d^3=1/27':          1/d**3,
    'pi/100':              pi/100,
    '1/(4d^2)=1/36':       1/(4*d**2),
    '0.03':                0.03,
    '0.035':               0.035,
    '0.04':                0.04,
}

for cname, cval in sorted(candidates.items(), key=lambda x: x[1]):
    errs = []
    enhanced_list = []
    for mol in molecules:
        D_tot, Dc, D_ion, q, ratio, enhanced, _ = compute(mol, threshold=cval)
        err = (D_tot - mol[2]) / mol[2] * 100
        errs.append(abs(err))
        if enhanced:
            enhanced_list.append(mol[0])

    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    print(f"  {cname:<25} = {cval:.5f}: avg={np.mean(errs):5.1f}%, w5={w5:2d}, w10={w10:2d}  [{', '.join(enhanced_list)}]")


# =============================================================================
# DETAILED RESULTS: Best threshold
# =============================================================================
print(f"\n{'='*80}")
print(f"  DETAILED: Node counting + Enhanced ionic (thr=1/d^3=1/27)")
print(f"{'='*80}")

threshold = 1/d**3  # 1/27 = 0.0370

print(f"\n  Node counting: |sin(phase)|/ceil(phase/pi) when has_nodes & phase > pi")
print(f"  Enhanced ionic: c = d/(2d+1) = 3/7 when q > 0.99 AND |D_cov|/delta < 1/d^3")
print(f"  Threshold = 1/{d}^3 = {threshold:.5f}")

header = f"{'Mol':<7} {'De_exp':>7} {'De_pred':>7} {'err%':>7}  {'D_cov':>7} {'D_ion':>6} {'q':>5} {'c':>5} {'ph/pi':>6} {'node':>4}"
print(f"\n{header}")
print("-" * 80)

errs = []
for mol in molecules:
    name, R, De_exp, bonds, orb1, orb2 = mol
    D_tot, Dc, D_ion, q, ratio, enhanced, phase = compute(mol, threshold=threshold)
    err = (D_tot - De_exp) / De_exp * 100
    errs.append(abs(err))

    h1, h2 = has_nodes(orb1), has_nodes(orb2)
    either_nodes = h1 > 0 or h2 > 0
    node_active = either_nodes and phase > pi
    c_str = "3/7" if enhanced else "1/7"
    node_str = "YES" if node_active else ""
    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''

    print(f"{name:<7} {De_exp:7.3f} {D_tot:7.3f} {err:+6.1f}%  "
          f"{Dc:7.3f} {D_ion:6.3f} {q:5.3f} {c_str:>5} {phase/pi:6.3f} {node_str:>4} {flag}")

print("-" * 80)
n = len(errs)
w2 = sum(1 for e in errs if e < 2)
w5 = sum(1 for e in errs if e < 5)
w10 = sum(1 for e in errs if e < 10)
print(f"  avg = {np.mean(errs):.1f}%, median = {np.median(errs):.1f}%")
print(f"  within 2%: {w2}/{n},  within 5%: {w5}/{n},  within 10%: {w10}/{n}")

# Count outliers
outliers = [(i, mol[0], errs[i]) for i, mol in enumerate(molecules) if errs[i] > 10]
print(f"\n  Remaining outliers (>10%):")
for i, name, e in outliers:
    print(f"    {name}: {e:.1f}%")


# =============================================================================
# TRY: Resonance symmetry ONLY for BF/CN (asymmetric triple bonds)
# =============================================================================
print(f"\n{'='*80}")
print(f"  ADDING RESONANCE SYMMETRY for asymmetric bonds")
print(f"{'='*80}")

for sym_power in [0.25, 0.5, 0.75, 1.0, 1.5, 2.0]:
    errs2 = []
    for mol in molecules:
        name, R, De_exp, bonds, orb1, orb2 = mol
        n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
        n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)

        a1 = 2+(1-2*l1)*alpha*h1; a2 = 2+(1-2*l2)*alpha*h2
        b1 = 1+beta*h1; b2 = 1+beta*h2
        E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)
        sigma_phase = R/n1**b1 + R/n2**b2
        either_nodes = (h1>0 or h2>0)

        # Resonance symmetry
        eps1, eps2 = orbital_energy(orb1), orbital_energy(orb2)
        if eps1 > 0 and eps2 > 0 and abs(eps1-eps2) > 0.01:
            sym = (2*np.sqrt(eps1*eps2)/(eps1+eps2)) ** sym_power
        else:
            sym = 1.0

        n_pb = sum(c for bt,c in bonds if 'pi' in bt and 'anti' not in bt)
        n_pa = sum(c for bt,c in bonds if 'pi' in bt and 'anti' in bt)
        is_full = (n_pa >= n_pb) if n_pb > 0 else True

        Dc = 0
        for bt, cnt in bonds:
            ph = sigma_phase if ('sigma' in bt or bt in ('ss','sp')) else sigma_phase*f_pi
            ov = abs(np.sin(ph))
            if either_nodes and ph > pi:
                ov /= int(np.ceil(ph/pi))
            ct = C_bond * E_scale * ov * sym
            if 'anti' in bt:
                fa = 1.0 if ('sigma' in bt or is_full) else f_anti
                Dc -= cnt*fa*ct
            else:
                Dc += cnt*ct

        delta = abs(eps1-eps2)
        V = max(abs(Dc), 0.01)
        q = delta/np.sqrt(delta**2+(2*V)**2) if delta > 0 else 0
        ratio = abs(Dc)/delta if delta > 0 else float('inf')

        if q > 0.99 and ratio < threshold:
            c_ion = c_ionic_full
        else:
            c_ion = c_ionic_base

        D_tot = Dc + c_ion*q**2*2*E_H/R
        errs2.append(abs((D_tot - De_exp)/De_exp*100))

    w5 = sum(1 for e in errs2 if e < 5)
    w10 = sum(1 for e in errs2 if e < 10)
    # Show BF and CN specifically
    bf_err = errs2[19]
    cn_err = errs2[20]
    co_err = errs2[10]
    print(f"  sym^{sym_power:.2f}: avg={np.mean(errs2):5.1f}%, w5={w5:2d}, w10={w10:2d}  "
          f"BF={bf_err:.1f}%, CN={cn_err:.1f}%, CO={co_err:.1f}%")


# =============================================================================
# FINAL BEST: Show what we have vs what remains
# =============================================================================
print(f"\n{'='*80}")
print(f"  SCORECARD: Baseline vs Best combo")
print(f"{'='*80}")

base_errs = []
best_errs = []
for mol in molecules:
    # Baseline (no corrections)
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
    n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)
    a1 = 2+(1-2*l1)*alpha*h1; a2 = 2+(1-2*l2)*alpha*h2
    b1 = 1+beta*h1; b2 = 1+beta*h2
    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)
    sp = R/n1**b1 + R/n2**b2
    n_pb = sum(c for bt,c in bonds if 'pi' in bt and 'anti' not in bt)
    n_pa = sum(c for bt,c in bonds if 'pi' in bt and 'anti' in bt)
    is_full = (n_pa >= n_pb) if n_pb > 0 else True
    Dc0 = 0
    for bt, cnt in bonds:
        ph = sp if ('sigma' in bt or bt in ('ss','sp')) else sp*f_pi
        ct = C_bond*E_scale*abs(np.sin(ph))
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or is_full) else f_anti
            Dc0 -= cnt*fa*ct
        else:
            Dc0 += cnt*ct
    eps1,eps2 = orbital_energy(orb1),orbital_energy(orb2)
    delta=abs(eps1-eps2); V=max(abs(Dc0),0.01)
    q0=delta/np.sqrt(delta**2+(2*V)**2) if delta>0 else 0
    D0 = Dc0 + c_ionic_base*q0**2*2*E_H/R
    base_errs.append(abs((D0-De_exp)/De_exp*100))

    # Best combo
    D_tot, _, _, _, _, _, _ = compute(mol, threshold=threshold)
    best_errs.append(abs((D_tot-De_exp)/De_exp*100))

print(f"\n  {'Mol':<7} {'baseline':>10} {'best':>10} {'delta':>8}")
print("  " + "-" * 40)
for i, mol in enumerate(molecules):
    diff = best_errs[i] - base_errs[i]
    marker = " <<<" if diff < -5 else ""
    print(f"  {mol[0]:<7} {base_errs[i]:9.1f}% {best_errs[i]:9.1f}% {diff:+7.1f}%{marker}")

print(f"\n  {'':15} {'Baseline':>10} {'Best':>10}")
print(f"  {'avg':15} {np.mean(base_errs):9.1f}% {np.mean(best_errs):9.1f}%")
print(f"  {'median':15} {np.median(base_errs):9.1f}% {np.median(best_errs):9.1f}%")
print(f"  {'w2':15} {sum(1 for e in base_errs if e<2):>9} {sum(1 for e in best_errs if e<2):>9}")
print(f"  {'w5':15} {sum(1 for e in base_errs if e<5):>9} {sum(1 for e in best_errs if e<5):>9}")
print(f"  {'w10':15} {sum(1 for e in base_errs if e<10):>9} {sum(1 for e in best_errs if e<10):>9}")
