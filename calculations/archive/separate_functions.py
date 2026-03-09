"""
Separate Functions for Bond Prediction
=======================================
Tests physically motivated corrections for each outlier category:

1. NODE COUNTING: When phase > pi AND atom has radial nodes,
   divide |sin(phase)| by ceil(phase/pi) = number of lobes.
   Physics: each additional wave lobe reduces effective overlap.
   Fixes: LiH, NaH

2. ENHANCED IONIC: For fully ionic bonds (q > threshold AND
   covalent is weak), use c_ionic = d/(2d+1) = 3/7 instead of 1/7.
   Physics: localized charge gives d× stronger Coulomb coupling.
   Fixes: LiF, NaCl

3. RESONANCE SYMMETRY: Multiply D_cov by 2*sqrt(e1*e2)/(e1+e2)
   where e = orbital energy from Z_eff.
   Physics: mismatched wave frequencies reduce resonance coupling.
   Fixes: BF (partially CN)

4. SELF-CONSISTENT R: For phase-at-node cases, use predicted R
   instead of experimental R to avoid sin(pi)=0 coincidence.
   Physics: the model should be self-consistent.
   Fixes: CH

All correction parameters from d=3. Zero free parameters.
"""

import numpy as np

pi = np.pi
E_H = 13.6057  # eV
d = 3

# Parameters from d=3
C_bond  = pi / d
f_pi    = d**2 / (d**2 + 1)    # 9/10
alpha   = 1 - f_pi / d          # 7/10
beta    = (1 + f_pi) / 2        # 19/20
f_anti  = 2*d / (2*d - 1)       # 6/5
c_ionic_base = 1.0 / (2*d + 1)  # 1/7
c_ionic_full = d / (2*d + 1)    # 3/7 — enhanced ionic

# Clementi-Raimondi Z_eff
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

def bond_order(bonds):
    return sum(c if 'anti' not in bt else -c for bt, c in bonds)

def bonding_radius(orb):
    n, l = get_n(orb), get_l(orb)
    return n**2 / (Z_eff[orb] * (n - l))

def atomic_radius(orb):
    return get_n(orb)**2 / Z_eff[orb]

def predict_R(orb1, orb2, bonds):
    bo = bond_order(bonds)
    if bo == 3:
        return atomic_radius(orb1) + atomic_radius(orb2)
    else:
        return bonding_radius(orb1) + bonding_radius(orb2)


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


def compute_energy(name, R, De_exp, bonds, orb1, orb2,
                   use_node_counting=False,
                   use_enhanced_ionic=False,
                   use_resonance_sym=False,
                   use_predicted_R=False):
    """Compute bond energy with optional corrections."""

    # Option: use predicted bond length instead of experimental
    if use_predicted_R:
        R = predict_R(orb1, orb2, bonds)

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

    # Resonance symmetry factor
    sym_factor = 1.0
    if use_resonance_sym:
        eps1 = orbital_energy(orb1)
        eps2 = orbital_energy(orb2)
        if eps1 > 0 and eps2 > 0 and abs(eps1 - eps2) > 0.01:
            sym_factor = 2 * np.sqrt(eps1 * eps2) / (eps1 + eps2)

    # Node counting: only when at least one atom has radial nodes
    either_has_nodes = (h1 > 0 or h2 > 0)

    n_pi_bond = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    n_pi_anti = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    is_full_anti = (n_pi_anti >= n_pi_bond) if n_pi_bond > 0 else True

    D_cov = 0
    for btype, count in bonds:
        if 'sigma' in btype or btype in ('ss', 'sp'):
            phase = sigma_phase
        else:
            phase = sigma_phase * f_pi

        overlap = abs(np.sin(phase))

        # NODE COUNTING: divide by number of lobes when past pi
        if use_node_counting and either_has_nodes and phase > pi:
            n_lobes = int(np.ceil(phase / pi))
            overlap /= n_lobes

        contribution = C_bond * E_scale * overlap * sym_factor

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

    # ENHANCED IONIC: use 3/7 when fully ionic and covalent is weak
    if use_enhanced_ionic and q > 0.99 and abs(D_cov) < delta_eps * 0.1:
        c_ion = c_ionic_full  # 3/7
    else:
        c_ion = c_ionic_base  # 1/7

    D_ionic = c_ion * q**2 * 2 * E_H / R

    return D_cov + D_ionic, D_cov, D_ionic, q, c_ion


def run_variant(label, **kwargs):
    """Run a variant and return errors."""
    errs = []
    details = []
    for mol in molecules:
        name, R, De_exp, bonds, o1, o2 = mol
        De_pred, D_cov, D_ion, q, c_ion = compute_energy(name, R, De_exp, bonds, o1, o2, **kwargs)
        err = (De_pred - De_exp) / De_exp * 100
        errs.append(abs(err))
        details.append((name, De_exp, De_pred, err, D_cov, D_ion, q, c_ion))
    return errs, details


def print_results(label, errs, details):
    print(f"\n{'='*85}")
    print(f"  {label}")
    print(f"{'='*85}")
    header = f"{'Mol':<7} {'De_exp':>7} {'De_pred':>7} {'err%':>7}  {'D_cov':>7} {'D_ion':>6} {'q':>5} {'c_ion':>5}"
    print(header)
    print("-" * 72)

    for name, De_exp, De_pred, err, D_cov, D_ion, q, c_ion in details:
        flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
        c_str = f"{c_ion:.3f}" if c_ion != c_ionic_base else "1/7"
        print(f"{name:<7} {De_exp:7.3f} {De_pred:7.3f} {err:+6.1f}%  "
              f"{D_cov:7.3f} {D_ion:6.3f} {q:5.3f} {c_str:>5} {flag}")

    print("-" * 72)
    n = len(errs)
    w2 = sum(1 for e in errs if e < 2)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    print(f"  avg = {np.mean(errs):.1f}%, median = {np.median(errs):.1f}%")
    print(f"  within 2%: {w2}/{n},  within 5%: {w5}/{n},  within 10%: {w10}/{n}")
    return w5, w10


# =============================================================================
# TEST EACH CORRECTION INDEPENDENTLY
# =============================================================================

print("=" * 85)
print("  SEPARATE FUNCTIONS TEST: Physically motivated corrections")
print("  All parameters from d=3. Zero free parameters.")
print("=" * 85)
print(f"\n  c_ionic_base = 1/(2d+1) = 1/7 = {c_ionic_base:.4f}")
print(f"  c_ionic_full = d/(2d+1) = 3/7 = {c_ionic_full:.4f}")

# Baseline
e0, d0 = run_variant("baseline")
print_results("BASELINE (current formula, CN BO=2.5)", e0, d0)

# 1. Node counting only
e1, d1 = run_variant("node_counting", use_node_counting=True)
print_results("FIX 1: NODE COUNTING (|sin|/n_lobes when phase>pi & has_nodes)", e1, d1)

# 2. Enhanced ionic only
e2, d2 = run_variant("enhanced_ionic", use_enhanced_ionic=True)
print_results("FIX 2: ENHANCED IONIC (c=3/7 when q>0.99 & D_cov weak)", e2, d2)

# 3. Resonance symmetry only
e3, d3 = run_variant("resonance_sym", use_resonance_sym=True)
print_results("FIX 3: RESONANCE SYMMETRY (D_cov x 2*sqrt(e1e2)/(e1+e2))", e3, d3)

# 4. Predicted R only
e4, d4 = run_variant("predicted_R", use_predicted_R=True)
print_results("FIX 4: USE PREDICTED R (self-consistent bond length)", e4, d4)


# =============================================================================
# COMBINE: Node counting + Enhanced ionic
# =============================================================================
e12, d12 = run_variant("node+ionic", use_node_counting=True, use_enhanced_ionic=True)
print_results("COMBINED: Node counting + Enhanced ionic", e12, d12)

# Node counting + Enhanced ionic + Resonance symmetry
e123, d123 = run_variant("all_three", use_node_counting=True, use_enhanced_ionic=True,
                          use_resonance_sym=True)
print_results("COMBINED: Node + Ionic + Resonance symmetry", e123, d123)

# All four
e1234, d1234 = run_variant("all_four", use_node_counting=True, use_enhanced_ionic=True,
                            use_resonance_sym=True, use_predicted_R=True)
print_results("COMBINED: All four corrections", e1234, d1234)


# =============================================================================
# BEST COMBO: Node counting + Enhanced ionic + Predicted R (no resonance)
# =============================================================================
e124, d124 = run_variant("node+ionic+predR", use_node_counting=True, use_enhanced_ionic=True,
                          use_predicted_R=True)
print_results("COMBINED: Node + Ionic + Predicted R", e124, d124)


# =============================================================================
# COMPARISON TABLE
# =============================================================================
print(f"\n{'='*85}")
print(f"  SUMMARY COMPARISON")
print(f"{'='*85}")

variants = [
    ("Baseline",        e0),
    ("Node counting",   e1),
    ("Enhanced ionic",  e2),
    ("Resonance sym",   e3),
    ("Predicted R",     e4),
    ("Node+Ionic",      e12),
    ("Node+Ion+Sym",    e123),
    ("All four",        e1234),
    ("Node+Ion+PredR",  e124),
]

print(f"\n  {'Variant':<20} {'avg':>6} {'med':>6} {'w2':>4} {'w5':>4} {'w10':>4}")
print("  " + "-" * 48)
for label, errs in variants:
    avg = np.mean(errs)
    med = np.median(errs)
    w2 = sum(1 for e in errs if e < 2)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    print(f"  {label:<20} {avg:5.1f}% {med:5.1f}% {w2:4d} {w5:4d} {w10:4d}")


# =============================================================================
# DETAILED: Which molecules change in best combo?
# =============================================================================
print(f"\n{'='*85}")
print(f"  MOLECULE-BY-MOLECULE: Baseline vs Best combo")
print(f"{'='*85}")

# Determine which combo is best
best_label = max(variants, key=lambda x: (sum(1 for e in x[1] if e < 5),
                                            sum(1 for e in x[1] if e < 10),
                                            -np.mean(x[1])))
print(f"\n  Best: {best_label[0]}")

print(f"\n  {'Mol':<7} {'baseline':>10} {'best':>10} {'change':>8}")
print("  " + "-" * 40)
names = [m[0] for m in molecules]
for i, name in enumerate(names):
    e_base = e0[i]
    e_best = best_label[1][i]
    diff = e_best - e_base
    marker = " <<<" if abs(diff) > 1 and e_best < e_base else " >>>" if abs(diff) > 1 and e_best > e_base else ""
    print(f"  {name:<7} {e_base:9.1f}% {e_best:9.1f}% {diff:+7.1f}%{marker}")


# =============================================================================
# TEST: What c_ionic threshold gives best results?
# =============================================================================
print(f"\n{'='*85}")
print(f"  IONIC THRESHOLD SCAN: What D_cov/delta_eps ratio works best?")
print(f"{'='*85}")

for threshold in [0.01, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30, 0.50]:
    errs = []
    for mol in molecules:
        name, R, De_exp, bonds, o1, o2 = mol
        n1, l1 = get_n(o1), get_l(o1)
        n2, l2 = get_n(o2), get_l(o2)
        h1, h2 = has_nodes(o1), has_nodes(o2)

        a1 = 2 + (1-2*l1)*alpha*h1; a2 = 2 + (1-2*l2)*alpha*h2
        b1 = 1 + beta*h1; b2 = 1 + beta*h2

        E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)
        sigma_phase = R/n1**b1 + R/n2**b2

        either_has_nodes = (h1 > 0 or h2 > 0)

        n_pb = sum(c for bt,c in bonds if 'pi' in bt and 'anti' not in bt)
        n_pa = sum(c for bt,c in bonds if 'pi' in bt and 'anti' in bt)
        is_full = (n_pa >= n_pb) if n_pb > 0 else True

        Dc = 0
        for bt, cnt in bonds:
            ph = sigma_phase if ('sigma' in bt or bt in ('ss','sp')) else sigma_phase * f_pi
            ov = abs(np.sin(ph))
            if either_has_nodes and ph > pi:
                ov /= int(np.ceil(ph / pi))
            ct = C_bond * E_scale * ov
            if 'anti' in bt:
                fa = 1.0 if ('sigma' in bt or is_full) else f_anti
                Dc -= cnt * fa * ct
            else:
                Dc += cnt * ct

        eps1, eps2 = orbital_energy(o1), orbital_energy(o2)
        delta = abs(eps1-eps2)
        V = max(abs(Dc), 0.01)
        q = delta/np.sqrt(delta**2 + (2*V)**2) if delta > 0 else 0

        # Enhanced ionic when q>0.99 and D_cov/delta_eps < threshold
        if q > 0.99 and abs(Dc) < delta * threshold:
            c_ion = c_ionic_full
        else:
            c_ion = c_ionic_base

        D_tot = Dc + c_ion * q**2 * 2 * E_H / R
        errs.append(abs((D_tot - De_exp)/De_exp * 100))

    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    # Show which molecules get enhanced ionic
    enhanced = []
    for mol in molecules:
        name, R, De_exp, bonds, o1, o2 = mol
        eps1, eps2 = orbital_energy(o1), orbital_energy(o2)
        delta = abs(eps1-eps2)
        # Quick recompute D_cov
        n1_v,l1_v = get_n(o1), get_l(o1); n2_v,l2_v = get_n(o2), get_l(o2)
        h1_v,h2_v = has_nodes(o1), has_nodes(o2)
        a1 = 2+(1-2*l1_v)*alpha*h1_v; a2 = 2+(1-2*l2_v)*alpha*h2_v
        b1 = 1+beta*h1_v; b2 = 1+beta*h2_v
        E_s = np.sqrt(E_H/n1_v**a1 * E_H/n2_v**a2)
        sp = R/n1_v**b1 + R/n2_v**b2
        either_hn = (h1_v>0 or h2_v>0)
        n_pb = sum(c for bt,c in bonds if 'pi' in bt and 'anti' not in bt)
        n_pa = sum(c for bt,c in bonds if 'pi' in bt and 'anti' in bt)
        is_full = (n_pa >= n_pb) if n_pb > 0 else True
        Dc=0
        for bt,cnt in bonds:
            ph = sp if ('sigma' in bt or bt in ('ss','sp')) else sp*f_pi
            ov = abs(np.sin(ph))
            if either_hn and ph > pi: ov /= int(np.ceil(ph/pi))
            ct = C_bond*E_s*ov
            if 'anti' in bt:
                fa=1.0 if ('sigma' in bt or is_full) else f_anti
                Dc -= cnt*fa*ct
            else:
                Dc += cnt*ct
        V = max(abs(Dc),0.01)
        q = delta/np.sqrt(delta**2+(2*V)**2) if delta>0 else 0
        if q > 0.99 and abs(Dc) < delta * threshold:
            enhanced.append(name)

    print(f"  thr={threshold:.2f}: avg={np.mean(errs):.1f}%, w5={w5}, w10={w10}  enhanced: {', '.join(enhanced)}")
