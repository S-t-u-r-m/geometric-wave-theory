"""
Systematically test modifications to fix the 7 outliers.
Goal: get as many as possible under 5%.

Current: 15/24 within 5%, avg 14.4%, median 4.1%
Outliers: LiH(+65%), NaH(+46%), CH(-40%), LiF(-41%), NaCl(-47%), BF(+27%), CN(+31%)

Modifications to test:
  A. sinc(phase) = sin(phase)/phase instead of |sin(phase)|
     Physical: 3D Green's function is sin(kR)/(kR)
  B. Different C_bond normalization to compensate
  C. Enhanced ionic coupling for full charge transfer
  D. Asymmetric triple bond damping
"""
import numpy as np
from itertools import product

pi = np.pi
E_H = 13.6057
d = 3

# Current parameters from d=3
C_bond = pi / d
f_pi = d**2 / (d**2 + 1)     # 9/10
alpha = 1 - f_pi / d          # 7/10
beta = (1 + f_pi) / 2         # 19/20
f_anti = 2*d / (2*d - 1)      # 6/5
c_ionic = 1.0 / (2*d + 1)     # 1/7

Z_eff = {
    'H_1s': 1.0000, 'Li_2s': 1.2792, 'B_2p': 2.4214,
    'C_2p': 3.1358, 'N_2p': 3.8340, 'O_2p': 4.4532,
    'F_2p': 5.0998, 'Na_3s': 2.5074, 'Cl_3p': 6.1161,
}

def get_n(orb): return int(orb.split('_')[1][0])
def get_l(orb): return {'s': 0, 'p': 1}[orb.split('_')[1][1]]
def has_nodes(orb): return min(get_n(orb) - get_l(orb) - 1, 1)

molecules = [
    ('H2',   1.401,  4.745, [('ss', 1)],           'H_1s',  'H_1s'),
    ('Li2',  5.051,  1.056, [('ss', 1)],            'Li_2s', 'Li_2s'),
    ('B2',   3.005,  3.02,  [('pi', 2)],            'B_2p',  'B_2p'),
    ('C2',   2.348,  6.32,  [('pi', 2), ('sp_sigma', 1), ('sp_sigma_anti', 1)], 'C_2p', 'C_2p'),
    ('N2',   2.074,  9.759, [('pp_sigma', 1), ('pi', 2)], 'N_2p', 'N_2p'),
    ('O2',   2.282,  5.213, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)], 'O_2p', 'O_2p'),
    ('F2',   2.668,  1.660, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)], 'F_2p', 'F_2p'),
    ('Na2',  5.818,  0.746, [('ss', 1)],            'Na_3s', 'Na_3s'),
    ('Cl2',  3.757,  2.514, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)], 'Cl_3p', 'Cl_3p'),
    ('HF',   1.733,  5.869, [('sp', 1)],            'H_1s',  'F_2p'),
    ('CO',   2.132, 11.225, [('pp_sigma', 1), ('pi', 2)], 'C_2p', 'O_2p'),
    ('NO',   2.175,  6.497, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)], 'N_2p', 'O_2p'),
    ('OH',   1.834,  4.392, [('sp', 1)],            'O_2p',  'H_1s'),
    ('HCl',  2.409,  4.434, [('sp', 1)],            'H_1s',  'Cl_3p'),
    ('LiH',  3.015,  2.515, [('ss', 1)],            'Li_2s', 'H_1s'),
    ('LiF',  2.955,  5.939, [('sp', 1)],            'Li_2s', 'F_2p'),
    ('BH',   2.329,  3.42,  [('sp', 1)],            'B_2p',  'H_1s'),
    ('CH',   2.116,  3.47,  [('sp', 1)],            'C_2p',  'H_1s'),
    ('NH',   1.958,  3.57,  [('sp', 1)],            'N_2p',  'H_1s'),
    ('BF',   2.386,  7.81,  [('pp_sigma', 1), ('pi', 2)], 'B_2p', 'F_2p'),
    ('CN',   2.214,  7.72,  [('pp_sigma', 1), ('pi', 2)], 'C_2p', 'N_2p'),
    ('NaH',  3.566,  1.97,  [('ss', 1)],            'Na_3s', 'H_1s'),
    ('NaCl', 4.461,  4.23,  [('sp', 1)],            'Na_3s', 'Cl_3p'),
    ('H2O',  1.809,  5.117, [('sp', 1)],            'O_2p',  'H_1s'),
]


def compute_all(overlap_func, C=C_bond, c_ion=c_ionic, f_a=f_anti,
                verbose=False, label=""):
    """Compute all bond energies with a given overlap function."""
    errs = []
    results = []
    for mol in molecules:
        name, R, De_exp, bonds, o1, o2 = mol
        n1, l1, n2, l2 = get_n(o1), get_l(o1), get_n(o2), get_l(o2)
        h1, h2 = has_nodes(o1), has_nodes(o2)

        a1 = 2 + (1-2*l1)*alpha*h1
        a2 = 2 + (1-2*l2)*alpha*h2
        b1 = 1 + beta*h1
        b2 = 1 + beta*h2

        E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)
        sigma_phase = R/n1**b1 + R/n2**b2

        npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
        npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
        ifa = (npa >= npb) if npb > 0 else True

        D_cov = 0
        for bt, cnt in bonds:
            ph = sigma_phase if ('sigma' in bt or bt in ('ss','sp')) else sigma_phase*f_pi
            cont = C * E_scale * overlap_func(ph)
            if 'anti' in bt:
                fa_use = 1.0 if ('sigma' in bt or ifa) else f_a
                D_cov -= cnt * fa_use * cont
            else:
                D_cov += cnt * cont

        eps1 = E_H*(Z_eff[o1]/n1)**2
        eps2 = E_H*(Z_eff[o2]/n2)**2
        dE = abs(eps1 - eps2)
        V = max(abs(D_cov), 0.01)
        q = dE/np.sqrt(dE**2 + (2*V)**2) if dE > 0 else 0
        D_ion = c_ion * q**2 * 2*E_H/R
        D_pred = D_cov + D_ion
        err = (D_pred - De_exp)/De_exp * 100
        errs.append(abs(err))
        results.append((name, De_exp, D_pred, D_cov, D_ion, err, sigma_phase))

    if verbose:
        print(f"\n{'='*80}")
        print(f"  {label}")
        print(f"{'='*80}")
        print(f"  {'Mol':<6} {'De_exp':>7} {'D_pred':>7} {'D_cov':>7} {'D_ion':>6} {'err%':>7} {'ph/pi':>6}")
        print("  " + "-"*55)
        for name, De, Dp, Dc, Di, err, ph in results:
            flag = '***' if abs(err)<2 else ' **' if abs(err)<5 else '  *' if abs(err)<10 else ''
            print(f"  {name:<6} {De:7.3f} {Dp:7.3f} {Dc:7.3f} {Di:6.3f} {err:+6.1f}% {ph/pi:6.3f} {flag}")

        w2 = sum(1 for e in errs if e<2)
        w5 = sum(1 for e in errs if e<5)
        w10 = sum(1 for e in errs if e<10)
        print(f"\n  avg={np.mean(errs):.1f}%, med={np.median(errs):.1f}%")
        print(f"  w2={w2}/24, w5={w5}/24, w10={w10}/24")

    return errs, results


# =============================================================================
# BASELINE
# =============================================================================
print("BASELINE (current formula):")
errs0, _ = compute_all(lambda ph: abs(np.sin(ph)), verbose=True,
                        label="CURRENT: |sin(phase)|, C=pi/3, c_ionic=1/7")


# =============================================================================
# TEST A: sinc-like overlap
# =============================================================================
# In 3D, Green's function ~ sin(kR)/(kR)
# sinc naturally decays for large phase, fixing wrapping

# A1: |sin(ph)/ph|  (raw sinc, unnormalized)
compute_all(lambda ph: abs(np.sin(ph)/ph), verbose=True,
            label="A1: |sin(ph)/ph| (raw sinc)")

# The sinc gives smaller values, so we might need larger C
# But C = pi/d is from first principles. Let's see what C gives best results
print("\n" + "="*80)
print("  SCANNING C for sinc overlap |sin(ph)/ph|")
print("="*80)
best_sinc = (999, 0)
for C_test in np.arange(0.5, 5.0, 0.01):
    errs, _ = compute_all(lambda ph: abs(np.sin(ph)/ph), C=C_test)
    avg = np.mean(errs)
    if avg < best_sinc[0]:
        best_sinc = (avg, C_test)
print(f"  Best C = {best_sinc[1]:.2f} (avg err = {best_sinc[0]:.1f}%)")
print(f"  Compare: pi/3 = {pi/3:.4f}, pi = {pi:.4f}, pi/2 = {pi/2:.4f}")
print(f"  Compare: pi^2/d = {pi**2/d:.4f}")

# Test with best C
compute_all(lambda ph: abs(np.sin(ph)/ph), C=best_sinc[1], verbose=True,
            label=f"A1 optimized: |sin(ph)/ph|, C={best_sinc[1]:.2f}")


# A2: |sin(ph)/(ph/pi)| = pi*|sin(ph)/ph|  (sinc normalized so max=pi at ph=0)
compute_all(lambda ph: pi*abs(np.sin(ph)/ph), verbose=True,
            label="A2: pi*|sin(ph)/ph| (sinc normalized to pi at origin)")


# A3: What about sin(ph)/sinh(ph)? Exponential decay
# This is the ratio of circular to hyperbolic sine - decays as exp(-ph) for large ph
def sin_over_sinh(ph):
    if abs(ph) < 1e-10: return 1.0
    return abs(np.sin(ph)/np.sinh(ph))

compute_all(sin_over_sinh, verbose=True,
            label="A3: |sin(ph)/sinh(ph)| (exponential decay)")

# Scan C for sin/sinh
best_ss = (999, 0)
for C_test in np.arange(0.5, 10.0, 0.01):
    errs, _ = compute_all(sin_over_sinh, C=C_test)
    avg = np.mean(errs)
    if avg < best_ss[0]:
        best_ss = (avg, C_test)
print(f"\n  Best C for sin/sinh = {best_ss[1]:.2f} (avg err = {best_ss[0]:.1f}%)")

compute_all(sin_over_sinh, C=best_ss[1], verbose=True,
            label=f"A3 optimized: |sin/sinh|, C={best_ss[1]:.2f}")


# =============================================================================
# TEST B: Different overlap functions with d=3 motivated C
# =============================================================================
print("\n" + "="*80)
print("  SYSTEMATIC OVERLAP FUNCTION SCAN")
print("="*80)

# Try various overlap functions with C = pi/d and with optimized C
overlap_funcs = {
    '|sin(ph)|':       lambda ph: abs(np.sin(ph)),
    '|sin(ph)/ph|':    lambda ph: abs(np.sin(ph)/ph) if abs(ph)>1e-10 else 1.0,
    '|sin(ph)|/ph':    lambda ph: abs(np.sin(ph))/max(ph, 0.01),
    'sin^2(ph/2)/ph':  lambda ph: np.sin(ph/2)**2/max(ph, 0.01),
    '|sin(ph)|*exp(-ph/pi)': lambda ph: abs(np.sin(ph))*np.exp(-ph/pi),
    '|sin(ph)|/(1+ph/pi)':   lambda ph: abs(np.sin(ph))/(1+ph/pi),
    '|sin(ph)|/(1+(ph/pi)^2)': lambda ph: abs(np.sin(ph))/(1+(ph/pi)**2),
    'sech(ph-pi/2)*|sin|':     lambda ph: abs(np.sin(ph))/np.cosh(ph-pi/2),
}

print(f"\n  {'Overlap func':>30} {'C_best':>6} {'avg%':>6} {'med%':>6} {'w5':>4} {'w10':>4}")
print("  " + "-"*60)

best_overall = (999, '', 0)
for name, func in overlap_funcs.items():
    # Find best C
    best = (999, 0)
    for C_test in np.arange(0.3, 8.0, 0.005):
        errs, _ = compute_all(func, C=C_test)
        avg = np.mean(errs)
        if avg < best[0]:
            best = (avg, C_test)

    errs, _ = compute_all(func, C=best[1])
    w5 = sum(1 for e in errs if e<5)
    w10 = sum(1 for e in errs if e<10)
    med = np.median(errs)
    print(f"  {name:>30} {best[1]:6.3f} {best[0]:6.1f} {med:6.1f} {w5:4d} {w10:4d}")

    if best[0] < best_overall[0]:
        best_overall = (best[0], name, best[1])

print(f"\n  BEST: {best_overall[1]} with C={best_overall[2]:.3f} (avg={best_overall[0]:.1f}%)")

# Show detailed results for the best
best_name = best_overall[1]
best_func = overlap_funcs[best_name]
best_C = best_overall[2]
compute_all(best_func, C=best_C, verbose=True,
            label=f"BEST OVERLAP: {best_name}, C={best_C:.3f}")


# =============================================================================
# TEST C: Different ionic coupling
# =============================================================================
print("\n" + "="*80)
print("  IONIC COUPLING SCAN (with best overlap)")
print("="*80)

# Scan c_ionic with the best overlap function
best_ion = (999, 0, 0)
for C_test in np.arange(0.3, 8.0, 0.02):
    for ci_test in np.arange(0.05, 0.6, 0.005):
        errs, _ = compute_all(best_func, C=C_test, c_ion=ci_test)
        avg = np.mean(errs)
        if avg < best_ion[0]:
            best_ion = (avg, C_test, ci_test)

print(f"  Best: C={best_ion[1]:.3f}, c_ionic={best_ion[2]:.4f} (avg={best_ion[0]:.1f}%)")
print(f"  Compare: 1/7={1/7:.4f}, 1/d={1/d:.4f}, pi/d^2={pi/d**2:.4f}")

compute_all(best_func, C=best_ion[1], c_ion=best_ion[2], verbose=True,
            label=f"BEST OVERLAP + IONIC: C={best_ion[1]:.3f}, c_ion={best_ion[2]:.4f}")


# =============================================================================
# TEST D: Keep C=pi/3 and c_ionic=1/7 but change overlap
# =============================================================================
# We want ZERO free parameters. So let's see which overlap with the
# d=3 constants works best.
print("\n" + "="*80)
print("  ZERO FREE PARAMETERS: Which overlap with C=pi/3, c=1/7?")
print("="*80)

# Key d=3 constant candidates for C:
C_candidates = {
    'pi/d = pi/3':   pi/3,
    'pi^2/(2d+1)':   pi**2/7,
    '2pi/d^2':       2*pi/9,
    'pi/(d-1)':      pi/2,
    '1':             1.0,
    'pi^2/d^2':      pi**2/9,
    'd':             3.0,
    'pi':            pi,
}

overlap_zero_param = {
    '|sin(ph)|':                  lambda ph: abs(np.sin(ph)),
    '|sin(ph)/ph|':               lambda ph: abs(np.sin(ph)/ph) if abs(ph)>1e-10 else 1.0,
    '|sin(ph)|/(1+ph/pi)':        lambda ph: abs(np.sin(ph))/(1+ph/pi),
    '|sin(ph)|/(1+(ph/pi)^2)':    lambda ph: abs(np.sin(ph))/(1+(ph/pi)**2),
    '|sin(ph)|*exp(-ph/pi)':      lambda ph: abs(np.sin(ph))*np.exp(-ph/pi),
    '|sin/sinh|':                 lambda ph: abs(np.sin(ph)/np.sinh(ph)) if abs(ph)>1e-10 else 1.0,
}

print(f"\n  {'Overlap':>28} {'C formula':>18} {'C val':>6} {'avg%':>6} {'med%':>6} {'w5':>4} {'w10':>4}")
print("  " + "-"*80)

best_zp = (999, '', '', 0)
for ov_name, ov_func in overlap_zero_param.items():
    for c_name, c_val in C_candidates.items():
        errs, _ = compute_all(ov_func, C=c_val)
        avg = np.mean(errs)
        w5 = sum(1 for e in errs if e<5)
        w10 = sum(1 for e in errs if e<10)
        med = np.median(errs)
        if avg < 20:  # only show reasonable ones
            print(f"  {ov_name:>28} {c_name:>18} {c_val:6.3f} {avg:6.1f} {med:6.1f} {w5:4d} {w10:4d}")
        if avg < best_zp[0]:
            best_zp = (avg, ov_name, c_name, c_val)

print(f"\n  BEST zero-param: {best_zp[1]} with C={best_zp[2]} = {best_zp[3]:.4f} (avg={best_zp[0]:.1f}%)")

# Show the best zero-param result
compute_all(overlap_zero_param[best_zp[1]], C=best_zp[3], verbose=True,
            label=f"BEST ZERO-PARAM: {best_zp[1]}, C={best_zp[2]}")


# =============================================================================
# TEST E: The SPECIFIC problem — can we fix just the phase formula?
# =============================================================================
print("\n" + "="*80)
print("  PHASE FORMULA VARIATIONS (keep |sin|, C=pi/3)")
print("="*80)

# What if the phase uses the BONDING radius instead of R?
# phase = r_bond_1/n1^b + r_bond_2/n2^b where r_bond = n^2/(Z*(n-l))
# This makes phase = n/(Z*(n-l)) for each atom

print("\nPhase = sum of n_i/(Z_eff_i * (n_i - l_i)):")
errs_br = []
results_br = []
for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    n1, l1, n2, l2 = get_n(o1), get_l(o1), get_n(o2), get_l(o2)
    h1, h2 = has_nodes(o1), has_nodes(o2)

    a1 = 2+(1-2*l1)*alpha*h1; a2 = 2+(1-2*l2)*alpha*h2
    b1 = 1+beta*h1; b2 = 1+beta*h2
    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)

    # Bond radius-based phase
    r1 = n1**2 / (Z_eff[o1] * (n1 - l1))
    r2 = n2**2 / (Z_eff[o2] * (n2 - l2))
    sigma_phase = r1/n1**b1 + r2/n2**b2

    npb = sum(c for bt,c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt,c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss','sp')) else sigma_phase*f_pi
        cont = C_bond * E_scale * abs(np.sin(ph))
        if 'anti' in bt:
            fa_use = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt*fa_use*cont
        else:
            D_cov += cnt*cont

    eps1 = E_H*(Z_eff[o1]/n1)**2; eps2 = E_H*(Z_eff[o2]/n2)**2
    dE = abs(eps1-eps2); V = max(abs(D_cov),0.01)
    q = dE/np.sqrt(dE**2+(2*V)**2) if dE>0 else 0
    D_ion = c_ionic*q**2*2*E_H/R
    D_pred = D_cov + D_ion
    err = (D_pred-De_exp)/De_exp*100
    errs_br.append(abs(err))

    flag = '***' if abs(err)<2 else ' **' if abs(err)<5 else '  *' if abs(err)<10 else ''
    print(f"  {name:<6} {De_exp:7.3f} {D_pred:7.3f} {err:+6.1f}% r1={r1:.3f} r2={r2:.3f} ph={sigma_phase:.3f} {flag}")

print(f"\n  avg={np.mean(errs_br):.1f}%, med={np.median(errs_br):.1f}%")


# What about a mixed phase: part R, part bonding radius?
# phase = (alpha_mix * R + (1-alpha_mix) * r_bond) * k
# With alpha_mix = some d=3 constant

print("\n\nPhase at CONTACT distance (r1 + r2) instead of R:")
errs_ct = []
for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    n1, l1, n2, l2 = get_n(o1), get_l(o1), get_n(o2), get_l(o2)
    h1, h2 = has_nodes(o1), has_nodes(o2)

    a1 = 2+(1-2*l1)*alpha*h1; a2 = 2+(1-2*l2)*alpha*h2
    b1 = 1+beta*h1; b2 = 1+beta*h2
    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)

    r1 = n1**2 / (Z_eff[o1] * (n1 - l1))
    r2 = n2**2 / (Z_eff[o2] * (n2 - l2))
    R_contact = r1 + r2

    sigma_phase = R_contact/n1**b1 + R_contact/n2**b2

    npb = sum(c for bt,c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt,c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss','sp')) else sigma_phase*f_pi
        cont = C_bond * E_scale * abs(np.sin(ph))
        if 'anti' in bt:
            fa_use = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt*fa_use*cont
        else:
            D_cov += cnt*cont

    eps1 = E_H*(Z_eff[o1]/n1)**2; eps2 = E_H*(Z_eff[o2]/n2)**2
    dE = abs(eps1-eps2); V = max(abs(D_cov),0.01)
    q = dE/np.sqrt(dE**2+(2*V)**2) if dE>0 else 0
    D_ion = c_ionic*q**2*2*E_H/R  # ionic still at actual R
    D_pred = D_cov + D_ion
    err = (D_pred-De_exp)/De_exp*100
    errs_ct.append(abs(err))

    flag = '***' if abs(err)<2 else ' **' if abs(err)<5 else '  *' if abs(err)<10 else ''
    print(f"  {name:<6} {De_exp:7.3f} {D_pred:7.3f} {err:+6.1f}% R_ct={R_contact:.3f} vs R={R:.3f} {flag}")

print(f"\n  avg={np.mean(errs_ct):.1f}%, med={np.median(errs_ct):.1f}%")
