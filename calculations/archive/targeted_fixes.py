"""
TARGETED FIXES for the 7 outliers.
The overlap scan proved |sin| with C=pi/3 is optimal globally.
The issues are structural — wrong bond types, missing physics.

Key observations:
1. CN has 9 valence electrons -> BO = 2.5, NOT 3
   We assigned [sigma, 2*pi] = BO 3. Wrong!
2. BF is isoelectronic to CO/N2 (BO=3) but overshoots.
   Maybe the large Z_eff asymmetry reduces coupling?
3. LiH/NaH: phase > pi, |sin| wraps back up
4. CH: phase = pi, sin = 0 (coincidence)
5. LiF/NaCl: ionic monopole too weak

Let's fix what we can.
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
def has_nodes(orb): return min(get_n(orb) - get_l(orb) - 1, 1)

def compute_bond(name, R, De_exp, bonds, o1, o2):
    n1, l1, n2, l2 = get_n(o1), get_l(o1), get_n(o2), get_l(o2)
    h1, h2 = has_nodes(o1), has_nodes(o2)
    a1 = 2 + (1-2*l1)*alpha*h1
    a2 = 2 + (1-2*l2)*alpha*h2
    b1 = 1 + beta*h1; b2 = 1 + beta*h2
    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)
    sigma_phase = R/n1**b1 + R/n2**b2

    npb = sum(c for bt,c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt,c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss','sp')) else sigma_phase*f_pi
        cont = C_bond * E_scale * abs(np.sin(ph))
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt*fa*cont
        else:
            D_cov += cnt*cont

    eps1 = E_H*(Z_eff[o1]/n1)**2
    eps2 = E_H*(Z_eff[o2]/n2)**2
    dE = abs(eps1 - eps2)
    V = max(abs(D_cov), 0.01)
    q = dE/np.sqrt(dE**2 + (2*V)**2) if dE > 0 else 0
    D_ion = c_ionic * q**2 * 2*E_H/R
    return D_cov + D_ion, D_cov, D_ion, q, sigma_phase, E_scale

# =============================================================================
# FIX 1: CN BOND ORDER
# =============================================================================
print("="*80)
print("  FIX 1: CN BOND ORDER")
print("="*80)
print()
print("  CN has 13 total electrons (C:6 + N:7):")
print("  Core: 1s(C)^2, 1s(N)^2 = 4")
print("  Valence MOs (9 electrons):")
print("    2sigma (bonding):     2 electrons")
print("    2sigma* (antibond):   2 electrons")
print("    1pi (bonding):        4 electrons")
print("    3sigma (bonding):     1 electron")
print("  BO = (2 + 4 + 1 - 2) / 2 = 2.5")
print()

# Current: [('pp_sigma', 1), ('pi', 2)] = BO 3
# Correct: [('pp_sigma', 1), ('pi', 2), ('pi_anti', 0.5)] = BO 2.5?
# Actually, the antibonding electron is in 2sigma*, not pi*
# Better: [('pp_sigma', 0.5), ('pi', 2)] = BO 2.5 (half sigma)
# Or: [('pp_sigma', 1), ('pi', 2), ('sp_sigma_anti', 1)] = BO 2
# Hmm, let me think about what the MO picture maps to in GWT

# In GWT, the bond_list represents which standing wave modes are occupied
# CN has one unpaired electron in a sigma bonding orbital
# Option A: half-filled sigma -> count = 0.5
# Option B: one sigma, two pi, one sigma anti (BO=2, but that's too low)
# Option C: keep as-is but note it should be 2.5

# Actually, looking at this more carefully:
# The sigma bond in CN has 1 electron, not 2. So it's a half-bond.
# In GWT, this means the standing wave pairing is only half-complete.

cn_options = [
    ("Current (BO=3)",           [('pp_sigma', 1), ('pi', 2)]),
    ("Half sigma (BO=2.5)",      [('pp_sigma', 0.5), ('pi', 2)]),
    ("Full sigma + half pi anti",[('pp_sigma', 1), ('pi', 2), ('pi_anti', 0.5)]),
    ("Double bond (BO=2)",       [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)]),
]

print(f"  {'Description':<35} {'D_pred':>7} {'De_exp':>7} {'err%':>7}")
print("  " + "-"*55)
for desc, bonds in cn_options:
    D, Dc, Di, q, ph, Es = compute_bond('CN', 2.214, 7.72, bonds, 'C_2p', 'N_2p')
    err = (D - 7.72)/7.72*100
    print(f"  {desc:<35} {D:7.3f} {7.72:7.3f} {err:+6.1f}%")


# =============================================================================
# FIX 1B: Also check NO bond order
# =============================================================================
print()
print("  Also check NO (currently BO=2):")
print("  NO has 15 electrons: N(7) + O(8)")
print("  Valence (11 electrons):")
print("    2sigma: 2, 2sigma*: 2, 1pi: 4, 3sigma: 2, 1pi*: 1")
print("  BO = (2+4+2-2-1)/2 = 2.5")
print()

no_options = [
    ("Current (BO=2)",           [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)]),
    ("BO=2.5 (half pi anti)",    [('pp_sigma', 1), ('pi', 2), ('pi_anti', 0.5)]),
]

print(f"  {'Description':<35} {'D_pred':>7} {'De_exp':>7} {'err%':>7}")
print("  " + "-"*55)
for desc, bonds in no_options:
    D, Dc, Di, q, ph, Es = compute_bond('NO', 2.175, 6.497, bonds, 'N_2p', 'O_2p')
    err = (D - 6.497)/6.497*100
    print(f"  {desc:<35} {D:7.3f} {6.497:7.3f} {err:+6.1f}%")


# =============================================================================
# FIX 2: BF TRIPLE BOND
# =============================================================================
print()
print("="*80)
print("  FIX 2: BF TRIPLE BOND ANALYSIS")
print("="*80)
print()
print("  BF: isoelectronic to N2, CO (10 valence e). BO=3 is correct.")
print("  But E_scale and phase differ from N2/CO because B and F have")
print("  very different Z_eff (2.42 vs 5.10).")
print()

# The issue: same E_scale formula gives 3.401 for N2, CO, BF, CN
# because they all use 2p orbitals with no radial nodes (h=0, a=2)
# E_scale = E_H/n^2 = E_H/4 = 3.401 regardless of Z_eff

# But physically, BF should be weaker than N2 because B's 2p is much
# less tightly bound than F's 2p. The AVERAGE makes sense for symmetric
# molecules but not for highly asymmetric ones.

# What if E_scale should be the HARMONIC mean instead of GEOMETRIC?
# Harmonic mean = 2*E1*E2/(E1+E2) penalizes asymmetry more
# Geometric mean = sqrt(E1*E2)

# With real orbital energies:
# N2: E_N = 13.606*(3.834/2)^2 = 49.999, E_scale_geom = 50.0, harm = 50.0
# CO: E_C = 33.45, E_O = 67.45, geom = 47.5, harm = 44.7
# BF: E_B = 19.94, E_F = 88.46, geom = 42.0, harm = 32.5
# CN: E_C = 33.45, E_N = 50.00, geom = 40.9, harm = 40.1

# But the formula uses E_H/n^a which is the SAME for all 2p (3.401)!
# The difference between molecules comes only from phase (different R)
# and the ionic term.

# CO works because the ionic term contributes 1.375 eV, boosting the total.
# BF's ionic term is 1.537 eV but covalent is too high at 8.415.
# N2 works because it's symmetric (no ionic term needed).

# The real question: why is BF's D_e only 7.81 when N2's is 9.76?
# They have the same E_scale and similar phase. The difference must be
# that the asymmetry reduces covalent coupling.

# Maybe the formula should multiply E_scale by a symmetry factor?
# symmetry_factor = 2*sqrt(E1*E2)/(E1+E2) = geometric/arithmetic mean ratio

# For N2: factor = 1 (symmetric)
# For CO: E_C/E_O = 33.45/67.45, factor = 2*sqrt(0.496)/(1+0.496) = 0.942
# For BF: E_B/E_F = 19.94/88.46, factor = 2*sqrt(0.225)/(1+0.225) = 0.774
# For CN: E_C/E_N = 33.45/50.00, factor = 2*sqrt(0.669)/(1+0.669) = 0.980

# This uses real orbital energies which we want to avoid...
# But what about using Z_eff to create the asymmetry penalty?
# ratio = min(Z1,Z2)/max(Z1,Z2)
# factor = 2*sqrt(ratio)/(1+ratio) = (same formula)

triples = [
    ('N2',  2.074,  9.759, [('pp_sigma', 1), ('pi', 2)], 'N_2p',  'N_2p'),
    ('CO',  2.132, 11.225, [('pp_sigma', 1), ('pi', 2)], 'C_2p',  'O_2p'),
    ('BF',  2.386,  7.81,  [('pp_sigma', 1), ('pi', 2)], 'B_2p',  'F_2p'),
    ('CN',  2.214,  7.72,  [('pp_sigma', 1), ('pi', 2)], 'C_2p',  'N_2p'),
]

print(f"  {'Mol':<5} {'Z1':>5} {'Z2':>5} {'ratio':>6} {'f_sym':>6} {'D_curr':>7} {'D_sym':>7} {'De_exp':>7}")
print("  " + "-"*55)
for name, R, De, bonds, o1, o2 in triples:
    z1, z2 = Z_eff[o1], Z_eff[o2]
    ratio = min(z1,z2)/max(z1,z2)
    f_sym = 2*np.sqrt(ratio)/(1+ratio)

    D_curr = compute_bond(name, R, De, bonds, o1, o2)[0]

    # Apply symmetry factor to covalent part
    n1, n2 = get_n(o1), get_n(o2)
    h1, h2 = has_nodes(o1), has_nodes(o2)
    b1, b2 = 1+beta*h1, 1+beta*h2
    ph = R/n1**b1 + R/n2**b2
    E_sc = np.sqrt(E_H/n1**2 * E_H/n2**2)

    # Recalculate with symmetry factor
    D_cov_sym = 0
    for bt, cnt in bonds:
        ph_use = ph if ('sigma' in bt or bt in ('ss','sp')) else ph*f_pi
        cont = C_bond * E_sc * f_sym * abs(np.sin(ph_use))
        if 'anti' not in bt:
            D_cov_sym += cnt*cont

    eps1, eps2 = E_H*(z1/n1)**2, E_H*(z2/n2)**2
    dE = abs(eps1-eps2); V = max(abs(D_cov_sym), 0.01)
    q = dE/np.sqrt(dE**2+(2*V)**2) if dE>0 else 0
    D_ion = c_ionic*q**2*2*E_H/R
    D_sym = D_cov_sym + D_ion

    print(f"  {name:<5} {z1:5.2f} {z2:5.02f} {ratio:6.3f} {f_sym:6.3f} "
          f"{D_curr:7.3f} {D_sym:7.3f} {De:7.3f}")


# =============================================================================
# FIX 3: PHASE WRAPPING — use signed sin for single bonds
# =============================================================================
print()
print("="*80)
print("  FIX 3: PHASE WRAPPING")
print("="*80)
print()

# When phase > pi, |sin| wraps back up. But physically the overlap
# should be WEAKER at larger distances.
#
# Idea: for SINGLE bonds (sigma or sp), if phase > pi,
# the overlapping region is past the first antinode.
# The amplitude should be reduced.
#
# Simple approach: clip phase to [0, pi]
# sin(min(phase, pi)) = 0 when phase >= pi
# This is too harsh — kills all covalent at large R
#
# Better: use the first-lobe approximation
# overlap = max(sin(phase), 0) for single bonds
# This means: bonding only in the first half-period

print("Test: max(sin(phase), 0) for sigma/sp bonds")
print("       (antibonding region gives zero, not reflected)")
print()

molecules_all = [
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


def compute_signed(mol, use_signed_sigma=False, cn_fix=False, f_sym_func=None):
    """Compute with optional fixes."""
    name, R, De_exp, bonds_orig, o1, o2 = mol
    n1, l1, n2, l2 = get_n(o1), get_l(o1), get_n(o2), get_l(o2)
    h1, h2 = has_nodes(o1), has_nodes(o2)
    a1 = 2+(1-2*l1)*alpha*h1; a2 = 2+(1-2*l2)*alpha*h2
    b1 = 1+beta*h1; b2 = 1+beta*h2
    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)
    sigma_phase = R/n1**b1 + R/n2**b2

    # Apply symmetry factor if provided
    if f_sym_func:
        z1, z2 = Z_eff[o1], Z_eff[o2]
        E_scale *= f_sym_func(z1, z2)

    # CN fix: adjust bonds
    bonds = list(bonds_orig)
    if cn_fix and name == 'CN':
        bonds = [('pp_sigma', 0.5), ('pi', 2)]

    npb = sum(c for bt,c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt,c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss','sp')) else sigma_phase*f_pi
        if use_signed_sigma and ('sigma' in bt or bt in ('ss','sp')) and 'anti' not in bt:
            overlap = max(np.sin(ph), 0)  # only first lobe
        else:
            overlap = abs(np.sin(ph))
        cont = C_bond * E_scale * overlap
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt*fa*cont
        else:
            D_cov += cnt*cont

    eps1 = E_H*(Z_eff[o1]/n1)**2; eps2 = E_H*(Z_eff[o2]/n2)**2
    dE = abs(eps1-eps2); V = max(abs(D_cov), 0.01)
    q = dE/np.sqrt(dE**2+(2*V)**2) if dE>0 else 0
    D_ion = c_ionic*q**2*2*E_H/R
    return D_cov+D_ion, D_cov, D_ion, q, sigma_phase


# Test signed sigma
print(f"  {'Mol':<6} {'De_exp':>7} {'D_curr':>7} {'D_sign':>7} {'err_c%':>7} {'err_s%':>7} {'ph/pi':>6}")
print("  " + "-"*60)

errs_curr = []; errs_sign = []
for mol in molecules_all:
    name = mol[0]; De = mol[2]
    D_c = compute_signed(mol, use_signed_sigma=False)[0]
    D_s = compute_signed(mol, use_signed_sigma=True)[0]
    ph = compute_signed(mol)[4]
    ec = (D_c-De)/De*100; es = (D_s-De)/De*100
    errs_curr.append(abs(ec)); errs_sign.append(abs(es))
    flag = '***' if abs(es)<2 else ' **' if abs(es)<5 else '  *' if abs(es)<10 else ''
    ch = '<--' if abs(es) < abs(ec) - 1 else ''
    print(f"  {name:<6} {De:7.3f} {D_c:7.3f} {D_s:7.3f} {ec:+6.1f}% {es:+6.1f}% {ph/pi:6.3f} {flag} {ch}")

print(f"\n  Current:  avg={np.mean(errs_curr):.1f}%, med={np.median(errs_curr):.1f}%")
print(f"  Signed:   avg={np.mean(errs_sign):.1f}%, med={np.median(errs_sign):.1f}%")
print(f"  w5: {sum(1 for e in errs_curr if e<5)} -> {sum(1 for e in errs_sign if e<5)}")
print(f"  w10: {sum(1 for e in errs_curr if e<10)} -> {sum(1 for e in errs_sign if e<10)}")


# =============================================================================
# FIX 4: ALL FIXES COMBINED
# =============================================================================
print()
print("="*80)
print("  COMBINED FIXES: signed sigma + CN BO=2.5")
print("="*80)
print()

# Symmetry factor function
def f_sym_zeff(z1, z2):
    ratio = min(z1,z2)/max(z1,z2)
    return 2*np.sqrt(ratio)/(1+ratio)

combos = [
    ("signed sigma only",               True,  False, None),
    ("CN BO=2.5 only",                  False, True,  None),
    ("signed + CN",                     True,  True,  None),
    ("signed + CN + Z_sym",             True,  True,  f_sym_zeff),
    ("Z_sym only",                      False, False, f_sym_zeff),
    ("signed + Z_sym",                  True,  False, f_sym_zeff),
]

for desc, sign, cn, fsym in combos:
    errs = []
    for mol in molecules_all:
        D = compute_signed(mol, use_signed_sigma=sign, cn_fix=cn, f_sym_func=fsym)[0]
        err = (D - mol[2])/mol[2]*100
        errs.append(abs(err))
    w5 = sum(1 for e in errs if e<5)
    w10 = sum(1 for e in errs if e<10)
    print(f"  {desc:<35} avg={np.mean(errs):5.1f}%, med={np.median(errs):5.1f}%, w5={w5}, w10={w10}")

# Show the best combo in detail
print()
print("Detailed: signed sigma + CN BO=2.5")
print(f"  {'Mol':<6} {'De_exp':>7} {'D_pred':>7} {'D_cov':>7} {'D_ion':>6} {'err%':>7} {'ph/pi':>6}")
print("  " + "-"*60)
errs_best = []
for mol in molecules_all:
    name = mol[0]; De = mol[2]
    D, Dc, Di, q, ph = compute_signed(mol, use_signed_sigma=True, cn_fix=True)
    err = (D-De)/De*100
    errs_best.append(abs(err))
    flag = '***' if abs(err)<2 else ' **' if abs(err)<5 else '  *' if abs(err)<10 else ''
    print(f"  {name:<6} {De:7.3f} {D:7.3f} {Dc:7.3f} {Di:6.3f} {err:+6.1f}% {ph/pi:6.3f} {flag}")

w5 = sum(1 for e in errs_best if e<5)
w10 = sum(1 for e in errs_best if e<10)
print(f"\n  avg={np.mean(errs_best):.1f}%, med={np.median(errs_best):.1f}%")
print(f"  w5={w5}/24, w10={w10}/24")


# =============================================================================
# REMAINING OUTLIERS ANALYSIS
# =============================================================================
print()
print("="*80)
print("  REMAINING OUTLIERS AFTER FIXES")
print("="*80)
print()
for i, mol in enumerate(molecules_all):
    if errs_best[i] > 10:
        name = mol[0]; De = mol[2]; R = mol[1]; o1 = mol[4]; o2 = mol[5]
        D, Dc, Di, q, ph = compute_signed(mol, use_signed_sigma=True, cn_fix=True)
        n1, l1, n2, l2 = get_n(o1), get_l(o1), get_n(o2), get_l(o2)

        print(f"  {name}: err={errs_best[i]:+.1f}%, D_cov={Dc:.3f}, D_ion={Di:.3f}, q={q:.3f}")
        print(f"    phase={ph:.3f} ({ph/pi:.3f}*pi), R={R:.3f}")
        print(f"    E_scale=sqrt(E_H/{n1}^a * E_H/{n2}^a)")

        if ph > pi:
            print(f"    PHASE > pi: signed sin -> max(sin,0) = 0, bond is purely ionic")
            print(f"    Need D_ion = {De:.3f} but have {Di:.3f}")
            print(f"    c_ionic needed = {(De-0)/(q**2*2*E_H/R):.4f} vs current {c_ionic:.4f}")
        elif q > 0.95:
            print(f"    FULLY IONIC: c_ionic needed = {(De-Dc)/(q**2*2*E_H/R):.4f}")
        else:
            print(f"    Covalent overshoot by {Dc-De:.2f} eV")
        print()
