"""
RECONCILIATION: Analyze all 7 outliers with correct Cl data.
Are these FORMULA problems or PHYSICS we're missing?

The 7 outliers:
  Phase wrapped:    LiH (+65%), NaH (+46%), CH (-40%)
  Ionic too weak:   LiF (-41%), NaCl (-47%)
  Cov overshoot:    BF (+27%), CN (+31%)

Question: did bad Cl data (4.8864 vs 6.1161) lead us to wrong assumptions?
"""
import numpy as np
pi = np.pi
E_H = 13.6057

d = 3
C_bond = pi / d
f_pi = d**2 / (d**2 + 1)    # 9/10
alpha = 1 - f_pi / d         # 7/10
beta = (1 + f_pi) / 2        # 19/20
f_anti = 2*d / (2*d - 1)     # 6/5
c_ionic = 1.0 / (2*d + 1)    # 1/7

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

print("="*80)
print("  CATEGORY 1: PHASE-WRAPPED MOLECULES")
print("  LiH, NaH, CH -- and BH which is borderline")
print("="*80)
print()

# All involve H (n=1, l=0, no nodes -> b=1, k=1)
# Bonded to atoms with radial nodes (n-l-1 >= 1)
# The node correction b = 1 + beta*h makes k = 1/n^(1+beta) very small
# But k_H = 1, so phase = R * (1 + 1/n^b) ~ R * 1 for large n

# Key question: should H really have k=1 (b=1)?
# H has no nodes, so h=0, b=1. That seems right.
# The PROBLEM is the asymmetry: H contributes phase = R*1 = R
# while the partner contributes R/n^1.95 which is tiny

print("Phase decomposition for all H-X bonds:")
print(f"  {'Mol':<5} {'R':>6} {'k_H':>6} {'k_X':>6} {'ph_H':>6} {'ph_X':>6} "
      f"{'ph_tot':>6} {'ph/pi':>6} {'|sin|':>6} {'De_exp':>6} {'err%':>6}")
print("-"*78)

for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    if 'H_1s' not in (o1, o2): continue

    h_orb = o1 if o1 == 'H_1s' else o2
    x_orb = o2 if o1 == 'H_1s' else o1

    n_h, l_h = 1, 0
    n_x, l_x = get_n(x_orb), get_l(x_orb)
    h_h, h_x = 0, has_nodes(x_orb)

    b_h = 1 + beta * h_h  # = 1
    b_x = 1 + beta * h_x
    k_h = 1.0 / n_h**b_h
    k_x = 1.0 / n_x**b_x

    ph_h = R * k_h
    ph_x = R * k_x
    ph = ph_h + ph_x

    a_h = 2 + (1-2*l_h)*alpha*h_h
    a_x = 2 + (1-2*l_x)*alpha*h_x
    E_scale = np.sqrt(E_H/n_h**a_h * E_H/n_x**a_x)
    D_cov = C_bond * E_scale * abs(np.sin(ph))

    # ionic
    eps_h = E_H * (Z_eff[h_orb]/n_h)**2
    eps_x = E_H * (Z_eff[x_orb]/n_x)**2
    dE = abs(eps_h - eps_x)
    V = max(abs(D_cov), 0.01)
    q = dE / np.sqrt(dE**2 + (2*V)**2) if dE > 0 else 0
    D_ion = c_ionic * q**2 * 2*E_H/R
    D_pred = D_cov + D_ion
    err = (D_pred - De_exp)/De_exp * 100

    print(f"  {name:<5} {R:6.3f} {k_h:6.3f} {k_x:6.3f} {ph_h:6.3f} {ph_x:6.3f} "
          f"{ph:6.3f} {ph/pi:6.3f} {abs(np.sin(ph)):6.4f} {De_exp:6.3f} {err:+5.1f}%")

print()
print("PATTERN: H's phase contribution (ph_H = R) dominates.")
print("When ph_H > pi (~3.14), the bond wraps regardless of partner.")
print(f"Critical R for wrapping: R > pi = {pi:.3f} bohr")
print()

# What if k_H isn't 1? What if there's a universal k = Z_eff/n^b?
print("="*80)
print("  TEST: What if k includes Z_eff?")
print("  k = Z_eff / n^b  (wavevector ~ momentum ~ Z/n)")
print("="*80)
print()
print("This makes physical sense: a tighter wavefunction (higher Z_eff)")
print("has a shorter wavelength (larger k), and the phase at distance R")
print("is k*R = (Z_eff/n^b)*R")
print()

print("Phase with k = Z_eff/n^b:")
print(f"  {'Mol':<5} {'R':>6} {'k1':>6} {'k2':>6} {'phase':>6} {'ph/pi':>6} "
      f"{'|sin|':>6} {'De_exp':>6} {'D_pred':>7} {'err%':>6}")
print("-"*78)

errs_zk = []
for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    n1, l1, n2, l2 = get_n(o1), get_l(o1), get_n(o2), get_l(o2)
    h1, h2 = has_nodes(o1), has_nodes(o2)

    b1 = 1 + beta*h1; b2 = 1 + beta*h2
    k1 = Z_eff[o1] / n1**b1
    k2 = Z_eff[o2] / n2**b2

    sigma_phase = R * (k1 + k2)

    a1 = 2 + (1-2*l1)*alpha*h1
    a2 = 2 + (1-2*l2)*alpha*h2
    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
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
    D_ion = c_ionic*q**2*2*E_H/R
    D_pred = D_cov + D_ion
    err = (D_pred - De_exp)/De_exp*100
    errs_zk.append(abs(err))

    flag = '***' if abs(err)<2 else ' **' if abs(err)<5 else '  *' if abs(err)<10 else ''
    print(f"  {name:<5} {R:6.3f} {k1:6.3f} {k2:6.3f} {sigma_phase:6.3f} {sigma_phase/pi:6.3f} "
          f"{abs(np.sin(sigma_phase)):6.4f} {De_exp:6.3f} {D_pred:7.3f} {err:+5.1f}% {flag}")

print(f"\n  avg={np.mean(errs_zk):.1f}%, med={np.median(errs_zk):.1f}%")
print(f"  within 5%: {sum(1 for e in errs_zk if e<5)}/24")
print(f"  within 10%: {sum(1 for e in errs_zk if e<10)}/24")


print()
print("="*80)
print("  TEST: k = 1/n^b (current formula, NO Z_eff in k)")
print("  vs k = Z_eff/n^b")
print("  vs k = sqrt(Z_eff)/n^b")
print("="*80)
print()

# Try different powers of Z_eff in k
for z_pow in [0, 0.25, 0.5, 0.75, 1.0]:
    errs = []
    for mol in molecules:
        name, R, De_exp, bonds, o1, o2 = mol
        n1, l1, n2, l2 = get_n(o1), get_l(o1), get_n(o2), get_l(o2)
        h1, h2 = has_nodes(o1), has_nodes(o2)

        b1 = 1 + beta*h1; b2 = 1 + beta*h2
        k1 = Z_eff[o1]**z_pow / n1**b1
        k2 = Z_eff[o2]**z_pow / n2**b2

        sigma_phase = R * (k1 + k2)

        a1 = 2 + (1-2*l1)*alpha*h1
        a2 = 2 + (1-2*l2)*alpha*h2
        E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)

        npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
        npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
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
        D_ion = c_ionic*q**2*2*E_H/R
        D_pred = D_cov + D_ion
        err = (D_pred - De_exp)/De_exp*100
        errs.append(abs(err))

    w5 = sum(1 for e in errs if e<5)
    w10 = sum(1 for e in errs if e<10)
    print(f"  Z_eff^{z_pow:.2f}: avg={np.mean(errs):5.1f}%, med={np.median(errs):5.1f}%, "
          f"w5={w5}/24, w10={w10}/24")


print()
print("="*80)
print("  CATEGORY 2: IONIC MONOPOLE TOO WEAK")
print("  LiF (-41%), NaCl (-47%)")
print("="*80)
print()

# These molecules are nearly 100% ionic (q ~ 1)
# The ionic formula: D_ionic = (1/7) * q^2 * 2*E_H/R
# For LiF: D_ionic = 1/7 * 1 * 2*13.606/2.955 = 1.312 eV
# Real De = 5.939, D_cov = 2.188, total = 3.500
# Missing: 5.939 - 3.500 = 2.439 eV

# What c_ionic would we need?
print("What c_ionic is needed to match experiment?")
print(f"  {'Mol':<6} {'De_exp':>7} {'D_cov':>7} {'D_ion_need':>9} {'c_need':>7} "
      f"{'q':>5} {'R':>6} {'c_need/c':>8}")
print("-"*65)

for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    n1, l1, n2, l2 = get_n(o1), get_l(o1), get_n(o2), get_l(o2)
    h1, h2 = has_nodes(o1), has_nodes(o2)

    a1 = 2 + (1-2*l1)*alpha*h1; a2 = 2 + (1-2*l2)*alpha*h2
    b1 = 1 + beta*h1; b2 = 1 + beta*h2
    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)
    sigma_phase = R/n1**b1 + R/n2**b2

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
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

    D_ion_needed = De_exp - D_cov
    if q > 0.01 and D_ion_needed > 0:
        c_needed = D_ion_needed / (q**2 * 2*E_H/R)
        if o1 != o2:  # heteronuclear only
            print(f"  {name:<6} {De_exp:7.3f} {D_cov:7.3f} {D_ion_needed:9.3f} "
                  f"{c_needed:7.4f} {q:5.3f} {R:6.3f} {c_needed/c_ionic:8.3f}")

print()
print(f"  Current c_ionic = 1/(2d+1) = {c_ionic:.4f}")
print()

# Is there a pattern in c_needed?
# Maybe ionic coupling scales with electronegativity difference?


print("="*80)
print("  CATEGORY 3: COVALENT OVERSHOOT")
print("  BF (+27%), CN (+31%)")
print("="*80)
print()

# These are TRIPLE bonds (BO=3) between atoms with very different Z_eff
# The formula gives too much covalent energy
# BF: B(2p, Ze=2.42) + F(2p, Ze=5.10) -- big asymmetry
# CN: C(2p, Ze=3.14) + N(2p, Ze=3.83) -- moderate asymmetry
# Compare with CO: C(2p) + O(2p) which works great at -1.2%
# And N2: N(2p) + N(2p) which works at +1.8%

print("Triple bond comparison:")
print(f"  {'Mol':<5} {'Ze1':>5} {'Ze2':>5} {'ratio':>6} {'De_exp':>7} {'D_pred':>7} {'err%':>6}")
print("-"*50)
for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    bo = sum(c for bt,c in bonds if 'anti' not in bt) - sum(c for bt,c in bonds if 'anti' in bt)
    if bo == 3:
        ze1 = Z_eff[o1]; ze2 = Z_eff[o2]
        ratio = max(ze1,ze2)/min(ze1,ze2)
        D_pred = 0  # recompute
        n1, l1, n2, l2 = get_n(o1), get_l(o1), get_n(o2), get_l(o2)
        h1, h2 = has_nodes(o1), has_nodes(o2)
        a1 = 2+(1-2*l1)*alpha*h1; a2 = 2+(1-2*l2)*alpha*h2
        b1 = 1+beta*h1; b2 = 1+beta*h2
        E_scale = np.sqrt(E_H/n1**a1*E_H/n2**a2)
        sigma_phase = R/n1**b1 + R/n2**b2

        npb = sum(c for bt,c in bonds if 'pi' in bt and 'anti' not in bt)
        npa = sum(c for bt,c in bonds if 'pi' in bt and 'anti' in bt)
        ifa = (npa >= npb) if npb > 0 else True
        D_cov = 0
        for bt,cnt in bonds:
            ph = sigma_phase if ('sigma' in bt or bt in ('ss','sp')) else sigma_phase*f_pi
            cont = C_bond*E_scale*abs(np.sin(ph))
            if 'anti' in bt:
                fa = 1.0 if ('sigma' in bt or ifa) else f_anti
                D_cov -= cnt*fa*cont
            else:
                D_cov += cnt*cont

        eps1 = E_H*(Z_eff[o1]/n1)**2; eps2 = E_H*(Z_eff[o2]/n2)**2
        dE = abs(eps1-eps2); V = max(abs(D_cov),0.01)
        q = dE/np.sqrt(dE**2+(2*V)**2) if dE>0 else 0
        D_ion = c_ionic*q**2*2*E_H/R
        D_pred = D_cov + D_ion
        err = (D_pred-De_exp)/De_exp*100

        print(f"  {name:<5} {ze1:5.2f} {ze2:5.2f} {ratio:6.3f} {De_exp:7.3f} {D_pred:7.3f} {err:+5.1f}%")

print()
print("  N2 and CO work, BF and CN overshoot.")
print("  BF has Z_eff ratio 2.11, CN has 1.22, CO has 1.42, N2 has 1.00")
print("  The overshoot doesn't correlate simply with Z_eff ratio.")


print()
print("="*80)
print("  DEEP DIVE: WHAT DOES EACH TERM CONTRIBUTE?")
print("="*80)
print()

print(f"  {'Mol':<6} {'E_scale':>7} {'phase':>6} {'|sin|':>6} {'D_sig':>7} {'D_pi':>7} "
      f"{'D_anti':>7} {'D_cov':>7} {'D_ion':>6} {'D_tot':>7} {'De_exp':>7} {'err':>6}")
print("-"*95)

for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    n1, l1, n2, l2 = get_n(o1), get_l(o1), get_n(o2), get_l(o2)
    h1, h2 = has_nodes(o1), has_nodes(o2)

    a1 = 2+(1-2*l1)*alpha*h1; a2 = 2+(1-2*l2)*alpha*h2
    b1 = 1+beta*h1; b2 = 1+beta*h2
    E_scale = np.sqrt(E_H/n1**a1*E_H/n2**a2)
    sigma_phase = R/n1**b1 + R/n2**b2

    npb = sum(c for bt,c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt,c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_sig = 0; D_pi_b = 0; D_anti = 0
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss','sp')) else sigma_phase*f_pi
        cont = C_bond*E_scale*abs(np.sin(ph))
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_anti += cnt*fa*cont
        elif 'pi' in bt:
            D_pi_b += cnt*cont
        else:
            D_sig += cnt*cont

    D_cov = D_sig + D_pi_b - D_anti

    eps1 = E_H*(Z_eff[o1]/n1)**2; eps2 = E_H*(Z_eff[o2]/n2)**2
    dE = abs(eps1-eps2); V = max(abs(D_cov),0.01)
    q = dE/np.sqrt(dE**2+(2*V)**2) if dE>0 else 0
    D_ion = c_ionic*q**2*2*E_H/R
    D_tot = D_cov + D_ion
    err = (D_tot-De_exp)/De_exp*100

    print(f"  {name:<6} {E_scale:7.3f} {sigma_phase:6.3f} {abs(np.sin(sigma_phase)):6.4f} "
          f"{D_sig:7.3f} {D_pi_b:7.3f} {D_anti:7.3f} {D_cov:7.3f} {D_ion:6.3f} "
          f"{D_tot:7.3f} {De_exp:7.3f} {err:+5.1f}%")


print()
print("="*80)
print("  KEY QUESTION: Is the energy exponent a = 2 + (1-2l)*alpha*h wrong?")
print("="*80)
print()

# The energy exponent determines E_scale = sqrt(E_H/n1^a1 * E_H/n2^a2)
# For s-orbitals with nodes (Li, Na): a = 2 + alpha = 2.7
#   E_scale = E_H / n^2.7
# For p-orbitals (B,C,N,O,F,Cl): a = 2 + (1-2)*alpha*1 = 2 - 0.7 = 1.3
#   E_scale = E_H / n^1.3  <-- much larger energy!
# For H (no nodes): a = 2, E_scale = E_H

# The (1-2l) factor means p-orbitals get LARGER E_scale than s-orbitals
# Is this physically right? p-orbitals have higher energy than s in the same n
# because they don't penetrate the core as much

print("Energy exponents and scales:")
orbs = ['H_1s', 'Li_2s', 'B_2p', 'C_2p', 'N_2p', 'O_2p', 'F_2p', 'Na_3s', 'Cl_3p']
print(f"  {'Orb':>6} {'n':>2} {'l':>2} {'h':>2} {'a':>5} {'b':>5} {'E_H/n^a':>8} {'E_real':>8}")
print("-"*50)
for orb in orbs:
    n, l = get_n(orb), get_l(orb)
    h = has_nodes(orb)
    a = 2 + (1-2*l)*alpha*h
    b = 1 + beta*h
    E_form = E_H / n**a
    E_real = E_H * (Z_eff[orb]/n)**2
    print(f"  {orb:>6} {n:2d} {l:2d} {h:2d} {a:5.2f} {b:5.2f} {E_form:8.4f} {E_real:8.4f}")

print()
print("  Note: the formula energy E_H/n^a is the 'free' energy scale")
print("  that doesn't use Z_eff. It's a UNIVERSAL function of (n, l).")
print("  The real orbital energy E_H*(Z_eff/n)^2 varies with atom.")
print()

# Check: does the formula energy ~ real energy for period 2?
print("  Ratio E_formula/E_real:")
for orb in orbs:
    n, l = get_n(orb), get_l(orb)
    h = has_nodes(orb)
    a = 2 + (1-2*l)*alpha*h
    E_form = E_H / n**a
    E_real = E_H * (Z_eff[orb]/n)**2
    ratio = E_form/E_real
    print(f"    {orb:>6}: E_form/E_real = {ratio:.4f}")


print()
print("="*80)
print("  WHAT IF: Use actual orbital energies instead of E_H/n^a?")
print("  E_scale = sqrt(E_real_1 * E_real_2) where E_real = E_H*(Z/n)^2")
print("="*80)
print()

errs_real_E = []
print(f"  {'Mol':<6} {'De_exp':>7} {'D_form':>7} {'D_realE':>7} {'err_f%':>7} {'err_rE%':>7}")
print("-"*55)

for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    n1, l1, n2, l2 = get_n(o1), get_l(o1), get_n(o2), get_l(o2)
    h1, h2 = has_nodes(o1), has_nodes(o2)

    # Formula energy
    a1 = 2+(1-2*l1)*alpha*h1; a2 = 2+(1-2*l2)*alpha*h2
    b1 = 1+beta*h1; b2 = 1+beta*h2
    E_scale_f = np.sqrt(E_H/n1**a1 * E_H/n2**a2)

    # Real orbital energy
    E1_r = E_H*(Z_eff[o1]/n1)**2; E2_r = E_H*(Z_eff[o2]/n2)**2
    E_scale_r = np.sqrt(E1_r * E2_r)

    sigma_phase = R/n1**b1 + R/n2**b2

    npb = sum(c for bt,c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt,c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    def calc_D(E_sc):
        Dc = 0
        for bt,cnt in bonds:
            ph = sigma_phase if ('sigma' in bt or bt in ('ss','sp')) else sigma_phase*f_pi
            cont = C_bond*E_sc*abs(np.sin(ph))
            if 'anti' in bt:
                fa = 1.0 if ('sigma' in bt or ifa) else f_anti
                Dc -= cnt*fa*cont
            else:
                Dc += cnt*cont
        eps1 = E_H*(Z_eff[o1]/n1)**2; eps2 = E_H*(Z_eff[o2]/n2)**2
        dE = abs(eps1-eps2); V = max(abs(Dc),0.01)
        q = dE/np.sqrt(dE**2+(2*V)**2) if dE>0 else 0
        Di = c_ionic*q**2*2*E_H/R
        return Dc+Di

    D_f = calc_D(E_scale_f)
    D_r = calc_D(E_scale_r)
    err_f = (D_f-De_exp)/De_exp*100
    err_r = (D_r-De_exp)/De_exp*100
    errs_real_E.append(abs(err_r))

    print(f"  {name:<6} {De_exp:7.3f} {D_f:7.3f} {D_r:7.3f} {err_f:+6.1f}% {err_r:+6.1f}%")

print(f"\n  Real E: avg={np.mean(errs_real_E):.1f}%, med={np.median(errs_real_E):.1f}%")
print(f"  within 5%: {sum(1 for e in errs_real_E if e<5)}/24")
print(f"  within 10%: {sum(1 for e in errs_real_E if e<10)}/24")


print()
print("="*80)
print("  WHAT IF: Phase uses Z_eff too?")
print("  k = Z_eff/n^b, E_scale still = E_H/n^a (formula)")
print("="*80)
print()

errs_zph = []
print(f"  {'Mol':<6} {'De_exp':>7} {'D_curr':>7} {'D_zph':>7} {'err_c%':>7} {'err_z%':>7}")
print("-"*55)

for mol in molecules:
    name, R, De_exp, bonds, o1, o2 = mol
    n1, l1, n2, l2 = get_n(o1), get_l(o1), get_n(o2), get_l(o2)
    h1, h2 = has_nodes(o1), has_nodes(o2)

    a1 = 2+(1-2*l1)*alpha*h1; a2 = 2+(1-2*l2)*alpha*h2
    b1 = 1+beta*h1; b2 = 1+beta*h2
    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)

    # Current phase: k = 1/n^b
    phase_curr = R/n1**b1 + R/n2**b2
    # Z_eff phase: k = Z_eff/n^b
    phase_z = R*Z_eff[o1]/n1**b1 + R*Z_eff[o2]/n2**b2

    def calc_D_ph(sigma_phase):
        npb = sum(c for bt,c in bonds if 'pi' in bt and 'anti' not in bt)
        npa = sum(c for bt,c in bonds if 'pi' in bt and 'anti' in bt)
        ifa = (npa >= npb) if npb > 0 else True
        Dc = 0
        for bt,cnt in bonds:
            ph = sigma_phase if ('sigma' in bt or bt in ('ss','sp')) else sigma_phase*f_pi
            cont = C_bond*E_scale*abs(np.sin(ph))
            if 'anti' in bt:
                fa = 1.0 if ('sigma' in bt or ifa) else f_anti
                Dc -= cnt*fa*cont
            else:
                Dc += cnt*cont
        eps1 = E_H*(Z_eff[o1]/n1)**2; eps2 = E_H*(Z_eff[o2]/n2)**2
        dE = abs(eps1-eps2); V = max(abs(Dc),0.01)
        q = dE/np.sqrt(dE**2+(2*V)**2) if dE>0 else 0
        Di = c_ionic*q**2*2*E_H/R
        return Dc+Di

    D_c = calc_D_ph(phase_curr)
    D_z = calc_D_ph(phase_z)
    err_c = (D_c-De_exp)/De_exp*100
    err_z = (D_z-De_exp)/De_exp*100
    errs_zph.append(abs(err_z))

    print(f"  {name:<6} {De_exp:7.3f} {D_c:7.3f} {D_z:7.3f} {err_c:+6.1f}% {err_z:+6.1f}%")

print(f"\n  Z_eff in phase: avg={np.mean(errs_zph):.1f}%, med={np.median(errs_zph):.1f}%")


print()
print("="*80)
print("  SUMMARY: What each modification does")
print("="*80)
print()
print(f"  Current formula (E_H/n^a, k=1/n^b):  avg=14.4%, med=4.1%, w5=15, w10=17")
print(f"  Z_eff in k only:                      avg={np.mean(errs_zk):.1f}%, med={np.median(errs_zk):.1f}%")
print(f"  Real E in E_scale only:               avg={np.mean(errs_real_E):.1f}%, med={np.median(errs_real_E):.1f}%")
print(f"  Z_eff in phase only:                  avg={np.mean(errs_zph):.1f}%, med={np.median(errs_zph):.1f}%")
