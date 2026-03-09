"""
DIRECT WAVE ENERGY PREDICTION
================================
Skip Z_eff entirely. Predict ionization potentials (IPs) directly
from wave geometry: Z, n, l, d=3.

IP = energy to remove the outermost standing wave quantum.
This is the REAL observable that drives chemical bonding.

In GWT: an atom with Z protons has Z energy quanta in standing
wave modes. The IP = energy of the outermost quantum.

For hydrogen (Z=1): IP = E_H = 13.606 eV (one quantum, one mode)

For multi-quantum atoms: the outermost mode's energy depends on
Z and the mode structure (n, l). The question: can we write
IP = E_H * f(Z, n, l, d)?
"""
import numpy as np

pi = np.pi
E_H = 13.6057
dd = 3
C_bond = pi / dd
f_pi = dd**2 / (dd**2 + 1)
alpha_n = 1 - f_pi / dd
beta_n = (1 + f_pi) / 2
f_anti = 2*dd / (2*dd - 1)
c_ionic = 1.0 / (2*dd + 1)

# Experimental first ionization potentials (eV)
atoms = {
    'H':  {'Z': 1,  'n': 1, 'l': 0, 'IP': 13.598},
    'He': {'Z': 2,  'n': 1, 'l': 0, 'IP': 24.587},
    'Li': {'Z': 3,  'n': 2, 'l': 0, 'IP': 5.392},
    'Be': {'Z': 4,  'n': 2, 'l': 0, 'IP': 9.323},
    'B':  {'Z': 5,  'n': 2, 'l': 1, 'IP': 8.298},
    'C':  {'Z': 6,  'n': 2, 'l': 1, 'IP': 11.260},
    'N':  {'Z': 7,  'n': 2, 'l': 1, 'IP': 14.534},
    'O':  {'Z': 8,  'n': 2, 'l': 1, 'IP': 13.618},
    'F':  {'Z': 9,  'n': 2, 'l': 1, 'IP': 17.423},
    'Ne': {'Z': 10, 'n': 2, 'l': 1, 'IP': 21.565},
    'Na': {'Z': 11, 'n': 3, 'l': 0, 'IP': 5.139},
    'Mg': {'Z': 12, 'n': 3, 'l': 0, 'IP': 7.646},
    'Al': {'Z': 13, 'n': 3, 'l': 1, 'IP': 5.986},
    'Si': {'Z': 14, 'n': 3, 'l': 1, 'IP': 8.152},
    'P':  {'Z': 15, 'n': 3, 'l': 1, 'IP': 10.487},
    'S':  {'Z': 16, 'n': 3, 'l': 1, 'IP': 10.360},
    'Cl': {'Z': 17, 'n': 3, 'l': 1, 'IP': 12.968},
    'Ar': {'Z': 18, 'n': 3, 'l': 1, 'IP': 15.760},
}

# =============================================================================
# LOOK AT THE PATTERN
# =============================================================================
print("=" * 80)
print("  IONIZATION POTENTIALS: The real energy differences")
print("=" * 80)

# Group by period and subshell
print("\n  Period 1 (1s):")
for name in ['H', 'He']:
    info = atoms[name]
    ratio = info['IP'] / E_H
    print(f"    {name:>2}: Z={info['Z']:2d}, IP={info['IP']:6.3f} eV, IP/E_H={ratio:.4f}")

print("\n  Period 2 (2s):")
for name in ['Li', 'Be']:
    info = atoms[name]
    ratio = info['IP'] / E_H
    print(f"    {name:>2}: Z={info['Z']:2d}, IP={info['IP']:6.3f} eV, IP/E_H={ratio:.4f}")

print("\n  Period 2 (2p):")
prev_ip = None
for name in ['B', 'C', 'N', 'O', 'F', 'Ne']:
    info = atoms[name]
    ratio = info['IP'] / E_H
    dip = info['IP'] - prev_ip if prev_ip else 0
    print(f"    {name:>2}: Z={info['Z']:2d}, IP={info['IP']:6.3f} eV, IP/E_H={ratio:.4f}, "
          f"dIP={dip:+6.3f}")
    prev_ip = info['IP']

print("\n  Period 3 (3s):")
for name in ['Na', 'Mg']:
    info = atoms[name]
    ratio = info['IP'] / E_H
    print(f"    {name:>2}: Z={info['Z']:2d}, IP={info['IP']:6.3f} eV, IP/E_H={ratio:.4f}")

print("\n  Period 3 (3p):")
prev_ip = None
for name in ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']:
    info = atoms[name]
    ratio = info['IP'] / E_H
    dip = info['IP'] - prev_ip if prev_ip else 0
    print(f"    {name:>2}: Z={info['Z']:2d}, IP={info['IP']:6.3f} eV, IP/E_H={ratio:.4f}, "
          f"dIP={dip:+6.3f}")
    prev_ip = info['IP']


# =============================================================================
# KEY PATTERNS
# =============================================================================
print()
print("=" * 80)
print("  KEY PATTERNS IN IP DATA")
print("=" * 80)

# 1. IP(He) / IP(H) = 24.587/13.598 = 1.808 ~ not 4 (Z^2)
#    Because of electron-electron repulsion
#    IP(He) = Z^2 * E_H - J = 4*13.606 - 29.837 = 24.587
#    So J(He) = 4*13.606 - 24.587 = 29.837... that's too large
#    Actually IP(He) = E(He+) - E(He) = Z^2*E_H - (Z^2*E_H - J_12) = J_12?
#    No: E(He) = 2*Z^2*E_H - J_12 = 2*4*13.606 - J = 108.848 - J
#    E(He+) = Z^2*E_H = 54.424
#    IP = E(He+) - E(He) = 54.424 - (108.848 - J) = J - 54.424
#    IP = 24.587, so J = 79.011... that's not right either.
#
#    Correct: E(He) = -(Z^2)*E_H/n^2 for each electron + repulsion
#    E(He) = -2*(4*13.606) + V_ee = -108.85 + V_ee
#    E(He+) = -4*13.606 = -54.42
#    IP = E(He+) - E(He) = -54.42 - (-108.85 + V_ee) = 54.42 - V_ee
#    24.59 = 54.42 - V_ee, so V_ee = 29.83 eV

# 2. Within a subshell, IP increases linearly (except pairing dips at O, S)
# The N->O dip is the HALF-FILLED SHELL effect
# The P->S dip is the same effect in period 3

# 3. The IP INCREMENT per Z within a subshell:
print("\n  IP increments per Z within 2p:")
ips_2p = [atoms[n]['IP'] for n in ['B', 'C', 'N', 'O', 'F', 'Ne']]
for i in range(1, len(ips_2p)):
    print(f"    Z={5+i}: dIP = {ips_2p[i]-ips_2p[i-1]:+6.3f} eV")

print("\n  IP increments per Z within 3p:")
ips_3p = [atoms[n]['IP'] for n in ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']]
for i in range(1, len(ips_3p)):
    print(f"    Z={13+i}: dIP = {ips_3p[i]-ips_3p[i-1]:+6.3f} eV")

# Average increment (excluding pairing dip):
# 2p: B->C: 2.96, C->N: 3.27, (N->O: -0.92 pairing), O->F: 3.81, F->Ne: 4.14
# Without pairing: avg = (2.96+3.27+3.81+4.14)/4 = 3.545
# 3p: similar pattern

print("\n  Average IP increment (no pairing dip):")
inc_2p_no_pair = [ips_2p[1]-ips_2p[0], ips_2p[2]-ips_2p[1],
                   ips_2p[4]-ips_2p[3], ips_2p[5]-ips_2p[4]]
inc_3p_no_pair = [ips_3p[1]-ips_3p[0], ips_3p[2]-ips_3p[1],
                   ips_3p[4]-ips_3p[3], ips_3p[5]-ips_3p[4]]
print(f"    2p: {np.mean(inc_2p_no_pair):.3f} eV")
print(f"    3p: {np.mean(inc_3p_no_pair):.3f} eV")
print(f"    Ratio: {np.mean(inc_2p_no_pair)/np.mean(inc_3p_no_pair):.3f}")
print(f"    Compare n2^2/n3^2 = 4/9 = {4/9:.3f}")
print(f"    Compare n2/n3 = 2/3 = {2/3:.3f}")

# What's the IP increment per Z in terms of E_H?
print(f"\n    2p increment / E_H = {np.mean(inc_2p_no_pair)/E_H:.4f}")
print(f"    3p increment / E_H = {np.mean(inc_3p_no_pair)/E_H:.4f}")
print(f"    Compare: 2/(d*n^2) for n=2: {2/(dd*4):.4f}")
print(f"    Compare: 2/(d*n^2) for n=3: {2/(dd*9):.4f}")


# =============================================================================
# PAIRING ENERGY
# =============================================================================
print()
print("=" * 80)
print("  PAIRING ENERGY: The half-filled shell effect")
print("=" * 80)

# When a p subshell goes from half-filled (l=1, 3 quanta) to half+1:
# the 4th quantum must pair with an existing one, costing energy

# N->O: IP drops by 0.916 eV relative to linear trend
# P->S: IP drops by 0.127 eV relative to linear trend

# Actually let me compute the pairing energy properly:
# Linear trend from unpaired increments, extrapolate
# 2p: increment = 3.545 eV
# Expected IP(O) = IP(N) + 3.545 = 14.534 + 3.545 = 18.079
# Actual IP(O) = 13.618
# Pairing energy = 18.079 - 13.618 = 4.461 eV

# 3p: increment = 2.481 eV
# Expected IP(S) = IP(P) + 2.481 = 10.487 + 2.481 = 12.968
# Actual IP(S) = 10.360
# Pairing energy = 12.968 - 10.360 = 2.608 eV

avg_inc_2p = np.mean(inc_2p_no_pair)
avg_inc_3p = np.mean(inc_3p_no_pair)
pair_2p = (atoms['N']['IP'] + avg_inc_2p) - atoms['O']['IP']
pair_3p = (atoms['P']['IP'] + avg_inc_3p) - atoms['S']['IP']

print(f"\n  2p pairing energy: {pair_2p:.3f} eV")
print(f"  3p pairing energy: {pair_3p:.3f} eV")
print(f"  Ratio: {pair_2p/pair_3p:.3f}")
print(f"  Compare: n3^2/n2^2 = 9/4 = {9/4:.3f}")

print(f"\n  2p pairing / E_H = {pair_2p/E_H:.4f}")
print(f"  3p pairing / E_H = {pair_3p/E_H:.4f}")
print(f"  Compare 1/d = {1/dd:.4f}")
print(f"  Compare 2/(d+1) = {2/(dd+1):.4f}")


# =============================================================================
# WAVE ENERGY MODEL
# =============================================================================
print()
print("=" * 80)
print("  WAVE ENERGY MODEL: IP from Z, n, l, d")
print("=" * 80)

# Model: IP = E_H/n^2 * (Z_eff_formula)
# where Z_eff_formula comes from the wave coupling law
#
# But we showed that Z_eff from wave coupling works for period 2.
# The issue is period 3.
#
# Alternative: model IP directly as a function of position in shell
#
# For the first electron in a new shell:
# IP(Li) = 5.39 = E_H/4 * Z_Li_eff^2
# Z_Li_eff = sqrt(5.39*4/13.606) = sqrt(1.585) = 1.259
# Compare: our harmonic model gives 1.286 (close!)
#
# For the first p electron:
# IP(B) = 8.30 = E_H/4 * Z_B_eff^2
# Z_B_eff = sqrt(8.30*4/13.606) = sqrt(2.441) = 1.562
# Compare: Clementi Z_eff = 2.421 (different because of Koopmans vs vertical IP)

# Actually, IP = -epsilon (Koopmans' theorem) = E_H * (Z_eff/n)^2
# But only approximately. The exact IP includes relaxation effects.

# Let me just use IPs directly for bond energy differences
# instead of going through Z_eff

# =============================================================================
# BOND TEST: Using IP differences directly
# =============================================================================
print()
print("=" * 80)
print("  BOND TEST: dE = |IP_1 - IP_2| for ionic correction")
print("=" * 80)

# Use only bond atoms
bond_atoms_IP = {
    'H': 13.598, 'Li': 5.392, 'B': 8.298, 'C': 11.260,
    'N': 14.534, 'O': 13.618, 'F': 17.423, 'Na': 5.139, 'Cl': 12.968,
}

bond_atoms_info = {
    'H':  {'n': 1, 'l': 0}, 'Li': {'n': 2, 'l': 0},
    'B':  {'n': 2, 'l': 1}, 'C':  {'n': 2, 'l': 1},
    'N':  {'n': 2, 'l': 1}, 'O':  {'n': 2, 'l': 1},
    'F':  {'n': 2, 'l': 1}, 'Na': {'n': 3, 'l': 0},
    'Cl': {'n': 3, 'l': 1},
}

# Clementi Z_eff for comparison
Z_eff_real = {
    'H': 1.0, 'Li': 1.2792, 'B': 2.4214, 'C': 3.1358, 'N': 3.8340,
    'O': 4.4532, 'F': 5.0998, 'Na': 2.5074, 'Cl': 6.1161
}

molecules = [
    ('H2',   1.401,  4.745, [('ss', 1)], 'H',  'H'),
    ('Li2',  5.051,  1.056, [('ss', 1)], 'Li', 'Li'),
    ('B2',   3.005,  3.02,  [('pi', 2)], 'B',  'B'),
    ('C2',   2.348,  6.32,  [('pi', 2), ('sp_sigma', 1), ('sp_sigma_anti', 1)], 'C',  'C'),
    ('N2',   2.074,  9.759, [('pp_sigma', 1), ('pi', 2)], 'N',  'N'),
    ('O2',   2.282,  5.213, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)], 'O',  'O'),
    ('F2',   2.668,  1.660, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)], 'F',  'F'),
    ('Na2',  5.818,  0.746, [('ss', 1)], 'Na', 'Na'),
    ('Cl2',  3.757,  2.514, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)], 'Cl', 'Cl'),
    ('HF',   1.733,  5.869, [('sp', 1)], 'H',  'F'),
    ('CO',   2.132, 11.225, [('pp_sigma', 1), ('pi', 2)], 'C',  'O'),
    ('NO',   2.175,  6.497, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)], 'N',  'O'),
    ('OH',   1.834,  4.392, [('sp', 1)], 'O',  'H'),
    ('HCl',  2.409,  4.434, [('sp', 1)], 'H',  'Cl'),
    ('LiH',  3.015,  2.515, [('ss', 1)], 'Li', 'H'),
    ('LiF',  2.955,  5.939, [('sp', 1)], 'Li', 'F'),
    ('BH',   2.329,  3.42,  [('sp', 1)], 'B',  'H'),
    ('CH',   2.116,  3.47,  [('sp', 1)], 'C',  'H'),
    ('NH',   1.958,  3.57,  [('sp', 1)], 'N',  'H'),
    ('BF',   2.386,  7.81,  [('pp_sigma', 1), ('pi', 2)], 'B',  'F'),
    ('CN',   2.214,  7.72,  [('pp_sigma', 1), ('pi', 2)], 'C',  'N'),
    ('NaH',  3.566,  1.97,  [('ss', 1)], 'Na', 'H'),
    ('NaCl', 4.461,  4.23,  [('sp', 1)], 'Na', 'Cl'),
    ('H2O',  1.809,  5.117, [('sp', 1)], 'O',  'H'),
]


def compute_bond(mol, dE_source='real'):
    """Compute bond energy with different dE sources."""
    name, R, De_exp, bonds, atom1, atom2 = mol
    info1 = bond_atoms_info[atom1]; info2 = bond_atoms_info[atom2]
    n1 = info1['n']; l1 = info1['l']; n2 = info2['n']; l2 = info2['l']
    h1 = min(n1 - l1 - 1, 1); h2 = min(n2 - l2 - 1, 1)

    a1 = 2 + (1 - 2*l1) * alpha_n * h1
    a2 = 2 + (1 - 2*l2) * alpha_n * h2
    b1 = 1 + beta_n * h1; b2 = 1 + beta_n * h2

    E_scale = np.sqrt(E_H / n1**a1 * E_H / n2**a2)
    k1 = 1.0 / n1**b1; k2 = 1.0 / n2**b2
    sigma_phase = R * (k1 + k2)

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss', 'sp')) else sigma_phase * f_pi
        cont = C_bond * E_scale * abs(np.sin(ph))
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else (2*dd/(2*dd-1))
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    # Different dE sources
    if dE_source == 'real':
        Ze1 = Z_eff_real[atom1]; Ze2 = Z_eff_real[atom2]
        dE = abs(E_H*(Ze1/n1)**2 - E_H*(Ze2/n2)**2)
    elif dE_source == 'IP':
        dE = abs(bond_atoms_IP[atom1] - bond_atoms_IP[atom2])
    elif dE_source == 'formula':
        dE = abs(E_H/n1**a1 - E_H/n2**a2)
    else:
        dE = 0

    V = max(abs(D_cov), 0.01)
    q = dE / np.sqrt(dE**2 + (2*V)**2) if dE > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / R
    return D_cov, D_cov + Di


# Compare all three dE sources
print()
print(f"{'Mol':<6} {'De_exp':>7} {'D_real':>7} {'D_IP':>7} {'D_form':>7} "
      f"{'err_r%':>7} {'err_IP%':>7} {'err_f%':>7}")
print("-" * 65)

errs_r = []; errs_ip = []; errs_f = []
for mol in molecules:
    _, D_r = compute_bond(mol, 'real')
    _, D_ip = compute_bond(mol, 'IP')
    _, D_f = compute_bond(mol, 'formula')
    De = mol[2]
    err_r = (D_r - De)/De*100
    err_ip = (D_ip - De)/De*100
    err_f = (D_f - De)/De*100
    errs_r.append(abs(err_r)); errs_ip.append(abs(err_ip)); errs_f.append(abs(err_f))
    print(f"{mol[0]:<6} {De:7.3f} {D_r:7.3f} {D_ip:7.3f} {D_f:7.3f} "
          f"{err_r:+6.1f}% {err_ip:+6.1f}% {err_f:+6.1f}%")

print(f"\n  Real Z_eff:  avg={np.mean(errs_r):.1f}%, med={np.median(errs_r):.1f}%, "
      f"<5%:{sum(1 for e in errs_r if e<5)}, <10%:{sum(1 for e in errs_r if e<10)}")
print(f"  IP-based:    avg={np.mean(errs_ip):.1f}%, med={np.median(errs_ip):.1f}%, "
      f"<5%:{sum(1 for e in errs_ip if e<5)}, <10%:{sum(1 for e in errs_ip if e<10)}")
print(f"  Formula:     avg={np.mean(errs_f):.1f}%, med={np.median(errs_f):.1f}%, "
      f"<5%:{sum(1 for e in errs_f if e<5)}, <10%:{sum(1 for e in errs_f if e<10)}")


# =============================================================================
# IP DIFFERENCES TABLE
# =============================================================================
print()
print("=" * 80)
print("  IP DIFFERENCES vs Z_eff ENERGY DIFFERENCES")
print("=" * 80)

pairs = set()
print(f"\n{'Pair':>8} {'dIP':>8} {'dE_Zeff':>8} {'dE_form':>8} {'dIP/dE_Z':>9}")
print("-" * 50)
for mol in molecules:
    a1, a2 = mol[4], mol[5]
    if a1 != a2:
        pair = tuple(sorted([a1, a2]))
        if pair not in pairs:
            pairs.add(pair)
            dIP = abs(bond_atoms_IP[pair[0]] - bond_atoms_IP[pair[1]])
            Ze1 = Z_eff_real[pair[0]]; Ze2 = Z_eff_real[pair[1]]
            n1 = bond_atoms_info[pair[0]]['n']; n2 = bond_atoms_info[pair[1]]['n']
            dE_Z = abs(E_H*(Ze1/n1)**2 - E_H*(Ze2/n2)**2)

            l1 = bond_atoms_info[pair[0]]['l']; l2 = bond_atoms_info[pair[1]]['l']
            h1 = min(n1-l1-1, 1); h2 = min(n2-l2-1, 1)
            a1_exp = 2 + (1-2*l1)*alpha_n*h1
            a2_exp = 2 + (1-2*l2)*alpha_n*h2
            dE_f = abs(E_H/n1**a1_exp - E_H/n2**a2_exp)

            ratio = dIP/dE_Z if dE_Z > 0.01 else float('inf')
            print(f"{pair[0]+'-'+pair[1]:>8} {dIP:8.3f} {dE_Z:8.3f} {dE_f:8.3f} {ratio:9.3f}")
