"""
DERIVING Z_eff FROM WAVE GEOMETRY
==================================
Goal: Find a mathematical formula that gives atom-specific energies
using only Z (atomic number), n, l, and geometric constants from d=3.

Z is just the number of standing wave modes in the atom.
Z_eff is how many of those modes the outermost electron "sees"
after inner modes screen the nuclear charge.

In wave terms: inner wave patterns partially cancel the nuclear
potential. The screening should follow from wave geometry.

Approach:
1. Map out Z vs Z_eff for all atoms we use
2. Look for patterns: Z_eff = f(Z, n, l) using GWT constants
3. Test whether bond predictions work with derived Z_eff
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

# Known Z_eff values (from Clementi-Raimondi, spectroscopic)
atoms = {
    'H':  {'Z': 1,  'n': 1, 'l': 0, 'Z_eff': 1.0000, 'orb': 'H_1s'},
    'Li': {'Z': 3,  'n': 2, 'l': 0, 'Z_eff': 1.2792, 'orb': 'Li_2s'},
    'B':  {'Z': 5,  'n': 2, 'l': 1, 'Z_eff': 2.4214, 'orb': 'B_2p'},
    'C':  {'Z': 6,  'n': 2, 'l': 1, 'Z_eff': 3.1358, 'orb': 'C_2p'},
    'N':  {'Z': 7,  'n': 2, 'l': 1, 'Z_eff': 3.8340, 'orb': 'N_2p'},
    'O':  {'Z': 8,  'n': 2, 'l': 1, 'Z_eff': 4.4532, 'orb': 'O_2p'},
    'F':  {'Z': 9,  'n': 2, 'l': 1, 'Z_eff': 5.0998, 'orb': 'F_2p'},
    'Na': {'Z': 11, 'n': 3, 'l': 0, 'Z_eff': 2.5074, 'orb': 'Na_3s'},
    'Cl': {'Z': 17, 'n': 3, 'l': 1, 'Z_eff': 6.1161, 'orb': 'Cl_3p'},
}

print("=" * 80)
print("  RAW DATA: Z vs Z_eff")
print("=" * 80)
print(f"{'Atom':>4} {'Z':>3} {'n':>2} {'l':>2} {'Z_eff':>7} {'Z-Z_eff':>8} {'Z_eff/Z':>8} "
      f"{'E_orb':>8} {'E_form':>8}")
print("-" * 65)

for name, info in atoms.items():
    Z = info['Z']; n = info['n']; l = info['l']; Ze = info['Z_eff']
    h = min(n - l - 1, 1)
    a = 2 + (1 - 2*l) * alpha_n * h
    E_orb = E_H * (Ze/n)**2
    E_form = E_H / n**a
    S = Z - Ze  # screening
    print(f"{name:>4} {Z:3d} {n:2d} {l:2d} {Ze:7.4f} {S:8.4f} {Ze/Z:8.4f} "
          f"{E_orb:8.4f} {E_form:8.4f}")


# =============================================================================
# PATTERN 1: Screening as function of Z, n, l
# =============================================================================
print()
print("=" * 80)
print("  SCREENING PATTERNS: S = Z - Z_eff")
print("=" * 80)

# Slater's rules (approximate):
# Same shell (n,l): screen by 0.35 each
# (n-1) shell: screen by 0.85 each
# (n-2) and below: screen by 1.00 each
# For 1s: other 1s screens by 0.30

# GWT approach: screening should come from wave geometry
# Inner wave modes (lower n) provide more screening than same-shell modes

# Let's look at the 2p series: B(Z=5) through F(Z=9)
print("\n2p series (B through F):")
print(f"{'Atom':>4} {'Z':>3} {'Z_eff':>7} {'S':>7} {'dZ_eff':>7} {'inner_e':>8} {'same_e':>8}")
print("-" * 55)

prev_Ze = None
for name in ['B', 'C', 'N', 'O', 'F']:
    info = atoms[name]
    Z = info['Z']; Ze = info['Z_eff']
    S = Z - Ze
    # For 2p: inner electrons = 1s^2 + 2s^2 = 4, same-shell = Z - 5 (for p electrons)
    inner_e = 4  # 1s^2 + 2s^2
    same_e = Z - 5  # other 2p electrons
    dZe = Ze - prev_Ze if prev_Ze else 0
    print(f"{name:>4} {Z:3d} {Ze:7.4f} {S:7.4f} {dZe:7.4f} {inner_e:8d} {same_e:8d}")
    prev_Ze = Ze

# Each additional proton adds +1 to Z, but screening by other 2p electrons
# reduces effective charge. The increment dZ_eff should tell us the screening per 2p electron.
print("\n  -> Each Z+1 adds ~0.67-0.72 to Z_eff (screens ~0.30 per same-shell)")

# =============================================================================
# PATTERN 2: Can we express Z_eff using GWT constants?
# =============================================================================
print()
print("=" * 80)
print("  TESTING Z_eff FORMULAS")
print("=" * 80)

# Formula candidates using d=3 constants:
# s_inner = screening per inner electron
# s_same = screening per same-shell electron

# Slater: s_inner ~ 0.85 (for n-1), s_same ~ 0.35
# GWT: what are the natural screening constants?

# Geometric candidates for s_same:
s_candidates = {
    '1/d = 1/3':       1/d,
    '1/pi':            1/pi,
    'alpha/d':         alpha_n/d,
    '1/(d+1) = 1/4':   1/(d+1),
    '1/(2d-1) = 1/5':  1/(2*d-1),
    'f_pi/d^2':        f_pi/d**2,
    '2/(2d+1) = 2/7':  2/(2*d+1),
}

# For inner shell screening, candidates:
s_inner_candidates = {
    '1-1/d = 2/3':     1 - 1/d,
    '1-1/pi':          1 - 1/pi,
    'f_pi = 9/10':     f_pi,
    '(d-1)/d = 2/3':   (d-1)/d,
    '1-1/(d+1) = 3/4': 1 - 1/(d+1),
    '1-1/(2d) = 5/6':  1 - 1/(2*d),
    'd/(d+1) = 3/4':   d/(d+1),
    '(2d-1)/(2d) = 5/6': (2*d-1)/(2*d),
}

# For each atom, compute: Z_eff_pred = Z - s_inner * N_inner - s_same * N_same
# Where N_inner = electrons in lower shells, N_same = other electrons in same (n,l) subshell

def count_electrons(name):
    """Return (N_inner, N_same) for valence orbital"""
    info = atoms[name]
    Z = info['Z']; n = info['n']; l = info['l']

    if name == 'H':
        return (0, 0)
    elif name == 'Li':
        # 2s: inner = 1s^2 = 2, same = 0 (first 2s electron)
        return (2, 0)
    elif name in ['B', 'C', 'N', 'O', 'F']:
        # 2p: inner = 1s^2 + 2s^2 = 4
        # same = other 2p electrons = Z - 5
        N_inner = 4
        N_same = Z - 5
        return (N_inner, N_same)
    elif name == 'Na':
        # 3s: inner = 1s^2 + 2s^2 + 2p^6 = 10, same = 0
        return (10, 0)
    elif name == 'Cl':
        # 3p: inner = 1s^2 + 2s^2 + 2p^6 + 3s^2 = 12
        # same = other 3p electrons = 17 - 13 = 4
        return (12, 4)

print(f"\n{'s_inner':>25} {'s_same':>25} {'rms_err':>8}")
print("-" * 65)

best = (999, '', '')
for si_name, si_val in s_inner_candidates.items():
    for ss_name, ss_val in s_candidates.items():
        errs = []
        for name, info in atoms.items():
            Z = info['Z']; Ze = info['Z_eff']
            N_inner, N_same = count_electrons(name)
            Ze_pred = Z - si_val * N_inner - ss_val * N_same
            errs.append((Ze_pred - Ze)**2)
        rms = np.sqrt(np.mean(errs))
        if rms < 0.5:  # only show good ones
            print(f"{si_name:>25} {ss_name:>25} {rms:8.4f}")
        if rms < best[0]:
            best = (rms, si_name, ss_name)

print(f"\n  Best: s_inner={best[1]}, s_same={best[2]}, RMS={best[0]:.4f}")


# =============================================================================
# DETAILED VIEW OF BEST SCREENING MODEL
# =============================================================================
# Use the best pair
si_best = s_inner_candidates[best[1]]
ss_best = s_candidates[best[2]]

print()
print("=" * 80)
print(f"  BEST SCREENING: s_inner={best[1]}={si_best:.6f}, s_same={best[2]}={ss_best:.6f}")
print("=" * 80)
print()

print(f"{'Atom':>4} {'Z':>3} {'N_in':>5} {'N_sm':>5} {'Z_eff':>7} {'Z_pred':>7} {'err':>7} "
      f"{'E_orb':>8} {'E_pred':>8} {'E_err%':>7}")
print("-" * 75)

for name, info in atoms.items():
    Z = info['Z']; n = info['n']; l = info['l']; Ze = info['Z_eff']
    N_inner, N_same = count_electrons(name)
    Ze_pred = Z - si_best * N_inner - ss_best * N_same

    E_orb = E_H * (Ze/n)**2
    E_pred = E_H * (Ze_pred/n)**2
    E_err = (E_pred - E_orb) / E_orb * 100

    print(f"{name:>4} {Z:3d} {N_inner:5d} {N_same:5d} {Ze:7.4f} {Ze_pred:7.4f} "
          f"{Ze_pred-Ze:7.4f} {E_orb:8.4f} {E_pred:8.4f} {E_err:+6.1f}%")


# =============================================================================
# TRY MORE NUANCED SCREENING: different s for different shell gaps
# =============================================================================
print()
print("=" * 80)
print("  NUANCED SCREENING: separate s for (n-1) vs (n-2) shells")
print("=" * 80)

# Slater distinguishes:
# s1 = screening by 1s electrons (for n>=2 valence)
# s2 = screening by (n-1) shell
# s3 = screening by same shell

# For 2p atoms: inner = 1s^2(deep) + 2s^2(n-1 same n, different l)
# Actually 2s screens 2p differently than 1s does

# Let's be more precise about electron configuration
def detailed_screening(name):
    """Return dict of electron groups and their counts"""
    configs = {
        'H':  {'1s': 0},
        'Li': {'1s': 2},
        'B':  {'1s': 2, '2s': 2, '2p': 0},
        'C':  {'1s': 2, '2s': 2, '2p': 1},
        'N':  {'1s': 2, '2s': 2, '2p': 2},
        'O':  {'1s': 2, '2s': 2, '2p': 3},
        'F':  {'1s': 2, '2s': 2, '2p': 4},
        'Na': {'1s': 2, '2s': 2, '2p': 6, '3s': 0},
        'Cl': {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 4},
    }
    return configs[name]

# Three screening constants: s_deep (n-2 and below), s_adj (n-1 or same n diff l), s_same
print("\nScanning 3-parameter model: Z_eff = Z - s_deep*N_deep - s_adj*N_adj - s_same*N_same")

# For valence orbital with quantum numbers (n, l):
# N_deep = electrons in shells n-2 and below
# N_adj = electrons in shell n-1, or same n but lower l
# N_same = other electrons in same (n, l) subshell

def get_groups(name):
    """Return (N_deep, N_adj, N_same) for the valence orbital"""
    info = atoms[name]
    n_val = info['n']; l_val = info['l']

    if name == 'H': return (0, 0, 0)
    elif name == 'Li': return (2, 0, 0)    # 1s^2 are n-1 for 2s
    elif name == 'B': return (2, 2, 0)     # 1s^2 deep, 2s^2 adj
    elif name == 'C': return (2, 2, 1)     # 1s^2 deep, 2s^2 adj, 1 other 2p
    elif name == 'N': return (2, 2, 2)
    elif name == 'O': return (2, 2, 3)
    elif name == 'F': return (2, 2, 4)
    elif name == 'Na': return (2, 8, 0)    # 1s^2 deep, 2s^2+2p^6 adj
    elif name == 'Cl': return (2+2+6, 2, 4)  # core deep, 3s adj, other 3p same

# Actually let me reconsider. For Li (2s valence):
# 1s^2 is n-1 shell, so adj
# For Na (3s valence):
# 1s^2 is n-2 (deep), 2s^2+2p^6 is n-1 (adj)
# For B-F (2p valence):
# 1s^2 is n-1 (deep? or adj?)
# 2s^2 is same n, lower l (adj)

# Let me use Slater's grouping more carefully:
# For 2p: (1s) screens as "n-1 and below", (2s,2p) screen as "same group"
# Slater groups: [1s] [2s,2p] [3s,3p] [3d] [4s,4p] ...
# Same group screens by 0.35, one group below by 0.85, two+ by 1.0

def slater_groups(name):
    """Return (N_below, N_same_group) using Slater's grouping"""
    if name == 'H':  return (0, 0)
    elif name == 'Li': return (0, 2)    # [2s,2p] group, but 1s is one group below?
    # Actually Slater: Li 2s -> same group = [2s,2p], groups below = [1s]
    # For Li: [1s]^2 below, [2s] has 0 others in same group
    elif name == 'Li': return (2, 0)
    elif name == 'B':  return (2, 2)    # [1s]^2 below, [2s,2p] group has 2s^2 = 2 others
    elif name == 'C':  return (2, 3)    # [1s]^2, [2s,2p] has 2s^2 + 1 2p = 3
    elif name == 'N':  return (2, 4)
    elif name == 'O':  return (2, 5)
    elif name == 'F':  return (2, 6)
    elif name == 'Na': return (10, 0)   # [1s]^2 + [2s,2p]^8 below, [3s,3p] has 0 others
    elif name == 'Cl': return (10, 6)   # below, [3s,3p] has 3s^2 + 4 3p = 6

# Scan 2 params: s_below and s_same_group
best2 = (999, 0, 0)
results = []
for s_below in np.arange(0.5, 1.05, 0.01):
    for s_same in np.arange(0.15, 0.55, 0.01):
        errs = []
        for name, info in atoms.items():
            Z = info['Z']; Ze = info['Z_eff']
            N_below, N_same = slater_groups(name)
            Ze_pred = Z - s_below * N_below - s_same * N_same
            errs.append((Ze_pred - Ze)**2)
        rms = np.sqrt(np.mean(errs))
        results.append((rms, s_below, s_same))
        if rms < best2[0]:
            best2 = (rms, s_below, s_same)

print(f"\n  Best 2-param Slater-style: s_below={best2[1]:.4f}, s_same={best2[2]:.4f}, RMS={best2[0]:.4f}")

# Check if best values align with GWT constants
sb = best2[1]
ss = best2[2]
print(f"\n  s_below = {sb:.4f}")
print(f"    compare: 1-1/(2d) = {1-1/(2*d):.4f}")
print(f"    compare: (2d-1)/(2d) = {(2*d-1)/(2*d):.4f}")
print(f"    compare: f_pi = {f_pi:.4f}")
print(f"    compare: d/(d+1) = {d/(d+1):.4f}")
print(f"    compare: 5/6 = {5/6:.4f}")
print(f"    compare: pi/4 = {pi/4:.4f}")

print(f"\n  s_same = {ss:.4f}")
print(f"    compare: 1/d = {1/d:.4f}")
print(f"    compare: 1/pi = {1/pi:.4f}")
print(f"    compare: 1/(d+1) = {1/(d+1):.4f}")
print(f"    compare: alpha = {alpha_n:.4f}")
print(f"    compare: 1/(2d-1) = {1/(2*d-1):.4f}")
print(f"    compare: 2/(2d+1) = {2/(2*d+1):.4f}")

# Show detailed predictions with best values
print()
print(f"{'Atom':>4} {'Z':>3} {'N_bel':>5} {'N_sam':>5} {'Z_eff':>7} {'Z_pred':>7} {'err':>7} "
      f"{'E_real':>8} {'E_pred':>8} {'E_err%':>7}")
print("-" * 75)

for name, info in atoms.items():
    Z = info['Z']; n = info['n']; Ze = info['Z_eff']
    N_below, N_same = slater_groups(name)
    Ze_pred = Z - sb * N_below - ss * N_same
    E_real = E_H * (Ze/n)**2
    E_pred = E_H * (Ze_pred/n)**2
    E_err = (E_pred - E_real) / E_real * 100
    print(f"{name:>4} {Z:3d} {N_below:5d} {N_same:5d} {Ze:7.4f} {Ze_pred:7.4f} "
          f"{Ze_pred-Ze:+7.4f} {E_real:8.4f} {E_pred:8.4f} {E_err:+6.1f}%")


# =============================================================================
# TEST: Use GWT-derived Z_eff in bond formula
# =============================================================================
print()
print("=" * 80)
print("  BOND TEST: Using GWT-derived Z_eff values")
print("=" * 80)

# Construct GWT Z_eff using best screening constants
def gwt_zeff(name):
    info = atoms[name]
    Z = info['Z']
    N_below, N_same = slater_groups(name)
    return Z - sb * N_below - ss * N_same

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

def compute_bond_gwt(mol):
    name, R, De_exp, bonds, atom1, atom2 = mol
    info1 = atoms[atom1]; info2 = atoms[atom2]
    n1 = info1['n']; l1 = info1['l']; n2 = info2['n']; l2 = info2['l']
    h1 = min(n1 - l1 - 1, 1); h2 = min(n2 - l2 - 1, 1)

    a1 = 2 + (1 - 2*l1) * alpha_n * h1
    a2 = 2 + (1 - 2*l2) * alpha_n * h2
    b1 = 1 + beta_n * h1
    b2 = 1 + beta_n * h2

    # Covalent part uses formula energy (no Z_eff)
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

    # Ionic part uses GWT-derived Z_eff
    Ze1 = gwt_zeff(atom1); Ze2 = gwt_zeff(atom2)
    E1_gwt = E_H * (Ze1/n1)**2
    E2_gwt = E_H * (Ze2/n2)**2
    dE_gwt = abs(E1_gwt - E2_gwt)

    V = max(abs(D_cov), 0.01)
    q = dE_gwt / np.sqrt(dE_gwt**2 + (2*V)**2) if dE_gwt > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / R
    D_total = D_cov + Di

    return D_cov, D_total, dE_gwt, Ze1, Ze2


# Also compute with real Z_eff for comparison
Z_eff_real = {
    'H': 1.0, 'Li': 1.2792, 'B': 2.4214, 'C': 3.1358, 'N': 3.8340,
    'O': 4.4532, 'F': 5.0998, 'Na': 2.5074, 'Cl': 6.1161
}

def compute_bond_real(mol):
    name, R, De_exp, bonds, atom1, atom2 = mol
    info1 = atoms[atom1]; info2 = atoms[atom2]
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
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    Ze1 = Z_eff_real[atom1]; Ze2 = Z_eff_real[atom2]
    E1 = E_H * (Ze1/n1)**2; E2 = E_H * (Ze2/n2)**2
    dE = abs(E1 - E2)
    V = max(abs(D_cov), 0.01)
    q = dE / np.sqrt(dE**2 + (2*V)**2) if dE > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / R
    return D_cov + Di

print()
print(f"{'Mol':<6} {'De_exp':>7} {'D_gwt':>7} {'D_real':>7} {'err_gwt%':>8} {'err_real%':>8} {'verdict':>8}")
print("-" * 65)

errs_gwt = []
errs_real = []
for mol in molecules:
    name = mol[0]; De_exp = mol[2]
    D_cov, D_gwt, dE_gwt, Ze1, Ze2 = compute_bond_gwt(mol)
    D_real = compute_bond_real(mol)

    err_gwt = (D_gwt - De_exp) / De_exp * 100
    err_real = (D_real - De_exp) / De_exp * 100
    errs_gwt.append(abs(err_gwt))
    errs_real.append(abs(err_real))

    v = 'BETTER' if abs(err_gwt) < abs(err_real) - 1 else 'WORSE' if abs(err_gwt) > abs(err_real) + 1 else '~'
    print(f"{name:<6} {De_exp:7.3f} {D_gwt:7.3f} {D_real:7.3f} {err_gwt:+7.1f}% {err_real:+7.1f}% {v}")

print(f"\n  GWT-derived Z_eff:  avg={np.mean(errs_gwt):.1f}%, med={np.median(errs_gwt):.1f}%")
print(f"  Real Z_eff:         avg={np.mean(errs_real):.1f}%, med={np.median(errs_real):.1f}%")

# Show GWT Z_eff vs real
print()
print("GWT-derived vs real Z_eff:")
print(f"{'Atom':>4} {'Z_real':>7} {'Z_gwt':>7} {'err':>7}")
print("-" * 30)
for name in atoms:
    Ze_real = Z_eff_real[name]
    Ze_gwt = gwt_zeff(name)
    print(f"{name:>4} {Ze_real:7.4f} {Ze_gwt:7.4f} {Ze_gwt-Ze_real:+7.4f}")
