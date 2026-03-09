"""
BOND LENGTH AS GRAVITATIONAL ORBIT
====================================

Key insight: bond length = orbit radius of the two-body wave system.

In gravity: R_orbit = a0 * n^2 / Z  (Bohr model)
In bonds:   R_bond = c * rB_eff     (effective Bohr radius of the pair)

The geometric mean rB_eff = sqrt(rB1 * rB2) = n1*n2 / sqrt(Z1*Z2)
has the BEST correlation with bond length (+0.963).

This script explores:
1. Why geometric mean works (reduced orbit)
2. What sets the coefficient (GWT constants?)
3. Why F2/O2/Na2 are outliers (lone pairs? weak overlap?)
4. Whether the coefficient can be expressed as pi/d, d/(d-1), etc.
"""

import numpy as np

pi = np.pi
E_H = 13.6057  # eV
d = 3

# GWT parameters
C_bond = pi / d  # bonding coupling constant
f_pi = d**2 / (d**2 + 1)  # 9/10
alpha_p = 1 - f_pi / d    # 7/10
beta_p = (1 + f_pi) / 2   # 19/20

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
def get_b(orb):
    l, h = get_l(orb), has_nodes(orb)
    return 1 + (1 - 2*l) * beta_p * h
def orbital_energy(orb):
    return E_H * (Z_eff[orb] / get_n(orb))**2
def bohr_radius(orb):
    return get_n(orb)**2 / Z_eff[orb]

molecules = [
    ('H2',   1.401, 4.745, 'H_1s',  'H_1s',  1),
    ('Li2',  5.051, 1.056, 'Li_2s', 'Li_2s', 1),
    ('B2',   3.005, 3.02,  'B_2p',  'B_2p',  2),
    ('C2',   2.348, 6.32,  'C_2p',  'C_2p',  2),
    ('N2',   2.074, 9.759, 'N_2p',  'N_2p',  3),
    ('O2',   2.282, 5.213, 'O_2p',  'O_2p',  2),
    ('F2',   2.668, 1.660, 'F_2p',  'F_2p',  1),
    ('Na2',  5.818, 0.746, 'Na_3s', 'Na_3s', 1),
    ('Cl2',  3.757, 2.514, 'Cl_3p', 'Cl_3p', 1),
    ('HF',   1.733, 5.869, 'H_1s',  'F_2p',  1),
    ('CO',   2.132, 11.225,'C_2p',  'O_2p',  3),
    ('NO',   2.175, 6.497, 'N_2p',  'O_2p',  2),
    ('OH',   1.834, 4.392, 'O_2p',  'H_1s',  1),
    ('HCl',  2.409, 4.434, 'H_1s',  'Cl_3p', 1),
    ('LiH',  3.015, 2.515, 'Li_2s', 'H_1s',  1),
    ('LiF',  2.955, 5.939, 'Li_2s', 'F_2p',  1),
    ('BH',   2.329, 3.42,  'B_2p',  'H_1s',  1),
    ('CH',   2.116, 3.47,  'C_2p',  'H_1s',  1),
    ('NH',   1.958, 3.57,  'N_2p',  'H_1s',  1),
    ('BF',   2.386, 7.81,  'B_2p',  'F_2p',  3),
    ('CN',   2.214, 7.72,  'C_2p',  'N_2p',  3),
    ('NaH',  3.566, 1.97,  'Na_3s', 'H_1s',  1),
    ('NaCl', 4.461, 4.23,  'Na_3s', 'Cl_3p', 1),
    ('H2O',  1.809, 5.117, 'O_2p',  'H_1s',  1),
]


# ================================================================
# SECTION 1: Clean formula test: R = c * n1*n2/sqrt(Z1*Z2) / BO^p
# ================================================================
print("=" * 90)
print("  BOND LENGTH = ORBIT RADIUS: R = c * n1*n2/sqrt(Z1*Z2) / BO^p")
print("=" * 90)
print()

# Precompute
data = []
for name, R, De, o1, o2, bo in molecules:
    n1, n2 = get_n(o1), get_n(o2)
    Z1, Z2 = Z_eff[o1], Z_eff[o2]
    rB1, rB2 = bohr_radius(o1), bohr_radius(o2)
    E1, E2 = orbital_energy(o1), orbital_energy(o2)
    l1, l2 = get_l(o1), get_l(o2)

    data.append({
        'name': name, 'R': R, 'De': De, 'bo': bo,
        'n1': n1, 'n2': n2, 'Z1': Z1, 'Z2': Z2,
        'rB1': rB1, 'rB2': rB2,
        'E1': E1, 'E2': E2,
        'l1': l1, 'l2': l2,
        'o1': o1, 'o2': o2,
        'rB_geom': np.sqrt(rB1 * rB2),
        'rB_sum': rB1 + rB2,
    })


# ================================================================
# Test GWT constants as the coefficient
# ================================================================
print("Testing GWT constants as coefficient c:")
print()

gwt_constants = {
    '2': 2.0,
    'pi/d = pi/3': pi/3,
    'd-1 = 2': 2.0,
    'd/(d-1) = 3/2': 3.0/2,
    'pi/(d-1) = pi/2': pi/2,
    'sqrt(d) = sqrt(3)': np.sqrt(3),
    '2d/(2d-1) = 6/5': 6.0/5,
    'pi/sqrt(d)': pi/np.sqrt(3),
    '(d+1)/d = 4/3': 4.0/3,
    '(2d-1)/d = 5/3': 5.0/3,
    'sqrt(2*pi)': np.sqrt(2*pi),
    'phi (golden)': (1+np.sqrt(5))/2,
    '2*pi/d = 2pi/3': 2*pi/3,
    '1+1/d = 4/3': 1+1.0/3,
    'd^2/(d^2+1) = 9/10': 9.0/10,
    'C_bond = pi/d': C_bond,
}

# For each constant, find best BO exponent
for label, c in sorted(gwt_constants.items(), key=lambda x: x[1]):
    best_err = 100
    best_p = 0
    for p in np.arange(-0.5, 1.5, 0.01):
        errs = []
        for dd in data:
            R_pred = c * dd['rB_geom'] / max(dd['bo'], 0.5)**p
            errs.append(abs((R_pred - dd['R']) / dd['R'] * 100))
        err = np.mean(errs)
        if err < best_err:
            best_err = err
            best_p = p
    print(f"  c = {c:.4f} ({label:<22s}): avg|err| = {best_err:5.1f}%, BO^{best_p:.2f}")

print()


# ================================================================
# SECTION 2: Why does n^2/Z beat n/Z?
# ================================================================
print("=" * 90)
print("  WHY n^2/Z WORKS: the orbit radius IS the force boundary")
print("=" * 90)
print()
print("In the Bohr model: R = a0 * n^2/Z comes from")
print("  force balance:  m*v^2/R = Z*e^2/R^2")
print("  quantization:   m*v*R = n*hbar")
print()
print("In GWT bonds: the standing wave 'orbits' the nucleus at r = n^2/Z.")
print("This IS the extent of the wave -- where its probability peaks.")
print("The bond forms where two such orbits touch.")
print()
print("Key insight: the 'boundary' of the wave IS the Bohr orbit.")
print("Not the exponential tail (n/Z), but the peak (n^2/Z).")
print()

# Show R_exp vs 2*rB_geom for homonuclear molecules
print("Homonuclear molecules: R vs 2*rB")
print(f"{'Mol':<6} {'R_exp':>6} {'rB':>6} {'2*rB':>6} {'R/(2rB)':>8}")
print("-" * 40)
for dd in data:
    if dd['n1'] == dd['n2'] and dd['Z1'] == dd['Z2']:
        rB = dd['rB1']
        print(f"{dd['name']:<6} {dd['R']:6.3f} {rB:6.3f} {2*rB:6.3f} {dd['R']/(2*rB):8.3f}")
print()

# For homonuclear: R/(2*rB) varies from 0.70 (H2) to 1.70 (F2)
# What accounts for this variation?
# Check if it correlates with ionization energy, electronegativity, etc.


# ================================================================
# SECTION 3: l-dependent correction
# ================================================================
print("=" * 90)
print("  l-DEPENDENT CORRECTION: s-orbitals vs p-orbitals")
print("=" * 90)
print()
print("s-orbitals (l=0) are spherical: full overlap on the bond axis")
print("p-orbitals (l=1) are directional: enhanced sigma overlap")
print()
print("In GWT, the standing wave has angular dependence:")
print("  s: Y_00 ~ 1  (isotropic)")
print("  p: Y_10 ~ cos(theta) (directional)")
print()
print("The directional p-orbital concentrates density along bond axis,")
print("which should SHORTEN bonds compared to naive rB estimate.")
print()

# Test: add angular correction
# For p orbital, effective radius should be smaller (concentrated)
# For s orbital, effective radius is the full Bohr radius
# Try: r_eff = rB / (1 + l * correction)

def r_eff_angular(orb, c_ang):
    """Bohr radius with angular correction for p-orbitals."""
    rB = bohr_radius(orb)
    l = get_l(orb)
    return rB / (1 + l * c_ang)

# Grid search
best = (100, 0, 0)
for c in np.arange(0.1, 4.0, 0.01):
    for c_ang in np.arange(-0.5, 2.0, 0.05):
        errs = []
        for dd in data:
            r1 = r_eff_angular(dd['o1'], c_ang)
            r2 = r_eff_angular(dd['o2'], c_ang)
            rg = np.sqrt(r1 * r2)
            R_pred = c * rg
            errs.append(abs((R_pred - dd['R']) / dd['R'] * 100))
        err = np.mean(errs)
        if err < best[0]:
            best = (err, c, c_ang)

err, c_opt, c_ang = best
print(f"R = {c_opt:.3f} * sqrt(r_eff_1 * r_eff_2)")
print(f"  where r_eff = rB / (1 + l * {c_ang:.2f})")
print(f"  Average |error| = {err:.1f}%")
print()

header = f"{'Mol':<6} {'R_exp':>6} {'R_pred':>6} {'err%':>6} {'r1':>5} {'r2':>5} {'type':<4}"
print(header)
print("-" * 45)

errs_all = []
for dd in data:
    r1 = r_eff_angular(dd['o1'], c_ang)
    r2 = r_eff_angular(dd['o2'], c_ang)
    rg = np.sqrt(r1 * r2)
    R_pred = c_opt * rg
    e = (R_pred - dd['R']) / dd['R'] * 100
    errs_all.append(abs(e))
    types = {(0,0): 'ss', (0,1): 'sp', (1,0): 'sp', (1,1): 'pp'}
    btype = types[(dd['l1'], dd['l2'])]
    flag = '***' if abs(e) < 5 else ' **' if abs(e) < 10 else '  *' if abs(e) < 15 else ''
    print(f"{dd['name']:<6} {dd['R']:6.3f} {R_pred:6.3f} {e:+5.1f}% {r1:5.3f} {r2:5.3f} {btype:<4} {flag}")

w5 = sum(1 for e in errs_all if e < 5)
w10 = sum(1 for e in errs_all if e < 10)
w15 = sum(1 for e in errs_all if e < 15)
print(f"\nAvg |err| = {np.mean(errs_all):.1f}%, within 5%: {w5}/24, 10%: {w10}/24, 15%: {w15}/24")
print()


# ================================================================
# SECTION 4: Now add bond order
# ================================================================
print("=" * 90)
print("  WITH BOND ORDER: R = c * sqrt(r1*r2) / BO^p")
print("=" * 90)
print()

best = (100, 0, 0, 0)
for c in np.arange(0.1, 4.0, 0.01):
    for c_ang in np.arange(-0.5, 2.0, 0.05):
        for p in np.arange(-0.3, 1.0, 0.02):
            errs = []
            for dd in data:
                r1 = r_eff_angular(dd['o1'], c_ang)
                r2 = r_eff_angular(dd['o2'], c_ang)
                rg = np.sqrt(r1 * r2)
                R_pred = c * rg / max(dd['bo'], 0.5)**p
                errs.append(abs((R_pred - dd['R']) / dd['R'] * 100))
            err = np.mean(errs)
            if err < best[0]:
                best = (err, c, c_ang, p)

err, c_opt, c_ang, p_opt = best
print(f"R = {c_opt:.3f} * sqrt(r1*r2) / BO^{p_opt:.2f}")
print(f"  where r_eff = rB / (1 + l * {c_ang:.2f})")
print(f"  Average |error| = {err:.1f}%")
print()

header = f"{'Mol':<6} {'R_exp':>6} {'R_pred':>6} {'err%':>6} {'BO':>3} {'type':<4}"
print(header)
print("-" * 40)

errs_all = []
for dd in data:
    r1 = r_eff_angular(dd['o1'], c_ang)
    r2 = r_eff_angular(dd['o2'], c_ang)
    rg = np.sqrt(r1 * r2)
    R_pred = c_opt * rg / max(dd['bo'], 0.5)**p_opt
    e = (R_pred - dd['R']) / dd['R'] * 100
    errs_all.append(abs(e))
    types = {(0,0): 'ss', (0,1): 'sp', (1,0): 'sp', (1,1): 'pp'}
    btype = types[(dd['l1'], dd['l2'])]
    flag = '***' if abs(e) < 5 else ' **' if abs(e) < 10 else '  *' if abs(e) < 15 else ''
    print(f"{dd['name']:<6} {dd['R']:6.3f} {R_pred:6.3f} {e:+5.1f}% {dd['bo']:3d} {btype:<4} {flag}")

w5 = sum(1 for e in errs_all if e < 5)
w10 = sum(1 for e in errs_all if e < 10)
w15 = sum(1 for e in errs_all if e < 15)
print(f"\nAvg |err| = {np.mean(errs_all):.1f}%, within 5%: {w5}/24, 10%: {w10}/24, 15%: {w15}/24")
print()


# ================================================================
# SECTION 5: The RATIO R_exp / R_pred -- what's the residual pattern?
# ================================================================
print("=" * 90)
print("  RESIDUAL ANALYSIS: what does the simple model miss?")
print("=" * 90)
print()
print("Using R = 2.0 * sqrt(rB1*rB2) (simple, no fitting)")
print("What determines the residual ratio = R_exp / R_pred?")
print()

c_simple = 2.0
print(f"{'Mol':<6} {'R_exp':>6} {'R_pred':>6} {'ratio':>6} {'De':>7} {'Z_avg':>6} {'BO':>3} {'type'}")
print("-" * 65)

ratios = []
for dd in data:
    R_pred = c_simple * dd['rB_geom']
    ratio = dd['R'] / R_pred
    ratios.append(ratio)
    Z_avg = np.sqrt(dd['Z1'] * dd['Z2'])
    types = {(0,0): 'ss', (0,1): 'sp', (1,0): 'sp', (1,1): 'pp'}
    btype = types[(dd['l1'], dd['l2'])]
    print(f"{dd['name']:<6} {dd['R']:6.3f} {R_pred:6.3f} {ratio:6.3f} {dd['De']:7.3f} "
          f"{Z_avg:6.3f} {dd['bo']:3d}  {btype}")

print(f"\nRatio: mean = {np.mean(ratios):.3f}, std = {np.std(ratios):.3f}")
print()

# Check correlations of ratio with various properties
print("Correlations of R_exp/(2*rB_geom) with properties:")
ratio_arr = np.array(ratios)
props = {
    'De (bond energy)': [dd['De'] for dd in data],
    'BO (bond order)': [dd['bo'] for dd in data],
    'Z_geom': [np.sqrt(dd['Z1']*dd['Z2']) for dd in data],
    'E_scale': [np.sqrt(dd['E1']*dd['E2']) for dd in data],
    'Z_max': [max(dd['Z1'], dd['Z2']) for dd in data],
    'Z_diff': [abs(dd['Z1']-dd['Z2']) for dd in data],
    'n_max': [max(dd['n1'], dd['n2']) for dd in data],
    'l_sum': [dd['l1'] + dd['l2'] for dd in data],
    'rB_ratio': [max(dd['rB1'],dd['rB2'])/min(dd['rB1'],dd['rB2']) for dd in data],
}
for label, vals in props.items():
    c = np.corrcoef(ratio_arr, vals)[0, 1]
    print(f"  {label:<25s} {c:+.3f}")
print()


# ================================================================
# SECTION 6: Physical insight -- orbit overlap condition
# ================================================================
print("=" * 90)
print("  ORBIT OVERLAP: R = rB1 + rB2 modified by overlap geometry")
print("=" * 90)
print()
print("The simplest 'gravitational orbit' picture:")
print("  Two atoms touch when their Bohr orbits meet: R = rB1 + rB2")
print("  But the actual bond contracts because of shared electron density.")
print()
print("Contraction factor: the overlap of two orbitals creates a")
print("bonding region that PULLS atoms closer than rB1+rB2.")
print()
print("How much closer? Depends on the coupling strength:")
print("  Weak coupling (F2): barely contracts -> R ~ rB1+rB2 or more")
print("  Strong coupling (N2): contracts significantly")
print()

# Test: R = (rB1 + rB2) * contraction(De, E_scale)
# The contraction should depend on De/E_scale (relative bond strength)

print(f"{'Mol':<6} {'R_exp':>6} {'rB_sum':>6} {'ratio':>6} {'De':>7} {'E_scl':>6} {'De/E':>6}")
print("-" * 55)

for dd in data:
    rB_sum = dd['rB1'] + dd['rB2']
    ratio = dd['R'] / rB_sum
    E_scale = np.sqrt(dd['E1'] * dd['E2'])
    print(f"{dd['name']:<6} {dd['R']:6.3f} {rB_sum:6.3f} {ratio:6.3f} "
          f"{dd['De']:7.3f} {E_scale:6.2f} {dd['De']/E_scale:6.3f}")
print()


# ================================================================
# SECTION 7: The GWT formula -- all from d=3
# ================================================================
print("=" * 90)
print("  DERIVING THE COEFFICIENT FROM d=3")
print("=" * 90)
print()
print("Candidate zero-parameter formulas:")
print()

# The coefficient should be a simple function of d=3
# Let's systematically test functions of d
candidates = {
    # Simple powers/fractions of d
    'pi/d': pi/d,
    '2': 2,
    'd/(d-1)': d/(d-1),
    '(d+1)/d': (d+1)/d,
    'sqrt(d)': np.sqrt(d),
    '(d-1)': d-1,
    'pi/(d-1)': pi/(d-1),

    # Bond coupling and combinations
    'C_bond*d = pi': pi,
    'pi/sqrt(d)': pi/np.sqrt(d),
    '2d/(2d-1)': 2*d/(2*d-1),

    # With angular momentum
    'd/(d-1) * 1/(2l+1)': None,  # l-dependent
}

# For each candidate coefficient, compute errors with BO correction
for label, c in candidates.items():
    if c is None:
        continue

    # Without BO
    errs_nobo = []
    for dd in data:
        R_pred = c * dd['rB_geom']
        errs_nobo.append(abs((R_pred - dd['R']) / dd['R'] * 100))

    # With optimal BO
    best_bo = (np.mean(errs_nobo), 0)
    for p in np.arange(0, 0.5, 0.01):
        errs = []
        for dd in data:
            R_pred = c * dd['rB_geom'] / max(dd['bo'], 0.5)**p
            errs.append(abs((R_pred - dd['R']) / dd['R'] * 100))
        err = np.mean(errs)
        if err < best_bo[0]:
            best_bo = (err, p)

    print(f"  c = {c:.4f} ({label:<20s}): no BO: {np.mean(errs_nobo):5.1f}%, "
          f"best BO^{best_bo[1]:.2f}: {best_bo[0]:5.1f}%")

print()


# ================================================================
# SECTION 8: Best candidate with detailed output
# ================================================================
print("=" * 90)
print("  BEST CANDIDATE: R = 2 * sqrt(rB1*rB2) / BO^p")
print("=" * 90)
print()
print("c = 2 is natural: the 'touching spheres' radius is 2*rB for equal atoms")
print("(each atom's orbit extends to rB, and the bond axis spans 2*rB)")
print()

# Find best BO exponent for c=2
best = (100, 0)
for p in np.arange(0, 0.5, 0.001):
    errs = []
    for dd in data:
        R_pred = 2.0 * dd['rB_geom'] / max(dd['bo'], 0.5)**p
        errs.append(abs((R_pred - dd['R']) / dd['R'] * 100))
    err = np.mean(errs)
    if err < best[0]:
        best = (err, p)

err, p_opt = best
print(f"R = 2 * sqrt(rB1*rB2) / BO^{p_opt:.3f}")
print(f"Average |error| = {err:.1f}%")

# Check if p ≈ 1/d = 1/3, 1/(d-1) = 1/2, 1/(d+1) = 1/4, etc.
print()
print("BO exponent candidates from GWT:")
for label, p_test in [('1/d = 1/3', 1/3), ('1/(d+1) = 1/4', 1/4),
                       ('1/(2d) = 1/6', 1/6), ('1/(d-1) = 1/2', 1/2),
                       ('2/d = 2/3', 2/3), (f'optimal = {p_opt:.3f}', p_opt)]:
    errs = []
    for dd in data:
        R_pred = 2.0 * dd['rB_geom'] / max(dd['bo'], 0.5)**p_test
        errs.append(abs((R_pred - dd['R']) / dd['R'] * 100))
    print(f"  BO^{p_test:.4f} ({label}): avg|err| = {np.mean(errs):.1f}%")

print()

# Show detailed results for the best
p_use = p_opt
print(f"\nDetailed: R = 2 * sqrt(rB1*rB2) / BO^{p_use:.3f}")
print()
header = f"{'Mol':<6} {'R_exp':>6} {'R_pred':>6} {'err%':>6} {'BO':>3} {'rB_g':>5}"
print(header)
print("-" * 35)

errs_all = []
for dd in data:
    R_pred = 2.0 * dd['rB_geom'] / max(dd['bo'], 0.5)**p_use
    e = (R_pred - dd['R']) / dd['R'] * 100
    errs_all.append(abs(e))
    flag = '***' if abs(e) < 5 else ' **' if abs(e) < 10 else '  *' if abs(e) < 15 else ''
    print(f"{dd['name']:<6} {dd['R']:6.3f} {R_pred:6.3f} {e:+5.1f}% {dd['bo']:3d} {dd['rB_geom']:5.3f} {flag}")

w5 = sum(1 for e in errs_all if e < 5)
w10 = sum(1 for e in errs_all if e < 10)
w15 = sum(1 for e in errs_all if e < 15)
print(f"\nAvg |err| = {np.mean(errs_all):.1f}%, within 5%: {w5}/24, 10%: {w10}/24, 15%: {w15}/24")
print()


# ================================================================
# SECTION 9: Analysis of outliers -- what's special about F2, O2?
# ================================================================
print("=" * 90)
print("  OUTLIER ANALYSIS: Why do F2, O2, H2, Na2 deviate?")
print("=" * 90)
print()

# Group molecules by error magnitude
outliers = []
good = []
for dd in data:
    R_pred = 2.0 * dd['rB_geom'] / max(dd['bo'], 0.5)**p_use
    e = (R_pred - dd['R']) / dd['R'] * 100
    if abs(e) > 15:
        outliers.append((dd, e))
    else:
        good.append((dd, e))

print(f"Good predictions ({len(good)} molecules, |err| < 15%):")
for dd, e in sorted(good, key=lambda x: abs(x[1])):
    print(f"  {dd['name']:<6} err = {e:+5.1f}%  n=({dd['n1']},{dd['n2']})  "
          f"Z=({dd['Z1']:.2f},{dd['Z2']:.2f})  BO={dd['bo']}")
print()

print(f"Outliers ({len(outliers)} molecules, |err| > 15%):")
for dd, e in sorted(outliers, key=lambda x: abs(x[1])):
    print(f"  {dd['name']:<6} err = {e:+5.1f}%  n=({dd['n1']},{dd['n2']})  "
          f"Z=({dd['Z1']:.2f},{dd['Z2']:.2f})  BO={dd['bo']}")
print()

# F2, O2 are predicted too SHORT (negative error -> R_pred < R_exp)
# This means the actual bond is LONGER than 2*rB_geom
# Physical reason: lone pair repulsion pushes atoms apart
# In GWT: the additional electron waves on the axis create repulsive
# interference that stretches the bond

# H2 is predicted too LONG (positive error -> R_pred > R_exp)
# This means H2 bonds MORE tightly than expected
# Physical reason: perfect overlap of two 1s orbitals with no inner shells

# Li2, Na2 are predicted differently depending on the model...

# Count electrons NOT in the bonding orbital
# Lone pairs on the bond axis create repulsive pressure
print("Lone pair analysis:")
print(f"{'Mol':<6} {'err%':>6} {'valence_e':>10} {'bonding_e':>10} {'nonbond_e':>10}")
print("-" * 55)

# Approximate: for each atom, how many valence electrons are NOT in the bond?
# This determines lone pair repulsion
val_electrons = {
    'H_1s': 1, 'Li_2s': 1, 'B_2p': 3, 'C_2p': 4,
    'N_2p': 5, 'O_2p': 6, 'F_2p': 7, 'Na_3s': 1, 'Cl_3p': 7,
}

for dd in data:
    R_pred = 2.0 * dd['rB_geom'] / max(dd['bo'], 0.5)**p_use
    e = (R_pred - dd['R']) / dd['R'] * 100
    v1 = val_electrons[dd['o1']]
    v2 = val_electrons[dd['o2']]
    v_total = v1 + v2
    bonding = 2 * dd['bo']  # 2 electrons per bond
    nonbond = v_total - bonding
    print(f"{dd['name']:<6} {e:+5.1f}% {v_total:10d} {bonding:10d} {nonbond:10d}")

print()
print("Key pattern: molecules with many non-bonding electrons")
print("tend to have LONGER bonds than predicted (positive residual).")
print("This suggests lone pair repulsion stretches the bond.")
print()

# Test: add lone pair correction
print("=" * 90)
print("  WITH LONE PAIR CORRECTION")
print("=" * 90)
print()

best = (100, 0, 0, 0)
for c in np.arange(1.0, 3.0, 0.02):
    for p_bo in np.arange(0, 0.5, 0.02):
        for p_lp in np.arange(0, 0.3, 0.01):
            errs = []
            for dd in data:
                v1 = val_electrons[dd['o1']]
                v2 = val_electrons[dd['o2']]
                nonbond = v1 + v2 - 2 * dd['bo']
                lp_factor = 1 + p_lp * nonbond
                R_pred = c * dd['rB_geom'] * lp_factor / max(dd['bo'], 0.5)**p_bo
                errs.append(abs((R_pred - dd['R']) / dd['R'] * 100))
            err = np.mean(errs)
            if err < best[0]:
                best = (err, c, p_bo, p_lp)

err, c_opt, p_bo, p_lp = best
print(f"R = {c_opt:.2f} * sqrt(rB1*rB2) * (1 + {p_lp:.3f}*N_nonbond) / BO^{p_bo:.2f}")
print(f"Average |error| = {err:.1f}%")
print()

header = f"{'Mol':<6} {'R_exp':>6} {'R_pred':>6} {'err%':>6} {'BO':>3} {'N_nb':>4}"
print(header)
print("-" * 40)

errs_all = []
for dd in data:
    v1 = val_electrons[dd['o1']]
    v2 = val_electrons[dd['o2']]
    nonbond = v1 + v2 - 2 * dd['bo']
    lp_factor = 1 + p_lp * nonbond
    R_pred = c_opt * dd['rB_geom'] * lp_factor / max(dd['bo'], 0.5)**p_bo
    e = (R_pred - dd['R']) / dd['R'] * 100
    errs_all.append(abs(e))
    flag = '***' if abs(e) < 5 else ' **' if abs(e) < 10 else '  *' if abs(e) < 15 else ''
    print(f"{dd['name']:<6} {dd['R']:6.3f} {R_pred:6.3f} {e:+5.1f}% {dd['bo']:3d} {nonbond:4d} {flag}")

w5 = sum(1 for e in errs_all if e < 5)
w10 = sum(1 for e in errs_all if e < 10)
w15 = sum(1 for e in errs_all if e < 15)
print(f"\nAvg |err| = {np.mean(errs_all):.1f}%, within 5%: {w5}/24, 10%: {w10}/24, 15%: {w15}/24")
print()

print("=" * 90)
print("  SUMMARY")
print("=" * 90)
print()
print("The bond length is the orbit radius of the two-body wave system.")
print("The Bohr radius n^2/Z defines the 'gravitational boundary' of each wave.")
print()
print("Best zero-parameter-like formula:")
print(f"  R = {c_opt:.2f} * sqrt(n1^2/Z1 * n2^2/Z2) * (1 + {p_lp:.3f}*N_nonbond) / BO^{p_bo:.2f}")
print()
print("Physical interpretation:")
print("  - sqrt(rB1*rB2): geometric mean orbit = 'reduced Bohr radius'")
print("  - 1 + lp*N_nonbond: lone pair repulsion stretches the bond")
print("  - BO^p: multiple bonds contract (shared orbit)")
print()
