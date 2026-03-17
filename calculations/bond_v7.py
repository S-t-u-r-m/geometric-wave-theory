#!/usr/bin/env python3
"""
Bond energy v7: n-dependent LP repulsion + ionic/dative correction.

v6 result: 9.3% covalent, zero params. Three remaining systematics:
  1. Cl₂ (-59%): LP repulsion overestimated for n=3 (diffuse LPs)
  2. CO (-16%), HF (-14%), C=O (-13%): polar bonds underpredicted
  3. CN (+23%), NH (+18%): overpredicted

Fix 1: LP repulsion scales with 1/(d*n²) instead of 1/(d*(d+1))
  - Period 2 (n=2): 1/(3*4) = 1/12 = 0.0833 (unchanged! since (d-1)²=d+1=4 for d=3)
  - Period 3 (n=3): 1/(3*9) = 1/27 = 1/d³ ≈ 0.0370

  Physical: LP cloud density at bond midpoint ~ 1/n per atom, overlap ~ 1/n²

Fix 2: Non-facing LP contribute POSITIVE energy (partial dative bonding)
  When atom B has LP but atom A has none, B can donate into A's empty orbital.
  This is a PARTIAL pi bond: weight = w_pi_bond / d = (d-1)/d² per dative LP.
  But only the LP on the MORE electronegative atom (higher E_ion) can donate.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np

d = 3
w_pi_bond = (d - 1) / d  # 2/3

E_ion_pred = {
    'H': 13.604, 'He': 25.192, 'Li': 5.616, 'Be': 8.920,
    'B': 8.345, 'C': 10.900, 'N': 13.981, 'O': 13.389,
    'F': 16.970, 'Ne': 22.000, 'Na': 5.179, 'Si': 7.971,
    'P': 10.178, 'S': 10.287, 'Cl': 12.391, 'K': 4.333,
}

val_info = {  # (n, lone_pairs, max_bonds)
    'H': (1, 0, 1), 'Li': (2, 0, 1), 'Be': (2, 0, 2),
    'B': (2, 0, 3), 'C': (2, 0, 4), 'N': (2, 1, 3),
    'O': (2, 2, 2), 'F': (2, 3, 1), 'Na': (3, 0, 1),
    'Si': (3, 0, 4), 'P': (3, 1, 3), 'S': (3, 2, 2),
    'Cl': (3, 3, 1), 'K': (4, 0, 1),
}

exp_bonds = [
    ('H', 'H', 4.478, 0.741, 1, 'H₂'),
    ('Li', 'Li', 1.046, 2.673, 1, 'Li₂'),
    ('N', 'N', 9.759, 1.098, 3, 'N₂'),
    ('O', 'O', 5.116, 1.208, 2, 'O₂'),
    ('F', 'F', 1.602, 1.412, 1, 'F₂'),
    ('H', 'F', 5.869, 0.917, 1, 'HF'),
    ('H', 'Cl', 4.434, 1.275, 1, 'HCl'),
    ('Na', 'Cl', 4.230, 2.361, 1, 'NaCl'),
    ('Li', 'H', 2.429, 1.596, 1, 'LiH'),
    ('H', 'O', 4.392, 0.970, 1, 'OH'),
    ('C', 'O', 11.09, 1.128, 3, 'CO'),
    ('N', 'O', 6.497, 1.151, 2, 'NO'),
    ('H', 'N', 3.910, 1.036, 1, 'NH'),
    ('C', 'H', 4.290, 1.089, 1, 'CH'),
    ('C', 'C', 3.600, 1.540, 1, 'C-C'),
    ('C', 'N', 7.760, 1.170, 3, 'CN'),
    ('C', 'C', 6.360, 1.340, 2, 'C=C'),
    ('C', 'O', 7.710, 1.200, 2, 'C=O'),
    ('C', 'C', 8.700, 1.200, 3, 'C≡C'),
    ('N', 'H', 4.513, 1.012, 1, 'NH(NH₃)'),
    ('O', 'H', 4.790, 0.958, 1, 'OH(H₂O)'),
    ('Cl', 'Cl', 2.514, 1.988, 1, 'Cl₂'),
    ('S', 'H', 3.78, 1.34, 1, 'SH'),
    ('S', 'S', 4.37, 1.89, 2, 'S₂'),
    ('P', 'H', 3.44, 1.42, 1, 'PH'),
]

seen = set()
bonds = []
for b in exp_bonds:
    a1, a2, De, re, bo, name = b
    key = (min(a1,a2), max(a1,a2), bo, name)
    if key not in seen:
        seen.add(key)
        bonds.append(b)

def E_harm(a1, a2):
    E1, E2 = E_ion_pred[a1], E_ion_pred[a2]
    return 2 * E1 * E2 / (E1 + E2)

def lp_facing(a1, a2):
    return min(val_info[a1][1], val_info[a2][1])

def eff_bo(bo):
    return 1 + (bo - 1) * w_pi_bond

def n_avg(a1, a2):
    """Average n for LP repulsion scaling."""
    return (val_info[a1][0] + val_info[a2][0]) / 2

# =====================================================================
# FIX 1: n-dependent LP repulsion
# C_lp(n) = 1/(d * n²)
# =====================================================================
print("=" * 100)
print("  FIX 1: LP repulsion ~ 1/(d*n²)")
print("  Period 2: 1/(3*4) = 1/12 = 0.0833 (same as v6)")
print("  Period 3: 1/(3*9) = 1/27 = 1/d³ = 0.0370")
print("=" * 100)

def De_v7a(a1, a2, bo):
    """n-dependent LP repulsion."""
    Eh = E_harm(a1, a2)
    lp = lp_facing(a1, a2)
    eb = eff_bo(bo)
    n = max(val_info[a1][0], val_info[a2][0])  # use larger n (more diffuse)
    C_lp = 1 / (d * n**2)
    return max(0, (eb / d - lp * C_lp) * Eh)

print(f"\n  {'Name':>10} {'bo':>3} {'eff':>5} {'lp':>3} {'n':>2} {'C_lp':>7} "
      f"{'De_pred':>8} {'De_obs':>8} {'err%':>7}")
errs_a = []
for a1, a2, De, re, bo, name in bonds:
    n = max(val_info[a1][0], val_info[a2][0])
    pred = De_v7a(a1, a2, bo)
    err = (pred - De) / De * 100
    errs_a.append(abs(err))
    eb = eff_bo(bo)
    lp = lp_facing(a1, a2)
    clp = 1/(d*n**2)
    print(f"  {name:>10} {bo:3.0f} {eb:5.3f} {lp:3d} {n:2d} {clp:7.4f} "
          f"{pred:8.3f} {De:8.3f} {err:+7.1f}%")

covalent = [i for i, (a1,a2,*_) in enumerate(bonds)
            if not (a1=='Li' and a2=='Li') and not ({a1,a2} == {'Na','Cl'})]
cov_errs_a = [errs_a[i] for i in covalent]
print(f"\n  All: {np.mean(errs_a):.1f}%")
print(f"  Covalent (excl Li₂, NaCl): {np.mean(cov_errs_a):.1f}%")

# =====================================================================
# FIX 2: Dative LP donation from electronegative atom
# When one atom has LP and the partner has empty capacity,
# the LP partially donates into the bond.
# =====================================================================
print("\n" + "=" * 100)
print("  FIX 2: Dative LP donation")
print("  Non-facing LP on more electronegative atom donate partially")
print("  Each dative LP adds C_dat * E_harm to the bond")
print("=" * 100)

def dative_lp(a1, a2, bo):
    """Count dative LP: LP on one atom that can donate to empty orbitals on the other."""
    lp1, lp2 = val_info[a1][1], val_info[a2][1]
    max_b1, max_b2 = val_info[a1][2], val_info[a2][2]
    E1, E2 = E_ion_pred[a1], E_ion_pred[a2]

    # Available empty orbitals: max_bonds - bonds_used
    # For atom making 'bo' bonds, remaining capacity = max_bonds - bo
    empty1 = max(0, max_b1 - bo)
    empty2 = max(0, max_b2 - bo)

    # LP from atom 2 can donate into empty orbitals on atom 1
    # But only if atom 2 is more electronegative (higher E_ion)
    dat = 0
    if E2 >= E1 and lp2 > 0 and empty1 > 0:
        dat += min(lp2, empty1)
    if E1 > E2 and lp1 > 0 and empty2 > 0:
        dat += min(lp1, empty2)
    return dat

print(f"\n  Dative LP analysis:")
print(f"  {'Name':>10} {'LP_A':>4} {'LP_B':>4} {'emptyA':>6} {'emptyB':>6} {'dative':>6}")
for a1, a2, De, re, bo, name in bonds:
    lp1, lp2 = val_info[a1][1], val_info[a2][1]
    empty1 = max(0, val_info[a1][2] - bo)
    empty2 = max(0, val_info[a2][2] - bo)
    dat = dative_lp(a1, a2, bo)
    print(f"  {name:>10} {lp1:4d} {lp2:4d} {empty1:6d} {empty2:6d} {dat:6d}")

# Sweep C_dative
print(f"\n  Sweep dative coefficient:")
for C_dat_trial in [0, 1/(d**2*d), 1/(d**2*(d+1)), w_pi_bond/d**2, 1/(d**2), w_pi_bond/d,
                     1/(d*(d+1)), 1/d**2 * w_pi_bond]:
    errs_d = []
    for a1, a2, De, re, bo, name in bonds:
        if (a1=='Li' and a2=='Li') or ({a1,a2} == {'Na','Cl'}):
            continue
        n = max(val_info[a1][0], val_info[a2][0])
        Eh = E_harm(a1, a2)
        lp = lp_facing(a1, a2)
        eb = eff_bo(bo)
        C_lp = 1 / (d * n**2)
        dat = dative_lp(a1, a2, bo)
        pred = max(0, (eb / d - lp * C_lp + dat * C_dat_trial) * Eh)
        errs_d.append(abs((pred - De) / De * 100))
    print(f"  C_dat={C_dat_trial:.5f}  cov_err={np.mean(errs_d):.2f}%")

# Fine sweep
best = None
for C_dat in np.arange(0.0, 0.2, 0.0005):
    sq = 0
    for a1, a2, De, re, bo, name in bonds:
        n = max(val_info[a1][0], val_info[a2][0])
        Eh = E_harm(a1, a2)
        lp = lp_facing(a1, a2)
        eb = eff_bo(bo)
        C_lp = 1 / (d * n**2)
        dat = dative_lp(a1, a2, bo)
        pred = max(0, (eb / d - lp * C_lp + dat * C_dat) * Eh)
        sq += ((pred - De) / De)**2
    if best is None or sq < best[1]:
        best = (C_dat, sq)

C_dat_best = best[0]
print(f"\n  Best fit C_dat = {C_dat_best:.4f}")
print(f"  d=3 candidates: 1/d²={1/d**2:.4f}, w_pi/d²={w_pi_bond/d**2:.4f}, "
      f"1/(d(d+1))={1/(d*(d+1)):.4f}")
print(f"                  (d-1)/d³={w_pi_bond/d**2:.4f}")

# =====================================================================
# COMBINED MODEL: n-LP + dative
# =====================================================================
print("\n" + "=" * 100)
print("  COMBINED MODEL v7: n-LP + dative")
print(f"  De = [(eff_bo + dat*C_dat)*d⁻¹ - lp/(d*n²)] × E_harm")
print(f"  C_dat_fit = {C_dat_best:.4f}")
print("=" * 100)

print(f"\n  {'Name':>10} {'bo':>3} {'lp':>3} {'dat':>3} {'De_cov':>7} {'De_dat':>7} "
      f"{'De_lp':>7} {'De_pred':>8} {'De_obs':>8} {'err%':>7}")
errs_comb = []
for a1, a2, De, re, bo, name in bonds:
    n = max(val_info[a1][0], val_info[a2][0])
    Eh = E_harm(a1, a2)
    lp = lp_facing(a1, a2)
    eb = eff_bo(bo)
    C_lp = 1 / (d * n**2)
    dat = dative_lp(a1, a2, bo)
    De_cov = eb / d * Eh
    De_dat = dat * C_dat_best / d * Eh  # dative adds to eff_bo
    De_lp = lp * C_lp * Eh
    pred = max(0, De_cov + De_dat - De_lp)
    err = (pred - De) / De * 100
    errs_comb.append(abs(err))
    print(f"  {name:>10} {bo:3.0f} {lp:3d} {dat:3d} {De_cov:7.2f} {De_dat:7.3f} "
          f"{De_lp:7.3f} {pred:8.3f} {De:8.3f} {err:+7.1f}%")

cov_comb = [errs_comb[i] for i in covalent]
print(f"\n  All: {np.mean(errs_comb):.1f}%")
print(f"  Covalent (excl Li₂, NaCl): {np.mean(cov_comb):.1f}%")

# =====================================================================
# Test with exact d=3 coefficient for dative
# =====================================================================
print("\n" + "=" * 100)
print("  TEST: d=3 coefficients for dative term")
print("=" * 100)

for C_dat_name, C_dat_val in [
    ('best_fit', C_dat_best),
    ('1/d²', 1/d**2),
    ('w_pi/d²', w_pi_bond/d**2),
    ('1/(d(d+1))', 1/(d*(d+1))),
    ('1/(2d)', 1/(2*d)),
    ('w_pi/d', w_pi_bond/d),
]:
    errs_test = []
    for a1, a2, De, re, bo, name in bonds:
        if (a1=='Li' and a2=='Li') or ({a1,a2} == {'Na','Cl'}):
            continue
        n = max(val_info[a1][0], val_info[a2][0])
        Eh = E_harm(a1, a2)
        lp = lp_facing(a1, a2)
        eb = eff_bo(bo)
        C_lp = 1 / (d * n**2)
        dat = dative_lp(a1, a2, bo)
        pred = max(0, (eb / d + dat * C_dat_val / d - lp * C_lp) * Eh)
        errs_test.append(abs((pred - De) / De * 100))
    print(f"  C_dat = {C_dat_name:>12} = {C_dat_val:.5f}  →  cov mean err = {np.mean(errs_test):.2f}%")

# =====================================================================
# Summary of what each fix contributes
# =====================================================================
print("\n" + "=" * 100)
print("  IMPROVEMENT BREAKDOWN")
print("=" * 100)

# v6 baseline (constant LP)
errs_v6 = []
for a1, a2, De, re, bo, name in bonds:
    if (a1=='Li' and a2=='Li') or ({a1,a2} == {'Na','Cl'}):
        continue
    Eh = E_harm(a1, a2)
    lp = lp_facing(a1, a2)
    eb = eff_bo(bo)
    pred = max(0, (eb / d - lp * 1/(d*(d+1))) * Eh)
    errs_v6.append(abs((pred - De) / De * 100))

# v7a: n-dependent LP only
errs_v7a_cov = []
for a1, a2, De, re, bo, name in bonds:
    if (a1=='Li' and a2=='Li') or ({a1,a2} == {'Na','Cl'}):
        continue
    pred = De_v7a(a1, a2, bo)
    errs_v7a_cov.append(abs((pred - De) / De * 100))

print(f"  v6 (constant LP):     {np.mean(errs_v6):.1f}%")
print(f"  + n-dependent LP:     {np.mean(errs_v7a_cov):.1f}%")
print(f"  + dative (best fit):  {np.mean(cov_comb):.1f}%  (excl Li₂, NaCl)")
