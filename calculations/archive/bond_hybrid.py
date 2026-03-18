#!/usr/bin/env python3
"""
Hybrid Bond Calculator — Oh Lookup + 3D-Calibrated Ratios
==========================================================
Angular structure: Oh tensor product lookup (exact, O(1))
Energy scale: E_H from alpha_em (exact)
Channel ratios: measured on 3D cubic lattice (from Lagrangian)

Three channel types from 3D GPU measurement:
  sigma (pz+pz half):  w_sigma = 1.0    (reference, bonding)
  pi (px+px half):     w_pi ≈ 0.15      (weaker bonding, perpendicular)
  LP_perp (px+px full): w_LP ≈ -7.0     (REPULSIVE, from 3D geometry)

Bond energy:
  D_e = (pi/d^2) * E_harm * [n_sigma*w_sigma + n_pi*w_pi + n_LP*w_LP_eff]

Where w_LP_eff = w_LP * (angular_weight) accounts for the reduced
angular overlap of LP modes (they point sideways, not at the partner).

The pi/sigma ratio from 3D was 0.15, but Oh predicts cos(pi/d) = 0.5.
The 3D value may be noisy — let's test both and see which matches experiment.
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial

d = 3
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6  # 13.604 eV
C_bond = np.pi / d

# ============================================================
# 3D-CALIBRATED CHANNEL WEIGHTS
# ============================================================
# From 3D GPU simulation at R=2.0 Bohr:
#   sigma (pz+pz half): V = -0.337  → reference = 1.0
#   pi (px+px half):    V = -0.051  → ratio = 0.151
#   LP_perp (full+full): V = +2.305 → ratio = -6.83 (REPULSIVE)
#
# The LP ratio means: each facing LP pair repels with ~7x the
# strength of a sigma bond, but the angular overlap is reduced.
#
# For the bond formula, the effective LP weight per facing pair:
#   w_LP_eff = |LP_ratio| * angular_reduction
#
# The angular reduction comes from LP orbitals pointing perpendicular
# to the bond. On the Oh cube, perpendicular coupling goes through
# the pi channel: cos(pi/d) = 0.5. But the 3D sim already includes
# this geometry. So w_LP_eff = LP_ratio directly.
#
# However, LP_ratio = 6.83 seems high. The ratio scales with amplitude^2
# since V ~ A^2 for the potential energy. full = 2*half, so
# LP_raw / sigma_raw = (2A)^2 / A^2 = 4 in the linear regime.
# The measured 6.83 includes nonlinear enhancement.
#
# For the EFFECTIVE LP repulsion in the bond formula:
# LP repulsion per pair = (LP_ratio / amplitude_ratio^2) * sigma_coupling
#                       = (6.83 / 4) * sigma = 1.71 * sigma
#
# But we also need the 1/n^2 radial dilution for period 3+.

# Amplitude-corrected LP ratio
LP_RAW_RATIO = 6.83       # from 3D sim
AMP_RATIO_SQ = 4.0        # (A_full/A_half)^2
LP_INTRINSIC = LP_RAW_RATIO / AMP_RATIO_SQ  # ~1.71 per LP pair

# Pi ratio: 3D gave 0.15, Oh predicts 0.5. Test both.
PI_3D = 0.151
PI_OH = 0.5  # cos(pi/d)


# ============================================================
# ATOM DATABASE
# ============================================================
atoms = {
    'H':  {'Z':1,  'IE':13.598, 'n':1, 'p':0, 'lp':0, 'max_bo':1},
    'Li': {'Z':3,  'IE':5.392,  'n':2, 'p':0, 'lp':0, 'max_bo':1},
    'B':  {'Z':5,  'IE':8.298,  'n':2, 'p':1, 'lp':0, 'max_bo':3},
    'C':  {'Z':6,  'IE':11.260, 'n':2, 'p':2, 'lp':0, 'max_bo':4},
    'N':  {'Z':7,  'IE':14.534, 'n':2, 'p':3, 'lp':1, 'max_bo':3},
    'O':  {'Z':8,  'IE':13.618, 'n':2, 'p':4, 'lp':2, 'max_bo':2},
    'F':  {'Z':9,  'IE':17.423, 'n':2, 'p':5, 'lp':3, 'max_bo':1},
    'Na': {'Z':11, 'IE':5.139,  'n':3, 'p':0, 'lp':0, 'max_bo':1},
    'Si': {'Z':14, 'IE':8.152,  'n':3, 'p':2, 'lp':0, 'max_bo':4},
    'P':  {'Z':15, 'IE':10.487, 'n':3, 'p':3, 'lp':1, 'max_bo':3},
    'S':  {'Z':16, 'IE':10.360, 'n':3, 'p':4, 'lp':2, 'max_bo':2},
    'Cl': {'Z':17, 'IE':12.968, 'n':3, 'p':5, 'lp':3, 'max_bo':1},
    'K':  {'Z':19, 'IE':4.341,  'n':4, 'p':0, 'lp':0, 'max_bo':1},
}


def E_harmonic(a, b):
    Ea, Eb = atoms[a]['IE'], atoms[b]['IE']
    return 2*Ea*Eb/(Ea+Eb)


def bond_energy(sym_a, sym_b, bo, w_pi_choice='oh'):
    """
    Bond energy from Oh lookup + 3D-calibrated weights.

    D_e = (pi/d^2) * E_harm * [sigma_term + pi_term - LP_term]

    sigma_term = 1.0 (always 1 sigma bond)
    pi_term = (bo-1) * w_pi (each pi bond)
    LP_term = n_LP_facing * LP_INTRINSIC * (2/n_max)^2

    The (2/n_max)^2 factor = radial dilution for period 3+:
      period 2 (n=2): (2/2)^2 = 1.0 (full strength)
      period 3 (n=3): (2/3)^2 = 0.44 (more diffuse LP)
    """
    a = atoms[sym_a]
    b = atoms[sym_b]
    E_harm = E_harmonic(sym_a, sym_b)
    n_max = max(a['n'], b['n'])

    w_pi = PI_OH if w_pi_choice == 'oh' else PI_3D

    # Sigma: always 1 channel
    sigma = 1.0

    # Pi: (bo-1) channels
    pi_term = (bo - 1) * w_pi

    # LP: facing lone pairs (perpendicular full+full)
    # Count facing LP: min of each atom's LP count
    # For LP counting: LP = p - d for p > d, else 0
    lp_a = max(0, a['p'] - d) if a['p'] > 0 else 0
    lp_b = max(0, b['p'] - d) if b['p'] > 0 else 0
    n_lp = min(lp_a, lp_b)

    # Radial dilution: LP overlap falls as (2/n)^2
    radial = (2.0 / n_max)**2 if n_max >= 2 else 1.0

    lp_term = n_lp * LP_INTRINSIC * radial

    # Ionic contribution (charge transfer)
    delta_E = abs(a['IE'] - b['IE'])
    E_avg = (a['IE'] + b['IE']) / 2
    asym = delta_E / E_avg
    c_ionic = 1/(2*d+1)  # 1/7
    D_ionic = c_ionic * delta_E if asym > 0.1 else 0

    # Total covalent
    coupling = sigma + pi_term - lp_term
    D_cov = (np.pi / d**2) * E_harm * max(coupling, 0)

    D_total = D_cov + D_ionic

    return {
        'D_total': D_total, 'D_cov': D_cov, 'D_ionic': D_ionic,
        'sigma': sigma, 'pi': pi_term, 'lp': lp_term, 'coupling': coupling,
        'E_harm': E_harm, 'n_lp': n_lp,
    }


# ============================================================
# EXPERIMENTAL DATA
# ============================================================
exp_bonds = [
    ('H', 'H', 1, 4.478, 'H2'),
    ('Li', 'Li', 1, 1.046, 'Li2'),
    ('N', 'N', 3, 9.759, 'N2'),
    ('O', 'O', 2, 5.116, 'O2'),
    ('F', 'F', 1, 1.602, 'F2'),
    ('H', 'F', 1, 5.869, 'HF'),
    ('H', 'Cl', 1, 4.434, 'HCl'),
    ('Na', 'Cl', 1, 4.230, 'NaCl'),
    ('Li', 'H', 1, 2.429, 'LiH'),
    ('H', 'O', 1, 4.392, 'OH'),
    ('C', 'O', 3, 11.09, 'CO'),
    ('N', 'O', 2, 6.497, 'NO'),
    ('H', 'N', 1, 3.910, 'NH'),
    ('C', 'H', 1, 4.290, 'CH'),
    ('C', 'C', 1, 3.600, 'C-C'),
    ('C', 'N', 3, 7.760, 'CN'),
    ('C', 'C', 2, 6.360, 'C=C'),
    ('C', 'O', 2, 7.710, 'C=O'),
    ('C', 'C', 3, 8.700, 'C≡C'),
    ('N', 'H', 1, 4.513, 'NH(NH3)'),
    ('O', 'H', 1, 4.790, 'OH(H2O)'),
    ('Cl', 'Cl', 1, 2.514, 'Cl2'),
    ('S', 'H', 1, 3.78, 'SH'),
    ('S', 'S', 2, 4.37, 'S2'),
    ('P', 'H', 1, 3.44, 'PH'),
]


# ============================================================
# RUN WITH Oh CHANNEL WEIGHTS
# ============================================================
for w_label, w_choice in [('Oh (w_pi=0.5)', 'oh'), ('3D (w_pi=0.15)', '3d')]:
    print("=" * 90)
    print(f"  HYBRID BOND CALCULATOR: Oh Lookup + 3D LP Ratio")
    print(f"  w_pi = {w_label}, LP_intrinsic = {LP_INTRINSIC:.3f}")
    print(f"  D = (pi/d^2) * E_harm * [1 + (bo-1)*w_pi - n_LP*{LP_INTRINSIC:.2f}*(2/n)^2]")
    print("=" * 90)
    print()
    print(f"  {'Name':>10} {'bo':>3} {'E_h':>6} {'sig':>5} {'pi':>5} {'lp':>5} "
          f"{'cpl':>6} {'D_cov':>7} {'D_ion':>5} {'D_pred':>7} {'D_obs':>7} {'err%':>7}")
    print("  " + "-" * 90)

    errs = []
    for sym_a, sym_b, bo, De_obs, name in exp_bonds:
        r = bond_energy(sym_a, sym_b, bo, w_choice)
        err = (r['D_total'] - De_obs) / De_obs * 100
        errs.append(abs(err))
        print(f"  {name:>10} {bo:>3} {r['E_harm']:>6.1f} {r['sigma']:>5.2f} {r['pi']:>5.2f} "
              f"{r['lp']:>5.2f} {r['coupling']:>6.3f} {r['D_cov']:>7.3f} {r['D_ionic']:>5.2f} "
              f"{r['D_total']:>7.3f} {De_obs:>7.3f} {err:>+7.1f}%")

    cov = [e for i, e in enumerate(errs) if exp_bonds[i][4] not in ('Li2', 'NaCl', 'LiH')]
    print()
    print(f"  ALL:      mean={np.mean(errs):.1f}%, median={np.median(errs):.1f}%")
    print(f"  Covalent: mean={np.mean(cov):.1f}%, median={np.median(cov):.1f}%")
    print(f"  Within 10%: {sum(1 for e in errs if e<10)}/{len(errs)}")
    print(f"  Within 20%: {sum(1 for e in errs if e<20)}/{len(errs)}")
    print(f"  Within 30%: {sum(1 for e in errs if e<30)}/{len(errs)}")
    print()
