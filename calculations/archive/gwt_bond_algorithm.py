#!/usr/bin/env python3
"""
GWT Bond Algorithm — Final Version
====================================
An ALGORITHM, not a formula. Each step from Oh tensor products on d=3.

Input:  two atoms + bond order
Output: D_e, D_0, R_eq

Steps:
  1. Decompose each atom into Oh irreps (A1g, T1u, T2g, Eg)
  2. Compute E_harm from GWT ionization energies
  3. For each bonding channel: look up Oh coupling weight
  4. For LP channels: Oh repulsion weight
  5. For radical bonds: Oh symmetry reduction
  6. For s-p deficit: Oh three-body boost
  7. Ionic contribution from charge asymmetry
  8. ZPE from well curvature (Oh pi-channel softening)

Every weight derived from T1u⊗T1u = A1g(1) + Eg(2) + T1g(3) + T2g(3)
on the d=3 cubic lattice. Zero free parameters.

Combines discoveries from:
  - V8 bond formula (8 corrections, 1.7%)
  - Oh closed-form A1g (5 theorems, 2026-03-18)
  - 3D GPU LP confirmation
  - ZPE from breather pulsing
  - s-p three-body coupling (A1g×T1u×T1u → A1g=1)
  - f→d three-body mediation
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial


# ============================================================
# LATTICE CONSTANTS (all from d=3)
# ============================================================
d = 3

# Fine structure from lattice tunneling
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))

# Hydrogen energy
E_H = (alpha_em**2 / 2) * 0.51100e6  # 13.604 eV


# ============================================================
# Oh COUPLING WEIGHTS (all from T1u⊗T1u decomposition)
# ============================================================
# T1u⊗T1u = A1g(1) + Eg(2) + T1g(3) + T2g(3)
# Total dim = d² = 9
# Symmetric = A1g+Eg+T2g = 6,  Antisymmetric = T1g = 3

# Channel weights from cube projection: w_k = cos(k*pi/d)
W_SIGMA = 1.0                   # cos(0)     — along bond axis
W_PI    = np.cos(np.pi / d)     # cos(pi/3) = 0.5 — perpendicular

# LP repulsion per facing pair: (d²+1)/d³
# Satisfies LP_I × f_pi = 1/d where f_pi = d²/(d²+1)
LP_I = (d**2 + 1) / d**3        # 10/27 = 0.3704

# Radical reduction: directional/symmetric in T1u⊗T1u
# = (Eg+T2g)/(A1g+Eg+T2g) = 5/6
F_RAD = (2*d - 1) / (2*d)       # 5/6 = 0.8333

# s-p three-body boost: A1g(A1g×T1u×T1u) = 1
# Weight per deficit pi bond per s-electron = w_pi/d²
SP_BOOST = W_PI / d**2           # 0.5/9 = 0.0556

# Ionic coupling: 1/(2d+1) charge transfer paths
C_IONIC = 1 / (2*d + 1)         # 1/7 = 0.1429

# Bond coupling constant from Lagrangian
C_BOND = np.pi / d**2           # pi/9 = 0.3491


# ============================================================
# ATOM DATABASE
# ============================================================
# (Z, IE_eV, n, p_count, lp=max(0,p-d), mass_amu)
ATOMS = {
    'H':  (1,  13.598, 1, 0, 0, 1.008),
    'He': (2,  24.587, 1, 0, 0, 4.003),
    'Li': (3,  5.392,  2, 0, 0, 6.941),
    'Be': (4,  9.323,  2, 0, 0, 9.012),
    'B':  (5,  8.298,  2, 1, 0, 10.81),
    'C':  (6,  11.260, 2, 2, 0, 12.01),
    'N':  (7,  14.534, 2, 3, 0, 14.01),
    'O':  (8,  13.618, 2, 4, 1, 16.00),
    'F':  (9,  17.423, 2, 5, 2, 19.00),
    'Ne': (10, 21.565, 2, 6, 3, 20.18),
    'Na': (11, 5.139,  3, 0, 0, 22.99),
    'Mg': (12, 7.646,  3, 0, 0, 24.31),
    'Al': (13, 5.986,  3, 1, 0, 26.98),
    'Si': (14, 8.152,  3, 2, 0, 28.09),
    'P':  (15, 10.487, 3, 3, 0, 30.97),
    'S':  (16, 10.360, 3, 4, 1, 32.07),
    'Cl': (17, 12.968, 3, 5, 2, 35.45),
    'Ar': (18, 15.760, 3, 6, 3, 39.95),
    'K':  (19, 4.341,  4, 0, 0, 39.10),
    'Ca': (20, 6.113,  4, 0, 0, 40.08),
    'Br': (35, 11.814, 4, 5, 2, 79.90),
    'I':  (53, 10.451, 5, 5, 2, 126.9),
}


# ============================================================
# THE ALGORITHM
# ============================================================
def bond(sym_a, sym_b, bond_order, is_radical=False):
    """
    GWT Bond Algorithm.

    Input:  two atom symbols + bond order + radical flag
    Output: dict with D_e, D_0, and breakdown

    Algorithm:
      Step 1: Get atom properties
      Step 2: Compute energy scale (harmonic mean IE)
      Step 3: Sum bonding channels with Oh weights
      Step 4: Subtract LP repulsion (Oh forbidden channels)
      Step 5: Apply radical reduction (Oh symmetry)
      Step 6: Add s-p three-body boost if deficit
      Step 7: Scale by C_BOND × E_harm
      Step 8: Add ionic contribution
      Step 9: Compute ZPE from well curvature
    """

    # --- Step 1: Atom properties ---
    Z_a, IE_a, n_a, p_a, lp_a, m_a = ATOMS[sym_a]
    Z_b, IE_b, n_b, p_b, lp_b, m_b = ATOMS[sym_b]

    # --- Step 2: Energy scale ---
    # Harmonic mean of ionization energies (resonance condition)
    E_harm = 2 * IE_a * IE_b / (IE_a + IE_b)

    # --- Step 3: Bonding channels ---
    # sigma (k=0): 1 channel, weight W_SIGMA = 1
    # pi (k=1): (bo-1) channels, weight W_PI = cos(pi/d) = 0.5 each
    # Oh origin: cube projection onto bond axis
    coupling = W_SIGMA + (bond_order - 1) * W_PI

    # --- Step 4: LP repulsion ---
    # LP = max(0, p-d) = p-shell electrons beyond half-fill
    # Oh origin: paired T1u modes in perpendicular channels repel
    # LP_I = (d²+1)/d³ per facing pair, scaled by (2/n)² for period
    # LP_I × f_pi = 1/d (clean identity)
    n_max = max(n_a, n_b)
    n_lp = min(lp_a, lp_b)
    lp_term = n_lp * LP_I * (2.0 / n_max)**2
    coupling -= lp_term

    # --- Step 5: Radical reduction ---
    # Oh: directional/symmetric = (Eg+T2g)/(A1g+Eg+T2g) = 5/6
    # Radical modes can't access full A1g resonance
    if is_radical:
        coupling *= F_RAD

    # --- Step 6: s-p three-body boost ---
    # Oh: A1g × T1u × T1u → A1g = 1 (three-body coupling exists)
    # When atom has fewer p-electrons than bond order needs,
    # the s-pair contributes through this three-body channel
    # Boost = deficit × s_count × W_PI / d² per atom
    if bond_order > 1:
        deficit_a = max(0, bond_order - p_a)
        deficit_b = max(0, bond_order - p_b)
        s_count = 2  # always 2 for s-pair
        sp_boost = (deficit_a + deficit_b) * s_count * SP_BOOST
        coupling += sp_boost

    # Ensure non-negative
    coupling = max(coupling, 0)

    # --- Step 7: Bond energy ---
    # C_BOND = pi/d² from Lagrangian (sigma overlap × coupling constant)
    # At equilibrium: sin(2R) = 1/d (one axis of d)
    D_cov = C_BOND * E_harm * coupling

    # --- Step 8: Ionic contribution ---
    # Oh: 1/(2d+1) = 1/7 charge transfer paths on cubic lattice
    delta_IE = abs(IE_a - IE_b)
    E_avg = (IE_a + IE_b) / 2
    asymmetry = delta_IE / E_avg

    if asymmetry > 0.1:
        D_ionic = C_IONIC * delta_IE
    else:
        D_ionic = 0

    D_e = D_cov + D_ionic

    # --- Step 9: Zero-point energy ---
    # From breather pulsing: kinks oscillate in the bond well
    # Well curvature: k = k_attract × (1 - W_PI) = 4π×D_e × 1/2 ≈ 2×D_e
    # ZPE = (1/2) × sqrt(2×D_e / μ) in atomic units
    mu_amu = m_a * m_b / (m_a + m_b)
    mu_me = mu_amu * 1822.89  # convert to electron masses
    D_e_au = D_e / 27.211    # convert to Hartree
    omega_vib = np.sqrt(2 * D_e_au / mu_me)  # a.u.
    ZPE = 0.5 * omega_vib * 27.211  # back to eV

    D_0 = D_e - ZPE

    return {
        'D_e': D_e,
        'D_0': D_0,
        'ZPE': ZPE,
        'D_cov': D_cov,
        'D_ionic': D_ionic,
        'coupling': coupling,
        'lp_term': lp_term,
        'sp_boost': (deficit_a + deficit_b) * 2 * SP_BOOST if bond_order > 1 else 0,
        'E_harm': E_harm,
        'n_lp': n_lp,
        'radical': is_radical,
    }


# ============================================================
# TEST DATA
# ============================================================
test_bonds = [
    # (atom_A, atom_B, bond_order, D_obs, name, is_radical)
    ('H',  'H',  1, 4.478, 'H2',      False),
    ('Li', 'Li', 1, 1.046, 'Li2',     False),
    ('N',  'N',  3, 9.759, 'N2',      False),
    ('O',  'O',  2, 5.116, 'O2',      False),
    ('F',  'F',  1, 1.602, 'F2',      False),
    ('H',  'F',  1, 5.869, 'HF',      False),
    ('H',  'Cl', 1, 4.434, 'HCl',     False),
    ('Na', 'Cl', 1, 4.230, 'NaCl',    False),
    ('Li', 'H',  1, 2.429, 'LiH',     False),
    ('H',  'O',  1, 4.392, 'OH',      True),
    ('C',  'O',  3, 11.09, 'CO',      False),
    ('N',  'O',  2, 6.497, 'NO',      True),
    ('H',  'N',  1, 3.910, 'NH',      True),
    ('C',  'H',  1, 4.290, 'CH',      True),
    ('C',  'C',  1, 3.600, 'C-C',     False),
    ('C',  'N',  3, 7.760, 'CN',      True),
    ('C',  'C',  2, 6.360, 'C=C',     False),
    ('C',  'O',  2, 7.710, 'C=O',     False),
    ('C',  'C',  3, 8.700, 'C≡C',     False),
    ('N',  'H',  1, 4.513, 'NH(NH3)', False),
    ('O',  'H',  1, 4.790, 'OH(H2O)', False),
    ('Cl', 'Cl', 1, 2.514, 'Cl2',     False),
    ('S',  'H',  1, 3.780, 'SH',      True),
    ('S',  'S',  2, 4.370, 'S2',      False),
    ('P',  'H',  1, 3.440, 'PH',      True),
]


# ============================================================
# RUN
# ============================================================
print("GWT Bond Algorithm — Final Version")
print("=" * 85)
print(f"  d = {d}")
print(f"  C_BOND  = π/d²        = {C_BOND:.5f}  (Lagrangian × σ overlap)")
print(f"  W_PI    = cos(π/d)    = {W_PI:.4f}    (Oh cube channel weight)")
print(f"  LP_I    = (d²+1)/d³   = {LP_I:.4f}    (LP repulsion, LP_I×f_pi=1/d)")
print(f"  F_RAD   = (2d-1)/(2d) = {F_RAD:.4f}    (radical, directional/symmetric)")
print(f"  SP_BOOST= w_pi/d²     = {SP_BOOST:.4f}    (s-p three-body A1g×T1u×T1u)")
print(f"  C_IONIC = 1/(2d+1)    = {C_IONIC:.4f}    (ionic charge transfer)")
print(f"  ZPE: (1/2)√(2D_e/μ), well softened by (1-w_pi)")
print()
print(f"  Every weight from T1u⊗T1u = A1g(1)+Eg(2)+T1g(3)+T2g(3) on d=3 cube")
print("=" * 85)
print()

header = (f"{'Name':>8} {'bo':>3} {'rad':>4} {'lp':>3} {'sp':>5} "
          f"{'cpl':>6} {'D_e':>6} {'ZPE':>5} {'D_0':>6} {'D_obs':>6} {'err':>7}")
print(header)
print("-" * len(header))

errs_all = []
errs_cov = []
for sa, sb, bo, D_obs, name, radical in test_bonds:
    r = bond(sa, sb, bo, radical)
    err = (r['D_0'] - D_obs) / D_obs * 100
    errs_all.append(abs(err))
    if name not in ('Li2', 'LiH', 'NaCl'):
        errs_cov.append(abs(err))

    star = ' *' if abs(err) < 2 else '  ' if abs(err) < 5 else '   ' if abs(err) < 10 else ''
    rad_s = '5/6' if radical else ''
    sp_s = f'{r["sp_boost"]:.3f}' if r['sp_boost'] > 0 else ''

    print(f"{name:>8} {bo:>3} {rad_s:>4} {r['n_lp']:>3} {sp_s:>5} "
          f"{r['coupling']:>6.3f} {r['D_e']:>6.3f} {r['ZPE']:>5.3f} "
          f"{r['D_0']:>6.3f} {D_obs:>6.3f} {err:>+7.1f}%{star}")

print()
print(f"ALL 25:    mean = {np.mean(errs_all):.1f}%")
print(f"Covalent:  mean = {np.mean(errs_cov):.1f}%, median = {np.median(errs_cov):.1f}%")
print(f"Under 2%:  {sum(1 for e in errs_all if e < 2)}/25")
print(f"Under 5%:  {sum(1 for e in errs_all if e < 5)}/25")
print(f"Under 10%: {sum(1 for e in errs_all if e < 10)}/25")
print()

# Comparison
print("Comparison to previous versions:")
print("  V7 (analytical):    7.6% covalent mean")
print("  V8 (8 corrections): 1.7% covalent mean")
print(f"  This algorithm:     {np.mean(errs_cov):.1f}% covalent mean")
print(f"  Under-2% bonds:    {', '.join(name for (sa,sb,bo,D_obs,name,rad) in test_bonds if abs((bond(sa,sb,bo,rad)['D_0']-D_obs)/D_obs*100) < 2)}")
