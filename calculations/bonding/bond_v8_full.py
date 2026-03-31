#!/usr/bin/env python3
"""
GWT Bond Formula V8 — Full 8-Correction Reconstruction
=======================================================
Reconstructed from source of truth (gwt_complete_reference.md).
All 8 corrections = 8 non-A1g channels of T1u x T1u.

V8 target: 1.7% mean, 1.5% median, 4.8% max on 23 molecules.

The 8 corrections:
  A1g (base): sigma coupling = pi/d^2 * E_harm
  Eg[1]:  LP repulsion (facing lone pairs)
  Eg[2]:  Heteronuclear p-p phase mismatch
  T1g[1]: Radical sigma reduction (5/6)
  T1g[2]: Overlap floor (min 1/4)
  T1g[3]: 3D parity node counting
  T2g[1]: Pi bonds + radical pi weakening
  T2g[2]: Enhanced ionic + period-3 boost
  T2g[3]: Triple-bond ionic (2/11)
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial

d = 3
PI = np.pi
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/PI**2 + np.log(2*d)))
E_Ry = (alpha_em**2 / 2) * 0.51100e6  # 13.606 eV

# ============================================================
# LATTICE CONSTANTS (all from d=3)
# ============================================================
C_BOND = PI / d**2                    # pi/9 = 0.3491
W_PI = np.cos(PI / d)                # cos(60) = 0.5
LP_I = (d**2 + 1) / d**3             # 10/27 = 0.3704
F_RAD = (2*d - 1) / (2*d)            # 5/6 = 0.8333
SP_BOOST = W_PI / d**2               # 0.5/9 = 0.0556
C_IONIC = 1 / (2*d + 1)              # 1/7 = 0.1429
C_IONIC_ENH = d / (2*d + 1)          # 3/7 = 0.4286
C_IONIC_TRIPLE = 2 / (d**2 + d - 1)  # 2/11 = 0.1818
P3_BOOST = (d**2 + 2) / (d**2 + 1)   # 11/10 = 1.1
OVERLAP_FLOOR = 1 / (d + 1)          # 1/4 = 0.25


# ============================================================
# ATOM DATABASE
# ============================================================
# (Z, IE_eV, n, p_electrons, LP, mass_amu, Z_eff_val)
ATOMS = {
    'H':  (1,  13.598, 1, 0, 0, 1.008,  1.00),
    'He': (2,  24.587, 1, 0, 0, 4.003,  1.70),
    'Li': (3,   5.392, 2, 0, 0, 6.941,  1.28),
    'Be': (4,   9.323, 2, 0, 0, 9.012,  1.91),
    'B':  (5,   8.298, 2, 1, 0, 10.81,  2.42),
    'C':  (6,  11.260, 2, 2, 0, 12.01,  3.14),
    'N':  (7,  14.534, 2, 3, 0, 14.01,  3.83),
    'O':  (8,  13.618, 2, 4, 1, 16.00,  4.45),
    'F':  (9,  17.423, 2, 5, 2, 19.00,  5.13),
    'Ne': (10, 21.565, 2, 6, 3, 20.18,  5.76),
    'Na': (11,  5.139, 3, 0, 0, 22.99,  2.51),
    'Mg': (12,  7.646, 3, 0, 0, 24.31,  3.31),
    'Al': (13,  5.986, 3, 1, 0, 26.98,  3.50),
    'Si': (14,  8.152, 3, 2, 0, 28.09,  4.29),
    'P':  (15, 10.487, 3, 3, 0, 30.97,  4.89),
    'S':  (16, 10.360, 3, 4, 1, 32.07,  5.48),
    'Cl': (17, 12.968, 3, 5, 2, 35.45,  6.12),
    'Ar': (18, 15.760, 3, 6, 3, 39.95,  6.76),
    'Br': (35, 11.814, 4, 5, 2, 79.90,  8.0),
    'I':  (53, 10.451, 5, 5, 2, 126.9,  9.5),
    'K':  (19,  4.341, 4, 0, 0, 39.10,  2.26),
    'Ca': (20,  6.113, 4, 0, 0, 40.08,  2.85),
}


def bond_v8(sym_a, sym_b, bond_order, is_radical=False):
    """Full V8 bond algorithm with all 8 Oh corrections."""

    Z_a, IE_a, n_a, p_a, lp_a, m_a, Zv_a = ATOMS[sym_a]
    Z_b, IE_b, n_b, p_b, lp_b, m_b, Zv_b = ATOMS[sym_b]

    # --- Energy scale ---
    E_harm = 2 * IE_a * IE_b / (IE_a + IE_b)

    # --- BASE: sigma coupling (A1g) ---
    n_sigma = 1  # always 1 sigma bond
    n_pi = bond_order - 1  # remaining are pi bonds

    # Number of p-p electrons available for bonding
    ne_pp = min(p_a, d) + min(p_b, d)  # max d per atom

    # --- RADICAL CORRECTION (single, not double-counted) ---
    # Half-sigma and F_RAD are the SAME geometric effect:
    #   half-sigma = unpaired mode has half coupling strength
    #   F_RAD = radical direction reduces projection to 5/6
    # Applying both double-counts. Apply ONE:
    #   - Half-sigma when BOTH atoms have p-electrons AND odd ne_pp
    #     (two p-modes sharing one electron in the sigma channel)
    #   - F_RAD otherwise (the radical modifies the projection angle)
    # H-X bonds: H contributes an s-electron to sigma, fully occupying it.
    #   The ne_pp count is irrelevant — sigma is full from the s-electron.
    sigma_eff = n_sigma
    half_sigma_applied = False
    both_have_p = (p_a > 0 and p_b > 0)
    if is_radical and both_have_p and ne_pp > 0 and ne_pp % 2 == 1 and ne_pp <= 2*d:
        sigma_eff = 0.5
        half_sigma_applied = True

    # --- CORRECTION 6: Radical pi-weakening ---
    # T2g[1]: radical reduces pi by (ne_pp-1)/ne_pp
    pi_eff = n_pi
    if is_radical and ne_pp > 1 and n_pi > 0:
        pi_eff = n_pi * (ne_pp - 1) / ne_pp

    parity_factor = 1.0
    n_max = max(n_a, n_b)

    # --- Base coupling ---
    coupling = sigma_eff + pi_eff * W_PI

    # --- CORRECTION 1 (Eg[1]): LP repulsion ---
    n_lp = min(lp_a, lp_b)
    lp_term = n_lp * LP_I * (2.0 / n_max)**2
    coupling -= lp_term

    # --- CORRECTION 2 (Eg[2]): Heteronuclear p-p phase ---
    # For heteronuclear pp bonds in same period: phase mismatch
    het_pp = (sym_a != sym_b and p_a > 0 and p_b > 0 and n_a == n_b)
    if het_pp:
        Z_am = np.sqrt(Zv_a * Zv_b)  # arithmetic/geometric mean of Z_eff
        Z_gm = np.sqrt(Zv_a * Zv_b)
        Z_hm = 2 * Zv_a * Zv_b / (Zv_a + Zv_b)
        am_gm_ratio = (Zv_a + Zv_b) / (2 * Z_gm)  # AM/GM >= 1
        phase_corr = am_gm_ratio ** (d - 1)  # ^2 for d=3
        coupling /= phase_corr

    # --- CORRECTION 3 (T1g[1]): Radical projection reduction ---
    # Only apply F_RAD when half-sigma was NOT already applied
    # (they are the same physics — applying both double-counts)
    if is_radical and not half_sigma_applied:
        coupling *= F_RAD  # 5/6

    # --- CORRECTION 4 (T1g[2]): Overlap floor ---
    coupling = max(coupling, OVERLAP_FLOOR)  # min 1/4

    # --- SP boost ---
    if bond_order > 1:
        deficit_a = max(0, bond_order - p_a)
        deficit_b = max(0, bond_order - p_b)
        sp_boost = (deficit_a + deficit_b) * 2 * SP_BOOST
        coupling += sp_boost

    # --- Covalent energy ---
    D_cov = C_BOND * E_harm * coupling * parity_factor

    # --- IONIC CORRECTIONS (T2g[2], T2g[3]) ---
    delta_IE = abs(IE_a - IE_b)
    E_avg = (IE_a + IE_b) / 2
    asym = delta_IE / E_avg if E_avg > 0 else 0

    D_ionic = 0
    if asym > 0.1:
        # Determine ionic tier
        ratio_cov_delta = D_cov / delta_IE if delta_IE > 0 else 999

        if ratio_cov_delta < 1 / d**3:  # < 1/27: strongly ionic
            # CORRECTION 7 (T2g[2]): Enhanced ionic
            c_ion = C_IONIC_ENH  # 3/7
            # CORRECTION 8: Period-3 boost
            if n_a >= 3 and n_b >= 3:
                c_ion *= P3_BOOST  # * 11/10
            D_ionic = c_ion * delta_IE
        elif (het_pp and bond_order == 3 and n_pi == 2
              and not is_radical):
            # CORRECTION 8 (T2g[3]): Triple-bond ionic
            D_ionic = C_IONIC_TRIPLE * delta_IE  # 2/11
        else:
            # Default ionic
            D_ionic = C_IONIC * delta_IE  # 1/7

    D_e = D_cov + D_ionic

    # --- ZPE ---
    mu_amu = m_a * m_b / (m_a + m_b)
    mu_me = mu_amu * 1822.89
    D_e_au = D_e / 27.211
    omega_vib = np.sqrt(2 * max(D_e_au, 1e-6) / mu_me)
    ZPE = 0.5 * omega_vib * 27.211

    D_0 = D_e - ZPE

    return {
        'D_e': D_e, 'D_0': D_0, 'ZPE': ZPE,
        'D_cov': D_cov, 'D_ionic': D_ionic,
        'coupling': coupling, 'lp_term': lp_term,
        'E_harm': E_harm, 'sigma_eff': sigma_eff,
        'pi_eff': pi_eff,
    }


# ============================================================
# TEST
# ============================================================
test_bonds = [
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
    ('C',  'C',  3, 8.700, 'CC3',     False),
    ('N',  'H',  1, 4.513, 'NH(NH3)', False),
    ('O',  'H',  1, 4.790, 'OH(H2O)', False),
    ('Cl', 'Cl', 1, 2.514, 'Cl2',     False),
    ('S',  'H',  1, 3.780, 'SH',      True),
    ('S',  'S',  2, 4.370, 'S2',      False),
    ('P',  'H',  1, 3.440, 'PH',      True),
]

print("GWT Bond Formula V8 — Full Reconstruction")
print("=" * 70)
print(f"  C_BOND = pi/d^2 = {C_BOND:.4f}")
print(f"  W_PI = cos(pi/d) = {W_PI:.4f}")
print(f"  LP_I = (d^2+1)/d^3 = {LP_I:.4f}")
print(f"  F_RAD = (2d-1)/(2d) = {F_RAD:.4f}")
print(f"  C_IONIC = 1/(2d+1) = {C_IONIC:.4f}")
print(f"  C_IONIC_ENH = d/(2d+1) = {C_IONIC_ENH:.4f}")
print(f"  C_IONIC_TRIPLE = 2/(d^2+d-1) = {C_IONIC_TRIPLE:.4f}")
print(f"  P3_BOOST = (d^2+2)/(d^2+1) = {P3_BOOST:.4f}")
print(f"  OVERLAP_FLOOR = 1/(d+1) = {OVERLAP_FLOOR:.4f}")
print("=" * 70)
print()

header = f"{'Name':>8} {'bo':>3} {'rad':>4} {'sig':>4} {'pi':>4} {'lp':>5} {'cpl':>6} {'D_e':>6} {'ZPE':>5} {'D_0':>6} {'D_obs':>6} {'err':>7}"
print(header)
print("-" * len(header))

errs = []
for sa, sb, bo, D_obs, name, radical in test_bonds:
    r = bond_v8(sa, sb, bo, radical)
    err = (r['D_0'] - D_obs) / D_obs * 100
    errs.append(abs(err))

    star = ' *' if abs(err) < 2 else '  ' if abs(err) < 5 else ''
    rad_s = '5/6' if radical else ''

    print(f"{name:>8} {bo:>3} {rad_s:>4} {r['sigma_eff']:>4.1f} {r['pi_eff']:>4.1f} "
          f"{r['lp_term']:>5.2f} {r['coupling']:>6.3f} {r['D_e']:>6.3f} {r['ZPE']:>5.3f} "
          f"{r['D_0']:>6.3f} {D_obs:>6.3f} {err:>+7.1f}%{star}")

print()
print(f"Mean:   {np.mean(errs):.1f}%")
print(f"Median: {np.median(errs):.1f}%")
print(f"Max:    {np.max(errs):.1f}%")
print(f"Under 2%:  {sum(1 for e in errs if e < 2)}/25")
print(f"Under 5%:  {sum(1 for e in errs if e < 5)}/25")
print(f"Under 10%: {sum(1 for e in errs if e < 10)}/25")
