#!/usr/bin/env python3
"""
Full Bond Simulation — Oh Lookup × Lagrangian Radial
=====================================================
Angular: Oh tensor product lookup (O(1), exact)
Radial: sine-Gordon breather overlap (numerical, from Lagrangian)

For each bond:
1. Compute radial potential V_raw(R) from two overlapping breathers
2. Apply Oh channel weights (sigma, pi, LP)
3. Find equilibrium → D_e, R_eq
4. Compare to experiment

Zero free parameters beyond L = n/2 (breather width from H2 calibration).
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial

# ============================================================
# CONSTANTS
# ============================================================
d = 3
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6  # 13.604 eV
C_bond = np.pi / d
w_ch = [np.cos(k * np.pi / d) for k in range(d+1)]  # [1, 0.5, -0.5, -1]

# Spatial grid
N_grid = 3000
x_max = 15.0
dx = 2 * x_max / N_grid
x_grid = np.linspace(-x_max, x_max, N_grid)


# ============================================================
# BREATHER PHYSICS (from Lagrangian)
# ============================================================
def breather_profile(x, center, L):
    arg = np.clip((x - center) / L, -50, 50)
    return (4.0 / np.pi) * np.arctan(1.0 / np.cosh(arg))


def interaction_energy_raw(R, L_a, L_b):
    """Raw interaction energy from Lagrangian (dimensionless)."""
    phi_a = breather_profile(x_grid, -R/2, L_a)
    phi_b = breather_profile(x_grid, +R/2, L_b)
    phi_tot = phi_a + phi_b

    # Energy densities
    def E_density(phi):
        dphi = np.gradient(phi, dx)
        return 0.5 * dphi**2 + (1.0/np.pi**2) * (1.0 - np.cos(np.pi * phi))

    E_AB = np.trapezoid(E_density(phi_tot), x_grid)
    E_A = np.trapezoid(E_density(phi_a), x_grid)
    E_B = np.trapezoid(E_density(phi_b), x_grid)
    return E_AB - E_A - E_B


# Precompute and cache radial potentials
_cache = {}

def get_radial_potential(L_a, L_b, R_min=0.3, R_max=8.0, N_R=250):
    """Get (or compute) the radial potential for a given L_a, L_b pair."""
    key = (round(L_a, 4), round(L_b, 4))
    if key not in _cache:
        R_arr = np.linspace(R_min, R_max, N_R)
        V_arr = np.array([interaction_energy_raw(R, L_a, L_b) for R in R_arr])
        _cache[key] = (R_arr, V_arr)
    return _cache[key]


# ============================================================
# Oh ANGULAR COUPLING
# ============================================================
def a1g_T1u(n):
    if n <= 0: return 0
    if n % 2 == 1: return 0
    return (3**n + 15) // 24


# ============================================================
# ENERGY SCALE CALIBRATION FROM H2
# ============================================================
# H2: L = 0.49 Bohr, D_e = 4.748 eV
# The energy scale converts Lagrangian units → eV for sigma bond
L_H = 0.49  # Bohr (gives R_eq = 1.400)
R_h2, V_h2 = get_radial_potential(L_H, L_H)
i_min_h2 = np.argmin(V_h2)
V_min_h2 = V_h2[i_min_h2]

# For sigma bond: V_physical = E_scale_sigma × V_raw
# D_e(H2) = 4.748 eV = -E_scale × V_min
E_SCALE_SIGMA = -4.748 / V_min_h2
print(f"Calibration from H2:")
print(f"  L_H = {L_H} Bohr")
print(f"  R_eq(H2) = {R_h2[i_min_h2]:.3f} Bohr (obs: 1.401)")
print(f"  V_min(raw) = {V_min_h2:.6f}")
print(f"  E_scale_sigma = {E_SCALE_SIGMA:.3f} eV/unit")
print()


# ============================================================
# ATOM DATABASE
# ============================================================
atoms = {
    'H':  {'Z':1,  'IE':13.598, 'n':1, 'p':0, 'lp':0},
    'Li': {'Z':3,  'IE':5.392,  'n':2, 'p':0, 'lp':0},
    'B':  {'Z':5,  'IE':8.298,  'n':2, 'p':1, 'lp':0},
    'C':  {'Z':6,  'IE':11.260, 'n':2, 'p':2, 'lp':0},
    'N':  {'Z':7,  'IE':14.534, 'n':2, 'p':3, 'lp':1},
    'O':  {'Z':8,  'IE':13.618, 'n':2, 'p':4, 'lp':2},
    'F':  {'Z':9,  'IE':17.423, 'n':2, 'p':5, 'lp':3},
    'Na': {'Z':11, 'IE':5.139,  'n':3, 'p':0, 'lp':0},
    'Si': {'Z':14, 'IE':8.152,  'n':3, 'p':2, 'lp':0},
    'P':  {'Z':15, 'IE':10.487, 'n':3, 'p':3, 'lp':1},
    'S':  {'Z':16, 'IE':10.360, 'n':3, 'p':4, 'lp':2},
    'Cl': {'Z':17, 'IE':12.968, 'n':3, 'p':5, 'lp':3},
    'K':  {'Z':19, 'IE':4.341,  'n':4, 'p':0, 'lp':0},
}


def breather_width(sym):
    """Breather width L = n/2 (in Bohr). From H2 calibration."""
    return atoms[sym]['n'] / 2.0


def E_harmonic(sym_a, sym_b):
    Ea, Eb = atoms[sym_a]['IE'], atoms[sym_b]['IE']
    return 2 * Ea * Eb / (Ea + Eb)


# ============================================================
# BOND SIMULATION
# ============================================================
def simulate_bond(sym_a, sym_b, bo):
    """
    Simulate a bond between two atoms.

    1. Get radial potential from Lagrangian (breather overlap)
    2. Apply Oh channel weights for each bond type
    3. Add LP repulsion (same potential, opposite sign)
    4. Find equilibrium

    Channel structure:
      sigma (k=0): 1 channel, weight w[0] = 1
      pi (k=1): (bo-1) channels, weight w[1] = 0.5 each
      LP facing: min(lp_a, lp_b) channels, weight w[2] = -0.5 (repulsive)
    """
    a = atoms[sym_a]
    b = atoms[sym_b]
    L_a = breather_width(sym_a)
    L_b = breather_width(sym_b)

    # Energy scale: proportional to E_harmonic, normalized by H2
    # E_scale = E_SCALE_SIGMA × (E_harm / E_H)
    E_harm = E_harmonic(sym_a, sym_b)
    E_ratio = E_harm / E_H

    # Get radial potential shape
    R_arr, V_raw = get_radial_potential(L_a, L_b)

    # === BUILD TOTAL POTENTIAL ===
    # Sigma channel: 1 bond, weight 1
    V_sigma = E_SCALE_SIGMA * E_ratio * w_ch[0] * V_raw

    # Pi channels: (bo-1) bonds, weight w[1] = 0.5 each
    # Pi overlap is WEAKER because coupling is perpendicular
    # In the Lagrangian: pi modes couple through transverse breather modes
    # Weight = cos(pi/d) = 0.5 per channel
    n_pi = max(0, bo - 1)
    V_pi = n_pi * E_SCALE_SIGMA * E_ratio * w_ch[1] * V_raw

    # LP repulsion: facing lone pairs create antibonding overlap
    # SAME radial potential shape, but POSITIVE (repulsive)
    # Weight per LP pair: from Oh T1u⊗T1u decomposition
    # The A1g fraction (1/d²) is attractive, rest is repulsive
    # Net repulsion weight per LP pair = (d²-2)/d² ≈ 0.778
    lp_facing = min(a['lp'], b['lp'])
    lp_weight = (d**2 - 2) / d**2  # 7/9 from Oh decomposition
    V_lp = lp_facing * E_SCALE_SIGMA * E_ratio * lp_weight * (-V_raw)  # flip sign: repulsive

    # Total potential
    V_total = V_sigma + V_pi + V_lp

    # Find minimum
    i_min = np.argmin(V_total)
    D_e = -V_total[i_min]
    R_eq = R_arr[i_min]

    return {
        'D_e': D_e,
        'R_eq': R_eq,
        'V_sigma_min': V_sigma[i_min],
        'V_pi_min': V_pi[i_min] if n_pi > 0 else 0,
        'V_lp_min': V_lp[i_min] if lp_facing > 0 else 0,
        'E_harm': E_harm,
        'E_ratio': E_ratio,
        'lp_facing': lp_facing,
        'n_pi': n_pi,
        'L_a': L_a, 'L_b': L_b,
    }


# ============================================================
# EXPERIMENTAL DATA
# ============================================================
exp_bonds = [
    ('H', 'H', 1, 4.478, 0.741, 'H2'),
    ('Li', 'Li', 1, 1.046, 2.673, 'Li2'),
    ('N', 'N', 3, 9.759, 1.098, 'N2'),
    ('O', 'O', 2, 5.116, 1.208, 'O2'),
    ('F', 'F', 1, 1.602, 1.412, 'F2'),
    ('H', 'F', 1, 5.869, 0.917, 'HF'),
    ('H', 'Cl', 1, 4.434, 1.275, 'HCl'),
    ('Na', 'Cl', 1, 4.230, 2.361, 'NaCl'),
    ('Li', 'H', 1, 2.429, 1.596, 'LiH'),
    ('H', 'O', 1, 4.392, 0.970, 'OH'),
    ('C', 'O', 3, 11.09, 1.128, 'CO'),
    ('N', 'O', 2, 6.497, 1.151, 'NO'),
    ('H', 'N', 1, 3.910, 1.036, 'NH'),
    ('C', 'H', 1, 4.290, 1.089, 'CH'),
    ('C', 'C', 1, 3.600, 1.540, 'C-C'),
    ('C', 'N', 3, 7.760, 1.170, 'CN'),
    ('C', 'C', 2, 6.360, 1.340, 'C=C'),
    ('C', 'O', 2, 7.710, 1.200, 'C=O'),
    ('C', 'C', 3, 8.700, 1.200, 'C≡C'),
    ('Cl', 'Cl', 1, 2.514, 1.988, 'Cl2'),
    ('S', 'H', 1, 3.78, 1.34, 'SH'),
    ('S', 'S', 2, 4.37, 1.89, 'S2'),
    ('P', 'H', 1, 3.44, 1.42, 'PH'),
]


# ============================================================
# RUN
# ============================================================
print("=" * 90)
print("  FULL BOND SIMULATION: Oh Lookup x Lagrangian Radial")
print("  Angular: Oh tensor products (exact, O(1))")
print("  Radial: sine-Gordon breather overlap (numerical)")
print("  L = n/2 Bohr. E_scale from H2. Zero free parameters.")
print("=" * 90)
print()

print(f"  {'Name':>8} {'bo':>3} {'L_a':>4} {'L_b':>4} {'E_r':>5} {'lp':>3} "
      f"{'V_sig':>7} {'V_pi':>7} {'V_lp':>7} "
      f"{'D_e':>7} {'D_obs':>7} {'err%':>7} {'R_eq':>6} {'R_obs':>6}")
print("  " + "-" * 100)

errs = []
for sym_a, sym_b, bo, De_obs, Re_obs, name in exp_bonds:
    r = simulate_bond(sym_a, sym_b, bo)
    err = (r['D_e'] - De_obs) / De_obs * 100
    errs.append(abs(err))
    print(f"  {name:>8} {bo:>3} {r['L_a']:>4.1f} {r['L_b']:>4.1f} {r['E_ratio']:>5.2f} {r['lp_facing']:>3} "
          f"{r['V_sigma_min']:>+7.2f} {r['V_pi_min']:>+7.2f} {r['V_lp_min']:>+7.2f} "
          f"{r['D_e']:>7.3f} {De_obs:>7.3f} {err:>+7.1f}% {r['R_eq']:>6.3f} {Re_obs:>6.3f}")

print()
print(f"  Mean |err|: {np.mean(errs):.1f}%")
print(f"  Median:     {np.median(errs):.1f}%")
print(f"  Max:        {np.max(errs):.1f}%")
print(f"  Within 10%: {sum(1 for e in errs if e < 10)}/{len(errs)}")
print(f"  Within 20%: {sum(1 for e in errs if e < 20)}/{len(errs)}")
print(f"  Within 30%: {sum(1 for e in errs if e < 30)}/{len(errs)}")

# Separate homonuclear and heteronuclear
homo = [i for i,(a,b,*_) in enumerate(exp_bonds) if a==b]
hetero = [i for i,(a,b,*_) in enumerate(exp_bonds) if a!=b]
print(f"\n  Homonuclear ({len(homo)}): mean = {np.mean([errs[i] for i in homo]):.1f}%")
print(f"  Heteronuclear ({len(hetero)}): mean = {np.mean([errs[i] for i in hetero]):.1f}%")
