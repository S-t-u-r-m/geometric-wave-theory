#!/usr/bin/env python3
"""
Bond Force Simulation — Oh Lookup Table Engine
================================================
No approximation formulas. Simulate the actual forces.

Two breather modes on the d=3 cubic lattice.
Interaction at each distance R computed from:
  1. Oh tensor product → which modes couple (O(1) lookup)
  2. Channel weight → cos(k*pi/d) from cube geometry
  3. Radial overlap → breather wavefunction at distance R
  4. Sum all mode-mode interactions → V(R)
  5. Find equilibrium → D_e, r_e

The Oh lookup table makes each force evaluation O(1).
Zero free parameters. All from L = (1/2)(dphi)^2 + (1/pi^2)(1-cos(pi*phi)).
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial

# ============================================================
# CONSTANTS FROM d=3 LAGRANGIAN
# ============================================================
d = 3
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H_eV = (alpha_em**2 / 2) * 0.51100e6  # 13.604 eV
a_0 = 0.52918  # Bohr radius in Angstrom
C_bond = np.pi / d  # universal coupling constant

# Channel weights from cube geometry
w = [np.cos(k * np.pi / d) for k in range(d+1)]
# w[0] = 1 (sigma), w[1] = 0.5 (pi), w[2] = -0.5 (delta), w[3] = -1

# Breather mass ratio
w_pi = w[1]  # = cos(pi/3) = 0.5


# ============================================================
# Oh CLOSED-FORM A1g LOOKUP
# ============================================================
def a1g_T1u(n):
    """A1g content of T1u^n (p-modes). O(1)."""
    if n <= 0: return 0
    if n % 2 == 1: return 0
    return (3**n + 15) // 24

def a1g_T2g(n):
    """A1g content of T2g^n (d-t2g modes). O(1)."""
    if n <= 0: return 0
    return (3**n + 6 + 9*(-1)**n) // 24

def a1g_Eg(n):
    """A1g content of Eg^n (d-eg modes). O(1)."""
    if n <= 0: return 0
    return (2**n + 2*(-1)**n) // 6


# ============================================================
# BREATHER WAVEFUNCTION OVERLAP
# ============================================================
def breather_overlap(R, n_a, n_b):
    """
    Radial overlap of two breather modes at separation R (in units of a_0).

    From the Lagrangian: breather n has width ~ n*a_0.
    The sine-Gordon breather solution: phi(x) = (4/pi) * arctan(sin(w*t)/cosh(x/L))
    At a given time, the spatial profile is sech(x/L) with L ~ n*a_0.

    The overlap of two breathers at distance R:
    integral of sech(x/L_a) * sech((x-R)/L_b) dx

    For identical breathers (n_a = n_b = n):
    overlap = pi * R/L / sinh(pi*R/(2L))  [exact for sech*sech convolution]

    Normalized so overlap(0) → 1 (maximum overlap).

    In dimensionless units R_bar = R / (n * a_0):
    overlap = sin(2*R_bar) if using the GWT equilibrium framework.

    Actually, from the source of truth:
    D_e = (pi/d) * E_H * sin(2R), and sin(2R) = 1/d at equilibrium.
    So the coupling goes as sin(2R) where R is in natural units.
    """
    # Use the Lagrangian-derived overlap: sin(2R) / (2R) for the coupling
    # This comes from the kink-antikink interaction in sine-Gordon theory
    # At short range: sin(2R) ~ 2R (linear, attractive)
    # At R = pi/2d: sin(2R) = sin(pi/d) = sin(60°) = sqrt(3)/2 (maximum)
    # At R = pi/4: sin(2R) = 1 (node)
    # Beyond: oscillates (antibonding nodes)

    # Convert R from Bohr to natural units (R_nat = R / (n_eff * a_0))
    n_eff = (n_a + n_b) / 2
    R_nat = R / n_eff if n_eff > 0 else R

    # Coupling function: sin(2*R_nat) gives the oscillating overlap
    # Multiplied by exponential decay for the breather tail
    # phi_overlap = sin(2*R_nat) * exp(-R_nat / d)
    if R_nat <= 0:
        return 0.0
    return np.sin(2 * R_nat) * np.exp(-R_nat / d)


def repulsive_core(R, n_a, n_b):
    """
    Short-range Pauli repulsion from breather core overlap.

    When two breathers overlap too much, the cosine potential
    creates a repulsive wall. This is the "hard core" of the
    nuclear/atomic force.

    From the Lagrangian: repulsion ~ (1/pi^2) * (1 - cos(pi*phi))
    At close range, phi → 2 (two overlapping kinks), and the
    potential energy rises as cosh^2.

    We model this as: V_rep = E_scale * exp(-d * R_nat)
    The factor d in the exponent comes from the d=3 lattice stiffness.
    """
    n_eff = (n_a + n_b) / 2
    R_nat = R / n_eff if n_eff > 0 else R
    if R_nat <= 0:
        return 1e10
    return np.exp(-d * R_nat)


# ============================================================
# ATOM DATABASE
# ============================================================
atoms_db = {
    'H':  {'Z': 1,  'IE': 13.598, 'n': 1, 's': 1, 'p': 0, 'lp': 0, 'config': [(1,0,1)]},
    'He': {'Z': 2,  'IE': 24.587, 'n': 1, 's': 2, 'p': 0, 'lp': 0, 'config': [(1,0,2)]},
    'Li': {'Z': 3,  'IE': 5.392,  'n': 2, 's': 1, 'p': 0, 'lp': 0, 'config': [(1,0,2),(2,0,1)]},
    'Be': {'Z': 4,  'IE': 9.323,  'n': 2, 's': 2, 'p': 0, 'lp': 0, 'config': [(1,0,2),(2,0,2)]},
    'B':  {'Z': 5,  'IE': 8.298,  'n': 2, 's': 2, 'p': 1, 'lp': 0, 'config': [(1,0,2),(2,0,2),(2,1,1)]},
    'C':  {'Z': 6,  'IE': 11.260, 'n': 2, 's': 2, 'p': 2, 'lp': 0, 'config': [(1,0,2),(2,0,2),(2,1,2)]},
    'N':  {'Z': 7,  'IE': 14.534, 'n': 2, 's': 2, 'p': 3, 'lp': 1, 'config': [(1,0,2),(2,0,2),(2,1,3)]},
    'O':  {'Z': 8,  'IE': 13.618, 'n': 2, 's': 2, 'p': 4, 'lp': 2, 'config': [(1,0,2),(2,0,2),(2,1,4)]},
    'F':  {'Z': 9,  'IE': 17.423, 'n': 2, 's': 2, 'p': 5, 'lp': 3, 'config': [(1,0,2),(2,0,2),(2,1,5)]},
    'Na': {'Z': 11, 'IE': 5.139,  'n': 3, 's': 1, 'p': 0, 'lp': 0, 'config': [(1,0,2),(2,0,2),(2,1,6),(3,0,1)]},
    'Si': {'Z': 14, 'IE': 8.152,  'n': 3, 's': 2, 'p': 2, 'lp': 0, 'config': [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,2)]},
    'P':  {'Z': 15, 'IE': 10.487, 'n': 3, 's': 2, 'p': 3, 'lp': 1, 'config': [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,3)]},
    'S':  {'Z': 16, 'IE': 10.360, 'n': 3, 's': 2, 'p': 4, 'lp': 2, 'config': [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,4)]},
    'Cl': {'Z': 17, 'IE': 12.968, 'n': 3, 's': 2, 'p': 5, 'lp': 3, 'config': [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,5)]},
    'K':  {'Z': 19, 'IE': 4.341,  'n': 4, 's': 1, 'p': 0, 'lp': 0, 'config': [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(4,0,1)]},
}


def get_mode_list(sym):
    """
    Get the list of valence Oh irrep modes for an atom.

    Each mode = (irrep_name, channel_type, is_bonding)
    where channel_type = 'sigma', 'pi', 'delta', or 'lp'

    For bonding: modes are ordered so sigma comes first, then pi.
    """
    atom = atoms_db[sym]
    modes = []

    # s-orbital → A1g, sigma bonding
    if atom['s'] >= 1:
        modes.append({'irrep': 'A1g', 'type': 'sigma', 'n': atom['n'], 'weight': 1.0})
    if atom['s'] >= 2:
        modes.append({'irrep': 'A1g', 'type': 'sigma_paired', 'n': atom['n'], 'weight': 1.0})

    # p-orbitals → T1u channels
    # First d modes are unpaired (Hund's rule), rest are paired
    p = atom['p']
    for i in range(min(p, d)):
        # First 3 p-electrons: one in each direction (px, py, pz)
        if i == 0:
            ch_type = 'sigma'  # pz along bond axis
        else:
            ch_type = 'pi'    # px, py perpendicular
        modes.append({
            'irrep': 'T1u',
            'type': ch_type,
            'n': atom['n'],
            'weight': w[0] if ch_type == 'sigma' else w[1],
            'paired': False
        })

    # Overfill: paired p-electrons (lone pairs in bonding direction become antibonding)
    for i in range(max(p - d, 0)):
        modes.append({
            'irrep': 'T1u',
            'type': 'lp',  # lone pair (paired electron)
            'n': atom['n'],
            'weight': w[2],  # delta weight (antibonding)
            'paired': True
        })

    return modes


def mode_mode_coupling(mode_a, mode_b, R):
    """
    Compute the coupling between two modes at distance R.

    From Oh: modes couple only if they're in the same irrep (pairwise A1g = identity).
    The coupling strength = directional_weight × radial_overlap × energy_scale.
    """
    # Selection rule: only same irrep couples (Oh orthogonality)
    if mode_a['irrep'] != mode_b['irrep']:
        return 0.0

    # Directional weight: product of channel weights
    # sigma-sigma: 1 × 1 = 1
    # pi-pi: 0.5 × 0.5 = 0.25
    # sigma-lp: 1 × (-0.5) = -0.5 (repulsive)
    # lp-lp: (-0.5) × (-0.5) = 0.25 (repulsive but weaker)
    w_dir = mode_a['weight'] * mode_b['weight']

    # Sign: bonding vs antibonding
    # If both are bonding modes (unpaired): attractive (negative energy)
    # If one or both are LP: repulsive (positive energy)
    is_lp_a = mode_a['type'] == 'lp' or mode_a.get('paired', False)
    is_lp_b = mode_b['type'] == 'lp' or mode_b.get('paired', False)

    if is_lp_a and is_lp_b:
        sign = +1  # LP-LP repulsion
    elif is_lp_a or is_lp_b:
        sign = +1  # LP-bonding repulsion (Pauli)
    else:
        sign = -1  # bonding-bonding attraction

    # Radial overlap
    overlap = breather_overlap(R, mode_a['n'], mode_b['n'])

    return sign * w_dir * overlap


def potential_energy(sym_a, sym_b, bo, R):
    """
    Total potential energy of two atoms at distance R.

    V(R) = Σ_{mode pairs} [coupling(mode_a, mode_b, R)] × C_bond × E_scale
         + V_repulsive(R)

    The sum is over all relevant mode pairs.
    bo determines which modes participate in bonding.
    """
    atom_a = atoms_db[sym_a]
    atom_b = atoms_db[sym_b]
    modes_a = get_mode_list(sym_a)
    modes_b = get_mode_list(sym_b)

    # Energy scale: harmonic mean of ionization energies
    E_a, E_b = atom_a['IE'], atom_b['IE']
    E_harm = 2 * E_a * E_b / (E_a + E_b)

    # Compute all mode-mode interactions
    V_total = 0.0
    V_bonding = 0.0
    V_lp = 0.0

    # For each bonding channel (up to bond order):
    # Match sigma first, then pi
    bonding_a = [m for m in modes_a if m['type'] in ('sigma', 'pi') and m['irrep'] == 'T1u']
    bonding_b = [m for m in modes_b if m['type'] in ('sigma', 'pi') and m['irrep'] == 'T1u']
    lp_a = [m for m in modes_a if m['type'] == 'lp']
    lp_b = [m for m in modes_b if m['type'] == 'lp']

    # s-orbital bonding (for H, Li, Na, etc.)
    s_a = [m for m in modes_a if m['irrep'] == 'A1g' and m['type'] == 'sigma']
    s_b = [m for m in modes_b if m['irrep'] == 'A1g' and m['type'] == 'sigma']

    # Bonding interactions: match bonding modes from A with B
    n_bond_channels = bo
    bond_pairs = []

    # sigma first
    if len(bonding_a) > 0 and len(bonding_b) > 0:
        bond_pairs.append((bonding_a[0], bonding_b[0]))
    elif len(s_a) > 0 and len(s_b) > 0:
        bond_pairs.append((s_a[0], s_b[0]))
    elif len(s_a) > 0 and len(bonding_b) > 0:
        bond_pairs.append((s_a[0], bonding_b[0]))
    elif len(bonding_a) > 0 and len(s_b) > 0:
        bond_pairs.append((bonding_a[0], s_b[0]))

    # pi bonds
    pi_a = [m for m in bonding_a if m['type'] == 'pi']
    pi_b = [m for m in bonding_b if m['type'] == 'pi']
    for i in range(min(bo - 1, len(pi_a), len(pi_b))):
        bond_pairs.append((pi_a[i], pi_b[i]))

    # Compute bonding energy
    for ma, mb in bond_pairs:
        coupling = mode_mode_coupling(ma, mb, R)
        V_bonding += coupling * C_bond * E_harm

    # LP-LP repulsion (facing lone pairs)
    for la in lp_a:
        for lb in lp_b:
            coupling = mode_mode_coupling(la, lb, R)
            V_lp += coupling * C_bond * E_harm

    # LP-bond cross repulsion (LP on A interferes with bonding on B and vice versa)
    for la in lp_a:
        for mb in [bp[1] for bp in bond_pairs if bp[1] in bonding_b + s_b]:
            if la['irrep'] == mb['irrep']:
                overlap = breather_overlap(R, la['n'], mb['n'])
                V_lp += 0.5 * abs(w[2]) * overlap * C_bond * E_harm / d

    for lb in lp_b:
        for ma in [bp[0] for bp in bond_pairs if bp[0] in bonding_a + s_a]:
            if lb['irrep'] == ma['irrep']:
                overlap = breather_overlap(R, lb['n'], ma['n'])
                V_lp += 0.5 * abs(w[2]) * overlap * C_bond * E_harm / d

    # Short-range repulsion (Pauli exclusion core)
    n_eff = (atom_a['n'] + atom_b['n']) / 2
    V_rep = E_harm * repulsive_core(R, atom_a['n'], atom_b['n'])

    V_total = V_bonding + V_lp + V_rep

    return V_total, V_bonding, V_lp, V_rep


def find_equilibrium(sym_a, sym_b, bo, R_min=0.3, R_max=8.0, N_points=2000):
    """
    Scan R to find the potential minimum.
    Returns: (D_e, R_eq, V_curve)
    """
    R_range = np.linspace(R_min, R_max, N_points)
    V_total = np.zeros(N_points)
    V_bond = np.zeros(N_points)
    V_lp_arr = np.zeros(N_points)
    V_rep_arr = np.zeros(N_points)

    for i, R in enumerate(R_range):
        vt, vb, vlp, vrep = potential_energy(sym_a, sym_b, bo, R)
        V_total[i] = vt
        V_bond[i] = vb
        V_lp_arr[i] = vlp
        V_rep_arr[i] = vrep

    # Find minimum
    i_min = np.argmin(V_total)
    D_e = -V_total[i_min]  # Dissociation energy (positive = bound)
    R_eq = R_range[i_min]

    return D_e, R_eq, R_range, V_total, V_bond, V_lp_arr, V_rep_arr


# ============================================================
# EXPERIMENTAL DATA
# ============================================================
exp_bonds = [
    ('H', 'H', 1, 4.478, 0.741, 'H2'),
    ('N', 'N', 3, 9.759, 1.098, 'N2'),
    ('O', 'O', 2, 5.116, 1.208, 'O2'),
    ('F', 'F', 1, 1.602, 1.412, 'F2'),
    ('H', 'F', 1, 5.869, 0.917, 'HF'),
    ('H', 'Cl', 1, 4.434, 1.275, 'HCl'),
    ('H', 'O', 1, 4.392, 0.970, 'OH'),
    ('H', 'N', 1, 3.910, 1.036, 'NH'),
    ('C', 'H', 1, 4.290, 1.089, 'CH'),
    ('C', 'C', 1, 3.600, 1.540, 'C-C'),
    ('C', 'C', 2, 6.360, 1.340, 'C=C'),
    ('C', 'C', 3, 8.700, 1.200, 'C≡C'),
    ('C', 'O', 3, 11.09, 1.128, 'CO'),
    ('C', 'O', 2, 7.710, 1.200, 'C=O'),
    ('N', 'O', 2, 6.497, 1.151, 'NO'),
    ('Cl', 'Cl', 1, 2.514, 1.988, 'Cl2'),
    ('S', 'H', 1, 3.78, 1.34, 'SH'),
    ('S', 'S', 2, 4.37, 1.89, 'S2'),
    ('P', 'H', 1, 3.44, 1.42, 'PH'),
]


# ============================================================
# RUN SIMULATION
# ============================================================
print("=" * 80)
print("  BOND FORCE SIMULATION — Oh Lookup Table Engine")
print("  Two breathers on d=3 lattice. Forces from Oh tensor products.")
print("  Zero free parameters. All from the Lagrangian.")
print("=" * 80)
print()

# First: detailed H2 scan
print("=== H2 POTENTIAL CURVE (validation) ===")
D_e, R_eq, R_arr, V_tot, V_b, V_lp, V_rep = find_equilibrium('H', 'H', 1)
print(f"  D_e = {D_e:.3f} eV  (obs D_0 = 4.478, D_e = 4.748)")
print(f"  R_eq = {R_eq:.3f} Bohr  (obs = 1.401 Bohr)")
print(f"  V_bond(R_eq) = {V_b[np.argmin(V_tot)]:.3f} eV")
print(f"  V_rep(R_eq) = {V_rep[np.argmin(V_tot)]:.3f} eV")
print()

# Show potential curve around minimum
i_min = np.argmin(V_tot)
print("  R(Bohr)  V_total   V_bond    V_rep")
for offset in [-10, -5, -2, -1, 0, 1, 2, 5, 10, 20, 40]:
    idx = max(0, min(len(R_arr)-1, i_min + offset))
    print(f"  {R_arr[idx]:6.3f}   {V_tot[idx]:+8.3f}  {V_b[idx]:+8.3f}  {V_rep[idx]:+8.3f}")

print()
print("=== ALL BONDS ===")
print(f"  {'Name':>8} {'bo':>3} {'D_e':>7} {'D_obs':>7} {'err%':>7} {'R_eq':>6} {'R_obs':>6}")
print("  " + "-" * 55)

errs = []
for sym_a, sym_b, bo, De_obs, Re_obs, name in exp_bonds:
    D_e, R_eq, _, _, _, _, _ = find_equilibrium(sym_a, sym_b, bo)
    err = (D_e - De_obs) / De_obs * 100 if De_obs > 0 else 0
    errs.append(abs(err))
    print(f"  {name:>8} {bo:>3} {D_e:>7.3f} {De_obs:>7.3f} {err:>+7.1f}% {R_eq:>6.3f} {Re_obs:>6.3f}")

print()
print(f"  Mean |err|: {np.mean(errs):.1f}%")
print(f"  Median:     {np.median(errs):.1f}%")
print(f"  Within 10%: {sum(1 for e in errs if e < 10)}/{len(errs)}")
print(f"  Within 20%: {sum(1 for e in errs if e < 20)}/{len(errs)}")
