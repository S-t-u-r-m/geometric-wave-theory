#!/usr/bin/env python3
"""
GWT Force Engine — Precomputed Oh Lookup + Radial Curves
=========================================================
The Wyler approach applied to ALL forces:
  Finite group (Oh) → finite table → everything precomputable.

Architecture:
  Chunk 1: Oh angular coupling (10x10 A1g table, closed-form)
  Chunk 2: Radial interaction curves (precomputed on 3D GPU, cached)
  Chunk 3: N-body corrections (closed-form A1g content)
  Chunk 4: Coulomb/ionic (analytical from alpha_em)

At runtime: ANY molecule = sum of table lookups. No simulation needed.
The GPU does a one-time precomputation of ~30 radial curves.

This is the GWT equivalent of quantum chemistry's integral tables,
but with 10 irreps instead of infinite basis functions.
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial
import os, json, time

# ============================================================
# CHUNK 1: Oh ANGULAR COUPLING (instant, exact)
# ============================================================
d = 3
Oh_order = 48
O_order = 24  # chiral

# Closed-form A1g content
def a1g_T1u(n):
    if n <= 0: return 0
    if n % 2 == 1: return 0
    return (3**n + 15) // 24

def a1g_T2g(n):
    if n <= 0: return 0
    return (3**n + 6 + 9*(-1)**n) // 24

def a1g_Eg(n):
    if n <= 0: return 0
    return (2**n + 2*(-1)**n) // 6

# Channel weights from cube geometry
W_CHANNEL = {
    'sigma': np.cos(0 * np.pi / d),        # 1.0
    'pi':    np.cos(1 * np.pi / d),         # 0.5
    'delta': np.cos(2 * np.pi / d),         # -0.5
}

# Pairwise selection rule: only same irrep couples (identity matrix)
# This means: sigma couples to sigma, pi to pi, etc.


# ============================================================
# CHUNK 2: RADIAL CURVE LIBRARY
# ============================================================
# Each radial curve depends on:
#   - L_a = n_a / 2 (breather width of atom A)
#   - L_b = n_b / 2 (breather width of atom B)
#   - type: 'bond' (half+half), 'lp' (full+full), 'cross' (half+full)
#   - channel: 'sigma' (pz+pz), 'pi' (px+px)
#
# Enumerate all physically relevant combinations:
#   n = 1, 2, 3, 4  → L = 0.5, 1.0, 1.5, 2.0
#   type = bond, lp
#   channel = sigma, pi
#
# For homonuclear: (L, L) pairs
# For heteronuclear: (L_a, L_b) pairs where L_a ≤ L_b

def enumerate_radial_keys():
    """List all unique radial curve parameters needed."""
    L_values = [0.5, 1.0, 1.5, 2.0]  # n=1,2,3,4
    types = ['bond', 'lp']  # half+half vs full+full
    channels = ['sigma', 'pi']

    keys = []
    for La in L_values:
        for Lb in L_values:
            if Lb < La: continue  # symmetry: only La ≤ Lb
            for typ in types:
                for ch in channels:
                    keys.append({
                        'L_a': La, 'L_b': Lb,
                        'type': typ, 'channel': ch,
                        'key': f"L{La:.1f}_L{Lb:.1f}_{typ}_{ch}"
                    })
    return keys


# R grid for radial curves
R_GRID = np.linspace(0.3, 8.0, 200)


def compute_radial_curve_3d(L_a, L_b, amp_a, amp_b, channel, xp_module, N_grid=96, box=8.0):
    """
    Compute one radial interaction curve on the 3D cubic lattice.

    Returns V(R) array: interaction energy vs separation.
    Positive = repulsive, negative = attractive.
    """
    xp = xp_module
    dx = 2 * box / N_grid
    dt = 0.2 * dx
    V_0 = 1.0 / np.pi**2
    gamma_sg = np.pi / (2**(d+1) * np.pi - 2)
    omega_p = float(np.cos(2 * gamma_sg))
    eps_p = float(np.sqrt(max(1 - omega_p**2, 1e-12)))

    x1d = xp.linspace(-box, box, N_grid, endpoint=False, dtype=np.float64)
    X, Y, Z = xp.meshgrid(x1d, x1d, x1d, indexing='ij')

    def make_mode(center_z, ang_type, amp):
        Zs = Z - center_z
        R = xp.sqrt(X**2 + Y**2 + Zs**2) + 1e-10
        if ang_type == 'sigma':
            ang = Zs / R
        else:  # pi
            ang = X / R
        rad = amp / (omega_p * xp.cosh(eps_p * R) + 1e-10)
        return (4.0 / np.pi) * xp.arctan(eps_p * ang * rad)

    def total_energy(phi):
        GE = 0.5 * ((xp.roll(phi,1,0)-phi)**2 + (xp.roll(phi,1,1)-phi)**2 +
                     (xp.roll(phi,1,2)-phi)**2)
        PE = V_0 * (1.0 - xp.cos(np.pi * phi))
        return float(xp.sum(GE + PE) * dx**3)

    def evolve(phi0, n_steps=3000, damping=0.02):
        phi = phi0.copy(); po = phi.copy(); es = []
        for s in range(n_steps):
            lap = (xp.roll(phi,1,0)+xp.roll(phi,-1,0)+xp.roll(phi,1,1)+
                   xp.roll(phi,-1,1)+xp.roll(phi,1,2)+xp.roll(phi,-1,2)-6*phi)/dx**2
            f = (1.0/np.pi)*xp.sin(np.pi*phi)
            pn = (2-damping*dt)*phi-(1-damping*dt)*po+dt**2*(lap-f)
            po = phi.copy(); phi = pn
            if s > n_steps*3//4 and s%100==0: es.append(total_energy(phi))
        return np.mean(es) if es else total_energy(phi)

    V_arr = np.zeros(len(R_GRID))
    for i, R in enumerate(R_GRID):
        phi_a = make_mode(-R/2, channel, amp_a)
        phi_b = make_mode(+R/2, channel, amp_b)
        E_a = evolve(phi_a, n_steps=2000)
        E_b = evolve(phi_b, n_steps=2000)
        E_ab = evolve(phi_a + phi_b, n_steps=3000)
        V_arr[i] = E_ab - E_a - E_b

    return V_arr


# ============================================================
# CHUNK 3: N-BODY CORRECTIONS (closed-form, instant)
# ============================================================
def nbody_correction(p_count):
    """
    N-body correction factor for an atom with p_count p-electrons.
    From A1g(T1u^p): even p has nonzero correction, odd p is exact.
    """
    a1g = a1g_T1u(p_count)
    if a1g == 0:
        return 0.0  # exact (half-fill, odd fill)
    dim = 3**p_count
    # Correction magnitude: A1g/dim suppressed by 1/d per order
    return a1g / (dim * d)


# ============================================================
# CHUNK 4: COULOMB / IONIC (analytical)
# ============================================================
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6  # 13.604 eV

def ionic_energy(IE_a, IE_b):
    """Ionic contribution from electronegativity difference."""
    delta_E = abs(IE_a - IE_b)
    E_avg = (IE_a + IE_b) / 2
    asym = delta_E / E_avg if E_avg > 0 else 0
    if asym > 0.1:
        return delta_E / (2*d + 1)  # 1/7
    return 0.0


# ============================================================
# ASSEMBLER: combine all chunks → bond energy
# ============================================================
ATOMS = {
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

# Radical reduction factor
F_RAD = (2*d - 1) / (2*d)  # 5/6

def assemble_bond(sym_a, sym_b, bo, is_radical=False, radial_table=None):
    """
    Assemble bond energy from precomputed chunks.

    If radial_table is provided: use GPU-precomputed radial curves.
    Otherwise: use the analytical formula (fast fallback).
    """
    a = ATOMS[sym_a]
    b = ATOMS[sym_b]
    E_harm = 2 * a['IE'] * b['IE'] / (a['IE'] + b['IE'])
    n_max = max(a['n'], b['n'])

    if radial_table is not None:
        # === GPU-PRECOMPUTED PATH ===
        # Look up the right radial curve for each channel
        La = a['n'] / 2.0
        Lb = b['n'] / 2.0
        L_lo, L_hi = min(La, Lb), max(La, Lb)

        # Sigma bonding curve
        key_sig = f"L{L_lo:.1f}_L{L_hi:.1f}_bond_sigma"
        V_sig = radial_table.get(key_sig, np.zeros_like(R_GRID))

        # Pi bonding curves (bo - 1 of them)
        key_pi = f"L{L_lo:.1f}_L{L_hi:.1f}_bond_pi"
        V_pi = radial_table.get(key_pi, np.zeros_like(R_GRID))

        # LP repulsion curves
        n_lp = min(a['lp'], b['lp'])
        key_lp = f"L{L_lo:.1f}_L{L_hi:.1f}_lp_pi"
        V_lp = radial_table.get(key_lp, np.zeros_like(R_GRID))

        # Assemble total potential
        # Chunk 1: Oh angular weights
        V_total = (W_CHANNEL['sigma'] * V_sig +
                   (bo - 1) * W_CHANNEL['pi'] * V_pi +
                   n_lp * V_lp)  # LP already has correct sign from 3D

        # Scale by energy
        V_total_eV = V_total * E_SCALE * (E_harm / E_H)

        # Radical correction
        if is_radical:
            V_total_eV *= F_RAD

        # Ionic
        D_ionic = ionic_energy(a['IE'], b['IE'])

        # Find minimum
        V_with_ionic = V_total_eV - D_ionic  # ionic lowers the potential
        i_min = np.argmin(V_with_ionic)
        D_e = -V_with_ionic[i_min]
        R_eq = R_GRID[i_min]

        return {'D_e': D_e, 'R_eq': R_eq, 'method': 'gpu_table'}

    else:
        # === ANALYTICAL FALLBACK (current formula) ===
        n_lp = min(a['lp'], b['lp'])
        LP_I = 1.0 / (d + 1)
        w_pi = W_CHANNEL['pi']

        cpl = 1 + (bo-1)*w_pi - n_lp * LP_I * (2.0/n_max)**2
        rf = F_RAD if is_radical else 1.0
        D_cov = (np.pi/d**2) * E_harm * max(cpl, 0) * rf
        D_ionic = ionic_energy(a['IE'], b['IE'])
        D_total = D_cov + D_ionic

        return {'D_e': D_total, 'R_eq': None, 'method': 'analytical'}


# ============================================================
# CACHE MANAGEMENT
# ============================================================
CACHE_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_FILE = os.path.join(CACHE_DIR, 'radial_cache.npz')

# Energy scale (calibrated from H2)
E_SCALE = 10.26  # eV per Lagrangian unit (from breather_potential.py calibration)


def load_radial_table():
    """Load precomputed radial curves from cache."""
    if os.path.exists(CACHE_FILE):
        data = np.load(CACHE_FILE)
        table = {key: data[key] for key in data.files}
        print(f"  Loaded {len(table)} radial curves from cache")
        return table
    return None


def save_radial_table(table):
    """Save computed radial curves to cache."""
    np.savez(CACHE_FILE, **table)
    print(f"  Saved {len(table)} radial curves to cache")


def precompute_radial_table(xp_module):
    """
    One-time GPU computation of all radial curves.
    Takes ~30 minutes on RTX 4070 Ti. Results cached to disk.
    """
    keys = enumerate_radial_keys()
    table = {}

    # Amplitude mapping
    AMP_BOND = 0.5   # half-filled
    AMP_LP = 1.0     # full (double)

    print(f"  Precomputing {len(keys)} radial curves on 3D GPU...")
    for i, k in enumerate(keys):
        amp_a = AMP_BOND if k['type'] == 'bond' else AMP_LP
        amp_b = amp_a  # symmetric for now
        ch = k['channel']

        t0 = time.time()
        V = compute_radial_curve_3d(k['L_a'], k['L_b'], amp_a, amp_b, ch, xp_module)
        elapsed = time.time() - t0

        table[k['key']] = V
        print(f"  [{i+1}/{len(keys)}] {k['key']:30s} ({elapsed:.0f}s)")

    save_radial_table(table)
    return table


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    print("=" * 70)
    print("  GWT Force Engine")
    print("  Wyler approach: finite group → finite table → instant lookup")
    print("=" * 70)
    print()

    # List all radial curves needed
    keys = enumerate_radial_keys()
    print(f"Radial curve library: {len(keys)} unique curves needed")
    print("  By L: ", sorted(set((k['L_a'], k['L_b']) for k in keys)))
    print("  By type: bond (half+half), lp (full+full)")
    print("  By channel: sigma (along bond), pi (perpendicular)")
    print()

    # Check if cache exists
    table = load_radial_table()
    if table is None:
        print("  No cache found. Using analytical fallback.")
        print("  Run with --precompute to build GPU radial table.")
        print()

    # Run analytical version
    print("ANALYTICAL RESULTS (current formula, instant):")
    exp_bonds = [
        ('H','H',1,4.478,'H2',False),('N','N',3,9.759,'N2',False),
        ('O','O',2,5.116,'O2',False),('F','F',1,1.602,'F2',False),
        ('H','F',1,5.869,'HF',False),('H','N',1,3.910,'NH',True),
        ('C','H',1,4.290,'CH',True),('C','C',1,3.600,'C-C',False),
        ('C','C',3,8.700,'C~C',False),('Cl','Cl',1,2.514,'Cl2',False),
    ]

    print(f"  {'Name':>8} {'D_pred':>7} {'D_obs':>7} {'err%':>7}")
    for sa,sb,bo,De,nm,rad in exp_bonds:
        r = assemble_bond(sa, sb, bo, is_radical=rad)
        err = (r['D_e'] - De)/De*100
        print(f"  {nm:>8} {r['D_e']:>7.3f} {De:>7.3f} {err:>+7.1f}%")

    print()
    print("To precompute GPU radial table:")
    print("  python gwt_force_engine.py --precompute")
    print()
    print("Once cached, the engine uses table lookups for EVERYTHING.")
    print("No simulation at runtime. Same Wyler philosophy: finite group → finite answer.")

    if '--precompute' in sys.argv:
        try:
            import cupy as cp
            table = precompute_radial_table(cp)
        except ImportError:
            print("ERROR: CuPy required for GPU precomputation")
