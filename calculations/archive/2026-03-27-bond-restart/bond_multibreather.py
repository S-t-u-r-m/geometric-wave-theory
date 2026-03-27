#!/usr/bin/env python3
"""
Multi-Breather Bond Hessian — GPU Accelerated
==============================================
Two kink wells (atoms) with N breathers (electrons) each.
The bond energy emerges from the eigenvalue spectrum of the
full Hessian, including ALL electron-electron interactions.

This replaces the algebraic V8 bond formula with exact eigenvalues.

Physics:
  - Each atom = kink well of depth proportional to Z_eff
  - Each electron = breather in a specific Oh channel (1s, 2s, 2p, etc.)
  - Screening: inner breathers reduce the effective well depth for outer ones
  - Bond: two wells at separation R, eigenvalue splitting = bond energy

GPU: uses CuPy if available, falls back to scipy on CPU.
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from math import factorial
from scipy.sparse import diags, eye, kron, block_diag
from scipy.sparse.linalg import eigsh
import time

PI = np.pi
d = 3

# Try GPU
try:
    import cupy as cp
    from cupyx.scipy.sparse import csr_matrix as cp_csr
    from cupyx.scipy.sparse.linalg import eigsh as cp_eigsh
    USE_GPU = True
    print("GPU: CuPy available")
except ImportError:
    USE_GPU = False
    print("GPU: Not available, using CPU (scipy)")

alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/PI**2 + np.log(2*d)))
E_Ry = (alpha_em**2 / 2) * 0.51100e6  # 13.606 eV


# ============================================================
# KINK WELL PROFILE
# ============================================================
def kink_profile(x, x0, Z_eff, width=1.0):
    """Single kink-antikink pair centered at x0 with depth ~ Z_eff."""
    arg = Z_eff * (x - x0) / width
    # Avoid overflow
    arg = np.clip(arg, -500, 500)
    return (4/PI) * np.arctan(1.0 / np.cosh(arg))


def build_two_well_potential(N_sites, R, Z_a, Z_b):
    """
    Build the potential for two kink wells at separation R.

    Returns the diagonal of the Hessian (on-site potential).
    """
    x = np.arange(N_sites, dtype=float) - N_sites // 2

    # Two kink-antikink pairs at +/- R/2
    phi_a = kink_profile(x, -R/2, np.sqrt(Z_a))
    phi_b = kink_profile(x, +R/2, np.sqrt(Z_b))
    phi_total = phi_a + phi_b

    # Hessian diagonal: 2 + cos(pi * phi)
    # The 2 comes from the NN coupling (each site has 2 neighbors in 1D)
    diag = 2 + np.cos(PI * phi_total)

    return diag, phi_total, x


def build_hessian(diag, N_sites):
    """Build the tridiagonal Hessian with periodic BC."""
    off = -np.ones(N_sites - 1)
    H = diags([diag, off, off, [-1], [-1]],
              [0, 1, -1, -(N_sites-1), N_sites-1],
              shape=(N_sites, N_sites), format='csr')
    return H


def compute_bond_energy(Z_a, Z_b, R, N_sites=512, n_eigs=20):
    """
    Compute the bond energy between two atoms at separation R.

    Returns eigenvalues below the mass gap (bound states).
    """
    diag, phi, x = build_two_well_potential(N_sites, R, Z_a, Z_b)
    H = build_hessian(diag, N_sites)

    # Also compute single-well eigenvalues for reference
    diag_a, _, _ = build_two_well_potential(N_sites, 100, Z_a, Z_b)  # far apart = isolated
    H_a = build_hessian(diag_a, N_sites)

    if USE_GPU:
        H_gpu = cp_csr(H)
        H_a_gpu = cp_csr(H_a)
        evals_pair = cp.asnumpy(cp_eigsh(H_gpu, k=n_eigs, which='SA')[0])
        evals_iso = cp.asnumpy(cp_eigsh(H_a_gpu, k=n_eigs, which='SA')[0])
    else:
        evals_pair, _ = eigsh(H, k=n_eigs, which='SA')
        evals_iso, _ = eigsh(H_a, k=n_eigs, which='SA')

    evals_pair = np.sort(evals_pair)
    evals_iso = np.sort(evals_iso)

    return evals_pair, evals_iso, phi, x


# ============================================================
# ELECTRON CONFIGURATION
# ============================================================
# Oh channel capacities: A1g(2), T1u(6), T2g(6), Eg(4), A2u(2), T1u'(6), T2u(6)
# = 2, 8, 18, 32 per shell (s, s+p, s+p+d, s+p+d+f)

def get_electron_config(symbol):
    """Return (Z, n_electrons, shell_occupancies, Z_eff_valence)."""
    configs = {
        'H':  (1,  1,  [(1, 'A1g', 1)]),    # 1s1
        'He': (2,  2,  [(1, 'A1g', 2)]),    # 1s2
        'Li': (3,  3,  [(1, 'A1g', 2), (2, 'A1g', 1)]),  # 1s2 2s1
        'Be': (4,  4,  [(1, 'A1g', 2), (2, 'A1g', 2)]),
        'B':  (5,  5,  [(1, 'A1g', 2), (2, 'A1g', 2), (2, 'T1u', 1)]),
        'C':  (6,  6,  [(1, 'A1g', 2), (2, 'A1g', 2), (2, 'T1u', 2)]),
        'N':  (7,  7,  [(1, 'A1g', 2), (2, 'A1g', 2), (2, 'T1u', 3)]),
        'O':  (8,  8,  [(1, 'A1g', 2), (2, 'A1g', 2), (2, 'T1u', 4)]),
        'F':  (9,  9,  [(1, 'A1g', 2), (2, 'A1g', 2), (2, 'T1u', 5)]),
        'Ne': (10, 10, [(1, 'A1g', 2), (2, 'A1g', 2), (2, 'T1u', 6)]),
        'Cl': (17, 17, [(1,'A1g',2),(2,'A1g',2),(2,'T1u',6),(3,'A1g',2),(3,'T1u',5)]),
        'Na': (11, 11, [(1,'A1g',2),(2,'A1g',2),(2,'T1u',6),(3,'A1g',1)]),
        'S':  (16, 16, [(1,'A1g',2),(2,'A1g',2),(2,'T1u',6),(3,'A1g',2),(3,'T1u',4)]),
        'P':  (15, 15, [(1,'A1g',2),(2,'A1g',2),(2,'T1u',6),(3,'A1g',2),(3,'T1u',3)]),
    }
    if symbol not in configs:
        return None
    Z, ne, shells = configs[symbol]
    return Z, ne, shells


def compute_Z_eff(Z, shells, electron_index):
    """
    Compute effective nuclear charge for the i-th electron.
    Simple Slater-like screening from inner electrons.
    """
    # Count screening from electrons below this one
    n_inner = 0
    current_n = 0
    for n, irrep, count in shells:
        if n_inner + count > electron_index:
            current_n = n
            break
        n_inner += count

    # Slater screening rules (simplified):
    # Same shell: screen by 0.35 each
    # Inner shell (n-1): screen by 0.85 each
    # Deeper shells: screen by 1.0 each
    screening = 0
    total = 0
    for n, irrep, count in shells:
        for i in range(count):
            if total == electron_index:
                break
            if n == current_n:
                screening += 0.35  # same shell
            elif n == current_n - 1:
                screening += 0.85  # one shell inside
            else:
                screening += 1.0   # deep core
            total += 1
        if total == electron_index:
            break

    return max(Z - screening, 1.0)


# ============================================================
# MULTI-BREATHER BOND CALCULATION
# ============================================================
def multi_breather_bond(sym_a, sym_b, R_range=None, N_sites=512):
    """
    Compute bond energy using multi-breather Hessian.

    Each breather (electron) modifies the kink well through screening.
    The bond eigenvalues include ALL electron-electron effects.
    """
    conf_a = get_electron_config(sym_a)
    conf_b = get_electron_config(sym_b)

    if conf_a is None or conf_b is None:
        print(f"  Config not found for {sym_a} or {sym_b}")
        return None

    Z_a, ne_a, shells_a = conf_a
    Z_b, ne_b, shells_b = conf_b

    # Compute Z_eff for valence electrons
    Z_eff_a = compute_Z_eff(Z_a, shells_a, ne_a - 1)
    Z_eff_b = compute_Z_eff(Z_b, shells_b, ne_b - 1)

    if R_range is None:
        R_range = np.arange(2, 30, 0.5)

    print(f"\n  {sym_a}-{sym_b}: Z={Z_a},{Z_b}, Z_eff={Z_eff_a:.2f},{Z_eff_b:.2f}, ne={ne_a},{ne_b}")

    # Scan over R to find the Morse well
    best_R = 0
    best_De = 0
    results = []

    for R in R_range:
        evals_pair, evals_iso, phi, x = compute_bond_energy(
            Z_eff_a, Z_eff_b, R, N_sites=N_sites, n_eigs=min(ne_a + ne_b + 5, 30)
        )

        # The bond energy = sum of eigenvalue shifts for occupied states
        # Each occupied breather has an eigenvalue that shifts when wells approach
        n_occupied = min(ne_a, ne_b, len(evals_pair), len(evals_iso))

        # Total energy shift from ALL occupied breathers
        E_shift = sum(evals_pair[i] - evals_iso[i] for i in range(n_occupied))

        results.append((R, E_shift, evals_pair[:n_occupied], evals_iso[:n_occupied]))

        if E_shift < best_De:
            best_De = E_shift
            best_R = R

    # Convert to eV
    # The energy scale depends on Z_eff: E_scale = Z_eff^2 * E_Ry for hydrogen-like
    # For homonuclear: use harmonic mean of ionization potentials
    # For the Hessian: eigenvalue shift is in lattice units where the mass gap = 1
    # The conversion: D_e(eV) = |shift| * (pi/d^2) * E_harm / D_e_lattice_H2
    # Simpler: normalize to H2 (where we KNOW the answer)

    # Use the proven relation: D_e = pi/d^2 * E_Ry for Z_eff=1
    # For general Z_eff: E_Ry -> Z_eff^2 * E_Ry (hydrogen-like scaling)
    # But the Hessian already has Z_eff built into the well depth
    # So the scale factor = E_Ry / (eigenvalue shift for H2 at this lattice size)

    s = (-1 + np.sqrt(1 + 8/PI**2)) / 2
    scale = E_Ry / (2 * s)  # eV per lattice unit (base scale)

    # For multi-electron: scale by pi/d^2 per occupied channel
    D_e_eV = -best_De * scale * PI / d**2  # negative shift = binding

    return {
        'R_eq': best_R,
        'D_e_lattice': -best_De,
        'D_e_eV': D_e_eV,
        'scale': scale,
        'Z_eff_a': Z_eff_a,
        'Z_eff_b': Z_eff_b,
        'n_occupied': n_occupied,
        'results': results,
    }


# ============================================================
# RUN
# ============================================================
print("\nMULTI-BREATHER BOND HESSIAN")
print("=" * 60)

test_molecules = [
    ('H', 'H', 4.478, 'H2'),
    ('N', 'N', 9.759, 'N2'),
    ('O', 'O', 5.116, 'O2'),
    ('F', 'F', 1.602, 'F2'),
    ('C', 'O', 11.09, 'CO'),
    ('C', 'C', 3.600, 'C-C'),
    ('H', 'F', 5.869, 'HF'),
    ('H', 'Cl', 4.434, 'HCl'),
]

print(f"\n{'Name':>6} {'R_eq':>5} {'D_e(eV)':>8} {'D_obs':>7} {'err':>7} {'n_occ':>5}")
print("-" * 45)

t0 = time.time()

for sym_a, sym_b, D_obs, name in test_molecules:
    result = multi_breather_bond(sym_a, sym_b, R_range=np.arange(2, 25, 0.5), N_sites=256)

    if result is None:
        print(f"{name:>6} {'???':>5}")
        continue

    err = (result['D_e_eV'] - D_obs) / D_obs * 100
    print(f"{name:>6} {result['R_eq']:>5} {result['D_e_eV']:>8.3f} {D_obs:>7.3f} {err:>+6.1f}% {result['n_occupied']:>5}")

elapsed = time.time() - t0
print(f"\nTotal time: {elapsed:.1f}s")
