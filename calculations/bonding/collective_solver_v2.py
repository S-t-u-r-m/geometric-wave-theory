#!/usr/bin/env python3
"""
Collective Field Solver — v2 (1D, multi-mode SCF with QM mean-field backreaction)
=================================================================================

What's new vs v1
-----------------
v1 was conceptually wrong: it added the bonding-mode amplitude back into
phi_bg as if it were a classical field shift. In QM, an occupied bound
mode has <delta_phi> = 0; the SCF correction lives in <delta_phi^2>,
which modifies the second derivative of V (the Hessian diagonal), not
the background field itself.

Correct mean-field SCF
----------------------
For sine-Gordon V(phi) = (1/pi^2)(1 - cos(pi*phi)):
  <V''(phi_bg + delta_phi)> = <cos(pi*(phi_bg+delta_phi))>
                            = cos(pi*phi_bg) * <cos(pi*delta_phi)>
                            ~= cos(pi*phi_bg) * exp(-pi^2 <delta_phi^2> / 2)

For bound mode n with occupation N_n and frequency omega_n=sqrt(w2_n):
  <delta_phi^2>(x) += (N_n + 1/2) * psi_n^2(x) / (2 omega_n)

Iterate:
  1. Build H with diag = 2 + cos(pi*phi_bg) * dressing_factor(x)
  2. Diagonalize -> w2_n, psi_n
  3. Update <delta_phi^2>(x) from occupied modes
  4. Repeat until lowest eigenvalues converge

Then bond energy uses bond_3d_emerge.py's convention:
  V(R) = w2_min_dimer(R) - w2_min_isolated_atom    (in lattice units)
  D_e_eV = -min_R V(R) * E_harm / (2*s_PT)

Test molecules
--------------
H2:  occupy mode 0 with 2 electrons (sigma).
     Bare predicts 4.79 eV (obs 4.748). v2 should not move much from
     bare since H2 has no LP electrons.

Cl2: occupy mode 0 (sigma, 2e), mode 1 (LP, 4e), mode 2 (LP, 4e).
     Bare predicts 4.51 eV; observed 2.51 eV. If LP backreaction
     creates repulsion, D_e should drop. Acid test.

F2:  same as Cl2 but E_harm higher.
     Bare predicts 13.72 eV; observed 1.60 eV. Most extreme test.
"""
import sys
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

PI = np.pi
d = 3
S_PT = (-1 + np.sqrt(1 + 8/PI**2)) / 2
KW = 3

# ============================================================
# FIELD CONSTRUCTION (same as v1)
# ============================================================
def kink(x, x0):
    return (4.0/PI) * np.arctan(np.exp(x - x0))
def antikink(x, x0):
    return 2.0 - (4.0/PI) * np.arctan(np.exp(x - x0))
def proton(x, center, w=KW):
    return kink(x, center - w/2) + antikink(x, center + w/2) - 2.0
def two_protons(x, R, w=KW):
    c = len(x)//2
    return proton(x, c - R/2, w) + proton(x, c + R/2, w)


# ============================================================
# DRESSED HESSIAN
# ============================================================
def build_dressed_hessian(phi_bg, delta_phi_sq):
    """H = -Laplacian + cos(pi*phi_bg) * exp(-pi^2 <delta_phi^2>/2)
    The dressing factor is <cos(pi*delta_phi)> for a Gaussian fluctuation
    with variance <delta_phi^2>.  Bounded in [0,1], equal to 1 when
    no occupation."""
    N = len(phi_bg)
    dress = np.exp(-PI**2 * delta_phi_sq / 2.0)
    diag = 2.0 + np.cos(PI * phi_bg) * dress
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1],
                     shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()


def low_eigs(phi_bg, delta_phi_sq, k=8):
    H = build_dressed_hessian(phi_bg, delta_phi_sq)
    w2, psi = eigsh(H, k=k, which='SA')   # SA = smallest algebraic
    idx = np.argsort(w2)
    return w2[idx], psi[:, idx]


# ============================================================
# SCF
# ============================================================
def occupation_dphi_sq(w2, psi, occupations):
    """<delta_phi^2>(x) = Sum_n (N_n + 1/2) * psi_n^2(x) / (2 omega_n)
    summed over bound modes (w2 > 0) with prescribed occupations."""
    dphi_sq = np.zeros(psi.shape[0])
    for mode_idx, N_e in occupations.items():
        if mode_idx >= len(w2):
            continue
        w2_n = w2[mode_idx]
        if w2_n <= 1e-6:
            continue  # skip tachyonic / zero modes
        omega = np.sqrt(w2_n)
        dphi_sq += (N_e + 0.5) * psi[:, mode_idx]**2 / (2.0 * omega)
    return dphi_sq


def scf(phi_bg, occupations, max_iter=30, mix=0.5, tol=1e-8, k=8, verbose=False):
    """Returns w2, psi at self-consistency."""
    N = len(phi_bg)
    dphi_sq = np.zeros(N)
    w2_prev = None
    for it in range(max_iter):
        w2, psi = low_eigs(phi_bg, dphi_sq, k=k)
        dphi_sq_new = occupation_dphi_sq(w2, psi, occupations)
        dphi_sq = mix * dphi_sq_new + (1 - mix) * dphi_sq

        # Convergence check: lowest eigenvalue stable
        if w2_prev is not None:
            change = abs(w2[0] - w2_prev[0])
            if verbose:
                print(f"    SCF iter {it:2d}: w2[0]={w2[0]:+.6f}  "
                      f"w2[1]={w2[1]:+.6f}  max<dphi^2>={dphi_sq.max():.4f}  "
                      f"d={change:.2e}")
            if change < tol:
                return {'w2': w2, 'psi': psi, 'dphi_sq': dphi_sq,
                        'iterations': it + 1, 'converged': True}
        w2_prev = w2.copy()
    return {'w2': w2, 'psi': psi, 'dphi_sq': dphi_sq,
            'iterations': max_iter, 'converged': False}


# ============================================================
# BOND ENERGY
# ============================================================
def bond_energy(N_lat, R, occ_dimer, occ_atom, verbose=False):
    """Run SCF for isolated atom and dimer at R, return D_e_lattice."""
    x = np.arange(N_lat, dtype=np.float64)

    # Isolated atom: occupy whatever occ_atom says (per-atom electrons)
    if verbose:
        print("  Isolated atom SCF:")
    iso = scf(proton(x, N_lat/2), occ_atom, mix=0.5, verbose=verbose)
    E_iso = iso['w2'][0]

    # Dimer: two atoms, doubled occupations
    if verbose:
        print(f"  Dimer SCF at R={R}:")
    dim = scf(two_protons(x, R), occ_dimer, mix=0.5, verbose=verbose)
    E_R = dim['w2'][0]

    return {
        'R': R,
        'E_iso': E_iso,
        'E_R': E_R,
        'V_R': E_R - E_iso,
        'D_e_lattice': -(E_R - E_iso),
        'iso_w2': iso['w2'],
        'dim_w2': dim['w2'],
        'iso_converged': iso['converged'],
        'dim_converged': dim['converged'],
        'iso_iter': iso['iterations'],
        'dim_iter': dim['iterations'],
    }


def scan_R(N_lat, R_values, occ_dimer, occ_atom):
    results = []
    for R in R_values:
        r = bond_energy(N_lat, R, occ_dimer, occ_atom)
        results.append(r)
    return results


# ============================================================
# RUN: H2, Cl2, F2 acid tests
# ============================================================
if __name__ == "__main__":
    N_lat = 256
    R_values = list(range(6, 13))

    print("=" * 78)
    print("COLLECTIVE FIELD SOLVER v2 — multi-mode SCF (QM mean-field backreaction)")
    print("=" * 78)
    print("Method: <delta_phi^2> dresses V'', not phi_bg.")
    print(f"Lattice N={N_lat}, R scan = {list(R_values)}")
    print()

    conv = lambda E_harm: E_harm / (2 * S_PT)

    # ------------------------------------------------------------
    # H2 — control: 2 electrons in sigma (mode 0), nothing else
    # ------------------------------------------------------------
    print("---- H2 (control: 1 sigma bond, no LP) ----")
    occ_atom_H = {0: 1}    # 1 electron per H atom
    occ_dimer_H = {0: 2}   # 2 electrons in sigma bond
    E_H = 13.598
    E_harm_H2 = E_H
    print(f"  Atom occupation:  {occ_atom_H}")
    print(f"  Dimer occupation: {occ_dimer_H}")
    print(f"  E_harm = {E_harm_H2:.3f} eV, conversion = {conv(E_harm_H2):.3f} eV/lat")
    res_H2 = scan_R(N_lat, R_values, occ_dimer_H, occ_atom_H)
    print(f"  {'R':>3} {'E_iso':>10} {'E_R':>10} {'V_R':>10} {'D_e_lat':>10} "
          f"{'D_e_eV':>10} {'iso_it':>7} {'dim_it':>7}")
    De_max_H2 = -1e9; R_eq_H2 = -1
    for r in res_H2:
        De_eV = r['D_e_lattice'] * conv(E_harm_H2)
        if r['D_e_lattice'] > De_max_H2:
            De_max_H2 = r['D_e_lattice']
            R_eq_H2 = r['R']
        print(f"  {r['R']:>3} {r['E_iso']:>+10.5f} {r['E_R']:>+10.5f} "
              f"{r['V_R']:>+10.5f} {r['D_e_lattice']:>+10.5f} "
              f"{De_eV:>+10.4f} {r['iso_iter']:>7d} {r['dim_iter']:>7d}")
    print(f"  >> H2 v2: D_e_eV = {De_max_H2 * conv(E_harm_H2):.3f}  "
          f"(bare 4.79, V8 4.48, obs 4.478)")
    print()

    # ------------------------------------------------------------
    # Cl2 — ACID TEST: 5 valence p-electrons per Cl, of which
    #       1 goes into sigma bond, 4 remain as LPs (in mode-1/2)
    # ------------------------------------------------------------
    print("---- Cl2 (4 LPs per atom; LP acid test) ----")
    # Per atom: 1 electron in mode 0 (sigma), 4 in higher modes (LPs)
    occ_atom_Cl = {0: 1, 1: 2, 2: 2}
    # Dimer: 2 electrons in sigma, 4 each in LP modes (4 from each atom)
    occ_dimer_Cl = {0: 2, 1: 4, 2: 4}
    E_Cl = 12.968
    E_harm_Cl2 = E_Cl
    print(f"  Atom occupation:  {occ_atom_Cl}")
    print(f"  Dimer occupation: {occ_dimer_Cl}")
    print(f"  E_harm = {E_harm_Cl2:.3f} eV, conversion = {conv(E_harm_Cl2):.3f} eV/lat")
    res_Cl2 = scan_R(N_lat, R_values, occ_dimer_Cl, occ_atom_Cl)
    print(f"  {'R':>3} {'E_iso':>10} {'E_R':>10} {'V_R':>10} {'D_e_lat':>10} "
          f"{'D_e_eV':>10} {'iso_it':>7} {'dim_it':>7}")
    De_max_Cl2 = -1e9; R_eq_Cl2 = -1
    for r in res_Cl2:
        De_eV = r['D_e_lattice'] * conv(E_harm_Cl2)
        if r['D_e_lattice'] > De_max_Cl2:
            De_max_Cl2 = r['D_e_lattice']
            R_eq_Cl2 = r['R']
        print(f"  {r['R']:>3} {r['E_iso']:>+10.5f} {r['E_R']:>+10.5f} "
              f"{r['V_R']:>+10.5f} {r['D_e_lattice']:>+10.5f} "
              f"{De_eV:>+10.4f} {r['iso_iter']:>7d} {r['dim_iter']:>7d}")
    print(f"  >> Cl2 v2: D_e_eV = {De_max_Cl2 * conv(E_harm_Cl2):.3f}  "
          f"(bare 4.51, V8 3.00, obs 2.514)")
    print()

    # ------------------------------------------------------------
    # F2 — extreme LP test (V8 already gets this within 4%; bare wildly off)
    # ------------------------------------------------------------
    print("---- F2 (most extreme LP case) ----")
    occ_atom_F = {0: 1, 1: 2, 2: 2}
    occ_dimer_F = {0: 2, 1: 4, 2: 4}
    E_F = 17.423
    E_harm_F2 = E_F
    print(f"  Atom occupation:  {occ_atom_F}")
    print(f"  Dimer occupation: {occ_dimer_F}")
    print(f"  E_harm = {E_harm_F2:.3f} eV, conversion = {conv(E_harm_F2):.3f} eV/lat")
    res_F2 = scan_R(N_lat, R_values, occ_dimer_F, occ_atom_F)
    print(f"  {'R':>3} {'E_iso':>10} {'E_R':>10} {'V_R':>10} {'D_e_lat':>10} "
          f"{'D_e_eV':>10} {'iso_it':>7} {'dim_it':>7}")
    De_max_F2 = -1e9; R_eq_F2 = -1
    for r in res_F2:
        De_eV = r['D_e_lattice'] * conv(E_harm_F2)
        if r['D_e_lattice'] > De_max_F2:
            De_max_F2 = r['D_e_lattice']
            R_eq_F2 = r['R']
        print(f"  {r['R']:>3} {r['E_iso']:>+10.5f} {r['E_R']:>+10.5f} "
              f"{r['V_R']:>+10.5f} {r['D_e_lattice']:>+10.5f} "
              f"{De_eV:>+10.4f} {r['iso_iter']:>7d} {r['dim_iter']:>7d}")
    print(f"  >> F2 v2: D_e_eV = {De_max_F2 * conv(E_harm_F2):.3f}  "
          f"(bare 13.72, V8 1.54, obs 1.602)")
    print()

    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print(f"{'Mol':>5} {'bare(eV)':>10} {'v2(eV)':>10} {'V8(eV)':>10} "
          f"{'obs(eV)':>10} {'v2_err':>9}")
    for label, bare, v2_val, V8_val, obs in [
            ('H2',  4.789, De_max_H2 * conv(E_harm_H2),  4.481, 4.478),
            ('Cl2', 4.512, De_max_Cl2 * conv(E_harm_Cl2), 3.001, 2.514),
            ('F2', 13.723, De_max_F2 * conv(E_harm_F2),  1.542, 1.602)]:
        err = (v2_val - obs)/obs * 100
        print(f"{label:>5} {bare:>10.3f} {v2_val:>10.3f} {V8_val:>10.3f} "
              f"{obs:>10.3f} {err:>+8.1f}%")

    print()
    print("Hypothesis: LP backreaction shrinks D_e for Cl2/F2 toward observed.")
    print("Pass: Cl2 v2 < 4 eV, F2 v2 < 6 eV.")
    print("Fail: v2 ~ bare, meaning multi-mode occupation isn't the LP knob.")
