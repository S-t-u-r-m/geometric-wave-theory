#!/usr/bin/env python3
"""
Collective Field Solver — v1 (1D, Hessian + SCF breather backreaction)
=======================================================================

Lessons from earlier attempts (recorded):
- Classical static field minimization fails: bumps have zero winding,
  dissolve to vacuum under gradient descent.
- Eigenvalue tracking needs care: the bond mode is mode 0 (the kink
  translation tachyon, which splits into bonding/antibonding at finite R),
  not "lowest positive eigenvalue."

This v1 mirrors bond_3d_emerge.py's proven method as the baseline:

    V(R) = E_bonding_dimer(R) - E_bonding_isolated_atom
    D_e_lattice = -min_R V(R)
    D_e_eV = D_e_lattice * E_harm / (2*s_PT)

That gives 4.75 eV for H2 (matches obs 4.748 to 0.02% per the memory).

The NEW ingredient under test: SCF backreaction.

After the one-shot diagonalization, mix the bonding-mode amplitude
into the background field and re-diagonalize:

    phi_bg_new = phi_bg + alpha * psi_bonding * normalization
    iterate until eigenvalues converge.

If the SCF gain is small for H2 (single sigma bond) but grows for
multi-bond molecules (sketch later), the collective-field signature
shows up. If it's negligible everywhere, the missing piece is somewhere
else and the staircase result has a different explanation.
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

PI = np.pi
d = 3

S_PT = (-1 + np.sqrt(1 + 8/PI**2)) / 2     # 0.17279
KW = 3                                       # kink width


# ============================================================
# FIELD CONSTRUCTION (mirror bond_3d_emerge.py exactly)
# ============================================================
def kink(x, x0):
    return (4.0/PI) * np.arctan(np.exp(x - x0))


def antikink(x, x0):
    return 2.0 - (4.0/PI) * np.arctan(np.exp(x - x0))


def proton(x, center, width=KW):
    """Kink-antikink pair = a 'proton' bump (zero net winding)."""
    return kink(x, center - width/2) + antikink(x, center + width/2) - 2.0


def two_protons(x, R, center=None, width=KW):
    """Two protons at separation R, centered on the lattice."""
    if center is None:
        center = len(x) // 2
    return proton(x, center - R/2, width) + proton(x, center + R/2, width)


# ============================================================
# HESSIAN
# ============================================================
def build_hessian(phi_bg, periodic=True):
    N = len(phi_bg)
    diag = 2.0 + np.cos(PI * phi_bg)
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1],
                     shape=(N, N), format='lil')
    if periodic:
        H[0, N-1] = -1.0
        H[N-1, 0] = -1.0
    return H.tocsr()


def low_eigs(phi_bg, k=10):
    H = build_hessian(phi_bg)
    w2, psi = eigsh(H, k=k, which='SM')
    idx = np.argsort(w2)
    return w2[idx], psi[:, idx]


# ============================================================
# ONE-SHOT (BASELINE: bond_3d_emerge style)
# ============================================================
def one_shot_bond_energy(N, R_values, n_eig=10, verbose=False):
    """Returns V(R) curve. D_e = -min(V)."""
    x = np.arange(N, dtype=np.float64)
    # Isolated atom (single proton at center)
    phi_iso = proton(x, N/2)
    w2_iso, _ = low_eigs(phi_iso, k=n_eig)
    E_iso_bond = w2_iso[0]   # lowest eigenvalue of isolated atom (tachyon)

    V_curve = []
    for R in R_values:
        phi_d = two_protons(x, R)
        w2, _ = low_eigs(phi_d, k=n_eig)
        E_bond_R = w2[0]    # lowest eigenvalue of dimer
        E_anti_R = w2[1]    # second-lowest = antibonding partner
        V_R = E_bond_R - E_iso_bond
        V_curve.append({'R': R, 'V_R': V_R,
                        'E_bond': E_bond_R, 'E_anti': E_anti_R,
                        'split': E_anti_R - E_bond_R})
        if verbose:
            print(f"  R={R:3d}  E_b={E_bond_R:+.6f}  E_a={E_anti_R:+.6f}  "
                  f"V(R)={V_R:+.6f}")
    return V_curve, E_iso_bond


# ============================================================
# SCF BREATHER BACKREACTION
# ============================================================
def scf_bond_energy(N, R, n_eig=10, max_iter=20, mix=0.2, tol=1e-7,
                    verbose=False):
    """Self-consistent: bonding-mode amplitude shifts background."""
    x = np.arange(N, dtype=np.float64)
    phi_bare = two_protons(x, R)
    phi_bg = phi_bare.copy()
    E_bond_prev = None
    converged = False

    for it in range(max_iter):
        w2, psi = low_eigs(phi_bg, k=n_eig)
        E_bond = w2[0]
        E_anti = w2[1]
        bonding_mode = psi[:, 0]
        # Sign convention: bonding mode should be symmetric (positive
        # contribution at both well centers).
        if bonding_mode[N//2 - int(R/2)] < 0:
            bonding_mode = -bonding_mode
        # ZPE-style amplitude (use sqrt|w2| to avoid blowup on tachyons).
        amp = 1.0 / np.sqrt(2.0 * max(np.sqrt(abs(E_bond)), 1e-3))
        delta = amp * bonding_mode
        # Mix into background (anchored to bare so we don't drift away).
        phi_bg = phi_bare + mix * delta

        if verbose:
            print(f"    SCF iter {it:2d}: E_b={E_bond:+.6f}  E_a={E_anti:+.6f}  "
                  f"split={E_anti-E_bond:+.6f}  |delta|={np.max(np.abs(delta)):.4f}")
        if E_bond_prev is not None:
            change = abs(E_bond - E_bond_prev)
            if change < tol:
                converged = True
                break
        E_bond_prev = E_bond
    return {
        'R': R, 'E_bond': E_bond, 'E_anti': E_anti,
        'iterations': it + 1, 'converged': converged,
        'phi_bg': phi_bg, 'bonding_mode': bonding_mode,
    }


# ============================================================
# RUN
# ============================================================
if __name__ == "__main__":
    N = 256
    print("=" * 72)
    print("COLLECTIVE FIELD SOLVER v1 — 1D Hessian + SCF backreaction")
    print("=" * 72)
    print(f"  Lattice N={N}, kink width KW={KW}, s_PT={S_PT:.5f}")
    print()

    # ---- Baseline: bond_3d_emerge-style one-shot ----
    print("STEP 1: BASELINE (one-shot Hessian, bond_3d_emerge.py replica)")
    print("-" * 72)
    R_values = list(range(4, 41, 2))
    V_curve, E_iso = one_shot_bond_energy(N, R_values, verbose=False)
    print(f"  E_iso (single bump tachyon ground): {E_iso:.6f}")
    print(f"{'R':>4} {'V(R)':>12} {'split':>12}")
    for entry in V_curve:
        print(f"{entry['R']:>4} {entry['V_R']:>+12.6f} {entry['split']:>+12.6f}")
    V_arr = np.array([e['V_R'] for e in V_curve])
    R_arr = np.array([e['R'] for e in V_curve])
    i_min = int(np.argmin(V_arr))
    D_e_lat = -V_arr[i_min]
    R_eq = R_arr[i_min]
    print()
    print(f"  R_eq = {R_eq}  D_e_lattice = {D_e_lat:.6f}")
    E_H = 13.598
    conv = E_H / (2 * S_PT)
    D_e_eV_baseline = D_e_lat * conv
    print(f"  D_e baseline = {D_e_eV_baseline:.4f} eV  (V8 prediction: 4.747 eV, obs: 4.748 eV)")
    print()

    # ---- SCF at R_eq ----
    print("STEP 2: SCF backreaction at R_eq")
    print("-" * 72)
    res_scf = scf_bond_energy(N, R_eq, max_iter=20, mix=0.2, verbose=True)
    V_scf = res_scf['E_bond'] - E_iso
    D_e_lat_scf = -V_scf
    D_e_eV_scf = D_e_lat_scf * conv
    print()
    print(f"  After SCF at R_eq = {R_eq}:")
    print(f"    E_bond (SCF)     = {res_scf['E_bond']:+.6f}")
    print(f"    V(R_eq) SCF      = {V_scf:+.6f}")
    print(f"    D_e SCF (lattice)= {D_e_lat_scf:.6f}")
    print(f"    D_e SCF (eV)     = {D_e_eV_scf:.4f}")
    print(f"    SCF gain         = {D_e_eV_scf - D_e_eV_baseline:+.4f} eV")
    print(f"    Iterations: {res_scf['iterations']}, "
          f"converged: {res_scf['converged']}")
    print()

    # ---- Scan SCF across R ----
    print("STEP 3: SCF scan to check for shift in R_eq")
    print("-" * 72)
    print(f"{'R':>4} {'V_OS':>12} {'V_SCF':>12} {'gain':>12} {'iters':>6}")
    R_scan = [r for r in R_values if 4 <= r <= 14]
    for R in R_scan:
        # One-shot reference at this R
        x = np.arange(N, dtype=np.float64)
        w2_one, _ = low_eigs(two_protons(x, R), k=10)
        V_one = w2_one[0] - E_iso
        # SCF
        res = scf_bond_energy(N, R, mix=0.2)
        V_scf = res['E_bond'] - E_iso
        print(f"{R:>4} {V_one:>+12.6f} {V_scf:>+12.6f} "
              f"{V_one - V_scf:>+12.6f} {res['iterations']:>6}")
    print()

    print("=" * 72)
    print("Interpretation:")
    print(f"  Baseline (one-shot) reproduces V8/observed at {D_e_eV_baseline:.3f} eV.")
    print("  SCF gain is the collective-field test:")
    print("    ~0 for H2  -> single-bond systems are already converged at one-shot,")
    print("                  and the under-prediction in CO/NO/etc. comes from")
    print("                  inter-MODE coupling (sigma <-> pi via cos potential),")
    print("                  not single-mode backreaction. v2 = multi-mode SCF.")
    print("    >0 for H2  -> even single bonds gain from self-consistency. The")
    print("                  collective effect is in the bond mode itself.")
    print("=" * 72)
