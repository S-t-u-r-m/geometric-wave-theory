#!/usr/bin/env python3
"""
Collective Field Solver — v3 (3D, sparse Hessian, optional SCF)
================================================================

Framing reminder (Jon, 2026-05-22):
GWT's "1D" refers ONLY to the scalar (single-component) character of phi.
Waves IN phi are fully 3D spatial objects on a 3D lattice. v1/v2 were
transverse-uniform projections of a 3D wave; v3 keeps the full 3D wave.

Architecture
------------
- 3D cubic lattice (N^3 sites, periodic BC).
- Spherical "atoms" = 3D radial bumps in phi.
- Sparse 3D Hessian H = -Laplacian + cos(pi*phi_bg) [diagonal].
- scipy.sparse.linalg.eigsh for lowest k eigenvalues+vectors.
- Bond energy = lowest-eigenvalue shift from isolated atom, converted via
  V8's E_harm/(2*s_PT) map (same convention as bond_3d_emerge.py).
- Mode classification: angular decomposition of each eigenmode at each
  atom center, projected onto Oh-irrep angular functions:
    A1g (sigma):   1
    T1u (p/LP):    x/r, y/r, z/r
    Eg  (d_z2,..): (3z^2 - r^2), (x^2 - y^2)
    T2g (d_xy,..): xy, yz, zx
  This tells us which eigenmodes carry transverse content where LP physics
  lives.

This pass: baseline only (no SCF). Two purposes:
  1. Validate H2 in 3D (should reproduce ~4.7 eV like 1D, since H2 sigma
     bond has minimal transverse structure).
  2. Check whether Cl2's BARE 3D baseline differs significantly from 1D
     bare — if yes, the 3D mode structure already captures something
     the 1D projection misses, and that's the SCF starting point. If no,
     SCF on top is the only knob and we need to be careful about the
     mean-field breakdown that killed v2 on Cl2.
"""
import sys
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

PI = np.pi
d = 3
S_PT = (-1 + np.sqrt(1 + 8/PI**2)) / 2
KW = 3   # kink width (lattice sites)


# ============================================================
# 3D LATTICE + SPHERICAL BUMP
# ============================================================
def make_grid_3d(N):
    """Centered 3D coords, integer lattice spacings."""
    x1d = np.arange(N, dtype=np.float64) - N/2
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing='ij')
    return X, Y, Z


def spherical_bump(X, Y, Z, center_xyz, w=KW):
    """3D radial extension of the 1D bump.
    Field peaks (~1.44) at center, decays to 0 at large r.
    Construction: phi(r) = (4/pi) * (arctan(exp(-r + w/2)) - arctan(exp(-r - w/2)))
    where r = distance from atom center."""
    cx, cy, cz = center_xyz
    r = np.sqrt((X - cx)**2 + (Y - cy)**2 + (Z - cz)**2)
    return (4.0/PI) * (np.arctan(np.exp(-r + w/2)) - np.arctan(np.exp(-r - w/2)))


def two_atoms_3d(X, Y, Z, R, w=KW):
    """Two spherical atoms separated by R along z-axis, centered at origin."""
    a = spherical_bump(X, Y, Z, (0, 0, -R/2), w)
    b = spherical_bump(X, Y, Z, (0, 0, +R/2), w)
    return a + b


# ============================================================
# 3D SPARSE HESSIAN (periodic BC)
# ============================================================
def build_hessian_3d(phi_3d, dress_3d=None):
    """Sparse 3D Hessian for sine-Gordon fluctuations.
    H = -Laplacian + cos(pi*phi) * dress
    Periodic BC. dress=None means dress=1 (no SCF backreaction)."""
    N = phi_3d.shape[0]
    N_tot = N**3
    phi_flat = phi_3d.flatten()

    # Diagonal: 6 (from -Laplacian on cubic) + cos(pi*phi)*dress
    if dress_3d is None:
        diag = 6.0 + np.cos(PI * phi_flat)
    else:
        diag = 6.0 + np.cos(PI * phi_flat) * dress_3d.flatten()

    # Build sparse matrix in COO format, then convert to CSR.
    # Indices: i = x*N*N + y*N + z
    idx_3d = np.arange(N_tot).reshape(N, N, N)
    rows_list = [np.arange(N_tot)]
    cols_list = [np.arange(N_tot)]
    vals_list = [diag]

    for axis in range(3):
        for shift in (-1, +1):
            neighbor = np.roll(idx_3d, shift, axis=axis).flatten()
            rows_list.append(np.arange(N_tot))
            cols_list.append(neighbor)
            vals_list.append(-np.ones(N_tot))

    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    vals = np.concatenate(vals_list)

    return sparse.coo_matrix((vals, (rows, cols)),
                             shape=(N_tot, N_tot)).tocsr()


def low_eigs_3d(phi_3d, dress_3d=None, k=8):
    """Lowest k eigenvalues + eigenvectors (algebraic order)."""
    H = build_hessian_3d(phi_3d, dress_3d)
    w2, psi = eigsh(H, k=k, which='SA')
    idx = np.argsort(w2)
    return w2[idx], psi[:, idx]


# ============================================================
# MODE ANGULAR DECOMPOSITION (Oh-irrep projection)
# ============================================================
def mode_angular_content(psi_flat, X, Y, Z, atom_center, r_max=KW+2):
    """Decompose psi(r) near an atom into Oh-irrep angular components.

    Restrict to a sphere of radius r_max around atom_center.
    Project onto angular basis functions:
      A1g (s, sigma):  1
      T1u (p):         x/r, y/r, z/r
      Eg  (d, eg):     (3z^2 - r^2)/r^2, (x^2 - y^2)/r^2
      T2g (d, t2g):    xy/r^2, yz/r^2, zx/r^2

    Returns dict of {irrep: total squared projection (= 'amount of this character')}.
    """
    cx, cy, cz = atom_center
    N = X.shape[0]
    psi = psi_flat.reshape(N, N, N)
    DX = X - cx
    DY = Y - cy
    DZ = Z - cz
    r = np.sqrt(DX**2 + DY**2 + DZ**2)
    mask = (r > 1e-3) & (r <= r_max)

    psi_loc = psi[mask]
    r_loc = r[mask]
    nx, ny, nz = DX[mask]/r_loc, DY[mask]/r_loc, DZ[mask]/r_loc

    # Build angular basis (real spherical harmonics)
    basis = {
        'A1g_s':  np.ones_like(r_loc),
        'T1u_x':  nx,
        'T1u_y':  ny,
        'T1u_z':  nz,
        'Eg_z2':  (3*nz**2 - 1),
        'Eg_x2y2':(nx**2 - ny**2),
        'T2g_xy': nx*ny,
        'T2g_yz': ny*nz,
        'T2g_zx': nz*nx,
    }
    # Project (without weighting by r^2 jacobian — we just want angular
    # character of the wavefunction near the atom).
    irrep_groups = {
        'A1g': ['A1g_s'],
        'T1u': ['T1u_x', 'T1u_y', 'T1u_z'],
        'Eg':  ['Eg_z2', 'Eg_x2y2'],
        'T2g': ['T2g_xy', 'T2g_yz', 'T2g_zx'],
    }
    out = {}
    for irrep, keys in irrep_groups.items():
        total = 0.0
        for k in keys:
            b = basis[k]
            # Normalize each basis function within the sphere
            b_norm2 = np.sum(b**2)
            if b_norm2 < 1e-12:
                continue
            coeff = np.sum(psi_loc * b) / np.sqrt(b_norm2)
            total += coeff**2
        out[irrep] = total
    # Total squared norm (for fractional report)
    total_norm2 = np.sum(psi_loc**2)
    if total_norm2 > 0:
        for k in out:
            out[k] /= total_norm2
    return out


# ============================================================
# BOND ENERGY (BARE BASELINE)
# ============================================================
def bond_energy_bare(N, R, k_eig=4, verbose=False):
    X, Y, Z = make_grid_3d(N)
    # Isolated atom (single bump at origin)
    phi_iso = spherical_bump(X, Y, Z, (0, 0, 0))
    t0 = time.time()
    w2_iso, psi_iso = low_eigs_3d(phi_iso, k=k_eig)
    t_iso = time.time() - t0

    # Dimer
    phi_dim = two_atoms_3d(X, Y, Z, R)
    t0 = time.time()
    w2_dim, psi_dim = low_eigs_3d(phi_dim, k=k_eig)
    t_dim = time.time() - t0

    E_iso = w2_iso[0]
    E_R = w2_dim[0]
    if verbose:
        print(f"    iso w2: {[f'{w:+.4f}' for w in w2_iso[:4]]}  ({t_iso:.1f}s)")
        print(f"    dim w2: {[f'{w:+.4f}' for w in w2_dim[:4]]}  ({t_dim:.1f}s)")

    return {
        'R': R, 'E_iso': E_iso, 'E_R': E_R,
        'V_R': E_R - E_iso,
        'D_e_lattice': -(E_R - E_iso),
        'w2_iso': w2_iso, 'w2_dim': w2_dim,
        'psi_iso': psi_iso, 'psi_dim': psi_dim,
        'X': X, 'Y': Y, 'Z': Z,
    }


# ============================================================
# RUN
# ============================================================
if __name__ == "__main__":
    N = 24    # 24^3 = 13824 sites (lightweight, ~seconds per eigsh)
    print("=" * 78)
    print("COLLECTIVE FIELD SOLVER v3 — 3D bare baseline")
    print("=" * 78)
    print(f"Lattice: N={N}, N^3={N**3:,}, periodic BC")
    print(f"Method: spherical 3D bumps, sparse Hessian, scipy eigsh (CPU)")
    print()

    conv = lambda E_harm: E_harm / (2 * S_PT)

    # ---- H2 baseline scan ----
    print("---- H2: 3D bare baseline scan ----")
    E_harm_H2 = 13.598
    print(f"  E_harm = {E_harm_H2:.3f} eV, conv = {conv(E_harm_H2):.3f} eV/lat")
    R_scan = [4, 5, 6, 7, 8]
    print(f"  {'R':>3} {'E_iso':>10} {'E_R':>10} {'V_R':>10} "
          f"{'D_e_lat':>10} {'D_e_eV':>10}")
    De_max = -1e9; R_eq = -1
    for R in R_scan:
        t0 = time.time()
        r = bond_energy_bare(N, R, k_eig=4, verbose=False)
        elapsed = time.time() - t0
        De_eV = r['D_e_lattice'] * conv(E_harm_H2)
        if r['D_e_lattice'] > De_max:
            De_max = r['D_e_lattice']
            R_eq = R
        print(f"  {R:>3} {r['E_iso']:>+10.5f} {r['E_R']:>+10.5f} "
              f"{r['V_R']:>+10.5f} {r['D_e_lattice']:>+10.5f} "
              f"{De_eV:>+10.4f}  ({elapsed:.1f}s)")
    print(f"  >> H2 3D bare: D_e = {De_max * conv(E_harm_H2):.3f} eV at R_eq={R_eq}")
    print(f"     1D bare: 4.789 eV; V8: 4.481; observed: 4.478")
    print()

    # ---- Compare angular content of low modes at R_eq ----
    print("---- Mode angular content at H2 R_eq ----")
    r = bond_energy_bare(N, R_eq, k_eig=6)
    print(f"  Eigenvalues at R={R_eq}: {[f'{w:+.4f}' for w in r['w2_dim']]}")
    print(f"  Angular content of each mode (around atom A at z=-R/2):")
    print(f"  {'mode':>5} {'w^2':>10}  {'A1g':>6} {'T1u':>6} {'Eg':>6} {'T2g':>6}")
    for i in range(min(6, len(r['w2_dim']))):
        ang = mode_angular_content(r['psi_dim'][:, i], r['X'], r['Y'], r['Z'],
                                    (0, 0, -R_eq/2))
        print(f"  {i:>5} {r['w2_dim'][i]:>+10.4f}  "
              f"{ang.get('A1g',0):>6.2f} {ang.get('T1u',0):>6.2f} "
              f"{ang.get('Eg',0):>6.2f} {ang.get('T2g',0):>6.2f}")
    print()

    # ============================================================
    # RATIO TEST: H2 vs Cl2 vs F2 in 3D bare
    # ============================================================
    # The LP hypothesis predicts that ratio D_e(Cl2)/D_e(H2) changes
    # significantly between 1D and 3D. In 1D bare ratio is 0.94 (basically
    # unchanged from E_harm scaling). Observed ratio is 0.53. If 3D bare
    # ratio is closer to observed (without any LP-correction term), the
    # 3D mode structure is doing real work.
    print("=" * 78)
    print("RATIO TEST: 3D bare D_e for H2 vs Cl2 vs F2")
    print("=" * 78)
    R_scan_ratio = [4, 5, 6, 7]   # short scan, take the max
    results_mol = {}
    for label, E_harm in [('H2', 13.598), ('Cl2', 12.968), ('F2', 17.423)]:
        De_max = -1e9
        R_best = -1
        for R in R_scan_ratio:
            r = bond_energy_bare(N, R, k_eig=4)
            if r['D_e_lattice'] > De_max:
                De_max = r['D_e_lattice']
                R_best = R
        De_eV = De_max * conv(E_harm)
        results_mol[label] = {'De_lat': De_max, 'De_eV': De_eV, 'R': R_best}
        print(f"  {label:>4}: D_e_lattice = {De_max:.4f} at R={R_best}, "
              f"E_harm={E_harm:.2f} -> D_e = {De_eV:.3f} eV")

    print()
    # The D_e_lattice value is Z-INDEPENDENT in this bare model — the bump
    # is the same shape regardless of "atom identity." So all molecules
    # give the same D_e_lattice; ratios are determined by E_harm alone.
    # That means 3D bare also says ratio = E_harm(Cl2)/E_harm(H2) = 0.95.
    # The 3D angular structure changes WHICH MODES are bound, but doesn't
    # touch the per-molecule scaling unless we add per-atom physics.
    print("Bare-baseline ratios (should match 1D bare if just E_harm scaling):")
    De_H2 = results_mol['H2']['De_eV']
    De_Cl2 = results_mol['Cl2']['De_eV']
    De_F2 = results_mol['F2']['De_eV']
    print(f"  Cl2/H2: 3D bare = {De_Cl2/De_H2:.3f}  |  observed = "
          f"{2.514/4.478:.3f}  |  1D bare = {4.512/4.789:.3f}")
    print(f"  F2/H2:  3D bare = {De_F2/De_H2:.3f}  |  observed = "
          f"{1.602/4.478:.3f}  |  1D bare = {13.723/4.789:.3f}")
    print()
    print("Reading: if 3D bare ratio = 1D bare ratio, the 3D structure doesn't")
    print("change anything at the bare level (because we have no per-atom")
    print("specificity in the bump). That's expected — LP physics requires")
    print("OCCUPATION of the 3D transverse modes (in SCF), not just their")
    print("existence in the Hessian spectrum. Next: multi-mode SCF.")
    print()

    # ============================================================
    # MULTI-MODE SCF — the actual LP test
    # ============================================================
    print("=" * 78)
    print("MULTI-MODE SCF in 3D — does occupying T1u (LP) modes shrink D_e?")
    print("=" * 78)
    print()
    print("Per-atom electron assignment (V8 convention):")
    print("  H:  1 electron  -> 0.5 in mode 0 (sigma)")
    print("  Cl: 7 valence e -> 1 in sigma, 4 in LPs, 2 in deep s (ignore)")
    print("  F:  7 valence e -> same as Cl")
    print()
    print("Dimer (occ_dimer):")
    print("  H2:  mode 0 = 2e (sigma only)")
    print("  Cl2: mode 0 = 2e (sigma), modes 2+3 = 2e each (T1u LP) "
          "[8 LP electrons want more modes but only 4 fit here]")

    def occupation_dphi_sq_3d(w2, psi, occupations):
        """<delta_phi^2>(x) on the 3D flat array."""
        N_tot = psi.shape[0]
        dphi_sq = np.zeros(N_tot)
        for idx, N_e in occupations.items():
            if idx >= len(w2):
                continue
            if w2[idx] <= 1e-6:
                continue
            omega = np.sqrt(w2[idx])
            dphi_sq += (N_e + 0.5) * psi[:, idx]**2 / (2.0 * omega)
        return dphi_sq

    def scf_3d(phi_bg_3d, occupations, max_iter=10, mix=0.3, tol=1e-6, k=8):
        N = phi_bg_3d.shape[0]
        N_tot = N**3
        dphi_sq = np.zeros(N_tot)
        w2_prev = None
        for it in range(max_iter):
            dress_3d = np.exp(-PI**2 * dphi_sq / 2.0).reshape(N, N, N)
            w2, psi = low_eigs_3d(phi_bg_3d, dress_3d, k=k)
            dphi_sq_new = occupation_dphi_sq_3d(w2, psi, occupations)
            dphi_sq = mix * dphi_sq_new + (1 - mix) * dphi_sq
            if w2_prev is not None and abs(w2[0] - w2_prev[0]) < tol:
                return {'w2': w2, 'psi': psi, 'iterations': it+1,
                        'converged': True, 'dphi_sq': dphi_sq}
            w2_prev = w2.copy()
        return {'w2': w2, 'psi': psi, 'iterations': max_iter,
                'converged': False, 'dphi_sq': dphi_sq}

    # ---- Test: H2 (control), Cl2 (LP test), F2 (extreme LP test) ----
    X, Y, Z = make_grid_3d(N)
    R_test = 4   # use R_eq from bare scan

    for label, occ_atom, occ_dimer, E_harm, obs in [
        ('H2',
         {0: 1},
         {0: 2},
         13.598, 4.478),
        ('Cl2',
         {0: 1, 1: 1, 2: 1},                  # atom: sigma + 2 transverse
         {0: 2, 2: 2, 3: 2},                  # dimer: sigma + 2 T1u (4 LP e)
         12.968, 2.514),
        ('F2',
         {0: 1, 1: 1, 2: 1},
         {0: 2, 2: 2, 3: 2},
         17.423, 1.602),
    ]:
        print(f"---- {label} (R={R_test}) ----")
        print(f"   occ_atom={occ_atom}  occ_dimer={occ_dimer}")
        # Isolated atom
        phi_iso = spherical_bump(X, Y, Z, (0, 0, 0))
        t0 = time.time()
        iso = scf_3d(phi_iso, occ_atom, max_iter=10, mix=0.3, k=8)
        t_iso = time.time() - t0
        E_iso = iso['w2'][0]
        # Dimer
        phi_dim = two_atoms_3d(X, Y, Z, R_test)
        t0 = time.time()
        dim = scf_3d(phi_dim, occ_dimer, max_iter=10, mix=0.3, k=8)
        t_dim = time.time() - t0
        E_R = dim['w2'][0]
        V_R = E_R - E_iso
        De_lat = -V_R
        De_eV = De_lat * conv(E_harm)

        # Also compute BARE (no SCF) for comparison at same R
        bare = bond_energy_bare(N, R_test, k_eig=4)
        De_bare_eV = bare['D_e_lattice'] * conv(E_harm)

        print(f"   E_iso(bare)={bare['E_iso']:+.4f}  E_iso(SCF)={E_iso:+.4f}  "
              f"(iso {iso['iterations']}/{'C' if iso['converged'] else 'N'}, {t_iso:.1f}s)")
        print(f"   E_R(bare)={bare['E_R']:+.4f}    E_R(SCF)={E_R:+.4f}    "
              f"(dim {dim['iterations']}/{'C' if dim['converged'] else 'N'}, {t_dim:.1f}s)")
        print(f"   D_e bare = {De_bare_eV:.3f} eV   D_e SCF = {De_eV:.3f} eV   "
              f"observed = {obs:.3f} eV")
        print(f"   SCF shift = {De_eV - De_bare_eV:+.3f} eV  "
              f"(negative = shrinks toward LP-repulsed observed)")
        print()
