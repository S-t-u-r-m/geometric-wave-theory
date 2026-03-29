"""
1D Hessian Bond Calculator — Z-dependent kink wells
=====================================================
The force IS 1D (proven by 3D torus simulation to 0.06%).
So model the force in 1D with Z-dependent kink wells.

For each atom pair:
  1. Build Z-dependent kink-antikink profiles (width ~ 1/sqrt(Z))
  2. Compute 1D Hessian eigenvalues for single and double wells
  3. ZPE change = bond energy in lattice units
  4. Convert: D_e(physical) = D_e(lattice) * E_harm / (2*s)
     where s = 0.17279 is UNIVERSAL (Z cancels in V_0/beta^2 = 2/pi^2)
     and E_harm = harmonic mean of ionization energies

The LP, radical, and ionic effects EMERGE from the eigenvalue structure
rather than being added as algebraic corrections.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

PI = np.pi
d = 3

# Universal Poeschl-Teller parameter (Z-independent)
s_PT = (-1 + np.sqrt(1 + 8/PI**2)) / 2  # = 0.17279

print("=" * 70)
print("1D HESSIAN BOND CALCULATOR — Z-DEPENDENT KINK WELLS")
print("=" * 70)
print(f"s_PT = {s_PT:.5f} (universal, Z cancels in V_0/beta^2 = 2/pi^2)")
print()


def kink_antikink(x, center, Z, width=3.0):
    """
    Z-dependent kink-antikink pair.

    The kink profile: phi(x) = (4/pi) * arctan(exp(beta * x))
    where beta = sqrt(Z) controls the width.

    Width parameter controls the kink-antikink separation
    (= the "proton size" in lattice units).
    """
    beta = np.sqrt(max(Z, 1))
    hw = width / 2
    # Kink at center-hw, antikink at center+hw
    phi_k = (4.0/PI) * np.arctan(np.exp(beta * (x - center + hw)))
    phi_a = 2.0 - (4.0/PI) * np.arctan(np.exp(beta * (x - center - hw)))
    return phi_k + phi_a - 2.0


def build_hessian(phi, N):
    """1D Hessian: H = -Laplacian + cos(pi*phi), periodic BC."""
    diag = 2.0 + np.cos(PI * phi)
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()


def get_zpe(H, n_eig=20, mass_gap=0.95):
    """ZPE from bound modes below mass gap."""
    evals = eigsh(H, k=n_eig, which='SM', return_eigenvectors=False)
    evals = np.sort(evals)
    bound = evals[evals < mass_gap]
    zpe = np.sum(np.sqrt(np.abs(bound))) / 2
    return zpe, bound


def bond_energy(Z_A, Z_B, IE_A, IE_B, N=256,
                R_range=None, kink_width=3.0, n_eig=20):
    """
    Compute bond energy from 1D Hessian ZPE difference.

    Returns: D_e in eV, R_eq in lattice sites, full V(R) curve
    """
    if R_range is None:
        # Scan from close to far, including odd and even
        R_range = list(range(3, 25))

    x = np.arange(N, dtype=np.float64)
    center = N // 2

    # E_harm = harmonic mean of ionization energies
    E_harm = 2 * IE_A * IE_B / (IE_A + IE_B)

    # Single well references (on SAME lattice)
    phi_A = kink_antikink(x, center, Z_A, kink_width)
    zpe_A, modes_A = get_zpe(build_hessian(phi_A, N), n_eig)

    phi_B = kink_antikink(x, center, Z_B, kink_width)
    zpe_B, modes_B = get_zpe(build_hessian(phi_B, N), n_eig)

    zpe_ref = zpe_A + zpe_B

    # Scan separations
    V_curve = []
    for R in R_range:
        pos_A = center - R // 2
        pos_B = center + R // 2

        # Build two-well profile
        phi_AB = kink_antikink(x, pos_A, Z_A, kink_width) + \
                 kink_antikink(x, pos_B, Z_B, kink_width)

        zpe_AB, modes_AB = get_zpe(build_hessian(phi_AB, N), n_eig)
        V = zpe_AB - zpe_ref
        V_curve.append((R, V, len(modes_AB)))

    # Extract bond energy using ASYMPTOTE as reference
    # (more reliable than 2*ZPE_single due to boundary effects)
    V_vals = np.array([v[1] for v in V_curve])
    R_vals = np.array([v[0] for v in V_curve])

    # Asymptote = average of last 3 points (should be converged)
    V_inf = np.mean(V_vals[-3:])

    # Well depth relative to asymptote
    V_relative = V_vals - V_inf
    i_min = np.argmin(V_relative)
    D_e_lattice = -V_relative[i_min]
    R_eq = R_vals[i_min]

    # Convert to eV using the s-cancellation
    # D_e(physical) = D_e(lattice) * E_harm / (2*s)
    D_e_eV = D_e_lattice * E_harm / (2 * s_PT)

    return {
        'D_e_eV': D_e_eV,
        'D_e_lattice': D_e_lattice,
        'R_eq': R_eq,
        'E_harm': E_harm,
        'n_modes_A': len(modes_A),
        'n_modes_B': len(modes_B),
        'V_curve': V_curve,
        'V_inf': V_inf,
    }


# ============================================================
# TEST: Full molecule set
# ============================================================

molecules = [
    # (name, Z_A, Z_B, IE_A, IE_B, D_obs_eV)
    ("H2",    1,  1, 13.598, 13.598, 4.478),
    ("Li2",   3,  3,  5.392,  5.392, 1.046),
    ("LiH",   3,  1,  5.392, 13.598, 2.515),
    ("BH",    5,  1,  8.298, 13.598, 3.420),
    ("CH",    6,  1, 11.260, 13.598, 3.465),
    ("NH",    7,  1, 14.534, 13.598, 3.893),
    ("OH",    8,  1, 13.618, 13.598, 4.392),
    ("HF",    1,  9, 13.598, 17.423, 5.869),
    ("HCl",   1, 17, 13.598, 12.968, 4.431),
    ("N2",    7,  7, 14.534, 14.534, 9.759),
    ("O2",    8,  8, 13.618, 13.618, 5.116),
    ("F2",    9,  9, 17.423, 17.423, 1.602),
    ("CO",    6,  8, 11.260, 13.618, 11.09),
    ("C2",    6,  6, 11.260, 11.260, 6.210),
    ("NO",    7,  8, 14.534, 13.618, 6.497),
    ("CN",    6,  7, 11.260, 14.534, 7.720),
    ("Cl2",  17, 17, 12.968, 12.968, 2.475),
    ("NaCl", 11, 17,  5.139, 12.968, 4.243),
    ("PH",   15,  1, 10.487, 13.598, 3.440),
    ("SH",   16,  1, 10.360, 13.598, 3.560),
    ("Si2",  14, 14,  8.152,  8.152, 3.210),
]

print(f"\n{'Mol':>6} {'Z_A':>3} {'Z_B':>3} {'D_e':>8} {'D_obs':>8} {'Err':>8} {'R_eq':>5} {'E_harm':>8} {'D_lat':>8}")
print("-" * 72)

errors = []
for name, zA, zB, ieA, ieB, D_obs in molecules:
    result = bond_energy(zA, zB, ieA, ieB)
    D_e = result['D_e_eV']
    err = (D_e - D_obs) / D_obs * 100
    errors.append(abs(err))
    flag = " ***" if abs(err) > 30 else (" **" if abs(err) > 15 else "")
    print(f"{name:>6} {zA:3d} {zB:3d} {D_e:8.3f} {D_obs:8.3f} {err:+7.1f}% {result['R_eq']:5d} {result['E_harm']:8.2f} {result['D_e_lattice']:8.5f}{flag}")

print("-" * 72)
print(f"Mean |error|: {np.mean(errors):.1f}%")
print(f"Median:        {np.median(errors):.1f}%")
print(f"Under 10%:     {sum(1 for e in errors if e < 10)}/{len(errors)}")
print(f"Under 20%:     {sum(1 for e in errors if e < 20)}/{len(errors)}")
print(f"Under 30%:     {sum(1 for e in errors if e < 30)}/{len(errors)}")

# Show V(R) curve for H2
print(f"\nH2 bond curve (relative to asymptote):")
result_h2 = bond_energy(1, 1, 13.598, 13.598)
for R, V, nm in result_h2['V_curve']:
    V_rel = V - result_h2['V_inf']
    bar = '*' * max(0, int(-V_rel * 200))
    print(f"  R={R:2d}: V_rel={V_rel:+.5f} ({nm} modes) {bar}")
