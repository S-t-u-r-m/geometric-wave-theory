"""
Particle-Point Bond Force Model
=================================
Two atoms as POINTS on a 1D axis. Three forces compete:

  1. Kink repulsion:     B * exp(-2*beta*R)     [short range, kink overlap]
  2. Breather attraction: A * exp(-s*beta*R)     [medium range, electron sharing]
  3. LP repulsion:       L * exp(-s*beta*R)      [medium range, occupied mode overlap]

where beta = sqrt(Z_eff), s = 0.17279 (universal PT parameter).

The bond energy = well depth of V(R) = V_rep + V_att + V_LP.
The equilibrium R_eq where dV/dR = 0.

All force parameters from the Lagrangian and Oh group theory.
No algebraic corrections — forces computed directly.
"""

import numpy as np
from scipy.optimize import minimize_scalar

PI = np.pi
d = 3

# Universal constants from the Lagrangian
s_PT = (-1 + np.sqrt(1 + 8/PI**2)) / 2    # 0.17279 Poeschl-Teller parameter
M_kink = 2**d / PI**2                       # 0.8106  kink mass
V_0 = 1.0 / PI**2                           # 0.1013  potential depth
alpha_bare = np.exp(-(2/6) * (2**7/PI**2 + np.log(6)))  # 1/137.042

# Oh channel fractions (all from d=3)
C_BOND = PI / d**2          # pi/9 = A1g fraction of T1u x T1u
W_PI = np.cos(PI / d)       # 1/2 = pi bond weight
F_RAD = (2*d - 1) / (2*d)   # 5/6 = radical coupling reduction
LP_I = (d**2 + 1) / d**3    # 10/27 = LP repulsion per pair
C_IONIC = 1 / (2*d + 1)     # 1/7 = ionic coupling


def bond_potential(R, Z_A, Z_B, IE_A, IE_B,
                   bond_order=1, n_lp=0, n_antibond=0,
                   is_radical=False, is_ionic=False):
    """
    Compute V(R) for two atoms at separation R.

    Three force components:
      V_rep:  kink-kink repulsion (always present)
      V_att:  breather tunneling attraction (bond order dependent)
      V_lp:   lone pair repulsion (LP count dependent)
      V_ion:  ionic attraction (IE asymmetry dependent)

    Returns: V_total, V_rep, V_att, V_lp, V_ion
    """
    if R < 0.5:
        return 1e6, 1e6, 0, 0, 0  # prevent collapse

    # Effective parameters
    Z_eff = np.sqrt(Z_A * Z_B)  # geometric mean
    beta = np.sqrt(Z_eff)
    E_harm = 2 * IE_A * IE_B / (IE_A + IE_B)

    # Decay rates (from Poeschl-Teller well structure)
    a_rep = 2 * beta       # kink overlap decay (fast)
    a_att = s_PT * beta    # breather tunneling decay (slow)

    # --- Force 1: Kink-kink repulsion ---
    # Two kink profiles overlap. Energy cost ~ M_kink per overlapping site.
    # Strength: proportional to kink mass and well depth
    B_rep = M_kink * Z_eff  # scales with Z (heavier = more repulsion)
    V_rep = B_rep * np.exp(-a_rep * R)

    # --- Force 2: Breather tunneling attraction ---
    # The electron (breather) tunnels between wells.
    # Bonding mode is lower energy by the tunneling splitting.
    # Strength: C_BOND * E_harm * bond_order
    # Each bond order adds one breather channel.
    A_att = C_BOND * E_harm * bond_order
    # Pi bonds are weaker by W_PI = cos(pi/d) = 1/2
    if bond_order > 1:
        # 1 sigma + (bo-1) pi bonds
        A_att = C_BOND * E_harm * (1 + (bond_order - 1) * W_PI)
    # Radical: directional coupling reduced by F_RAD = 5/6
    if is_radical:
        A_att *= F_RAD
    V_att = -A_att * np.exp(-a_att * R)

    # --- Force 3: Lone pair repulsion ---
    # Occupied breather modes resist overlap.
    # Each facing LP pair contributes LP_I = 10/27 of repulsion
    # relative to the bond attraction.
    # The LP force has the SAME range as the attraction (same tunneling decay)
    # but OPPOSITE sign (destructive interference).
    V_lp = 0
    if n_lp > 0:
        L_rep = C_BOND * E_harm * n_lp * LP_I
        V_lp = L_rep * np.exp(-a_att * R)

    # --- Force 4: Antibonding ---
    # Antibonding orbitals ADD repulsion at the attraction range
    if n_antibond > 0:
        f_anti = 2*d / (2*d - 1)  # 6/5 antibonding enhancement
        V_anti = C_BOND * E_harm * n_antibond * W_PI * f_anti * np.exp(-a_att * R)
        V_lp += V_anti

    # --- Force 5: Ionic attraction ---
    # Charge transfer between asymmetric wells.
    # Ionic energy ~ delta_IE^2 / E_H, screened by distance.
    V_ion = 0
    if is_ionic or abs(IE_A - IE_B) > 2.0:
        delta_IE = abs(IE_A - IE_B)
        E_H = alpha_bare**2 * 0.51100e6 / 2
        # Ionic decay: slower than covalent (Coulomb-like at medium range)
        a_ion = s_PT * beta / 2  # longer range than covalent
        c_ion = C_IONIC
        if delta_IE / E_harm > 0.3:  # strongly ionic
            c_ion = d / (2*d + 1)  # 3/7 enhanced
        V_ion = -c_ion * delta_IE * np.exp(-a_ion * R)

    V_total = V_rep + V_att + V_lp + V_ion
    return V_total, V_rep, V_att, V_lp, V_ion


def compute_bond_energy(Z_A, Z_B, IE_A, IE_B,
                        bond_order=1, n_lp=0, n_antibond=0,
                        is_radical=False, is_ionic=False):
    """
    Find equilibrium and bond energy from the force model.
    """
    # Find minimum of V(R)
    def V_func(R):
        return bond_potential(R, Z_A, Z_B, IE_A, IE_B,
                             bond_order, n_lp, n_antibond,
                             is_radical, is_ionic)[0]

    # Scan to find approximate minimum
    R_scan = np.linspace(0.5, 15.0, 300)
    V_scan = np.array([V_func(R) for R in R_scan])

    # Asymptote
    V_inf = V_func(15.0)

    # Find minimum
    i_min = np.argmin(V_scan)
    if V_scan[i_min] >= V_inf - 1e-10:
        # No bound state
        return 0, 0, V_inf, {}

    # Refine with optimizer
    result = minimize_scalar(V_func, bounds=(0.5, 12.0), method='bounded')
    R_eq = result.x
    V_min = result.fun

    D_e = -(V_min - V_inf)

    # Get force breakdown at equilibrium
    V_tot, V_rep, V_att, V_lp, V_ion = bond_potential(
        R_eq, Z_A, Z_B, IE_A, IE_B,
        bond_order, n_lp, n_antibond, is_radical, is_ionic)

    return D_e, R_eq, V_inf, {
        'V_rep': V_rep, 'V_att': V_att, 'V_lp': V_lp, 'V_ion': V_ion
    }


# ============================================================
# MOLECULE DATABASE
# ============================================================
# (name, Z_A, Z_B, IE_A, IE_B, bond_order, n_lp, n_antibond, radical, ionic, D_obs)

molecules = [
    ("H2",    1,  1, 13.598, 13.598, 1, 0, 0, False, False, 4.478),
    ("Li2",   3,  3,  5.392,  5.392, 1, 0, 0, False, False, 1.046),
    ("LiH",   3,  1,  5.392, 13.598, 1, 0, 0, False, True,  2.515),
    ("BH",    5,  1,  8.298, 13.598, 1, 0, 0, False, False, 3.420),
    ("CH",    6,  1, 11.260, 13.598, 1, 0, 0, True,  False, 3.465),
    ("NH",    7,  1, 14.534, 13.598, 1, 1, 0, True,  False, 3.893),
    ("OH",    8,  1, 13.618, 13.598, 1, 2, 0, True,  False, 4.392),
    ("HF",    1,  9, 13.598, 17.423, 1, 3, 0, False, True,  5.869),
    ("HCl",   1, 17, 13.598, 12.968, 1, 3, 0, False, True,  4.431),
    ("N2",    7,  7, 14.534, 14.534, 3, 0, 0, False, False, 9.759),
    ("O2",    8,  8, 13.618, 13.618, 2, 2, 1, True,  False, 5.116),
    ("F2",    9,  9, 17.423, 17.423, 1, 3, 0, False, False, 1.602),
    ("CO",    6,  8, 11.260, 13.618, 3, 1, 0, False, True, 11.090),
    ("C2",    6,  6, 11.260, 11.260, 2, 0, 0, False, False, 6.210),
    ("NO",    7,  8, 14.534, 13.618, 2, 1, 0, True,  False, 6.497),
    ("CN",    6,  7, 11.260, 14.534, 3, 0, 0, True,  False, 7.720),
    ("Cl2",  17, 17, 12.968, 12.968, 1, 2, 0, False, False, 2.475),
    ("NaCl", 11, 17,  5.139, 12.968, 1, 3, 0, False, True,  4.243),
    ("PH",   15,  1, 10.487, 13.598, 1, 1, 0, True,  False, 3.440),
    ("SH",   16,  1, 10.360, 13.598, 1, 2, 0, True,  False, 3.560),
    ("Si2",  14, 14,  8.152,  8.152, 2, 0, 0, False, False, 3.210),
    ("HBr",   1, 35, 13.598, 11.814, 1, 3, 0, False, True,  3.758),
    ("NH3",   1,  7, 13.598, 14.534, 1, 1, 0, False, False, 4.514),
    ("H2O",   1,  8, 13.598, 13.618, 1, 2, 0, False, False, 5.110),
]

# ============================================================
# RUN
# ============================================================
print("=" * 75)
print("PARTICLE-POINT BOND FORCE MODEL")
print("=" * 75)
print(f"Forces: kink repulsion + breather attraction + LP repulsion + ionic")
print(f"All parameters from d={d} Lagrangian. Zero fitting.")
print()

print(f"{'Mol':>6} {'D_e':>7} {'D_obs':>7} {'Err':>7} {'R_eq':>5} "
      f"{'V_rep':>7} {'V_att':>7} {'V_lp':>7} {'V_ion':>7}")
print("-" * 75)

errors = []
for (name, zA, zB, ieA, ieB, bo, nlp, nab, rad, ionic, D_obs) in molecules:
    D_e, R_eq, V_inf, forces = compute_bond_energy(
        zA, zB, ieA, ieB, bo, nlp, nab, rad, ionic)

    if D_e > 0:
        err = (D_e - D_obs) / D_obs * 100
        errors.append(abs(err))
        flag = " ***" if abs(err) > 30 else ""
        print(f"{name:>6} {D_e:7.3f} {D_obs:7.3f} {err:+6.1f}% {R_eq:5.2f} "
              f"{forces['V_rep']:7.4f} {forces['V_att']:7.3f} "
              f"{forces['V_lp']:7.3f} {forces['V_ion']:7.3f}{flag}")
    else:
        errors.append(100)
        print(f"{name:>6}   NONE {D_obs:7.3f}  -100% {'---':>5} "
              f"  ---     ---     ---     ---  ***")

print("-" * 75)
print(f"Mean |error|: {np.mean(errors):.1f}%")
print(f"Median:        {np.median(errors):.1f}%")
print(f"Under 5%:      {sum(1 for e in errors if e < 5)}/{len(errors)}")
print(f"Under 10%:     {sum(1 for e in errors if e < 10)}/{len(errors)}")
print(f"Under 20%:     {sum(1 for e in errors if e < 20)}/{len(errors)}")
print(f"Under 30%:     {sum(1 for e in errors if e < 30)}/{len(errors)}")

# Show force profile for H2
print()
print("H2 FORCE PROFILE:")
print(f"{'R':>5} {'V_tot':>8} {'V_rep':>8} {'V_att':>8} {'V_lp':>8}")
for R in np.arange(0.5, 8.0, 0.5):
    vt, vr, va, vl, vi = bond_potential(R, 1, 1, 13.598, 13.598)
    print(f"{R:5.1f} {vt:8.4f} {vr:8.4f} {va:8.4f} {vl:8.4f}")
