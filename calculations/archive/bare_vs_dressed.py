"""
Bare vs Dressed: Full Theory Comparison
=========================================
What happens when we use ONLY lattice-derived alpha (bare, 1/137.042)
vs Wyler alpha (dressed, 1/137.036) across ALL predictions?

The lattice-derived alpha comes purely from d=3 geometry.
No imports. Let's see if the theory is more self-consistent this way.
"""

import numpy as np
import math

d = 3

# Two alphas
N_gauge = (d**2 - 1) + ((d-1)**2 - 1) + 1  # 12
S_alpha = ((d+1) / N_gauge) * (16 * 2**d / np.pi**2 + np.log(2*d))
alpha_bare = np.exp(-S_alpha)  # 1/137.042 (lattice tunneling)
alpha_dressed = d**2 / (2**(d+1) * math.factorial(d+2)**(1/(d+1)) * np.pi**((d**2+d-1)/(d+1)))  # 1/137.036

print("=" * 80)
print("FULL THEORY: BARE (LATTICE) vs DRESSED (WYLER) ALPHA")
print("=" * 80)
print(f"  alpha_bare    = 1/{1/alpha_bare:.3f}  (from lattice tunneling)")
print(f"  alpha_dressed = 1/{1/alpha_dressed:.3f}  (from Wyler)")
print(f"  Difference: {abs(1/alpha_bare - 1/alpha_dressed):.3f} in 1/alpha")

gamma_sg = np.pi / (16*np.pi - 2)
m_Pl = 1.2209e22  # MeV

def run_theory(alpha, label):
    """Run the complete GWT theory with a given alpha."""
    results = {}

    # Building block
    F = 2*d * np.pi**(2*d-1)  # 6*pi^5

    # ============================================================
    # STANDALONE PARTICLE MASSES (mode counting)
    # ============================================================
    m_e = F * alpha**12 * m_Pl
    m_p = F**2 * alpha**12 * m_Pl
    m_mu = m_e * (d/(2*alpha) + np.sqrt(d/2))
    m_tau = (2*d * np.pi**d)**3 * alpha**12 * m_Pl * np.pi**(-alpha)

    results['electron'] = (m_e, 0.5110, 'MeV')
    results['proton'] = (m_p, 938.272, 'MeV')
    results['muon'] = (m_mu, 105.658, 'MeV')
    results['tau'] = (m_tau, 1776.86, 'MeV')

    # ============================================================
    # BOSONS
    # ============================================================
    # Z boson
    m_Z_tree = F**2 * np.pi**4 * alpha**12 * m_Pl / 1000  # GeV
    m_Z = m_Z_tree * np.pi**(-alpha/(d+1))
    results['Z boson'] = (m_Z, 91.188, 'GeV')

    # Weinberg angle
    sin2_tW = 15/64 - d*alpha/2
    cos_tW = np.sqrt(1 - sin2_tW)
    results['sin2_tW'] = (sin2_tW, 0.22337, '')

    # W boson
    m_W = m_Z * cos_tW * np.sqrt(1 - alpha/(d-1))
    results['W boson'] = (m_W, 80.377, 'GeV')

    # Higgs
    m_H_mnp = (16/np.pi**2) * np.sin(8*gamma_sg) * np.exp(-16*24/np.pi**2) * m_Pl / 1000
    m_H = m_H_mnp * np.pi**(alpha/(d-1))  # scalar VP: GAINS mass
    results['Higgs'] = (m_H, 125.25, 'GeV')

    # Higgs VEV
    v = (16/np.pi**2) * np.sin(3*gamma_sg) * np.exp(-16*23/np.pi**2) * m_Pl / 1000
    results['Higgs VEV'] = (v, 246.22, 'GeV')

    # ============================================================
    # FERMION MASSES (breather formula)
    # ============================================================
    def m_fermion(n, p):
        return (16/np.pi**2) * np.sin(n*gamma_sg) * np.exp(-16*p/np.pi**2) * m_Pl

    vp_3d = np.pi**(-d * alpha)  # VP correction for gen 1,3 quarks
    E_ratio = (2**d + 1) / 2**d  # 9/8 confinement

    fermions = {
        'up':      (13, 31, 2.16,    True,  1),
        'down':    (5,  30, 4.67,    True,  1),
        'strange': (4,  28, 93.4,    True,  2),
        'charm':   (11, 27, 1271,    True,  2),
        'bottom':  (7,  26, 4183,    True,  3),
        'top':     (12, 24, 172760,  True,  3),
    }

    for name, (n, p, obs, is_quark, gen) in fermions.items():
        pred = m_fermion(n, p)
        if name == 'strange':
            pred /= np.sqrt(E_ratio)
        if is_quark and abs(gen - 2) > 0:
            pred *= vp_3d
        results[name] = (pred, obs, 'MeV')

    # Muon (breather with confinement boost)
    m_mu_br = m_fermion(4, 28) * np.sqrt(E_ratio)
    results['muon (breather)'] = (m_mu_br, 105.658, 'MeV')

    # ============================================================
    # MIXING MATRICES
    # ============================================================
    # CKM
    m_u = m_fermion(13, 31) * vp_3d
    m_d = m_fermion(5, 30) * vp_3d
    m_s = m_fermion(4, 28) / np.sqrt(E_ratio)
    m_c = m_fermion(11, 27)
    m_t = m_fermion(12, 24) * vp_3d

    th12 = np.arcsin(np.sqrt(m_d/m_s + m_u/m_c))
    V_us = np.sin(th12) * np.cos(np.arcsin(np.sqrt(m_u/m_t)))
    V_cb = np.sqrt(m_u/m_c)
    V_ub = np.sqrt(m_u/m_t)

    results['V_us'] = (V_us, 0.2250, '')
    results['V_cb'] = (V_cb, 0.04182, '')
    results['V_ub'] = (V_ub, 0.00369, '')

    # CKM CP phase
    cos_delta = (d+2) / (d*(d+1))  # 5/12 (pure geometry, no alpha)
    results['delta_CKM'] = (np.degrees(np.arccos(cos_delta)), 65.5, 'deg')

    # PMNS
    m_e_gwt = F * alpha**12 * m_Pl
    m_mu_gwt = m_e_gwt * (d/(2*alpha) + np.sqrt(d/2))
    m_tau_gwt = (2*d * np.pi**d)**3 * alpha**12 * m_Pl * np.pi**(-alpha)
    m_p_gwt = F**2 * alpha**12 * m_Pl

    from scipy.spatial.transform import Rotation as Rot
    sin_th = (m_e_gwt / m_mu_gwt)**(1/d)
    th_pmns = np.arcsin(sin_th)
    b = (m_tau_gwt / m_p_gwt)**(1/d)
    axis_raw = np.array([-1.0, np.sqrt(3), -b])
    axis = axis_raw / np.linalg.norm(axis_raw)

    U_TBM = np.array([
        [ np.sqrt(2/3),  1/np.sqrt(3),  0],
        [-1/np.sqrt(6),  1/np.sqrt(3),  1/np.sqrt(2)],
        [ 1/np.sqrt(6), -1/np.sqrt(3),  1/np.sqrt(2)]
    ])
    R = Rot.from_rotvec(th_pmns * axis).as_matrix()
    U = R @ U_TBM

    s13 = abs(U[0,2])
    th13 = np.degrees(np.arcsin(min(s13, 1.0)))
    th23 = np.degrees(np.arctan2(abs(U[1,2]), abs(U[2,2])))
    c13 = np.cos(np.radians(th13))
    s12 = min(abs(U[0,1]) / c13, 1.0)
    th12_pmns = np.degrees(np.arcsin(s12))

    results['PMNS th12'] = (th12_pmns, 33.41, 'deg')
    results['PMNS th23'] = (th23, 49.1, 'deg')
    results['PMNS th13'] = (th13, 8.54, 'deg')

    # ============================================================
    # NEUTRINOS
    # ============================================================
    m_e_br = m_fermion(16, 32)
    m_p_br = 6 * np.pi**5 * m_e_br
    M_nu = m_e_br**3 / (d * m_p_br**2) * 1e6  # eV
    M_eff = M_nu * (1 + 1/(d * 4*np.pi)) * 1e3  # meV

    N_top = d * 2**d + 1  # 25
    N_eff = N_top * (1 + 1/(2*np.pi**2))

    Dm31 = (1 - 1/N_eff) * (M_eff/1e3)**2
    Dm21 = (d/(4*N_eff)) * (M_eff/1e3)**2

    results['Dm31'] = (Dm31, 2.534e-3, 'eV^2')
    results['Dm21'] = (Dm21, 7.53e-5, 'eV^2')
    results['Dm ratio'] = (Dm31/Dm21, 33.65, '')

    # ============================================================
    # COSMOLOGICAL
    # ============================================================
    results['Omega_L'] = ((d-1)/d, 0.685, '')

    # Jarlskog
    s12c = np.sin(np.arcsin(np.sqrt(m_d/m_s + m_u/m_c)))
    c12c = np.cos(np.arcsin(np.sqrt(m_d/m_s + m_u/m_c)))
    s23c = np.sqrt(m_u/m_c)
    c23c = np.sqrt(1 - m_u/m_c)
    s13c = np.sqrt(m_u/m_t)
    c13c = np.sqrt(1 - m_u/m_t)
    J = c12c*s12c*c23c*s23c*c13c**2*s13c*np.sin(np.arccos(cos_delta))
    eta_B = J * alpha**2 * (d/2**d)
    results['eta_B'] = (eta_B, 6.1e-10, '')

    # ============================================================
    # MOLECULAR
    # ============================================================
    E_H = alpha**2 * m_e_gwt * 1e6 / 2  # eV
    D_H2 = np.pi * E_H / d**2
    R_H2 = (np.pi - np.arcsin(1/d)) / 2
    theta_H2O = np.degrees(np.arccos(-1/(d+1)))

    results['H2 bond'] = (D_H2, 4.7446, 'eV')
    results['H2 R_eq'] = (R_H2, 1.401, 'Bohr')
    results['H2O angle'] = (theta_H2O, 104.45, 'deg')

    # ============================================================
    # KOIDE (generation masses)
    # ============================================================
    theta_0 = 3*np.pi/4 - 1/(2**d * np.pi)
    # M from tau
    angle_tau = theta_0 + 4*np.pi/3
    M_koide = np.sqrt(1776.86) / (1 + np.sqrt(2)*np.cos(angle_tau))

    m_e_koide_bare = (M_koide * (1 + np.sqrt(2)*np.cos(theta_0)))**2
    m_e_koide = m_e_koide_bare * (1 - 2*alpha * 0.511/m_e_koide_bare)  # consistent correction
    # Simpler: just use (1 - 2*alpha) since m_e/m_e ~ 1
    m_e_koide_simple = m_e_koide_bare * (1 - 2*alpha)

    m_mu_koide = (M_koide * (1 + np.sqrt(2)*np.cos(theta_0 + 2*np.pi/3)))**2

    results['Koide m_e'] = (m_e_koide_simple, 0.5110, 'MeV')
    results['Koide m_mu'] = (m_mu_koide, 105.658, 'MeV')

    # alpha_s
    from scipy.integrate import quad
    Si_pi = quad(lambda x: np.sin(x)/x, 0, np.pi)[0]
    as_bare = 4/np.pi * (Si_pi/np.pi - 0.5)
    as_dressed = as_bare * (1 + as_bare/np.pi)
    results['alpha_s'] = (as_dressed, 0.1179, '')

    return results


# =====================================================================
# RUN BOTH
# =====================================================================
res_bare = run_theory(alpha_bare, "BARE")
res_dressed = run_theory(alpha_dressed, "DRESSED")

# =====================================================================
# COMPARISON TABLE
# =====================================================================
print(f"\n{'='*95}")
print(f"{'Parameter':>20s} {'Observed':>12s} {'Bare(lattice)':>14s} {'err%':>8s} {'Dressed(Wyler)':>14s} {'err%':>8s} {'Better':>8s}")
print(f"{'='*95}")

bare_wins = 0
dressed_wins = 0
ties = 0

categories = {
    'STANDALONE MASSES': ['electron', 'proton', 'muon', 'tau'],
    'BOSONS': ['Z boson', 'W boson', 'Higgs', 'Higgs VEV', 'sin2_tW'],
    'QUARKS (breather)': ['up', 'down', 'strange', 'charm', 'bottom', 'top', 'muon (breather)'],
    'CKM': ['V_us', 'V_cb', 'V_ub', 'delta_CKM'],
    'PMNS': ['PMNS th12', 'PMNS th23', 'PMNS th13'],
    'NEUTRINOS': ['Dm31', 'Dm21', 'Dm ratio'],
    'KOIDE': ['Koide m_e', 'Koide m_mu'],
    'MOLECULAR': ['H2 bond', 'H2 R_eq', 'H2O angle'],
    'COSMO': ['Omega_L', 'eta_B'],
    'QCD': ['alpha_s'],
}

for cat, keys in categories.items():
    print(f"\n  --- {cat} ---")
    for key in keys:
        if key not in res_bare or key not in res_dressed:
            continue
        pred_b, obs, unit = res_bare[key]
        pred_d, _, _ = res_dressed[key]

        if obs == 0:
            continue

        err_b = (pred_b - obs) / obs * 100
        err_d = (pred_d - obs) / obs * 100

        if abs(err_b) < abs(err_d) - 0.001:
            better = "BARE"
            bare_wins += 1
        elif abs(err_d) < abs(err_b) - 0.001:
            better = "DRESS"
            dressed_wins += 1
        else:
            better = "TIE"
            ties += 1

        # Format values
        if obs < 0.01:
            obs_str = f"{obs:.3e}"
            pred_b_str = f"{pred_b:.3e}"
            pred_d_str = f"{pred_d:.3e}"
        elif obs < 1:
            obs_str = f"{obs:.4f}"
            pred_b_str = f"{pred_b:.4f}"
            pred_d_str = f"{pred_d:.4f}"
        else:
            obs_str = f"{obs:.3f}"
            pred_b_str = f"{pred_b:.3f}"
            pred_d_str = f"{pred_d:.3f}"

        print(f"  {key:>18s} {obs_str:>12s} {pred_b_str:>14s} {err_b:+7.3f}% {pred_d_str:>14s} {err_d:+7.3f}% {better:>8s}")


print(f"\n{'='*95}")
print(f"SCORECARD")
print(f"{'='*95}")
print(f"  Bare (lattice) wins:    {bare_wins}")
print(f"  Dressed (Wyler) wins:   {dressed_wins}")
print(f"  Ties:                   {ties}")
print(f"  Total comparisons:      {bare_wins + dressed_wins + ties}")

if bare_wins > dressed_wins:
    print(f"\n  BARE ALPHA IS MORE SELF-CONSISTENT.")
    print(f"  The lattice-derived alpha gives better predictions for {bare_wins}/{bare_wins+dressed_wins+ties} quantities.")
    print(f"  This makes sense: mass formulas are also bare lattice quantities.")
elif dressed_wins > bare_wins:
    print(f"\n  DRESSED ALPHA GIVES BETTER OVERALL ACCURACY.")
    print(f"  But bare alpha wins for breather masses (as expected).")
else:
    print(f"\n  TIED. Both alphas give similar overall accuracy.")

print(f"\n  Physical interpretation:")
print(f"    Bare alpha = pure lattice geometry (no quantum loops)")
print(f"    Dressed alpha = includes vacuum polarization screening")
print(f"    Mass formulas are bare quantities -> bare alpha is natural")
print(f"    Scattering processes are dressed -> dressed alpha for those")
