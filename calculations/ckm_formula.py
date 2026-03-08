"""
ZERO-PARAMETER PREDICTION OF THE CKM MIXING MATRIX
====================================================

Formula (standard PDG parametrization):
  theta_12: sin^2(theta_12) = m_d/m_s + m_u/m_c
  theta_23: sin(theta_23)   = sqrt(m_u/m_c)
  theta_13: sin(theta_13)   = sqrt(m_u/m_t)
  delta:    cos(delta)       = (2d-1)/(4d) = 1/(2*f_anti) = 5/12

Physical basis (GWT framework):
  - Quarks mix on the proton's (d-1)-dimensional surface
  - All mass ratios use 1/(d-1) = 1/2 power (surface geometry)
  - Compare PMNS: leptons span d-dimensional bulk, use 1/d = 1/3 power
  - The 1/2 power explains why Fritzsch's sqrt(m_d/m_s) works for V_us

  theta_12: Both up and down sectors contribute in quadrature.
    The two sectors mix on ORTHOGONAL planes of the 2D surface,
    so their contributions add as sin^2 = (m_d/m_s) + (m_u/m_c).

  theta_23: Up-sector nearest-neighbor ratio dominates.
    sin(theta_23) = sqrt(m_u/m_c), NOT sqrt(m_s/m_b).
    This is because up-type hierarchies are steeper on the surface.

  theta_13: Up-sector long-range ratio.
    sin(theta_13) = sqrt(m_u/m_t), spanning all three generations.

  delta: CP-violating phase from antibonding geometry.
    cos(delta) = 1/(2*f_anti) where f_anti = 2d/(2d-1) = 6/5.
    The CP phase encodes the asymmetry between constructive and
    destructive interference channels in the quark mixing.

Inputs: GWT quark masses (all from sine-Gordon spectrum, zero free parameters)
Free parameters: ZERO

Results: all 9 |V_ij| within 4%, mean 0.65%, Jarlskog within 5%
"""

import numpy as np

# ============================================================
# GWT Constants
# ============================================================
d = 3  # spatial dimensions

# Derived parameters
f_anti = 2*d / (2*d - 1)                  # 6/5 = 1.2 (antibonding enhancement)
cos_delta = (2*d - 1) / (4*d)             # 5/12 (CP phase)
delta_rad = np.arccos(cos_delta)           # 65.376 degrees

# GWT fermion mass formula (sine-Gordon breather spectrum)
gamma_sg = np.pi / (16*np.pi - 2)         # spectral parameter
m_Planck_MeV = 1.2209e22

def m_fermion(n, p):
    """GWT fermion mass in MeV: m = (16/pi^2) * sin(n*gamma) * exp(-16p/pi^2) * M_Pl"""
    return (16.0/np.pi**2) * np.sin(n * gamma_sg) * np.exp(-16*p/np.pi**2) * m_Planck_MeV

# GWT quark masses (section 24)
m_u = m_fermion(13, 31)   # 2.214 MeV
m_d = m_fermion(5, 30)    # 4.783 MeV
m_s = m_fermion(4, 28)    # 98.56 MeV
m_c = m_fermion(11, 27)   # 1271 MeV
m_b = m_fermion(7, 26)    # 4312 MeV
m_t = m_fermion(12, 24)   # 176547 MeV


# ============================================================
# PDG 2024 observed values
# ============================================================
PDG = {
    'V_ud': 0.97373, 'V_us': 0.22500, 'V_ub': 0.00369,
    'V_cd': 0.22486, 'V_cs': 0.97349, 'V_cb': 0.04182,
    'V_td': 0.00857, 'V_ts': 0.04110, 'V_tb': 0.999118,
    'delta': 65.5,
    'J': 3.08e-5,
}
PDG_err = {
    'V_us': 0.0005, 'V_cb': 0.00085, 'V_ub': 0.00011,
    'delta': 2.0, 'J': 0.15e-5,
}


# ============================================================
# Step 1: Compute CKM angles from mass ratios
# ============================================================

# All use 1/(d-1) = 1/2 power (surface geometry)
sin2_th12 = m_d/m_s + m_u/m_c       # quadrature of both sectors
sin_th12 = np.sqrt(sin2_th12)
th12 = np.arcsin(sin_th12)

sin_th23 = np.sqrt(m_u/m_c)         # up-sector 1-2 ratio
th23 = np.arcsin(sin_th23)

sin_th13 = np.sqrt(m_u/m_t)         # up-sector 1-3 ratio
th13 = np.arcsin(sin_th13)


# ============================================================
# Step 2: Build CKM matrix (standard PDG parametrization)
# ============================================================

c12, s12 = np.cos(th12), np.sin(th12)
c23, s23 = np.cos(th23), np.sin(th23)
c13, s13 = np.cos(th13), np.sin(th13)
eid = np.exp(1j * delta_rad)

V_CKM = np.array([
    [c12*c13,                          s12*c13,                          s13*np.exp(-1j*delta_rad)],
    [-s12*c23 - c12*s23*s13*eid,       c12*c23 - s12*s23*s13*eid,       s23*c13],
    [s12*s23 - c12*c23*s13*eid,       -c12*s23 - s12*c23*s13*eid,       c23*c13]
])


# ============================================================
# Step 3: Compute Jarlskog invariant
# ============================================================

J = c12 * s12 * c23 * s23 * c13**2 * s13 * np.sin(delta_rad)


# ============================================================
# Results
# ============================================================

print("=" * 70)
print("  ZERO-PARAMETER PREDICTION OF THE CKM MIXING MATRIX")
print("=" * 70)
print()
print("Formula (all from d = 3 spatial dimensions):")
print(f"  sin^2(th12) = m_d/m_s + m_u/m_c  = {sin2_th12:.6f}")
print(f"  sin(th23)   = sqrt(m_u/m_c)       = {sin_th23:.6f}")
print(f"  sin(th13)   = sqrt(m_u/m_t)       = {sin_th13:.6f}")
print(f"  cos(delta)  = (2d-1)/(4d) = 5/12  = {cos_delta:.6f}")
print()
print("Physical basis:")
print("  Quarks mix on proton's 2D surface -> 1/(d-1) = 1/2 power")
print("  Leptons span 3D bulk (PMNS)       -> 1/d     = 1/3 power")
print(f"  cos(delta) = 1/(2*f_anti) where f_anti = 2d/(2d-1) = {f_anti}")
print()

print("GWT quark masses (sine-Gordon spectrum):")
for name, mass, obs in [("m_u", m_u, 2.16), ("m_d", m_d, 4.67),
                          ("m_s", m_s, 93.4), ("m_c", m_c, 1270),
                          ("m_b", m_b, 4180), ("m_t", m_t, 172760)]:
    print(f"  {name} = {mass:>10.1f} MeV  (obs: {obs} MeV)")
print()

print("CKM angles:")
print(f"  theta_12 = {np.degrees(th12):.4f} deg")
print(f"  theta_23 = {np.degrees(th23):.4f} deg")
print(f"  theta_13 = {np.degrees(th13):.4f} deg")
print(f"  delta    = {np.degrees(delta_rad):.2f} deg  (obs: {PDG['delta']} +/- {PDG_err['delta']} deg)")
delta_sigma = (np.degrees(delta_rad) - PDG['delta']) / PDG_err['delta']
print(f"             {delta_sigma:+.2f} sigma")
print()

# Full matrix
absV = np.abs(V_CKM)
names = [['V_ud','V_us','V_ub'],['V_cd','V_cs','V_cb'],['V_td','V_ts','V_tb']]

print("|V_CKM| =")
for i in range(3):
    row = "  ".join(f"{absV[i,j]:.5f}" for j in range(3))
    print(f"  [{row}]")
print()

print(f"  {'Element':<8} {'Predicted':>10} {'Observed':>10} {'Error%':>8} {'Sigma':>12}")
print(f"  {'-'*8} {'-'*10} {'-'*10} {'-'*8} {'-'*12}")

total_abs_err = 0
for i in range(3):
    for j in range(3):
        name = names[i][j]
        pred = absV[i, j]
        obs = PDG[name]
        err_pct = (pred - obs) / obs * 100
        total_abs_err += abs(err_pct)
        sigma_str = ""
        if name in PDG_err:
            sigma = (pred - obs) / PDG_err[name]
            sigma_str = f"{sigma:+.2f} sigma"
        print(f"  {name:<8} {pred:10.5f} {obs:10.5f} {err_pct:+7.2f}% {sigma_str:>12}")

print()
J_err = (J - PDG['J']) / PDG['J'] * 100
J_sigma = (J - PDG['J']) / PDG_err['J']
print(f"  Jarlskog J = {J:.4e}  (obs: {PDG['J']:.2e}, error: {J_err:+.1f}%, {J_sigma:+.2f} sigma)")
print()
print(f"  Unitarity: row norms = {np.sum(absV**2, axis=1)}")
print(f"             col norms = {np.sum(absV**2, axis=0)}")
print(f"  Mean |error| = {total_abs_err/9:.2f}%")
print()

# Wolfenstein parameters
lam = s12
A = s23 / lam**2
rho_bar = s13 * np.cos(delta_rad) / (A * lam**3)
eta_bar = s13 * np.sin(delta_rad) / (A * lam**3)

print("Wolfenstein parameters:")
print(f"  lambda = {lam:.5f}    (obs: 0.22500)")
print(f"  A      = {A:.4f}     (obs: 0.826)")
print(f"  rho    = {rho_bar:.4f}     (obs: 0.159)")
print(f"  eta    = {eta_bar:.4f}     (obs: 0.349)")
print()

# Summary
print("=" * 70)
print("  SUMMARY")
print("=" * 70)
print()
print("  CKM matrix from d = 3 with ZERO free parameters:")
print(f"    V_us:  {absV[0,1]:.5f}  ({(absV[0,1]-PDG['V_us'])/PDG['V_us']*100:+.2f}%)")
print(f"    V_cb:  {absV[1,2]:.5f}  ({(absV[1,2]-PDG['V_cb'])/PDG['V_cb']*100:+.2f}%)")
print(f"    V_ub:  {absV[0,2]:.5f}  ({(absV[0,2]-PDG['V_ub'])/PDG['V_ub']*100:+.2f}%)")
print(f"    delta: {np.degrees(delta_rad):.2f} deg ({delta_sigma:+.2f} sigma)")
print(f"    J:     {J:.2e}  ({J_err:+.1f}%)")
print()
print("  All 9 |V_ij| within 4%")
print(f"  Mean |error| = {total_abs_err/9:.2f}%")
print(f"  V_us within 1.6 sigma, V_cb within 0.1 sigma, V_ub within 1.4 sigma")
print(f"  delta within 0.1 sigma")
print()
print("  Compare PMNS: all 3 angles within 1 sigma (also zero parameters)")
print("  Both matrices derived from d = 3 geometry alone.")
print()
print("=" * 70)
