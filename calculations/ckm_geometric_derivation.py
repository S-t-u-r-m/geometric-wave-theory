"""
CKM Matrix — Unified Geometric Derivation
============================================
Goal: Derive the full CKM matrix from zero free parameters,
mirroring the PMNS derivation structure.

PMNS pattern (worked, all 3 angles within 1σ):
  U_PMNS = R(θ_corr, axis) × U_TBM
  θ_corr = arcsin((m_e/m_μ)^(1/d))    [1/d = 1/3 power, bulk geometry]
  U_TBM = tribimaximal matrix          [democratic base]
  axis = (−1, √3, −1) × wrapping      [Wyler S³ correction]

CKM key insight (from plan):
  Quarks use 1/(d−1) = 1/2 power (surface geometry, all confined in proton)
  Leptons use 1/d = 1/3 power (bulk geometry, cross proton boundary)

Current CKM predictions (ad hoc formulas):
  V_us = √(m_d/m_s + m_u/m_c) = 0.2242  (obs 0.2243, −0.04%)  <-- excellent
  V_cb = √(2/d) × λ²           = 0.0410  (obs 0.0408, +0.6%)   <-- good
  V_ub = √(m_u/m_t)            = 0.00354 (obs 0.00369, −4.0%)   <-- weakest
  δ_CKM = arccos(5/12)         = 65.38°  (obs 65.5°, −0.2%)     <-- excellent

Question: Can we unify these into a SINGLE rotation construction?
  V_CKM = R(θ_corr, axis) × V_Cabibbo
  or
  V_CKM = U_up† × U_down  (mismatch between up and down sectors)
"""

import numpy as np
import math

# ============================================================
# GWT CONSTANTS AND MASSES
# ============================================================
d = 3
gamma_gwt = np.pi / (16 * np.pi - 2)
m_Planck_MeV = 1.2209e22

def m_gwt(n, p):
    return (16.0/np.pi**2) * np.sin(n * gamma_gwt) * np.exp(-16*p/np.pi**2) * m_Planck_MeV

# GWT §24 quark masses
m_u = m_gwt(13, 31)  # 2.214 MeV
m_d = m_gwt(5, 30)   # 4.783 MeV
m_s = m_gwt(4, 28)   # 98.56 MeV
m_c = m_gwt(11, 27)  # 1271 MeV
m_b = m_gwt(7, 26)   # 4312 MeV
m_t = m_gwt(12, 24)  # 176547 MeV

print("GWT Quark Masses:")
print(f"  u = {m_u:.3f}, d = {m_d:.3f}, s = {m_s:.3f}")
print(f"  c = {m_c:.1f}, b = {m_b:.1f}, t = {m_t:.0f}")

# PDG observed CKM
V_us_obs = 0.2243
V_cb_obs = 0.0408
V_ub_obs = 0.00369
V_td_obs = 0.0082
V_ts_obs = 0.0394
V_tb_obs = 0.9991
delta_CKM_obs = 65.5  # degrees

# ============================================================
# APPROACH 1: V_CKM = U_up† × U_down (sector mismatch)
# ============================================================
print("\n" + "=" * 80)
print("APPROACH 1: V_CKM = U_up† × U_down")
print("=" * 80)

def rotation_12(theta):
    """Rotation in the 1-2 plane."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])

def rotation_23(theta):
    """Rotation in the 2-3 plane."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])

def rotation_13(theta):
    """Rotation in the 1-3 plane."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

def extract_ckm_angles(V):
    """Extract standard parametrization angles from CKM matrix."""
    s13 = abs(V[0, 2])
    c13 = np.sqrt(1 - s13**2)
    s12 = abs(V[0, 1]) / c13 if c13 > 0 else 0
    s23 = abs(V[1, 2]) / c13 if c13 > 0 else 0
    theta12 = np.degrees(np.arcsin(s12))
    theta23 = np.degrees(np.arcsin(s23))
    theta13 = np.degrees(np.arcsin(s13))
    return theta12, theta23, theta13

# Down-type rotation angle: arcsin(√(m_d/m_s))
theta_down_12 = np.arcsin(np.sqrt(m_d / m_s))
theta_down_23 = np.arcsin(np.sqrt(m_s / m_b))
theta_down_13 = np.arcsin(np.sqrt(m_d / m_b))

# Up-type rotation angle: arcsin(√(m_u/m_c))
theta_up_12 = np.arcsin(np.sqrt(m_u / m_c))
theta_up_23 = np.arcsin(np.sqrt(m_c / m_t))
theta_up_13 = np.arcsin(np.sqrt(m_u / m_t))

print(f"\nDown-type angles: th12={np.degrees(theta_down_12):.3f} deg, "
      f"th23={np.degrees(theta_down_23):.3f} deg, th13={np.degrees(theta_down_13):.3f} deg")
print(f"Up-type angles:   th12={np.degrees(theta_up_12):.3f} deg, "
      f"th23={np.degrees(theta_up_23):.3f} deg, th13={np.degrees(theta_up_13):.3f} deg")

# Test 1a: Simple 12-plane rotations
U_down_12 = rotation_12(theta_down_12)
U_up_12 = rotation_12(theta_up_12)
V_test = U_up_12.T @ U_down_12

print(f"\n  Test 1a: V = R12_up† × R12_down")
print(f"  |V_us| = {abs(V_test[0,1]):.5f} (obs {V_us_obs})")
print(f"  |V_cb| = {abs(V_test[1,2]):.5f} (obs {V_cb_obs})")
print(f"  |V_ub| = {abs(V_test[0,2]):.5f} (obs {V_ub_obs})")

# Test 1b: Full 3-angle rotations per sector
U_down_full = rotation_23(theta_down_23) @ rotation_13(theta_down_13) @ rotation_12(theta_down_12)
U_up_full = rotation_23(theta_up_23) @ rotation_13(theta_up_13) @ rotation_12(theta_up_12)
V_full = U_up_full.T @ U_down_full

print(f"\n  Test 1b: V = U_up†(full) × U_down(full)")
print(f"  |V| = ")
for row in abs(V_full):
    print(f"    [{row[0]:.6f}  {row[1]:.6f}  {row[2]:.6f}]")
print(f"  |V_us| = {abs(V_full[0,1]):.5f} (obs {V_us_obs}, err {(abs(V_full[0,1])-V_us_obs)/V_us_obs*100:+.2f}%)")
print(f"  |V_cb| = {abs(V_full[1,2]):.5f} (obs {V_cb_obs}, err {(abs(V_full[1,2])-V_cb_obs)/V_cb_obs*100:+.2f}%)")
print(f"  |V_ub| = {abs(V_full[0,2]):.5f} (obs {V_ub_obs}, err {(abs(V_full[0,2])-V_ub_obs)/V_ub_obs*100:+.2f}%)")


# ============================================================
# APPROACH 2: V_CKM = R(θ_corr, axis) × V_Cabibbo
# (mirror PMNS structure)
# ============================================================
print("\n" + "=" * 80)
print("APPROACH 2: V_CKM = R(θ_corr, axis) × V_Cabibbo")
print("=" * 80)

# Leading-order Cabibbo rotation
theta_C = np.arcsin(np.sqrt(m_d / m_s))  # Fritzsch-like
V_Cabibbo = rotation_12(theta_C)

print(f"  Cabibbo angle: θ_C = arcsin(√(m_d/m_s)) = {np.degrees(theta_C):.3f}°")
print(f"  Leading V_us = sin(θ_C) = {np.sin(theta_C):.5f} (obs {V_us_obs})")

# Correction angle candidates (mirroring PMNS θ_corr)
correction_candidates = {
    "arcsin((m_u/m_c)^(1/2))": np.arcsin(np.sqrt(m_u / m_c)),
    "arcsin((m_u/m_c)^(1/3))": np.arcsin((m_u / m_c)**(1.0/3)),
    "arcsin((m_d/m_b)^(1/2))": np.arcsin(np.sqrt(m_d / m_b)),
    "arcsin((m_s/m_b)^(1/2))": np.arcsin(np.sqrt(m_s / m_b)),
    "arcsin((m_u/m_t)^(1/2))": np.arcsin(np.sqrt(m_u / m_t)),
    "arcsin((m_d/m_s)^(1/d))": np.arcsin((m_d / m_s)**(1.0/d)),
    "arcsin((m_u/m_c)^(1/(d-1)))": np.arcsin((m_u / m_c)**(1.0/(d-1))),
    "arcsin((m_u/m_t)^(1/d))": np.arcsin((m_u / m_t)**(1.0/d)),
}

print(f"\n  Correction angle candidates:")
for desc, theta in sorted(correction_candidates.items(), key=lambda x: x[1]):
    print(f"    {desc:.<45s} = {np.degrees(theta):8.4f}°")


# ============================================================
# APPROACH 3: Quadrature formula as a rotation
# ============================================================
print("\n" + "=" * 80)
print("APPROACH 3: Understand the quadrature formula geometrically")
print("=" * 80)

# The current V_us formula: √(m_d/m_s + m_u/m_c)
# This is actually: √(sin²θ_d + sin²θ_u) where θ_d and θ_u are
# the individual sector rotation angles in the 12-plane.
# Geometrically: two perpendicular rotations add in quadrature.

sin_theta_d = np.sqrt(m_d / m_s)  # = 0.2203
sin_theta_u = np.sqrt(m_u / m_c)  # = 0.04174
lambda_quad = np.sqrt(sin_theta_d**2 + sin_theta_u**2)

print(f"  sin(θ_d) = √(m_d/m_s) = {sin_theta_d:.5f}")
print(f"  sin(θ_u) = √(m_u/m_c) = {sin_theta_u:.5f}")
print(f"  λ = √(sin²θ_d + sin²θ_u) = {lambda_quad:.5f} (obs V_us = {V_us_obs})")
print(f"  This is Pythagoras: two perpendicular rotations in flavor space")

# Now for V_cb and V_ub, can we use the same quadrature structure?
# V_cb should involve 23 rotations from both sectors
sin_theta_d_23 = np.sqrt(m_s / m_b)  # = 0.1512
sin_theta_u_23 = np.sqrt(m_c / m_t)  # = 0.08484

V_cb_quad = np.sqrt(sin_theta_d_23**2 + sin_theta_u_23**2)
print(f"\n  V_cb quadrature test:")
print(f"    sin(θ_d_23) = √(m_s/m_b) = {sin_theta_d_23:.5f}")
print(f"    sin(θ_u_23) = √(m_c/m_t) = {sin_theta_u_23:.5f}")
print(f"    √(sin²θ_d_23 + sin²θ_u_23) = {V_cb_quad:.5f} (obs {V_cb_obs})")
print(f"    Error: {(V_cb_quad - V_cb_obs)/V_cb_obs*100:+.1f}%  *** TOO LARGE ***")

# The quadrature doesn't work for V_cb directly because it's 2nd order
# V_cb ≈ A·λ² in Wolfenstein. Let's try: V_cb = √(m_s/m_b · m_u/m_c)
V_cb_geom = np.sqrt((m_s/m_b) * (m_u/m_c))
print(f"\n  V_cb = √((m_s/m_b)·(m_u/m_c)) = {V_cb_geom:.5f} (obs {V_cb_obs}, err {(V_cb_geom-V_cb_obs)/V_cb_obs*100:+.1f}%)")

# V_cb = √(m_s/m_b) · √(m_u/m_c) — product not quadrature
V_cb_prod = np.sqrt(m_s/m_b) * np.sqrt(m_u/m_c)
print(f"  V_cb = √(m_s/m_b)·√(m_u/m_c) = {V_cb_prod:.5f} (obs {V_cb_obs}, err {(V_cb_prod-V_cb_obs)/V_cb_obs*100:+.1f}%)")

# Try: V_cb from the mismatch approach (approach 1b)
print(f"  V_cb from U_up†·U_down = {abs(V_full[1,2]):.5f} (obs {V_cb_obs})")


# ============================================================
# APPROACH 4: Systematic scan of (m_i/m_j)^(1/p) formulas
# ============================================================
print("\n" + "=" * 80)
print("APPROACH 4: Systematic scan of mass ratio powers")
print("=" * 80)

masses_down = {"d": m_d, "s": m_s, "b": m_b}
masses_up = {"u": m_u, "c": m_c, "t": m_t}

targets = {
    "V_us": (V_us_obs, 0.0005),
    "V_cb": (V_cb_obs, 0.0014),
    "V_ub": (V_ub_obs, 0.00011),
    "V_td": (V_td_obs, 0.0003),
}

print(f"\n  Scanning (m_i/m_j)^(1/p) for each CKM element...")
print(f"  {'Target':>6s}  {'Formula':>35s}  {'Value':>10s}  {'Observed':>10s}  {'Error':>8s}")
print("  " + "-" * 80)

for target_name, (obs_val, obs_err) in targets.items():
    best_formulas = []

    # Test all mass ratios with powers 1/2, 1/3, 1/4, 1/6
    all_masses = {**masses_down, **masses_up}
    for n1, m1 in all_masses.items():
        for n2, m2 in all_masses.items():
            if m1 >= m2:
                continue
            ratio = m1 / m2
            for p_desc, p in [("1/2", 0.5), ("1/3", 1/3), ("1/4", 0.25), ("1/6", 1/6)]:
                val = ratio**p
                err_pct = abs(val - obs_val) / obs_val * 100
                if err_pct < 10:
                    best_formulas.append((f"(m_{n1}/m_{n2})^({p_desc})", val, err_pct))
            # Also test sqrt(ratio1 + ratio2) for quadrature

    # Also test products
    for n1, m1 in masses_down.items():
        for n2, m2 in masses_down.items():
            if m1 >= m2:
                continue
            for n3, m3 in masses_up.items():
                for n4, m4 in masses_up.items():
                    if m3 >= m4:
                        continue
                    # sqrt(r1 + r2) quadrature
                    val = np.sqrt(m1/m2 + m3/m4)
                    err_pct = abs(val - obs_val) / obs_val * 100
                    if err_pct < 5:
                        best_formulas.append((f"√(m_{n1}/m_{n2} + m_{n3}/m_{n4})", val, err_pct))
                    # sqrt(r1 * r2) geometric
                    val2 = np.sqrt((m1/m2) * (m3/m4))
                    err_pct2 = abs(val2 - obs_val) / obs_val * 100
                    if err_pct2 < 5:
                        best_formulas.append((f"√(m_{n1}/m_{n2} · m_{n3}/m_{n4})", val2, err_pct2))
                    # r1 * r2 direct product
                    val3 = np.sqrt(m1/m2) * np.sqrt(m3/m4)
                    err_pct3 = abs(val3 - obs_val) / obs_val * 100
                    if err_pct3 < 5:
                        best_formulas.append((f"√(m_{n1}/m_{n2})·√(m_{n3}/m_{n4})", val3, err_pct3))

    best_formulas.sort(key=lambda x: x[2])
    for desc, val, err in best_formulas[:5]:
        sigma = abs(val - obs_val) / obs_err if obs_err > 0 else 0
        print(f"  {target_name:>6s}  {desc:>35s}  {val:10.6f}  {obs_val:10.6f}  {err:+7.2f}% ({sigma:.1f}σ)")
    print()


# ============================================================
# APPROACH 5: Build CKM from the "frequency ratio" picture
# ============================================================
print("=" * 80)
print("APPROACH 5: Frequency ratio construction")
print("=" * 80)
print("CKM elements as coupling amplitudes between wave modes.")
print("Two modes couple proportionally to the overlap of their")
print("frequency distributions on the (d-1)-dimensional surface.")
print()

# The key formula: V_ij ~ (m_i/m_j)^(1/(d-1)) for surface coupling
# where i is the lighter quark and j is the heavier.
# The 1/(d-1) = 1/2 power comes from d-1 = 2 dimensional surface.

print(f"Surface coupling: V_ij = (m_light/m_heavy)^(1/(d-1)) = (m_i/m_j)^(1/2)")
print(f"Bulk coupling:    V_ij = (m_light/m_heavy)^(1/d) = (m_i/m_j)^(1/3)")
print()

# For CKM 12: two perpendicular contributions (down-type and up-type)
# V_us² = (m_d/m_s)^(2/(d-1)) + (m_u/m_c)^(2/(d-1))
# At d=3: V_us² = m_d/m_s + m_u/m_c  ← this is exactly the quadrature formula!

V_us_surface = np.sqrt((m_d/m_s) + (m_u/m_c))
print(f"V_us = √((m_d/m_s)^(2/(d-1)) + (m_u/m_c)^(2/(d-1)))")
print(f"     = √(m_d/m_s + m_u/m_c) = {V_us_surface:.5f} (obs {V_us_obs})")
print(f"  This IS the existing formula — validated as surface geometry.")
print()

# For CKM 23: same structure but 2-3 generation
# V_cb² = (m_s/m_b)^(2/(d-1)) + (m_c/m_t)^(2/(d-1))
# At d=3: V_cb² = m_s/m_b + m_c/m_t
V_cb_surface = np.sqrt((m_s/m_b) + (m_c/m_t))
print(f"V_cb = √((m_s/m_b)^(2/(d-1)) + (m_c/m_t)^(2/(d-1)))")
print(f"     = √(m_s/m_b + m_c/m_t) = {V_cb_surface:.5f} (obs {V_cb_obs})")
print(f"  Error: {(V_cb_surface-V_cb_obs)/V_cb_obs*100:+.1f}% — TOO LARGE (factor ~4)")
print()

# Hmm, the direct quadrature for V_cb gives ~0.17, way too large.
# V_cb is SECOND order in the hierarchy. The structure must be:
# V_cb ~ V_us² × (geometric factor)
# OR: V_cb is the 23-element of the MISMATCH matrix, not a direct ratio.

# Let me try the mismatch more carefully.
# Each sector has a principal rotation in the 12 plane.
# The 23 and 13 elements come from the commutator of the two rotations.

# For the mismatch V = U_up† × U_down with ONLY 12 rotations:
# V_us = sin(θ_d - θ_u) ≈ sin(θ_d) - cos(θ_d)sin(θ_u)
# V_cb = 0 (no 23 mixing from 12 rotations alone)
# So we NEED the 23 rotation in at least one sector.

# Key: V_cb comes from the 23-plane mismatch.
# θ_down_23 = arcsin(√(m_s/m_b)) and θ_up_23 = arcsin(√(m_c/m_t))
# V_cb ≈ sin(θ_down_23 - θ_up_23) for small angles

V_cb_mismatch = np.sin(theta_down_23 - theta_up_23)
print(f"V_cb = sin(θ_d_23 - θ_u_23)")
print(f"     = sin({np.degrees(theta_down_23):.3f}° - {np.degrees(theta_up_23):.3f}°)")
print(f"     = sin({np.degrees(theta_down_23 - theta_up_23):.3f}°)")
print(f"     = {V_cb_mismatch:.5f} (obs {V_cb_obs})")
print(f"  Error: {(V_cb_mismatch - V_cb_obs)/V_cb_obs*100:+.1f}%")
print()

# V_ub from 13-plane mismatch
V_ub_mismatch = np.sin(theta_down_13 - theta_up_13)
print(f"V_ub = sin(θ_d_13 - θ_u_13)")
print(f"     = sin({np.degrees(theta_down_13):.3f}° - {np.degrees(theta_up_13):.3f}°)")
print(f"     = {V_ub_mismatch:.5f} (obs {V_ub_obs})")
print(f"  Error: {(V_ub_mismatch - V_ub_obs)/V_ub_obs*100:+.1f}%")
print()


# ============================================================
# BEST UNIFIED CONSTRUCTION
# ============================================================
print("=" * 80)
print("UNIFIED CKM CONSTRUCTION")
print("=" * 80)

# The mismatch V = U_up† × U_down using √(m_i/m_j) angles seems
# to be the cleanest. Let me build the full matrix.

print("\nBuilding V_CKM = U_up†(θ12, θ23, θ13) × U_down(θ12, θ23, θ13)")
print(f"where θ_ij = arcsin(√(m_light/m_heavy)) in each sector\n")

# Add CP phase: the mismatch of two real rotations is real (no CP).
# CP phase must come from the GEOMETRY — the tetrahedral angle.
# Insert it as a phase in the 13 rotation of the down sector.
delta = np.arccos(5.0/12)  # = 65.38° from cos(δ) = (d+2)/(d(d+1))

def rotation_13_cp(theta, delta_phase):
    """Rotation in 1-3 plane with CP phase."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([
        [c, 0, s * np.exp(-1j * delta_phase)],
        [0, 1, 0],
        [-s * np.exp(1j * delta_phase), 0, c]
    ])

# Build U_down with CP phase
U_down_cp = rotation_23(theta_down_23) @ rotation_13_cp(theta_down_13, delta) @ rotation_12(theta_down_12)
U_up_real = rotation_23(theta_up_23) @ rotation_13(theta_up_13) @ rotation_12(theta_up_12)

V_unified = U_up_real.conj().T @ U_down_cp

print("|V_CKM| (unified geometric construction):")
V_abs = np.abs(V_unified)
labels = [["V_ud", "V_us", "V_ub"], ["V_cd", "V_cs", "V_cb"], ["V_td", "V_ts", "V_tb"]]
obs_matrix = [
    [0.9742, 0.2243, 0.00369],
    [0.218, 0.997, 0.0408],
    [0.0082, 0.0394, 0.9991]
]

print(f"\n  {'Element':>6s}  {'Predicted':>10s}  {'Observed':>10s}  {'Error':>8s}")
print("  " + "-" * 40)

total_err = 0
n_elem = 0
for i in range(3):
    for j in range(3):
        pred = V_abs[i, j]
        obs = obs_matrix[i][j]
        err = (pred - obs) / obs * 100
        print(f"  {labels[i][j]:>6s}  {pred:10.6f}  {obs:10.6f}  {err:+7.2f}%")
        total_err += abs(err)
        n_elem += 1

mean_err = total_err / n_elem
print(f"\n  Mean |error| across all 9 elements: {mean_err:.2f}%")

# Check unitarity
print(f"\n  Unitarity check:")
VV = V_unified @ V_unified.conj().T
for i in range(3):
    print(f"    Row {i+1} norm: {np.abs(VV[i,i]):.8f}")

# Jarlskog invariant
J = np.imag(V_unified[0,0] * V_unified[1,1] * np.conj(V_unified[0,1]) * np.conj(V_unified[1,0]))
print(f"\n  Jarlskog invariant: J = {J:.2e} (obs 3.08e-5)")

# Extract CP phase
if abs(V_unified[0,2]) > 0:
    delta_extracted = np.degrees(np.angle(V_unified[0,2]))
    print(f"  Extracted CP phase: {delta_extracted:.1f}° (input: {np.degrees(delta):.1f}°)")


# ============================================================
# COMPARISON: Current vs Unified
# ============================================================
print("\n" + "=" * 80)
print("COMPARISON: CURRENT AD HOC vs UNIFIED GEOMETRIC")
print("=" * 80)

current = {
    "V_us": (np.sqrt(m_d/m_s + m_u/m_c), V_us_obs, "√(m_d/m_s + m_u/m_c)"),
    "V_cb": (np.sqrt(2.0/d) * (m_d/m_s + m_u/m_c), V_cb_obs, "√(2/d)·λ²"),
    "V_ub": (np.sqrt(m_u/m_t), V_ub_obs, "√(m_u/m_t)"),
}

unified = {
    "V_us": (V_abs[0, 1], V_us_obs, "U_up†·U_down [0,1]"),
    "V_cb": (V_abs[1, 2], V_cb_obs, "U_up†·U_down [1,2]"),
    "V_ub": (V_abs[0, 2], V_ub_obs, "U_up†·U_down [0,2]"),
}

print(f"\n  {'Element':>6s}  {'Current':>10s}  {'Unified':>10s}  {'Observed':>10s}  {'Err_curr':>9s}  {'Err_unif':>9s}")
print("  " + "-" * 65)

for elem in ["V_us", "V_cb", "V_ub"]:
    c_val, c_obs, c_desc = current[elem]
    u_val, u_obs, u_desc = unified[elem]
    err_c = (c_val - c_obs) / c_obs * 100
    err_u = (u_val - u_obs) / u_obs * 100
    better = "UNIFIED" if abs(err_u) < abs(err_c) else "CURRENT"
    print(f"  {elem:>6s}  {c_val:10.6f}  {u_val:10.6f}  {c_obs:10.6f}  {err_c:+8.2f}%  {err_u:+8.2f}%  {better}")
