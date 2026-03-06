"""
WHY sigma_p/sigma_tau AND NOT sigma_tau/sigma_p?
Deep investigation of the tau wrapping ratio geometry.
"""

import numpy as np
from scipy.spatial.transform import Rotation as Rot

m_e, m_mu, m_tau, m_p = 0.51100, 105.658, 1776.86, 938.272
d = 3
sigmas = np.array([m_e, m_mu, m_tau])**(-1./d)
s_p = m_p**(-1./d)

OBS = (33.41, 8.57, 49.20)
ERRS = (0.78, 0.12, 1.05)
U_tbm = np.array([
    [ np.sqrt(2./3),  1/np.sqrt(3),  0],
    [-1/np.sqrt(6),   1/np.sqrt(3),  1/np.sqrt(2)],
    [ 1/np.sqrt(6),  -1/np.sqrt(3),  1/np.sqrt(2)]
])
theta = np.arcsin((m_e/m_mu)**(1./d))

def predict(axis_raw):
    axis = axis_raw / np.linalg.norm(axis_raw)
    rot = Rot.from_rotvec(theta * axis)
    U = rot.as_matrix() @ U_tbm
    s13 = abs(U[0,2])
    t13 = np.degrees(np.arcsin(min(s13,1.0)))
    t23 = np.degrees(np.arctan2(abs(U[1,2]), abs(U[2,2])))
    c13 = np.cos(np.radians(t13))
    s12 = min(abs(U[0,1])/c13, 1.0)
    t12 = np.degrees(np.arcsin(s12))
    return t12, t13, t23

print("="*70)
print("  WHY sigma_p/sigma_tau AND NOT sigma_tau/sigma_p?")
print("="*70)
print()

# Scan the tau component
best_c = 0
best_score = 999
for c_100 in range(50, 200):
    c = c_100 * 0.01
    t12, t13, t23 = predict(np.array([-1, np.sqrt(3), -c]))
    score = sum(((t-o)/e)**2 for t,o,e in zip((t12,t13,t23), OBS, ERRS))
    if score < best_score:
        best_score = score
        best_c = c

print(f"Best-fit c = {best_c:.4f}")
print(f"sigma_p/sigma_tau = {s_p/sigmas[2]:.4f}")
print(f"(m_tau/m_p)^(1/d) = {(m_tau/m_p)**(1./d):.4f}")
print()

# Compare candidates
candidates = {
    "sigma_tau/sigma_p (overlap)": sigmas[2]/s_p,
    "sigma_p/sigma_tau (wrapping)": s_p/sigmas[2],
    "(m_tau/m_p)^(1/d)": (m_tau/m_p)**(1./d),
    "(m_p/m_tau)^(1/d)": (m_p/m_tau)**(1./d),
    "(m_tau/m_p)^(2/d)": (m_tau/m_p)**(2./d),
    "sigma_mu/sigma_tau": sigmas[1]/sigmas[2],
    "1 (pure geometric)": 1.0,
}

print(f"Best c = {best_c:.4f}. Candidate values:")
for name, val in candidates.items():
    err = abs(val/best_c - 1)*100
    marker = " <<<<" if err < 0.5 else (" <<<" if err < 2 else "")
    t12, t13, t23 = predict(np.array([-1, np.sqrt(3), -val]))
    sigs = [(t-o)/e for t,o,e in zip((t12,t13,t23), OBS, ERRS)]
    print(f"  {name:35s} = {val:.4f}  err={err:.2f}%  [{sigs[0]:+.2f}s {sigs[1]:+.2f}s {sigs[2]:+.2f}s]{marker}")

print()
print("="*70)
print("  THE KEY INSIGHT: AXIS = MU DIRECTION * TAU RESCALING")
print("="*70)
print()

mu_dir = np.array([-1, np.sqrt(3), -1])
scale = np.array([1, 1, s_p/sigmas[2]])
axis_test = mu_dir * scale
print("The axis = mu_direction * (1, 1, sigma_p/sigma_tau)")
print(f"  mu direction: ({mu_dir[0]:.4f}, {mu_dir[1]:.4f}, {mu_dir[2]:.4f})")
print(f"  rescaling:    (1, 1, {s_p/sigmas[2]:.4f})")
print(f"  result:       ({axis_test[0]:.4f}, {axis_test[1]:.4f}, {axis_test[2]:.4f})")
print(f"  exact:        (-1, 1.7321, -{(m_tau/m_p)**(1./d):.4f})")
print()

t12, t13, t23 = predict(axis_test)
sigs = [(t-o)/e for t,o,e in zip((t12,t13,t23), OBS, ERRS)]
print(f"  Prediction: t12={t12:.2f}({sigs[0]:+.2f}s) t13={t13:.2f}({sigs[1]:+.2f}s) t23={t23:.2f}({sigs[2]:+.2f}s)")
print()

print("="*70)
print("  PHYSICAL INTERPRETATION")
print("="*70)
print()
print("In the far-field regime (e, mu both extend beyond proton):")
print("  The axis component is GEOMETRIC - determined by the TBM")
print("  equilateral triangle structure alone. The mass ratios")
print("  are irrelevant because both knots are so large that the")
print("  proton sees them as effectively infinite backgrounds.")
print()
print("In the near-field regime (tau inside proton):")
print("  The proton wraps around the tau. The effective coupling")
print("  is ENHANCED by how much the proton extends beyond:")
print(f"  sigma_p/sigma_tau = {s_p/sigmas[2]:.4f}")
print()
print("  This rescales the tau component of the mu direction from")
print(f"  -1 to -1 * {s_p/sigmas[2]:.4f} = {-s_p/sigmas[2]:.4f}")
print()
print("  The physical reason: when the proton wraps the tau, the")
print("  tau-axis coupling is amplified by the proton-to-tau size")
print("  ratio. The tau 'sees' more proton than itself, increasing")
print("  its contribution to the rotation.")
print()

# ================================================================
# DEEPER: What does the (1,1,1) component of the axis mean?
# ================================================================
print("="*70)
print("  THE (1,1,1) COMPONENT AND THETA_13")
print("="*70)
print()

a1 = np.array([1, -1, 0]) / np.sqrt(2)
a2 = np.array([1, 1, -2]) / np.sqrt(6)
e3 = np.array([1, 1, 1]) / np.sqrt(3)

b = (m_tau/m_p)**(1./d)
axis_emp = np.array([-1, np.sqrt(3), -b])
axis_geo = np.array([-1, np.sqrt(3), -1])

c3_geo = np.dot(axis_geo, e3)
c3_emp = np.dot(axis_emp, e3)

print(f"(1,1,1) component of geometric axis: {c3_geo:.4f}")
print(f"(1,1,1) component of empirical axis:  {c3_emp:.4f}")
print(f"Ratio: {c3_emp/c3_geo:.3f} (the tau correction nearly DOUBLES it)")
print()

# What if the (1,1,1) component is proportional to theta_13?
# In perturbation theory, the (1,1,1) component mixes the
# non-degenerate eigenvector with the degenerate subspace.
# The mixing angle is proportional to the (1,1,1) component
# of the perturbation direction times sin(theta_rot).

# For a rotation R(axis, theta) applied to U_TBM:
# The amount of (1,1,1) mixing generated is:
# ~ sin(theta) * cos(angle between axis and degenerate subspace)
# = sin(theta) * sin(tilt angle from degenerate subspace)

# The tilt angle of the axis from the degenerate subspace:
norm_full = np.linalg.norm(axis_emp)
tilt_geo = np.arcsin(abs(c3_geo)/np.linalg.norm(axis_geo))
tilt_emp = np.arcsin(abs(c3_emp)/norm_full)
print(f"Tilt from degenerate subspace:")
print(f"  Geometric: {np.degrees(tilt_geo):.2f} deg")
print(f"  Empirical: {np.degrees(tilt_emp):.2f} deg")
print()

# The effective theta_13 is approximately:
# sin(theta_13) ~ sin(theta_rot) * sin(tilt) * geometric_factor
sin_th = np.sin(theta)
print(f"sin(theta_rot) = {sin_th:.4f}")
print(f"sin(theta_rot) * sin(tilt_emp) = {sin_th * np.sin(tilt_emp):.4f}")
print(f"Observed sin(theta_13) = {np.sin(np.radians(8.57)):.4f}")
print()

# ================================================================
# Now explore: what if we can write the FULL axis as a function
# that works in all regimes?
# ================================================================
print("="*70)
print("  UNIFIED FORMULA FOR THE AXIS?")
print("="*70)
print()

# axis_i = (-1)^f(i) * max(1, sigma_p/sigma_i) * sign_i
# where sign_i comes from the mu direction

# Actually: the mu direction in flavor space is (-1, 2, -1)/sqrt(6)
# (in the degenerate subspace). The FULL mu direction including the
# (1,1,1) component depends on how we embed it.

# The empirical axis is: mu_direction * diag(1, 1, sigma_p/sigma_tau)
# This is an ANISOTROPIC scaling of the mu direction.

# What if we generalize: axis_i = mu_dir_i * max(1, sigma_p/sigma_i)?
# For e: max(1, 0.082) = 1 -> axis_e = -1*1 = -1  (correct!)
# For mu: max(1, 0.483) = 1 -> axis_mu = sqrt(3)*1 = sqrt(3)  (correct!)
# For tau: max(1, 1.237) = 1.237 -> axis_tau = -1*1.237 = -1.237  (correct!)

print("UNIFIED FORMULA: axis_i = mu_dir_i * max(1, sigma_p/sigma_i)")
print()
for i, name in enumerate(["e", "mu", "tau"]):
    ratio = s_p/sigmas[i]
    factor = max(1, ratio)
    mu_comp = [-1, np.sqrt(3), -1][i]
    result = mu_comp * factor
    emp_comp = [-1, np.sqrt(3), -b][i]
    print(f"  {name}: sigma_p/sigma = {ratio:.4f}, max(1,.) = {factor:.4f}, "
          f"mu_dir = {mu_comp:+.4f} -> axis = {result:+.4f} (exact: {emp_comp:+.4f})")

print()
print("THIS IS EXACT! The formula is:")
print()
print("  axis_i = mu_direction_i * max(1, sigma_p / sigma_i)")
print()
print("  where mu_direction = (-1, sqrt(d), -1)")
print("  and sigma_i = m_i^(-1/d), sigma_p = m_p^(-1/d)")
print()
print("Physical meaning: the rotation axis is the mu direction in")
print("the TBM flavor triangle, with each component rescaled by")
print("how much the proton wraps around that lepton. For leptons")
print("larger than the proton (e, mu), the factor is 1 (no wrapping).")
print("For the tau (inside the proton), the factor is sigma_p/sigma_tau > 1.")
print()

# Verify this is exactly the same as the empirical formula
axis_unified = np.array([-1, np.sqrt(3), -1]) * np.array([
    max(1, s_p/sigmas[0]),
    max(1, s_p/sigmas[1]),
    max(1, s_p/sigmas[2])
])
axis_empirical = np.array([-1, np.sqrt(3), -(m_tau/m_p)**(1./d)])
print(f"Unified:   ({axis_unified[0]:.6f}, {axis_unified[1]:.6f}, {axis_unified[2]:.6f})")
print(f"Empirical: ({axis_empirical[0]:.6f}, {axis_empirical[1]:.6f}, {axis_empirical[2]:.6f})")
print(f"Difference: {np.linalg.norm(axis_unified - axis_empirical):.2e}")
