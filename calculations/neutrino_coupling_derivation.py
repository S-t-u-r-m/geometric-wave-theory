import math

print("="*70)
print("DERIVING g(m) FROM WAVE OVERLAP INTEGRALS")
print("="*70)

print("""
The coupling g(m_CL) in the 3rd-order process nu -> CL -> p -> CL -> nu
comes from the OVERLAP of three standing waves on the lattice:
  1. The neutrino wave (very long wavelength, ~10 um)
  2. The charged lepton wave (wavelength depends on mass)
  3. The proton wave (localized breather, ~0.84 fm)

The coupling vertex is the integral over space of the product
of the wave amplitudes at each lattice site.
""")

# Physical constants
hbar_c = 197.327  # MeV*fm (hbar*c)
m_e = 0.51100     # MeV
m_mu = 105.658    # MeV
m_tau = 1776.86   # MeV
m_p = 938.272     # MeV
m_nu = 0.05e-3    # MeV (~50 meV)

# Compton wavelengths (characteristic size of each wave)
lambda_e = hbar_c / m_e      # electron
lambda_mu = hbar_c / m_mu    # muon
lambda_tau = hbar_c / m_tau  # tau
lambda_p = hbar_c / m_p      # proton
lambda_nu = hbar_c / m_nu    # neutrino

print("Compton wavelengths (characteristic wave sizes):")
print(f"  Neutrino:  {lambda_nu:.2e} fm = {lambda_nu*1e-9:.2f} mm")
print(f"  Electron:  {lambda_e:.2f} fm = {lambda_e*1e-3:.2f} pm")
print(f"  Muon:      {lambda_mu:.4f} fm")
print(f"  Tau:       {lambda_tau:.4f} fm")
print(f"  Proton:    {lambda_p:.4f} fm")

print(f"""
Size hierarchy:
  lambda_nu >> lambda_e >> lambda_p > lambda_mu > lambda_tau

The neutrino wave is essentially constant over the proton's extent.
The electron wave is ~1800x larger than the proton.
The muon wave is comparable to the proton (lambda_mu/lambda_p = {lambda_mu/lambda_p:.2f}).
The tau wave is smaller than the proton (lambda_tau/lambda_p = {lambda_tau/lambda_p:.4f}).
""")

print("="*70)
print("THE WAVE OVERLAP INTEGRAL")
print("="*70)

print("""
The coupling vertex at a lattice site is proportional to:

  V ~ integral psi_nu(x) * psi_CL(x) * psi_p(x) dx

Each wave is a standing wave solution to the GWT Hamiltonian:
  psi(x) ~ A * exp(-x/lambda) * cos(k*x)

where lambda is the Compton wavelength (localization scale)
and k is the wave vector.

For the neutrino: psi_nu is essentially constant (lambda_nu >> everything)
  -> psi_nu ~ A_nu (constant over the interaction region)

For the proton: psi_p is localized within ~lambda_p ~ 0.21 fm
  -> psi_p ~ A_p * exp(-|x|/lambda_p)

The coupling reduces to:
  V ~ A_nu * integral psi_CL(x) * psi_p(x) dx

The key integral is:
  I(m_CL) = integral psi_CL(x) * psi_p(x) dx
""")

print("="*70)
print("EVALUATING THE OVERLAP INTEGRAL")
print("="*70)

print("""
For exponentially localized waves:
  psi_CL(x) ~ exp(-|x| / lambda_CL)
  psi_p(x)  ~ exp(-|x| / lambda_p)

The overlap integral is:
  I = integral_{-inf}^{inf} exp(-|x|/lambda_CL) * exp(-|x|/lambda_p) dx
    = integral_{-inf}^{inf} exp(-|x| * (1/lambda_CL + 1/lambda_p)) dx
    = 2 / (1/lambda_CL + 1/lambda_p)
    = 2 * lambda_CL * lambda_p / (lambda_CL + lambda_p)

This is the HARMONIC MEAN of the two wavelengths!

  I(m_CL) = 2 * lambda_CL * lambda_p / (lambda_CL + lambda_p)

In terms of masses (lambda = hbar_c / m):
  I(m_CL) = 2 * (hbar_c/m_CL) * (hbar_c/m_p) / (hbar_c/m_CL + hbar_c/m_p)
           = 2 * hbar_c / (m_CL + m_p)

So: g(m_CL) ~ 1 / (m_CL + m_p)
""")

# Compute g(m) = 1/(m + m_p)
m_CL = [m_e, m_mu, m_tau]
labels = ["electron", "muon", "tau"]

g_harmonic = [1.0 / (m + m_p) for m in m_CL]
g_h_norm = max(g_harmonic)
g_h = [g / g_h_norm for g in g_harmonic]

print("g(m) = 1/(m + m_p) [from harmonic mean overlap]:\n")
print(f"  {'Lepton':<12} {'m (MeV)':>10} {'g(m)':>12} {'g/g_max':>10} {'g/g_e':>10}")
print(f"  {'-'*55}")
for i, (lab, m) in enumerate(zip(labels, m_CL)):
    print(f"  {lab:<12} {m:10.3f} {g_harmonic[i]:12.6e} {g_h[i]:10.6f} {g_harmonic[i]/g_harmonic[0]:10.4f}")

print(f"\n  g_mu/g_e = {g_harmonic[1]/g_harmonic[0]:.4f}")
print(f"  g_tau/g_e = {g_harmonic[2]/g_harmonic[0]:.4f}")
print(f"  g_tau/g_mu = {g_harmonic[2]/g_harmonic[1]:.4f}")

print(f"""
PROBLEM: g = 1/(m+m_p) gives g_mu/g_e = {g_harmonic[1]/g_harmonic[0]:.4f}
and g_tau/g_e = {g_harmonic[2]/g_harmonic[0]:.4f}

All three are nearly equal! The electron's g is only slightly larger
because m_e << m_p, so 1/(m_e + m_p) ~ 1/m_p for all.

This is TOO democratic — essentially gives back TBM.
""")

print("="*70)
print("WAIT — THE OVERLAP INCLUDES NORMALIZATION")
print("="*70)

print("""
The wave amplitude A depends on the mass too!
A normalized wave has: integral |psi|^2 dx = 1

For psi ~ A * exp(-|x|/lambda):
  A^2 * integral exp(-2|x|/lambda) dx = A^2 * lambda = 1
  -> A = 1/sqrt(lambda) = sqrt(m / hbar_c)

The coupling vertex includes the AMPLITUDES:
  V ~ A_CL * A_p * I
    = sqrt(m_CL/hbar_c) * sqrt(m_p/hbar_c) * 2*hbar_c/(m_CL + m_p)
    = 2 * sqrt(m_CL * m_p) / (m_CL + m_p)

So: g(m_CL) = 2 * sqrt(m_CL * m_p) / (m_CL + m_p)

This is the GEOMETRIC MEAN divided by the ARITHMETIC MEAN!
It equals 1 when m_CL = m_p (perfect matching) and falls off
for both lighter and heavier leptons.

Note: 2*sqrt(a*b)/(a+b) is exactly the ratio of geometric to
arithmetic mean. Maximum at a = b, always <= 1.
""")

# Compute g(m) = 2*sqrt(m*m_p)/(m+m_p)
g_overlap = [2*math.sqrt(m * m_p) / (m + m_p) for m in m_CL]

print("g(m) = 2*sqrt(m*m_p) / (m+m_p)  [normalized overlap integral]:\n")
print(f"  {'Lepton':<12} {'m (MeV)':>10} {'m/m_p':>10} {'g(m)':>10} {'g/g_max':>10}")
print(f"  {'-'*55}")
for i, (lab, m) in enumerate(zip(labels, m_CL)):
    print(f"  {lab:<12} {m:10.3f} {m/m_p:10.5f} {g_overlap[i]:10.6f} {g_overlap[i]/max(g_overlap):10.6f}")

print(f"\nRatios:")
print(f"  g_mu/g_e   = {g_overlap[1]/g_overlap[0]:.4f}")
print(f"  g_tau/g_e  = {g_overlap[2]/g_overlap[0]:.4f}")
print(f"  g_tau/g_mu = {g_overlap[2]/g_overlap[1]:.4f}")
print(f"  g_tau/g_mu mass ratio was: {m_tau/m_mu:.1f}")

print(f"""
NOW we have real differentiation:
  g_e   = {g_overlap[0]:.6f}  (tiny — electron is 1800x smaller than proton)
  g_mu  = {g_overlap[1]:.6f}  (moderate — muon is ~9x smaller than proton)
  g_tau = {g_overlap[2]:.6f}  (large — tau is ~2x the proton, good overlap)

The mu-tau ratio is {g_overlap[2]/g_overlap[1]:.3f} — close to equal!
(Mass ratio was {m_tau/m_mu:.1f}, coupling ratio is only {g_overlap[2]/g_overlap[1]:.1f})

The proton resonance at m = m_p compresses the hierarchy:
  g(m_p) = 2*sqrt(m_p^2)/(2*m_p) = 1.000 (perfect overlap)

Geometric mean / arithmetic mean:
  This ratio measures how "matched" the two waves are.
  Equal waves: ratio = 1 (perfect)
  Very different: ratio -> 0
""")

# Now compute the mass matrix and mixing angles
print("="*70)
print("MASS MATRIX FROM DERIVED g(m)")
print("="*70)

# M_ij = g_i * G_p * g_j where G_p is the proton propagator
# The proton propagator connects axes through its 24 breather modes

# For different proton propagator structures:
print("M_ij = g_i * P_ij * g_j, varying proton propagator structure:\n")

def find_eigvec(mat, lam):
    B = [[mat[i][j] - (lam if i==j else 0) for j in range(3)] for i in range(3)]
    best_v = [1, 0, 0]
    best_norm = 0
    for r1 in range(3):
        for r2 in range(r1+1, 3):
            v = [
                B[r1][1]*B[r2][2] - B[r1][2]*B[r2][1],
                B[r1][2]*B[r2][0] - B[r1][0]*B[r2][2],
                B[r1][0]*B[r2][1] - B[r1][1]*B[r2][0]
            ]
            norm = math.sqrt(sum(x**2 for x in v))
            if norm > best_norm:
                best_norm = norm
                best_v = v
    if best_norm < 1e-15:
        return [1, 0, 0]
    return [x/best_norm for x in best_v]

def compute_angles(A):
    tr = A[0][0] + A[1][1] + A[2][2]
    q = (A[0][0]*A[1][1] - A[0][1]**2) + \
        (A[0][0]*A[2][2] - A[0][2]**2) + \
        (A[1][1]*A[2][2] - A[1][2]**2)
    det = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) - \
          A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]) + \
          A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0])

    p = tr/3
    qq_v = (tr**2 - 3*q) / 9
    r_v = (2*tr**3 - 9*tr*q + 27*det) / 54

    eigenvalues = [p, p, p]
    if qq_v > 1e-30:
        disc = r_v**2 - qq_v**3
        if disc < 0:
            theta_eig = math.acos(max(-1, min(1, r_v / qq_v**1.5)))
            sq = math.sqrt(qq_v)
            eigenvalues = sorted([
                p + 2*sq*math.cos(theta_eig/3),
                p + 2*sq*math.cos((theta_eig - 2*math.pi)/3),
                p + 2*sq*math.cos((theta_eig + 2*math.pi)/3)
            ])

    eigvecs = [find_eigvec(A, lam) for lam in eigenvalues]
    U = [[eigvecs[j][i] for j in range(3)] for i in range(3)]

    s13_2 = abs(U[0][2])**2
    s12_2 = abs(U[0][1])**2 / (1 - s13_2) if s13_2 < 0.999 else 0
    s23_2 = abs(U[1][2])**2 / (1 - s13_2) if s13_2 < 0.999 else 0

    t13 = math.degrees(math.asin(math.sqrt(max(0, min(1, s13_2)))))
    t23 = math.degrees(math.asin(math.sqrt(max(0, min(1, s23_2)))))
    t12 = math.degrees(math.asin(math.sqrt(max(0, min(1, s12_2)))))

    return s12_2, s13_2, s23_2, t12, t13, t23, eigenvalues

g = g_overlap  # Use the derived coupling

# Model 1: Pure democratic proton propagator P_ij = 1 for all i,j
print("--- Model 1: Democratic proton propagator P_ij = 1 ---")
A = [[g[i]*1.0*g[j] for j in range(3)] for i in range(3)]
print("  Matrix (normalized):")
norm = A[2][2]  # normalize to largest
for i in range(3):
    row = [A[i][j]/norm for j in range(3)]
    print(f"    [{row[0]:8.5f} {row[1]:8.5f} {row[2]:8.5f}]")
s12, s13, s23, t12, t13, t23, eigs = compute_angles(A)
print(f"  Rank-1 matrix: only one nonzero eigenvalue")
print(f"  Eigenvector ~ (g_e, g_mu, g_tau) = ({g[0]:.4f}, {g[1]:.4f}, {g[2]:.4f})")
print(f"  This is tau-dominated, NOT democratic")

# Model 2: Democratic + identity proton propagator
# P_ij = D*delta_ij + C  (democratic + diagonal)
print(f"\n--- Model 2: P_ij = D*delta_ij + C (mixed propagator) ---")
print(f"  D controls diagonal preference, C controls cross-axis coupling")
print(f"  r = C/D = cross/diagonal ratio")
print()

print(f"  {'r':>6}  {'s12^2':>7} {'s13^2':>8} {'s23^2':>7} {'t12':>6} {'t13':>6} {'t23':>6}")
print(f"  {'-'*50}")

for r in [0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0]:
    D = 1.0
    C = r * D
    A = [[g[i]*(D*(1 if i==j else 0) + C)*g[j] for j in range(3)] for i in range(3)]
    s12, s13, s23, t12, t13, t23, eigs = compute_angles(A)
    marker = ""
    if abs(t13 - 8.57) < 2:
        marker = " <-- t13!"
    if abs(t13 - 8.57) < 2 and abs(t23 - 49.1) < 5:
        marker = " <== CLOSE!"
    print(f"  {r:6.3f}  {s12:7.4f} {s13:8.5f} {s23:7.4f} {t12:6.1f} {t13:6.1f} {t23:6.1f}{marker}")

print(f"  {'OBS':>6}  {'0.304':>7} {'0.02219':>8} {'0.573':>7} {'33.4':>6} {'8.57':>6} {'49.1':>6}")

# The r->infinity limit is pure democratic (TBM)
# The r->0 limit is pure diagonal (no mixing)
# Somewhere in between should give the right angles!

print(f"""
KEY INSIGHT:
  r -> 0 (no cross-axis coupling): diagonal matrix, no mixing
  r -> inf (all cross-axis): approaches democratic -> TBM

  At intermediate r, the matrix interpolates between:
  - Pure axis-dependent (diagonal, from g_i^2)
  - Pure democratic (off-diagonal, from g_i * g_j)

  The mixing angles depend on the RATIO r = C/D.

  In GWT, this ratio comes from the proton's internal structure:
  - D = on-axis coupling (8 breather modes per axis)
  - C = cross-axis coupling (through the proton's cross-axis modes)
  - r = C/D is determined by the proton's mode structure
""")

# What r gives the best fit?
print("="*70)
print("FINDING BEST r (SINGLE PARAMETER FIT)")
print("="*70)

best_score = 999
best_r = 0
best_result = None

for i_r in range(1, 10000):
    r = i_r * 0.001
    D = 1.0
    C = r * D
    A = [[g[i]*(D*(1 if i==j else 0) + C)*g[j] for j in range(3)] for i in range(3)]
    s12, s13, s23, t12, t13, t23, eigs = compute_angles(A)
    score = ((s12 - 0.304)/0.304)**2 + ((s13 - 0.02219)/0.02219)**2 + ((s23 - 0.573)/0.573)**2
    if score < best_score:
        best_score = score
        best_r = r
        best_result = (s12, s13, s23, t12, t13, t23)

print(f"Best fit: r = C/D = {best_r:.4f}")
print(f"  sin2(t12) = {best_result[0]:.4f}  (obs: 0.304)")
print(f"  sin2(t13) = {best_result[1]:.5f}  (obs: 0.02219)")
print(f"  sin2(t23) = {best_result[2]:.4f}  (obs: 0.573)")
print(f"  t12 = {best_result[3]:.2f}  t13 = {best_result[4]:.2f}  t23 = {best_result[5]:.2f}")
print(f"  Score = {best_score:.6f}")

# What does this r correspond to physically?
print(f"""
Physical interpretation:
  r = C/D = {best_r:.4f}

  If the proton has 8 modes per axis:
    D = 8 (on-axis modes)
    C = r * D = {best_r * 8:.2f} (effective cross-axis modes)

  Cross-axis coupling epsilon = C/(2*D) = {best_r/2:.4f}
  (factor of 2 because there are 2 other axes)
""")

# Also try with the full 24-mode proton propagator
print("="*70)
print("WITH PROTON MODE STRUCTURE: P_ii = 8+16*eps^2, P_ij = 16*eps+8*eps^2")
print("="*70)

best_score = 999
best_eps = 0
best_result = None

for i_eps in range(1, 1000):
    eps = i_eps * 0.001
    P_diag = 8*(1 + 2*eps**2)
    P_off = 16*eps + 8*eps**2

    A = [[g[i]*(P_diag if i==j else P_off)*g[j] for j in range(3)] for i in range(3)]
    s12, s13, s23, t12, t13, t23, eigs = compute_angles(A)
    score = ((s12 - 0.304)/0.304)**2 + ((s13 - 0.02219)/0.02219)**2 + ((s23 - 0.573)/0.573)**2
    if score < best_score:
        best_score = score
        best_eps = eps
        best_result = (s12, s13, s23, t12, t13, t23)

print(f"Best fit: eps = {best_eps:.4f}")
P_d = 8*(1 + 2*best_eps**2)
P_o = 16*best_eps + 8*best_eps**2
print(f"  P_diag = {P_d:.3f}, P_off = {P_o:.3f}, ratio = {P_o/P_d:.4f}")
print(f"  sin2(t12) = {best_result[0]:.4f}  (obs: 0.304)")
print(f"  sin2(t13) = {best_result[1]:.5f}  (obs: 0.02219)")
print(f"  sin2(t23) = {best_result[2]:.4f}  (obs: 0.573)")
print(f"  t12 = {best_result[3]:.2f}  t13 = {best_result[4]:.2f}  t23 = {best_result[5]:.2f}")
print(f"  Score = {best_score:.6f}")


print(f"\n{'='*70}")
print("SUMMARY: DERIVED COUPLING FUNCTION")
print(f"{'='*70}")

print(f"""
The wave overlap integral gives:

  g(m) = 2 * sqrt(m * m_p) / (m + m_p)

This is the ratio of geometric mean to arithmetic mean of the
charged lepton and proton masses. It measures how well the two
standing waves "match" in spatial extent.

  g(m_e)   = {g_overlap[0]:.6f}  (electron: poorly matched, 1800x smaller)
  g(m_mu)  = {g_overlap[1]:.6f}  (muon: moderately matched, 9x smaller)
  g(m_tau) = {g_overlap[2]:.6f}  (tau: well matched, 2x larger)
  g(m_p)   = 1.000000  (proton: perfect self-matching)

Key ratios:
  g_tau/g_mu = {g_overlap[2]/g_overlap[1]:.4f} (vs mass ratio {m_tau/m_mu:.1f})
  g_mu/g_e   = {g_overlap[1]/g_overlap[0]:.2f}
  g_tau/g_e  = {g_overlap[2]/g_overlap[0]:.2f}

The coupling function is DERIVED, not guessed:
  - Exponential localization of standing waves
  - Overlap integral of two exponentials = harmonic mean
  - Normalization of wave amplitudes = sqrt(m)
  - Combined: g = geometric_mean / arithmetic_mean
""")
