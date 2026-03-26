"""
Angular D_e(θ) → Oh Irrep Decomposition
==========================================
The V8 bond model uses Oh irrep fractions as correction weights:
  A1g (base):  σ coupling strength
  T2g:         π bonds (W_PI = cos(π/d) = 0.5)
  Eg:          LP repulsion ((d²+1)/d³ = 10/27)
  T1g:         Radical corrections (5/6)

The angular D_e(θ) from the torus simulation contains the SAME
information, encoded as an angular function. The Oh harmonics
on the sphere decompose any angular function into irreps:

  D_e(θ) = Σ_Γ c_Γ × χ_Γ(θ)

where χ_Γ(θ) are the Oh-symmetric angular basis functions and
c_Γ are the coefficients. If c_σ/c_π matches W_PI = 0.5, we've
connected the simulation to the bond model.

The Oh irreps and their angular content on the sphere:
  A1g: l=0 (isotropic)
  Eg:  l=2, m=0 and l=2, m=±2 (quadrupolar)
  T2g: l=2, m=±1 and l=2, m=±2 (another quadrupolar combo)
  T1g: l=1 (dipolar, but Oh forbids — appears at l=3)

For axial symmetry (torus around z-axis), only m=0 survives:
  D_e(θ) = Σ_l a_l × P_l(cos θ)

where P_l are Legendre polynomials. The mapping to Oh irreps:
  l=0 → A1g
  l=2 → Eg (the cos²θ part)
  l=4 → A1g + Eg + T2g
  l=6 → A1g + Eg + T1g + T2g
"""
import sys, io, os, time
import numpy as np
from scipy.special import legendre
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.interpolate import interp1d

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "angular_oh_decomposition_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("ANGULAR D_e(θ) → Oh IRREP DECOMPOSITION")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("")

# ============================================================
# RECOMPUTE D_e(θ) — FINE ANGULAR SCAN
# ============================================================
N3 = 64
R_maj = 8
kw = 3

ix = np.arange(N3, dtype=np.float64) - N3/2
X, Y, Z = np.meshgrid(ix, ix, ix, indexing='ij')

def make_torus(X, Y, Z, z_center, R_major, kw):
    rho_xy = np.sqrt(X**2 + Y**2 + 1e-30)
    rho_tube = np.sqrt((rho_xy - R_major)**2 + (Z - z_center)**2)
    return (4.0/PI) * (np.arctan(np.exp(rho_tube + kw/2.0))
                      - np.arctan(np.exp(rho_tube - kw/2.0)))

phi_torus = make_torus(X, Y, Z, 0.0, R_maj, kw)

def extract_slice(phi_3d, direction, N3):
    dx, dy, dz = direction
    norm = np.sqrt(dx**2 + dy**2 + dz**2)
    dx, dy, dz = dx/norm, dy/norm, dz/norm
    r_max = N3 // 2 - 2
    n_samples = N3 * 4
    r_values = np.linspace(-r_max, r_max, n_samples)
    phi_values = np.zeros(n_samples)
    for i, r in enumerate(r_values):
        px, py, pz = r*dx, r*dy, r*dz
        ix_f, iy_f, iz_f = px + N3/2, py + N3/2, pz + N3/2
        ix0 = int(np.floor(ix_f)) % N3
        iy0 = int(np.floor(iy_f)) % N3
        iz0 = int(np.floor(iz_f)) % N3
        ix1, iy1, iz1 = (ix0+1)%N3, (iy0+1)%N3, (iz0+1)%N3
        fx = ix_f - np.floor(ix_f)
        fy = iy_f - np.floor(iy_f)
        fz = iz_f - np.floor(iz_f)
        phi_values[i] = (
            phi_3d[ix0,iy0,iz0]*(1-fx)*(1-fy)*(1-fz) +
            phi_3d[ix1,iy0,iz0]*fx*(1-fy)*(1-fz) +
            phi_3d[ix0,iy1,iz0]*(1-fx)*fy*(1-fz) +
            phi_3d[ix0,iy0,iz1]*(1-fx)*(1-fy)*fz +
            phi_3d[ix1,iy1,iz0]*fx*fy*(1-fz) +
            phi_3d[ix1,iy0,iz1]*fx*(1-fy)*fz +
            phi_3d[ix0,iy1,iz1]*(1-fx)*fy*fz +
            phi_3d[ix1,iy1,iz1]*fx*fy*fz)
    return r_values, phi_values

N1d = 512

def build_1d_hessian(phi_1d, N):
    diag = 2.0 + np.cos(PI * phi_1d)
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()

def get_De(r_vals, phi_vals, R_scan):
    x_1d = np.arange(N1d, dtype=np.float64) - N1d/2
    center_1d = N1d // 2
    f_interp = interp1d(r_vals, phi_vals, kind='cubic', fill_value=0.0, bounds_error=False)
    phi_single = f_interp(x_1d - center_1d)
    H_single = build_1d_hessian(phi_single, N1d)
    ev_single, _ = eigsh(H_single, k=5, which='SM')
    E0 = np.sort(ev_single)[0]
    V_min = 0
    for R in R_scan:
        phi_A = f_interp(x_1d - center_1d + R/2)
        phi_B = f_interp(x_1d - center_1d - R/2)
        H_d = build_1d_hessian(phi_A + phi_B, N1d)
        try:
            ev_d, _ = eigsh(H_d, k=3, which='SM')
            V = np.sort(ev_d)[0] - E0
            if V < V_min:
                V_min = V
        except Exception:
            pass
    return -V_min  # D_e (positive)

# Fine angular scan
report("Computing D_e(θ) with fine angular resolution...")
theta_degrees = list(range(0, 91, 3))  # every 3 degrees
R_scan = list(range(4, 30, 2))

De_theta = {}
for theta_deg in theta_degrees:
    theta = theta_deg * PI / 180
    if theta_deg == 0:
        direction = (0, 0, 1)
    elif theta_deg == 90:
        direction = (1, 0, 0)
    else:
        direction = (np.sin(theta), 0, np.cos(theta))

    r_vals, phi_vals = extract_slice(phi_torus, direction, N3)
    De = get_De(r_vals, phi_vals, R_scan)
    De_theta[theta_deg] = De

report("Done.")
report("")

# ============================================================
# PART 1: LEGENDRE POLYNOMIAL DECOMPOSITION
# ============================================================
report("PART 1: LEGENDRE POLYNOMIAL DECOMPOSITION")
report("-" * 55)
report("D_e(θ) = Σ_l a_l × P_l(cos θ)")
report("")

theta_rad = np.array([t * PI / 180 for t in theta_degrees])
cos_theta = np.cos(theta_rad)
De_arr = np.array([De_theta[t] for t in theta_degrees])

# Fit Legendre coefficients up to l=10
l_max = 10
a_l = np.zeros(l_max + 1)

for l in range(l_max + 1):
    Pl = legendre(l)
    Pl_values = Pl(cos_theta)
    # Coefficient: a_l = (2l+1)/2 × ∫ D_e(θ) × P_l(cosθ) × sinθ dθ
    integrand = De_arr * Pl_values * np.sin(theta_rad)
    a_l[l] = (2*l + 1) / 2.0 * np.trapezoid(integrand, theta_rad)

report(f"{'l':>3} {'a_l':>14} {'a_l/a_0':>10} {'Oh irrep':>15}")
report("-" * 45)

# Oh irrep assignments for each l
oh_map = {
    0: 'A1g',
    1: '(forbidden)',
    2: 'Eg',
    3: '(T1g+T2g)',
    4: 'A1g+Eg+T2g',
    5: '(T1g+T2g)',
    6: 'A1g+Eg+T1g+T2g',
    7: '(A2g+T1g+T2g)',
    8: 'A1g+Eg+T2g+...',
    9: '(...)',
    10: 'A1g+Eg+...'
}

for l in range(l_max + 1):
    ratio = a_l[l] / a_l[0] if abs(a_l[0]) > 1e-12 else 0
    oh = oh_map.get(l, '...')
    report(f"  {l:3d} {a_l[l]:14.8f} {ratio:+10.4f} {oh:>15}")

report("")

# Quality check: reconstruct D_e(θ) from the fit
theta_check = np.linspace(0, PI/2, 100)
De_reconstructed = np.zeros_like(theta_check)
for l in range(l_max + 1):
    Pl = legendre(l)
    De_reconstructed += a_l[l] * Pl(np.cos(theta_check))

# Compare at measured points
report("Fit quality:")
De_fit_at_measured = np.zeros(len(theta_degrees))
for l in range(l_max + 1):
    Pl = legendre(l)
    De_fit_at_measured += a_l[l] * Pl(cos_theta)

residual = np.sqrt(np.mean((De_fit_at_measured - De_arr)**2))
rms_De = np.sqrt(np.mean(De_arr**2))
report(f"  RMS residual: {residual:.6f}")
report(f"  RMS D_e: {rms_De:.6f}")
report(f"  Fit quality: {(1 - residual/rms_De)*100:.1f}%")
report("")

# ============================================================
# PART 2: MAP TO Oh IRREPS
# ============================================================
report("PART 2: Oh IRREP CONTENT")
report("-" * 55)
report("")
report("On the sphere with Oh symmetry, the angular function D_e(θ,φ)")
report("decomposes into irreps. For axial symmetry (m=0 only):")
report("")
report("  A1g content: l=0, 4, 6, 8, 10, ...")
report("  Eg content:  l=2, 4, 6, 8, 10, ...")
report("  T2g content: l=4, 6, 8, ...")
report("  T1g content: l=6, 8, ...")
report("")

# A1g: isotropic part = a_0 (l=0) + contributions from l=4,6,...
# Eg: anisotropic quadrupolar = a_2 (l=2) + contributions from l=4,6,...
# The l=2 component is purely Eg on the Oh sphere.

A1g_strength = a_l[0]  # leading A1g contribution
Eg_strength = a_l[2]   # leading Eg contribution (l=2 → pure Eg in Oh)

report(f"Leading contributions:")
report(f"  A1g (l=0):  {A1g_strength:.8f} — isotropic bond strength")
report(f"  Eg  (l=2):  {Eg_strength:.8f} — quadrupolar anisotropy")
report(f"  l=4:        {a_l[4]:.8f} — A1g + Eg + T2g mix")
report(f"  l=6:        {a_l[6]:.8f} — A1g + Eg + T1g + T2g mix")
report("")

# ============================================================
# PART 3: COMPARE WITH V8 BOND MODEL FRACTIONS
# ============================================================
report("PART 3: COMPARISON WITH V8 BOND MODEL")
report("-" * 55)
report("")

# V8 uses these Oh-derived constants:
# C_BOND = π/d² = 0.349 (A1g base coupling)
# W_PI = cos(π/d) = 0.5 (pi bond weight)
# LP_I = (d²+1)/d³ = 0.370 (LP repulsion)
# F_RAD = (2d-1)/(2d) = 5/6 (radical reduction)

C_BOND = PI / d**2
W_PI = np.cos(PI / d)
LP_I = (d**2 + 1) / d**3
F_RAD = (2*d - 1) / (2*d)

report("V8 constants from Oh group theory:")
report(f"  C_BOND = π/d² = {C_BOND:.6f}")
report(f"  W_PI   = cos(π/d) = {W_PI:.6f}  (π bond / σ bond ratio)")
report(f"  LP_I   = (d²+1)/d³ = {LP_I:.6f}")
report(f"  F_RAD  = (2d-1)/(2d) = {F_RAD:.6f}")
report("")

# The σ/π ratio from the simulation
# σ bond ≈ D_e at θ=75° (strongest single-peak)
# π bond ≈ D_e at θ=55-60° (oblique tube crossing)

De_sigma = De_theta.get(75, 0)
De_pi_60 = De_theta.get(60, 0)
De_pi_57 = De_theta.get(57, De_theta.get(54, De_theta.get(55, 0)))

# But more precisely: the σ/π ratio should come from the
# angular decomposition. In the bond model:
#   D_σ = C_BOND × E_harm × 1 (σ contribution)
#   D_π = C_BOND × E_harm × W_PI (π contribution)
#   So W_PI = D_π / D_σ = cos(π/d) = 0.5

# From the Legendre decomposition:
# The "σ part" is the isotropic component (A1g, l=0)
# The "π part" is the quadrupolar component (Eg, l=2)
# The ratio Eg/A1g should relate to W_PI

if abs(A1g_strength) > 1e-10:
    ratio_Eg_A1g = abs(Eg_strength) / A1g_strength
    report(f"From angular decomposition:")
    report(f"  A1g (isotropic):   {A1g_strength:.6f}")
    report(f"  Eg (quadrupolar):  {Eg_strength:.6f}")
    report(f"  |Eg/A1g| ratio:    {ratio_Eg_A1g:.4f}")
    report(f"  V8 W_PI = cos(π/d): {W_PI:.4f}")
    report("")

# Alternative: use the angular integral directly
# σ-type: integral of D_e(θ)sin(θ) for θ > 60° (near-equatorial)
# π-type: integral of D_e(θ)sin(θ) for 40° < θ < 60° (oblique)

theta_rad_fine = np.array([t * PI / 180 for t in theta_degrees])
sin_t = np.sin(theta_rad_fine)
De_fine = np.array([De_theta[t] for t in theta_degrees])

# Define angular regions
sigma_mask = np.array([t >= 65 for t in theta_degrees])
pi_mask = np.array([(t >= 45 and t < 65) for t in theta_degrees])
delta_mask = np.array([(t >= 30 and t < 45) for t in theta_degrees])

D_sigma_int = np.trapezoid((De_fine * sin_t)[sigma_mask], theta_rad_fine[sigma_mask]) if np.sum(sigma_mask) > 1 else 0
D_pi_int = np.trapezoid((De_fine * sin_t)[pi_mask], theta_rad_fine[pi_mask]) if np.sum(pi_mask) > 1 else 0
D_delta_int = np.trapezoid((De_fine * sin_t)[delta_mask], theta_rad_fine[delta_mask]) if np.sum(delta_mask) > 1 else 0
D_total_int = np.trapezoid(De_fine * sin_t, theta_rad_fine)

report("Angular-integrated bond strengths:")
report(f"  σ region (θ≥65°): {D_sigma_int:.6f} ({D_sigma_int/D_total_int*100:.1f}%)")
report(f"  π region (45-65°): {D_pi_int:.6f} ({D_pi_int/D_total_int*100:.1f}%)")
report(f"  δ region (30-45°): {D_delta_int:.6f} ({D_delta_int/D_total_int*100:.1f}%)")
report(f"  Total: {D_total_int:.6f}")
report("")

if D_sigma_int > 0:
    sim_pi_sigma = D_pi_int / D_sigma_int
    report(f"  π/σ ratio from simulation: {sim_pi_sigma:.4f}")
    report(f"  V8 W_PI = cos(π/d):        {W_PI:.4f}")
    report(f"  Match: {abs(sim_pi_sigma - W_PI)/W_PI*100:.1f}% error")
    report("")

# ============================================================
# PART 4: THE CRITICAL QUESTION — DOES A1g FRACTION MATCH C_BOND?
# ============================================================
report("PART 4: DOES THE A1g FRACTION MATCH C_BOND = π/d²?")
report("-" * 55)
report("")

# The A1g fraction of the angular function D_e(θ):
# A1g_frac = a_0 / max(D_e) or a_0 / ∫D_e×sinθ dθ

A1g_frac_of_total = A1g_strength / D_total_int if D_total_int > 0 else 0
report(f"  a_0 / ∫D_e sinθ dθ = {A1g_frac_of_total:.6f}")
report(f"  C_BOND = π/d² = {C_BOND:.6f}")
report(f"  π/d² / (a_0/∫) = {C_BOND/A1g_frac_of_total:.4f}" if A1g_frac_of_total > 0 else "  (zero)")
report("")

# What fraction of the PEAK D_e is the isotropic part?
De_max = max(De_fine)
report(f"  a_0 / D_e_max = {A1g_strength / De_max:.6f}" if De_max > 0 else "  (zero)")
report(f"  C_BOND = {C_BOND:.6f}")
report("")

# The ratio of σ integral to total should be the A1g/T1u² fraction
# In Oh: T1u ⊗ T1u = A1g + Eg + T1g + T2g
# The A1g fraction of T1u² = 1/9 on the Oh group
# But the A1g fraction of the 3D solid angle might be different
A1g_of_T1u2 = 1 / d**2  # 1/9
report(f"  Oh prediction: A1g fraction of T1u⊗T1u = 1/d² = {A1g_of_T1u2:.6f}")
report(f"  C_BOND = π/d² = π × (1/d²) = π × A1g_fraction")
report(f"  The π factor comes from the potential periodicity!")
report("")

# ============================================================
# PART 5: SIGMA BOND = A1g, PI BOND = T2g?
# ============================================================
report("PART 5: BOND TYPE ↔ Oh IRREP MAPPING")
report("-" * 55)
report("")
report("Hypothesis:")
report("  σ bond → A1g channel (isotropic, strongest)")
report("  π bond → Eg or T2g channel (quadrupolar, weaker)")
report("  δ bond → higher-l channels (octupolar, weakest)")
report("")
report("In V8:")
report("  coupling = σ_eff + π_eff × W_PI")
report("         = A1g_part + T2g_part × cos(π/d)")
report("")
report("From our simulation:")
report(f"  Isotropic (A1g) strength:    a_0 = {A1g_strength:.6f}")
report(f"  Quadrupolar (Eg) strength:  |a_2| = {abs(Eg_strength):.6f}")
report(f"  Octupolar (l=3) strength:   |a_3| = {abs(a_l[3]):.6f}")
report(f"  l=4 strength:               |a_4| = {abs(a_l[4]):.6f}")
report("")

# The coupling formula in V8 is:
# D_e = (π/d²) × E_harm × (1 + n_pi × cos(π/d))
# The (1 + n_pi × W_PI) factor encodes the angular structure
# For a single bond (n_pi=0): D_e = (π/d²) × E_harm × 1
# For a double bond (n_pi=1): D_e = (π/d²) × E_harm × 1.5
# For a triple bond (n_pi=2): D_e = (π/d²) × E_harm × 2.0

# Does our D_e(θ) have the same ratio structure?
# At the σ peak (75°): D_e ≈ 0.46
# Add one π bond (this means adding contribution from 55-60° angular range)
# The "double bond" would be σ + π channel

De_at_75 = De_theta.get(75, 0)
De_at_60 = De_theta.get(60, 0)

if De_at_75 > 0:
    # Double bond ≈ σ + π
    # The ratio should be (1 + W_PI) = 1.5
    double_single_ratio = (De_at_75 + De_at_60) / De_at_75
    report(f"Single bond (θ=75°): D_e = {De_at_75:.6f}")
    report(f"Adding π channel (θ=60°): D_e = {De_at_60:.6f}")
    report(f"Double/single ratio: {double_single_ratio:.4f}")
    report(f"V8 prediction (1+W_PI): {1+W_PI:.4f}")
    report(f"Match: {abs(double_single_ratio - (1+W_PI))/(1+W_PI)*100:.1f}%")
    report("")

    # Triple bond = σ + 2π
    triple_single = (De_at_75 + 2 * De_at_60) / De_at_75
    report(f"Triple/single ratio: {triple_single:.4f}")
    report(f"V8 prediction (1+2×W_PI): {1+2*W_PI:.4f}")
    report(f"Match: {abs(triple_single - (1+2*W_PI))/(1+2*W_PI)*100:.1f}%")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report("Legendre decomposition of D_e(θ) from torus simulation:")
for l in range(7):
    report(f"  l={l}: a_l = {a_l[l]:+.6f} (a_l/a_0 = {a_l[l]/a_l[0]:+.4f})")
report("")

report("KEY COMPARISONS WITH V8:")
report(f"  π/σ angular integral ratio: {sim_pi_sigma:.4f} vs V8 W_PI = {W_PI:.4f}")
if D_sigma_int > 0:
    report(f"  Error: {abs(sim_pi_sigma - W_PI)/W_PI*100:.1f}%")
report("")

report(f"  The isotropic component a_0 captures the σ-bond strength.")
report(f"  The l=2 (Eg) component captures the π-bond anisotropy.")
report(f"  Higher l encode δ bonds and LP repulsion geometry.")
report("")

if abs(sim_pi_sigma - W_PI)/W_PI < 0.3:
    report("  *** π/σ RATIO MATCHES V8 within 30% ***")
    report("  The torus geometry PRODUCES the cos(π/d) weight.")
elif abs(sim_pi_sigma - W_PI)/W_PI < 0.5:
    report("  π/σ ratio is in the right ballpark but needs refinement.")
else:
    report("  π/σ ratio does not match — the angular regions may need")
    report("  different boundaries, or the mapping is more subtle.")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
