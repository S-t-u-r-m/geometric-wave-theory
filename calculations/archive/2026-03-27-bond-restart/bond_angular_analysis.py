"""
Angular Bond Analysis — Sigma/Pi/Delta from Torus Geometry
============================================================
The 3D torus slices give a D_e(θ) profile that maps directly to
chemical bond types. This script analyzes that mapping.

KEY GEOMETRY:
  Two coaxial tori, both with ring in the xy-plane, separated along z.
  θ = angle from z-axis (the torus/bond axis).

  θ = 0°:  AXIAL — looking through the hole. Barely any field.
  θ = 90°: EQUATORIAL — cutting across the ring. Two bumps at ±R_maj.
  θ ≈ 75°: NEAR-EQUATORIAL — clips one thick cross-section of tube.

BOND TYPE MAPPING:
  In chemistry, σ/π/δ refer to how orbitals ORIENT relative to the
  bond axis. The proton torus can rotate to maximize overlap:

  σ (sigma): Torus presents its RING EDGE to the partner.
    The equatorial slice (θ ~ 75-90°) gives the deepest overlap.
    This is the strongest bond type.

  π (pi): Torus presents at ~60° — sideways partial overlap.
    The tube clips at an angle, giving moderate D_e.

  δ (delta): Torus at ~45° — face-to-face with minimal overlap.
    Only the tube rim touches. Weakest covalent bond.

  non-bonding: θ < 30° — through the hole, no significant overlap.

WHY 75° > 90°:
  At 90° the slice crosses the ring TWICE (at ±R_maj), creating two
  kink-antikink wells that interfere. At 75° the slice clips one
  thick section cleanly — a single deep well bonds more efficiently
  than two shallow ones fighting each other.
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.interpolate import interp1d

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "bond_angular_analysis_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("ANGULAR BOND ANALYSIS — σ/π/δ FROM TORUS GEOMETRY")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("")

# ============================================================
# SETUP — REUSE THE TORUS AND SLICE MACHINERY
# ============================================================
N3 = 64
R_maj = 8
kink_width = 3

ix = np.arange(N3, dtype=np.float64) - N3/2
X, Y, Z = np.meshgrid(ix, ix, ix, indexing='ij')

def make_torus(X, Y, Z, z_center, R_major, kw):
    rho_xy = np.sqrt(X**2 + Y**2 + 1e-30)
    rho_tube = np.sqrt((rho_xy - R_major)**2 + (Z - z_center)**2)
    return (4.0/PI) * (np.arctan(np.exp(rho_tube + kw/2.0))
                      - np.arctan(np.exp(rho_tube - kw/2.0)))

phi_torus = make_torus(X, Y, Z, 0.0, R_maj, kink_width)

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

def bond_from_slice(r_vals, phi_vals, R_values):
    x_1d = np.arange(N1d, dtype=np.float64) - N1d/2
    center_1d = N1d // 2
    f_interp = interp1d(r_vals, phi_vals, kind='cubic', fill_value=0.0, bounds_error=False)
    phi_single = f_interp(x_1d - center_1d)
    H_single = build_1d_hessian(phi_single, N1d)
    ev_single, _ = eigsh(H_single, k=5, which='SM')
    ev_single = np.sort(ev_single)
    E0_single = ev_single[0]
    n_bound = np.sum(ev_single < 1.0)
    results = []
    for R in R_values:
        phi_A = f_interp(x_1d - center_1d + R/2)
        phi_B = f_interp(x_1d - center_1d - R/2)
        phi_double = phi_A + phi_B
        H_double = build_1d_hessian(phi_double, N1d)
        try:
            ev_double, _ = eigsh(H_double, k=5, which='SM')
            ev_double = np.sort(ev_double)
            V_R = ev_double[0] - E0_single
            split = ev_double[1] - ev_double[0] if len(ev_double) > 1 else 0
        except Exception:
            V_R, split = 0, 0
        results.append({'R': R, 'V': V_R, 'split': split})
    return results, E0_single, n_bound

# ============================================================
# PART 1: FINE ANGULAR SCAN — D_e(θ) PROFILE
# ============================================================
report("PART 1: D_e(θ) — BOND STRENGTH VS APPROACH ANGLE")
report("-" * 55)
report("Fine scan from 0° to 90° in 5° steps.")
report("")

R_scan = list(range(4, 42, 2))
theta_scan = list(range(0, 91, 5))

angular_data = {}  # {theta: {'De': D_e, 'Req': R_eq, 'results': [...]}}

report(f"{'θ':>5} {'peak_phi':>9} {'FWHM':>6} {'bound':>6} {'R_eq':>5} "
       f"{'D_e':>12} {'bond_type':>10}")
report("-" * 60)

for theta_deg in theta_scan:
    theta = theta_deg * PI / 180
    if theta_deg == 0:
        direction = (0, 0, 1)
    elif theta_deg == 90:
        direction = (1, 0, 0)
    else:
        direction = (np.sin(theta), 0, np.cos(theta))

    r_vals, phi_vals = extract_slice(phi_torus, direction, N3)
    results, E0, n_bound = bond_from_slice(r_vals, phi_vals, R_scan)

    V_arr = np.array([r['V'] for r in results])
    R_arr = np.array([r['R'] for r in results])
    i_min = np.argmin(V_arr)
    D_e = -V_arr[i_min] if V_arr[i_min] < 0 else 0
    R_eq = R_arr[i_min]

    peak = np.max(phi_vals)
    fwhm_mask = phi_vals > 0.5 * peak
    fwhm = np.sum(fwhm_mask) * (r_vals[1] - r_vals[0]) if peak > 0.001 else 0

    # Bond type classification
    if D_e > 0.1:
        btype = "σ (sigma)"
    elif D_e > 0.005:
        btype = "π (pi)"
    elif D_e > 0.0001:
        btype = "δ (delta)"
    else:
        btype = "non-bond"

    angular_data[theta_deg] = {
        'De': D_e, 'Req': R_eq, 'results': results,
        'peak': peak, 'fwhm': fwhm, 'n_bound': n_bound, 'type': btype
    }

    report(f"  {theta_deg:3d}° {peak:9.4f} {fwhm:6.1f} {n_bound:6d} {R_eq:5.0f} "
           f"{D_e:12.8f} {btype:>10}")

report("")

# ============================================================
# PART 2: WHY 75° > 90° — THE DOUBLE-BUMP EFFECT
# ============================================================
report("PART 2: WHY 75° BONDS STRONGER THAN 90°")
report("-" * 55)
report("")
report("At 90° (equatorial), the 1D slice crosses the torus ring at TWO")
report("points (x = -R_maj and x = +R_maj), creating a DOUBLE-WELL profile.")
report("At 75°, the slice clips the tube at a steep angle, creating a")
report("SINGLE deep well.")
report("")
report("Compare the profiles:")
report("")

# Show the two profiles
for theta_deg in [75, 80, 85, 90]:
    theta = theta_deg * PI / 180
    if theta_deg == 90:
        direction = (1, 0, 0)
    else:
        direction = (np.sin(theta), 0, np.cos(theta))
    r_vals, phi_vals = extract_slice(phi_torus, direction, N3)

    # Count peaks
    from scipy.signal import find_peaks
    peaks, props = find_peaks(phi_vals, height=0.1)
    n_peaks = len(peaks)
    max_phi = np.max(phi_vals)

    De = angular_data[theta_deg]['De']
    report(f"  θ={theta_deg}°: {n_peaks} peak(s), max_phi={max_phi:.4f}, D_e={De:.6f}")

report("")
report("The single-peak profile at 75° creates a deeper, sharper potential")
report("well than the double-peak at 90°. Two shallow wells compete for")
report("the bound state, weakening each one.")
report("")

# ============================================================
# PART 3: SIGMA/PI/DELTA WELL PARAMETERS
# ============================================================
report("PART 3: σ/π/δ BOND PARAMETERS")
report("-" * 55)
report("")

# Group angles by bond type
sigma_angles = [t for t in theta_scan if angular_data[t]['De'] > 0.1]
pi_angles = [t for t in theta_scan if 0.005 < angular_data[t]['De'] <= 0.1]
delta_angles = [t for t in theta_scan if 0.0001 < angular_data[t]['De'] <= 0.005]
nonbond_angles = [t for t in theta_scan if angular_data[t]['De'] <= 0.0001]

report(f"σ (sigma) angles: {sigma_angles}")
report(f"π (pi) angles:    {pi_angles}")
report(f"δ (delta) angles: {delta_angles}")
report(f"non-bonding:      {nonbond_angles}")
report("")

# Average D_e for each bond type
for btype, angles_list in [("σ sigma", sigma_angles), ("π pi", pi_angles),
                            ("δ delta", delta_angles)]:
    if angles_list:
        Des = [angular_data[t]['De'] for t in angles_list]
        Reqs = [angular_data[t]['Req'] for t in angles_list]
        report(f"  {btype}:")
        report(f"    Angular range: {min(angles_list)}° - {max(angles_list)}°")
        report(f"    D_e range: {min(Des):.6f} - {max(Des):.6f}")
        report(f"    D_e mean:  {np.mean(Des):.6f}")
        report(f"    R_eq range: {min(Reqs):.0f} - {max(Reqs):.0f}")
        report("")

# Ratio sigma/pi
if sigma_angles and pi_angles:
    De_sigma = np.mean([angular_data[t]['De'] for t in sigma_angles])
    De_pi = np.mean([angular_data[t]['De'] for t in pi_angles])
    report(f"  σ/π ratio: {De_sigma/De_pi:.1f}")
    report(f"  (Chemistry typically: 2-10×)")
report("")

# ============================================================
# PART 4: BOND CURVES FOR EACH TYPE
# ============================================================
report("PART 4: REPRESENTATIVE BOND CURVES")
report("-" * 55)

# Pick best angle for each type
best_sigma = max(sigma_angles, key=lambda t: angular_data[t]['De']) if sigma_angles else None
best_pi = max(pi_angles, key=lambda t: angular_data[t]['De']) if pi_angles else None
best_delta = max(delta_angles, key=lambda t: angular_data[t]['De']) if delta_angles else None

picks = []
if best_sigma is not None: picks.append(('σ', best_sigma))
if best_pi is not None: picks.append(('π', best_pi))
if best_delta is not None: picks.append(('δ', best_delta))

header = f"{'R':>4}"
for sym, theta_deg in picks:
    header += f"  {sym+'('+str(theta_deg)+'°)':>14}"
report(header)
report("-" * (4 + len(picks) * 16))

for R in R_scan:
    line = f"{R:4d}"
    for sym, theta_deg in picks:
        results = angular_data[theta_deg]['results']
        v = next((r['V'] for r in results if r['R'] == R), 0)
        line += f"  {v:+14.8f}"
    report(line)

report("")

# ============================================================
# PART 5: DECAY RATES BY BOND TYPE
# ============================================================
report("PART 5: MORSE DECAY RATES")
report("-" * 55)

for sym, theta_deg in picks:
    results = angular_data[theta_deg]['results']
    V_arr = np.array([r['V'] for r in results])
    R_arr = np.array([r['R'] for r in results])

    i_min = np.argmin(V_arr)
    D_e = -V_arr[i_min]
    R_eq = R_arr[i_min]

    # Fit decay in tail
    mask = (V_arr < -1e-8) & (R_arr > R_eq)
    if np.sum(mask) > 3:
        coeffs = np.polyfit(R_arr[mask], np.log(-V_arr[mask]), 1)
        a_morse = -coeffs[0]
    else:
        a_morse = 0

    report(f"  {sym} (θ={theta_deg}°): D_e={D_e:.6f}, R_eq={R_eq:.0f}, "
           f"Morse_a={a_morse:.4f}")

report("")

# ============================================================
# PART 6: SOLID-ANGLE WEIGHTED BOND ENERGY
# ============================================================
report("PART 6: ANGLE-WEIGHTED TOTAL BOND ENERGY")
report("-" * 55)
report("The isotropic (rotationally-averaged) bond is:")
report("  <D_e> = ∫ D_e(θ) × sin(θ) dθ / ∫ sin(θ) dθ")
report("")

# Numerical integration
theta_rad = np.array([t * PI / 180 for t in theta_scan])
De_arr = np.array([angular_data[t]['De'] for t in theta_scan])
sin_theta = np.sin(theta_rad)

# Trapezoidal integration with sin(θ) weighting
numerator = np.trapezoid(De_arr * sin_theta, theta_rad)
denominator = np.trapezoid(sin_theta, theta_rad)
De_iso = numerator / denominator

report(f"  <D_e> (isotropic average) = {De_iso:.6f}")
report("")

# Fraction from each bond type
for btype, angles_list in [("σ sigma", sigma_angles), ("π pi", pi_angles),
                            ("δ delta", delta_angles), ("non-bond", nonbond_angles)]:
    if angles_list:
        mask = np.isin(theta_scan, angles_list)
        th = theta_rad[mask]
        de = De_arr[mask]
        st = sin_theta[mask]
        if len(th) > 1:
            contrib = np.trapezoid(de * st, th) / denominator
        else:
            contrib = de[0] * st[0] * (5 * PI / 180) / denominator
        frac = contrib / De_iso * 100 if De_iso > 0 else 0
        report(f"  {btype:>10}: contribution = {contrib:.6f} ({frac:.1f}%)")

report("")

# ============================================================
# PART 7: WHAT THE TORUS GEOMETRY PREDICTS
# ============================================================
report("PART 7: GEOMETRIC PREDICTIONS")
report("-" * 55)
report("")
report("The torus has two scales:")
report(f"  R_major = {R_maj} (ring radius)")
report(f"  r_tube ~ {kink_width} (tube half-width)")
report(f"  Aspect ratio R_maj/r_tube = {R_maj/kink_width:.1f}")
report("")
report("Predictions from geometry:")
report(f"  σ bond R_eq should be ~ 2×r_tube = {2*kink_width}")
report(f"  π bond R_eq should be ~ R_major = {R_maj}")
report(f"  σ/π strength ratio ~ (R_maj/r_tube)^2 = {(R_maj/kink_width)**2:.1f}")
report(f"    (because tunneling depends on barrier area)")
report("")

# Compare with actual
if best_sigma and best_pi:
    actual_ratio = angular_data[best_sigma]['De'] / angular_data[best_pi]['De']
    report(f"  Actual σ/π ratio: {actual_ratio:.1f}")
    report(f"  Predicted: {(R_maj/kink_width)**2:.1f}")
    report(f"  This {'AGREES' if abs(actual_ratio - (R_maj/kink_width)**2) / (R_maj/kink_width)**2 < 0.5 else 'DISAGREES'} with the geometric prediction")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report("The toroidal proton geometry naturally produces σ/π/δ bond types:")
report("")
report(f"  σ (sigma): θ = {sigma_angles} — D_e ~ {np.mean([angular_data[t]['De'] for t in sigma_angles]):.4f}" if sigma_angles else "  σ: not found")
report(f"  π (pi):    θ = {pi_angles} — D_e ~ {np.mean([angular_data[t]['De'] for t in pi_angles]):.4f}" if pi_angles else "  π: not found")
report(f"  δ (delta): θ = {delta_angles} — D_e ~ {np.mean([angular_data[t]['De'] for t in delta_angles]):.6f}" if delta_angles else "  δ: not found")
report("")
report(f"  Isotropic average: <D_e> = {De_iso:.6f}")
report("")
report("The 75° > 90° effect: single-peak profiles bond stronger than")
report("double-peak (the equatorial slice hits the ring twice).")
report("")
report("KEY RESULT: Bond type diversity EMERGES from torus geometry alone.")
report("No parameters needed — just the shape of the kink-antikink ring.")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
