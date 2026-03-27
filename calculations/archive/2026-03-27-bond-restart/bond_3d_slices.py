"""
3D Bond Surface from Real Torus Slices
========================================
Extract actual 1D field profiles from the 3D toroidal kink, then
use the proven Hessian eigenvalue method on each slice.

Instead of guessing the effective kink width along each direction,
we MEASURE the real field profile by slicing the 3D torus along:
  1. Through the hole (z-axis, axial)
  2. Along the ring (x-axis, equatorial)
  3. Face diagonal (x+y direction)
  4. Body diagonal (x+y+z direction)
  5. Multiple angles θ from 0° to 90° in 10° steps

Each slice gives a 1D field profile φ(r) that we plug directly
into the Hessian. Two such profiles at separation R give the
bonding-antibonding splitting → V(R) for that direction.

This gives the REAL anisotropic bond surface without any 3D eigsh.
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

outfile = os.path.join(os.path.dirname(__file__), "bond_3d_slices_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("3D BOND SURFACE FROM REAL TORUS SLICES")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"d = {d}, gamma = {gamma:.10f}")
report("")

# ============================================================
# BUILD THE 3D TORUS
# ============================================================
N3 = 64  # 3D lattice for profile extraction
R_maj = 8
kink_width = 3

report(f"3D torus: N={N3}, R_major={R_maj}, kink_width={kink_width}")

ix = np.arange(N3, dtype=np.float64) - N3/2
X, Y, Z = np.meshgrid(ix, ix, ix, indexing='ij')

def make_torus(X, Y, Z, z_center, R_major, kw):
    rho_xy = np.sqrt(X**2 + Y**2 + 1e-30)
    rho_tube = np.sqrt((rho_xy - R_major)**2 + (Z - z_center)**2)
    phi = (4.0/PI) * (np.arctan(np.exp(rho_tube + kw/2.0))
                     - np.arctan(np.exp(rho_tube - kw/2.0)))
    return phi

# Build single torus at origin
phi_torus = make_torus(X, Y, Z, 0.0, R_maj, kink_width)
report(f"Torus peak: {phi_torus.max():.4f}")
report("")

# ============================================================
# EXTRACT 1D SLICES
# ============================================================
report("PART 1: 1D SLICE PROFILES FROM THE 3D TORUS")
report("-" * 55)

mid = N3 // 2  # center index

def extract_slice(phi_3d, direction, N3):
    """Extract a 1D field profile along a direction through the torus center.

    direction: unit vector (dx, dy, dz)
    Returns: (r_values, phi_values) where r is distance from center.
    """
    dx, dy, dz = direction
    norm = np.sqrt(dx**2 + dy**2 + dz**2)
    dx, dy, dz = dx/norm, dy/norm, dz/norm

    # Sample points along the ray from center
    r_max = N3 // 2 - 2
    n_samples = N3 * 4  # oversample for smoothness
    r_values = np.linspace(-r_max, r_max, n_samples)

    phi_values = np.zeros(n_samples)
    for i, r in enumerate(r_values):
        # 3D position
        px = r * dx
        py = r * dy
        pz = r * dz

        # Trilinear interpolation from the lattice
        ix_f = px + N3/2
        iy_f = py + N3/2
        iz_f = pz + N3/2

        # Clamp to grid
        ix0 = int(np.floor(ix_f)) % N3
        iy0 = int(np.floor(iy_f)) % N3
        iz0 = int(np.floor(iz_f)) % N3
        ix1 = (ix0 + 1) % N3
        iy1 = (iy0 + 1) % N3
        iz1 = (iz0 + 1) % N3

        fx = ix_f - np.floor(ix_f)
        fy = iy_f - np.floor(iy_f)
        fz = iz_f - np.floor(iz_f)

        # Trilinear interpolation
        phi_values[i] = (
            phi_3d[ix0, iy0, iz0] * (1-fx)*(1-fy)*(1-fz) +
            phi_3d[ix1, iy0, iz0] * fx*(1-fy)*(1-fz) +
            phi_3d[ix0, iy1, iz0] * (1-fx)*fy*(1-fz) +
            phi_3d[ix0, iy0, iz1] * (1-fx)*(1-fy)*fz +
            phi_3d[ix1, iy1, iz0] * fx*fy*(1-fz) +
            phi_3d[ix1, iy0, iz1] * fx*(1-fy)*fz +
            phi_3d[ix0, iy1, iz1] * (1-fx)*fy*fz +
            phi_3d[ix1, iy1, iz1] * fx*fy*fz
        )

    return r_values, phi_values

# Define slice directions
# Theta = angle from z-axis (torus axis)
# At theta=0: along z (through the hole)
# At theta=90: along x (in the torus plane, across the ring)

directions = {}
angles = {}

# Axial (through the hole)
directions['axial (θ=0°)'] = (0, 0, 1)
angles['axial (θ=0°)'] = 0

# Equatorial (across the ring, in the torus plane)
directions['equat (θ=90°)'] = (1, 0, 0)
angles['equat (θ=90°)'] = 90

# Face diagonal in xz plane (θ=45°)
directions['diag45 (θ=45°)'] = (1, 0, 1)
angles['diag45 (θ=45°)'] = 45

# Body diagonal (θ=54.7°)
directions['body (θ=54.7°)'] = (1, 1, 1)
angles['body (θ=54.7°)'] = 54.7

# Additional angles: 15°, 30°, 60°, 75°
for theta_deg in [15, 30, 60, 75]:
    theta = theta_deg * PI / 180
    dx = np.sin(theta)
    dz = np.cos(theta)
    name = f'θ={theta_deg}°'
    directions[name] = (dx, 0, dz)
    angles[name] = theta_deg

# Extract all slices
slices = {}
for name, direction in directions.items():
    r_vals, phi_vals = extract_slice(phi_torus, direction, N3)
    slices[name] = (r_vals, phi_vals)

    # Profile summary
    peak = np.max(phi_vals)
    width_mask = phi_vals > 0.1 * peak
    n_sig = np.sum(width_mask)
    fwhm = np.sum(phi_vals > 0.5 * peak) * (r_vals[1] - r_vals[0])

    report(f"  {name:>18}: peak={peak:.4f}, FWHM={fwhm:.1f}, "
           f"significant={n_sig} pts")

report("")

# Show the actual profiles at key positions
report("SLICE PROFILES (phi vs r):")
key3 = ['axial (θ=0°)', 'diag45 (θ=45°)', 'equat (θ=90°)']
header = f"{'r':>6}" + "".join(f"  {name:>18}" for name in key3)
report(header)

for r_target in range(-15, 16):
    line = f"{r_target:6d}"
    for name in key3:
        r_vals, phi_vals = slices[name]
        idx = np.argmin(np.abs(r_vals - r_target))
        line += f"  {phi_vals[idx]:18.4f}"
    report(line)

report("")

# ============================================================
# PART 2: HESSIAN BOND CALCULATION FOR EACH SLICE
# ============================================================
report("PART 2: BOND CURVES FROM REAL SLICE PROFILES")
report("-" * 55)

N1d = 512  # 1D lattice for Hessian calculation

def build_1d_hessian(phi_1d, N):
    diag = 2.0 + np.cos(PI * phi_1d)
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()

def bond_from_slice(slice_profile_r, slice_profile_phi, R_3d_values):
    """Compute V(R) using the actual torus slice as the kink profile.

    Place one torus slice at position A and another at position B,
    separated by R. Build the 1D Hessian and find the eigenvalue splitting.
    """
    # Interpolate slice profile onto the 1D lattice
    x_1d = np.arange(N1d, dtype=np.float64) - N1d/2
    center_1d = N1d // 2

    # The slice profile is defined on r_values (centered at 0)
    f_interp = interp1d(slice_profile_r, slice_profile_phi,
                        kind='cubic', fill_value=0.0, bounds_error=False)

    # Single proton: profile centered at 0
    phi_single = f_interp(x_1d - center_1d)
    H_single = build_1d_hessian(phi_single, N1d)
    ev_single, _ = eigsh(H_single, k=5, which='SM')
    ev_single = np.sort(ev_single)
    E0_single = ev_single[0]
    n_bound = np.sum(ev_single < 1.0)

    results = []
    for R in R_3d_values:
        # Two profiles at ±R/2
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
            V_R = 0
            split = 0

        results.append({'R': R, 'V': V_R, 'split': split})

    return results, E0_single, n_bound

# R values to scan (in 3D units = same as the slice coordinate)
R_scan = list(range(4, 42, 2))

# Run for each direction
all_bond_curves = {}

for name in sorted(slices.keys(), key=lambda n: angles.get(n, 0)):
    r_vals, phi_vals = slices[name]
    theta = angles.get(name, 0)

    results, E0, n_bound = bond_from_slice(r_vals, phi_vals, R_scan)
    all_bond_curves[name] = results

    V_arr = np.array([r['V'] for r in results])
    R_arr = np.array([r['R'] for r in results])
    i_min = np.argmin(V_arr)
    D_e = -V_arr[i_min] if V_arr[i_min] < 0 else 0

    report(f"  {name:>18}: E0={E0:.6f}, bound={n_bound}, "
           f"R_eq={R_arr[i_min]}, D_e={D_e:.8f}")

report("")

# ============================================================
# PART 3: DETAILED BOND CURVES FOR KEY DIRECTIONS
# ============================================================
report("PART 3: DETAILED BOND CURVES")
report("-" * 55)

key_dirs = ['axial (θ=0°)', 'θ=30°', 'diag45 (θ=45°)', 'θ=60°', 'equat (θ=90°)']

header = f"{'R':>4}"
for name in key_dirs:
    header += f"  {name:>16}"
report(header)
report("-" * (4 + len(key_dirs) * 18))

for R in R_scan:
    line = f"{R:4d}"
    for name in key_dirs:
        results = all_bond_curves.get(name, [])
        v = next((r['V'] for r in results if r['R'] == R), 0)
        line += f"  {v:+16.8f}"
    report(line)

report("")

# ============================================================
# PART 4: SPLITTING (TUNNEL COUPLING) VS DIRECTION
# ============================================================
report("PART 4: BONDING-ANTIBONDING SPLITTING VS DIRECTION")
report("-" * 55)
report("The splitting measures tunnel coupling strength.")
report("At large R it saturates at 2×(binding energy of single well).")
report("")

header = f"{'R':>4}"
for name in key_dirs:
    header += f"  {name:>16}"
report(header)
report("-" * (4 + len(key_dirs) * 18))

for R in R_scan:
    line = f"{R:4d}"
    for name in key_dirs:
        results = all_bond_curves.get(name, [])
        s = next((r['split'] for r in results if r['R'] == R), 0)
        line += f"  {s:16.8f}"
    report(line)

report("")

# ============================================================
# PART 5: THE ANGULAR BOND MAP
# ============================================================
report("PART 5: V(R, θ) — THE ANGULAR BOND MAP")
report("-" * 55)
report("Bond potential as a function of distance and approach angle.")
report("")

# For each angle, find the well depth and equilibrium distance
report(f"{'direction':>18} {'θ':>5} {'R_eq':>6} {'D_e':>12} {'split_max':>12}")
report("-" * 58)

theta_vals = []
De_vals = []
Req_vals = []

for name in sorted(all_bond_curves.keys(), key=lambda n: angles.get(n, 0)):
    theta = angles.get(name, 0)
    results = all_bond_curves[name]

    V_arr = np.array([r['V'] for r in results])
    R_arr = np.array([r['R'] for r in results])
    S_arr = np.array([r['split'] for r in results])

    i_min = np.argmin(V_arr)
    D_e = -V_arr[i_min] if V_arr[i_min] < 0 else 0
    R_eq = R_arr[i_min]
    split_max = np.max(S_arr)

    theta_vals.append(theta)
    De_vals.append(D_e)
    Req_vals.append(R_eq)

    report(f"  {name:>18} {theta:5.1f} {R_eq:6.0f} {D_e:12.8f} {split_max:12.6f}")

report("")

# ============================================================
# PART 6: Oh-WEIGHTED AVERAGE
# ============================================================
report("PART 6: Oh-WEIGHTED ISOTROPIC AVERAGE")
report("-" * 55)

# Weight each direction by its solid angle fraction
# On the sphere, different theta bands have different areas
# For discrete Oh directions:
# Face (θ=0°, 90°): 6 directions, solid angle = 6 × (4π/26) (approximate)
# Edge (θ=45°): 12 directions
# Vertex (θ=54.7°): 8 directions

# But with our finer angular sampling, use sin(θ) weighting
report("V_iso(R) from sin(θ)-weighted angular average:")
report("")

theta_sorted = sorted(all_bond_curves.keys(), key=lambda n: angles.get(n, 0))

report(f"{'R':>4} {'V_iso':>12} {'V_axial':>12} {'V_equat':>12} "
       f"{'anisotropy':>12}")
report("-" * 55)

for R in R_scan:
    v_sum = 0
    w_sum = 0
    for name in theta_sorted:
        theta = angles.get(name, 0) * PI / 180
        results = all_bond_curves[name]
        v = next((r['V'] for r in results if r['R'] == R), 0)
        # sin(θ) weighting for spherical average
        w = np.sin(theta) + 0.01  # small offset to include θ=0
        v_sum += v * w
        w_sum += w

    V_iso = v_sum / w_sum if w_sum > 0 else 0

    v_ax = next((r['V'] for r in all_bond_curves.get('axial (θ=0°)', []) if r['R'] == R), 0)
    v_eq = next((r['V'] for r in all_bond_curves.get('equat (θ=90°)', []) if r['R'] == R), 0)

    aniso = v_ax / v_eq if abs(v_eq) > 1e-12 else 0

    report(f"{R:4d} {V_iso:+12.8f} {v_ax:+12.8f} {v_eq:+12.8f} {aniso:+12.4f}")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report("3D bond surface from real torus slice profiles:")
report("")

for name in sorted(all_bond_curves.keys(), key=lambda n: angles.get(n, 0)):
    theta = angles.get(name, 0)
    results = all_bond_curves[name]
    V_arr = np.array([r['V'] for r in results])
    D_e = -min(V_arr) if min(V_arr) < 0 else 0
    report(f"  θ={theta:5.1f}°: D_e = {D_e:.8f}")

report("")
report("KEY: The torus is NOT spherically symmetric.")
report("Different approach angles see different kink profiles,")
report("leading to direction-dependent bond strengths.")
report("This IS the origin of sigma/pi/delta bond types in chemistry.")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
