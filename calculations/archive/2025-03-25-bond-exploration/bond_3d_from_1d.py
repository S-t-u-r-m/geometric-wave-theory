"""
3D Bond Surface from 1D Hessian Slices
========================================
Build the 3D bond potential cheaply by running the proven Hessian
eigenvalue method along multiple 1D directions through 3D space.

A toroidal proton has different interaction profiles depending on
the approach angle. The Oh group gives 3 unique direction classes:

  FACE (6 directions):   along a cube axis (e.g., z-axis)
    → approaches through the torus hole or along the ring
  EDGE (12 directions):  along a face diagonal (e.g., x+y)
    → approaches at 45° to the torus plane
  VERTEX (8 directions): along the body diagonal (x+y+z)
    → approaches at 54.7° to the torus plane

For a spherical proton, all 3 would give the same V(R).
For a toroidal proton, they give different V(R) — the ANISOTROPY
is the 3D information we're after.

The 3D bond surface is then:
  V(R, θ, φ) = V_face(R) × f_face(θ,φ) + V_edge(R) × f_edge(θ,φ) + V_vertex(R) × f_vertex(θ,φ)

where f are angular weight functions from Oh symmetry (known, not fit).

Method: Hessian eigenvalue splitting (proven in bond_3d_emerge.py).
Each 1D slice takes ~1 second. Total: ~3 seconds for full 3D picture.
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "bond_3d_from_1d_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("3D BOND SURFACE FROM 1D HESSIAN SLICES")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"d = {d}, gamma = {gamma:.10f}")
report("")

# ============================================================
# 1D HESSIAN BOND CALCULATION (proven method)
# ============================================================
N = 512
x = np.arange(N, dtype=np.float64)
center = N // 2

def kink_profile(x, x_center):
    return (4.0/PI) * np.arctan(np.exp(x - x_center))

def antikink_profile(x, x_center):
    return 2.0 - (4.0/PI) * np.arctan(np.exp(x - x_center))

def make_proton(x, center, width):
    return kink_profile(x, center - width/2) + antikink_profile(x, center + width/2) - 2.0

def build_hessian(phi, N):
    diag = 2.0 + np.cos(PI * phi)
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()

def bond_curve_1d(kink_width, R_values, label=""):
    """Compute V(R) using Hessian eigenvalue splitting.

    Two kink-antikink pairs at separation R.
    V(R) = lowest eigenvalue of double well - lowest eigenvalue of single well.
    """
    # Single well reference
    phi_single = make_proton(x, center, kink_width)
    H_single = build_hessian(phi_single, N)
    ev_single, _ = eigsh(H_single, k=3, which='SM')
    E0_single = np.sort(ev_single)[0]

    results = []
    for R in R_values:
        pos_A = center - R // 2
        pos_B = center + R // 2
        phi_double = make_proton(x, pos_A, kink_width) + make_proton(x, pos_B, kink_width)
        H_double = build_hessian(phi_double, N)
        ev_double, _ = eigsh(H_double, k=3, which='SM')
        ev_double = np.sort(ev_double)

        V_R = ev_double[0] - E0_single
        split = ev_double[1] - ev_double[0] if len(ev_double) > 1 else 0
        results.append({'R': R, 'V': V_R, 'split': split, 'E0': ev_double[0]})

    return results, E0_single

# ============================================================
# THE THREE DIRECTION CLASSES
# ============================================================
# On the d=3 cubic lattice with a toroidal proton:
#
# The kink-antikink pair has an intrinsic WIDTH that depends on
# the approach direction relative to the torus orientation.
#
# For a torus with major radius R_maj and tube radius r_tube:
#   FACE approach (along torus axis): sees the tube cross-section
#     → effective kink width ~ 2 × r_tube
#   EDGE approach (45° to axis): sees an elliptical cross-section
#     → effective kink width ~ sqrt(2) × r_tube (geometric mean)
#   VERTEX approach (54.7° to axis): sees the most oblique section
#     → effective kink width ~ sqrt(3) × r_tube
#
# But for the 1D effective model, the key parameter is the
# EFFECTIVE KINK WIDTH along each direction.
#
# The torus has two characteristic lengths:
#   R_maj = 8 (major radius, from the 3D sim)
#   r_tube ~ kink_width = 3 (tube radius)
#
# Along the torus axis (FACE): two protons approach through the holes.
#   The kink profiles overlap along the axis. Effective width = kink_width.
#   This is the standard 1D calculation.
#
# Perpendicular to axis (FACE, in-plane): approach around the ring.
#   The overlap is through the ring circumference.
#   Effective width ~ kink_width (same radial profile).
#   BUT the effective separation is the chord distance, not R.
#
# For the 1D model, the DIRECTION enters through:
#   1. The effective kink width (how broad the proton looks from that angle)
#   2. The effective potential depth (how strong the interaction along that line)
#
# On the cubic lattice, these are determined by the number of nearest
# neighbors along each direction:
#   FACE: 1 NN at distance 1
#   EDGE: √2 NN at distance √2
#   VERTEX: √3 NN at distance √3

report("THREE DIRECTION CLASSES (Oh symmetry):")
report("-" * 55)
report("")
report("  FACE   (6 dirs): along cube axis")
report("    Approach through torus hole or along ring")
report("    Lattice spacing: a = 1")
report("    Multiplicity: 6")
report("")
report("  EDGE   (12 dirs): along face diagonal")
report("    Approach at 45° to torus plane")
report("    Lattice spacing: a√2 = 1.414")
report("    Multiplicity: 12")
report("")
report("  VERTEX (8 dirs): along body diagonal")
report("    Approach at 54.7° to torus plane")
report("    Lattice spacing: a√3 = 1.732")
report("    Multiplicity: 8")
report("")

# For the effective 1D model along each direction:
# The SEPARATION R in 3D maps to an effective 1D coordinate r:
#   FACE:   r = R (direct mapping)
#   EDGE:   r = R/√2 (each lattice step covers √2 in 3D)
#   VERTEX: r = R/√3 (each lattice step covers √3 in 3D)
#
# Also, the coupling strength scales with the number of NNs:
# The discrete Laplacian has coupling 1 per NN.
# Along FACE: 1 neighbor per step
# Along EDGE: 2 NNs contribute (from the 2 axes), coupling per step is 2
# Along VERTEX: 3 NNs contribute, coupling per step is 3
#
# BUT in the 1D Hessian the coupling is already 1 per bond.
# The direction changes the EFFECTIVE KINK SHAPE (how the 3D kink
# projects onto the 1D line) and the EFFECTIVE SPACING.

# ============================================================
# PART 1: STANDARD FACE DIRECTION (baseline)
# ============================================================
report("PART 1: FACE DIRECTION (standard, along cube axis)")
report("-" * 55)

R_scan = list(range(4, 52, 2))
kw_face = 3  # standard kink width

face_results, E0_face = bond_curve_1d(kw_face, R_scan, "face")

report(f"Single well E0 = {E0_face:.8f}")
report(f"{'R':>4} {'V(R)':>12} {'split':>12}")
report("-" * 30)
for r in face_results:
    marker = ""
    if r['V'] == min(rr['V'] for rr in face_results):
        marker = " <-- min"
    report(f"{r['R']:4d} {r['V']:+12.8f} {r['split']:12.8f}{marker}")

# Find minimum
V_arr = np.array([r['V'] for r in face_results])
R_arr = np.array([r['R'] for r in face_results])
i_min = np.argmin(V_arr)
report(f"\nFACE: R_eq = {R_arr[i_min]}, D_e = {-V_arr[i_min]:.8f}")
report("")

# ============================================================
# PART 2: EDGE DIRECTION (face diagonal, effective spacing √2)
# ============================================================
report("PART 2: EDGE DIRECTION (face diagonal)")
report("-" * 55)
report("Along the face diagonal, each 1D step = √2 in 3D.")
report("The kink appears WIDER (projected onto the longer axis).")
report("")

# The kink width along the diagonal is wider by √2
# because the profile is spread over more lattice sites when
# projected onto the face diagonal.
kw_edge = kw_face * np.sqrt(2)

# The R values in 3D map to r = R/√2 in the 1D effective model
# But we keep R as the 3D separation for comparison
# So we use R_1d = R / sqrt(2) as the 1D separation
R_3d_scan = list(range(4, 52, 2))

edge_results, E0_edge = bond_curve_1d(kw_edge, R_3d_scan, "edge")

report(f"Effective kink width (edge): {kw_edge:.3f}")
report(f"{'R_3d':>5} {'R_1d':>6} {'V(R)':>12} {'split':>12}")
report("-" * 40)
for r in edge_results:
    R_3d = r['R'] * np.sqrt(2)  # convert back to 3D
    report(f"{R_3d:5.1f} {r['R']:6d} {r['V']:+12.8f} {r['split']:12.8f}")

V_edge = np.array([r['V'] for r in edge_results])
R_edge = np.array([r['R'] for r in edge_results])
i_min_e = np.argmin(V_edge)
report(f"\nEDGE: R_eq(1d) = {R_edge[i_min_e]}, D_e = {-V_edge[i_min_e]:.8f}")
report(f"       R_eq(3d) = {R_edge[i_min_e] * np.sqrt(2):.1f}")
report("")

# ============================================================
# PART 3: VERTEX DIRECTION (body diagonal, effective spacing √3)
# ============================================================
report("PART 3: VERTEX DIRECTION (body diagonal)")
report("-" * 55)

kw_vertex = kw_face * np.sqrt(3)

vertex_results, E0_vertex = bond_curve_1d(kw_vertex, R_3d_scan, "vertex")

report(f"Effective kink width (vertex): {kw_vertex:.3f}")
report(f"{'R_3d':>5} {'R_1d':>6} {'V(R)':>12} {'split':>12}")
report("-" * 40)
for r in vertex_results:
    R_3d = r['R'] * np.sqrt(3)
    report(f"{R_3d:5.1f} {r['R']:6d} {r['V']:+12.8f} {r['split']:12.8f}")

V_vertex = np.array([r['V'] for r in vertex_results])
R_vertex = np.array([r['R'] for r in vertex_results])
i_min_v = np.argmin(V_vertex)
report(f"\nVERTEX: R_eq(1d) = {R_vertex[i_min_v]}, D_e = {-V_vertex[i_min_v]:.8f}")
report(f"         R_eq(3d) = {R_vertex[i_min_v] * np.sqrt(3):.1f}")
report("")

# ============================================================
# PART 4: COMPARISON — THE ANGULAR ANISOTROPY
# ============================================================
report("PART 4: ANGULAR ANISOTROPY")
report("-" * 55)
report("At the same 3D distance R, how does V depend on direction?")
report("")

# For each 3D distance, interpolate V from each direction
from scipy.interpolate import interp1d

# Build interpolators (1D R values → V)
R_face_arr = np.array([r['R'] for r in face_results], dtype=float)
V_face_arr = np.array([r['V'] for r in face_results])
f_face = interp1d(R_face_arr, V_face_arr, kind='cubic', fill_value=0, bounds_error=False)

R_edge_1d = np.array([r['R'] for r in edge_results], dtype=float)
V_edge_arr = np.array([r['V'] for r in edge_results])
# Edge: R_3d = R_1d * sqrt(2), so R_1d = R_3d / sqrt(2)
f_edge = interp1d(R_edge_1d * np.sqrt(2), V_edge_arr, kind='cubic', fill_value=0, bounds_error=False)

R_vert_1d = np.array([r['R'] for r in vertex_results], dtype=float)
V_vert_arr = np.array([r['V'] for r in vertex_results])
f_vertex = interp1d(R_vert_1d * np.sqrt(3), V_vert_arr, kind='cubic', fill_value=0, bounds_error=False)

report(f"{'R_3d':>5} {'V_face':>12} {'V_edge':>12} {'V_vertex':>12} "
       f"{'anisotropy':>12}")
report("-" * 58)

for R3d in [6, 8, 10, 12, 14, 16, 20, 24, 30, 40]:
    vf = float(f_face(R3d))
    ve = float(f_edge(R3d))
    vv = float(f_vertex(R3d))

    # Anisotropy: max/min ratio
    vals = [vf, ve, vv]
    neg_vals = [v for v in vals if v < 0]
    if len(neg_vals) >= 2:
        aniso = max(neg_vals) / min(neg_vals) if min(neg_vals) != 0 else 0
    else:
        aniso = 0

    report(f"{R3d:5d} {vf:+12.8f} {ve:+12.8f} {vv:+12.8f} {aniso:12.4f}")

report("")

# ============================================================
# PART 5: ISOTROPIC AVERAGE (Oh-weighted)
# ============================================================
report("PART 5: Oh-WEIGHTED ISOTROPIC AVERAGE")
report("-" * 55)
report("V_iso(R) = (6·V_face + 12·V_edge + 8·V_vertex) / 26")
report("This is the spherically-averaged potential (A1g component).")
report("")

report(f"{'R_3d':>5} {'V_face':>10} {'V_edge':>10} {'V_vertex':>10} "
       f"{'V_iso':>10} {'V_face/V_iso':>12}")
report("-" * 60)

for R3d in [6, 8, 10, 12, 14, 16, 20, 24, 30, 40]:
    vf = float(f_face(R3d))
    ve = float(f_edge(R3d))
    vv = float(f_vertex(R3d))

    V_iso = (6 * vf + 12 * ve + 8 * vv) / 26.0
    ratio = vf / V_iso if abs(V_iso) > 1e-12 else 0

    report(f"{R3d:5d} {vf:+10.6f} {ve:+10.6f} {vv:+10.6f} "
           f"{V_iso:+10.6f} {ratio:12.4f}")

report("")

# ============================================================
# PART 6: BOND WELL PARAMETERS FOR EACH DIRECTION
# ============================================================
report("PART 6: MORSE WELL PARAMETERS BY DIRECTION")
report("-" * 55)

for name, results, scale in [("FACE", face_results, 1.0),
                               ("EDGE", edge_results, np.sqrt(2)),
                               ("VERTEX", vertex_results, np.sqrt(3))]:
    V = np.array([r['V'] for r in results])
    R = np.array([r['R'] for r in results]) * scale  # convert to 3D

    i_min = np.argmin(V)
    D_e = -V[i_min]
    R_eq = R[i_min]

    # Fit decay rate from the tail
    mask = (V < -1e-8) & (R > R_eq)
    if np.sum(mask) > 3:
        coeffs = np.polyfit(R[mask], np.log(-V[mask]), 1)
        a_morse = -coeffs[0]
    else:
        a_morse = 0

    report(f"  {name:8s}: R_eq = {R_eq:5.1f}, D_e = {D_e:.8f}, decay = {a_morse:.4f}")

report("")

# ============================================================
# PART 7: THE 3D PICTURE — WHAT DOES THE SURFACE LOOK LIKE?
# ============================================================
report("PART 7: 3D BOND SURFACE CROSS-SECTIONS")
report("-" * 55)
report("V(R, θ) where θ = angle from torus axis:")
report("  θ = 0°:    FACE direction (through the hole)")
report("  θ = 45°:   EDGE direction (face diagonal)")
report("  θ = 54.7°: VERTEX direction (body diagonal)")
report("  θ = 90°:   FACE direction (in the torus plane)")
report("")
report("(θ = 0° and θ = 90° both map to FACE — torus has axial symmetry)")
report("")

# Angular interpolation between directions
# θ = 0 → face, θ = π/4 → edge, θ = arctan(√2) ≈ 54.7° → vertex
theta_face = 0.0
theta_edge = PI / 4.0
theta_vertex = np.arctan(np.sqrt(2))  # 54.7°

report(f"{'R':>4} {'θ=0°':>10} {'θ=22.5°':>10} {'θ=45°':>10} {'θ=54.7°':>10}")
report("-" * 48)

for R3d in [6, 8, 10, 12, 14, 16, 20]:
    vf = float(f_face(R3d))
    ve = float(f_edge(R3d))
    vv = float(f_vertex(R3d))

    # Simple angular interpolation: V(θ) between known directions
    # θ = 0 → face, θ = 45° → edge, θ = 54.7° → vertex
    # Linear interpolation for intermediate angles
    v_22 = (vf + ve) / 2  # halfway between face and edge

    report(f"{R3d:4d} {vf:+10.6f} {v_22:+10.6f} {ve:+10.6f} {vv:+10.6f}")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report("3D bond surface reconstructed from 3 unique 1D Hessian calculations.")
report("")

# Get the minima
V_face_min = min(r['V'] for r in face_results)
V_edge_min = min(r['V'] for r in edge_results)
V_vertex_min = min(r['V'] for r in vertex_results)

report(f"Well depths:")
report(f"  FACE   (6 dirs):  D_e = {-V_face_min:.8f}")
report(f"  EDGE   (12 dirs): D_e = {-V_edge_min:.8f}")
report(f"  VERTEX (8 dirs):  D_e = {-V_vertex_min:.8f}")
report("")

if V_face_min < 0 and V_edge_min < 0 and V_vertex_min < 0:
    # Anisotropy
    depths = [-V_face_min, -V_edge_min, -V_vertex_min]
    report(f"Anisotropy (max/min depth): {max(depths)/min(depths):.2f}")
    report(f"The interaction is {'ISOTROPIC' if max(depths)/min(depths) < 1.5 else 'ANISOTROPIC'}")
    report("")

    # Which direction bonds strongest?
    strongest = ["FACE", "EDGE", "VERTEX"][np.argmax(depths)]
    report(f"Strongest bonding direction: {strongest}")
    report(f"  FACE → through torus hole / along ring")
    report(f"  EDGE → at 45° to torus plane")
    report(f"  VERTEX → along body diagonal")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
