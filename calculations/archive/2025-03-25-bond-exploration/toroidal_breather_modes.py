"""
Toroidal Breather Modes — 8 Vertex Modes on the Torus
======================================================
The cube has 8 vertices, and each vertex maps to a breather orientation.
On the torus, these correspond to specific (m, n) Fourier modes.

The 24 breathers decompose under Oh:
  8 vertices × 3 rotations = 24 (vertex class)
  12 edges × 2 rotations = 24 (edge class)
  6 faces × 4 rotations = 24 (face class)

On the torus cross-section, the 8 vertex directions project as:
  Vertices of cube at (±1, ±1, ±1):
    - Toroidal projection (xy-plane): ±1, ±1 → m = ±1 (4 combinations)
    - Poloidal projection (z): ±1 → n = ±1
    - Combined: (m, n) = (±1, ±1) — these are TWIST modes!

So the 8 vertex breathers should show up as:
  4 modes at (m=+1, n=+1), (m=+1, n=-1), (m=-1, n=+1), (m=-1, n=-1)
  × 2 (time-reversal pairs) = 8

This predicts that the 8 vertex modes are ALL twist modes.
The 12 edge modes should be m=±1,n=0 or m=0,n=±1 (toroidal/poloidal).
The 6 face modes should be m=0,n=0 with different radial quantum numbers.

Let's test this by examining the single-torus eigenspectrum in detail.

Also: repeat all analysis with the 8 lowest-lying BREATHER profiles
explicitly placed on the torus, to see if they produce the predicted
shift pattern.
"""
import sys, io, os, time
import numpy as np
import cupy as cp
import cupyx.scipy.sparse as csp
import cupyx.scipy.sparse.linalg as csla
from scipy import sparse

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "toroidal_breather_modes_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("TOROIDAL BREATHER MODES — 8 VERTEX MODES")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("")

# ============================================================
# SETUP
# ============================================================
N = 64
Ntot = N**3
R_maj = 8
kink_width = 3

ix_1d = np.arange(N, dtype=np.float64) - N/2
X, Y, Z = np.meshgrid(ix_1d, ix_1d, ix_1d, indexing='ij')

def make_torus_proton(X, Y, Z, z_center, R_major, kw):
    rho_xy = np.sqrt(X**2 + Y**2 + 1e-30)
    rho_tube = np.sqrt((rho_xy - R_major)**2 + (Z - z_center)**2)
    phi = (4.0/PI) * (np.arctan(np.exp(rho_tube + kw/2.0))
                     - np.arctan(np.exp(rho_tube - kw/2.0)))
    return phi

def torus_angles(X, Y, Z, z_center, R_major):
    rho_xy = np.sqrt(X**2 + Y**2 + 1e-30)
    theta = np.arctan2(Y, X)
    phi_pol = np.arctan2(Z - z_center, rho_xy - R_major)
    rho_tube = np.sqrt((rho_xy - R_major)**2 + (Z - z_center)**2)
    return theta, phi_pol, rho_tube

def build_3d_hessian_gpu(phi_3d_flat):
    diag = 2.0 * d + np.cos(PI * phi_3d_flat)
    idx = np.arange(Ntot, dtype=np.int32)
    ix = idx % N
    iy = (idx // N) % N
    iz = idx // (N * N)
    neighbors = []
    for dx, dy, dz in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
        jx = (ix + dx) % N
        jy = (iy + dy) % N
        jz = (iz + dz) % N
        j = (jz * N * N + jy * N + jx).astype(np.int32)
        neighbors.append(j)
    all_rows = np.concatenate([idx] * 7)
    all_cols = np.concatenate([idx] + neighbors)
    all_vals = np.concatenate([diag] + [-np.ones(Ntot, dtype=np.float64)] * 6)
    H_cpu = sparse.coo_matrix((all_vals, (all_rows, all_cols)),
                               shape=(Ntot, Ntot)).tocsr()
    return csp.csr_matrix(H_cpu)

# ============================================================
# PART 1: FULL EIGENSPECTRUM — MAP TO Oh DECOMPOSITION
# ============================================================
report("PART 1: FULL EIGENSPECTRUM — 24 BREATHER MAPPING")
report("-" * 55)
report("")
report("Oh group has 24 elements = 8 vertices × 3 + 12 edges × 2 + 6 faces × 4")
report("On the torus:")
report("  VERTEX modes (8): (m, n) with |m|>0 AND |n|>0 → TWIST")
report("  EDGE modes (12): (m, n) with |m|>0,n=0 or m=0,|n|>0 → TOROIDAL or POLOIDAL")
report("  FACE modes (6): (m=0, n=0) with radial quantum number → RADIAL")
report("")

z_center = 0.0
phi_single = make_torus_proton(X, Y, Z, z_center, R_maj, kink_width)
theta_s, phi_pol_s, rho_tube_s = torus_angles(X, Y, Z, z_center, R_maj)

rho_flat = rho_tube_s.ravel()
theta_flat = theta_s.ravel()
phipol_flat = phi_pol_s.ravel()
torus_mask = rho_flat < 8.0

# Get many eigenvalues to capture all bound states
H_gpu = build_3d_hessian_gpu(phi_single.ravel())
n_eig = 30  # get more than 16 to be thorough
evals_gpu, evecs_gpu = csla.eigsh(H_gpu, k=n_eig, which='SA')
cp.cuda.Stream.null.synchronize()
evals = cp.asnumpy(evals_gpu)
evecs = cp.asnumpy(evecs_gpu)
idx_sort = np.argsort(evals)
evals = evals[idx_sort]
evecs = evecs[:, idx_sort]

E0_ref = evals[0]

# Count bound states
n_bound = np.sum(evals < 1.0)
report(f"Bound states found: {n_bound} (below mass gap ω²=1)")
report("")

# Full Fourier decomposition with high resolution
def full_decompose(evec, m_max=8, n_max=6):
    mask = torus_mask
    psi2 = evec[mask]**2
    th = theta_flat[mask]
    ph = phipol_flat[mask]

    power = {}
    for m in range(-m_max, m_max + 1):
        for n in range(-n_max, n_max + 1):
            c = np.sum(psi2 * np.exp(-1j * m * th) * np.exp(-1j * n * ph))
            power[(m, n)] = np.abs(c)**2
    total = sum(power.values()) + 1e-30

    # Classify by Oh mapping
    radial = power.get((0,0), 0) / total  # FACE type
    toroidal = sum(p for (m,n),p in power.items() if abs(m)>0 and n==0) / total  # EDGE (toroidal)
    poloidal = sum(p for (m,n),p in power.items() if m==0 and abs(n)>0) / total  # EDGE (poloidal)
    twist = sum(p for (m,n),p in power.items() if abs(m)>0 and abs(n)>0) / total  # VERTEX

    # Get dominant (m,n) pair
    sorted_p = sorted(power.items(), key=lambda x: -x[1])
    dom = sorted_p[0][0]
    sub = sorted_p[1][0] if len(sorted_p) > 1 else (0,0)

    return radial, toroidal, poloidal, twist, dom, sub, power, total

report(f"{'i':>3} {'ω²':>10} {'ω':>8} {'Oh_type':>10} {'rad':>6} {'tor':>6} "
       f"{'pol':>6} {'twi':>6} {'dom(m,n)':>10} {'sub(m,n)':>10}")
report("-" * 88)

# Count by type
n_radial = 0
n_toroidal = 0
n_poloidal = 0
n_twist = 0

eigenmode_data = []

for i in range(n_bound):
    ev = evals[i]
    rad, tor, pol, twi, dom, sub, power, total = full_decompose(evecs[:, i])

    if rad > max(tor, pol, twi):
        oh_type = "FACE"
        n_radial += 1
    elif tor > max(rad, pol, twi):
        oh_type = "EDGE-tor"
        n_toroidal += 1
    elif pol > max(rad, tor, twi):
        oh_type = "EDGE-pol"
        n_poloidal += 1
    else:
        oh_type = "VERTEX"
        n_twist += 1

    omega = np.sqrt(abs(ev))
    dom_str = f"({dom[0]:+d},{dom[1]:+d})"
    sub_str = f"({sub[0]:+d},{sub[1]:+d})"

    eigenmode_data.append({
        'eval': ev, 'omega': omega, 'type': oh_type,
        'rad': rad, 'tor': tor, 'pol': pol, 'twi': twi,
        'dom': dom, 'sub': sub
    })

    report(f"  {i:3d} {ev:10.6f} {omega:8.4f} {oh_type:>10} {rad:6.3f} {tor:6.3f} "
           f"{pol:6.3f} {twi:6.3f} {dom_str:>10} {sub_str:>10}")

report("")
report(f"Oh decomposition of {n_bound} bound states:")
report(f"  FACE (radial, m=0 n=0):    {n_radial}")
report(f"  EDGE-toroidal (|m|>0 n=0): {n_toroidal}")
report(f"  EDGE-poloidal (m=0 |n|>0): {n_poloidal}")
report(f"  VERTEX (twist, |m|>0 |n|>0): {n_twist}")
report(f"  Total: {n_radial + n_toroidal + n_poloidal + n_twist}")
report("")
report(f"Prediction: 6 face + 12 edge + 8 vertex = 26")
report(f"Or: modes count should reflect Oh irrep dimensions")
report("")

# ============================================================
# PART 2: DEGENERACY STRUCTURE
# ============================================================
report("PART 2: DEGENERACY STRUCTURE")
report("-" * 55)
report("Group modes by near-degenerate eigenvalues (tolerance 0.001).")
report("")

groups = []
current_group = [0]
for i in range(1, n_bound):
    if abs(evals[i] - evals[current_group[0]]) < 0.003:
        current_group.append(i)
    else:
        groups.append(current_group)
        current_group = [i]
groups.append(current_group)

report(f"{'group':>5} {'ω²':>10} {'deg':>4} {'modes':>20} {'types':>25}")
report("-" * 70)

for gi, group in enumerate(groups):
    ev = evals[group[0]]
    deg = len(group)
    mode_str = ",".join(str(g) for g in group)
    types = [eigenmode_data[g]['type'] for g in group]
    type_str = ",".join(types)

    # What toroidal harmonic m?
    dom_ms = [eigenmode_data[g]['dom'][0] for g in group]
    dom_ns = [eigenmode_data[g]['dom'][1] for g in group]

    report(f"  {gi:5d} {ev:10.6f} {deg:4d} [{mode_str:>18}] {type_str:>25}")

report("")

# ============================================================
# PART 3: BREATHER MODES ON TORUS — GWT FREQUENCY COMPARISON
# ============================================================
report("PART 3: TORUS EIGENVALUES VS GWT BREATHER FREQUENCIES")
report("-" * 55)
report("GWT predicts: ω_n = cos(n·γ), n = 1..24")
report("The torus eigenvalues should correspond to specific breather modes.")
report("")

report(f"{'torus_i':>7} {'ω²_torus':>10} {'ω_torus':>10} "
       f"{'best_n':>6} {'ω_gwt':>10} {'error%':>8}")
report("-" * 55)

gwt_omegas = [np.cos(n * gamma) for n in range(1, 25)]

for i in range(n_bound):
    ev = evals[i]
    omega = np.sqrt(abs(ev)) if ev > 0 else np.sqrt(abs(ev))

    # Find closest GWT mode
    errors = [abs(omega - gw) / gw * 100 for gw in gwt_omegas]
    best_n = np.argmin(errors) + 1
    best_err = errors[best_n - 1]

    report(f"  {i:7d} {ev:10.6f} {omega:10.6f} "
           f"{best_n:6d} {gwt_omegas[best_n-1]:10.6f} {best_err:8.2f}%")

report("")

# ============================================================
# PART 4: EXPLICIT 8-MODE ANALYSIS
# ============================================================
report("PART 4: THE 8 BREATHERS (first 8 bound modes)")
report("-" * 55)
report("For each of the first 8 modes, compute:")
report("  - Fourier structure (m, n)")
report("  - Spatial extent (how far does it reach?)")
report("  - Shift when neighbor torus present (R-dependence)")
report("")

# Get the 8 breather modes
first_8 = []
for i in range(min(8, n_bound)):
    ev = evals[i]
    evec = evecs[:, i]

    # Spatial extent: 90% containment radius
    psi2 = evec**2
    rho_bins = np.linspace(0, 12, 100)
    cum = 0
    total = np.sum(psi2[torus_mask])
    r90 = 12.0
    for j in range(len(rho_bins) - 1):
        mask_r = (rho_flat >= rho_bins[j]) & (rho_flat < rho_bins[j+1]) & torus_mask
        cum += np.sum(psi2[mask_r])
        if cum > 0.90 * total and r90 == 12.0:
            r90 = rho_bins[j+1]
            break

    data = eigenmode_data[i]
    first_8.append({
        'idx': i, 'eval': ev, 'omega': np.sqrt(abs(ev)),
        'r90': r90, **data
    })

    report(f"  Mode {i}: ω²={ev:+.6f}, ω={np.sqrt(abs(ev)):.4f}, "
           f"type={data['type']}, r90={r90:.1f} sites")

report("")

# Now scan R for each of the first 8 modes — track how each shifts
report("SHIFT OF FIRST 8 MODES VS R (two-torus system):")
report("")

R_values = [8, 10, 12, 14, 16, 20]

header = f"{'R':>4} "
for i in range(8):
    header += f"  {'dw2_'+str(i):>10}"
report(header)
report("-" * (4 + 8 * 12))

for R in R_values:
    z_A = -R / 2.0
    z_B = +R / 2.0
    phi_double = (make_torus_proton(X, Y, Z, z_A, R_maj, kink_width)
                + make_torus_proton(X, Y, Z, z_B, R_maj, kink_width))

    H_gpu = build_3d_hessian_gpu(phi_double.ravel())
    ev_d, _ = csla.eigsh(H_gpu, k=20, which='SA')
    cp.cuda.Stream.null.synchronize()
    ev_d = np.sort(cp.asnumpy(ev_d))

    # The two-torus system doubles each mode (bonding + antibonding)
    # So mode i maps to modes 2i (bonding) and 2i+1 (antibonding)
    line = f"{R:4d} "
    for i in range(8):
        if 2*i < len(ev_d):
            # Shift of bonding mode relative to single-torus
            shift = ev_d[2*i] - evals[i]
            line += f"  {shift:+10.6f}"
        else:
            line += f"  {'---':>10}"
    report(line)

report("")

# ============================================================
# PART 5: WHICH MODES COUPLE MOST TO NEIGHBOR?
# ============================================================
report("PART 5: RELATIVE COUPLING STRENGTH")
report("-" * 55)
report("Which of the 8 modes is most affected by a neighbor torus?")
report("Stronger shift = more exposed to neighbor = stronger bonding channel.")
report("")

report("At R=10 (moderate separation):")
z_A = -5.0
z_B = +5.0
phi_double = (make_torus_proton(X, Y, Z, z_A, R_maj, kink_width)
            + make_torus_proton(X, Y, Z, z_B, R_maj, kink_width))

H_gpu = build_3d_hessian_gpu(phi_double.ravel())
ev_d, _ = csla.eigsh(H_gpu, k=20, which='SA')
cp.cuda.Stream.null.synchronize()
ev_d = np.sort(cp.asnumpy(ev_d))

report(f"{'mode':>4} {'ω²_single':>12} {'ω²_bond':>12} {'ω²_anti':>12} "
       f"{'split':>12} {'shift':>12} {'type':>10}")
report("-" * 80)

for i in range(min(8, n_bound)):
    ev_single = evals[i]
    ev_bond = ev_d[2*i] if 2*i < len(ev_d) else ev_single
    ev_anti = ev_d[2*i+1] if 2*i+1 < len(ev_d) else ev_single
    split = ev_anti - ev_bond
    shift = ev_bond - ev_single
    oh_type = eigenmode_data[i]['type']

    report(f"  {i:4d} {ev_single:12.6f} {ev_bond:12.6f} {ev_anti:12.6f} "
           f"{split:12.8f} {shift:+12.8f} {oh_type:>10}")

report("")

# ============================================================
# PART 6: MAPPING TO CUBE VERTICES/EDGES/FACES
# ============================================================
report("PART 6: CUBE ELEMENT MAPPING")
report("-" * 55)
report("")
report("The d=3 cube has:")
report(f"  2^d = {2**d} vertices at (±1,±1,±1)")
report(f"  2d(d-1) = {2*d*(d-1)} edges")
report(f"  2d = {2*d} faces")
report(f"  Total orientation elements: |O| = 24")
report("")
report("Torus Fourier modes → cube elements:")
report("  (m=0, n=0): radial breathing → FACE normal direction (6)")
report("  (|m|>0, n=0): toroidal → EDGE along ring (toroidal edges)")
report("  (m=0, |n|>0): poloidal → EDGE through hole (poloidal edges)")
report("  (|m|>0, |n|>0): twist → VERTEX diagonal direction (8)")
report("")
report("Expected eigenvalue hierarchy (from torus geometry):")
report("  FACE modes: lowest (pure radial, most bound)")
report("  TOROIDAL edges: next (m²/R² excitation cost)")
report("  POLOIDAL edges: higher (n²/r² excitation cost, r < R)")
report("  VERTEX modes: highest (both m² and n² costs)")
report("")

# Check if the observed hierarchy matches
report("Observed hierarchy:")
for data in eigenmode_data[:min(n_bound, 20)]:
    report(f"  ω²={data['eval']:+.6f}, type={data['type']:>10}, "
           f"dom=({data['dom'][0]:+d},{data['dom'][1]:+d})")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report(f"Torus eigenspectrum: {n_bound} bound states")
report(f"  FACE (radial):      {n_radial:2d}   (expect: related to 2d=6)")
report(f"  EDGE (toroidal):    {n_toroidal:2d}   (expect: part of 2d(d-1)=12)")
report(f"  EDGE (poloidal):    {n_poloidal:2d}   (expect: part of 2d(d-1)=12)")
report(f"  VERTEX (twist):     {n_twist:2d}   (expect: related to 2^d=8)")
report("")
report("Degeneracy groups:")
for gi, group in enumerate(groups):
    deg = len(group)
    types = set(eigenmode_data[g]['type'] for g in group)
    report(f"  Group {gi}: deg={deg}, types={types}")
report("")
report("Key finding: the first 8 bound states decompose as:")
for i in range(min(8, n_bound)):
    d = eigenmode_data[i]
    report(f"  [{i}] {d['type']:>10}: ω²={d['eval']:+.6f}, "
           f"tor={d['tor']:.3f}, pol={d['pol']:.3f}, twi={d['twi']:.3f}")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
