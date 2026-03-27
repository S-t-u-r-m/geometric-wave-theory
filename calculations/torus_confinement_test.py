"""
Torus Confinement Test — Does the Tachyon Stabilize in 3D?
============================================================
In 1D, the kink-antikink pair has omega^2 = -0.372 (TACHYON).
This is the annihilation channel — the kink can slide into the antikink.

On the 3D torus, topology should PREVENT annihilation.
The tachyon mode should stabilize to omega^2 >= 0.

If it does: confinement emerges from topology.
If it doesn't: the torus model needs revision.

Method:
  1. Build a single toroidal kink on 64^3 lattice
  2. Compute ALL eigenvalues near and below the mass gap
  3. Check: are ANY eigenvalues negative?
  4. If all >= 0: tachyon is GONE (topology stabilizes)
  5. Compare the mode structure to the 1D case

GPU: RTX 4070 Ti, CuPy. N=64 eigsh takes ~3s.
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

outfile = os.path.join(os.path.dirname(__file__), "torus_confinement_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("TORUS CONFINEMENT TEST — DOES THE TACHYON STABILIZE?")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("")

# ============================================================
# 3D TORUS SETUP
# ============================================================
N = 64
Ntot = N**3
R_maj = 8
kink_width = 3

report(f"Lattice: {N}^3 = {Ntot:,} sites")
report(f"Torus: R_major = {R_maj}, kink_width = {kink_width}")
report("")

ix = np.arange(N, dtype=np.float64) - N/2
X, Y, Z = np.meshgrid(ix, ix, ix, indexing='ij')

def make_torus(X, Y, Z, z_center, R_major, kw):
    rho_xy = np.sqrt(X**2 + Y**2 + 1e-30)
    rho_tube = np.sqrt((rho_xy - R_major)**2 + (Z - z_center)**2)
    return (4.0/PI) * (np.arctan(np.exp(rho_tube + kw/2.0))
                      - np.arctan(np.exp(rho_tube - kw/2.0)))

def build_3d_hessian_gpu(phi_flat):
    diag = 2.0 * d + np.cos(PI * phi_flat)
    idx = np.arange(Ntot, dtype=np.int32)
    ix_arr = idx % N
    iy_arr = (idx // N) % N
    iz_arr = idx // (N * N)
    neighbors = []
    for dx, dy, dz in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
        jx = (ix_arr + dx) % N
        jy = (iy_arr + dy) % N
        jz = (iz_arr + dz) % N
        j = (jz * N * N + jy * N + jx).astype(np.int32)
        neighbors.append(j)
    all_rows = np.concatenate([idx] * 7)
    all_cols = np.concatenate([idx] + neighbors)
    all_vals = np.concatenate([diag] + [-np.ones(Ntot, dtype=np.float64)] * 6)
    H_cpu = sparse.coo_matrix((all_vals, (all_rows, all_cols)),
                               shape=(Ntot, Ntot)).tocsr()
    return csp.csr_matrix(H_cpu)

# ============================================================
# PART 1: SINGLE TORUS — THE CONFINEMENT TEST
# ============================================================
report("PART 1: SINGLE TORUS EIGENSPECTRUM")
report("-" * 55)

phi_torus = make_torus(X, Y, Z, 0.0, R_maj, kink_width)
report(f"Torus profile: phi_max = {phi_torus.max():.4f}")
report(f"Sites with phi > 0.1: {np.sum(phi_torus > 0.1):,}")
report("")

# Build Hessian
report("Building 3D Hessian...")
t0 = time.time()
H_gpu = build_3d_hessian_gpu(phi_torus.ravel())
report(f"  Built in {time.time()-t0:.2f}s")

# Find eigenvalues — request enough to see all sub-gap modes
report("Finding eigenvalues on GPU...")
n_eig = 40
t0 = time.time()
evals_gpu, evecs_gpu = csla.eigsh(H_gpu, k=n_eig, which='SA')
cp.cuda.Stream.null.synchronize()
t_eigsh = time.time() - t0

evals = cp.asnumpy(evals_gpu)
evecs = cp.asnumpy(evecs_gpu)
idx_sort = np.argsort(evals)
evals = evals[idx_sort]
evecs = evecs[:, idx_sort]

report(f"  Eigsh completed in {t_eigsh:.2f}s")
report("")

# THE KEY QUESTION: any negative eigenvalues?
n_negative = np.sum(evals < 0)
n_below_gap = np.sum(evals < 1.0)

report("=" * 50)
if n_negative == 0:
    report("*** NO NEGATIVE EIGENVALUES ***")
    report("*** THE TACHYON IS GONE ***")
    report("*** TOPOLOGY STABILIZES THE KINK-ANTIKINK PAIR ***")
    report("*** THIS IS CONFINEMENT ***")
else:
    report(f"*** {n_negative} NEGATIVE EIGENVALUE(S) FOUND ***")
    report("*** THE TACHYON PERSISTS EVEN ON THE TORUS ***")
report("=" * 50)
report("")

report(f"Eigenvalues below mass gap: {n_below_gap}")
report(f"Negative eigenvalues: {n_negative}")
report("")

# 3D mass gap: on the d=3 cubic lattice with the cosine potential,
# the phonon band starts at omega^2 = 1 (same as 1D, since the
# potential is the same; the kinetic term changes from 2 to 2d
# for the diagonal but the mass gap is still at omega^2 = 1
# because the 2d from the Laplacian is already in the diagonal)
# Actually: in 3D, diagonal = 2d + cos(pi*phi). At phi=0: 2d+1 = 7.
# The phonon dispersion: omega^2 = 1 + 4(sin^2(kx/2)+sin^2(ky/2)+sin^2(kz/2))
# Band bottom: omega^2 = 1 (at k=0). Band top: omega^2 = 1+4*3 = 13.
# So mass gap is still at omega^2 = 1.

report(f"3D mass gap: omega^2 = 1 (phonon band: [1, 13])")
report("")

report(f"{'idx':>4} {'omega^2':>12} {'omega':>10} {'binding':>10} {'type':>10}")
report("-" * 50)
for i in range(min(n_eig, 40)):
    ev = evals[i]
    omega = np.sqrt(abs(ev))
    binding = 1.0 - ev
    if ev < -0.01:
        t = "TACHYON"
    elif ev < 0.01:
        t = "~ZERO"
    elif ev < 1.0:
        t = "BOUND"
    elif ev < 1.01:
        t = "band edge"
    else:
        t = "band"
    report(f"  {i:3d} {ev:+12.6f} {omega:10.6f} {binding:+10.6f} {t:>10}")
    if ev > 1.1 and i > 10:
        break

report("")

# ============================================================
# PART 2: COMPARE 1D vs 3D
# ============================================================
report("PART 2: 1D vs 3D COMPARISON")
report("-" * 55)
report("")
report("1D kink-antikink (from bond_eigenvalue_flow.py):")
report("  mode 0: omega^2 = -0.372  TACHYON")
report("  mode 1: omega^2 = +0.125  BOUND")
report("  mode 2: omega^2 = +0.938  BOUND")
report("")
report("3D torus (this calculation):")
report(f"  Lowest: omega^2 = {evals[0]:+.6f}  {'TACHYON' if evals[0]<0 else 'STABLE'}")
if n_below_gap > 1:
    report(f"  Second: omega^2 = {evals[1]:+.6f}")
if n_below_gap > 2:
    report(f"  Third:  omega^2 = {evals[2]:+.6f}")
report("")

if n_negative == 0:
    report("The 3D torus topology REMOVES the tachyon.")
    report("The kink-antikink annihilation channel is BLOCKED.")
    report("")
    report("Physical interpretation:")
    report("  In 1D: kink can slide into antikink -> annihilation")
    report("  In 3D torus: kink wraps around ring -> can't reach antikink")
    report("  The negative eigenvalue becomes a confined oscillation")
    report("  omega^2 goes from -0.372 (unstable) to >= 0 (stable)")
    report("")
    report("THIS IS CONFINEMENT emerging from lattice topology.")
else:
    report("The tachyon PERSISTS on the torus.")
    report("Possible reasons:")
    report("  - The torus R_maj is too small (kink can still reach antikink)")
    report("  - The kink_width is too large relative to R_maj")
    report("  - The topology doesn't help for this potential shape")
    report("  - The 3D model needs a different kink profile")

report("")

# ============================================================
# PART 3: EIGENVALUE DEGENERACIES
# ============================================================
report("PART 3: DEGENERACY PATTERN")
report("-" * 55)
report("Group eigenvalues by near-degeneracy (within 0.001):")
report("")

groups = []
i = 0
while i < len(evals) and evals[i] < 1.1:
    group = [evals[i]]
    j = i + 1
    while j < len(evals) and abs(evals[j] - evals[i]) < 0.001:
        group.append(evals[j])
        j += 1
    groups.append(group)
    i = j

for g_idx, group in enumerate(groups[:15]):
    avg = np.mean(group)
    t = "TACHYON" if avg < -0.01 else "~ZERO" if avg < 0.01 else \
        "BOUND" if avg < 1.0 else "band"
    report(f"  Group {g_idx}: avg={avg:+.6f}, degeneracy={len(group)}, [{t}]")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report(f"3D Torus (N={N}, R_maj={R_maj}, kw={kink_width}):")
report(f"  Negative eigenvalues: {n_negative}")
report(f"  Below mass gap: {n_below_gap}")
report(f"  Lowest eigenvalue: {evals[0]:+.8f}")
report("")

if n_negative == 0:
    report("RESULT: TACHYON STABILIZED BY TORUS TOPOLOGY")
    report("  The kink-antikink pair is CONFINED on the torus.")
    report("  This is the mechanism for proton stability in GWT.")
else:
    report("RESULT: TACHYON PERSISTS ON TORUS")
    report("  Further investigation needed.")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
