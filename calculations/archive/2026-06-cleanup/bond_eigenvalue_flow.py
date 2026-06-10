"""
Bond Eigenvalue Flow — Phase 1: Raw Physics, No Assumptions
=============================================================
Two kink-antikink pairs on the 1D lattice.
Compute ALL eigenvalues below and near the mass gap at every integer R.
No interpretation. Just the data.

Questions this answers:
  - How many eigenvalues exist below the mass gap?
  - How do they move as R changes?
  - Do any cross, repel, or merge?
  - Does mode 0 (tachyon in single well) remain tachyonic in the double well?
  - At what R does each mode's splitting become significant?
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3

outfile = os.path.join(os.path.dirname(__file__), "bond_eigenvalue_flow_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("BOND EIGENVALUE FLOW — RAW PHYSICS")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("")

N = 512
x = np.arange(N, dtype=np.float64)
center = N // 2
kw = 3

def kink(x, pos):
    return (4.0/PI) * (np.arctan(np.exp(x - pos + kw/2))
                      - np.arctan(np.exp(x - pos - kw/2)))

def get_eigenvalues(phi, n_eig=20):
    diag = 2.0 + np.cos(PI * phi)
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    ev, evec = eigsh(H.tocsr(), k=n_eig, which='SM')
    idx = np.argsort(ev)
    return ev[idx], evec[:, idx]

# ============================================================
# PART 1: SINGLE WELL — REFERENCE SPECTRUM
# ============================================================
report("PART 1: SINGLE KINK-ANTIKINK WELL")
report("-" * 55)

phi_single = kink(x, center)
ev_single, evec_single = get_eigenvalues(phi_single, n_eig=20)

n_below_gap = np.sum(ev_single < 1.0)
report(f"Eigenvalues below mass gap (omega^2 < 1): {n_below_gap}")
report("")

report(f"{'idx':>4} {'omega^2':>12} {'omega':>10} {'binding':>10} {'type':>10}")
report("-" * 50)
for i in range(min(12, len(ev_single))):
    ev = ev_single[i]
    omega = np.sqrt(abs(ev)) * (1 if ev >= 0 else -1)  # signed
    binding = 1.0 - ev
    if ev < 0:
        t = "TACHYON"
    elif ev < 1.0:
        t = "BOUND"
    elif ev < 1.01:
        t = "band edge"
    else:
        t = "band"
    report(f"  {i:3d} {ev:+12.6f} {omega:+10.6f} {binding:+10.6f} {t:>10}")

report("")

# ============================================================
# PART 2: DOUBLE WELL — FULL EIGENVALUE FLOW
# ============================================================
report("PART 2: TWO KINK-ANTIKINK WELLS — EIGENVALUE FLOW")
report("-" * 55)
report("Every integer R from 4 to 60. All eigenvalues below omega^2 = 1.1")
report("")

R_range = list(range(4, 61))
n_eig = 20  # enough to capture all sub-gap modes

# Store everything
all_evals = {}  # {R: array of eigenvalues}
all_evecs = {}  # {R: matrix of eigenvectors}
n_below = {}    # {R: count below gap}

for R in R_range:
    pA = center - R // 2
    pB = center + R // 2
    phi_double = kink(x, pA) + kink(x, pB)

    ev, evec = get_eigenvalues(phi_double, n_eig=n_eig)
    all_evals[R] = ev
    all_evecs[R] = evec
    n_below[R] = np.sum(ev < 1.0)

# Print the flow table
# How many eigenvalues below gap changes with R?
header = f"{'R':>3} {'n<gap':>5} "
for i in range(10):
    header += f"{'ev_'+str(i):>11} "
report(header)
report("-" * (8 + 10 * 12))

for R in R_range:
    ev = all_evals[R]
    line = f"{R:3d} {n_below[R]:5d} "
    for i in range(min(10, len(ev))):
        line += f"{ev[i]:+11.6f} "
    report(line)

report("")

# ============================================================
# PART 3: TRACK SPECIFIC EIGENVALUE PAIRS
# ============================================================
report("PART 3: EIGENVALUE PAIR TRACKING")
report("-" * 55)
report("At large R, each single-well eigenvalue appears twice (degenerate).")
report("As R decreases, they split into bonding + antibonding.")
report("Track the splitting for each pair.")
report("")

# At R=60 (large), eigenvalues should be ~doubly degenerate
ev_large = all_evals[60]
report(f"At R=60 (isolated wells):")
for i in range(min(8, len(ev_large))):
    report(f"  ev[{i}] = {ev_large[i]:+.8f}")
report("")

# Identify pairs at large R
pairs = []
i = 0
while i < len(ev_large) - 1:
    if abs(ev_large[i+1] - ev_large[i]) < 0.01:
        pairs.append((i, i+1))
        i += 2
    else:
        pairs.append((i, None))  # unpaired
        i += 1

report(f"Identified {len(pairs)} eigenvalue pairs at R=60:")
for p_idx, (i, j) in enumerate(pairs):
    if j is not None:
        avg = (ev_large[i] + ev_large[j]) / 2
        split = ev_large[j] - ev_large[i]
        report(f"  Pair {p_idx}: ev[{i},{j}] avg={avg:+.6f} split={split:.2e}")
    else:
        report(f"  Pair {p_idx}: ev[{i}] unpaired = {ev_large[i]:+.6f}")
    if p_idx >= 5:
        break

report("")

# Track splittings vs R for each pair
report("SPLITTING vs R for each pair:")
header2 = f"{'R':>3} "
for p_idx in range(min(4, len(pairs))):
    header2 += f"{'pair'+str(p_idx)+'_avg':>12} {'split'+str(p_idx):>12} "
report(header2)
report("-" * (3 + 4 * 24))

for R in R_range:
    ev = all_evals[R]
    line = f"{R:3d} "
    for p_idx, (i, j) in enumerate(pairs[:4]):
        if j is not None and j < len(ev):
            avg = (ev[i] + ev[j]) / 2
            split = ev[j] - ev[i]
            line += f"{avg:+12.6f} {split:12.8f} "
        elif i < len(ev):
            line += f"{ev[i]:+12.6f} {'---':>12} "
        else:
            line += f"{'---':>12} {'---':>12} "
    report(line)

report("")

# ============================================================
# PART 4: THE TACHYON QUESTION
# ============================================================
report("PART 4: DOES THE TACHYON STABILIZE?")
report("-" * 55)
report("Single well: mode 0 has omega^2 = -0.372 (TACHYON).")
report("In the double well, does the lowest eigenvalue remain negative?")
report("If it becomes positive at some R, the tachyon 'heals' in pairs.")
report("")

report(f"{'R':>3} {'ev_lowest':>12} {'status':>10}")
report("-" * 28)

for R in R_range:
    ev = all_evals[R]
    lowest = ev[0]
    if lowest < -0.01:
        status = "TACHYON"
    elif lowest < 0.01:
        status = "~ZERO"
    else:
        status = "STABLE"
    if R <= 20 or R % 10 == 0 or abs(lowest) < 0.05:
        report(f"{R:3d} {lowest:+12.6f} {status:>10}")

report("")

# Is there a crossing where the lowest eigenvalue changes sign?
for i in range(len(R_range) - 1):
    R1, R2 = R_range[i], R_range[i+1]
    ev1 = all_evals[R1][0]
    ev2 = all_evals[R2][0]
    if ev1 * ev2 < 0:  # sign change
        report(f"*** SIGN CHANGE: ev[0] goes from {ev1:+.6f} at R={R1} "
               f"to {ev2:+.6f} at R={R2} ***")

report("")

# ============================================================
# PART 5: EIGENVECTOR LOCALIZATION
# ============================================================
report("PART 5: WHERE DO THE EIGENVECTORS LIVE?")
report("-" * 55)
report("For each eigenvalue at R=10, compute the weight on well A vs well B.")
report("")

R_test = 10
pA = center - R_test // 2
pB = center + R_test // 2
ev_test = all_evals[R_test]
evec_test = all_evecs[R_test]

# Define well regions: well A centered at pA, well B at pB
hw = 10  # half-width of well region
mask_A = np.zeros(N, dtype=bool)
mask_A[max(0, pA-hw):min(N, pA+hw+1)] = True
mask_B = np.zeros(N, dtype=bool)
mask_B[max(0, pB-hw):min(N, pB+hw+1)] = True
mask_mid = ~mask_A & ~mask_B  # between and outside wells

report(f"R = {R_test}, wells at {pA} and {pB}")
report(f"Well A region: [{pA-hw}, {pA+hw}]")
report(f"Well B region: [{pB-hw}, {pB+hw}]")
report("")

report(f"{'idx':>4} {'omega^2':>10} {'wt_A':>8} {'wt_B':>8} {'wt_mid':>8} "
       f"{'A-B':>8} {'type':>12}")
report("-" * 60)

for i in range(min(10, len(ev_test))):
    psi = evec_test[:, i]
    psi2 = psi**2
    wt_A = np.sum(psi2[mask_A])
    wt_B = np.sum(psi2[mask_B])
    wt_mid = np.sum(psi2[mask_mid])
    asym = wt_A - wt_B

    if abs(asym) < 0.1:
        if wt_A + wt_B > 0.5:
            loc = "BONDING" if i % 2 == 0 else "ANTIBOND"
        else:
            loc = "delocalized"
    else:
        loc = "well_A" if asym > 0 else "well_B"

    report(f"  {i:3d} {ev_test[i]:+10.6f} {wt_A:8.4f} {wt_B:8.4f} {wt_mid:8.4f} "
           f"{asym:+8.4f} {loc:>12}")

report("")

# ============================================================
# PART 6: SAME FOR R=6 (close range)
# ============================================================
report("PART 6: EIGENVECTOR LOCALIZATION AT R=6 (close range)")
report("-" * 55)

R_test2 = 6
pA2 = center - R_test2 // 2
pB2 = center + R_test2 // 2
ev_test2 = all_evals[R_test2]
evec_test2 = all_evecs[R_test2]

mask_A2 = np.zeros(N, dtype=bool)
mask_A2[max(0, pA2-hw):min(N, pA2+hw+1)] = True
mask_B2 = np.zeros(N, dtype=bool)
mask_B2[max(0, pB2-hw):min(N, pB2+hw+1)] = True

report(f"R = {R_test2}, wells at {pA2} and {pB2}")
report(f"{'idx':>4} {'omega^2':>10} {'wt_A':>8} {'wt_B':>8} {'A-B':>8} {'type':>12}")
report("-" * 55)

for i in range(min(10, len(ev_test2))):
    psi = evec_test2[:, i]
    psi2 = psi**2
    wt_A = np.sum(psi2[mask_A2])
    wt_B = np.sum(psi2[mask_B2])
    asym = wt_A - wt_B

    if abs(asym) < 0.1:
        loc = "symmetric" if i % 2 == 0 else "antisym"
    else:
        loc = "well_A" if asym > 0 else "well_B"

    report(f"  {i:3d} {ev_test2[i]:+10.6f} {wt_A:8.4f} {wt_B:8.4f} "
           f"{asym:+8.4f} {loc:>12}")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY — RAW OBSERVATIONS (NO INTERPRETATION)")
report("=" * 70)
report("")
report(f"Single well: {n_below_gap} eigenvalues below mass gap")
report(f"  Lowest: omega^2 = {ev_single[0]:+.6f} ({'TACHYON' if ev_single[0] < 0 else 'STABLE'})")
report("")

# Does the tachyon persist at all R?
all_lowest = [all_evals[R][0] for R in R_range]
min_lowest = min(all_lowest)
max_lowest = max(all_lowest)
report(f"Double well lowest eigenvalue range: [{min_lowest:+.6f}, {max_lowest:+.6f}]")
if all(v < 0 for v in all_lowest):
    report("  TACHYON persists at ALL separations R=4..60")
elif all(v >= 0 for v in all_lowest):
    report("  TACHYON is GONE in the double well (stabilized by pairing)")
else:
    report("  TACHYON transitions: negative at some R, positive at others")

report("")

# How many sub-gap eigenvalues at different R?
report("Sub-gap eigenvalue count vs R:")
for R in [4, 6, 8, 10, 15, 20, 30, 60]:
    report(f"  R={R:3d}: {n_below[R]} eigenvalues below mass gap")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
