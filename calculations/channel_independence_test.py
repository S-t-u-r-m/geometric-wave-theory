"""
Channel Independence Test — Are Covalent, Ionic, and Occupancy Separable?
===========================================================================
Three channels affect bond energy:
  1. GEOMETRIC: kink-kink overlap at distance R
  2. OCCUPANCY: number of breather modes in each well (electron config)
  3. ASYMMETRY: different well depths (different IE → ionic character)

Test: vary each independently and check if effects are additive.
  If D_e(all) = D_e(geom) + D_e(occ) + D_e(asym) → SEPARABLE
  If D_e(all) ≠ sum → COUPLED (and by how much?)

Method: 1D Hessian with two kink wells. Modify:
  - R: separation between wells
  - n_occ: number of breather perturbations in each well (0 to 7)
  - delta_V: asymmetric potential depth (shifts cos(pi*phi) in one well)
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)

outfile = os.path.join(os.path.dirname(__file__), "channel_independence_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("CHANNEL INDEPENDENCE TEST")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("")

# ============================================================
# SETUP
# ============================================================
N = 512
x = np.arange(N, dtype=np.float64)
center = N // 2

def make_proton(x, pos, kw=3):
    return (4.0/PI) * (np.arctan(np.exp(x - pos + kw/2))
                      - np.arctan(np.exp(x - pos - kw/2)))

def breather_sq(x, pos, eps):
    phi = (4.0/PI) * np.arctan(1.0 / np.cosh(eps * (x - pos) + 1e-30))
    return phi**2

def build_hessian(phi_bg, extra_diag=None):
    diag = 2.0 + np.cos(PI * phi_bg)
    if extra_diag is not None:
        diag = diag + extra_diag
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()

def get_E0(phi_bg, extra_diag=None):
    H = build_hessian(phi_bg, extra_diag)
    ev, _ = eigsh(H, k=1, which='SM')
    return ev[0]

def breather_pert(x, pos, mode_list):
    """Combined perturbation from multiple breather modes at a position."""
    pert = np.zeros(N)
    for n in mode_list:
        eps_n = np.sin(n * gamma)
        pert += -(PI**2 / 6.0) * breather_sq(x, pos, eps_n)
    return pert

def asymmetry_pert(x, pos, strength, width=8):
    """Asymmetric potential: shifts the well depth at one position.
    Adds 'strength' to the Hessian diagonal near 'pos'.
    Positive = stiffer well (higher IE), negative = softer (lower IE).
    """
    return strength * np.exp(-(x - pos)**2 / (2 * width**2))

# ============================================================
# REFERENCE: symmetric, no breathers, no asymmetry
# ============================================================
report("REFERENCE CONFIGURATIONS")
report("-" * 55)

# Single proton reference
phi_single = make_proton(x, center)
E0_single_bare = get_E0(phi_single)
report(f"Single proton (bare): E0 = {E0_single_bare:.8f}")

# Single proton with 3 breather modes
modes_3 = [1, 2, 3]
pert_3 = breather_pert(x, center, modes_3)
E0_single_3 = get_E0(phi_single, pert_3)
report(f"Single proton (3 modes): E0 = {E0_single_3:.8f}")

# Single proton with 7 breather modes
modes_7 = [1, 2, 3, 4, 5, 6, 7]
pert_7 = breather_pert(x, center, modes_7)
E0_single_7 = get_E0(phi_single, pert_7)
report(f"Single proton (7 modes): E0 = {E0_single_7:.8f}")
report("")

# ============================================================
# THE THREE CHANNELS
# ============================================================
R_test = 10
pos_A = center - R_test // 2
pos_B = center + R_test // 2
phi_double = make_proton(x, pos_A) + make_proton(x, pos_B)

report(f"Test separation: R = {R_test}")
report("")

# --- CHANNEL 1: GEOMETRIC (just the bare kink overlap) ---
E0_geom = get_E0(phi_double)
V_geom = E0_geom - E0_single_bare
report(f"CHANNEL 1 — GEOMETRIC (bare kinks at R={R_test}):")
report(f"  V_geom = {V_geom:+.8f}")
report("")

# --- CHANNEL 2: OCCUPANCY (add breather modes symmetrically) ---
# Both wells have the same number of modes
for n_modes_test in [1, 3, 5, 7]:
    mode_list = list(range(1, n_modes_test + 1))
    pert_A = breather_pert(x, pos_A, mode_list)
    pert_B = breather_pert(x, pos_B, mode_list)
    pert_occ = pert_A + pert_B

    # Also need single-well reference with same modes
    pert_single_occ = breather_pert(x, center, mode_list)
    E0_ref_occ = get_E0(phi_single, pert_single_occ)

    E0_occ = get_E0(phi_double, pert_occ)
    V_occ = E0_occ - E0_ref_occ

    report(f"CHANNEL 2 — OCCUPANCY ({n_modes_test} modes per well):")
    report(f"  V_occ = {V_occ:+.8f}, V_geom = {V_geom:+.8f}, "
           f"difference = {V_occ - V_geom:+.8f}")

report("")

# --- CHANNEL 3: ASYMMETRY (different well depths) ---
# Shift one well's potential while keeping the other the same
for delta_V_strength in [0.1, 0.3, 0.5, 1.0]:
    asym_A = asymmetry_pert(x, pos_A, +delta_V_strength)  # stiffer
    asym_B = asymmetry_pert(x, pos_B, -delta_V_strength)  # softer

    E0_asym = get_E0(phi_double, asym_A + asym_B)
    V_asym = E0_asym - E0_single_bare  # rough reference

    # Symmetric version for comparison
    E0_sym = get_E0(phi_double, asymmetry_pert(x, pos_A, 0) + asymmetry_pert(x, pos_B, 0))

    report(f"CHANNEL 3 — ASYMMETRY (delta_V = ±{delta_V_strength}):")
    report(f"  V_asym = {V_asym:+.8f}, V_sym = {E0_sym - E0_single_bare:+.8f}, "
           f"ionic_shift = {V_asym - (E0_sym - E0_single_bare):+.8f}")

report("")

# ============================================================
# THE INDEPENDENCE TEST — ALL COMBINATIONS
# ============================================================
report("INDEPENDENCE TEST — FACTORIAL DESIGN")
report("-" * 55)
report("Compare: V(all three) vs V(geom) + V(occ) + V(asym)")
report("")

# Use moderate values for each channel
modes_occ = [1, 2, 3]  # 3 breather modes
delta_V = 0.3           # moderate asymmetry

# A: geometric only (bare kinks)
E0_A = get_E0(phi_double)
V_A = E0_A - E0_single_bare

# B: geometric + occupancy
pert_occ_AB = breather_pert(x, pos_A, modes_occ) + breather_pert(x, pos_B, modes_occ)
E0_ref_occ = get_E0(phi_single, breather_pert(x, center, modes_occ))
E0_B = get_E0(phi_double, pert_occ_AB)
V_B = E0_B - E0_ref_occ

# C: geometric + asymmetry
pert_asym = asymmetry_pert(x, pos_A, +delta_V) + asymmetry_pert(x, pos_B, -delta_V)
E0_C = get_E0(phi_double, pert_asym)
# Reference: single well with average asymmetry (none)
V_C = E0_C - E0_single_bare

# D: occupancy only (no geometric — infinite R approximation)
# Use same-well occupancy difference: one well has modes, other doesn't
# This tests whether occupancy creates a force even without geometric overlap
pert_occ_A_only = breather_pert(x, pos_A, modes_occ)
E0_D = get_E0(phi_double, pert_occ_A_only)
V_D = E0_D - E0_single_bare

# E: asymmetry only (at the base geometric R, no extra occupancy)
V_E = V_C - V_A  # subtract the geometric part

# F: ALL THREE — geometric + occupancy + asymmetry
pert_all = pert_occ_AB + pert_asym
E0_F = get_E0(phi_double, pert_all)
V_F = E0_F - E0_ref_occ

# Additive prediction: V_F ≈ V_B + V_E?
# Or: V_F ≈ V_A + (V_B - V_A) + (V_C - V_A)?
V_additive = V_A + (V_B - V_A) + (V_C - V_A)
nonadd = V_F - V_additive
nonadd_pct = nonadd / abs(V_F) * 100 if abs(V_F) > 1e-12 else 0

report(f"  V(geometric only):          {V_A:+.8f}")
report(f"  V(geometric + occupancy):   {V_B:+.8f}")
report(f"  V(geometric + asymmetry):   {V_C:+.8f}")
report(f"  V(ALL THREE):               {V_F:+.8f}")
report(f"  V(additive prediction):     {V_additive:+.8f}")
report(f"  Non-additive residual:      {nonadd:+.8f} ({nonadd_pct:+.1f}%)")
report("")

if abs(nonadd_pct) < 5:
    report("  → Channels are approximately INDEPENDENT (<5% coupling)")
elif abs(nonadd_pct) < 20:
    report("  → Channels have WEAK coupling (5-20%)")
else:
    report(f"  → Channels are COUPLED ({abs(nonadd_pct):.0f}% non-additive)")
report("")

# ============================================================
# SCAN R — HOW DOES COUPLING CHANGE WITH DISTANCE?
# ============================================================
report("R-DEPENDENCE OF CHANNEL COUPLING")
report("-" * 55)
report("Does the coupling between channels change with separation?")
report("")

report(f"{'R':>4} {'V_geom':>12} {'V_g+occ':>12} {'V_g+asym':>12} "
       f"{'V_all':>12} {'V_add':>12} {'nonadd%':>9}")
report("-" * 75)

for R in [6, 8, 10, 12, 14, 16, 20, 30]:
    pA = center - R // 2
    pB = center + R // 2
    phi2 = make_proton(x, pA) + make_proton(x, pB)

    # Geometric only
    ea = get_E0(phi2)
    va = ea - E0_single_bare

    # Geometric + occupancy
    po = breather_pert(x, pA, modes_occ) + breather_pert(x, pB, modes_occ)
    eb = get_E0(phi2, po)
    vb = eb - E0_ref_occ

    # Geometric + asymmetry
    pa = asymmetry_pert(x, pA, +delta_V) + asymmetry_pert(x, pB, -delta_V)
    ec = get_E0(phi2, pa)
    vc = ec - E0_single_bare

    # All three
    ef = get_E0(phi2, po + pa)
    vf = ef - E0_ref_occ

    v_add = va + (vb - va) + (vc - va)
    na = vf - v_add
    na_pct = na / abs(vf) * 100 if abs(vf) > 1e-10 else 0

    report(f"{R:4d} {va:+12.8f} {vb:+12.8f} {vc:+12.8f} "
           f"{vf:+12.8f} {v_add:+12.8f} {na_pct:+9.1f}%")

report("")

# ============================================================
# ASYMMETRIC OCCUPANCY — THE IONIC-COVALENT COUPLING
# ============================================================
report("ASYMMETRIC OCCUPANCY — IONIC-COVALENT COUPLING")
report("-" * 55)
report("What happens when one well has more modes than the other?")
report("This is the analog of ionic bonding: one atom has more electrons.")
report("")

R_ionic = 10
pA = center - R_ionic // 2
pB = center + R_ionic // 2
phi2 = make_proton(x, pA) + make_proton(x, pB)

report(f"R = {R_ionic}")
report(f"{'n_A':>4} {'n_B':>4} {'V(R)':>12} {'vs_sym':>12} {'note':>20}")
report("-" * 58)

# Symmetric reference: both have 3 modes
pert_sym = breather_pert(x, pA, [1,2,3]) + breather_pert(x, pB, [1,2,3])
E0_ref_sym = get_E0(phi_single, breather_pert(x, center, [1,2,3]))
E0_sym_ref = get_E0(phi2, pert_sym)
V_sym = E0_sym_ref - E0_ref_sym

for n_A, n_B in [(0,0), (1,1), (2,2), (3,3), (4,4), (5,5), (7,7),
                  (1,3), (2,4), (3,5), (1,5), (0,3), (0,7),
                  (3,1), (5,1), (7,0)]:
    modes_A = list(range(1, n_A + 1))
    modes_B = list(range(1, n_B + 1))

    pert = breather_pert(x, pA, modes_A) + breather_pert(x, pB, modes_B)

    # Reference: average of the two single-well energies
    if n_A > 0:
        E0_refA = get_E0(phi_single, breather_pert(x, center, modes_A))
    else:
        E0_refA = E0_single_bare
    if n_B > 0:
        E0_refB = get_E0(phi_single, breather_pert(x, center, modes_B))
    else:
        E0_refB = E0_single_bare
    E0_ref_avg = (E0_refA + E0_refB) / 2

    E0_double = get_E0(phi2, pert)
    V_R = E0_double - E0_ref_avg

    vs_sym = V_R - V_sym
    note = "symmetric" if n_A == n_B else f"asym {n_A}-{n_B}"
    if n_A == 3 and n_B == 3:
        note = "REFERENCE"
    report(f"  {n_A:3d} {n_B:3d} {V_R:+12.8f} {vs_sym:+12.8f} {note:>20}")

report("")

# ============================================================
# SUMMARY
# ============================================================
report("SUMMARY")
report("=" * 70)
report("")
report("Three bond channels tested for independence:")
report("  1. GEOMETRIC (kink overlap at distance R)")
report("  2. OCCUPANCY (breather modes in each well)")
report("  3. ASYMMETRY (different well depths / IE)")
report("")
report("KEY QUESTION: V(all three) = sum of individual channels?")
report("")
report("If YES → can model each channel separately (V8 approach is valid)")
report("If NO → must model them as a coupled system")
report("        (explains why V8 needs molecule-specific corrections)")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f"\nResults saved to: {outfile}")
