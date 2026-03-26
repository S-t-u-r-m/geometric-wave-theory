"""
R-Dependent W_PI Bond Model — Replacing 8 Corrections with 1 Function
=======================================================================
From the 2D Hessian simulation, the pi/sigma tunnel ratio is:

  W_PI(R) = 0.286 × exp(0.411 × R)    (R in Angstroms)

This is NOT a constant — it depends on bond length. The V8 model
uses W_PI = cos(π/d) = 0.5 as a fixed constant plus 7 other corrections.

KEY FINDING: All V8 correction constants (0.349, 0.5, 0.667, 0.833)
lie on the SAME exponential curve at different bond lengths R.

Test: can W_PI(R) replace all 8 corrections with one function?

Models to test:
  V8:     coupling = σ + π×0.5 - LP + radical + ionic (8 corrections)
  Model A: coupling = 1 + (bo-1) × W(R)  (1 function, uses bo and R)
  Model B: coupling = W(R) × bo  (pure R-dependent)
  Model C: V8 but with W(R) replacing fixed 0.5
  Model D: coupling = W(R)^(1/bo)  (R-dependent, self-consistent)
"""
import sys, io, os, time
import numpy as np
from math import factorial

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3

outfile = os.path.join(os.path.dirname(__file__), "bond_wpi_r_dependent_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("R-DEPENDENT W_PI BOND MODEL")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("")

# ============================================================
# THE W_PI(R) FUNCTION
# ============================================================
# From 2D Hessian: sigma decays at 1.142, pi at 0.925
# W(R_latt) = 0.0068 × exp(0.2175 × R_latt)
# Mapping to Angstroms: R_latt = 17.19 + R_bond/0.529
# Simplifies to: W(R) = 0.286 × exp(0.411 × R)

A_w = 0.286
B_w = 0.411  # per Angstrom

def W_PI_of_R(R_angstrom):
    """Pi/sigma tunnel ratio as a function of bond length."""
    return A_w * np.exp(B_w * R_angstrom)

report("W_PI(R) = 0.286 × exp(0.411 × R)")
report("")
report(f"{'R(A)':>6} {'W_PI(R)':>8} {'GWT const':>12}")
report("-" * 30)
for R in [0.74, 0.92, 1.00, 1.10, 1.20, 1.34, 1.54, 1.89, 2.00, 2.36, 2.67]:
    W = W_PI_of_R(R)
    # Find nearest GWT constant
    consts = {'pi/d^2': PI/d**2, 'cos(pi/d)': 0.5, '(d-1)/d': 2/3,
              '5/6': 5/6, '1': 1.0, '10/27': 10/27, '1/d': 1/3}
    nearest = min(consts, key=lambda k: abs(consts[k] - W))
    report(f"  {R:5.2f}  {W:8.4f}  {nearest} = {consts[nearest]:.4f}")
report("")

# ============================================================
# FUNDAMENTAL CONSTANTS
# ============================================================
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/PI**2 + np.log(2*d)))
E_Ry = (alpha_em**2 / 2) * 0.51100e6 / 1000  # in eV: 13.606 eV
C_BOND = PI / d**2  # 0.3491

# V8 constants for comparison
W_PI_FIXED = np.cos(PI / d)  # 0.5
LP_I = (d**2 + 1) / d**3     # 10/27
F_RAD = (2*d - 1) / (2*d)    # 5/6
C_IONIC = 1 / (2*d + 1)      # 1/7

# ============================================================
# MOLECULE DATABASE
# ============================================================
# (atom_a, atom_b, bo, D_obs_eV, name, radical, R_angstrom)
ATOMS = {
    'H':  (1,  13.598, 1, 0, 0),
    'Li': (3,   5.392, 2, 0, 0),
    'Be': (4,   9.323, 2, 0, 0),
    'B':  (5,   8.298, 2, 1, 0),
    'C':  (6,  11.260, 2, 2, 0),
    'N':  (7,  14.534, 2, 3, 0),
    'O':  (8,  13.618, 2, 4, 1),
    'F':  (9,  17.423, 2, 5, 2),
    'Na': (11,  5.139, 3, 0, 0),
    'S':  (16, 10.360, 3, 4, 1),
    'Cl': (17, 12.968, 3, 5, 2),
    'P':  (15, 10.487, 3, 3, 0),
}
# (Z, IE, n, p_electrons, LP)

molecules = [
    ('H',  'H',  1, 4.478, 'H2',      False, 0.741),
    ('Li', 'Li', 1, 1.046, 'Li2',     False, 2.673),
    ('N',  'N',  3, 9.759, 'N2',      False, 1.098),
    ('O',  'O',  2, 5.116, 'O2',      False, 1.208),
    ('F',  'F',  1, 1.602, 'F2',      False, 1.412),
    ('H',  'F',  1, 5.869, 'HF',      False, 0.917),
    ('H',  'Cl', 1, 4.434, 'HCl',     False, 1.275),
    ('Na', 'Cl', 1, 4.230, 'NaCl',    False, 2.361),
    ('Li', 'H',  1, 2.429, 'LiH',     False, 1.596),
    ('H',  'O',  1, 4.392, 'OH',      True,  0.970),
    ('C',  'O',  3, 11.09, 'CO',      False, 1.128),
    ('N',  'O',  2, 6.497, 'NO',      True,  1.151),
    ('H',  'N',  1, 3.910, 'NH',      True,  1.036),
    ('C',  'H',  1, 4.290, 'CH',      True,  1.120),
    ('C',  'C',  1, 3.600, 'C-C',     False, 1.540),
    ('C',  'N',  3, 7.760, 'CN',      True,  1.172),
    ('C',  'C',  2, 6.360, 'C=C',     False, 1.339),
    ('C',  'O',  2, 7.710, 'C=O',     False, 1.200),
    ('C',  'C',  3, 8.700, 'CC3',     False, 1.203),
    ('N',  'H',  1, 4.513, 'NH(NH3)', False, 1.012),
    ('O',  'H',  1, 4.790, 'OH(H2O)', False, 0.958),
    ('Cl', 'Cl', 1, 2.514, 'Cl2',     False, 1.988),
    ('S',  'H',  1, 3.780, 'SH',      True,  1.340),
    ('S',  'S',  2, 4.370, 'S2',      False, 1.889),
    ('P',  'H',  1, 3.440, 'PH',      True,  1.422),
]

def E_harm(sa, sb):
    IE_a = ATOMS[sa][1]
    IE_b = ATOMS[sb][1]
    return 2 * IE_a * IE_b / (IE_a + IE_b)

def ZPE(D_e, m_a_amu, m_b_amu):
    mu_amu = m_a_amu * m_b_amu / (m_a_amu + m_b_amu)
    mu_me = mu_amu * 1822.89
    D_e_au = D_e / 27.211
    omega = np.sqrt(2 * max(D_e_au, 1e-6) / mu_me)
    return 0.5 * omega * 27.211

mass = {'H': 1.008, 'Li': 6.941, 'C': 12.01, 'N': 14.01, 'O': 16.00,
        'F': 19.00, 'Na': 22.99, 'Cl': 35.45, 'S': 32.07, 'P': 30.97}

# ============================================================
# MODEL A: coupling = 1 + (bo-1) × W(R)
# ============================================================
report("MODEL A: coupling = 1 + (bo-1) × W(R)")
report("-" * 55)
report("W(R) only matters for double/triple bonds (bo > 1).")
report("Single bonds get coupling = 1 regardless of R.")
report("")

header = f"{'Name':>8} {'bo':>3} {'R':>5} {'W(R)':>6} {'cpl':>6} {'D_e':>6} {'D_0':>6} {'obs':>6} {'err':>7}"
report(header)
report("-" * len(header))

errs_A = []
for sa, sb, bo, D_obs, name, rad, R_ang in molecules:
    W = W_PI_of_R(R_ang)
    coupling = 1.0 + (bo - 1) * W
    Eh = E_harm(sa, sb)
    De = C_BOND * Eh * coupling
    zpe = ZPE(De, mass[sa], mass[sb])
    D0 = De - zpe
    err = (D0 - D_obs) / D_obs * 100
    errs_A.append(abs(err))
    report(f"{name:>8} {bo:3d} {R_ang:5.3f} {W:6.3f} {coupling:6.3f} {De:6.3f} {D0:6.3f} {D_obs:6.3f} {err:+7.1f}%")

report(f"\nModel A: mean={np.mean(errs_A):.1f}%, median={np.median(errs_A):.1f}%, max={np.max(errs_A):.1f}%")
report("")

# ============================================================
# MODEL B: coupling = W(R) × bo
# ============================================================
report("MODEL B: coupling = W(R) × bo")
report("-" * 55)

errs_B = []
header = f"{'Name':>8} {'bo':>3} {'R':>5} {'W(R)':>6} {'cpl':>6} {'D_e':>6} {'D_0':>6} {'obs':>6} {'err':>7}"
report(header)
report("-" * len(header))

for sa, sb, bo, D_obs, name, rad, R_ang in molecules:
    W = W_PI_of_R(R_ang)
    coupling = W * bo
    Eh = E_harm(sa, sb)
    De = C_BOND * Eh * coupling
    zpe = ZPE(De, mass[sa], mass[sb])
    D0 = De - zpe
    err = (D0 - D_obs) / D_obs * 100
    errs_B.append(abs(err))
    report(f"{name:>8} {bo:3d} {R_ang:5.3f} {W:6.3f} {coupling:6.3f} {De:6.3f} {D0:6.3f} {D_obs:6.3f} {err:+7.1f}%")

report(f"\nModel B: mean={np.mean(errs_B):.1f}%, median={np.median(errs_B):.1f}%, max={np.max(errs_B):.1f}%")
report("")

# ============================================================
# MODEL C: V8 with W(R) replacing fixed W_PI
# ============================================================
report("MODEL C: V8 formula but W_PI → W(R)")
report("-" * 55)
report("Same 8 corrections as V8, but cos(pi/d) replaced by W(R).")
report("")

errs_C = []
header = f"{'Name':>8} {'bo':>3} {'R':>5} {'W(R)':>6} {'cpl':>6} {'D_e':>6} {'D_0':>6} {'obs':>6} {'err':>7}"
report(header)
report("-" * len(header))

for sa, sb, bo, D_obs, name, rad, R_ang in molecules:
    W = W_PI_of_R(R_ang)
    _, IE_a, n_a, p_a, lp_a = ATOMS[sa]
    _, IE_b, n_b, p_b, lp_b = ATOMS[sb]
    Eh = E_harm(sa, sb)

    n_sigma = 1
    n_pi = bo - 1
    sigma_eff = n_sigma
    pi_eff = n_pi

    # Use W(R) instead of fixed 0.5
    coupling = sigma_eff + pi_eff * W

    # LP repulsion (same as V8)
    n_lp = min(lp_a, lp_b)
    n_max = max(n_a, n_b)
    lp_term = n_lp * LP_I * (2.0 / n_max)**2
    coupling -= lp_term

    # Radical
    if rad:
        coupling *= F_RAD

    # Floor
    coupling = max(coupling, 1/(d+1))

    De = C_BOND * Eh * coupling

    # Ionic
    delta_IE = abs(IE_a - IE_b)
    E_avg = (IE_a + IE_b) / 2
    asym = delta_IE / E_avg if E_avg > 0 else 0
    D_ionic = 0
    if asym > 0.1:
        D_ionic = C_IONIC * delta_IE
    De += D_ionic

    zpe = ZPE(De, mass[sa], mass[sb])
    D0 = De - zpe
    err = (D0 - D_obs) / D_obs * 100
    errs_C.append(abs(err))
    report(f"{name:>8} {bo:3d} {R_ang:5.3f} {W:6.3f} {coupling:6.3f} {De:6.3f} {D0:6.3f} {D_obs:6.3f} {err:+7.1f}%")

report(f"\nModel C: mean={np.mean(errs_C):.1f}%, median={np.median(errs_C):.1f}%, max={np.max(errs_C):.1f}%")
report("")

# ============================================================
# MODEL D: ALL corrections from W(R) at different scales
# ============================================================
report("MODEL D: ALL corrections as W(R) at different R-scales")
report("-" * 55)
report("Hypothesis: each V8 correction IS W_PI(R) at a different scale.")
report("LP repulsion = W(R - delta_LP) where delta_LP shifts to shorter R.")
report("Radical = W(R + delta_rad) where delta_rad shifts to longer R.")
report("")

# The idea: the 8 corrections are the SAME function sampled at
# different effective separations. The molecule's actual bond length
# determines the base R, and each physical effect (LP, radical, etc.)
# shifts R by a characteristic amount.
#
# From our data:
# LP_I = 0.370 corresponds to R_latt = 18.37, gap = 2.37
# W_PI = 0.500 corresponds to R_latt = 19.75, gap = 3.75
# F_RAD = 0.833 corresponds to R_latt = 22.10, gap = 6.10
#
# Delta from the W_PI=0.5 point:
# LP:  gap = 2.37 vs 3.75 → delta = -1.38 lattice units = -0.73 Angstrom
# RAD: gap = 6.10 vs 3.75 → delta = +2.35 lattice units = +1.24 Angstrom

delta_LP = -0.73   # LP is a shorter-range effect
delta_RAD = +1.24  # radical is a longer-range effect

errs_D = []
header = f"{'Name':>8} {'bo':>3} {'R':>5} {'W_base':>6} {'W_lp':>6} {'W_rad':>6} {'D_0':>6} {'obs':>6} {'err':>7}"
report(header)
report("-" * len(header))

for sa, sb, bo, D_obs, name, rad, R_ang in molecules:
    _, IE_a, n_a, p_a, lp_a = ATOMS[sa]
    _, IE_b, n_b, p_b, lp_b = ATOMS[sb]
    Eh = E_harm(sa, sb)

    W_base = W_PI_of_R(R_ang)
    W_lp = W_PI_of_R(R_ang + delta_LP)    # shorter range → smaller W
    W_rad = W_PI_of_R(R_ang + delta_RAD)   # longer range → larger W

    n_sigma = 1
    n_pi = bo - 1

    coupling = 1.0 + n_pi * W_base

    # LP repulsion uses W at shorter range
    n_lp = min(lp_a, lp_b)
    n_max = max(n_a, n_b)
    lp_term = n_lp * W_lp * (2.0 / n_max)**2
    coupling -= lp_term

    # Radical uses W at longer range
    if rad:
        coupling *= W_rad

    coupling = max(coupling, 1/(d+1))

    De = C_BOND * Eh * coupling

    # Ionic
    delta_IE = abs(IE_a - IE_b)
    E_avg = (IE_a + IE_b) / 2
    asym = delta_IE / E_avg if E_avg > 0 else 0
    if asym > 0.1:
        De += C_IONIC * delta_IE

    zpe = ZPE(De, mass[sa], mass[sb])
    D0 = De - zpe
    err = (D0 - D_obs) / D_obs * 100
    errs_D.append(abs(err))
    report(f"{name:>8} {bo:3d} {R_ang:5.3f} {W_base:6.3f} {W_lp:6.3f} {W_rad:6.3f} {D0:6.3f} {D_obs:6.3f} {err:+7.1f}%")

report(f"\nModel D: mean={np.mean(errs_D):.1f}%, median={np.median(errs_D):.1f}%, max={np.max(errs_D):.1f}%")
report("")

# ============================================================
# COMPARISON
# ============================================================
report("COMPARISON")
report("=" * 70)
report("")
report(f"{'Model':>10} {'mean%':>7} {'median%':>8} {'max%':>7} {'<5%':>5} {'<10%':>5}")
report("-" * 45)
for name, errs in [('A (simple)', errs_A), ('B (W×bo)', errs_B),
                    ('C (V8+W(R))', errs_C), ('D (all W(R))', errs_D)]:
    n5 = sum(1 for e in errs if e < 5)
    n10 = sum(1 for e in errs if e < 10)
    report(f"{name:>10} {np.mean(errs):7.1f} {np.median(errs):8.1f} {np.max(errs):7.1f} "
           f"{n5:3d}/25 {n10:3d}/25")

report("")
report("V8 reference: mean=1.7%, median=1.5%, max=4.8%, 25/25 under 5%")
report("")

# ============================================================
# WHICH MOLECULES BENEFIT MOST FROM R-DEPENDENCE?
# ============================================================
report("MOLECULE-BY-MOLECULE COMPARISON (Model C vs V8 target):")
report("-" * 55)

report(f"{'Name':>8} {'R':>5} {'W(R)':>6} {'err_C':>8} {'note':>20}")
report("-" * 50)
for i, (sa, sb, bo, D_obs, name, rad, R_ang) in enumerate(molecules):
    W = W_PI_of_R(R_ang)
    err = errs_C[i]
    note = ""
    if W < 0.45: note = "W < cos(pi/d)"
    elif W > 0.55: note = "W > cos(pi/d)"
    if bo > 1: note += " pi-bond"
    report(f"{name:>8} {R_ang:5.3f} {W:6.3f} {err:+8.1f}% {note:>20}")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
log.close()
print(f'\nResults saved to: {outfile}')
