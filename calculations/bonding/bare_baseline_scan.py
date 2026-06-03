#!/usr/bin/env python3
"""
Bare-baseline scan on V8's full set.

The hypothesis to test: is bond_3d_emerge.py's 4.73 eV match on H2 a
universal feature of the bare 1D Hessian method, or is H2 a special case
and V8's 8-correction stack is what makes everything else work?

Three predictions per bond:
  D_bare    = bare 1D sigma-only Hessian (no corrections, no multi-bond)
  D_bare_BO = bare + naive bond-order scaling (1 + (bo-1)*W_PI)
  D_V8      = V8's full 8-correction prediction
  D_obs     = observed

We then look at signed errors. The diagnostic question:
  - Bare 1D matches H2 to 0.6%. Does it match other single bonds equally well?
  - If yes -> bond_3d_emerge.py is universally accurate at the sigma-channel
    level; V8's corrections are addressing multi-bond + LP + ionic, not a
    flaw in the base.
  - If no -> the bare match on H2 is partly coincidence, and the collective
    solver has to do work even for simple single bonds.

This guides v2's design: which V8 corrections are essential (must be
preserved in the collective solver) vs which are residual modeling errors
that the collective solver replaces from scratch.
"""
import sys
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

PI = np.pi
d = 3
S_PT = (-1 + np.sqrt(1 + 8/PI**2)) / 2     # 0.17279
KW = 3
W_PI = np.cos(PI/d)                          # 0.5

# ----- Reuse V8 atom database and V8 predictor -----
import importlib.util
spec = importlib.util.spec_from_file_location(
    "v8", "calculations/bonding/bond_v8_full.py")
v8 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v8)
ATOMS = v8.ATOMS
bond_v8 = v8.bond_v8

# Mirror V8's test_bonds verbatim
TEST_BONDS = [
    ('H','H',1,4.478,'H2',False),
    ('Li','Li',1,1.046,'Li2',False),
    ('N','N',3,9.759,'N2',False),
    ('O','O',2,5.116,'O2',False),
    ('F','F',1,1.602,'F2',False),
    ('H','F',1,5.869,'HF',False),
    ('H','Cl',1,4.434,'HCl',False),
    ('Na','Cl',1,4.230,'NaCl',False),
    ('Li','H',1,2.429,'LiH',False),
    ('H','O',1,4.392,'OH',True),
    ('C','O',3,11.09,'CO',False),
    ('N','O',2,6.497,'NO',True),
    ('H','N',1,3.910,'NH',True),
    ('C','H',1,4.290,'CH',True),
    ('C','C',1,3.600,'C-C',False),
    ('C','N',3,7.760,'CN',True),
    ('C','C',2,6.360,'C=C',False),
    ('C','O',2,7.710,'C=O',False),
    ('C','C',3,8.700,'CC3',False),
    ('N','H',1,4.513,'NH(NH3)',False),
    ('O','H',1,4.790,'OH(H2O)',False),
    ('Cl','Cl',1,2.514,'Cl2',False),
    ('S','H',1,3.780,'SH',True),
    ('S','S',2,4.370,'S2',False),
    ('P','H',1,3.440,'PH',True),
]

# ============================================================
# 1D BARE HESSIAN — universal D_e_lattice (Z-independent kinks)
# ============================================================
def kink(x, x0):
    return (4.0/PI) * np.arctan(np.exp(x - x0))
def antikink(x, x0):
    return 2.0 - (4.0/PI) * np.arctan(np.exp(x - x0))
def proton(x, center, w=KW):
    return kink(x, center - w/2) + antikink(x, center + w/2) - 2.0
def two_protons(x, R, w=KW):
    c = len(x)//2
    return proton(x, c - R/2, w) + proton(x, c + R/2, w)

def build_hess(phi):
    N = len(phi)
    diag = 2.0 + np.cos(PI*phi)
    off = -np.ones(N-1)
    H = sparse.diags([off, diag, off], [-1,0,1], shape=(N,N), format='lil')
    H[0,N-1] = -1.0
    H[N-1,0] = -1.0
    return H.tocsr()

def lowest_eig(phi, k=4):
    H = build_hess(phi)
    w2, _ = eigsh(H, k=k, which='SM')
    return float(np.min(w2))

def bare_De_lattice(N=256, R_scan=range(6,15)):
    """One-shot bare 1D Hessian baseline (Z-independent).
    R>=6 avoids the bump-merger artifact at R<=5 where the two protons
    coalesce into a single big soliton (gives unphysically large D_e)."""
    x = np.arange(N, dtype=np.float64)
    E_iso = lowest_eig(proton(x, N/2))
    V = []
    for R in R_scan:
        E_R = lowest_eig(two_protons(x, R))
        V.append((R, E_R - E_iso))
    R_arr = np.array([v[0] for v in V])
    V_arr = np.array([v[1] for v in V])
    i_min = int(np.argmin(V_arr))
    return -float(V_arr[i_min]), int(R_arr[i_min])


# ============================================================
# SCAN
# ============================================================
print("=" * 78)
print("BARE 1D BASELINE SCAN — is bond_3d_emerge.py's H2 match universal?")
print("=" * 78)

# Precompute the universal sigma-channel lattice value (the Hessian is
# Z-independent; only E_harm changes per molecule).
print("Computing universal bare D_e_lattice (one-shot)...")
De_lat_universal, R_eq_lat = bare_De_lattice()
print(f"  D_e_lattice = {De_lat_universal:.6f}  at R_eq = {R_eq_lat}")
conv_per_E_harm = 1.0 / (2 * S_PT)
print(f"  D_e (eV) = D_e_lattice * E_harm / (2*s_PT)")
print(f"           = {De_lat_universal:.4f} * E_harm * {conv_per_E_harm:.4f}")
print(f"           = E_harm * {De_lat_universal * conv_per_E_harm:.4f}")
print()
print(f"For pure sigma bonds, V8 base = pi/d^2 * E_harm = {PI/d**2:.4f} * E_harm")
print(f"Bare 1D Hessian      = {De_lat_universal * conv_per_E_harm:.4f} * E_harm")
print(f"Ratio bare/V8base    = {De_lat_universal*conv_per_E_harm/(PI/d**2):.4f}  "
      f"(should be ~1)")
print()

# ============================================================
# Per-bond table
# ============================================================
print(f"{'Name':>8} {'bo':>3} {'rad':>3} "
      f"{'E_harm':>7} {'D_bare':>8} {'D_bareBO':>9} "
      f"{'D_V8':>7} {'D_obs':>7} "
      f"{'bare_err':>9} {'bareBO_err':>11} {'V8_err':>8}")
print("-" * 110)

rows = []
for sa, sb, bo, D_obs, name, rad in TEST_BONDS:
    Z_a, IE_a, n_a, p_a, lp_a, m_a, _ = ATOMS[sa]
    Z_b, IE_b, n_b, p_b, lp_b, m_b, _ = ATOMS[sb]
    E_harm = 2*IE_a*IE_b/(IE_a+IE_b)

    D_bare = De_lat_universal * conv_per_E_harm * E_harm
    D_bare_BO = D_bare * (1 + (bo-1)*W_PI)

    v8r = bond_v8(sa, sb, bo, rad)
    D_V8 = v8r['D_0']

    bare_err = (D_bare - D_obs)/D_obs * 100
    bareBO_err = (D_bare_BO - D_obs)/D_obs * 100
    V8_err = (D_V8 - D_obs)/D_obs * 100

    rad_s = 'rad' if rad else ''
    print(f"{name:>8} {bo:>3} {rad_s:>3} "
          f"{E_harm:>7.2f} {D_bare:>8.3f} {D_bare_BO:>9.3f} "
          f"{D_V8:>7.3f} {D_obs:>7.3f} "
          f"{bare_err:>+8.1f}% {bareBO_err:>+10.1f}% {V8_err:>+7.1f}%")
    rows.append({'name': name, 'bo': bo, 'rad': rad,
                 'lp': lp_a + lp_b, 'het': sa != sb,
                 'D_bare': D_bare, 'D_bare_BO': D_bare_BO,
                 'D_V8': D_V8, 'D_obs': D_obs,
                 'bare_err': bare_err, 'bareBO_err': bareBO_err,
                 'V8_err': V8_err})

print()

# ============================================================
# Summary by category
# ============================================================
def stats(label, filtered, key='bare_err'):
    if not filtered:
        return
    errs = [r[key] for r in filtered]
    abs_errs = [abs(e) for e in errs]
    print(f"  {label:<40} n={len(filtered):>2}  "
          f"mean|err|={np.mean(abs_errs):>5.1f}%  "
          f"mean(err)={np.mean(errs):>+6.1f}%  "
          f"min={min(errs):>+6.1f}%  max={max(errs):>+6.1f}%")

print("=" * 78)
print("BARE 1D ERRORS BY CATEGORY (positive = over-predict, neg = under)")
print("=" * 78)
no_Li = [r for r in rows if r['name'] not in ('Li2', 'LiH')]
print("Bare (sigma-only, no bond-order, no corrections):")
stats("ALL (incl Li)", rows, 'bare_err')
stats("ALL (excl Li)", no_Li, 'bare_err')
stats("single-bond closed-shell, no LP",
      [r for r in no_Li if r['bo']==1 and not r['rad'] and r['lp']==0], 'bare_err')
stats("single-bond closed-shell, has LP",
      [r for r in no_Li if r['bo']==1 and not r['rad'] and r['lp']>0], 'bare_err')
stats("single-bond radical",
      [r for r in no_Li if r['bo']==1 and r['rad']], 'bare_err')
stats("bond order 2",
      [r for r in no_Li if r['bo']==2], 'bare_err')
stats("bond order 3",
      [r for r in no_Li if r['bo']==3], 'bare_err')

print()
print("Bare + naive bond-order scaling (1 + (bo-1)*W_PI):")
stats("ALL (excl Li)", no_Li, 'bareBO_err')
stats("single-bond closed-shell, no LP",
      [r for r in no_Li if r['bo']==1 and not r['rad'] and r['lp']==0], 'bareBO_err')
stats("single-bond closed-shell, has LP",
      [r for r in no_Li if r['bo']==1 and not r['rad'] and r['lp']>0], 'bareBO_err')
stats("bond order 2",
      [r for r in no_Li if r['bo']==2], 'bareBO_err')
stats("bond order 3",
      [r for r in no_Li if r['bo']==3], 'bareBO_err')

print()
print("V8 full (8 corrections):")
stats("ALL (excl Li)", no_Li, 'V8_err')

print()
print("=" * 78)
print("DIAGNOSTIC INTERPRETATION")
print("=" * 78)
print("""
The question: is the bare 1D Hessian's H2 match (0.6%) universal?

Look at 'single-bond closed-shell, no LP' bare error:
  - If small and unbiased -> bare 1D captures all single-sigma physics
    correctly. V8's job is only multi-bond + LP + ionic. Collective solver
    inherits the bare baseline as already-good for sigma bonds.
  - If systematically biased -> H2 is special, and even simple sigma bonds
    show the collective gap. v2 must do work even for trivial cases.

Look at 'bond order 2' and 'bond order 3' bare error vs bareBO error:
  - If bareBO closes most of the gap -> multi-bond is well-modeled by
    naive (1 + (bo-1)*W_PI) scaling. V8's extra multi-bond corrections
    are residual not load-bearing.
  - If bareBO still wildly off -> multi-bond requires real mode-mode
    coupling, and v2's multi-mode SCF is exactly the right addition.
""")
