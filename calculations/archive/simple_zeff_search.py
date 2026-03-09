"""
SIMPLE Z_eff SEARCH
====================
Strip it down to basics. Z_eff must follow from Z, n, l
through some simple wave relationship. Electrons are energy
in a standing wave — not particles being screened.

What is the SIMPLEST formula?
"""
import numpy as np

pi = np.pi
E_H = 13.6057
dd = 3
phi = (1 + np.sqrt(5)) / 2  # golden ratio

# Data: Z_eff from Clementi-Raimondi
data = [
    # name, Z, n, l, Z_eff, N_inner (electrons below valence shell)
    ('H',   1,  1, 0, 1.0000, 0),
    ('Li',  3,  2, 0, 1.2792, 2),
    ('B',   5,  2, 1, 2.4214, 4),
    ('C',   6,  2, 1, 3.1358, 5),
    ('N',   7,  2, 1, 3.8340, 6),
    ('O',   8,  2, 1, 4.4532, 7),
    ('F',   9,  2, 1, 5.0998, 8),
    ('Na', 11,  3, 0, 2.5074, 10),
    ('Cl', 17,  3, 1, 6.1161, 16),
]

# N_core = electrons in COMPLETE shells below (noble gas core)
# For Li-F: N_core = 2 (He core)
# For Na-Cl: N_core = 10 (Ne core)
cores = {'H': 0, 'Li': 2, 'B': 2, 'C': 2, 'N': 2, 'O': 2, 'F': 2,
         'Na': 10, 'Cl': 10}

# N_same = other electrons in same Slater group [ns, np]
sames = {'H': 0, 'Li': 0, 'B': 2, 'C': 3, 'N': 4, 'O': 5, 'F': 6,
         'Na': 0, 'Cl': 6}

# =============================================================================
# JUST LOOK AT THE NUMBERS
# =============================================================================
print("=" * 80)
print("  RAW PATTERNS")
print("=" * 80)

print(f"\n{'name':>4} {'Z':>3} {'n':>2} {'l':>2} {'Nin':>4} {'Nc':>3} {'Ns':>3} "
      f"{'Ze':>7} {'Z-Nin':>6} {'Z-Nc':>5} {'Ze/(Z-Nc)':>9} {'Ze*n/Z':>8}")
print("-" * 75)

for name, Z, n, l, Ze, Nin in data:
    Nc = cores[name]
    Ns = sames[name]
    Znet = Z - Nc
    ratio = Ze / Znet if Znet > 0 else 0
    print(f"{name:>4} {Z:3d} {n:2d} {l:2d} {Nin:4d} {Nc:3d} {Ns:3d} "
          f"{Ze:7.4f} {Z-Nin:6d} {Znet:5d} {ratio:9.4f} {Ze*n/Z:8.4f}")

# Key observation: Ze/(Z - N_core) is the fraction of valence charge retained
# B-F: 0.807, 0.784, 0.767, 0.742, 0.729 — decreasing
# These are close to 1 - Ns/(something)

print("\n  2p series: Ze/(Z-2) and 1 - Ns*s_same:")
for name, Z, n, l, Ze, Nin in data:
    if l == 1 and n == 2:
        Nc = cores[name]; Ns = sames[name]
        ratio = Ze / (Z - Nc)
        s_implied = (1 - ratio) * (Z - Nc) / Ns if Ns > 0 else 0
        print(f"    {name}: Ze/(Z-2)={ratio:.4f}, implied s_same = {s_implied:.4f}")


# =============================================================================
# SYSTEMATIC FORMULA SEARCH
# =============================================================================
print()
print("=" * 80)
print("  FORMULA SEARCH: Z_eff = f(Z, n, l, N_core, N_same)")
print("=" * 80)

def test_formula(formula_func, label):
    errs = []
    for name, Z, n, l, Ze, Nin in data:
        Nc = cores[name]; Ns = sames[name]
        Ze_pred = formula_func(Z, n, l, Nc, Ns)
        errs.append((Ze_pred - Ze)**2)
    rms = np.sqrt(np.mean(errs))
    return rms

# Test many formulas
formulas = {}

# Slater-style: Z - Nc*s_core - Ns*s_same
for sc_name, sc in [('6/7', 6/7), ('5/6', 5/6), ('f_pi', dd**2/(dd**2+1)),
                     ('(2d-1)/2d', (2*dd-1)/(2*dd))]:
    for ss_name, ss in [('1/3', 1/dd), ('1/pi', 1/pi), ('2/7', 2/7),
                         ('(d-1)/2d', (dd-1)/(2*dd)), ('1/(d+1)', 1/(dd+1))]:
        label = f"Z - {sc_name}*Nc - {ss_name}*Ns"
        f = lambda Z, n, l, Nc, Ns, _sc=sc, _ss=ss: Z - _sc*Nc - _ss*Ns
        formulas[label] = f

# With n-dependence: Z_eff = Z - Nc*(1-a/n^2) - Ns*b
for a_name, a in [('2/7', 2/7), ('4/7', 4/7), ('1/d', 1/dd), ('2/d', 2/dd)]:
    for b_name, b in [('1/3', 1/dd), ('2/7', 2/7), ('1/pi', 1/pi)]:
        label = f"Z - Nc*(1-{a_name}/n^2) - Ns*{b_name}"
        f = lambda Z, n, l, Nc, Ns, _a=a, _b=b: Z - Nc*(1-_a/n**2) - Ns*_b
        formulas[label] = f

# Z_eff = Z - Z*f(n) + correction
# Z_eff = (Z-Nc)*(1 - Ns/(2l+1+x))
for x_name, x in [('d', dd), ('2d', 2*dd), ('d+1', dd+1), ('2d+1', 2*dd+1),
                   ('d^2', dd**2), ('d*n', 0)]:
    label = f"(Z-Nc)*(1-Ns/({x_name}+Ns))" if x != 0 else f"(Z-Nc)*(1-Ns/(d*n+Ns))"
    if x != 0:
        f = lambda Z, n, l, Nc, Ns, _x=x: (Z-Nc)*(1-Ns/(_x+Ns))
    else:
        f = lambda Z, n, l, Nc, Ns: (Z-Nc)*(1-Ns/(dd*n+Ns))
    formulas[label] = f

# Power law: Z_eff = Z^a / n^b
for a in np.arange(0.5, 1.5, 0.1):
    for b in np.arange(0.0, 2.0, 0.1):
        label = f"Z^{a:.1f}/n^{b:.1f}"
        f = lambda Z, n, l, Nc, Ns, _a=a, _b=b: Z**_a / n**_b
        formulas[label] = f

# With l-dependence
for a_name, a in [('1/d', 1/dd), ('2/(2d+1)', 2/(2*dd+1)), ('1/(d+1)', 1/(dd+1))]:
    label = f"Z - (6/7+{a_name}*l)*Nc - (1/3)*Ns"
    f = lambda Z, n, l, Nc, Ns, _a=a: Z - (6/7+_a*l)*Nc - (1/dd)*Ns
    formulas[label] = f

# My harmonic model for comparison
label = "harmonic(g_s=2/3,g_d=4/7)"
def harmonic(Z, n, l, Nc, Ns):
    configs = {
        0: [],  # H
        2: [(1, 0, 2)],  # Li
        4: [(1, 0, 2), (2, 0, 2)],  # B (Nc=2, Ns=2, total inner=4)
        5: [(1, 0, 2), (2, 0, 2), (2, 1, 1)],
        6: [(1, 0, 2), (2, 0, 2), (2, 1, 2)],
        7: [(1, 0, 2), (2, 0, 2), (2, 1, 3)],
        8: [(1, 0, 2), (2, 0, 2), (2, 1, 4)],
        10: [(1, 0, 2), (2, 0, 2), (2, 1, 6)],
        16: [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 4)],
    }
    Nin = Nc + Ns
    if Nin not in configs:
        return Z  # fallback
    S = 0
    for (ni, li, cnt) in configs[Nin]:
        if ni == n and li == l:
            g = 2/dd
        else:
            g = 4/(2*dd+1)
        S += cnt * (1 - g*(ni/n)**2)
    return Z - S
formulas[label] = harmonic

# Sort by RMS and show top 20
results = []
for label, f in formulas.items():
    try:
        rms = test_formula(f, label)
        results.append((rms, label))
    except:
        pass

results.sort()
print(f"\n{'RMS':>8}  Formula")
print("-" * 70)
for rms, label in results[:25]:
    print(f"{rms:8.4f}  {label}")


# =============================================================================
# DEEP DIVE: Best formulas
# =============================================================================
print()
print("=" * 80)
print("  DETAIL: Top formulas")
print("=" * 80)

for rms, label in results[:5]:
    f = formulas[label]
    print(f"\n  {label} (RMS={rms:.4f}):")
    print(f"  {'name':>4} {'Z_pred':>7} {'Z_real':>7} {'err':>7}")
    for name, Z, n, l, Ze, Nin in data:
        Nc = cores[name]; Ns = sames[name]
        Ze_pred = f(Z, n, l, Nc, Ns)
        print(f"  {name:>4} {Ze_pred:7.4f} {Ze:7.4f} {Ze_pred-Ze:+7.4f}")


# =============================================================================
# REALLY SIMPLE: scan Z_eff = Z - a*Nc - b*Ns (2 params from d)
# =============================================================================
print()
print("=" * 80)
print("  EXHAUSTIVE 2-PARAM SCAN: Z_eff = Z - a*Nc - b*Ns")
print("=" * 80)

# Where a, b are GWT constants or simple fractions of d
candidates = {}
for p in range(1, 30):
    for q in range(1, 30):
        candidates[f"{p}/{q}"] = p/q

# Also add GWT constants
gwt = {'pi/d': pi/dd, 'd/pi': dd/pi, 'sqrt(d)': np.sqrt(dd),
       '(d+1)/d': (dd+1)/dd, 'sqrt(2)': np.sqrt(2),
       'e/d': np.e/dd, 'pi/(d+1)': pi/(dd+1)}
candidates.update(gwt)

best = (999, '', '', 0, 0)
for a_name, a_val in candidates.items():
    if a_val < 0.5 or a_val > 1.2:
        continue
    for b_name, b_val in candidates.items():
        if b_val < 0.1 or b_val > 0.6:
            continue
        errs = []
        for name, Z, n, l, Ze, Nin in data:
            Nc = cores[name]; Ns = sames[name]
            Ze_pred = Z - a_val*Nc - b_val*Ns
            errs.append((Ze_pred - Ze)**2)
        rms = np.sqrt(np.mean(errs))
        if rms < best[0]:
            best = (rms, a_name, b_name, a_val, b_val)

print(f"\n  Best: a={best[1]}={best[3]:.6f}, b={best[2]}={best[4]:.6f}, RMS={best[0]:.4f}")
a_best = best[3]; b_best = best[4]

print(f"\n  {'name':>4} {'Z_pred':>7} {'Z_real':>7} {'err':>7}")
print("  " + "-" * 30)
for name, Z, n, l, Ze, Nin in data:
    Nc = cores[name]; Ns = sames[name]
    Ze_pred = Z - a_best*Nc - b_best*Ns
    print(f"  {name:>4} {Ze_pred:7.4f} {Ze:7.4f} {Ze_pred-Ze:+7.4f}")


# =============================================================================
# 3-PARAM: Z_eff = Z - a*N_deep - b*N_adj - c*Ns (Slater groups)
# =============================================================================
print()
print("=" * 80)
print("  3-PARAM SLATER: Z_eff = Z - a*N_deep - b*N_adj - c*Ns")
print("=" * 80)

# N_deep = electrons 2+ groups below (screen most)
# N_adj = electrons 1 group below (screen somewhat)
# Ns = same group
# For 2p atoms: N_deep = 0 (no groups 2+ below), N_adj = 2 (1s^2), Ns = 2s^2 + other 2p

# Actually, Slater groups: [1s] [2s,2p] [3s,3p]
# For Li (2s): N_deep=0, N_adj=2(1s), Ns=0
# For B (2p): N_deep=0, N_adj=2(1s), Ns=2(2s)
# For Na (3s): N_deep=2(1s), N_adj=8(2s2p), Ns=0
# For Cl (3p): N_deep=2(1s), N_adj=8(2s2p), Ns=6(3s+3p)

slater = {
    'H':  (0, 0, 0),
    'Li': (0, 2, 0),
    'B':  (0, 2, 2),
    'C':  (0, 2, 3),
    'N':  (0, 2, 4),
    'O':  (0, 2, 5),
    'F':  (0, 2, 6),
    'Na': (2, 8, 0),
    'Cl': (2, 8, 6),
}

best3 = (999, 0, 0, 0)
for a in np.arange(0.8, 1.05, 0.005):
    for b in np.arange(0.7, 1.0, 0.005):
        for c in np.arange(0.2, 0.5, 0.005):
            errs = []
            for name, Z, n, l, Ze, Nin in data:
                Nd, Na, Ns = slater[name]
                Ze_pred = Z - a*Nd - b*Na - c*Ns
                errs.append((Ze_pred - Ze)**2)
            rms = np.sqrt(np.mean(errs))
            if rms < best3[0]:
                best3 = (rms, a, b, c)

a3, b3, c3 = best3[1], best3[2], best3[3]
print(f"\n  Best: s_deep={a3:.4f}, s_adj={b3:.4f}, s_same={c3:.4f}, RMS={best3[0]:.4f}")

# Compare to GWT
print(f"\n  s_deep = {a3:.4f}")
print(f"    1 = {1.0:.4f}")
print(f"    (2d+1)/(2d) = {(2*dd+1)/(2*dd):.4f}")

print(f"\n  s_adj = {b3:.4f}")
print(f"    6/7 = 2d/(2d+1) = {2*dd/(2*dd+1):.4f}")
print(f"    5/6 = (2d-1)/(2d) = {(2*dd-1)/(2*dd):.4f}")
print(f"    pi/d - 1/(2d) = {pi/dd - 1/(2*dd):.4f}")

print(f"\n  s_same = {c3:.4f}")
print(f"    1/3 = 1/d = {1/dd:.4f}")
print(f"    1/pi = {1/pi:.4f}")
print(f"    2/7 = 2/(2d+1) = {2/(2*dd+1):.4f}")
print(f"    5/14 = {5/14:.4f}")

# Predictions
print(f"\n  {'name':>4} {'Nd':>3} {'Na':>3} {'Ns':>3} {'Z_pred':>7} {'Z_real':>7} {'err':>7}")
print("  " + "-" * 40)
for name, Z, n, l, Ze, Nin in data:
    Nd, Na_val, Ns = slater[name]
    Ze_pred = Z - a3*Nd - b3*Na_val - c3*Ns
    print(f"  {name:>4} {Nd:3d} {Na_val:3d} {Ns:3d} {Ze_pred:7.4f} {Ze:7.4f} {Ze_pred-Ze:+7.4f}")
