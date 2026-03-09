"""Deep pattern hunt for Z_eff - looking for the simplest formula."""
import numpy as np
from math import pi, sqrt

atoms = {
    'H':  {'Z': 1,  'n': 1, 'l': 0, 'Ze': 1.0000},
    'Li': {'Z': 3,  'n': 2, 'l': 0, 'Ze': 1.2792},
    'B':  {'Z': 5,  'n': 2, 'l': 1, 'Ze': 2.4214},
    'C':  {'Z': 6,  'n': 2, 'l': 1, 'Ze': 3.1358},
    'N':  {'Z': 7,  'n': 2, 'l': 1, 'Ze': 3.8340},
    'O':  {'Z': 8,  'n': 2, 'l': 1, 'Ze': 4.4532},
    'F':  {'Z': 9,  'n': 2, 'l': 1, 'Ze': 5.0998},
    'Na': {'Z': 11, 'n': 3, 'l': 0, 'Ze': 2.5074},
    'Cl': {'Z': 17, 'n': 3, 'l': 1, 'Ze': 6.1161},
}

d = 3
g_same = 2/3
g_diff = 4/7

configs = {
    'H':  [(1, 0, 0)],
    'Li': [(1, 0, 2), (2, 0, 0)],
    'B':  [(1, 0, 2), (2, 0, 2), (2, 1, 0)],
    'C':  [(1, 0, 2), (2, 0, 2), (2, 1, 1)],
    'N':  [(1, 0, 2), (2, 0, 2), (2, 1, 2)],
    'O':  [(1, 0, 2), (2, 0, 2), (2, 1, 3)],
    'F':  [(1, 0, 2), (2, 0, 2), (2, 1, 4)],
    'Na': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 0)],
    'Cl': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 4)],
}

# Reference: harmonic model predictions
print("="*70)
print("  HARMONIC MODEL (reference)")
print("="*70)
harmonic_pred = {}
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs[name]:
        if ni == nv and li == lv:
            s_per = 1 - g_same * (ni/nv)**2
        else:
            s_per = 1 - g_diff * (ni/nv)**2
        total_s += s_per * Ne
    pred = Z - total_s
    harmonic_pred[name] = pred
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# Key: Na error = +0.76, Cl error = +0.54. Both overestimate.
# Period 2 errors are small. What's different about period 3?

print("\n" + "="*70)
print("  SECOND-ORDER COUPLING: s = 1 - g*(ni/nv)^2 + g2*(ni/nv)^4")
print("="*70)
print("  Idea: higher-order term increases screening for deeper inner shells")

# Single g2 scan
best_rms = 999
best_g2 = 0
for g2_1000 in range(-500, 500):
    g2 = g2_1000 / 1000
    rms = 0
    for name, at in atoms.items():
        Z, nv, lv = at['Z'], at['n'], at['l']
        total_s = 0
        for (ni, li, Ne) in configs[name]:
            r2 = (ni/nv)**2
            r4 = (ni/nv)**4
            g = g_same if (ni == nv and li == lv) else g_diff
            s_per = 1 - g * r2 + g2 * r4
            total_s += s_per * Ne
        pred = Z - total_s
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_g2 = g2

print(f"  Best g2 = {best_g2:.3f}, RMS = {best_rms:.4f}")
print(f"    Compare: 2/7 = {2/7:.4f}, 1/3 = {1/3:.4f}, 2/9 = {2/9:.4f}")
print(f"    g_same*g_diff = {g_same*g_diff:.4f}")

for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs[name]:
        r2 = (ni/nv)**2
        r4 = (ni/nv)**4
        g = g_same if (ni == nv and li == lv) else g_diff
        s_per = 1 - g * r2 + best_g2 * r4
        total_s += s_per * Ne
    pred = Z - total_s
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# Separate g2 for same/diff
print("\n--- Separate g2_same, g2_diff ---")
best_rms = 999
best_params = (0, 0)
for g2s_100 in range(-100, 100):
    g2s = g2s_100 / 100
    for g2d_100 in range(-100, 100):
        g2d = g2d_100 / 100
        rms = 0
        for name, at in atoms.items():
            Z, nv, lv = at['Z'], at['n'], at['l']
            total_s = 0
            for (ni, li, Ne) in configs[name]:
                r2 = (ni/nv)**2
                r4 = (ni/nv)**4
                if ni == nv and li == lv:
                    s_per = 1 - g_same * r2 + g2s * r4
                else:
                    s_per = 1 - g_diff * r2 + g2d * r4
                total_s += s_per * Ne
            pred = Z - total_s
            rms += (pred - at['Ze'])**2
        rms = sqrt(rms / len(atoms))
        if rms < best_rms:
            best_rms = rms
            best_params = (g2s, g2d)

g2s, g2d = best_params
print(f"  Best: g2_same={g2s:.2f}, g2_diff={g2d:.2f}, RMS={best_rms:.4f}")
print(f"  Compare: g_same^2={g_same**2:.4f}, g_diff^2={g_diff**2:.4f}")
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs[name]:
        r2 = (ni/nv)**2
        r4 = (ni/nv)**4
        if ni == nv and li == lv:
            s_per = 1 - g_same * r2 + g2s * r4
        else:
            s_per = 1 - g_diff * r2 + g2d * r4
        total_s += s_per * Ne
    pred = Z - total_s
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# RADICAL SIMPLICITY SEARCH
print("\n" + "="*70)
print("  RADICAL SIMPLICITY")
print("="*70)

N_inner = {'H': 0, 'Li': 2, 'B': 4, 'C': 5, 'N': 6, 'O': 7, 'F': 8, 'Na': 10, 'Cl': 16}

# Ze = Z - N*(1 - c/n)
print("\n--- Ze = Z - N_inner*(1 - c/n) ---")
best_rms = 999
for c1000 in range(1, 2000):
    c = c1000 / 1000
    rms = 0
    for name, at in atoms.items():
        pred = at['Z'] - N_inner[name] * (1 - c / at['n'])
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_c = c
print(f"  Best c = {best_c:.3f}, RMS = {best_rms:.4f}")
for name, at in atoms.items():
    pred = at['Z'] - N_inner[name] * (1 - best_c / at['n'])
    print(f"    {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# Ze = Z - N*(1 - c/n^2)
print("\n--- Ze = Z - N_inner*(1 - c/n^2) ---")
best_rms = 999
for c1000 in range(1, 3000):
    c = c1000 / 1000
    rms = 0
    for name, at in atoms.items():
        pred = at['Z'] - N_inner[name] * (1 - c / at['n']**2)
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_c = c
print(f"  Best c = {best_c:.3f}, RMS = {best_rms:.4f}")

# Ze = Z^(n/(n+c))
print("\n--- Ze = Z^(n/(n+c)) ---")
best_rms = 999
for c1000 in range(1, 3000):
    c = c1000 / 1000
    rms = 0
    for name, at in atoms.items():
        pred = at['Z'] ** (at['n'] / (at['n'] + c))
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_c = c
print(f"  Best c = {best_c:.3f}, RMS = {best_rms:.4f}")
for name, at in atoms.items():
    pred = at['Z'] ** (at['n'] / (at['n'] + best_c))
    print(f"    {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# WHAT IF: Ze = (Z - N_core) * (1 - s_same * N_same / (Z - N_core))
# with s_same that depends on n?
# Try: s_same = a / n
N_core = {'H': 0, 'Li': 2, 'B': 2, 'C': 2, 'N': 2, 'O': 2, 'F': 2, 'Na': 10, 'Cl': 10}
N_same = {'H': 0, 'Li': 0, 'B': 2, 'C': 3, 'N': 4, 'O': 5, 'F': 6, 'Na': 0, 'Cl': 6}

print("\n--- 3-param: Ze = Z - a*Nc - b(n)*Ns where b = c/n + d ---")
best_rms = 999
for a100 in range(50, 120):
    a = a100 / 100
    for c100 in range(0, 100):
        c = c100 / 100
        for d100 in range(0, 80):
            dd = d100 / 100
            rms = 0
            for name, at in atoms.items():
                b = c / at['n'] + dd
                pred = at['Z'] - a * N_core[name] - b * N_same[name]
                rms += (pred - at['Ze'])**2
            rms = sqrt(rms / len(atoms))
            if rms < best_rms:
                best_rms = rms
                best_params = (a, c, dd)

a, c, dd = best_params
print(f"  Best: s_core={a:.2f}, c={c:.2f}, d={dd:.2f}")
print(f"  s_same(n=2) = {c/2+dd:.4f}, s_same(n=3) = {c/3+dd:.4f}")
print(f"  RMS = {best_rms:.4f}")
for name, at in atoms.items():
    b = c / at['n'] + dd
    pred = at['Z'] - a * N_core[name] - b * N_same[name]
    print(f"    {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# KEY QUESTION: What does the harmonic model get WRONG about period 3?
# Let me decompose the Na error
print("\n" + "="*70)
print("  DECOMPOSING HARMONIC MODEL ERRORS")
print("="*70)

print("\n  Na (real Ze = 2.507):")
nv = 3
pred_s = 0
for (ni, li, Ne) in configs['Na']:
    if Ne == 0: continue
    r2 = (ni/nv)**2
    g = g_same if (ni == nv and li == 0) else g_diff
    s_per = 1 - g * r2
    contrib = s_per * Ne
    pred_s += contrib
    print(f"    ({ni},{li}) x{Ne}: r2={r2:.4f}, g={g:.4f}, s={s_per:.4f}, total_s={contrib:.4f}")
print(f"    Sum screening = {pred_s:.4f}")
print(f"    Pred Ze = 11 - {pred_s:.4f} = {11-pred_s:.4f}")
print(f"    Need extra screening of {11-pred_s - 2.5074:.4f}")

# What fraction of the n=2 shell screening is "missing"?
# n=2 shell: 2s^2 + 2p^6 = 8 electrons, each at (2/3)^2 = 4/9 ratio
# With g_diff (2p->3s): s = 1 - (4/7)(4/9) = 1 - 16/63 = 47/63 = 0.746
# With g_same (2s->3s): s = 1 - (2/3)(4/9) = 1 - 8/27 = 19/27 = 0.704
missing = 11 - pred_s - 2.5074
print(f"\n    Missing screening = {missing:.4f}")
print(f"    Per n=2 electron (8): {missing/8:.4f}")
print(f"    Per all inner (10): {missing/10:.4f}")

# WHAT IF: the coupling g should be LARGER for inner shells that are
# more than 1 quantum number away?
# Currently: g is the same whether n_i=1 or n_i=2 for outer n=3
# What if g_far = g * (nv - ni) / (nv - 1)?
print("\n\n--- DISTANCE-DEPENDENT g ---")
print("  g_eff = g * (nv - ni + 1) / nv  (farther = stronger coupling)")
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        g_base = g_same if (ni == nv and li == lv) else g_diff
        # Distance factor: how far is this shell?
        dist_factor = (nv - ni + 1) / nv if ni < nv else 1
        g_eff = g_base * dist_factor
        s_per = 1 - g_eff * r2
        total_s += s_per * Ne
    pred = Z - total_s
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# Try: g_eff = g * n_v / (n_v - (n_v - n_i) * c)
print("\n--- SCAN: g_eff = g * (1 + c*(nv-ni)/nv) ---")
best_rms = 999
for c1000 in range(-2000, 2000):
    c = c1000 / 1000
    rms = 0
    for name, at in atoms.items():
        Z, nv, lv = at['Z'], at['n'], at['l']
        total_s = 0
        for (ni, li, Ne) in configs[name]:
            if Ne == 0: continue
            r2 = (ni/nv)**2
            g_base = g_same if (ni == nv and li == lv) else g_diff
            dn = nv - ni
            g_eff = g_base * (1 + c * dn / nv)
            s_per = 1 - g_eff * r2
            total_s += s_per * Ne
        pred = Z - total_s
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_c = c

print(f"  Best c = {best_c:.3f}, RMS = {best_rms:.4f}")
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        g_base = g_same if (ni == nv and li == lv) else g_diff
        dn = nv - ni
        g_eff = g_base * (1 + best_c * dn / nv)
        s_per = 1 - g_eff * r2
        total_s += s_per * Ne
    pred = Z - total_s
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# ALTERNATIVE: What if the coupling uses (ni/nv)^p where p varies?
# We found p=2 works for period 2. What if p should be DIFFERENT?
print("\n--- SCAN: s = 1 - g * (ni/nv)^p for different p ---")
for p10 in [15, 18, 20, 22, 25, 30]:
    p = p10 / 10
    rms = 0
    for name, at in atoms.items():
        Z, nv, lv = at['Z'], at['n'], at['l']
        total_s = 0
        for (ni, li, Ne) in configs[name]:
            rp = (ni/nv)**p
            g = g_same if (ni == nv and li == lv) else g_diff
            s_per = 1 - g * rp
            total_s += s_per * Ne
        pred = Z - total_s
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    print(f"  p={p:.1f}: RMS={rms:.4f}")

# Fine scan of p
best_rms = 999
for p100 in range(100, 400):
    p = p100 / 100
    rms = 0
    for name, at in atoms.items():
        Z, nv, lv = at['Z'], at['n'], at['l']
        total_s = 0
        for (ni, li, Ne) in configs[name]:
            rp = (ni/nv)**p
            g = g_same if (ni == nv and li == lv) else g_diff
            s_per = 1 - g * rp
            total_s += s_per * Ne
        pred = Z - total_s
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_p = p

print(f"\n  Best p = {best_p:.2f}, RMS = {best_rms:.4f}")
print(f"  Compare: p=2 gives harmonic model")
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs[name]:
        rp = (ni/nv)**best_p
        g = g_same if (ni == nv and li == lv) else g_diff
        s_per = 1 - g * rp
        total_s += s_per * Ne
    pred = Z - total_s
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# What if g itself depends on the ratio? g_eff = g * (ni/nv)?
# Then s = 1 - g * (ni/nv)^3 effectively
print("\n--- s = 1 - g * (ni/nv)^p, best p with g_same, g_diff optimized ---")
for p10 in [15, 18, 20, 22, 25, 28, 30]:
    p = p10 / 10
    best_rms = 999
    for gs100 in range(10, 100):
        gs = gs100 / 100
        for gd100 in range(10, 100):
            gd = gd100 / 100
            rms = 0
            for name, at in atoms.items():
                Z, nv, lv = at['Z'], at['n'], at['l']
                total_s = 0
                for (ni, li, Ne) in configs[name]:
                    rp = (ni/nv)**p
                    g = gs if (ni == nv and li == lv) else gd
                    s_per = 1 - g * rp
                    total_s += s_per * Ne
                pred = Z - total_s
                rms += (pred - at['Ze'])**2
            rms = sqrt(rms / len(atoms))
            if rms < best_rms:
                best_rms = rms
                best_gs, best_gd = gs, gd
    print(f"  p={p:.1f}: gs={best_gs:.2f}, gd={best_gd:.2f}, RMS={best_rms:.4f}")

# COMPLETELY NEW IDEA: What if the wave coupling depends on the
# OVERLAP INTEGRAL between harmonics?
# For spherical harmonics, the overlap between n1,l1 and n2,l2
# goes as something related to their spatial overlap
# In a hydrogen atom, <r> = n^2 * a0 / Z
# Spatial overlap ~ exp(-|r1-r2|) or some power law

# The key physical insight: inner wave modes DON'T just partially
# screen. They shift the effective potential. The AMOUNT of shift
# depends on how much the inner wave's energy density overlaps
# with the outer wave.

# For hydrogenic wavefunctions:
# Probability of finding r < r0: depends on n,l
# The inner electrons create a charge distribution, and the outer
# electron sees Z_eff(r) that varies with r.

# SIMPLEST VERSION: What if Z_eff = Z - sum_i [1 - overlap(i, outer)]
# where overlap = (n_i/n_v)^2 * angular_factor ?
# This IS the harmonic model! So the issue must be in the angular factor.

# What if g depends on whether the inner shell is COMPLETE?
print("\n\n--- COMPLETENESS-DEPENDENT g ---")
print("  Complete subshells screen differently than partial ones?")

# In Na, the 2p^6 is a complete subshell (filled)
# In F, the 2p^4 is NOT complete
# What if a filled subshell screens with a DIFFERENT g?
# Filled subshells: 1s^2 always; 2s^2 in B-F,Na,Cl; 2p^6 in Na,Cl; 3s^2 in Cl
# For Cl outer 3p: 3s^2 is filled, 3p^4 is partial

configs_filled = {
    'H':  [(1, 0, 0, False)],
    'Li': [(1, 0, 2, True), (2, 0, 0, False)],
    'B':  [(1, 0, 2, True), (2, 0, 2, True), (2, 1, 0, False)],
    'C':  [(1, 0, 2, True), (2, 0, 2, True), (2, 1, 1, False)],
    'N':  [(1, 0, 2, True), (2, 0, 2, True), (2, 1, 2, False)],
    'O':  [(1, 0, 2, True), (2, 0, 2, True), (2, 1, 3, False)],
    'F':  [(1, 0, 2, True), (2, 0, 2, True), (2, 1, 4, False)],
    'Na': [(1, 0, 2, True), (2, 0, 2, True), (2, 1, 6, True), (3, 0, 0, False)],
    'Cl': [(1, 0, 2, True), (2, 0, 2, True), (2, 1, 6, True), (3, 0, 2, True), (3, 1, 4, False)],
}
# Scan: filled shells get g_filled instead of g_diff
best_rms = 999
for gf100 in range(10, 100):
    gf = gf100 / 100
    rms = 0
    for name, at in atoms.items():
        Z, nv, lv = at['Z'], at['n'], at['l']
        total_s = 0
        for (ni, li, Ne, filled) in configs_filled[name]:
            r2 = (ni/nv)**2
            if ni == nv and li == lv:
                g = g_same
            elif filled and ni < nv:
                g = gf  # filled inner subshell
            else:
                g = g_diff
            s_per = 1 - g * r2
            total_s += s_per * Ne
        pred = Z - total_s
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_gf = gf

print(f"  Best g_filled = {best_gf:.2f}, RMS = {best_rms:.4f}")
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne, filled) in configs_filled[name]:
        r2 = (ni/nv)**2
        if ni == nv and li == lv:
            g = g_same
        elif filled and ni < nv:
            g = best_gf
        else:
            g = g_diff
        s_per = 1 - g * r2
        total_s += s_per * Ne
    pred = Z - total_s
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# What if there is a COLLECTIVE effect: filled shells screen MORE
# because they form a complete spherical shield?
# Extra screening = delta * N_filled_shells (number of complete SHELLS, not electrons)
print("\n--- Collective shell correction ---")
n_complete_shells = {'H': 0, 'Li': 0, 'B': 0, 'C': 0, 'N': 0, 'O': 0, 'F': 0, 'Na': 1, 'Cl': 1}
# Na has one complete shell below (n=2: 2s^2 2p^6 = 8)
# Cl has one complete shell below (n=2), but also n=1 is complete
# Actually: n=1 is always complete for Z>=2. n=2 is complete for Z>=10.
n_complete_inner = {'H': 0, 'Li': 1, 'B': 1, 'C': 1, 'N': 1, 'O': 1, 'F': 1, 'Na': 2, 'Cl': 2}
# Li-F: n=1 shell (1s^2) is complete = 1 complete shell
# Na,Cl: n=1 (1s^2) AND n=2 (2s^2 2p^6) = 2 complete shells

best_rms = 999
for delta100 in range(0, 200):
    delta = delta100 / 100
    rms = 0
    for name, at in atoms.items():
        pred = harmonic_pred[name] - delta * n_complete_inner[name]
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_delta = delta

print(f"  Best delta = {best_delta:.2f} per complete shell, RMS = {best_rms:.4f}")
print(f"    1/d = {1/d:.4f}, 2/7 = {2/7:.4f}, 1/pi = {1/pi:.4f}")
for name, at in atoms.items():
    pred = harmonic_pred[name] - best_delta * n_complete_inner[name]
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# What about total electrons in complete inner shells?
print("\n--- delta per electron in complete inner shells ---")
n_e_complete = {'H': 0, 'Li': 2, 'B': 2, 'C': 2, 'N': 2, 'O': 2, 'F': 2, 'Na': 10, 'Cl': 10}
best_rms = 999
for delta1000 in range(0, 200):
    delta = delta1000 / 1000
    rms = 0
    for name, at in atoms.items():
        pred = harmonic_pred[name] - delta * n_e_complete[name]
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_delta = delta

print(f"  Best delta = {best_delta:.3f} per electron, RMS = {best_rms:.4f}")
for name, at in atoms.items():
    pred = harmonic_pred[name] - best_delta * n_e_complete[name]
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# FINAL SCAN: What if the harmonic model is RIGHT but needs one
# simple additive correction?
# Ze_corrected = Ze_harmonic - c * f(atom properties)
print("\n" + "="*70)
print("  CORRECTIONS TO HARMONIC MODEL")
print("="*70)
print("  Harmonic errors:")
errs = {}
for name in atoms:
    errs[name] = harmonic_pred[name] - atoms[name]['Ze']
    print(f"    {name}: {errs[name]:+.4f}")

# What correlates with the error?
print("\n  Error vs various quantities:")
for name in atoms:
    at = atoms[name]
    Ni = N_inner[name]  # was at['Z'] - 1 for total inner, let me use actual
    n_val = at['n']
    print(f"    {name}: err={errs[name]:+.4f}  Z/n={at['Z']/n_val:.2f}  Ni/n^2={Ni/n_val**2:.3f}  Z*Ni/n^3={at['Z']*Ni/n_val**3:.3f}")
