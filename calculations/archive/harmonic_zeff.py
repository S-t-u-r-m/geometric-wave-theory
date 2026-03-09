"""
HARMONIC PATTERN IN Z_eff
===========================
The user's insight: Z_eff should follow a harmonic pattern because
electrons are wave modes. The screening by inner electrons should
depend on harmonic relationships between their wave patterns.

Key observation from data:
  2p series: Z_eff increments are ~0.67 = 2/d per Z step
  This means each additional proton adds 1 charge, but each
  additional same-shell electron screens by ~1/d.

Let's look for the harmonic structure.
"""
import numpy as np

pi = np.pi
E_H = 13.6057
dd = 3
f_pi = dd**2 / (dd**2 + 1)
alpha_n = 1 - f_pi / dd
beta_n = (1 + f_pi) / 2

# Known data
atoms = {
    'H':  {'Z': 1,  'n': 1, 'l': 0, 'Z_eff': 1.0000},
    'Li': {'Z': 3,  'n': 2, 'l': 0, 'Z_eff': 1.2792},
    'B':  {'Z': 5,  'n': 2, 'l': 1, 'Z_eff': 2.4214},
    'C':  {'Z': 6,  'n': 2, 'l': 1, 'Z_eff': 3.1358},
    'N':  {'Z': 7,  'n': 2, 'l': 1, 'Z_eff': 3.8340},
    'O':  {'Z': 8,  'n': 2, 'l': 1, 'Z_eff': 4.4532},
    'F':  {'Z': 9,  'n': 2, 'l': 1, 'Z_eff': 5.0998},
    'Na': {'Z': 11, 'n': 3, 'l': 0, 'Z_eff': 2.5074},
    'Cl': {'Z': 17, 'n': 3, 'l': 1, 'Z_eff': 6.1161},
}

# =============================================================================
# OBSERVATION 1: The 2p increment pattern
# =============================================================================
print("=" * 80)
print("  2p SERIES: Increment pattern")
print("=" * 80)

series_2p = ['B', 'C', 'N', 'O', 'F']
print(f"\n{'Atom':>4} {'Z':>3} {'Z_eff':>7} {'dZ_eff':>7} {'1-dZ_eff':>8} {'screen':>7}")
print("-" * 45)

prev = None
for name in series_2p:
    Z = atoms[name]['Z']; Ze = atoms[name]['Z_eff']
    if prev is not None:
        dZe = Ze - prev
        screen = 1 - dZe  # how much the new same-shell electron screens
        print(f"{name:>4} {Z:3d} {Ze:7.4f} {dZe:7.4f} {screen:8.4f} {screen:7.4f}")
    else:
        print(f"{name:>4} {Z:3d} {Ze:7.4f}       -        -       -")
    prev = Ze

print(f"\n  Average increment: {np.mean([atoms[s]['Z_eff'] for s in series_2p[1:]]) - np.mean([atoms[s]['Z_eff'] for s in series_2p[:-1]]):.4f}")

# The screening per 2p electron varies: 0.286, 0.302, 0.381, 0.354
# Average: ~0.33 = 1/d !
screens = []
for i in range(1, len(series_2p)):
    dZe = atoms[series_2p[i]]['Z_eff'] - atoms[series_2p[i-1]]['Z_eff']
    screens.append(1 - dZe)
print(f"  Average same-shell screening: {np.mean(screens):.4f}")
print(f"  Compare 1/d = {1/dd:.4f}")
print(f"  Compare 1/pi = {1/pi:.4f}")

# =============================================================================
# OBSERVATION 2: Boron as the "fundamental" 2p atom
# =============================================================================
print()
print("=" * 80)
print("  BORON: The fundamental 2p (first electron in shell)")
print("=" * 80)

# B has Z=5: 1s^2, 2s^2, 2p^1
# Z_eff(B) = 5 - S(1s^2) - S(2s^2)
# = 5 - 2*s_1s - 2*s_2s = 2.4214
# So 2*s_1s + 2*s_2s = 2.5786, s_1s + s_2s = 1.2893

# Li has Z=3: 1s^2, 2s^1
# Z_eff(Li) = 3 - 2*s_1s_for_2s = 1.2792
# s_1s_for_2s = (3 - 1.2792) / 2 = 0.8604

s_1s_for_2s = (3 - atoms['Li']['Z_eff']) / 2
print(f"\n  Li: s_1s(for 2s) = (3 - {atoms['Li']['Z_eff']}) / 2 = {s_1s_for_2s:.4f}")
print(f"    compare: (2d-1)/(2d) = {(2*dd-1)/(2*dd):.4f}")
print(f"    compare: f_pi = {f_pi:.4f}")
print(f"    compare: pi/d - 1/(2d) = {pi/dd - 1/(2*dd):.4f}")

# For B, 1s screens 2p: could be different from 1s screening 2s
# s_1s_for_2p + s_2s_for_2p = 1.2893
s_sum = (5 - atoms['B']['Z_eff']) / 2  # total screening per pair
print(f"\n  B: (s_1s + s_2s for 2p) / 2 = {s_sum:.4f}")

# If s_1s_for_2p = s_1s_for_2s (same screening):
s_2s_for_2p = s_sum - s_1s_for_2s
print(f"  If s_1s same: s_2s(for 2p) = {s_2s_for_2p:.4f}")
print(f"    compare: 1/(2d) = {1/(2*dd):.4f}")
print(f"    compare: 1/(d+1) = {1/(dd+1):.4f}")

# =============================================================================
# HARMONIC SCREENING MODEL
# =============================================================================
print()
print("=" * 80)
print("  HARMONIC SCREENING MODEL")
print("=" * 80)
print()
print("  Hypothesis: screening depends on HARMONIC RATIO of wave modes")
print("  s(n_i, l_i -> n_v, l_v) = f(n_i/n_v, l_i, l_v)")

# For hydrogen wavefunctions, the screening integral depends on
# how much of the inner electron's density lies inside the outer electron's orbit
# This is related to the ratio of their principal quantum numbers

# Harmonic model: s = 1 - amplitude of harmonic coupling
# If inner mode n_i couples to outer mode n_v with strength ~1/n_v
# then screening = 1 - 1/n_v (for n_i < n_v)
# But this gives same screening for all inner electrons regardless of n_i

# Better: s depends on the RATIO n_i/n_v
# At resonance (n_i = n_v, same shell), coupling is maximum
# Off-resonance (n_i << n_v), inner electron fully screens

# Try: s = 1 - |sin(pi * n_i / n_v)| / something
# Or: s = 1 - (n_i/n_v)^p for some power p

# For same shell (n_i = n_v): s should be small (~1/d)
# For different shell (n_i < n_v): s should be large (~5/6)

# What if: s = 1 - f(l) where f depends on orbital type matching?

# Let's try a systematic approach: compute what screening each shell provides
print("Computing implied screening constants from data:")
print()

# For each atom, decompose Z_eff = Z - sum(s_i)
# Group electrons by (n, l)

electron_configs = {
    'H':  [(1, 0, 0)],  # (n, l, count of OTHER electrons in that subshell)
    'Li': [(1, 0, 2)],  # 1s^2 screens 2s
    'B':  [(1, 0, 2), (2, 0, 2)],  # 1s^2 and 2s^2 screen 2p^1
    'C':  [(1, 0, 2), (2, 0, 2), (2, 1, 1)],  # + 1 other 2p
    'N':  [(1, 0, 2), (2, 0, 2), (2, 1, 2)],
    'O':  [(1, 0, 2), (2, 0, 2), (2, 1, 3)],
    'F':  [(1, 0, 2), (2, 0, 2), (2, 1, 4)],
    'Na': [(1, 0, 2), (2, 0, 2), (2, 1, 6)],  # screens 3s
    'Cl': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 4)],  # screens 3p
}

# Total screening S = Z - Z_eff must equal sum of all screening contributions
# S = sum_over_subshells(count_i * s_i)

# We have these unknowns:
# s_1s_2s: 1s electron screening 2s (used in Li)
# s_1s_2p: 1s electron screening 2p (used in B-F)
# s_2s_2p: 2s electron screening 2p (used in B-F)
# s_2p_2p: 2p electron screening 2p (used in C-F)
# s_1s_3s: 1s screening 3s (used in Na)
# s_2s_3s: 2s screening 3s
# s_2p_3s: 2p screening 3s
# s_1s_3p: 1s screening 3p (used in Cl)
# etc.

# Let's simplify: assume screening depends only on the SHELL GAP
# s(n_i -> n_v) with same-shell having different l counted separately

# From Li: 2*s(1->2) = 3 - 1.2792 = 1.7208, s(1->2) = 0.8604
s_12 = (3 - atoms['Li']['Z_eff']) / 2
print(f"  s(1s->2s): {s_12:.4f}  [from Li]")

# From B: 2*s(1->2p) + 2*s(2s->2p) = 5 - 2.4214 = 2.5786
# We have 2 unknowns. Assume s(1->2p) = s(1->2s) = s_12:
s_12p = s_12  # assume same
s_2s2p = (5 - atoms['B']['Z_eff'] - 2*s_12p) / 2
print(f"  s(1s->2p): {s_12p:.4f}  [assumed = s(1s->2s)]")
print(f"  s(2s->2p): {s_2s2p:.4f}  [from B]")

# Same-shell screening from C-F increments
s_2p2p_vals = []
prev_Ze = atoms['B']['Z_eff']
for name in ['C', 'N', 'O', 'F']:
    Ze = atoms[name]['Z_eff']
    s_pp = 1 - (Ze - prev_Ze)
    s_2p2p_vals.append(s_pp)
    prev_Ze = Ze

s_2p2p = np.mean(s_2p2p_vals)
print(f"  s(2p->2p): {s_2p2p:.4f}  [from C-F increments, avg of {[f'{v:.4f}' for v in s_2p2p_vals]}]")

# Now predict all 2p atoms:
print(f"\n  2p predictions with s(1->2)={s_12:.4f}, s(2s->2p)={s_2s2p:.4f}, s(2p->2p)={s_2p2p:.4f}:")
print(f"  {'Atom':>4} {'Z':>3} {'Z_pred':>7} {'Z_real':>7} {'err':>7}")
for name in ['B', 'C', 'N', 'O', 'F']:
    Z = atoms[name]['Z']
    n_2p_others = Z - 5  # other 2p electrons
    Ze_pred = Z - 2*s_12 - 2*s_2s2p - n_2p_others * s_2p2p
    Ze_real = atoms[name]['Z_eff']
    print(f"  {name:>4} {Z:3d} {Ze_pred:7.4f} {Ze_real:7.4f} {Ze_pred - Ze_real:+7.4f}")

# For Na (3s): 1s^2 + 2s^2 + 2p^6 all screen 3s
# Need s(1->3), s(2->3) for both s and p
# From Na: Z_eff = 11 - 2*s(1->3) - 2*s(2s->3) - 6*s(2p->3) = 2.5074
# 3 unknowns, 1 equation. Need more constraints.

# Harmonic hypothesis: s(n_i -> n_v) = 1 - (n_i/n_v)^p * geometric_factor
print()
print("=" * 80)
print("  HARMONIC HYPOTHESIS: s = 1 - g * (n_i/n_v)^p")
print("=" * 80)

# For s(1->2) = 0.8604: 1 - g*(1/2)^p = 0.8604, so g*(1/2)^p = 0.1396
# For s(2s->2p) = 0.4289: 1 - g*(2/2)^p = 0.4289, so g = 0.5711 (when n_i = n_v)
# But wait, 2s->2p is same n, different l. The ratio n_i/n_v = 1.

# So g = 0.5711 and (1/2)^p = 0.1396/0.5711 = 0.2444
# (1/2)^p = 0.2444, p = -log(0.2444)/log(2) = 2.033 !
# p ~ 2 !

g = 1 - s_2s2p  # g when n_i/n_v = 1 (but different l)
# Wait, s_2s_2p is for same-n, different-l. Same-l would be s_2p_2p.
# For same subshell: s_2p_2p = 0.3305, so g_same = 1 - 0.3305 = 0.6695
# For different l in same n: s_2s_2p = 0.4289, so g_diff_l = 1 - 0.4289 = 0.5711

# Let's separate: for same n, the coupling depends on l matching
# For different n, the coupling falls off as (n_i/n_v)^p

# From s(1->2) = 0.8604 and g ~ 0.57 (diff l, using s_2s->2p):
# g * (1/2)^p = 1 - 0.8604 = 0.1396
# 0.5711 * (1/2)^p = 0.1396
# (1/2)^p = 0.2444
# p = log(1/0.2444) / log(2) = 2.033

ratio_12 = 0.5  # n_i/n_v = 1/2
coupling_12 = 1 - s_12  # = 0.1396
coupling_same_n_diff_l = 1 - s_2s2p  # = 0.5711
coupling_same_subshell = 1 - s_2p2p  # = 0.6695

p_harmonic = np.log(coupling_12 / coupling_same_n_diff_l) / np.log(ratio_12)

print(f"\n  Coupling at different distances:")
print(f"    Same subshell (2p->2p): g = 1 - s = {coupling_same_subshell:.4f}")
print(f"    Same n, diff l (2s->2p): g = 1 - s = {coupling_same_n_diff_l:.4f}")
print(f"    Different n (1s->2p): g = 1 - s = {coupling_12:.4f}")
print(f"\n  Harmonic falloff: g(1->2) / g(2s->2p) = (1/2)^p")
print(f"    p = {p_harmonic:.4f}  (very close to 2!)")
print(f"\n  This means: coupling falls off as (n_i/n_v)^2")
print(f"  Which is the SQUARE of the harmonic ratio!")
print(f"  In wave terms: overlap integral ~ (wavelength_inner/wavelength_outer)^2")

# =============================================================================
# FULL HARMONIC SCREENING MODEL
# =============================================================================
print()
print("=" * 80)
print("  FULL HARMONIC SCREENING MODEL")
print("=" * 80)

# Model: s(n_i, l_i -> n_v, l_v) = 1 - g(l_i, l_v) * (n_i/n_v)^2
# Where g depends on orbital type matching:
#   g_same = coupling for same subshell (l_i = l_v, same n)
#   g_diff = coupling for different l (same or different n)

# From data:
# g_same = 1 - s_2p2p = 0.6695
# g_diff_l = 1 - s_2s2p = 0.5711
# Check: g_diff_l * (1/2)^2 = 0.5711 * 0.25 = 0.1428
#         1 - s_12 = 0.1396 (close!)

# Let's look for GWT constants:
print(f"\n  g_same_subshell = {coupling_same_subshell:.4f}")
print(f"    compare: 2/d = {2/dd:.4f}")
print(f"    compare: (d-1)/(d-1+1/d) = {(dd-1)/(dd-1+1/dd):.4f}")
print(f"    compare: 2/(d+1/(d-1)) = hmm")

print(f"\n  g_diff_l = {coupling_same_n_diff_l:.4f}")
print(f"    compare: 1/(d-1) = {1/(dd-1):.4f}")
print(f"    compare: (d-1)/d = {(dd-1)/dd:.4f}")
print(f"    compare: (2d-1)/(2d+1) = {(2*dd-1)/(2*dd+1):.4f}")
print(f"    compare: 4/(2d+1) = {4/(2*dd+1):.4f}")

# Let me try: g_same = 2/3, g_diff = 4/7
# Then s(1->2, diff l) = 1 - 4/7 * (1/2)^2 = 1 - 4/28 = 1 - 1/7 = 6/7 = 0.857
# Actual: 0.860 - close!

# s(2s->2p, same n diff l) = 1 - 4/7 * 1 = 3/7 = 0.429
# Actual: 0.429 - exact!

# s(2p->2p, same subshell) = 1 - 2/3 * 1 = 1/3 = 0.333
# Actual: 0.331 - close!

g_same_test = 2/dd  # 2/3
g_diff_test = 4/(2*dd+1)  # 4/7

print(f"\n  TEST: g_same = 2/d = {g_same_test:.4f}, g_diff = 4/(2d+1) = {g_diff_test:.4f}")
print(f"  Formula: s(n_i -> n_v) = 1 - g * (n_i/n_v)^2")

# Predict all atoms
print(f"\n  {'Atom':>4} {'Z':>3} {'Z_pred':>7} {'Z_real':>7} {'err':>7} {'E_pred':>8} {'E_real':>8} {'E_err%':>7}")
print("  " + "-" * 70)

def harmonic_zeff(name):
    info = atoms[name]
    Z = info['Z']; n_v = info['n']; l_v = info['l']

    config = electron_configs[name]
    S = 0  # total screening
    for (n_i, l_i, count) in config:
        if n_i == n_v and l_i == l_v:
            # Same subshell
            g = g_same_test
        else:
            # Different subshell (same or different n)
            g = g_diff_test
        s = 1 - g * (n_i / n_v)**2
        S += count * s
    return Z - S

for name in atoms:
    info = atoms[name]
    Ze_pred = harmonic_zeff(name)
    Ze_real = info['Z_eff']
    n = info['n']
    E_pred = E_H * (Ze_pred/n)**2
    E_real = E_H * (Ze_real/n)**2
    print(f"  {name:>4} {info['Z']:3d} {Ze_pred:7.4f} {Ze_real:7.4f} {Ze_pred-Ze_real:+7.4f} "
          f"{E_pred:8.4f} {E_real:8.4f} {(E_pred-E_real)/E_real*100:+6.1f}%")


# =============================================================================
# SCAN: What are the best g_same and g_diff from d?
# =============================================================================
print()
print("=" * 80)
print("  SCANNING g_same and g_diff for best fit")
print("=" * 80)

# Test all combinations of GWT-motivated constants
candidates = {
    '1/d': 1/dd,
    '2/d': 2/dd,
    '1/(d-1)': 1/(dd-1),
    '(d-1)/d': (dd-1)/dd,
    '2/(2d+1)': 2/(2*dd+1),
    '4/(2d+1)': 4/(2*dd+1),
    '1/(d+1)': 1/(dd+1),
    '2/(d+1)': 2/(dd+1),
    'd/(d+1)': dd/(dd+1),
    'f_pi': f_pi,
    '1-1/d': 1-1/dd,
    '1/pi': 1/pi,
    '2/pi': 2/pi,
    '(d-1)/(d+1)': (dd-1)/(dd+1),
}

best = (999, '', '', 0, 0)
print(f"\n{'g_same':>15} {'g_diff':>15} {'rms_Ze':>8} {'rms_E%':>8}")
print("-" * 55)

for gs_name, gs_val in candidates.items():
    for gd_name, gd_val in candidates.items():
        errs_ze = []
        errs_E = []
        for name in atoms:
            info = atoms[name]
            Z = info['Z']; n_v = info['n']; l_v = info['l']
            config = electron_configs[name]
            S = 0
            for (n_i, l_i, count) in config:
                if n_i == n_v and l_i == l_v:
                    g = gs_val
                else:
                    g = gd_val
                s = 1 - g * (n_i / n_v)**2
                S += count * s
            Ze_pred = Z - S
            Ze_real = info['Z_eff']
            errs_ze.append((Ze_pred - Ze_real)**2)
            if Ze_pred > 0:
                E_pred = E_H * (Ze_pred/n_v)**2
                E_real = E_H * (Ze_real/n_v)**2
                errs_E.append(((E_pred - E_real)/E_real * 100)**2)
            else:
                errs_E.append(1e6)

        rms_ze = np.sqrt(np.mean(errs_ze))
        rms_E = np.sqrt(np.mean(errs_E))

        if rms_ze < 0.35:
            print(f"{gs_name:>15} {gd_name:>15} {rms_ze:8.4f} {rms_E:8.2f}%")

        if rms_ze < best[0]:
            best = (rms_ze, gs_name, gd_name, gs_val, gd_val)

print(f"\n  BEST: g_same={best[1]}={best[3]:.6f}, g_diff={best[2]}={best[4]:.6f}, RMS(Z_eff)={best[0]:.4f}")


# =============================================================================
# ALSO TRY: p != 2 (scan the harmonic power)
# =============================================================================
print()
print("=" * 80)
print("  SCANNING harmonic power p in s = 1 - g * (n_i/n_v)^p")
print("=" * 80)

best_overall = (999, 0, 0, 0)
for p in np.arange(1.0, 4.1, 0.1):
    for gs in np.arange(0.3, 0.9, 0.02):
        for gd in np.arange(0.3, 0.9, 0.02):
            errs = []
            for name in atoms:
                info = atoms[name]
                Z = info['Z']; n_v = info['n']; l_v = info['l']
                config = electron_configs[name]
                S = 0
                for (n_i, l_i, count) in config:
                    if n_i == n_v and l_i == l_v:
                        g = gs
                    else:
                        g = gd
                    s = 1 - g * (n_i / n_v)**p
                    S += count * s
                Ze_pred = Z - S
                errs.append((Ze_pred - atoms[name]['Z_eff'])**2)
            rms = np.sqrt(np.mean(errs))
            if rms < best_overall[0]:
                best_overall = (rms, p, gs, gd)

p_best, gs_best, gd_best = best_overall[1], best_overall[2], best_overall[3]
print(f"\n  Best fit: p={p_best:.1f}, g_same={gs_best:.4f}, g_diff={gd_best:.4f}, RMS={best_overall[0]:.4f}")

# Check against GWT constants
print(f"\n  p = {p_best:.1f}")
print(f"    p=2 means (n_i/n_v)^2 = energy ratio!")
print(f"\n  g_same = {gs_best:.4f}")
for cname, cval in candidates.items():
    if abs(cval - gs_best) < 0.03:
        print(f"    ~ {cname} = {cval:.4f}")
print(f"\n  g_diff = {gd_best:.4f}")
for cname, cval in candidates.items():
    if abs(cval - gd_best) < 0.03:
        print(f"    ~ {cname} = {cval:.4f}")

# Show predictions with best values
print(f"\n  Predictions with best (p={p_best}, g_s={gs_best}, g_d={gd_best}):")
print(f"  {'Atom':>4} {'Z':>3} {'Z_pred':>7} {'Z_real':>7} {'err':>7}")
print("  " + "-" * 35)
for name in atoms:
    info = atoms[name]
    Z = info['Z']; n_v = info['n']; l_v = info['l']
    config = electron_configs[name]
    S = 0
    for (n_i, l_i, count) in config:
        if n_i == n_v and l_i == l_v:
            g = gs_best
        else:
            g = gd_best
        s = 1 - g * (n_i / n_v)**p_best
        S += count * s
    Ze_pred = Z - S
    print(f"  {name:>4} {Z:3d} {Ze_pred:7.4f} {info['Z_eff']:7.4f} {Ze_pred-info['Z_eff']:+7.4f}")
