"""
Z-dependent screening: inner electrons screen MORE in heavier atoms
because more protons pull them inward (more compact = better screening).

The harmonic model s = 1 - g*(ni/nv)^2 doesn't depend on Z at all.
But screening MUST depend on Z: in Cl (Z=17), the 10 inner electrons
are pulled much tighter than in Na (Z=11), so they screen better.

Key idea: the atom's wave structure is ONE pattern. The inner modes'
spatial extent depends on the TOTAL wave pattern, not just their
quantum numbers. Higher Z = tighter inner modes = better screening.
"""
import numpy as np
from math import pi, sqrt

d = 3
g_same = 2/d       # 2/3
g_diff = 4/(2*d+1) # 4/7

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

print("="*70)
print("  REQUIRED SCREENING PER INNER ELECTRON")
print("="*70)

# First: what screening per inner electron does reality REQUIRE?
# For B-F: inner = 1s^2, same-shell = rest of n=2
# For Na: inner = 1s^2 + 2s^2 + 2p^6, same-shell = none
# For Cl: inner = 1s^2 + 2s^2 + 2p^6, same-shell = 3s^2 + 3p^4

# Same-shell screening (trust harmonic for now):
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    inner_s = 0
    same_s = 0
    n_inner = 0
    for (ni, li, Ne) in configs[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        if ni == nv and li == lv:
            g = g_same
            same_s += (1 - g * r2) * Ne
        elif ni == nv:
            g = g_diff
            same_s += (1 - g * r2) * Ne
        else:
            n_inner += Ne
            # We want to know what screening is NEEDED, not predicted

    total_s_needed = Z - at['Ze']
    inner_s_needed = total_s_needed - same_s
    s_per_inner = inner_s_needed / n_inner if n_inner > 0 else 0

    # Compare with harmonic prediction
    inner_harmonic = 0
    for (ni, li, Ne) in configs[name]:
        if Ne == 0 or ni == nv: continue
        r2 = (ni/nv)**2
        g = g_same if (ni == nv and li == lv) else g_diff
        inner_harmonic += (1 - g * r2) * Ne

    s_per_harmonic = inner_harmonic / n_inner if n_inner > 0 else 0

    print(f"  {name:3s}: Z={Z:2d}, N_inner={n_inner:2d}, "
          f"s_needed/e={s_per_inner:.4f}, s_harmonic/e={s_per_harmonic:.4f}, "
          f"gap={s_per_inner - s_per_harmonic:+.4f}")


print("\n" + "="*70)
print("  WHAT DETERMINES THE GAP?")
print("="*70)
print("  The gap between needed and harmonic screening per inner electron:")
print("  Grows with Z? With N_inner/Z? With n?")

# For 2p atoms: gap varies from +0.004 (B) to -0.036 (N)
# These are SMALL. Period 2 is essentially solved.
# For Na: gap = +0.065. Moderate.
# For Cl: gap = +0.208. Large.
# What correlates?

# Na: Z=11, N_inner=10, Z/N_inner=1.1
# Cl: Z=17, N_inner=10, Z/N_inner=1.7
# Gap ratio: 0.208/0.065 = 3.2
# Z ratio: 17/11 = 1.55
# (Z/N)^2 ratio: (1.7/1.1)^2 = 2.39
# Z-N ratio: 7/1 = 7... nope

# What if the gap is proportional to (Z-N_inner-1)/N_inner * something?
# For Na: (11-10-1)/10 = 0 -> nah
# For Cl: (17-10-1)/10 = 0.6

# Actually, the gap might be proportional to Z_eff_inner - Z_eff_inner_hydrogen
# In heavy atoms, inner electrons see MORE nuclear charge than in light ones
# This makes them more compact, increasing screening

# What if we compute Z_eff for each INNER subshell, then use that?
print("\n--- Inner Z_eff analysis ---")
print("  Z_eff seen by each inner subshell (using harmonic for same-level screening):")

def compute_inner_zeff(name, at):
    """Compute Z_eff for each subshell from inside out."""
    Z = at['Z']
    results = {}

    cfg = configs[name]
    for idx, (n_out, l_out, N_out) in enumerate(cfg):
        if N_out == 0: continue
        total_screen = 0
        for jdx, (ni, li, Ne) in enumerate(cfg):
            if jdx == idx:
                # Same subshell: screen by N-1 electrons
                total_screen += (1 - g_same) * (Ne - 1) if Ne > 1 else 0
            elif Ne > 0:
                if ni == n_out and li == l_out:
                    continue  # skip self (shouldn't happen, handled above)
                elif ni == n_out:
                    # Same shell, different subshell
                    total_screen += (1 - g_diff) * Ne
                elif ni < n_out:
                    # Inner shell screens this subshell
                    r2 = (ni/n_out)**2
                    total_screen += (1 - g_diff * r2) * Ne
                else:
                    # Outer shell - doesn't screen
                    pass
        ze = Z - total_screen
        results[(n_out, l_out)] = ze
        print(f"    {name:3s} ({n_out},{l_out}): Z_eff = {ze:.4f} (N={N_out})")

    return results

for name, at in atoms.items():
    compute_inner_zeff(name, at)
    print()


print("="*70)
print("  MODEL A: Z-CORRECTED HARMONIC")
print("  s = 1 - g*(ni/nv)^2 * (Z_ref/Z)^a")
print("  where Z_ref normalizes so period 2 is unchanged")
print("="*70)

# The idea: in the harmonic model, the leakage g*(ni/nv)^2 represents
# the fraction of inner wave extending beyond the outer wave.
# In heavier atoms, the inner wave is more compact (higher Z pulls it in).
# So the leakage should DECREASE with Z.
# leakage = g * (ni/nv)^2 * (Z_ref/Z)^a
# For period 2, Z = 3-9. For period 3, Z = 11-17.

best_rms = 999
best_params = (0, 0)
for zref100 in range(100, 1500):
    zref = zref100 / 100
    for a100 in range(10, 300):
        a = a100 / 100
        rms = 0
        for name, at in atoms.items():
            Z, nv, lv = at['Z'], at['n'], at['l']
            total_s = 0
            for (ni, li, Ne) in configs[name]:
                if Ne == 0: continue
                r2 = (ni/nv)**2
                if ni == nv and li == lv:
                    g = g_same
                    leak = g * r2  # same shell: no Z correction
                elif ni == nv:
                    g = g_diff
                    leak = g * r2  # same shell: no Z correction
                else:
                    g = g_diff
                    leak = g * r2 * (zref/Z)**a
                s_per = 1 - leak
                total_s += s_per * Ne
            pred = Z - total_s
            rms += (pred - at['Ze'])**2
        rms = sqrt(rms / len(atoms))
        if rms < best_rms:
            best_rms = rms
            best_params = (zref, a)

zref, a = best_params
print(f"  Best: Z_ref={zref:.2f}, a={a:.2f}, RMS={best_rms:.4f}")
print(f"  Leakage scaled by (Z_ref/Z)^a")
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        if ni == nv and li == lv:
            leak = g_same * r2
        elif ni == nv:
            leak = g_diff * r2
        else:
            leak = g_diff * r2 * (zref/Z)**a
        total_s += (1 - leak) * Ne
    pred = Z - total_s
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")


print("\n" + "="*70)
print("  MODEL B: SAME-SHELL SCREENING ALSO DEPENDS ON Z")
print("  s_same = 1 - g*(Z_ref/Z)^b for same-shell")
print("="*70)

best_rms = 999
for zref100 in range(100, 1500):
    zref = zref100 / 100
    for a100 in range(10, 300):
        a = a100 / 100
        for b100 in range(0, 200):
            b = b100 / 100
            rms = 0
            for name, at in atoms.items():
                Z, nv, lv = at['Z'], at['n'], at['l']
                total_s = 0
                for (ni, li, Ne) in configs[name]:
                    if Ne == 0: continue
                    r2 = (ni/nv)**2
                    if ni == nv and li == lv:
                        leak = g_same * r2 * (zref/Z)**b
                    elif ni == nv:
                        leak = g_diff * r2 * (zref/Z)**b
                    else:
                        leak = g_diff * r2 * (zref/Z)**a
                    total_s += (1 - leak) * Ne
                pred = Z - total_s
                rms += (pred - at['Ze'])**2
            rms = sqrt(rms / len(atoms))
            if rms < best_rms:
                best_rms = rms
                best_params = (zref, a, b)

zref, a, b = best_params
print(f"  Best: Z_ref={zref:.2f}, a={a:.2f}, b={b:.2f}, RMS={best_rms:.4f}")
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        if ni == nv and li == lv:
            leak = g_same * r2 * (zref/Z)**b
        elif ni == nv:
            leak = g_diff * r2 * (zref/Z)**b
        else:
            leak = g_diff * r2 * (zref/Z)**a
        total_s += (1 - leak) * Ne
    pred = Z - total_s
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")


print("\n" + "="*70)
print("  MODEL C: RECURSIVE INNER Z_eff")
print("  Inner electrons are more compact when they see higher Z_eff")
print("  More compact = more screening of outer shell")
print("="*70)

# Build Z_eff from inside out:
# 1. Compute Z_eff for n=1 electrons
# 2. Use that to determine how compact they are
# 3. Their compactness determines their screening of n=2
# 4. Compute Z_eff for n=2 electrons
# 5. Use that for screening of n=3
# etc.

# Screening per electron: s = 1 - g*(ni/nv)^2 * correction
# correction = (ni^2 / Z_eff_inner)^c — less compact inner = less screening
# When Z_eff_inner is large: inner electron is tight, correction is small, s is large
# When Z_eff_inner is small: inner electron is diffuse, correction is large, s is small

def recursive_model(c_param):
    results = {}
    for name, at in atoms.items():
        Z, nv, lv = at['Z'], at['n'], at['l']

        # Group by shell
        shell_data = {}  # n -> list of (l, N)
        for (ni, li, Ne) in configs[name]:
            if ni not in shell_data:
                shell_data[ni] = []
            shell_data[ni].append((li, Ne))

        sorted_shells = sorted(shell_data.keys())
        ze_inner = {}  # Z_eff per shell
        inner_screen_total = 0

        for idx, n_sh in enumerate(sorted_shells):
            # Compute Z_eff for this shell
            # It sees: Z - screening from all inner shells - self screening

            self_screen = 0
            N_in_shell = sum(Ne for (li, Ne) in shell_data[n_sh])
            for (li, Ne) in shell_data[n_sh]:
                # Self-screening within shell
                self_screen += (1 - g_same) * max(Ne - 1, 0)  # simplified
                # Cross-subshell in same shell
                for (lj, Nj) in shell_data[n_sh]:
                    if lj != li:
                        self_screen += (1 - g_diff) * Nj
            # Avoid double counting
            self_screen = self_screen / max(N_in_shell, 1)  # per electron average

            ze_sh = Z - inner_screen_total - self_screen * (N_in_shell - 1)
            ze_inner[n_sh] = max(ze_sh, 0.1)

            if n_sh < nv:
                # This shell screens the outer shell
                for (li, Ne) in shell_data[n_sh]:
                    r2 = (n_sh / nv) ** 2
                    g = g_diff  # inner to outer is always different shell

                    # Compactness correction: inner orbital size ~ n_i^2 / Z_eff_inner
                    # Outer orbital size ~ n_v^2 / Z_outer (but Z_outer is what we want)
                    # Leakage ~ (inner_size / outer_size) = (n_i/n_v)^2 * (Z_outer/Z_eff_inner)
                    # But we don't know Z_outer yet. Use Z as proxy.

                    # Simple: leakage scales as (n_i^2 / Z_eff_inner)^c_param
                    # Higher Z_eff_inner = smaller orbital = less leakage = more screening
                    if c_param == 0:
                        leak = g * r2
                    else:
                        # Normalize so that for hydrogen-like (Ze=Z), we get standard result
                        leak = g * r2 * (n_sh**2 / ze_inner[n_sh]) ** c_param

                    s_per = 1 - leak
                    inner_screen_total += s_per * Ne

        # Same-shell screening
        same_screen = 0
        if nv in shell_data:
            for (li, Ne) in shell_data[nv]:
                if li == lv:
                    same_screen += (1 - g_same) * Ne
                else:
                    same_screen += (1 - g_diff) * Ne

        pred = Z - inner_screen_total - same_screen
        results[name] = pred

    return results

best_rms = 999
best_c = 0
for c1000 in range(0, 2000):
    c = c1000 / 1000
    results = recursive_model(c)
    rms = 0
    for name in atoms:
        rms += (results[name] - atoms[name]['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_c = c

print(f"  Best c = {best_c:.3f}, RMS = {best_rms:.4f}")
results = recursive_model(best_c)
for name, at in atoms.items():
    print(f"  {name:3s}: pred={results[name]:.4f}, real={at['Ze']:.4f}, err={results[name]-at['Ze']:+.4f}")

# Also show with c=0 (standard harmonic) for comparison
print(f"\n  For reference, c=0 (standard harmonic-like):")
results0 = recursive_model(0)
for name, at in atoms.items():
    print(f"  {name:3s}: pred={results0[name]:.4f}, real={at['Ze']:.4f}, err={results0[name]-at['Ze']:+.4f}")


print("\n" + "="*70)
print("  MODEL D: SIMPLEST POSSIBLE Z-DEPENDENT")
print("  What if leakage = g * (ni/(nv*Z))^2 * Z_ref^2 ?")
print("  i.e., the coupling ratio uses Z-SCALED quantum numbers")
print("="*70)

# In GWT, the wave mode's spatial extent ~ n^2 / Z (Bohr model)
# The RATIO of extents is (n_i^2/Z_i) / (n_v^2/Z_v)
# If inner and outer electrons see DIFFERENT Z_eff values,
# the ratio changes from (ni/nv)^2 to (ni^2/Zi) / (nv^2/Zv) = (ni/nv)^2 * (Zv/Zi)

# For simplicity, what if we use:
# Zi = Z for inner (they see nearly full nuclear charge)
# Zv = Z_eff for outer (what we're trying to calculate)

# This is self-referential! But we can solve self-consistently:
# Z_eff = Z - sum_i [s_i * N_i]
# s_i = 1 - g * (ni/nv)^2 * (Z_eff / Z)   for inner shells

# Rearranging: Z_eff = Z - sum_i [(1 - g*(ni/nv)^2*(Z_eff/Z)) * N_i]
# Let x = Z_eff:
# x = Z - sum_inner[(1 - g*r2_i*x/Z)*N_i] - sum_same[(1-g_s)*N_s]
# x = Z - sum_inner[N_i - g*r2_i*N_i*x/Z] - S_same
# x = Z - N_inner + (g/Z)*x*sum_inner[r2_i*N_i] - S_same
# x = (Z - N_inner - S_same) + x * (g/Z) * sum_inner[r2_i * N_i]
# x * (1 - (g/Z)*sum[r2*N]) = Z - N_inner - S_same
# x = (Z - N_inner - S_same) / (1 - (g/Z)*sum[r2*N])

print("  Self-consistent solution:")
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']

    N_inner = 0
    sum_r2N = 0
    S_same = 0

    for (ni, li, Ne) in configs[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        if ni == nv and li == lv:
            S_same += (1 - g_same) * Ne
        elif ni == nv:
            S_same += (1 - g_diff) * Ne
        else:
            N_inner += Ne
            sum_r2N += r2 * Ne

    numerator = Z - N_inner - S_same
    denominator = 1 - (g_diff / Z) * sum_r2N

    if abs(denominator) > 0.001:
        pred = numerator / denominator
    else:
        pred = float('inf')

    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# What about different powers?
print("\n--- Self-consistent with (Z_eff/Z)^a ---")
# s_i = 1 - g*(ni/nv)^2 * (x/Z)^a
# Harder to solve analytically for a != 1, iterate instead

for a_try in [0.5, 1.0, 1.5, 2.0]:
    print(f"\n  a = {a_try}:")
    for name, at in atoms.items():
        Z, nv, lv = at['Z'], at['n'], at['l']

        S_same = 0
        inner_configs = []
        for (ni, li, Ne) in configs[name]:
            if Ne == 0: continue
            r2 = (ni/nv)**2
            if ni == nv and li == lv:
                S_same += (1 - g_same) * Ne
            elif ni == nv:
                S_same += (1 - g_diff) * Ne
            else:
                inner_configs.append((ni, li, Ne, r2))

        # Iterate to self-consistency
        x = Z / 2  # initial guess
        for _ in range(100):
            inner_screen = 0
            for (ni, li, Ne, r2) in inner_configs:
                leak = g_diff * r2 * (max(x, 0.1) / Z) ** a_try
                inner_screen += (1 - leak) * Ne
            x_new = Z - inner_screen - S_same
            if abs(x_new - x) < 1e-8:
                break
            x = x_new

        print(f"    {name:3s}: pred={x:.4f}, real={at['Ze']:.4f}, err={x-at['Ze']:+.4f}")


print("\n" + "="*70)
print("  MODEL E: WAVE OVERLAP (WHOLE ATOM)")
print("  The outer mode's leakage into inner region creates screening")
print("  But inner modes' leakage OUT reduces screening")
print("  Net screening = 1 - (inner leakage out) + (outer leakage in)")
print("="*70)

# In the standard model: s_inner = 1 - g*(ni/nv)^2
# The g*(ni/nv)^2 term is the inner mode's "tail" extending into the
# outer region. But the outer mode ALSO has a "tail" extending into
# the inner region. This INCREASES the charge the outer electron
# sees (anti-screening).

# However, in GWT, the outer mode IS the electron whose Z_eff we want.
# The inner modes are what create the screening potential.

# What if screening depends on BOTH the inner mode's leakage AND
# the outer mode's penetration?
#
# Inner leakage out: g*(ni/nv)^2 (reduces screening)
# Outer penetration in: some function of (ni/nv)
#
# For the outer electron in the inner region, it sees FULL nuclear charge
# minus NOTHING (the inner electrons are also there, but they're a negative
# charge that partially cancels). This is the penetration effect.

# In quantum mechanics, s-orbitals penetrate more than p-orbitals.
# This is why 2s has lower energy than 2p in multi-electron atoms.
# Penetration depends on l of the OUTER electron.

# For the outer electron with l_v:
# Penetration probability ~ (Z*a0/n_v)^(2*l_v+1) or something
# More precisely, for l=0 (s): significant penetration
# For l=1 (p): less penetration
# For l=2 (d): even less

# This could explain why Na (3s, l=0) has different screening than
# Cl (3p, l=1) even with same inner shells!

print("  The outer electron's angular momentum affects how much it")
print("  penetrates inner shells:")
print("  l=0 (s): maximum penetration (sees more nuclear charge)")
print("  l=1 (p): less penetration")
print("  l=2 (d): even less")
print()

# What if: s = 1 - g*(ni/nv)^2 + penetration_correction(l_v)
# where penetration makes Z_eff LARGER (more nuclear charge seen)
# penetration ~ c / (2*l_v + 1)  (s-orbitals penetrate most)

# For Na (l=0): correction = c/1 = c (large extra Z_eff)
# For Cl (l=1): correction = c/3 (smaller extra Z_eff)

# This would make Na Z_eff larger than Cl for same inner structure

# But wait: Na real Z_eff = 2.51, Cl real Z_eff = 4.89
# Na has Z=11 with 10 inner, Cl has Z=17 with 10 inner + 6 same-shell
# So it's not directly comparable.

# Let's just scan: Z_eff = harmonic + c/(2*l_v+1)
print("--- Z_eff = harmonic_pred + c/(2*l_v+1) ---")
harmonic = {}
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        if ni == nv and li == lv:
            total_s += (1 - g_same * r2) * Ne
        elif ni == nv:
            total_s += (1 - g_diff * r2) * Ne
        else:
            total_s += (1 - g_diff * r2) * Ne
    harmonic[name] = Z - total_s

best_rms = 999
for c100 in range(-500, 500):
    c = c100 / 100
    rms = 0
    for name, at in atoms.items():
        pred = harmonic[name] + c / (2*at['l'] + 1)
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_c = c

print(f"  Best c = {best_c:.2f}, RMS = {best_rms:.4f}")
for name, at in atoms.items():
    pred = harmonic[name] + best_c / (2*at['l'] + 1)
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")


print("\n" + "="*70)
print("  SUMMARY OF ALL MODELS")
print("="*70)
print(f"  Harmonic (g=2/d, 4/(2d+1), p=2):     RMS ~ 0.73")
print(f"  Best whole-wave models above will be compared here")
