"""
Whole-wave Z_eff model: treat the atom's standing wave as ONE structure.

Instead of summing per-electron screening, treat each SHELL as a resonant
unit. A complete shell forms a coherent spherical wave pattern that screens
differently than the sum of individual modes.

Key idea: The screening of a complete shell depends on its WAVE STRUCTURE
as a whole, not on counting individual electrons.
"""
import numpy as np
from math import pi, sqrt

d = 3
g_same = 2/d       # 2/3
g_diff = 4/(2*d+1) # 4/7

# Clementi-Raimondi Z_eff values
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

# Shell structure: each shell as a UNIT
# Shell capacity: n=1 → 2 electrons, n=2 → 8, n=3 → 18
# A "complete shell" has all angular modes filled
shells = {
    'H':  [(1, 1, 2)],    # (n_shell, N_electrons, capacity)
    'Li': [(1, 2, 2), (2, 1, 8)],
    'B':  [(1, 2, 2), (2, 3, 8)],  # 2 core + 3 valence in n=2
    'C':  [(1, 2, 2), (2, 4, 8)],
    'N':  [(1, 2, 2), (2, 5, 8)],
    'O':  [(1, 2, 2), (2, 6, 8)],
    'F':  [(1, 2, 2), (2, 7, 8)],
    'Na': [(1, 2, 2), (2, 8, 8), (3, 1, 18)],
    'Cl': [(1, 2, 2), (2, 8, 8), (3, 7, 18)],
}

print("="*70)
print("  MODEL 1: SHELL-BASED SCREENING")
print("  Each shell screens as a unit based on its completeness fraction")
print("="*70)

# Idea: a shell with N/N_max electrons screens with effectiveness
# that depends on its filling fraction f = N/N_max
# Complete shell (f=1) screens almost perfectly
# Partial shell screens less
# s_shell = 1 - g * (n_shell/n_outer)^2 * h(f)
# where h(f) modifies the coupling based on completeness

# h(f=1) should be small (complete shell → strong screening → small leakage)
# h(f→0) should be 1 (single electron → normal harmonic model)

# Try: h(f) = 1 - c*f (linear in filling)
# or: h(f) = (1-f)^a
# or: h(f) = 1/(1 + c*f)

print("\n--- Model 1a: leakage = g*(n/nv)^2 * (1 - c*f_fill) ---")
print("  where f_fill = (N_in_shell - 1) / (capacity - 1) for same shell")
print("  or f_fill = N_in_shell / capacity for inner shells")

best_rms = 999
best_c = 0
for c100 in range(0, 200):
    c = c100 / 100
    rms = 0
    for name, at in atoms.items():
        Z, nv = at['Z'], at['n']
        total_s = 0
        for (ns, Ne, cap) in shells[name]:
            if ns == nv:
                # Same shell: self is counted, so N_other = Ne - 1
                if Ne <= 1:
                    continue
                f = (Ne - 1) / (cap - 1) if cap > 1 else 0
                # Same shell has mixed same/diff subshell
                # Average g for the shell
                g_avg = (g_same + g_diff) / 2  # rough average
                leak = g_avg * 1 * (1 - c * f)  # (ns/nv)^2 = 1
                s_per = 1 - leak
                total_s += s_per * (Ne - 1)
            else:
                # Inner shell
                f = Ne / cap
                r2 = (ns / nv) ** 2
                leak = g_diff * r2 * (1 - c * f)
                s_per = 1 - leak
                total_s += s_per * Ne
        pred = Z - total_s
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_c = c

print(f"  Best c = {best_c:.2f}, RMS = {best_rms:.4f}")
for name, at in atoms.items():
    Z, nv = at['Z'], at['n']
    total_s = 0
    for (ns, Ne, cap) in shells[name]:
        if ns == nv:
            if Ne <= 1: continue
            f = (Ne - 1) / (cap - 1) if cap > 1 else 0
            g_avg = (g_same + g_diff) / 2
            leak = g_avg * 1 * (1 - best_c * f)
            s_per = 1 - leak
            total_s += s_per * (Ne - 1)
        else:
            f = Ne / cap
            r2 = (ns / nv) ** 2
            leak = g_diff * r2 * (1 - best_c * f)
            s_per = 1 - leak
            total_s += s_per * Ne
    pred = Z - total_s
    print(f"  {name:3s}: f_inner={'-'}, pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")


print("\n" + "="*70)
print("  MODEL 2: GAUSS'S LAW IN WAVE SPACE")
print("  Complete inner wave = perfect shielding (s=1)")
print("  Incomplete wave = partial shielding (s < 1)")
print("="*70)

# In classical EM, a spherically symmetric charge distribution inside
# radius r screens perfectly at r (Gauss's law).
# For a wave function, the screening is imperfect because the wave
# extends beyond its "classical" radius.

# What if: s = 1 - leakage, where leakage is the fraction of the
# inner mode's wave that extends BEYOND the outer mode?
# For a hydrogen-like atom, this goes as (n_i/n_v)^(2l+2) or similar.

# But in GWT: the wave modes are STANDING WAVES in the atom.
# The "leakage" of inner mode n_i past the outer mode n_v
# is related to their spatial overlap.

# Key reframe: Instead of g*(n/nv)^2, what if the leakage depends
# on the NUMBER OF NODES between the two modes?
# n_v - n_i = number of radial nodes separating them

# More nodes = less leakage (more complete screening)

print("\n--- Model 2a: s = 1 - g * (1/delta_n)^p for inner shells ---")
print("  where delta_n = n_outer - n_inner")

# For delta_n = 1: s = 1 - g (maximal leakage)
# For delta_n = 2: s = 1 - g/2^p (less leakage)
# Same shell: delta_n = 0 → need separate treatment

# Subshell-level configs for same-shell screening
configs_sub = {
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

best_rms = 999
best_p = 0
for p100 in range(50, 400):
    p = p100 / 100
    rms = 0
    for name, at in atoms.items():
        Z, nv, lv = at['Z'], at['n'], at['l']
        total_s = 0
        for (ni, li, Ne) in configs_sub[name]:
            if ni == nv and li == lv:
                # Same subshell
                s_per = 1 - g_same  # (ni/nv)^2 = 1
                total_s += s_per * Ne
            elif ni == nv:
                # Same shell, different subshell
                s_per = 1 - g_diff
                total_s += s_per * Ne
            else:
                # Inner shell
                dn = nv - ni
                leak = g_diff / dn**p
                s_per = 1 - leak
                total_s += s_per * Ne
        pred = Z - total_s
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_p = p

print(f"  Best p = {best_p:.2f}, RMS = {best_rms:.4f}")
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs_sub[name]:
        if ni == nv and li == lv:
            s_per = 1 - g_same
            total_s += s_per * Ne
        elif ni == nv:
            s_per = 1 - g_diff
            total_s += s_per * Ne
        else:
            dn = nv - ni
            leak = g_diff / dn**best_p
            s_per = 1 - leak
            total_s += s_per * Ne
    pred = Z - total_s
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")


print("\n" + "="*70)
print("  MODEL 3: WHOLE-ATOM WAVE ENERGY")
print("  Z_eff emerges from the total wave pattern, not per-electron sums")
print("="*70)

# The atom's standing wave has total quantum numbers.
# What if Z_eff = Z * f(total_wave_structure)?
#
# For a hydrogenic atom (1 electron): Z_eff = Z
# For adding more electrons, each one modifies the TOTAL wave pattern.
#
# Simplest whole-atom idea:
# Z_eff = Z - (Z-1) * (1 - 1/n_v^2)
# This gives Z_eff → 1 as n → inf (all charge screened)
# and Z_eff = Z when n=1 (no screening)

print("\n--- Z_eff = Z - (Z-1)*(1 - a/n^b) ---")
best_rms = 999
for a100 in range(10, 300):
    a = a100 / 100
    for b100 in range(50, 300):
        b = b100 / 100
        rms = 0
        for name, at in atoms.items():
            pred = at['Z'] - (at['Z']-1) * (1 - a/at['n']**b)
            rms += (pred - at['Ze'])**2
        rms = sqrt(rms / len(atoms))
        if rms < best_rms:
            best_rms = rms
            best_a, best_b = a, b

print(f"  Best: a={best_a:.2f}, b={best_b:.2f}, RMS={best_rms:.4f}")
for name, at in atoms.items():
    pred = at['Z'] - (at['Z']-1) * (1 - best_a/at['n']**best_b)
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")

# What if Z_eff depends on the WAVE ENERGY RATIO between the outer mode
# and the total atom?
# Outer mode energy ∝ Z_eff^2 / n^2
# Total atom energy ∝ sum of Z_eff_i^2 / n_i^2
# This is self-referential... but what if we use bare Z for total energy?

print("\n--- Z_eff from energy balance ---")
# Total nuclear attraction: Z^2 (units of E_H)
# Each shell at n has energy ∝ Z^2/n^2 per electron
# The outer electron sees Z_eff based on its share of the total potential

# What if Z_eff = Z * sqrt(n_outer^2 / sum(N_i * n_outer^2 / n_i^2))?
# i.e., the outer electron's share of potential weighted by n ratios

for name, at in atoms.items():
    Z, nv = at['Z'], at['n']
    # Total "potential weight": sum of N_i * (n_v/n_i)^2
    total_weight = 0
    for (ns, Ne, cap) in shells[name]:
        total_weight += Ne * (nv / ns) ** 2
    # Outer electron sees: Z * 1/total_weight
    pred = Z / total_weight * nv**2 / nv**2  # simplifies
    pred2 = Z * 1 / total_weight
    print(f"  {name:3s}: weight={total_weight:.2f}, Ze=Z/w={pred2:.4f}, real={at['Ze']:.4f}")

print("\n--- Normalized: Z_eff = Z * n_v^2 / sum_i(N_i * n_v^2/n_i^2) ---")
# This doesn't account for self-screening correctly. Let me think differently.

# In a standing wave, each mode has a FREQUENCY proportional to energy.
# The coupling between modes depends on frequency ratio.
# For a hydrogen atom, mode frequencies are: f_n ∝ 1/n^2
# Frequency ratio: f_outer/f_inner = (n_i/n_v)^2

# What if the outer mode is "screened" by the BEAT FREQUENCY between it
# and each inner mode?

print("\n" + "="*70)
print("  MODEL 4: BEAT FREQUENCY SCREENING")
print("="*70)
# Beat between modes n_i and n_v: |1/n_i^2 - 1/n_v^2|
# For n_i < n_v: beat = 1/n_i^2 - 1/n_v^2
# Screening proportional to beat frequency?
# Higher beat = modes more separated = better screening

print("  Beat frequencies (1/ni^2 - 1/nv^2):")
for name, at in atoms.items():
    nv = at['n']
    beats = []
    for (ns, Ne, cap) in shells[name]:
        if ns < nv:
            beat = 1/ns**2 - 1/nv**2
            beats.append((ns, Ne, beat))
            print(f"    {name}: n={ns}, N={Ne}, beat={beat:.4f}")
    if not beats:
        print(f"    {name}: (no inner shells)")

# What if s = c * beat_frequency?
# For n_i=1, n_v=2: beat = 1 - 1/4 = 3/4, s = c*3/4
# For n_i=1, n_v=3: beat = 1 - 1/9 = 8/9, s = c*8/9
# For n_i=2, n_v=3: beat = 1/4 - 1/9 = 5/36, s = c*5/36
# Hmm, 5/36 is small — 2→3 screening would be weak. But actually the
# 2→3 electrons are the ones that need STRONG screening for Na.

# Actually beat frequency means DIFFERENT modes separate more cleanly
# → they interact LESS → less leakage → MORE screening
# So: s = 1 - g / beat^a ? or s = 1 - g * exp(-c*beat)?

print("\n--- s = 1 - g / (beat + epsilon) for inner shells ---")
best_rms = 999
for g100 in range(1, 100):
    g = g100 / 100
    for eps100 in range(100, 500):
        eps = eps100 / 100
        rms = 0
        for name, at in atoms.items():
            Z, nv, lv = at['Z'], at['n'], at['l']
            total_s = 0
            for (ni, li, Ne) in configs_sub[name]:
                if ni == nv and li == lv:
                    s_per = 1 - g_same
                elif ni == nv:
                    s_per = 1 - g_diff
                else:
                    beat = 1/ni**2 - 1/nv**2
                    s_per = 1 - g / (beat + eps)
                total_s += s_per * Ne
            pred = Z - total_s
            rms += (pred - at['Ze'])**2
        rms = sqrt(rms / len(atoms))
        if rms < best_rms:
            best_rms = rms
            best_g, best_eps = g, eps

print(f"  Best: g={best_g:.2f}, eps={best_eps:.2f}, RMS={best_rms:.4f}")
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs_sub[name]:
        if ni == nv and li == lv:
            s_per = 1 - g_same
        elif ni == nv:
            s_per = 1 - g_diff
        else:
            beat = 1/ni**2 - 1/nv**2
            s_per = 1 - best_g / (beat + best_eps)
        total_s += s_per * Ne
    pred = Z - total_s
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")


print("\n" + "="*70)
print("  MODEL 5: WAVE AMPLITUDE OVERLAP")
print("="*70)
# Key GWT insight: screening = wave amplitude overlap
# The AMPLITUDE of mode n at distance r goes as r^l * exp(-Z*r/(n*a0)) * L(r)
# At the outer mode's characteristic radius r_v ~ n_v^2/Z:
# The inner mode's amplitude ∝ exp(-Z*r_v/(n_i*a0)) = exp(-n_v^2/n_i)
#
# Leakage = fraction of inner mode beyond r_v ∝ exp(-n_v^2/n_i) approximately
# Or more precisely: exp(-(n_v/n_i)^2) — the exponential of the energy ratio!

print("  Leakage = exp(-(nv/ni)^a) — exponential suppression of overlap")
for a_try in [1.0, 1.5, 2.0]:
    print(f"\n  a = {a_try}:")
    for name, at in atoms.items():
        Z, nv, lv = at['Z'], at['n'], at['l']
        total_s = 0
        for (ni, li, Ne) in configs_sub[name]:
            if ni == nv and li == lv:
                s_per = 1 - g_same
            elif ni == nv:
                s_per = 1 - g_diff
            else:
                leak = np.exp(-(nv/ni)**a_try)
                s_per = 1 - leak
            total_s += s_per * Ne
        pred = Z - total_s
        err = pred - at['Ze']
        print(f"    {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={err:+.4f}")

# Scan a for best fit
best_rms = 999
for a100 in range(50, 500):
    a = a100 / 100
    rms = 0
    for name, at in atoms.items():
        Z, nv, lv = at['Z'], at['n'], at['l']
        total_s = 0
        for (ni, li, Ne) in configs_sub[name]:
            if ni == nv and li == lv:
                s_per = 1 - g_same
            elif ni == nv:
                s_per = 1 - g_diff
            else:
                leak = np.exp(-(nv/ni)**a)
                s_per = 1 - leak
            total_s += s_per * Ne
        pred = Z - total_s
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_a = a

print(f"\n  Best a = {best_a:.2f}, RMS = {best_rms:.4f}")
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs_sub[name]:
        if ni == nv and li == lv:
            s_per = 1 - g_same
        elif ni == nv:
            s_per = 1 - g_diff
        else:
            leak = np.exp(-(nv/ni)**best_a)
            s_per = 1 - leak
        total_s += s_per * Ne
    pred = Z - total_s
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")


print("\n" + "="*70)
print("  MODEL 6: THE ATOM AS A SINGLE RESONANT CAVITY")
print("="*70)
# Total wave function of the atom has:
# - Z protons (Z units of charge at center)
# - N_total = Z electrons (neutral atom) distributed in shells
# - The outer mode sees an EFFECTIVE charge = Z_eff
#
# In a resonant cavity, the mode frequencies shift due to boundary conditions.
# What if Z_eff is determined by the RESONANCE CONDITION of the whole atom?
#
# For a standing wave in a sphere, the resonance condition involves
# the total phase accumulated across the atom.
#
# Phase = integral of k(r) dr from 0 to infinity
# k(r) = sqrt(2m(E - V(r))) / hbar
# V(r) = -Z_eff(r) * e^2 / r
#
# But in GWT, we can think of it more simply:
# The outer mode's frequency = Z_eff^2 / n^2 (in Rydberg units)
# This frequency must RESONATE with the cavity formed by inner modes.

# Simplest resonance: Z_eff^2/n_v^2 = Z^2/n_v^2 - sum_inner(Z^2/n_i^2) * coupling
# i.e., the outer mode's energy is the nuclear energy minus the coupled inner energies

# What if: Z_eff^2 = Z^2 - c * sum_i(N_i * Z^2 * n_v^2/n_i^2) ?
# Nah, too many free parameters.

# SIMPLEST: What if Z_eff = Z - sum_i[N_i * s_i] where s_i depends on
# the ENERGY RATIO in a nonlinear way?

# We've been using s = 1 - g*(ni/nv)^2 (linear coupling)
# What if it's: s = 1 - g*(ni/nv)^2 / sqrt(1 + N_shell*(ni/nv)^2) ?
# i.e., the screening PER ELECTRON decreases as more electrons pile up
# No wait, it should INCREASE (collective effect)

# What if: s = 1 - g*(ni/nv)^2 / (1 + c*N_in_shell*(ni/nv)^2) ?
# More electrons in the shell → denominator grows → leakage shrinks → more screening

print("  Testing: collective screening enhancement")
print("  s = 1 - g*(ni/nv)^2 / (1 + c*N_shell*(ni/nv)^2)")

# For this we need to know which shell each subshell belongs to
shell_N = {
    'H':  {1: 1},
    'Li': {1: 2, 2: 1},
    'B':  {1: 2, 2: 5},
    'C':  {1: 2, 2: 6},
    'N':  {1: 2, 2: 7},
    'O':  {1: 2, 2: 8},
    'F':  {1: 2, 2: 9},
    'Na': {1: 2, 2: 10, 3: 1},
    'Cl': {1: 2, 2: 10, 3: 7},
}
# Wait, shell_N should be total electrons with n <= n_shell
# No, it should be electrons IN that shell
shell_pop = {
    'H':  {1: 1},
    'Li': {1: 2, 2: 1},
    'B':  {1: 2, 2: 3},  # 2s^2 + 2p^1
    'C':  {1: 2, 2: 4},
    'N':  {1: 2, 2: 5},
    'O':  {1: 2, 2: 6},
    'F':  {1: 2, 2: 7},
    'Na': {1: 2, 2: 8, 3: 1},
    'Cl': {1: 2, 2: 8, 3: 7},
}

best_rms = 999
best_c = 0
for c100 in range(0, 200):
    c = c100 / 100
    rms = 0
    for name, at in atoms.items():
        Z, nv, lv = at['Z'], at['n'], at['l']
        total_s = 0
        for (ni, li, Ne) in configs_sub[name]:
            if Ne == 0: continue
            r2 = (ni/nv)**2
            if ni == nv and li == lv:
                g = g_same
                N_sh = shell_pop[name].get(ni, 0)
            elif ni == nv:
                g = g_diff
                N_sh = shell_pop[name].get(ni, 0)
            else:
                g = g_diff
                N_sh = shell_pop[name].get(ni, 0)

            leak = g * r2 / (1 + c * N_sh * r2)
            s_per = 1 - leak
            total_s += s_per * Ne
        pred = Z - total_s
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_c = c

print(f"  Best c = {best_c:.2f}, RMS = {best_rms:.4f}")
for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs_sub[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        if ni == nv and li == lv:
            g = g_same
        elif ni == nv:
            g = g_diff
        else:
            g = g_diff
        N_sh = shell_pop[name].get(ni, 0)
        leak = g * r2 / (1 + best_c * N_sh * r2)
        s_per = 1 - leak
        total_s += s_per * Ne
    pred = Z - total_s
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")


print("\n" + "="*70)
print("  MODEL 7: RECURSIVE WAVE COUPLING")
print("="*70)
# Each shell sees Z_eff from the shells INSIDE it.
# Shell 1 sees Z_nuclear directly.
# Shell 2 sees Z_nuclear screened by shell 1's effect.
# Shell 3 sees Z_nuclear screened by shells 1+2.
#
# But the screening by shell k depends on Z_eff SEEN BY shell k,
# because a tighter-bound inner shell screens more effectively.

# Z_eff(n=1) = Z - s_11 * (N_1 - 1)
# Z_eff(n=2) = Z - N_1 * s_12(Z_eff_1) - s_22 * (N_2 - 1)
# Z_eff(n=3) = Z - N_1 * s_13 - N_2 * s_23(Z_eff_2) - s_33 * (N_3 - 1)

# The key idea: inner electrons that are MORE tightly bound (higher Z_eff_inner)
# are more compact → screen better

print("  Z_eff builds recursively from inside out")
print("  Inner Z_eff → tighter binding → better screening of outer")

# Simple version: s_inner->outer = 1 - g * (n_i/n_v)^2 * (Z_bare/Z_eff_inner)^a
# When Z_eff_inner ≈ Z_bare, screening is as before.
# When Z_eff_inner < Z_bare (screened inner electron), it's more diffuse → screens LESS
# Wait, that's backwards. A more tightly bound electron is MORE compact,
# so it's MORE likely to be INSIDE the outer electron → BETTER screening.
# So: s should INCREASE with Z_eff_inner.

# s = 1 - g * (n_i/n_v)^2 / (Z_eff_inner / n_i^2)^a
# or simpler: s = 1 - g * (n_i/n_v)^2 * (n_i^2 / Z_eff_inner)^a

# Let's try: compute Z_eff for each subshell from inside out
print("\n  Recursive Z_eff calculation:")

# Z_eff for 1s electrons (all atoms with Z>=2):
# Z_eff_1s = Z - s_same * 1 = Z - 1/3
# (one other 1s electron, s = 1 - g_same = 1/3)

for name, at in atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']

    # Build Z_eff from inside
    ze_1s = Z - g_same * 0 if Z == 1 else Z - (1 - g_same) * 1  # 1 other 1s
    # 1s orbital "size" ∝ 1/ze_1s

    total_s = 0
    for (ni, li, Ne) in configs_sub[name]:
        if Ne == 0: continue
        if ni == nv and li == lv:
            s_per = 1 - g_same
        elif ni == nv:
            s_per = 1 - g_diff
        else:
            # Use standard harmonic for now
            s_per = 1 - g_diff * (ni/nv)**2
        total_s += s_per * Ne

    pred = Z - total_s
    # (Same as harmonic for now - placeholder for recursive version)

# Actually, let me implement a proper recursive model
def recursive_zeff(name, at_data, configs, c_recursive):
    Z = at_data['Z']
    nv = at_data['n']
    lv = at_data['l']

    # Group configs by shell
    shell_configs = {}
    for (ni, li, Ne) in configs:
        if ni not in shell_configs:
            shell_configs[ni] = []
        shell_configs[ni].append((li, Ne))

    # Compute Z_eff for each shell from inside out
    ze_by_shell = {}
    screening_so_far = 0

    for n_sh in sorted(shell_configs.keys()):
        # This shell sees: Z minus screening from all inner shells
        # Plus self-screening within the shell
        N_in_shell = sum(Ne for (li, Ne) in shell_configs[n_sh])

        # Self-screening (within same shell)
        self_screen = 0
        for (li, Ne) in shell_configs[n_sh]:
            if n_sh == nv and li == lv:
                self_screen += (1 - g_same) * Ne  # same subshell
            elif n_sh == nv:
                self_screen += (1 - g_diff) * Ne   # same shell diff sub
            else:
                # For inner shells: self-screening among their own electrons
                self_screen += (1 - g_same) * Ne  # approximate

        ze_sh = Z - screening_so_far - self_screen
        ze_by_shell[n_sh] = max(ze_sh, 0.1)

        if n_sh < nv:
            # This shell screens the outer shell
            # Tighter inner electrons (higher ze_sh) screen MORE effectively
            # s = 1 - g*(ni/nv)^2 * correction(ze_inner)
            r2 = (n_sh / nv)**2

            # correction: inner electrons are more compact when ze_inner is larger
            # relative compactness = ze_inner / Z (fraction of nuclear charge seen)
            # When ze_inner/Z → 1, electrons are tight, screen well
            # When ze_inner/Z → 0, electrons are diffuse, screen poorly
            compact = ze_sh / Z

            for (li, Ne) in shell_configs[n_sh]:
                g = g_diff  # cross-shell is always different
                leak = g * r2 * (1 - c_recursive * compact)
                leak = max(leak, 0)
                screening_so_far += (1 - leak) * Ne

    # Outer Z_eff
    # Need to add same-shell screening
    same_shell_screen = 0
    if nv in shell_configs:
        for (li, Ne) in shell_configs[nv]:
            if li == lv:
                same_shell_screen += (1 - g_same) * Ne
            else:
                same_shell_screen += (1 - g_diff) * Ne

    return Z - screening_so_far - same_shell_screen

best_rms = 999
best_c = 0
for c100 in range(-200, 200):
    c = c100 / 100
    rms = 0
    for name, at in atoms.items():
        pred = recursive_zeff(name, at, configs_sub[name], c)
        rms += (pred - at['Ze'])**2
    rms = sqrt(rms / len(atoms))
    if rms < best_rms:
        best_rms = rms
        best_c = c

print(f"\n  Best c_recursive = {best_c:.2f}, RMS = {best_rms:.4f}")
for name, at in atoms.items():
    pred = recursive_zeff(name, at, configs_sub[name], best_c)
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")


print("\n" + "="*70)
print("  SUMMARY: Best models comparison")
print("="*70)
print(f"  Harmonic (reference):         RMS ≈ 0.73")
print(f"  Harmonic errors concentrated in Cl (+2.08) and Na (+0.65)")
print(f"  Key question: what 'whole wave' effect are we missing?")
