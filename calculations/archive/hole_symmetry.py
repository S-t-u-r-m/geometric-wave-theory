"""
Hole symmetry in GWT: Cl as a 'flipped' harmonic pattern.

Key insight from user: Cl's 3p^5 is ONE HOLE away from a complete 3p^6 shell.
A complete shell is a 'perfect harmonic' (all angular modes filled = spherically
symmetric wave). If we 'flip' the picture and treat Cl as having a HOLE instead
of 5 electrons, the wave pattern might be simpler.

Ar (Z=18) = perfect harmonic (all shells complete through n=3)
Cl (Z=17) = Ar - 1 proton - 1 electron
Na (Z=11) = Ne + 1 proton + 1 electron

Both Na and Cl are ONE STEP away from a noble gas.
Na = Ne + 1 (one extra)
Cl = Ar - 1 (one missing)

The 'particle-hole symmetry' in wave terms:
- Na's outer electron is ONE mode above a complete wave
- Cl's outer hole is ONE mode below a complete wave

What if Z_eff for near-noble-gas atoms can be computed from the
noble gas structure + a perturbation?
"""
import numpy as np
from math import pi, sqrt

d = 3
g_same = 2/d
g_diff = 4/(2*d+1)

atoms = {
    'H':  {'Z': 1,  'n': 1, 'l': 0, 'Ze': 1.0000},
    'He': {'Z': 2,  'n': 1, 'l': 0, 'Ze': 1.6875},  # noble gas
    'Li': {'Z': 3,  'n': 2, 'l': 0, 'Ze': 1.2792},
    'Be': {'Z': 4,  'n': 2, 'l': 0, 'Ze': 1.9120},
    'B':  {'Z': 5,  'n': 2, 'l': 1, 'Ze': 2.4214},
    'C':  {'Z': 6,  'n': 2, 'l': 1, 'Ze': 3.1358},
    'N':  {'Z': 7,  'n': 2, 'l': 1, 'Ze': 3.8340},
    'O':  {'Z': 8,  'n': 2, 'l': 1, 'Ze': 4.4532},
    'F':  {'Z': 9,  'n': 2, 'l': 1, 'Ze': 5.0998},
    'Ne': {'Z': 10, 'n': 2, 'l': 1, 'Ze': 5.7584},  # noble gas
    'Na': {'Z': 11, 'n': 3, 'l': 0, 'Ze': 2.5074},
    'Mg': {'Z': 12, 'n': 3, 'l': 0, 'Ze': 3.3075},
    'Al': {'Z': 13, 'n': 3, 'l': 1, 'Ze': 3.4959},
    'Si': {'Z': 14, 'n': 3, 'l': 1, 'Ze': 4.1259},
    'P':  {'Z': 15, 'n': 3, 'l': 1, 'Ze': 4.8864},  # wait, same as Cl?
    'S':  {'Z': 16, 'n': 3, 'l': 1, 'Ze': 5.4819},
    'Cl': {'Z': 17, 'n': 3, 'l': 1, 'Ze': 6.1161},
    'Ar': {'Z': 18, 'n': 3, 'l': 1, 'Ze': 6.7641},  # noble gas
}

# Wait, I had Cl Z_eff = 6.1161 before, but that might have been P!
# Let me check: Clementi-Raimondi Z_eff for 3p orbital:
# P (Z=15): 4.8864
# Cl (Z=17): 6.1161
# Need to verify which we've been using in bonds...

print("="*70)
print("  CLEMENTI-RAIMONDI Z_eff VALUES (FULL PERIODS 1-3)")
print("="*70)
for name, at in sorted(atoms.items(), key=lambda x: x[1]['Z']):
    print(f"  {name:3s}: Z={at['Z']:2d}, n={at['n']}, l={at['l']}, Ze={at['Ze']:.4f}")

print("\n" + "="*70)
print("  INCREMENT PATTERNS")
print("="*70)
print("\n  Period 2 (2p): B through Ne")
p2_names = ['B', 'C', 'N', 'O', 'F', 'Ne']
for i in range(1, len(p2_names)):
    a1 = atoms[p2_names[i-1]]
    a2 = atoms[p2_names[i]]
    dZe = a2['Ze'] - a1['Ze']
    print(f"    {p2_names[i-1]}->{p2_names[i]}: dZe = {dZe:.4f}")

print("\n  Period 3 (3p): Al through Ar")
p3_names = ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']
for i in range(1, len(p3_names)):
    a1 = atoms[p3_names[i-1]]
    a2 = atoms[p3_names[i]]
    dZe = a2['Ze'] - a1['Ze']
    print(f"    {p3_names[i-1]}->{p3_names[i]}: dZe = {dZe:.4f}")

print("\n  Period 2 (2s): Li, Be")
print(f"    Li->Be: dZe = {atoms['Be']['Ze'] - atoms['Li']['Ze']:.4f}")

print("\n  Period 3 (3s): Na, Mg")
print(f"    Na->Mg: dZe = {atoms['Mg']['Ze'] - atoms['Na']['Ze']:.4f}")


print("\n" + "="*70)
print("  HOLE SYMMETRY: Distance from noble gas")
print("="*70)

# For each atom, compute distance from nearest noble gas
noble_gases = {
    'He': {'Z': 2,  'n': 1, 'Ze': 1.6875, 'shell_cap': 2},
    'Ne': {'Z': 10, 'n': 2, 'Ze': 5.7584, 'shell_cap': 8},
    'Ar': {'Z': 18, 'n': 3, 'Ze': 6.7641, 'shell_cap': 8},  # 3s+3p cap = 8
}

print("  Atom   Z   d_from_NG   NG   Ze_atom   Ze_NG")
for name, at in sorted(atoms.items(), key=lambda x: x[1]['Z']):
    Z = at['Z']
    # Find nearest noble gas
    if Z <= 2:
        ng_name, ng = 'He', noble_gases['He']
        d_ng = Z - ng['Z']  # negative = before, positive = after
    elif Z <= 10:
        # Could be after He or before Ne
        d_he = Z - 2
        d_ne = Z - 10
        if abs(d_he) <= abs(d_ne):
            ng_name, ng = 'He', noble_gases['He']
            d_ng = d_he
        else:
            ng_name, ng = 'Ne', noble_gases['Ne']
            d_ng = d_ne
    else:
        d_ne = Z - 10
        d_ar = Z - 18
        if abs(d_ne) <= abs(d_ar):
            ng_name, ng = 'Ne', noble_gases['Ne']
            d_ng = d_ne
        else:
            ng_name, ng = 'Ar', noble_gases['Ar']
            d_ng = d_ar

    print(f"  {name:3s}   {Z:2d}   {d_ng:+3d}         {ng_name}   {at['Ze']:.4f}   {ng['Ze']:.4f}")


print("\n" + "="*70)
print("  PARTICLE-HOLE PAIRS")
print("="*70)
# Na (+1 from Ne) and F (-1 from Ne): same distance, opposite sides
# Mg (+2 from Ne) and O (-2 from Ne)
# Li (+1 from He) and H (-1 from He)
# etc.

pairs = [
    ('H', 'He', 'Li', 'He'),   # 1 before/after He (not great, H isn't really -1 from He)
    ('F', 'Ne', 'Na', 'Ne'),   # 1 before/after Ne
    ('O', 'Ne', 'Mg', 'Ne'),   # 2 before/after Ne
    ('N', 'Ne', 'Al', 'Ne'),   # 3 before/after Ne
    ('Cl', 'Ar', 'Na', 'Ne'),  # 1 before Ar vs 1 after Ne
]

print("\n  SYMMETRIC PAIRS around noble gas:")
print("  atom(-d)  Ze      atom(+d)  Ze      Ze_sum   Ze_NG*2")
for neg, ng1, pos, ng2 in pairs:
    ze_sum = atoms[neg]['Ze'] + atoms[pos]['Ze']
    ng_ze = noble_gases[ng1]['Ze']
    if ng1 == ng2:
        print(f"  {neg}(-d)    {atoms[neg]['Ze']:.4f}  {pos}(+d)    {atoms[pos]['Ze']:.4f}  sum={ze_sum:.4f}  2*Ze_{ng1}={2*ng_ze:.4f}")

# What about just within period 3?
# Ar is the complete wave. Cl = Ar - 1, S = Ar - 2, etc.
# Na = Ne + 1, Mg = Ne + 2, etc.
# Is there a symmetry: Ze(Ar-k) + Ze(Ne+k) = const?
print("\n  Period 3 symmetry: Ze(Ne+k) + Ze(Ar-k)")
p3_atoms = [('Na',11), ('Mg',12), ('Al',13), ('Si',14), ('P',15), ('S',16), ('Cl',17)]
for name, Z in p3_atoms:
    k_from_Ne = Z - 10
    k_from_Ar = 18 - Z
    # Find the symmetric partner
    partner_Z = 18 - k_from_Ne + 10  # same distance from Ar as this is from Ne
    # Actually: if this is Ne+k, partner is Ar-k
    partner_Z2 = 18 - k_from_Ne
    for pname, pZ in p3_atoms:
        if pZ == partner_Z2:
            ze_sum = atoms[name]['Ze'] + atoms[pname]['Ze']
            ze_ar = noble_gases['Ar']['Ze']
            print(f"    {name}(Ne+{k_from_Ne}) + {pname}(Ar-{k_from_Ne}) = "
                  f"{atoms[name]['Ze']:.4f} + {atoms[pname]['Ze']:.4f} = {ze_sum:.4f}  "
                  f"(2*Ar_Ze = {2*ze_ar:.4f})")

# Same for period 2
print("\n  Period 2 symmetry: Ze(He+k) + Ze(Ne-k)")
p2_atoms = [('Li',3), ('Be',4), ('B',5), ('C',6), ('N',7), ('O',8), ('F',9)]
for name, Z in p2_atoms:
    k_from_He = Z - 2
    partner_Z = 10 - k_from_He
    for pname, pZ in p2_atoms:
        if pZ == partner_Z:
            ze_sum = atoms[name]['Ze'] + atoms[pname]['Ze']
            ze_ne = noble_gases['Ne']['Ze']
            print(f"    {name}(He+{k_from_He}) + {pname}(Ne-{k_from_He}) = "
                  f"{atoms[name]['Ze']:.4f} + {atoms[pname]['Ze']:.4f} = {ze_sum:.4f}  "
                  f"(2*Ne_Ze = {2*ze_ne:.4f})")


print("\n" + "="*70)
print("  HARMONIC MODEL WITH CORRECT Cl Z_eff")
print("="*70)
# CRITICAL: Check if we had the wrong Z_eff for Cl!
# The value 4.8864 is actually for PHOSPHORUS (Z=15), not Cl (Z=17)!
# Real Cl 3p Z_eff = 6.1161

print("  IMPORTANT: We may have been using the WRONG Z_eff for Cl!")
print(f"  Previous Cl Ze = 4.8864 -- this is actually P (Z=15)")
print(f"  Correct  Cl Ze = 6.1161 -- Clementi-Raimondi for Z=17")
print()

# Rerun harmonic model with correct values
test_atoms = {
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

print("  Harmonic model with CORRECT Cl Z_eff:")
for name, at in test_atoms.items():
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        if ni == nv and li == lv:
            g = g_same
        elif ni == nv:
            g = g_diff
        else:
            g = g_diff
        total_s += (1 - g * r2) * Ne
    pred = Z - total_s
    print(f"  {name:3s}: pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}")


# Now let's check: what Z_eff does the BOND calculation actually use?
print("\n" + "="*70)
print("  CHECK: What Z_eff is in the bond calculation source?")
print("="*70)
# The source of truth is gwt_lagrangian.py
# Let's just print what matters: if Cl is 4.89 or 6.12
print("  Need to check gwt_lagrangian.py for the Cl Z_eff value!")


print("\n" + "="*70)
print("  PERIOD 3 HARMONIC WITH ALL ATOMS")
print("="*70)

# Full period 3 analysis
p3_full = {
    'Na': {'Z': 11, 'n': 3, 'l': 0, 'Ze': 2.5074},
    'Mg': {'Z': 12, 'n': 3, 'l': 0, 'Ze': 3.3075},
    'Al': {'Z': 13, 'n': 3, 'l': 1, 'Ze': 3.4959},
    'Si': {'Z': 14, 'n': 3, 'l': 1, 'Ze': 4.1259},
    'P':  {'Z': 15, 'n': 3, 'l': 1, 'Ze': 4.8864},
    'S':  {'Z': 16, 'n': 3, 'l': 1, 'Ze': 5.4819},
    'Cl': {'Z': 17, 'n': 3, 'l': 1, 'Ze': 6.1161},
    'Ar': {'Z': 18, 'n': 3, 'l': 1, 'Ze': 6.7641},
}

p3_configs = {
    'Na': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 0)],
    'Mg': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 1)],
    'Al': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 0)],
    'Si': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 1)],
    'P':  [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 2)],
    'S':  [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 3)],
    'Cl': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 4)],
    'Ar': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 5)],
}

print("  Harmonic model for full period 3:")
for name in ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']:
    at = p3_full[name]
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in p3_configs[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        if ni == nv and li == lv:
            g = g_same
        elif ni == nv:
            g = g_diff
        else:
            g = g_diff
        total_s += (1 - g * r2) * Ne
    pred = Z - total_s
    print(f"  {name:3s}: Z={Z:2d}, pred={pred:.4f}, real={at['Ze']:.4f}, err={pred-at['Ze']:+.4f}, err%={100*(pred-at['Ze'])/at['Ze']:+.1f}%")

# Now check the INCREMENT pattern: is the harmonic increment correct?
print("\n  Harmonic increment per Z (should match real):")
print("  For 3p: each new Z adds 1 proton and 1 same-shell electron")
print("  dZe/dZ = 1 - s_same = 1 - (1 - g_same) = g_same = 2/3")
print(f"  Predicted increment = {g_same:.4f}")

# Real increment:
for i in range(1, len(p3_names)):
    a1 = p3_full[p3_names[i-1]]
    a2 = p3_full[p3_names[i]]
    print(f"    {p3_names[i-1]}->{p3_names[i]}: real dZe = {a2['Ze']-a1['Ze']:.4f}")

print(f"\n  Average 3p increment: {(p3_full['Ar']['Ze'] - p3_full['Al']['Ze'])/5:.4f}")
print(f"  Expected from g_same = 2/3 = {2/3:.4f}")

# For 2p:
print(f"\n  Average 2p increment: {(atoms['Ne']['Ze'] - atoms['B']['Ze'])/5:.4f}")

# KEY: The increment IS g_same = 2/3! The harmonic model gets the SLOPE right!
# The error is only in the OFFSET (intercept)!

print("\n" + "="*70)
print("  THE ERROR IS A CONSTANT OFFSET!")
print("="*70)
print("  If the increment (slope) is correct at 2/3, the error is in the")
print("  BASE Z_eff of the first atom in each subshell.")
print("  The base is determined by INNER shell screening.")
print()

# For 3p atoms: Ze = Z - inner_screening - (N_same-1)*s_same_subshell - N_diff*s_diff_subshell
# inner_screening = 2*s(1s->3p) + 2*s(2s->3p) + 6*s(2p->3p)
# With harmonic: s = 1 - g_diff*(ni/nv)^2

# What offset correction fixes all 3p atoms?
offsets = []
for name in ['Al', 'Si', 'P', 'S', 'Cl', 'Ar']:
    at = p3_full[name]
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in p3_configs[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        if ni == nv and li == lv:
            g = g_same
        elif ni == nv:
            g = g_diff
        else:
            g = g_diff
        total_s += (1 - g * r2) * Ne
    pred = Z - total_s
    offset = at['Ze'] - pred
    offsets.append(offset)
    print(f"  {name}: pred={pred:.4f}, real={at['Ze']:.4f}, offset={offset:+.4f}")

print(f"\n  Average offset for 3p: {np.mean(offsets):.4f}")
print(f"  Std of offset: {np.std(offsets):.4f}")
print(f"  This offset = extra screening of inner shells needed")

# The offset should be the same for all 3p atoms if the slope is right!
# It represents the FIXED error in inner shell screening.

# What about 2p?
print("\n  2p offsets:")
offsets_2p = []
for name in ['B', 'C', 'N', 'O', 'F', 'Ne']:
    at = atoms[name]
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        if ni == nv and li == lv:
            g = g_same
        elif ni == nv:
            g = g_diff
        else:
            g = g_diff
        total_s += (1 - g * r2) * Ne
    pred = Z - total_s
    offset = at['Ze'] - pred
    offsets_2p.append(offset)
    print(f"  {name}: pred={pred:.4f}, real={at['Ze']:.4f}, offset={offset:+.4f}")

print(f"\n  Average offset for 2p: {np.mean(offsets_2p):.4f}")
print(f"  Std of offset: {np.std(offsets_2p):.4f}")

# 3s offsets:
print("\n  3s offsets:")
for name in ['Na', 'Mg']:
    at = p3_full[name]
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in p3_configs[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        if ni == nv and li == lv:
            g = g_same
        elif ni == nv:
            g = g_diff
        else:
            g = g_diff
        total_s += (1 - g * r2) * Ne
    pred = Z - total_s
    offset = at['Ze'] - pred
    print(f"  {name}: pred={pred:.4f}, real={at['Ze']:.4f}, offset={offset:+.4f}")

# 2s offsets:
print("\n  2s offsets:")
for name in ['Li', 'Be']:
    at = atoms[name]
    Z, nv, lv = at['Z'], at['n'], at['l']
    total_s = 0
    for (ni, li, Ne) in configs[name]:
        if Ne == 0: continue
        r2 = (ni/nv)**2
        if ni == nv and li == lv:
            g = g_same
        elif ni == nv:
            g = g_diff
        else:
            g = g_diff
        total_s += (1 - g * r2) * Ne
    pred = Z - total_s
    offset = at['Ze'] - pred
    print(f"  {name}: pred={pred:.4f}, real={at['Ze']:.4f}, offset={offset:+.4f}")
