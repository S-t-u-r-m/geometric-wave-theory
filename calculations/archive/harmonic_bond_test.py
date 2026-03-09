"""
HARMONIC Z_eff IN BOND FORMULA
================================
Use s = 1 - g * (n_i/n_v)^2 with g_same=2/d, g_diff=4/(2d+1)
to derive Z_eff from first principles, then use in bond formula.

Also: investigate WHY Na/Cl fail and whether filled-shell effects
follow a harmonic pattern too.
"""
import numpy as np

pi = np.pi
E_H = 13.6057
dd = 3
C_bond = pi / dd
f_pi = dd**2 / (dd**2 + 1)
alpha_n = 1 - f_pi / dd
beta_n = (1 + f_pi) / 2
f_anti = 2*dd / (2*dd - 1)
c_ionic = 1.0 / (2*dd + 1)

g_same = 2.0 / dd          # 2/3
g_diff = 4.0 / (2*dd + 1)  # 4/7

# Known Z_eff values (from Clementi-Raimondi, spectroscopic)
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

electron_configs = {
    'H':  [],
    'Li': [(1, 0, 2)],
    'B':  [(1, 0, 2), (2, 0, 2)],
    'C':  [(1, 0, 2), (2, 0, 2), (2, 1, 1)],
    'N':  [(1, 0, 2), (2, 0, 2), (2, 1, 2)],
    'O':  [(1, 0, 2), (2, 0, 2), (2, 1, 3)],
    'F':  [(1, 0, 2), (2, 0, 2), (2, 1, 4)],
    'Na': [(1, 0, 2), (2, 0, 2), (2, 1, 6)],
    'Cl': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 4)],
}

# =============================================================================
# FILLED SHELL ANALYSIS
# =============================================================================
print("=" * 80)
print("  WHY DO Na AND Cl FAIL?")
print("=" * 80)

# The 2p subshell holds 2*(2l+1) = 6 electrons
# When FULL (2p^6 in Na), it forms a complete spherical shell
# A complete shell screens more effectively than individual electrons
# because the density is uniform in angle -> perfect Gauss law screening

# In wave terms: a FILLED harmonic = complete standing wave pattern
# This provides TOTAL screening (s=1) for that shell

# Hypothesis: filled subshell screening is DIFFERENT
# s_filled = 1 (complete screening, like Gauss's law for spherical shell)
# s_partial = 1 - g*(n_i/n_v)^2 (partial screening)

# For Na (3s^1): inner shells are 1s^2 (filled), 2s^2 (filled), 2p^6 (filled)
# If filled subshells screen completely (s=1 each):
# S = 2*1 + 2*1 + 6*1 = 10
# Z_eff = 11 - 10 = 1.0
# But real Z_eff = 2.507, so that's too much screening

# Partial filling effect: a filled shell IS a complete wave,
# but the outer electron still penetrates it on some orbits
# The penetration probability depends on (n_outer/n_inner)^something

# For Na: the 3s electron penetrates into the n=2 core
# Penetration ~ (n_inner/n_outer)^2 (same harmonic ratio!)
# So effective screening by filled n=2 shell:
# s_filled(2->3) = 1 - g_penetrate * (n_inner/n_outer)^2... no wait, opposite

# Actually: penetration REDUCES screening
# s = 1 - penetration_fraction
# penetration ~ g * (n_i/n_v)^p  BUT this is what we already have!

# The problem: we use g_diff for 2p electrons screening 3s
# s(2p->3s) = 1 - 4/7 * (2/3)^2 = 1 - 4/7 * 4/9 = 1 - 16/63 = 0.746
# With 6 electrons: 6 * 0.746 = 4.476
# With 2 electrons (2s): 2 * 0.746 = 1.492
# With 2 electrons (1s): s = 1 - 4/7 * (1/3)^2 = 1 - 4/63 = 0.937, 2*0.937 = 1.873
# Total: 1.873 + 1.492 + 4.476 = 7.841
# Z_eff = 11 - 7.841 = 3.159 (real: 2.507)

# We're under-screening by 0.65. The 2p^6 fills 6 slots but each only screens
# by 0.746 when it should screen more.

# What if FILLED subshells get g_filled instead of g_diff?
# A filled subshell has 2*(2l+1) electrons forming a complete harmonic
# The complete harmonic screens more effectively

# For Na: need s(2p_filled -> 3s) such that:
# 2*s(1s->3s) + 2*s(2s->3s) + 6*s(2p->3s) = 11 - 2.507 = 8.493

# With s(1s->3s) = 1 - g_diff*(1/3)^2 = 59/63 = 0.937
# And s(2s->3s) = 1 - g_diff*(2/3)^2 = 47/63 = 0.746
# We need: 1.873 + 1.492 + 6*s_2p = 8.493
# 6*s_2p = 5.128, s_2p = 0.855

# Compare with the unfilled prediction: s=0.746
# The filled shell gives s=0.855, an INCREASE of 0.109
# Relative increase: 0.855/0.746 = 1.146

# What if filled shells screen by s_filled = 1 - g_filled * (n_i/n_v)^2?
# 0.855 = 1 - g_filled * 4/9
# g_filled = (1-0.855) * 9/4 = 0.145 * 2.25 = 0.326

# Compare: g_diff = 4/7 = 0.571, g_filled = 0.326 ~ 1/d = 0.333!
# Hmm, or it could be that filled shells use g_same = 2/d = 0.667
# 1 - 2/3 * 4/9 = 1 - 8/27 = 19/27 = 0.704 * 6 = 4.222
# Total = 1.873 + 1.492 + 4.222 = 7.587, Z_eff = 3.413 (still too high)

# Let's be more systematic. Check if there's a FILLING FACTOR
print("\nNa analysis (Z=11, 3s^1, Z_eff=2.5074):")
print(f"  Need total screening S = {11 - 2.5074:.4f}")

s_1s_3 = 1 - g_diff * (1/3)**2
s_2_3 = 1 - g_diff * (2/3)**2
print(f"  s(1s->3s) = 1 - 4/7*(1/3)^2 = {s_1s_3:.4f}")
print(f"  s(2x->3s) = 1 - 4/7*(2/3)^2 = {s_2_3:.4f}")
print(f"  With 2(1s) + 2(2s) + 6(2p):")
S_naive = 2*s_1s_3 + 2*s_2_3 + 6*s_2_3
print(f"  S_naive = {S_naive:.4f}, Z_eff_naive = {11-S_naive:.4f}")

S_needed = 11 - 2.5074
s_2p_needed = (S_needed - 2*s_1s_3 - 2*s_2_3) / 6
print(f"  s(2p->3s) needed = {s_2p_needed:.4f}")
g_2p_needed = (1 - s_2p_needed) / (2/3)**2
print(f"  g(2p->3s) needed = {g_2p_needed:.4f}")

# =============================================================================
# OBSERVATION: Maybe the screening depends on ANGULAR MOMENTUM
# =============================================================================
print()
print("=" * 80)
print("  ANGULAR MOMENTUM DEPENDENCE")
print("=" * 80)

# p orbitals have l=1, s orbitals have l=0
# p orbitals have angular momentum that pushes density AWAY from nucleus
# So p electrons screen LESS effectively than s electrons at the same n
# Wait no — p electrons have their density further out, so they screen
# the outer electron MORE (they're between it and the nucleus)

# Actually: s electrons have more density NEAR the nucleus (penetration)
# so they're LESS effective at screening an outer electron
# p electrons have density at intermediate radii, better screening

# For Na: 2s^2 should screen LESS than 2p^6 per electron because
# 2s has more penetration (reaches closer to nucleus, less in the way)

# Let's try: screening depends on l of the screening electron
# s_electrons: g_s = g_diff (less screening = more penetration = larger g)
# p_electrons: g_p = g_less (more screening = less penetration = smaller g)

# From the 2p series, same-shell p->p screening: s = 1/3 (g=2/3)
# From B, 2s->2p screening: s = 3/7 (g = 4/7)
# So 2s screens 2p MORE than 2p screens 2p? That seems backwards.
# Wait: s = 3/7 > s = 1/3, so yes, 2s screens MORE.
# This is because 2s has more penetration toward the nucleus
# but also more density between the nucleus and the 2p electron's main lobe

# Hmm, this is getting complicated. Let me instead try a data-driven approach
# for the key question: can we get ALL Z_eff from a formula with d=3 constants?

# =============================================================================
# DATA-DRIVEN: What screening per (n_i, l_i) -> (n_v, l_v) pair?
# =============================================================================
print()
print("=" * 80)
print("  DATA-DRIVEN: Screening per orbital pair")
print("=" * 80)

# We have 9 atoms but many screening pairs
# Let's solve for what we can

# KNOWN (from period 2):
# s(1s->2s) = 0.860  [from Li]
# s(1s->2p) = assume same = 0.860  [or solve from B]
# s(2s->2p) = 0.429  [from B, assuming s(1s->2p)=s(1s->2s)]
# s(2p->2p) = 0.330  [from C-F average]

# FOR PERIOD 3, we need:
# s(1s->3s), s(2s->3s), s(2p->3s)  [from Na]
# s(1s->3p), s(2s->3p), s(2p->3p), s(3s->3p)  [from Cl, 4 unknowns]

# Model: s(n_i, l_i -> n_v, l_v) = 1 - g(l_match) * (n_i/n_v)^2
# With three g values:
# g1: same subshell (same n, same l)
# g2: same n, different l (or adjacent, same l)
# g3: different n entirely

# Actually, let me try the simplest thing: ONE g value for everything
# s = 1 - g * (n_i/n_v)^2
# This gives the same screening to s and p electrons at the same n

print("\nTesting single-g model: s = 1 - g * (n_i/n_v)^2")
best_g1 = (999, 0)
for g in np.arange(0.3, 0.9, 0.001):
    errs = []
    for aname in atoms:
        info = atoms[aname]
        Z = info['Z']; n_v = info['n']; l_v = info['l']
        # Build electron config
        configs = {
            'H': [],
            'Li': [(1, 0, 2)],
            'B': [(1, 0, 2), (2, 0, 2)],
            'C': [(1, 0, 2), (2, 0, 2), (2, 1, 1)],
            'N': [(1, 0, 2), (2, 0, 2), (2, 1, 2)],
            'O': [(1, 0, 2), (2, 0, 2), (2, 1, 3)],
            'F': [(1, 0, 2), (2, 0, 2), (2, 1, 4)],
            'Na': [(1, 0, 2), (2, 0, 2), (2, 1, 6)],
            'Cl': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 4)],
        }
        S = sum(cnt * (1 - g * (n_i/n_v)**2) for n_i, l_i, cnt in configs[aname])
        Ze_pred = Z - S
        errs.append((Ze_pred - info['Z_eff'])**2)
    rms = np.sqrt(np.mean(errs))
    if rms < best_g1[0]:
        best_g1 = (rms, g)

print(f"  Best single g = {best_g1[1]:.4f}, RMS = {best_g1[0]:.4f}")
print(f"  Compare: 4/(2d+1) = {4/(2*dd+1):.4f}")
print(f"  Compare: 1/(dd-1) = {1/(dd-1):.4f}")
print(f"  Compare: (d-1)/d = {(dd-1)/dd:.4f}")

g_single = best_g1[1]
configs_all = {
    'H': [],
    'Li': [(1, 0, 2)],
    'B': [(1, 0, 2), (2, 0, 2)],
    'C': [(1, 0, 2), (2, 0, 2), (2, 1, 1)],
    'N': [(1, 0, 2), (2, 0, 2), (2, 1, 2)],
    'O': [(1, 0, 2), (2, 0, 2), (2, 1, 3)],
    'F': [(1, 0, 2), (2, 0, 2), (2, 1, 4)],
    'Na': [(1, 0, 2), (2, 0, 2), (2, 1, 6)],
    'Cl': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 4)],
}

print(f"\n  {'Atom':>4} {'Z':>3} {'Z_pred':>7} {'Z_real':>7} {'err':>7} {'E_pred':>8} {'E_real':>8}")
print("  " + "-" * 55)
for name in atoms:
    info = atoms[name]
    Z = info['Z']; n_v = info['n']
    S = sum(cnt * (1 - g_single * (n_i/n_v)**2) for n_i, l_i, cnt in configs_all[name])
    Ze_pred = Z - S
    Ze_real = info['Z_eff']
    E_pred = E_H * (Ze_pred/n_v)**2
    E_real = E_H * (Ze_real/n_v)**2
    print(f"  {name:>4} {Z:3d} {Ze_pred:7.4f} {Ze_real:7.4f} {Ze_pred-Ze_real:+7.4f} "
          f"{E_pred:8.4f} {E_real:8.4f}")


# =============================================================================
# KEY INSIGHT: n_i^2/n_v^2 IS the energy ratio!
# =============================================================================
print()
print("=" * 80)
print("  KEY INSIGHT: (n_i/n_v)^2 = E_outer/E_inner (energy ratio)")
print("=" * 80)
print()
print("  For hydrogen-like atoms:")
print("    E_n = -E_H/n^2")
print("    E_outer/E_inner = (n_i/n_v)^2")
print()
print("  So screening = 1 - g * (E_outer/E_inner)")
print("  Electrons whose energy is CLOSE to the valence energy")
print("  screen LESS (they're at similar radii, more 'resonant')")
print("  Electrons whose energy is FAR from valence (deep core)")
print("  screen MORE (they're definitely inside)")
print()
print("  This IS the harmonic pattern: resonance reduces screening!")


# =============================================================================
# BOND PREDICTIONS with harmonic Z_eff
# =============================================================================
print()
print("=" * 80)
print("  BOND PREDICTIONS: harmonic Z_eff")
print("=" * 80)

Z_eff_real = {
    'H': 1.0, 'Li': 1.2792, 'B': 2.4214, 'C': 3.1358, 'N': 3.8340,
    'O': 4.4532, 'F': 5.0998, 'Na': 2.5074, 'Cl': 6.1161
}

def harmonic_zeff_2g(name):
    """Z_eff with g_same=2/d, g_diff=4/(2d+1), p=2"""
    info = atoms[name]
    Z = info['Z']; n_v = info['n']; l_v = info['l']
    S = 0
    for n_i, l_i, cnt in configs_all[name]:
        if n_i == n_v and l_i == l_v:
            g = g_same  # 2/3
        else:
            g = g_diff  # 4/7
        s = 1 - g * (n_i/n_v)**2
        S += cnt * s
    return Z - S

def harmonic_zeff_1g(name):
    """Z_eff with single g"""
    info = atoms[name]
    Z = info['Z']; n_v = info['n']
    S = sum(cnt * (1 - g_single * (n_i/n_v)**2) for n_i, l_i, cnt in configs_all[name])
    return Z - S

molecules = [
    ('H2',   1.401,  4.745, [('ss', 1)], 'H',  'H'),
    ('Li2',  5.051,  1.056, [('ss', 1)], 'Li', 'Li'),
    ('B2',   3.005,  3.02,  [('pi', 2)], 'B',  'B'),
    ('C2',   2.348,  6.32,  [('pi', 2), ('sp_sigma', 1), ('sp_sigma_anti', 1)], 'C',  'C'),
    ('N2',   2.074,  9.759, [('pp_sigma', 1), ('pi', 2)], 'N',  'N'),
    ('O2',   2.282,  5.213, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)], 'O',  'O'),
    ('F2',   2.668,  1.660, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)], 'F',  'F'),
    ('Na2',  5.818,  0.746, [('ss', 1)], 'Na', 'Na'),
    ('Cl2',  3.757,  2.514, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)], 'Cl', 'Cl'),
    ('HF',   1.733,  5.869, [('sp', 1)], 'H',  'F'),
    ('CO',   2.132, 11.225, [('pp_sigma', 1), ('pi', 2)], 'C',  'O'),
    ('NO',   2.175,  6.497, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)], 'N',  'O'),
    ('OH',   1.834,  4.392, [('sp', 1)], 'O',  'H'),
    ('HCl',  2.409,  4.434, [('sp', 1)], 'H',  'Cl'),
    ('LiH',  3.015,  2.515, [('ss', 1)], 'Li', 'H'),
    ('LiF',  2.955,  5.939, [('sp', 1)], 'Li', 'F'),
    ('BH',   2.329,  3.42,  [('sp', 1)], 'B',  'H'),
    ('CH',   2.116,  3.47,  [('sp', 1)], 'C',  'H'),
    ('NH',   1.958,  3.57,  [('sp', 1)], 'N',  'H'),
    ('BF',   2.386,  7.81,  [('pp_sigma', 1), ('pi', 2)], 'B',  'F'),
    ('CN',   2.214,  7.72,  [('pp_sigma', 1), ('pi', 2)], 'C',  'N'),
    ('NaH',  3.566,  1.97,  [('ss', 1)], 'Na', 'H'),
    ('NaCl', 4.461,  4.23,  [('sp', 1)], 'Na', 'Cl'),
    ('H2O',  1.809,  5.117, [('sp', 1)], 'O',  'H'),
]


def compute_bond(mol, zeff_func):
    name, R, De_exp, bonds, atom1, atom2 = mol
    info1 = atoms[atom1]; info2 = atoms[atom2]
    n1 = info1['n']; l1 = info1['l']; n2 = info2['n']; l2 = info2['l']
    h1 = min(n1 - l1 - 1, 1); h2 = min(n2 - l2 - 1, 1)

    a1 = 2 + (1 - 2*l1) * alpha_n * h1
    a2 = 2 + (1 - 2*l2) * alpha_n * h2
    b1 = 1 + beta_n * h1
    b2 = 1 + beta_n * h2

    E_scale = np.sqrt(E_H / n1**a1 * E_H / n2**a2)
    k1 = 1.0 / n1**b1; k2 = 1.0 / n2**b2
    sigma_phase = R * (k1 + k2)

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss', 'sp')) else sigma_phase * f_pi
        cont = C_bond * E_scale * abs(np.sin(ph))
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    # Ionic with derived Z_eff
    Ze1 = zeff_func(atom1); Ze2 = zeff_func(atom2)
    E1 = E_H * (Ze1/n1)**2; E2 = E_H * (Ze2/n2)**2
    dE = abs(E1 - E2)
    V = max(abs(D_cov), 0.01)
    q = dE / np.sqrt(dE**2 + (2*V)**2) if dE > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / R
    return D_cov + Di


print()
print("  Using 2-g model (g_same=2/d, g_diff=4/(2d+1)):")
print(f"  {'Mol':<6} {'De_exp':>7} {'D_harm':>7} {'D_real':>7} {'err_h%':>7} {'err_r%':>7}")
print("  " + "-" * 50)

errs_h2 = []; errs_r = []
for mol in molecules:
    D_h = compute_bond(mol, harmonic_zeff_2g)
    D_r = compute_bond(mol, lambda n: Z_eff_real[n])
    err_h = (D_h - mol[2]) / mol[2] * 100
    err_r = (D_r - mol[2]) / mol[2] * 100
    errs_h2.append(abs(err_h)); errs_r.append(abs(err_r))
    print(f"  {mol[0]:<6} {mol[2]:7.3f} {D_h:7.3f} {D_r:7.3f} {err_h:+6.1f}% {err_r:+6.1f}%")

print(f"\n  Harmonic 2g: avg={np.mean(errs_h2):.1f}%, med={np.median(errs_h2):.1f}%")
print(f"  Real Z_eff:  avg={np.mean(errs_r):.1f}%, med={np.median(errs_r):.1f}%")


# =============================================================================
# ENERGY DIFFERENCE COMPARISON
# =============================================================================
print()
print("=" * 80)
print("  ENERGY DIFFERENCES: Harmonic vs Real Z_eff")
print("=" * 80)
print(f"\n  {'Pair':>8} {'dE_harm':>8} {'dE_real':>8} {'dE_form':>8} {'ratio':>7}")
print("  " + "-" * 50)

# For each heteronuclear pair, show energy differences
pairs = set()
for mol in molecules:
    a1, a2 = mol[4], mol[5]
    if a1 != a2:
        pair = tuple(sorted([a1, a2]))
        if pair not in pairs:
            pairs.add(pair)
            info1 = atoms[pair[0]]; info2 = atoms[pair[1]]
            n1 = info1['n']; n2 = info2['n']

            Ze1_h = harmonic_zeff_2g(pair[0]); Ze2_h = harmonic_zeff_2g(pair[1])
            Ze1_r = Z_eff_real[pair[0]]; Ze2_r = Z_eff_real[pair[1]]

            E1_h = E_H * (Ze1_h/n1)**2; E2_h = E_H * (Ze2_h/n2)**2
            E1_r = E_H * (Ze1_r/n1)**2; E2_r = E_H * (Ze2_r/n2)**2

            l1 = info1['l']; l2 = info2['l']
            h1 = min(n1-l1-1, 1); h2 = min(n2-l2-1, 1)
            a1 = 2 + (1-2*l1)*alpha_n*h1
            a2 = 2 + (1-2*l2)*alpha_n*h2
            E1_f = E_H/n1**a1; E2_f = E_H/n2**a2

            dE_h = abs(E1_h - E2_h)
            dE_r = abs(E1_r - E2_r)
            dE_f = abs(E1_f - E2_f)

            ratio = dE_h / dE_r if dE_r > 0.01 else float('inf')
            print(f"  {pair[0]+'-'+pair[1]:>8} {dE_h:8.3f} {dE_r:8.3f} {dE_f:8.3f} {ratio:7.3f}")


# =============================================================================
# WHAT IF: formula energy + harmonic correction?
# =============================================================================
print()
print("=" * 80)
print("  HYBRID: Formula energy + harmonic Z_eff for ionic only")
print("=" * 80)
print()
print("  Covalent uses E_H/n^a (no Z_eff)")
print("  Ionic uses harmonic Z_eff for dE")
print("  -> This keeps covalent self-consistent while giving")
print("     atom-specific ionic corrections")

# The ionic term only needs the ENERGY DIFFERENCE
# Use harmonic Z_eff just for that
print(f"\n  Harmonic Z_eff values:")
for name in atoms:
    Ze_h = harmonic_zeff_2g(name)
    Ze_r = Z_eff_real[name]
    n = atoms[name]['n']
    E_h = E_H * (Ze_h/n)**2
    E_r = E_H * (Ze_r/n)**2
    print(f"    {name:>4}: Z_h={Ze_h:.4f} (real={Ze_r:.4f}), E_h={E_h:.3f} (real={E_r:.3f})")
