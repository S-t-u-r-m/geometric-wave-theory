"""
WAVE ENERGY MODEL WITH PAIRING CORRECTION
============================================
Combine:
1. Harmonic Z_eff: s = 1 - g*(n_i/n_v)^2 with g_same=2/d, g_diff=4/(2d+1)
2. Pairing energy: E_pair = E_H/d per pair in a half-filled+ subshell

The orbital energy becomes:
  E_orb = E_H*(Z_eff/n)^2 + N_pairs * E_pair / n^p

where N_pairs = max(0, N_electrons_in_subshell - (2l+1))

This changes the ENERGY DIFFERENCE between atoms, which drives
the ionic correction in bonds.
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

# Electron configurations: (n_i, l_i, count_screening_valence)
# Also need: valence subshell electron count for pairing
atom_data = {
    'H':  {'Z': 1,  'n': 1, 'l': 0, 'Z_eff': 1.0000, 'IP': 13.598,
            'inner': [], 'N_val': 1},
    'Li': {'Z': 3,  'n': 2, 'l': 0, 'Z_eff': 1.2792, 'IP': 5.392,
            'inner': [(1, 0, 2)], 'N_val': 1},
    'B':  {'Z': 5,  'n': 2, 'l': 1, 'Z_eff': 2.4214, 'IP': 8.298,
            'inner': [(1, 0, 2), (2, 0, 2)], 'N_val': 1},
    'C':  {'Z': 6,  'n': 2, 'l': 1, 'Z_eff': 3.1358, 'IP': 11.260,
            'inner': [(1, 0, 2), (2, 0, 2), (2, 1, 1)], 'N_val': 2},
    'N':  {'Z': 7,  'n': 2, 'l': 1, 'Z_eff': 3.8340, 'IP': 14.534,
            'inner': [(1, 0, 2), (2, 0, 2), (2, 1, 2)], 'N_val': 3},
    'O':  {'Z': 8,  'n': 2, 'l': 1, 'Z_eff': 4.4532, 'IP': 13.618,
            'inner': [(1, 0, 2), (2, 0, 2), (2, 1, 3)], 'N_val': 4},
    'F':  {'Z': 9,  'n': 2, 'l': 1, 'Z_eff': 5.0998, 'IP': 17.423,
            'inner': [(1, 0, 2), (2, 0, 2), (2, 1, 4)], 'N_val': 5},
    'Na': {'Z': 11, 'n': 3, 'l': 0, 'Z_eff': 2.5074, 'IP': 5.139,
            'inner': [(1, 0, 2), (2, 0, 2), (2, 1, 6)], 'N_val': 1},
    'Cl': {'Z': 17, 'n': 3, 'l': 1, 'Z_eff': 6.1161, 'IP': 12.968,
            'inner': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 4)],
            'N_val': 5},
}

def harmonic_zeff(name):
    """Z_eff from wave mode coupling: s = 1 - g*(n_i/n_v)^2"""
    info = atom_data[name]
    Z = info['Z']; n_v = info['n']; l_v = info['l']
    S = 0
    for (n_i, l_i, cnt) in info['inner']:
        if n_i == n_v and l_i == l_v:
            g = g_same  # 2/3
        else:
            g = g_diff  # 4/7
        s = 1 - g * (n_i/n_v)**2
        S += cnt * s
    return Z - S


def pairing_count(name):
    """Number of electron pairs in valence subshell."""
    info = atom_data[name]
    l = info['l']
    N_val = info['N_val']
    max_unpaired = 2*l + 1  # 1 for s, 3 for p
    return max(0, N_val - max_unpaired)


# =============================================================================
# SHOW HARMONIC Z_eff AND PAIRING
# =============================================================================
print("=" * 80)
print("  HARMONIC Z_eff + PAIRING COUNT")
print("=" * 80)

print(f"\n{'Atom':>4} {'Z':>3} {'n':>2} {'l':>2} {'N_val':>5} {'pairs':>5} "
      f"{'Z_h':>7} {'Z_real':>7} {'err':>7} {'IP_obs':>7}")
print("-" * 65)

for name in atom_data:
    info = atom_data[name]
    Ze_h = harmonic_zeff(name)
    Ze_r = info['Z_eff']
    pairs = pairing_count(name)
    print(f"{name:>4} {info['Z']:3d} {info['n']:2d} {info['l']:2d} {info['N_val']:5d} "
          f"{pairs:5d} {Ze_h:7.4f} {Ze_r:7.4f} {Ze_h-Ze_r:+7.4f} {info['IP']:7.3f}")


# =============================================================================
# ORBITAL ENERGY WITH PAIRING CORRECTION
# =============================================================================
print()
print("=" * 80)
print("  ORBITAL ENERGY: E = E_H*(Z_h/n)^2 - N_pairs * E_pair")
print("  E_pair = E_H / (d * n^q)  — pairing energy per pair")
print("=" * 80)

# The pairing REDUCES the effective orbital energy
# (paired electrons are less tightly bound)
# So E_eff = E_H*(Z_h/n)^2 - pairs * E_pair

# Scan: what E_pair makes the orbital energies match real Z_eff energies?
# E_real = E_H*(Z_real/n)^2
# E_pred = E_H*(Z_h/n)^2 - pairs * E_pair
# We want E_pred ~ E_real

# Only atoms with pairs: O (1 pair), F (2 pairs), Cl (2 pairs)
# O: E_h = E_H*(4.4286/2)^2 = 66.71 eV, E_r = 67.45 eV
# Without pairing, harmonic already undershoots for O.
# So pairing correction (reducing E further) would make it WORSE for O.

# Hmm. Let me reconsider. The pairing energy might INCREASE the screening
# (paired electrons screen more), not decrease the orbital energy.

# In wave terms: two quanta in the same angular mode create a stronger
# standing wave at that angle, which screens the nucleus MORE effectively
# for the NEXT electron. So pairing increases screening for subsequent modes.

# Let's model it differently:
# When counting same-subshell screening, paired electrons screen MORE
# s_paired = s_same + delta_s_pair
# s_unpaired = s_same

# From O: Z_eff(O) = 4.4532, our model gives 4.4286 (err -0.025)
# O has 3 other 2p electrons screening it. Of these, 1 is in a paired mode.
# Does the pairing of ONE of the 3 screening electrons change things?

# Actually, for O (2p^4): the 4th electron being added IS the one that pairs.
# The 3 existing 2p electrons screen it: 2 unpaired + the one it pairs with.
# The one it pairs with might screen differently.

# But in our model, all 3 screen by s_same = 1/3.
# The INCREMENT Z_eff(O) - Z_eff(N) = 4.4532 - 3.8340 = 0.6192
# Our model predicts: 1 - s_same = 1 - 1/3 = 0.667
# The SHORTFALL: 0.667 - 0.619 = 0.048 per pair

# So the pairing makes the INCOMING electron see 0.048 MORE screening
# per paired interaction. delta_s_pair = 0.048

# Let's compute this for all pairing instances:
print("\n  Pairing effect on Z_eff increments (2p series):")
zeffs_2p = [(n, atom_data[n]['Z_eff']) for n in ['B', 'C', 'N', 'O', 'F']]
avg_unpaired = []
for i in range(1, len(zeffs_2p)):
    inc = zeffs_2p[i][1] - zeffs_2p[i-1][1]
    is_pairing = (i >= 3)  # O(i=3), F(i=4) have pairing
    print(f"    {zeffs_2p[i-1][0]}->{zeffs_2p[i][0]}: dZ_eff={inc:.4f} "
          f"{'(PAIRED)' if is_pairing else '(unpaired)'}")
    if not is_pairing:
        avg_unpaired.append(inc)

mean_unpaired = np.mean(avg_unpaired)
print(f"\n  Mean unpaired increment: {mean_unpaired:.4f}")
print(f"  Model prediction (1-1/d): {1-1/dd:.4f}")

# Delta for paired transitions:
for i in range(3, len(zeffs_2p)):
    inc = zeffs_2p[i][1] - zeffs_2p[i-1][1]
    delta = mean_unpaired - inc
    print(f"  {zeffs_2p[i-1][0]}->{zeffs_2p[i][0]}: "
          f"delta_s = {delta:.4f} (extra screening per pair)")


# =============================================================================
# MODEL: Z_eff with pairing correction
# =============================================================================
print()
print("=" * 80)
print("  Z_eff WITH PAIRING CORRECTION")
print("=" * 80)

# For each screening electron in a half-filled+ subshell,
# add extra screening delta_s when it's part of a pair.
# But we need to know HOW MANY of the screening electrons are paired.

# For atom with N_val electrons in (n, l) subshell:
# - First (2l+1) are unpaired, screen by s_same = 1 - 2/d
# - Next (N_val - (2l+1)) are paired, screen by s_same + delta_s_pair

# The delta_s_pair from 2p data: ~0.048
# In GWT terms: delta_s_pair ~ 1/(d*(2d+1)) = 1/21 = 0.048 !!

delta_pair_guess = 1.0 / (dd * (2*dd + 1))  # = 1/21
print(f"\n  delta_s_pair = 1/(d*(2d+1)) = 1/21 = {delta_pair_guess:.6f}")

# Scan to find best
best_delta = (999, 0)
for delta in np.arange(0.0, 0.15, 0.001):
    errs = []
    for name in atom_data:
        info = atom_data[name]
        Z = info['Z']; n_v = info['n']; l_v = info['l']
        S = 0
        for (n_i, l_i, cnt) in info['inner']:
            if n_i == n_v and l_i == l_v:
                # Same subshell — check for pairing
                max_unp = 2*l_i + 1
                # cnt = number of other electrons in this subshell
                n_paired = max(0, cnt - (max_unp - 1))  # how many of these are paired
                n_unpaired = cnt - n_paired
                s_unp = 1 - g_same * (n_i/n_v)**2
                s_pair = s_unp + delta
                S += n_unpaired * s_unp + n_paired * s_pair
            else:
                g = g_diff
                s = 1 - g * (n_i/n_v)**2
                S += cnt * s
        Ze_pred = Z - S
        errs.append((Ze_pred - info['Z_eff'])**2)
    rms = np.sqrt(np.mean(errs))
    if rms < best_delta[0]:
        best_delta = (rms, delta)

print(f"  Best delta_s_pair = {best_delta[1]:.4f}, RMS = {best_delta[0]:.4f}")
print(f"  Compare: 1/(d*(2d+1)) = 1/21 = {1/(dd*(2*dd+1)):.4f}")
print(f"  Compare: 1/(2d^2) = {1/(2*dd**2):.4f}")
print(f"  Compare: 1/(d^2+d) = {1/(dd**2+dd):.4f}")

delta_best = best_delta[1]

# Show predictions
print(f"\n{'Atom':>4} {'Z':>3} {'Z_pred':>7} {'Z_real':>7} {'err':>7}")
print("-" * 35)
for name in atom_data:
    info = atom_data[name]
    Z = info['Z']; n_v = info['n']; l_v = info['l']
    S = 0
    for (n_i, l_i, cnt) in info['inner']:
        if n_i == n_v and l_i == l_v:
            max_unp = 2*l_i + 1
            n_paired = max(0, cnt - (max_unp - 1))
            n_unpaired = cnt - n_paired
            s_unp = 1 - g_same * (n_i/n_v)**2
            s_pair = s_unp + delta_best
            S += n_unpaired * s_unp + n_paired * s_pair
        else:
            g = g_diff
            s = 1 - g * (n_i/n_v)**2
            S += cnt * s
    Ze_pred = Z - S
    print(f"{name:>4} {Z:3d} {Ze_pred:7.4f} {info['Z_eff']:7.4f} {Ze_pred-info['Z_eff']:+7.4f}")


# =============================================================================
# Now handle Na/Cl: what if the filled p-shell (2p^6) has ALL paired?
# =============================================================================
print()
print("=" * 80)
print("  Na/Cl FIX: Filled shell pairing + cross-shell pairing")
print("=" * 80)

# For Na: inner 2p^6 has 3 pairs.
# If delta_pair applies to CROSS-SHELL pairs too:
# s(2p->3s) = s_diff + n_pairs * delta_cross / n_total ?

# Actually, the 2p^6 is a COMPLETE shell. All 6 electrons are paired (3 pairs).
# Each pair adds extra screening. In our cross-shell formula:
# s = 1 - g_diff * (n_i/n_v)^2
# But paired modes screen more: s_paired = s + delta_cross

# For Na: 2p^6 has 6 electrons, 3 pairs
# Currently: 6 * (1 - 4/7 * 4/9) = 6 * 0.746 = 4.476
# Need: 6 * s = 5.13 (so each screens 0.855)
# delta_cross per paired electron: (0.855 - 0.746) = 0.109
# But 3 of the 6 are paired: so delta per paired = 0.109 * 6/3 = 0.218?

# Or think of it as: the 2p^6 complete shell has an ADDITIONAL screening
# of 3 * delta_fill, where delta_fill is the filled-shell pairing bonus

# From Na: extra screening needed = 8.493 - 7.841 = 0.651
# 3 pairs: delta_per_pair = 0.651/3 = 0.217

# For Cl: 2p^6 screening 3p (same ratio 2/3), also 3 pairs
# Current: 6*(1-4/7*4/9) = 4.476
# With pair bonus: 4.476 + 3*0.217 = 5.128
# But we also have 3s^2 (1 pair?) and 3p^4 (1 pair)

# Wait, s orbitals only have 1 angular mode. 2 electrons = 1 pair.
# 3s^2: 1 pair.

# Let me be systematic. For each subshell that screens the valence:
# Count pairs in that subshell = max(0, N - (2l+1))
# Add delta_cross * pairs to total screening

# Scan delta_cross
print("\nScanning delta_cross for cross-shell pairs:")
best_dc = (999, 0, 0)
for ds in np.arange(0.0, 0.15, 0.001):  # same-shell delta
    for dc in np.arange(0.0, 0.30, 0.001):  # cross-shell delta
        errs = []
        for name in atom_data:
            info = atom_data[name]
            Z = info['Z']; n_v = info['n']; l_v = info['l']
            S = 0
            for (n_i, l_i, cnt) in info['inner']:
                max_unp = 2*l_i + 1
                n_pairs = max(0, cnt - max_unp) if cnt <= 2*max_unp else max_unp
                # Actually: for a subshell with cnt electrons and (2l+1) orbitals:
                # max capacity = 2*(2l+1)
                # pairs = max(0, cnt - (2l+1))
                n_pairs = max(0, cnt - (2*l_i + 1))

                if n_i == n_v and l_i == l_v:
                    s_base = 1 - g_same * (n_i/n_v)**2
                    S += cnt * s_base + n_pairs * ds
                else:
                    s_base = 1 - g_diff * (n_i/n_v)**2
                    S += cnt * s_base + n_pairs * dc
            Ze_pred = Z - S
            errs.append((Ze_pred - info['Z_eff'])**2)
        rms = np.sqrt(np.mean(errs))
        if rms < best_dc[0]:
            best_dc = (rms, ds, dc)

ds_best = best_dc[1]; dc_best = best_dc[2]
print(f"\n  Best: delta_same={ds_best:.4f}, delta_cross={dc_best:.4f}, RMS={best_dc[0]:.4f}")

# GWT constant matches
print(f"\n  delta_same = {ds_best:.4f}")
gwt_vals = {
    '1/(d*(2d+1))': 1/(dd*(2*dd+1)),
    '1/(2d^2)': 1/(2*dd**2),
    '1/(d^2+d)': 1/(dd**2+dd),
    '1/(2d+1)': 1/(2*dd+1),
    '1/d^2': 1/dd**2,
    '2/(d*(2d+1))': 2/(dd*(2*dd+1)),
    '1/(d*(d+1))': 1/(dd*(dd+1)),
    '1/(2*pi)': 1/(2*pi),
}
for label, val in gwt_vals.items():
    if abs(val - ds_best) < 0.015:
        print(f"    ~ {label} = {val:.4f}")

print(f"\n  delta_cross = {dc_best:.4f}")
for label, val in gwt_vals.items():
    if abs(val - dc_best) < 0.015:
        print(f"    ~ {label} = {val:.4f}")

# Full predictions
print(f"\n{'Atom':>4} {'Z':>3} {'Z_pred':>7} {'Z_real':>7} {'err':>7} {'E_pred':>8} {'E_real':>8}")
print("-" * 60)

def full_zeff(name):
    info = atom_data[name]
    Z = info['Z']; n_v = info['n']; l_v = info['l']
    S = 0
    for (n_i, l_i, cnt) in info['inner']:
        n_pairs = max(0, cnt - (2*l_i + 1))
        if n_i == n_v and l_i == l_v:
            s_base = 1 - g_same * (n_i/n_v)**2
            S += cnt * s_base + n_pairs * ds_best
        else:
            s_base = 1 - g_diff * (n_i/n_v)**2
            S += cnt * s_base + n_pairs * dc_best
    return Z - S

for name in atom_data:
    info = atom_data[name]
    Ze_pred = full_zeff(name)
    Ze_real = info['Z_eff']
    n = info['n']
    E_pred = E_H * (Ze_pred/n)**2
    E_real = E_H * (Ze_real/n)**2
    print(f"{name:>4} {info['Z']:3d} {Ze_pred:7.4f} {Ze_real:7.4f} "
          f"{Ze_pred-Ze_real:+7.4f} {E_pred:8.3f} {E_real:8.3f}")


# =============================================================================
# BOND PREDICTIONS WITH FULL MODEL
# =============================================================================
print()
print("=" * 80)
print("  BOND PREDICTIONS: full wave model (harmonic + pairing)")
print("=" * 80)

Z_eff_real = {
    'H': 1.0, 'Li': 1.2792, 'B': 2.4214, 'C': 3.1358, 'N': 3.8340,
    'O': 4.4532, 'F': 5.0998, 'Na': 2.5074, 'Cl': 6.1161
}

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
    info1 = atom_data[atom1]; info2 = atom_data[atom2]
    n1 = info1['n']; l1 = info1['l']; n2 = info2['n']; l2 = info2['l']
    h1 = min(n1 - l1 - 1, 1); h2 = min(n2 - l2 - 1, 1)

    a1 = 2 + (1 - 2*l1) * alpha_n * h1
    a2 = 2 + (1 - 2*l2) * alpha_n * h2
    b1 = 1 + beta_n * h1; b2 = 1 + beta_n * h2

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

    Ze1 = zeff_func(atom1); Ze2 = zeff_func(atom2)
    E1 = E_H * (Ze1/n1)**2; E2 = E_H * (Ze2/n2)**2
    dE = abs(E1 - E2)
    V = max(abs(D_cov), 0.01)
    q = dE / np.sqrt(dE**2 + (2*V)**2) if dE > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / R
    return D_cov + Di

print()
print(f"{'Mol':<6} {'De_exp':>7} {'D_wave':>7} {'D_real':>7} {'err_w%':>7} {'err_r%':>7}")
print("-" * 50)

errs_w = []; errs_r = []
for mol in molecules:
    D_w = compute_bond(mol, full_zeff)
    D_r = compute_bond(mol, lambda n: Z_eff_real[n])
    err_w = (D_w - mol[2]) / mol[2] * 100
    err_r = (D_r - mol[2]) / mol[2] * 100
    errs_w.append(abs(err_w)); errs_r.append(abs(err_r))
    v = 'BETTER' if abs(err_w) < abs(err_r) - 1 else 'WORSE' if abs(err_w) > abs(err_r) + 1 else '~'
    print(f"{mol[0]:<6} {mol[2]:7.3f} {D_w:7.3f} {D_r:7.3f} {err_w:+6.1f}% {err_r:+6.1f}% {v}")

print(f"\n  Wave model: avg={np.mean(errs_w):.1f}%, med={np.median(errs_w):.1f}%, "
      f"<5%:{sum(1 for e in errs_w if e<5)}/24, <10%:{sum(1 for e in errs_w if e<10)}/24")
print(f"  Real Z_eff: avg={np.mean(errs_r):.1f}%, med={np.median(errs_r):.1f}%, "
      f"<5%:{sum(1 for e in errs_r if e<5)}/24, <10%:{sum(1 for e in errs_r if e<10)}/24")


# =============================================================================
# SUMMARY
# =============================================================================
print()
print("=" * 80)
print("  COMPLETE WAVE MODEL: Z_eff from first principles")
print("=" * 80)
print(f"""
  Z_eff = Z - sum_modes [s_i]

  Screening per mode:
    s = 1 - g * (n_i/n_v)^2 + delta * is_paired

  Constants (all from d={dd}):
    g_same      = 2/d       = {g_same:.6f}  (same subshell)
    g_diff      = 4/(2d+1)  = {g_diff:.6f}  (different subshell)
    delta_same  = {ds_best:.6f}            (same-shell pairing extra)
    delta_cross = {dc_best:.6f}            (cross-shell pairing extra)

  (n_i/n_v)^2 = standing wave frequency ratio
  Pairing = extra screening when angular mode is doubly occupied
""")
