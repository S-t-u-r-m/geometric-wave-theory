"""
WAVE MODE COUPLING: Z_eff from standing wave interference
==========================================================
In GWT, atoms are standing waves with Z energy quanta distributed
among wave modes. There are no "electrons screening electrons."

Instead: Z energy modes couple to each other through their
wave overlap integrals. The outermost mode's effective energy
depends on how it couples to all inner modes.

KEY FINDING (from harmonic_zeff.py):
  coupling(n_i -> n_v) = g * (n_i/n_v)^2
  where (n_i/n_v)^2 = E_outer/E_inner = frequency ratio

  g_same = 2/d = 2/3  (same angular pattern, same n)
  g_diff = 4/(2d+1) = 4/7  (different angular pattern)

This gives the effective charge as:
  Z_eff = Z - sum_i [ (1 - g_i * (n_i/n_v)^2) ]
        = Z - N_inner + sum_i [ g_i * (n_i/n_v)^2 ]

Rewriting: Z_eff = Z_bare + sum_of_couplings
  where Z_bare = Z - N_inner (nuclear charge minus inner quanta count)
  and each inner mode adds back g*(n_i/n_v)^2 of coupling

Physical meaning: each inner energy mode reduces the nuclear charge
by 1 (it occupies potential), but then COUPLES BACK some fraction
through wave resonance. The coupling fraction = g * (n_i/n_v)^2.

FILLED MODE PATTERN:
When all 2(2l+1) angular modes at the same (n,l) are occupied,
they form a COMPLETE spherical harmonic with no angular nodes.
This is like a standing wave sphere — it couples differently to
outer modes because its radial profile is uniform in angle.

In Gauss's law terms: a complete spherical shell looks like a
point charge from outside. In wave terms: a complete harmonic
has only monopole coupling (no higher multipoles).

This should REDUCE the coupling back (g_filled < g_diff),
meaning MORE effective screening.
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

# Inner mode configurations: (n_i, l_i, count, is_subshell_filled)
# A subshell is filled when it has 2*(2l+1) energy quanta
mode_configs = {
    'H':  [],  # no inner modes
    'Li': [(1, 0, 2, True)],  # 1s^2 filled (2 of 2)
    'B':  [(1, 0, 2, True), (2, 0, 2, True)],
    'C':  [(1, 0, 2, True), (2, 0, 2, True), (2, 1, 1, False)],
    'N':  [(1, 0, 2, True), (2, 0, 2, True), (2, 1, 2, False)],
    'O':  [(1, 0, 2, True), (2, 0, 2, True), (2, 1, 3, False)],
    'F':  [(1, 0, 2, True), (2, 0, 2, True), (2, 1, 4, False)],
    'Na': [(1, 0, 2, True), (2, 0, 2, True), (2, 1, 6, True)],
    'Cl': [(1, 0, 2, True), (2, 0, 2, True), (2, 1, 6, True),
           (3, 0, 2, True), (3, 1, 4, False)],
}

# =============================================================================
# COUPLING CONSTANTS FROM d=3
# =============================================================================
g_same = 2.0 / dd          # 2/3: same subshell coupling
g_diff = 4.0 / (2*dd+1)    # 4/7: different angular pattern coupling

print("=" * 80)
print("  WAVE MODE COUPLING MODEL")
print("=" * 80)
print(f"\n  g_same = 2/d = {g_same:.6f}")
print(f"  g_diff = 4/(2d+1) = {g_diff:.6f}")
print(f"  Coupling law: c(n_i -> n_v) = g * (n_i/n_v)^2")
print(f"  Z_eff = Z - N_inner + sum(g_i * (n_i/n_v)^2)")

# =============================================================================
# TEST: g_filled for complete spherical harmonics
# =============================================================================
print()
print("=" * 80)
print("  SEARCHING: g_filled for complete angular mode sets")
print("=" * 80)

# When a subshell (n,l) has all 2(2l+1) quanta, it forms a complete
# spherical harmonic. Test different g_filled values.

# Only Na and Cl have filled p-shells screening the valence
# Na: 2p^6 (filled) screens 3s
# Cl: 2p^6 (filled) screens 3p, plus 1s^2, 2s^2, 3s^2 (filled s-shells)

# GWT constant candidates for g_filled:
candidates = [
    ('1/(2d+1) = 1/7', 1/(2*dd+1)),
    ('1/(d+1) = 1/4', 1/(dd+1)),
    ('1/d = 1/3', 1/dd),
    ('2/(2d+1) = 2/7', 2/(2*dd+1)),
    ('(d-1)/(2d) = 1/3', (dd-1)/(2*dd)),
    ('1/pi', 1/pi),
    ('2/pi^2', 2/pi**2),
    ('alpha = 7/10', alpha_n),
    ('1/(d^2-1) = 1/8', 1/(dd**2-1)),
]

print(f"\n  Testing: filled p-subshells use g_filled instead of g_diff")
print(f"  (Filled s-subshells still use g_diff since s-shell is always 'complete')")
print(f"\n{'g_filled':>25} {'Na_pred':>8} {'Na_real':>8} {'Cl_pred':>8} {'Cl_real':>8} {'rms':>7}")
print("-" * 65)

def compute_zeff(name, g_filled_p):
    info = atoms[name]
    Z = info['Z']; n_v = info['n']; l_v = info['l']
    S = 0
    for (n_i, l_i, cnt, filled) in mode_configs[name]:
        if n_i == n_v and l_i == l_v:
            g = g_same  # same subshell
        elif filled and l_i >= 1:
            # Complete angular mode set (p or higher)
            g = g_filled_p
        else:
            g = g_diff
        coupling = g * (n_i / n_v)**2
        s = 1 - coupling  # screening per quantum
        S += cnt * s
    return Z - S

for label, gf in candidates:
    Na_pred = compute_zeff('Na', gf)
    Cl_pred = compute_zeff('Cl', gf)
    rms = np.sqrt(((Na_pred - 2.5074)**2 + (Cl_pred - 6.1161)**2) / 2)
    print(f"{label:>25} {Na_pred:8.4f} {2.5074:8.4f} {Cl_pred:8.4f} {6.1161:8.4f} {rms:7.4f}")

# Fine scan
print("\nFine scan:")
best = (999, 0)
for gf in np.arange(0.0, 0.8, 0.001):
    Na_p = compute_zeff('Na', gf)
    Cl_p = compute_zeff('Cl', gf)
    errs_all = []
    for name in atoms:
        Ze_p = compute_zeff(name, gf)
        errs_all.append((Ze_p - atoms[name]['Z_eff'])**2)
    rms = np.sqrt(np.mean(errs_all))
    if rms < best[0]:
        best = (rms, gf)

gf_best = best[1]
print(f"  Best g_filled = {gf_best:.4f}, RMS(all atoms) = {best[0]:.4f}")

# Check against GWT constants
print(f"  Compare:")
for label, val in candidates:
    if abs(val - gf_best) < 0.05:
        print(f"    {label} = {val:.4f} (diff: {abs(val-gf_best):.4f})")

# Also check: 2/(d*(2d+1)) = 2/21
print(f"    2/(d*(2d+1)) = 2/21 = {2/(dd*(2*dd+1)):.4f}")
print(f"    1/(d*pi) = {1/(dd*pi):.4f}")
print(f"    (d-1)/(d*(d+1)) = {(dd-1)/(dd*(dd+1)):.4f}")
print(f"    2/(d^2+d) = {2/(dd**2+dd):.4f}")


# =============================================================================
# FULL PREDICTIONS WITH BEST g_filled
# =============================================================================
print()
print("=" * 80)
print(f"  FULL PREDICTIONS: g_same=2/d, g_diff=4/(2d+1), g_filled_p={gf_best:.4f}")
print("=" * 80)

print(f"\n{'Atom':>4} {'Z':>3} {'Z_pred':>7} {'Z_real':>7} {'err':>7} {'E_pred':>8} {'E_real':>8} {'E_err%':>7}")
print("-" * 65)

for name in atoms:
    info = atoms[name]
    Ze_pred = compute_zeff(name, gf_best)
    Ze_real = info['Z_eff']
    n = info['n']
    E_pred = E_H * (Ze_pred/n)**2
    E_real = E_H * (Ze_real/n)**2
    E_err = (E_pred - E_real) / E_real * 100
    print(f"{name:>4} {info['Z']:3d} {Ze_pred:7.4f} {Ze_real:7.4f} {Ze_pred-Ze_real:+7.4f} "
          f"{E_pred:8.4f} {E_real:8.4f} {E_err:+6.1f}%")


# =============================================================================
# COUPLING INTERPRETATION
# =============================================================================
print()
print("=" * 80)
print("  PHYSICAL INTERPRETATION: Wave mode coupling")
print("=" * 80)
print()
print("  Each energy quantum occupies a wave mode (n, l, m).")
print("  The mode's coupling to the nuclear potential is Z_eff/n^2.")
print("  Inner modes REDUCE the effective nuclear coupling by occupying")
print("  potential, but COUPLE BACK through wave resonance.")
print()
print("  Coupling back = g * (n_i/n_v)^2 = g * (frequency ratio)")
print()
print("  g depends on the MODE TYPE MATCH:")
print(f"    g_same = 2/d = {g_same:.4f}: same angular mode, maximum resonance")
print(f"    g_diff = 4/(2d+1) = {g_diff:.4f}: different angular mode, partial resonance")
print(f"    g_filled_p = {gf_best:.4f}: complete angular set, minimal resonance")
print()
print("  WHY complete angular sets couple LESS:")
print("    A complete set of angular modes (all m for given l)")
print("    forms a SPHERICALLY SYMMETRIC energy distribution.")
print("    From outside, this looks like a point source (Gauss's law).")
print("    It has NO angular structure to resonate with the outer mode.")
print("    Only the radial overlap remains, which is weaker.")


# =============================================================================
# BOND FORMULA TEST WITH WAVE-DERIVED Z_eff
# =============================================================================
print()
print("=" * 80)
print("  BOND PREDICTIONS: wave-derived Z_eff vs empirical Z_eff")
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
    info1 = atoms[atom1]; info2 = atoms[atom2]
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

    # Ionic with Z_eff from zeff_func
    Ze1 = zeff_func(atom1); Ze2 = zeff_func(atom2)
    E1 = E_H * (Ze1/info1['n'])**2
    E2 = E_H * (Ze2/info2['n'])**2
    dE = abs(E1 - E2)
    V = max(abs(D_cov), 0.01)
    q = dE / np.sqrt(dE**2 + (2*V)**2) if dE > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / R
    return D_cov, D_cov + Di

wave_zeff = lambda name: compute_zeff(name, gf_best)

print()
print(f"{'Mol':<6} {'De_exp':>7} {'D_wave':>7} {'D_real':>7} {'err_w%':>7} {'err_r%':>7} {'verdict':>8}")
print("-" * 60)

errs_w = []; errs_r = []
for mol in molecules:
    _, D_w = compute_bond(mol, wave_zeff)
    _, D_r = compute_bond(mol, lambda n: Z_eff_real[n])
    err_w = (D_w - mol[2]) / mol[2] * 100
    err_r = (D_r - mol[2]) / mol[2] * 100
    errs_w.append(abs(err_w)); errs_r.append(abs(err_r))
    v = 'BETTER' if abs(err_w) < abs(err_r) - 1 else 'WORSE' if abs(err_w) > abs(err_r) + 1 else '~'
    print(f"{mol[0]:<6} {mol[2]:7.3f} {D_w:7.3f} {D_r:7.3f} {err_w:+6.1f}% {err_r:+6.1f}% {v}")

print(f"\n  Wave-derived: avg={np.mean(errs_w):.1f}%, med={np.median(errs_w):.1f}%, "
      f"<5%:{sum(1 for e in errs_w if e<5)}/24, <10%:{sum(1 for e in errs_w if e<10)}/24")
print(f"  Real Z_eff:   avg={np.mean(errs_r):.1f}%, med={np.median(errs_r):.1f}%, "
      f"<5%:{sum(1 for e in errs_r if e<5)}/24, <10%:{sum(1 for e in errs_r if e<10)}/24")

# Show Z_eff comparison
print()
print("Wave-derived Z_eff:")
print(f"{'Atom':>4} {'Z':>3} {'Z_wave':>7} {'Z_real':>7} {'E_wave':>8} {'E_real':>8}")
print("-" * 45)
for name in atoms:
    info = atoms[name]
    Ze_w = wave_zeff(name)
    Ze_r = Z_eff_real[name]
    n = info['n']
    print(f"{name:>4} {info['Z']:3d} {Ze_w:7.4f} {Ze_r:7.4f} "
          f"{E_H*(Ze_w/n)**2:8.3f} {E_H*(Ze_r/n)**2:8.3f}")


# =============================================================================
# SUMMARY: All constants from d=3
# =============================================================================
print()
print("=" * 80)
print("  SUMMARY: Wave mode coupling constants from d=3")
print("=" * 80)
print(f"""
  Dimension: d = 3

  Coupling law: c(n_i -> n_v) = g * (n_i/n_v)^2
    where (n_i/n_v)^2 is the standing wave frequency ratio

  Three coupling constants:
    g_same    = 2/d       = {g_same:.6f}  (same angular mode)
    g_diff    = 4/(2d+1)  = {g_diff:.6f}  (different angular mode)
    g_filled  = {gf_best:.6f}            (complete angular set, l>=1)

  Z_eff = Z - sum_modes [ 1 - g * (n_i/n_v)^2 ]
        = (Z - N_inner) + sum_modes [ g * (n_i/n_v)^2 ]

  Bond formula uses this Z_eff for ionic energy difference:
    dE = |E_H*(Z1/n1)^2 - E_H*(Z2/n2)^2|
""")
