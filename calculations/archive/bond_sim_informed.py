"""
Bond Predictions — Simulation-Informed Overlap
===============================================
Uses 3D breather simulation results to inform the overlap function.

KEY FINDINGS from breather_3d_gpu simulation:
  - s-wave (l=0): interaction is MONOTONIC, no oscillation
    -> overlap ~ exp(-eta*R) smooth decay, NOT sin(kR)
  - p-wave (l=1): interaction OSCILLATES with sign reversal
    -> overlap ~ sin(k*R + phi0) * exp(-alpha*R)
    -> k ~ eta = sin(n*gamma), half-period ~ pi/eta ~ 8.5 lattice units
    -> first zero at R ~ 3.3, peak bonding at R ~ 7
  - mixed s+p: always repulsive in simulation (angular mismatch)

APPROACH: Replace |sin(phase)| with simulation-informed functions:
  - For s-wave bonds (ss): use smooth overlap S(R) = exp(-decay*R)
  - For p-wave bonds (pp, pi): keep sin-like oscillation but with envelope
  - For mixed (sp): hybrid

All GWT parameters still from d=3, but overlap function changes.
"""

import numpy as np

pi = np.pi
E_H = 13.6057  # eV

# =============================================================================
# GWT PARAMETERS (d=3, unchanged)
# =============================================================================
d = 3
C_bond = pi / d
f_pi   = d**2 / (d**2 + 1)
alpha  = 1 - f_pi / d
beta   = (1 + f_pi) / 2
f_anti = 2*d / (2*d - 1)
c_ionic = 1.0 / (2*d + 1)

Z_eff = {
    'H_1s':  1.0000, 'Li_2s': 1.2792, 'B_2p':  2.4214,
    'C_2p':  3.1358, 'N_2p':  3.8340, 'O_2p':  4.4532,
    'F_2p':  5.0998, 'Na_3s': 2.5074, 'Cl_3p': 6.1161,
}

molecules = [
    ('H2',   1.401,  4.745, [('ss', 1)],                                           'H_1s',  'H_1s'),
    ('Li2',  5.051,  1.056, [('ss', 1)],                                           'Li_2s', 'Li_2s'),
    ('B2',   3.005,  3.02,  [('pi', 2)],                                           'B_2p',  'B_2p'),
    ('C2',   2.348,  6.32,  [('pi', 2), ('sp_sigma', 1), ('sp_sigma_anti', 1)],    'C_2p',  'C_2p'),
    ('N2',   2.074,  9.759, [('pp_sigma', 1), ('pi', 2)],                          'N_2p',  'N_2p'),
    ('O2',   2.282,  5.213, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)],          'O_2p',  'O_2p'),
    ('F2',   2.668,  1.660, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)],          'F_2p',  'F_2p'),
    ('Na2',  5.818,  0.746, [('ss', 1)],                                           'Na_3s', 'Na_3s'),
    ('Cl2',  3.757,  2.514, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)],          'Cl_3p', 'Cl_3p'),
    ('HF',   1.733,  5.869, [('sp', 1)],                                           'H_1s',  'F_2p'),
    ('CO',   2.132, 11.225, [('pp_sigma', 1), ('pi', 2)],                          'C_2p',  'O_2p'),
    ('NO',   2.175,  6.497, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)],          'N_2p',  'O_2p'),
    ('OH',   1.834,  4.392, [('sp', 1)],                                           'O_2p',  'H_1s'),
    ('HCl',  2.409,  4.434, [('sp', 1)],                                           'H_1s',  'Cl_3p'),
    ('LiH',  3.015,  2.515, [('ss', 1)],                                           'Li_2s', 'H_1s'),
    ('LiF',  2.955,  5.939, [('sp', 1)],                                           'Li_2s', 'F_2p'),
    ('BH',   2.329,  3.42,  [('sp', 1)],                                           'B_2p',  'H_1s'),
    ('CH',   2.116,  3.47,  [('sp', 1)],                                           'C_2p',  'H_1s'),
    ('NH',   1.958,  3.57,  [('sp', 1)],                                           'N_2p',  'H_1s'),
    ('BF',   2.386,  7.81,  [('pp_sigma', 1), ('pi', 2)],                          'B_2p',  'F_2p'),
    ('CN',   2.214,  7.72,  [('pp_sigma', 1), ('pi', 2)],                          'C_2p',  'N_2p'),
    ('NaH',  3.566,  1.97,  [('ss', 1)],                                           'Na_3s', 'H_1s'),
    ('NaCl', 4.461,  4.23,  [('sp', 1)],                                           'Na_3s', 'Cl_3p'),
    ('H2O',  1.809,  5.117, [('sp', 1)],                                           'O_2p',  'H_1s'),
]

def get_n(orb):
    return int(orb.split('_')[1][0])

def get_l(orb):
    return {'s': 0, 'p': 1}[orb.split('_')[1][1]]

def has_nodes(orb):
    n, l = get_n(orb), get_l(orb)
    return min(n - l - 1, 1)

def orbital_energy(orb):
    n, z = get_n(orb), Z_eff[orb]
    return E_H * (z / n)**2

def bond_order(bonds):
    bo = 0
    for bt, cnt in bonds:
        bo += cnt if 'anti' not in bt else -cnt
    return bo


# =============================================================================
# OVERLAP FUNCTIONS — from simulation
# =============================================================================

def overlap_original(phase):
    """Original GWT: |sin(phase)|"""
    return abs(np.sin(phase))


def overlap_sim_v1(phase, l1, l2):
    """Simulation-informed overlap V1.

    s-wave (l=0): simulation shows monotonic decay, no oscillation.
      Use 1/(1 + phase^2) — Lorentzian-like decay.
      At small phase this ~ 1, at large phase -> 0.

    p-wave (l=1): simulation shows oscillation with sign reversal.
      Keep sin-like behavior but with exponential envelope.
      Use |sin(phase)| * exp(-phase/pi) (damped oscillation).

    sp mixed: simulation shows always repulsive with simple superposition.
      But real sp bonds work — use geometric mean of s and p overlaps.
    """
    if l1 == 0 and l2 == 0:
        # s-wave: smooth monotonic decay
        # Normalize: at phase=0 overlap=1, decays smoothly
        return 1.0 / (1.0 + (phase / pi)**2)
    elif l1 == 1 and l2 == 1:
        # p-wave: damped oscillation
        return abs(np.sin(phase)) * np.exp(-0.1 * phase)
    else:
        # sp mixed: geometric mean of s and p behavior
        s_part = 1.0 / (1.0 + (phase / pi)**2)
        p_part = abs(np.sin(phase)) * np.exp(-0.1 * phase)
        return np.sqrt(s_part * p_part)


def overlap_sim_v2(phase, l1, l2):
    """V2: Use simulation shape more directly.

    s-wave: exp(-phase/pi) — simple exponential decay
    p-wave: |sin(phase)| unchanged (it works well for pp bonds!)
    sp: the p-wave part oscillates, s-wave part decays
        -> use p-wave sin(phase) but modulated by s-wave decay
    """
    if l1 == 0 and l2 == 0:
        # s-wave: exponential decay
        return np.exp(-phase / pi)
    elif l1 == 1 and l2 == 1:
        # p-wave: keep original (works for N2, O2, F2, etc.)
        return abs(np.sin(phase))
    else:
        # sp mixed: sin modulated by decay
        return abs(np.sin(phase)) * np.exp(-0.5 * phase / pi)


def overlap_sim_v3(phase, l1, l2, h1, h2):
    """V3: Node-aware overlap.

    Key insight from simulation: the breather width ~ 1/sin(n*gamma).
    Higher modes are narrower, so their overlap decays faster.
    Nodes (h=1) reduce the effective overlap range.

    s-wave with nodes (Li_2s, Na_3s): faster decay, like exp(-phase)
    s-wave no nodes (H_1s): slower decay, like exp(-phase/pi)
    p-wave: |sin(phase)| with decay for nodal orbitals
    sp: blend

    Also: use |sin(phase)| / ceil(phase/pi) for multi-lobe averaging
    when phase > pi (node counting from previous analysis).
    """
    has_any_nodes = (h1 + h2 > 0)

    if l1 == 0 and l2 == 0:
        # ss bonds
        if has_any_nodes:
            # Nodal s-wave: faster decay
            return np.exp(-phase * 0.5)
        else:
            # Nodeless s-wave (H-H): moderate decay
            return np.exp(-phase / pi)
    elif l1 == 1 and l2 == 1:
        # pp bonds: oscillatory
        base = abs(np.sin(phase))
        if has_any_nodes and phase > pi:
            # Node counting: average over lobes
            n_lobes = int(np.ceil(phase / pi))
            base = base / n_lobes
        return base
    else:
        # sp bonds: sin with decay
        base = abs(np.sin(phase))
        if has_any_nodes and phase > pi:
            n_lobes = int(np.ceil(phase / pi))
            base = base / n_lobes
        decay = np.exp(-0.3 * phase / pi) if has_any_nodes else 1.0
        return base * decay


def compute_bond_energy(name, R, De_exp, bonds, orb1, orb2, overlap_fn='original'):
    """Compute bond energy with selectable overlap function."""
    n1, l1 = get_n(orb1), get_l(orb1)
    n2, l2 = get_n(orb2), get_l(orb2)
    h1, h2 = has_nodes(orb1), has_nodes(orb2)

    a1 = 2 + (1 - 2*l1) * alpha * h1
    a2 = 2 + (1 - 2*l2) * alpha * h2
    b1 = 1 + beta * h1
    b2 = 1 + beta * h2

    E1 = E_H / n1**a1
    E2 = E_H / n2**a2
    E_scale = np.sqrt(E1 * E2)

    sigma_phase = R / n1**b1 + R / n2**b2

    n_pi_bond = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    n_pi_anti = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    is_full_anti = (n_pi_anti >= n_pi_bond) if n_pi_bond > 0 else True

    D_cov = 0
    for btype, count in bonds:
        if 'sigma' in btype or btype in ('ss', 'sp'):
            phase = sigma_phase
        else:
            phase = sigma_phase * f_pi

        # Select overlap function
        if overlap_fn == 'original':
            S = overlap_original(phase)
        elif overlap_fn == 'v1':
            # Determine l for this bond type
            bl1 = l1 if btype in ('ss', 'sp') else l1
            bl2 = l2 if btype in ('ss', 'sp') else l2
            if btype == 'ss':
                bl1, bl2 = 0, 0
            elif btype == 'sp':
                bl1, bl2 = min(l1, l2), max(l1, l2)
            elif 'pi' in btype or 'pp' in btype:
                bl1, bl2 = 1, 1
            S = overlap_sim_v1(phase, bl1, bl2)
        elif overlap_fn == 'v2':
            if btype == 'ss':
                bl1, bl2 = 0, 0
            elif btype == 'sp':
                bl1, bl2 = min(l1, l2), max(l1, l2)
            else:
                bl1, bl2 = 1, 1
            S = overlap_sim_v2(phase, bl1, bl2)
        elif overlap_fn == 'v3':
            if btype == 'ss':
                bl1, bl2 = 0, 0
            elif btype == 'sp':
                bl1, bl2 = min(l1, l2), max(l1, l2)
            else:
                bl1, bl2 = 1, 1
            S = overlap_sim_v3(phase, bl1, bl2, h1, h2)
        else:
            S = overlap_original(phase)

        contribution = C_bond * E_scale * S

        if 'anti' in btype:
            if 'sigma' in btype or is_full_anti:
                f_a = 1.0
            else:
                f_a = f_anti
            D_cov -= count * f_a * contribution
        else:
            D_cov += count * contribution

    # Ionic correction (unchanged)
    eps1 = orbital_energy(orb1)
    eps2 = orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    D_ionic = c_ionic * q**2 * 2 * E_H / R

    return D_cov + D_ionic, D_cov, D_ionic, q


# =============================================================================
# MAIN: Compare all overlap models
# =============================================================================

print("=" * 95)
print("  SIMULATION-INFORMED BOND PREDICTIONS")
print("=" * 95)
print()
print("Comparing overlap functions:")
print("  Original: |sin(phase)|  for all bonds")
print("  V1: s-wave -> Lorentzian decay, p-wave -> damped sin, sp -> geometric mean")
print("  V2: s-wave -> exp decay, p-wave -> |sin| (unchanged), sp -> sin*exp")
print("  V3: Node-aware: nodal s -> fast decay, pp node counting, sp node decay")
print()

for version in ['original', 'v1', 'v2', 'v3']:
    print(f"\n{'='*95}")
    print(f"  VERSION: {version.upper()}")
    print(f"{'='*95}")

    header = f"{'Mol':<7} {'De_exp':>7} {'De_pred':>7} {'err%':>7}  {'D_cov':>7} {'D_ion':>6} {'q':>5}"
    print(header)
    print("-" * 60)

    errs = []
    for mol in molecules:
        name, R, De_exp, bonds, o1, o2 = mol
        De_pred, D_cov, D_ion, q = compute_bond_energy(*mol, overlap_fn=version)
        err = (De_pred - De_exp) / De_exp * 100
        errs.append(abs(err))
        flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
        print(f"{name:<7} {De_exp:7.3f} {De_pred:7.3f} {err:+6.1f}%  "
              f"{D_cov:7.3f} {D_ion:6.3f} {q:5.3f} {flag}")

    print("-" * 60)
    n = len(errs)
    w2 = sum(1 for e in errs if e < 2)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"  avg={np.mean(errs):.1f}%, median={np.median(errs):.1f}%")
    print(f"  w2={w2}/{n}, w5={w5}/{n}, w10={w10}/{n}, w20={w20}/{n}")

    # Show top outliers
    indexed = sorted(enumerate(errs), key=lambda x: -x[1])
    print(f"  Top 5 outliers:")
    for idx, err_val in indexed[:5]:
        mol = molecules[idx]
        print(f"    {mol[0]:<7} {err_val:.1f}%")


# =============================================================================
# DETAILED COMPARISON: Show how each outlier changes
# =============================================================================
print(f"\n\n{'='*95}")
print(f"  OUTLIER COMPARISON: How each version affects the 7 problem molecules")
print(f"{'='*95}")

outliers = ['CH', 'LiH', 'NaH', 'LiF', 'NaCl', 'BF', 'CN']
print(f"\n{'Mol':<7} {'De_exp':>7} ", end='')
for v in ['original', 'v1', 'v2', 'v3']:
    print(f"  {v:>10} {'err':>7}", end='')
print()
print("-" * 90)

for mol in molecules:
    name = mol[0]
    if name not in outliers:
        continue
    De_exp = mol[2]
    print(f"{name:<7} {De_exp:7.3f} ", end='')
    for v in ['original', 'v1', 'v2', 'v3']:
        De_pred = compute_bond_energy(*mol, overlap_fn=v)[0]
        err = (De_pred - De_exp) / De_exp * 100
        print(f"  {De_pred:10.3f} {err:+6.1f}%", end='')
    print()


# =============================================================================
# Check: Do the good molecules stay good?
# =============================================================================
print(f"\n\n{'='*95}")
print(f"  STABILITY CHECK: Do good molecules (within 5%) stay good?")
print(f"{'='*95}")

good_mols = []
for i, mol in enumerate(molecules):
    De_pred_orig = compute_bond_energy(*mol, overlap_fn='original')[0]
    err_orig = abs((De_pred_orig - mol[2]) / mol[2] * 100)
    if err_orig < 5:
        good_mols.append(mol[0])

print(f"\nGood molecules in original: {', '.join(good_mols)}")
print(f"\n{'Mol':<7} {'De_exp':>7} {'orig':>10} {'orig%':>7}", end='')
for v in ['v1', 'v2', 'v3']:
    print(f"  {v:>10} {v+'%':>7}", end='')
print()
print("-" * 85)

for mol in molecules:
    if mol[0] not in good_mols:
        continue
    De_exp = mol[2]
    De_orig = compute_bond_energy(*mol, overlap_fn='original')[0]
    err_orig = (De_orig - De_exp) / De_exp * 100
    print(f"{mol[0]:<7} {De_exp:7.3f} {De_orig:10.3f} {err_orig:+6.1f}%", end='')
    for v in ['v1', 'v2', 'v3']:
        De_v = compute_bond_energy(*mol, overlap_fn=v)[0]
        err_v = (De_v - De_exp) / De_exp * 100
        broke = ' BROKE' if abs(err_v) > 5 and abs(err_orig) < 5 else ''
        print(f"  {De_v:10.3f} {err_v:+6.1f}%{broke}", end='')
    print()
