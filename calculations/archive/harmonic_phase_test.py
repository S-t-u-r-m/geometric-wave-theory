"""
HARMONIC MEAN PHASE + ONE FORCE
================================
Key findings from one_force_test.py:
1. The "ionic" enhancement correlates at +0.979 with 1/|sin(phase)|
   -> The ionic term is compensating for sin hitting nodes, not a separate force
2. Harmonic mean phase is the best alternative: avg 20.3%, 13/24 within 10%
   -> The overlap is limited by the SLOWER wavevector (larger orbital)

Physics of harmonic mean:
  k_harm = 2*k1*k2/(k1+k2) = harmonic mean wavevector
  phase = 2*R*k_harm = 4*R*k1*k2/(k1+k2)

  For homonuclear (k1=k2=k): phase = 2kR (same as current)
  For heteronuclear: gives MORE weight to smaller k (larger orbital)

  Physical meaning: the overlap integral is dominated by the orbital
  that extends less far. The smaller k (larger orbital) determines
  the effective wavelength of the overlap region.

Also test: 3 coupling modes per bond (d=3 axes)
  Each axis contributes independently to the bonding
  sigma: 1 axis (along bond), pi: 2 axes (perpendicular)
  What if each has its own coupling strength?
"""
import numpy as np

pi = np.pi
E_H = 13.6057
d = 3
C_bond = pi / d
f_pi = d**2 / (d**2 + 1)
alpha_n = 1 - f_pi / d
beta_n = (1 + f_pi) / 2
f_anti = 2*d / (2*d - 1)
c_ionic = 1.0 / (2*d + 1)

Z_eff = {
    'H_1s': 1.0, 'Li_2s': 1.2792, 'B_2p': 2.4214,
    'C_2p': 3.1358, 'N_2p': 3.8340, 'O_2p': 4.4532,
    'F_2p': 5.0998, 'Na_3s': 2.5074, 'Cl_3p': 6.1161,
}

molecules = [
    ('H2',   1.401,  4.745, [('ss', 1)], 'H_1s',  'H_1s'),
    ('Li2',  5.051,  1.056, [('ss', 1)], 'Li_2s', 'Li_2s'),
    ('B2',   3.005,  3.02,  [('pi', 2)], 'B_2p',  'B_2p'),
    ('C2',   2.348,  6.32,  [('pi', 2), ('sp_sigma', 1), ('sp_sigma_anti', 1)], 'C_2p',  'C_2p'),
    ('N2',   2.074,  9.759, [('pp_sigma', 1), ('pi', 2)], 'N_2p',  'N_2p'),
    ('O2',   2.282,  5.213, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)], 'O_2p',  'O_2p'),
    ('F2',   2.668,  1.660, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)], 'F_2p',  'F_2p'),
    ('Na2',  5.818,  0.746, [('ss', 1)], 'Na_3s', 'Na_3s'),
    ('Cl2',  3.757,  2.514, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)], 'Cl_3p', 'Cl_3p'),
    ('HF',   1.733,  5.869, [('sp', 1)], 'H_1s',  'F_2p'),
    ('CO',   2.132, 11.225, [('pp_sigma', 1), ('pi', 2)], 'C_2p',  'O_2p'),
    ('NO',   2.175,  6.497, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)], 'N_2p',  'O_2p'),
    ('OH',   1.834,  4.392, [('sp', 1)], 'O_2p',  'H_1s'),
    ('HCl',  2.409,  4.434, [('sp', 1)], 'H_1s',  'Cl_3p'),
    ('LiH',  3.015,  2.515, [('ss', 1)], 'Li_2s', 'H_1s'),
    ('LiF',  2.955,  5.939, [('sp', 1)], 'Li_2s', 'F_2p'),
    ('BH',   2.329,  3.42,  [('sp', 1)], 'B_2p',  'H_1s'),
    ('CH',   2.116,  3.47,  [('sp', 1)], 'C_2p',  'H_1s'),
    ('NH',   1.958,  3.57,  [('sp', 1)], 'N_2p',  'H_1s'),
    ('BF',   2.386,  7.81,  [('pp_sigma', 1), ('pi', 2)], 'B_2p',  'F_2p'),
    ('CN',   2.214,  7.72,  [('pp_sigma', 1), ('pi', 2)], 'C_2p',  'N_2p'),
    ('NaH',  3.566,  1.97,  [('ss', 1)], 'Na_3s', 'H_1s'),
    ('NaCl', 4.461,  4.23,  [('sp', 1)], 'Na_3s', 'Cl_3p'),
    ('H2O',  1.809,  5.117, [('sp', 1)], 'O_2p',  'H_1s'),
]


def get_n(orb): return int(orb.split('_')[1][0])
def get_l(orb): return {'s': 0, 'p': 1}[orb.split('_')[1][1]]
def has_nodes(orb): return min(get_n(orb) - get_l(orb) - 1, 1)
def orb_energy(orb): return E_H * (Z_eff[orb] / get_n(orb))**2


def compute_model(mol, phase_type='current', C=None, include_ionic=True,
                  use_phase_fix=True):
    """
    Flexible bond computation.
    phase_type: 'current' (k1+k2)*R, 'harmonic' 4k1k2R/(k1+k2), 'geometric' 2R*sqrt(k1k2)
    """
    if C is None:
        C = C_bond

    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
    n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)

    a1 = 2 + (1 - 2*l1) * alpha_n * h1
    a2 = 2 + (1 - 2*l2) * alpha_n * h2

    if use_phase_fix:
        b1 = 1 + beta_n * h1
        b2 = 1 + beta_n * h2
    else:
        b1 = 1 + (1 - 2*l1) * beta_n * h1
        b2 = 1 + (1 - 2*l2) * beta_n * h2

    E_scale = np.sqrt(E_H / n1**a1 * E_H / n2**a2)
    k1 = 1.0 / n1**b1
    k2 = 1.0 / n2**b2

    if phase_type == 'current':
        sigma_phase = R * (k1 + k2)
    elif phase_type == 'harmonic':
        sigma_phase = 4 * R * k1 * k2 / (k1 + k2) if (k1 + k2) > 0 else 0
    elif phase_type == 'geometric':
        sigma_phase = 2 * R * np.sqrt(k1 * k2)
    else:
        sigma_phase = R * (k1 + k2)

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss', 'sp')) else sigma_phase * f_pi
        cont = C * E_scale * abs(np.sin(ph))
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    D_ion = 0
    if include_ionic:
        eps1 = E_H * (Z_eff[orb1] / n1)**2
        eps2 = E_H * (Z_eff[orb2] / n2)**2
        de = abs(eps1 - eps2)
        Vc = max(abs(D_cov), 0.01)
        q = de / np.sqrt(de**2 + (2*Vc)**2) if de > 0 else 0
        D_ion = c_ionic * q**2 * 2 * E_H / R

    return D_cov + D_ion, D_cov, D_ion, sigma_phase


# =============================================================================
# HEAD TO HEAD: current vs harmonic, with and without ionic
# =============================================================================
configs = [
    ('Current + ionic',    'current',  True,  True),
    ('Current, NO ionic',  'current',  False, True),
    ('Harmonic + ionic',   'harmonic', True,  True),
    ('Harmonic, NO ionic', 'harmonic', False, True),
    ('Geometric + ionic',  'geometric', True, True),
    ('Geometric, NO ionic','geometric', False, True),
    # Also with old phase convention
    ('Current+ion (old b)', 'current', True,  False),
    ('Harmonic+ion (old b)','harmonic', True,  False),
]

print("=" * 100)
print("  COMPARISON: Phase models x ionic on/off x phase fix on/off")
print("=" * 100)
print()
print(f"{'Config':<25} {'avg%':>6} {'med%':>6} {'<5%':>5} {'<10%':>5} {'<20%':>5}")
print("-" * 60)

best_config = None
best_avg = 999

for label, pt, ionic, pfix in configs:
    errs = []
    for mol in molecules:
        D_tot = compute_model(mol, phase_type=pt, include_ionic=ionic,
                              use_phase_fix=pfix)[0]
        err = abs((D_tot - mol[2]) / mol[2] * 100)
        errs.append(err)

    avg = np.mean(errs)
    med = np.median(errs)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"{label:<25} {avg:6.1f} {med:6.1f} {w5:5d} {w10:5d} {w20:5d}")

    if avg < best_avg:
        best_avg = avg
        best_config = (label, pt, ionic, pfix)


# =============================================================================
# DETAILED VIEW: harmonic, no ionic (one force)
# =============================================================================
print()
print("=" * 100)
print("  HARMONIC MEAN PHASE, NO IONIC (one force)")
print("  phase = 4*R*k1*k2/(k1+k2),  b = 1 + beta*h (fixed)")
print("=" * 100)
print()

print(f"{'Mol':<6} {'De_exp':>7} {'D_pred':>7} {'err%':>7} {'ph/pi':>6} {'k1':>6} {'k2':>6} {'k_h':>6}")
print("-" * 65)

errs_harm = []
for mol in molecules:
    name, R = mol[0], mol[1]
    De_exp = mol[2]
    orb1, orb2 = mol[4], mol[5]

    D_tot, D_cov, D_ion, phase = compute_model(mol, 'harmonic', include_ionic=False)
    err = (D_tot - De_exp) / De_exp * 100
    errs_harm.append(abs(err))

    n1 = get_n(orb1); h1 = has_nodes(orb1)
    n2 = get_n(orb2); h2 = has_nodes(orb2)
    b1 = 1 + beta_n * h1
    b2 = 1 + beta_n * h2
    k1 = 1.0 / n1**b1
    k2 = 1.0 / n2**b2
    k_h = 2*k1*k2/(k1+k2) if (k1+k2) > 0 else 0

    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
    print(f"{name:<6} {De_exp:7.3f} {D_tot:7.3f} {err:+6.1f}% {phase/pi:6.3f} "
          f"{k1:6.3f} {k2:6.3f} {k_h:6.3f} {flag}")

w5 = sum(1 for e in errs_harm if e < 5)
w10 = sum(1 for e in errs_harm if e < 10)
w20 = sum(1 for e in errs_harm if e < 20)
print(f"\nAvg={np.mean(errs_harm):.1f}%, Med={np.median(errs_harm):.1f}%, "
      f"<5%:{w5}/24, <10%:{w10}/24, <20%:{w20}/24")


# =============================================================================
# DETAILED VIEW: harmonic + ionic (best of both worlds?)
# =============================================================================
print()
print("=" * 100)
print("  HARMONIC MEAN PHASE + IONIC")
print("=" * 100)
print()

print(f"{'Mol':<6} {'De_exp':>7} {'D_cov':>7} {'D_ion':>7} {'D_tot':>7} {'err%':>7} {'ph/pi':>6}")
print("-" * 60)

errs_hi = []
for mol in molecules:
    name = mol[0]
    De_exp = mol[2]
    D_tot, D_cov, D_ion, phase = compute_model(mol, 'harmonic', include_ionic=True)
    err = (D_tot - De_exp) / De_exp * 100
    errs_hi.append(abs(err))
    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
    print(f"{name:<6} {De_exp:7.3f} {D_cov:7.3f} {D_ion:7.3f} {D_tot:7.3f} {err:+6.1f}% "
          f"{phase/pi:6.3f} {flag}")

w5 = sum(1 for e in errs_hi if e < 5)
w10 = sum(1 for e in errs_hi if e < 10)
w20 = sum(1 for e in errs_hi if e < 20)
print(f"\nAvg={np.mean(errs_hi):.1f}%, Med={np.median(errs_hi):.1f}%, "
      f"<5%:{w5}/24, <10%:{w10}/24, <20%:{w20}/24")


# =============================================================================
# HEAD TO HEAD: current+ionic vs harmonic no ionic
# =============================================================================
print()
print("=" * 100)
print("  HEAD TO HEAD: Current+ionic vs Harmonic (no ionic)")
print("=" * 100)
print()

print(f"{'Mol':<6} {'De_exp':>7} {'curr+ion':>9} {'e_ci%':>7} {'harm':>9} {'e_h%':>7} {'better':>8}")
print("-" * 65)

nb = 0; nw = 0; ns = 0
for i, mol in enumerate(molecules):
    name = mol[0]
    De_exp = mol[2]
    D_ci = compute_model(mol, 'current', include_ionic=True, use_phase_fix=True)[0]
    D_h = compute_model(mol, 'harmonic', include_ionic=False, use_phase_fix=True)[0]
    e_ci = abs((D_ci - De_exp) / De_exp * 100)
    e_h = abs((D_h - De_exp) / De_exp * 100)

    if e_h < e_ci - 1:
        v = 'HARM'; nb += 1
    elif e_h > e_ci + 1:
        v = 'CURR'; nw += 1
    else:
        v = '~'; ns += 1

    print(f"{name:<6} {De_exp:7.3f} {D_ci:9.3f} {e_ci:6.1f}% {D_h:9.3f} {e_h:6.1f}% {v:>8}")

print(f"\nHarmonic wins: {nb}, Current wins: {nw}, Ties: {ns}")


# =============================================================================
# WHAT ABOUT 3 COUPLING MODES? (d=3 axes)
# =============================================================================
print()
print("=" * 100)
print("  3-AXIS MODEL: Each spatial axis contributes a coupling mode")
print("=" * 100)
print()
print("In d=3, a bond involves 3 spatial axes:")
print("  - 1 axis along bond direction (sigma)")
print("  - 2 axes perpendicular (pi)")
print("Current: sigma uses full phase, pi uses f_pi*phase")
print()
print("What if EACH axis has its own resonance condition?")
print("  sigma axis: sin(k*R)   -- wave along bond")
print("  pi axis 1:  sin(k*R*f_pi) -- wave perpendicular")
print("  pi axis 2:  sin(k*R*f_pi) -- wave perpendicular")
print()
print("For a single bond, the total resonance = sum of all 3 axes?")
print("  D = C * E_scale * [sin(ph_sigma) + 2*sin(ph_pi)] / 3")
print("  = C * E_scale * [sin(ph) + 2*sin(f_pi*ph)] / 3")
print()
print("For homonuclear:")
print("  H2 (ss): sin(2.802) + 2*sin(0.9*2.802) = sin(2.802) + 2*sin(2.522)")
print(f"         = {np.sin(2.802):.4f} + 2*{np.sin(0.9*2.802):.4f} = {np.sin(2.802) + 2*np.sin(0.9*2.802):.4f}")
print(f"  Current: sin(2.802) = {np.sin(2.802):.4f}")
print(f"  Ratio: {(np.sin(2.802) + 2*np.sin(0.9*2.802))/(3*np.sin(2.802)):.4f}")
print()

# Test: 3-axis model with various normalizations
def compute_3axis(mol, use_phase_fix=True, include_ionic=False):
    """3-axis coupling: each axis contributes independently."""
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
    n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)

    a1 = 2 + (1 - 2*l1) * alpha_n * h1
    a2 = 2 + (1 - 2*l2) * alpha_n * h2
    b1 = 1 + beta_n * h1 if use_phase_fix else 1 + (1-2*l1)*beta_n*h1
    b2 = 1 + beta_n * h2 if use_phase_fix else 1 + (1-2*l2)*beta_n*h2

    E_scale = np.sqrt(E_H / n1**a1 * E_H / n2**a2)
    k1 = 1.0 / n1**b1
    k2 = 1.0 / n2**b2
    sigma_phase = R * (k1 + k2)
    pi_phase = sigma_phase * f_pi

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        if 'sigma' in bt or bt in ('ss', 'sp'):
            # Sigma bond: resonance on all 3 axes
            # Along bond: sin(sigma_phase)
            # Perpendicular: sin(pi_phase) x 2
            overlap = (abs(np.sin(sigma_phase)) + 2*abs(np.sin(pi_phase))) / d
        elif 'pi' in bt:
            # Pi bond: primarily perpendicular
            overlap = abs(np.sin(pi_phase))
        else:
            overlap = abs(np.sin(sigma_phase))

        cont = C_bond * E_scale * overlap
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    D_ion = 0
    if include_ionic:
        eps1 = E_H * (Z_eff[orb1] / n1)**2
        eps2 = E_H * (Z_eff[orb2] / n2)**2
        de = abs(eps1 - eps2)
        Vc = max(abs(D_cov), 0.01)
        q = de / np.sqrt(de**2 + (2*Vc)**2) if de > 0 else 0
        D_ion = c_ionic * q**2 * 2 * E_H / R

    return D_cov + D_ion, sigma_phase

print("--- 3-axis model (no ionic) ---")
errs_3a = []
for mol in molecules:
    D_tot, phase = compute_3axis(mol, include_ionic=False)
    err = abs((D_tot - mol[2]) / mol[2] * 100)
    errs_3a.append(err)

w5 = sum(1 for e in errs_3a if e < 5)
w10 = sum(1 for e in errs_3a if e < 10)
w20 = sum(1 for e in errs_3a if e < 20)
print(f"  avg={np.mean(errs_3a):.1f}%, med={np.median(errs_3a):.1f}%, "
      f"<5%:{w5}/24, <10%:{w10}/24, <20%:{w20}/24")

print("--- 3-axis model + ionic ---")
errs_3ai = []
for mol in molecules:
    D_tot, phase = compute_3axis(mol, include_ionic=True)
    err = abs((D_tot - mol[2]) / mol[2] * 100)
    errs_3ai.append(err)

w5 = sum(1 for e in errs_3ai if e < 5)
w10 = sum(1 for e in errs_3ai if e < 10)
w20 = sum(1 for e in errs_3ai if e < 20)
print(f"  avg={np.mean(errs_3ai):.1f}%, med={np.median(errs_3ai):.1f}%, "
      f"<5%:{w5}/24, <10%:{w10}/24, <20%:{w20}/24")


# =============================================================================
# DEEPER: What if sigma bonds get contribution from ALL 3 axes?
# =============================================================================
print()
print("=" * 100)
print("  SIGMA BONDS AS 3-AXIS RESONANCE")
print("  D_sigma = C * E_scale * (1/d) * sum_i |sin(phase_i)|")
print("  where phase_sigma = phase, phase_pi = f_pi * phase")
print("=" * 100)
print()

# What this changes for the problematic molecules:
# For HF (phase=0.827*pi):
#   current: |sin(0.827*pi)| = 0.507
#   3-axis: (|sin(0.827*pi)| + 2*|sin(0.744*pi)|)/3 = (0.507 + 2*0.716)/3 = 0.647
# Bigger! About 28% enhancement. That's the right direction for HF.

for mol in molecules:
    name, R = mol[0], mol[1]
    De_exp = mol[2]
    orb1, orb2 = mol[4], mol[5]

    # Check if this is a sigma-type bond molecule
    is_sigma = any(bt in ('ss', 'sp') or 'sigma' in bt
                   for bt, _ in mol[3] if 'anti' not in bt)
    if not is_sigma:
        continue

    D_curr = compute_model(mol, 'current', include_ionic=False)[0]
    D_3ax, phase = compute_3axis(mol, include_ionic=False)

    ph_s = phase
    ph_p = phase * f_pi
    sin_s = abs(np.sin(ph_s))
    sin_p = abs(np.sin(ph_p))
    avg_3 = (sin_s + 2*sin_p) / 3

    e_curr = (D_curr - De_exp) / De_exp * 100
    e_3ax = (D_3ax - De_exp) / De_exp * 100

    print(f"  {name:<6} ph/pi={ph_s/pi:.3f}  sin_s={sin_s:.3f}  sin_pi={sin_p:.3f}  "
          f"avg3={avg_3:.3f}  e_curr={e_curr:+.1f}%  e_3ax={e_3ax:+.1f}%")


# =============================================================================
# COMBINED: harmonic phase + 3 axis
# =============================================================================
print()
print("=" * 100)
print("  COMBINED: harmonic phase + 3-axis sigma coupling")
print("=" * 100)
print()

def compute_combined(mol, include_ionic=False):
    """Harmonic mean phase + 3-axis sigma coupling."""
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
    n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)

    a1 = 2 + (1 - 2*l1) * alpha_n * h1
    a2 = 2 + (1 - 2*l2) * alpha_n * h2
    b1 = 1 + beta_n * h1
    b2 = 1 + beta_n * h2

    E_scale = np.sqrt(E_H / n1**a1 * E_H / n2**a2)
    k1 = 1.0 / n1**b1
    k2 = 1.0 / n2**b2

    # Harmonic mean phase
    sigma_phase = 4 * R * k1 * k2 / (k1 + k2) if (k1+k2) > 0 else 0
    pi_phase = sigma_phase * f_pi

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        if 'sigma' in bt or bt in ('ss', 'sp'):
            # 3-axis resonance for sigma
            overlap = (abs(np.sin(sigma_phase)) + 2*abs(np.sin(pi_phase))) / d
        else:
            overlap = abs(np.sin(pi_phase))

        cont = C_bond * E_scale * overlap
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    D_ion = 0
    if include_ionic:
        eps1 = E_H * (Z_eff[orb1] / n1)**2
        eps2 = E_H * (Z_eff[orb2] / n2)**2
        de = abs(eps1 - eps2)
        Vc = max(abs(D_cov), 0.01)
        q = de / np.sqrt(de**2 + (2*Vc)**2) if de > 0 else 0
        D_ion = c_ionic * q**2 * 2 * E_H / R

    return D_cov + D_ion, D_cov, D_ion, sigma_phase

print(f"{'Mol':<6} {'De_exp':>7} {'D_pred':>7} {'err%':>7} {'ph/pi':>6}")
print("-" * 40)

errs_comb = []
for mol in molecules:
    D_tot, D_cov, D_ion, phase = compute_combined(mol, include_ionic=False)
    err = (D_tot - mol[2]) / mol[2] * 100
    errs_comb.append(abs(err))
    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
    print(f"{mol[0]:<6} {mol[2]:7.3f} {D_tot:7.3f} {err:+6.1f}% {phase/pi:6.3f} {flag}")

w5 = sum(1 for e in errs_comb if e < 5)
w10 = sum(1 for e in errs_comb if e < 10)
w20 = sum(1 for e in errs_comb if e < 20)
print(f"\nAvg={np.mean(errs_comb):.1f}%, Med={np.median(errs_comb):.1f}%, "
      f"<5%:{w5}/24, <10%:{w10}/24, <20%:{w20}/24")

# With ionic too
print()
print("Same but with ionic:")
errs_combi = []
for mol in molecules:
    D_tot, D_cov, D_ion, phase = compute_combined(mol, include_ionic=True)
    err = abs((D_tot - mol[2]) / mol[2] * 100)
    errs_combi.append(err)

w5 = sum(1 for e in errs_combi if e < 5)
w10 = sum(1 for e in errs_combi if e < 10)
print(f"  avg={np.mean(errs_combi):.1f}%, med={np.median(errs_combi):.1f}%, "
      f"<5%:{w5}/24, <10%:{w10}/24")
