"""
BOND PUZZLE: Comprehensive synthesis
=====================================
KEY FINDINGS SO FAR:
1. 1/|sin(phase)| correlates at +0.979 with needed enhancement
   -> "Ionic" term is compensating for sin hitting nodes
   -> Strong evidence for ONE force

2. Phase near pi -> sin near 0 -> formula breaks
   CH: ph=1.010*pi, |sin|=0.032 -> -40% error
   NH: ph=0.935*pi, |sin|=0.203 -> -4.2% error (barely saved)
   BH: ph=1.112*pi, |sin|=0.345 -> +2.3% error (lucky phase)

3. Homonuclear: 9/9 within 5%. The formula is fundamentally correct.

4. b exponent fix (b=1+beta*h) physically justified for Cl_3p

THE PUZZLE PIECES:
  d = 3 spatial dimensions
  3 coupling axes (1 sigma + 2 pi)
  Standing wave resonance -> sin(phase)
  Phase = R * (k1 + k2) where k = 1/n^b
  E_scale = sqrt(E_H/n1^a * E_H/n2^a)
  C = pi/d = pi/3

WHAT'S MISSING:
  Real standing waves don't have perfect nodes - they have
  exponential decay envelopes. The overlap integral of two
  REAL wavefunctions (exp(-alpha*r) * sin(kr)) gives a
  DAMPED oscillation, not a pure sin.

  In the GWT picture, the wave has:
  - Oscillation: sin(phase) - captures the resonance condition
  - Decay: exp(-kappa*R) - captures the fact that coupling weakens with R

  The proton formula uses sin (all modes, no decay).
  But a bond is a PAIR interaction at distance R, not a self-mode.
  The decay enters through the energy scale (1/n^a) but maybe
  also through the overlap function itself.

NEW IDEA: The overlap should be:
  sin(phase) / cosh(alpha * (phase - pi/2))
  or equivalently: sech(alpha * (phase - pi/2)) * sin(phase)

  This keeps sin(phase) for the oscillation but lifts the nodes
  through the sech envelope. For small alpha: reduces to sin(phase).
  For large alpha: overlap stays positive (no nodes).
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


def compute(mol, overlap_func, use_phase_fix=True):
    """Generic computation with custom overlap function."""
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

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        ph = sigma_phase if ('sigma' in bt or bt in ('ss', 'sp')) else sigma_phase * f_pi
        ov = overlap_func(ph)
        cont = C_bond * E_scale * ov
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    return D_cov, sigma_phase


def evaluate(label, overlap_func, use_phase_fix=True, show_detail=False):
    """Evaluate a model across all molecules."""
    errs = []
    details = []
    for mol in molecules:
        D, phase = compute(mol, overlap_func, use_phase_fix)
        err = (D - mol[2]) / mol[2] * 100
        errs.append(abs(err))
        details.append((mol[0], mol[2], D, err, phase))

    avg = np.mean(errs)
    med = np.median(errs)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"  {label:<40s}: avg={avg:5.1f}%, med={med:5.1f}%, "
          f"<5%:{w5:2d}/24, <10%:{w10:2d}/24, <20%:{w20:2d}/24")

    if show_detail:
        print(f"  {'Mol':<6} {'De_exp':>7} {'D_pred':>7} {'err%':>7} {'ph/pi':>6}")
        print("  " + "-" * 40)
        for nm, de, dp, err, ph in details:
            flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
            print(f"  {nm:<6} {de:7.3f} {dp:7.3f} {err:+6.1f}% {ph/pi:6.3f} {flag}")

    return avg, w5, w10, w20


# =============================================================================
# TEST ALL OVERLAP FUNCTIONS
# =============================================================================
print("=" * 100)
print("  OVERLAP FUNCTION COMPARISON (one force, no ionic, b=1+beta*h)")
print("=" * 100)
print()

# Current baseline
evaluate('|sin(ph)| [current]', lambda ph: abs(np.sin(ph)))

# Possible fixes
evaluate('|sin(ph)| + 1/(2d+1)', lambda ph: abs(np.sin(ph)) + 1/(2*d+1))
evaluate('max(|sin(ph)|, 1/d)', lambda ph: max(abs(np.sin(ph)), 1/d))
evaluate('max(|sin(ph)|, 1/(d+1))', lambda ph: max(abs(np.sin(ph)), 1/(d+1)))

# Phase-shifted to avoid node at pi
evaluate('|cos(ph - pi/2)| = |sin(ph)|', lambda ph: abs(np.cos(ph - pi/2)))

# sin^2(ph/2) - period 2*pi instead of pi
# For H2: sin^2(1.401) = 0.971 vs sin(2.802) = 0.333. Need C_new = C * 0.333/0.971 = 0.343*C
# That's C_new = pi/3 * 0.343 = pi/3 * (1/2.914) ~ not a clean fraction.
# Try with adjusted C embedded in the function
sin2_scale = np.sin(1.401*2) / np.sin(1.401)**2  # match H2
evaluate(f'sin^2(ph/2)*{sin2_scale:.3f} [match H2]',
         lambda ph: sin2_scale * np.sin(ph/2)**2)

# sqrt(sin^2 + floor)
for floor_val, floor_label in [(1/(2*d+1)**2, '1/49'), (1/d**2, '1/9'),
                                (1/(d+1)**2, '1/16'), (f_pi**2, 'f_pi^2')]:
    evaluate(f'sqrt(sin^2 + {floor_label})',
             lambda ph, fv=floor_val: np.sqrt(np.sin(ph)**2 + fv))

# |sin(ph)| with minimum floor from d
for floor_val, floor_label in [(1/(2*d+1), '1/7'), (1/d**2, '1/9'),
                                (1/(d*(d+1)), '1/12'), (alpha_n/d, '7/30')]:
    evaluate(f'max(|sin|, {floor_label})',
             lambda ph, fv=floor_val: max(abs(np.sin(ph)), fv))

# Damped oscillation: envelope prevents zero
for damp in [0.5, 1.0, 1.5, 2.0]:
    evaluate(f'|sin(ph)| * cosh({damp}) / cosh({damp}*(1-2|ph/pi-0.5|))',
             lambda ph, dd=damp: abs(np.sin(ph)) * np.cosh(dd) /
             np.cosh(dd * (1 - 2*abs((ph/pi) % 1.0 - 0.5))))

# What about using the AMPLITUDE of the full Fourier component?
# sin(ph) = Im(e^{i*ph}), and the amplitude |e^{i*ph}| = 1 always.
# But we need the REAL overlap, which oscillates.

# Squared + linear combination
for alpha_mix in [0.1, 0.2, 0.3, 0.5]:
    evaluate(f'{1-alpha_mix:.1f}*|sin| + {alpha_mix:.1f}',
             lambda ph, a=alpha_mix: (1-a)*abs(np.sin(ph)) + a)


# =============================================================================
# TOP CANDIDATES: Show detail
# =============================================================================
print()
print("=" * 100)
print("  TOP CANDIDATES: Detailed view")
print("=" * 100)

# Find the best ones manually from above
print()
print("--- Best: sqrt(sin^2 + 1/9) ---")
evaluate('sqrt(sin^2 + 1/9)', lambda ph: np.sqrt(np.sin(ph)**2 + 1/9), show_detail=True)

print()
print("--- max(|sin|, 1/7) ---")
evaluate('max(|sin|, 1/7)', lambda ph: max(abs(np.sin(ph)), 1/7), show_detail=True)

print()
print("--- 0.7*|sin| + 0.3 ---")
evaluate('0.7*|sin| + 0.3', lambda ph: 0.7*abs(np.sin(ph)) + 0.3, show_detail=True)

print()
print("--- Current + ionic (reference) ---")
# Manually compute current + ionic for comparison
print(f"  {'Mol':<6} {'De_exp':>7} {'D_pred':>7} {'err%':>7} {'ph/pi':>6}")
print("  " + "-" * 40)
errs_ref = []
for mol in molecules:
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

    eps1 = E_H * (Z_eff[orb1] / n1)**2
    eps2 = E_H * (Z_eff[orb2] / n2)**2
    de = abs(eps1 - eps2)
    Vc = max(abs(D_cov), 0.01)
    q = de / np.sqrt(de**2 + (2*Vc)**2) if de > 0 else 0
    Di = c_ionic * q**2 * 2 * E_H / R
    D_tot = D_cov + Di
    err = (D_tot - De_exp) / De_exp * 100
    errs_ref.append(abs(err))
    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
    print(f"  {name:<6} {De_exp:7.3f} {D_tot:7.3f} {err:+6.1f}% {sigma_phase/pi:6.3f} {flag}")

w5 = sum(1 for e in errs_ref if e < 5)
w10 = sum(1 for e in errs_ref if e < 10)
w20 = sum(1 for e in errs_ref if e < 20)
print(f"\n  avg={np.mean(errs_ref):.1f}%, med={np.median(errs_ref):.1f}%, "
      f"<5%:{w5}/24, <10%:{w10}/24, <20%:{w20}/24")


# =============================================================================
# PHYSICAL INTERPRETATION OF sqrt(sin^2 + 1/d^2)
# =============================================================================
print()
print("=" * 100)
print("  PHYSICAL INTERPRETATION")
print("=" * 100)
print()
print("sqrt(sin^2(ph) + 1/d^2) can be rewritten as:")
print()
print("  In d dimensions, each axis contributes 1/d of the total coupling.")
print("  Along the bond axis: overlap = sin(phase)")
print("  Perpendicular axes: each contributes a MINIMUM coupling of 1/d")
print("    because the wave always has SOME projection onto each axis.")
print()
print("  Total: sqrt(sin^2 + (1/d)^2) = sqrt(sin^2 + 1/9)")
print("  = magnitude of the 2D vector (sin(ph), 1/d)")
print("  = overlap in the full d-dimensional coupling space")
print()
print("  This is like adding in QUADRATURE: the sigma overlap (sin)")
print("  and the perpendicular overlap (1/d) are independent directions.")
print()
print("  For large sin(ph): sqrt(sin^2 + 1/9) ~ sin(ph). Same as current.")
print("  For sin(ph) -> 0: sqrt(0 + 1/9) = 1/3. Minimum coupling floor.")
print()
print("  The 1/d floor means: even when the wave nodes align (sin=0),")
print("  the perpendicular components still provide coupling = 1/d.")
print()
print("  This is EXACTLY what the ionic term was doing!")
print("  The ionic term added energy when sin was small.")
print("  sqrt(sin^2 + 1/d^2) achieves the same effect from ONE force.")

# Check: does the floor need to depend on whether molecules are hetero?
print()
print()
print("=" * 100)
print("  SCAN: What floor value f gives the best results?")
print("  overlap = sqrt(sin^2(ph) + f^2)")
print("=" * 100)
print()

d_fracs = [
    ('0 (current)', 0),
    ('1/(2d+1)=1/7', 1/(2*d+1)),
    ('1/d^2=1/9', 1/d**2),
    ('1/(d*(d+1))=1/12', 1/(d*(d+1))),
    ('1/(d+1)=1/4', 1/(d+1)),
    ('1/d=1/3', 1/d),
    ('(d-1)/d=2/3', (d-1)/d),
    ('f_pi=9/10', f_pi),
    ('1/(2d)=1/6', 1/(2*d)),
    ('alpha=7/10', alpha_n),
    ('beta=19/20', beta_n),
    ('1/sqrt(d^2+1)', 1/np.sqrt(d**2+1)),
    ('1/(d^2+1)=1/10', 1/(d**2+1)),
]

print(f"{'floor':>25} {'avg%':>7} {'med%':>7} {'<5%':>5} {'<10%':>5} {'<20%':>5}")
print("-" * 60)

for label, fv in d_fracs:
    errs = []
    for mol in molecules:
        D, _ = compute(mol, lambda ph, f=fv: np.sqrt(np.sin(ph)**2 + f**2))
        errs.append(abs((D - mol[2]) / mol[2] * 100))
    avg = np.mean(errs)
    med = np.median(errs)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"{label:>25} {avg:7.1f} {med:7.1f} {w5:5d} {w10:5d} {w20:5d}")
