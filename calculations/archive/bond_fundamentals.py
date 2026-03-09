"""
BOND FORMULA FROM FUNDAMENTALS
================================
List every physical ingredient, then build the formula step by step.

FUNDAMENTALS:
=============
1. d=3 spatial dimensions
2. Each dimension is an independent coupling axis
3. Two standing waves overlap to form a bond
4. The overlap depends on:
   a. Wavevectors k1, k2 (de Broglie: spatial frequency of each wave)
   b. Separation R (distance between wave centers)
   c. Energy scale (amplitude of each wave)
   d. Relative phase (constructive vs destructive)

DERIVED CONSTANTS (all from d=3):
  C = pi/d = pi/3                    (single-mode coupling)
  f_pi = d^2/(d^2+1) = 9/10          (pi vs sigma phase coupling)
  alpha = (d^2-d+1)/(d^2+1) = 7/10   (energy node correction)
  beta = (2d^2+1)/(2(d^2+1)) = 19/20 (phase node correction)
  f_anti = 2d/(2d-1) = 6/5           (antibonding enhancement)

KEY QUESTION: What is the correct overlap integral for two waves
with DIFFERENT wavevectors?

For IDENTICAL waves (k1 = k2 = k):
  Overlap ~ sin(kR + kR) = sin(2kR)  [phase from both centers]
  This WORKS perfectly for all 9 homonuclear molecules.

For DIFFERENT waves (k1 != k2):
  Current: sin((k1+k2)*R)  [sum of phases]
  Problem: this can hit sin=0 at the wrong R

  The physics: wave 1 has wavelength lambda1 = 2*pi/k1
                wave 2 has wavelength lambda2 = 2*pi/k2
  The OVERLAP depends on how well they "fit" together.
  Two waves with different wavelengths create BEATS.
  The beating pattern has:
    - Fast oscillation at average frequency: (k1+k2)/2
    - Slow envelope at difference frequency: (k1-k2)/2

  Overlap amplitude = cos((k1-k2)*R/2) * sin((k1+k2)*R/2)
                     = [sin(k1*R) + sin(k2*R)] / 2

  Wait - this is just the AVERAGE of the individual overlaps!
  For homonuclear: [sin(kR) + sin(kR)] / 2 = sin(kR)
  Current formula: sin(2kR) = 2*sin(kR)*cos(kR)
  These are DIFFERENT functions!

  So if the correct form is sin(k1*R) + sin(k2*R),
  then for homonuclear we get 2*sin(kR), not sin(2kR).
  Since sin(2kR) works, the overlap is NOT a simple sum.

  The current formula sin((k1+k2)*R) = sin(k1*R+k2*R) uses the
  angle addition: sin(A+B) = sinA*cosB + cosA*sinB.
  This represents a PRODUCT-like coupling, not a sum.

  Physical interpretation:
    sin(k1*R) = wave 1's amplitude at wave 2's center
    cos(k2*R) = wave 2's phase alignment at wave 2's center
    And vice versa.
    Total = sin(k1*R)*cos(k2*R) + cos(k1*R)*sin(k2*R)

  For homonuclear: 2*sin(kR)*cos(kR) = sin(2kR). Checks out.

  So the current formula IS physically motivated. But it assumes
  the phase from each wave adds LINEARLY to a total phase.
  This is correct for identical waves. For different waves,
  the question is whether the phases add equally.

  What if the phase contribution is WEIGHTED by the coupling strength?
  Heavier atom contributes less phase because its wave decays faster.

  Weighted phase: w1*k1*R + w2*k2*R where w1 + w2 = 2
  For homonuclear: w1 = w2 = 1, so phase = 2kR. Same.

  What should the weights be? Options:
    - Equal: w = (1, 1) -> phase = (k1+k2)*R [current]
    - By E: w proportional to sqrt(E_i) [heavier wave weighted more]
    - By k: w proportional to k_i [faster wave weighted more]
    - Inverse k: w proportional to 1/k_i [slower wave weighted more]

ENERGY CONSIDERATIONS:
  Current: E_scale = sqrt(E_H/n1^a * E_H/n2^a) = E_H / sqrt(n1^a * n2^a)
  This gives ALL 2p orbitals (B, C, N, O, F) the same energy = E_H/4 = 3.40 eV
  But actual orbital energies range from 19.9 (B) to 88.5 (F) eV!

  For homonuclear, this doesn't matter - the E_scale * sin(phase) product
  works because E_H/4 * sin(phase) happens to give the right numbers.
  The Z_eff information is encoded in the phase (via observed R) and the
  bond list (which MOs are filled).

  For heteronuclear, E_scale = sqrt(3.40 * 3.40) = 3.40 for ALL 2p-2p pairs.
  But CO (C=3.40, O=3.40) and BF (B=3.40, F=3.40) need very different energies.
  The formula can't distinguish them through E_scale or phase alone.

THREE-AXIS COUPLING:
  In d=3, a sigma bond couples through all 3 axes:
    - Along bond (sigma): phase_sigma = (k1+k2)*R
    - Perpendicular 1 (pi): phase_pi = f_pi * (k1+k2)*R
    - Perpendicular 2 (pi): phase_pi = f_pi * (k1+k2)*R

  If each axis contributes independently:
    D_sigma_bond = C/d * E_scale * [|sin(ph_s)| + 2*|sin(ph_pi)|]

  For homonuclear, this changes the functional form and breaks the
  perfect results. So the 3 axes don't couple independently for a
  single sigma bond mode.

  Instead: the 3 axes determine the NUMBER of bond types:
    - 1 sigma direction
    - 2 pi directions
    - Total: d bond modes in d dimensions

WHAT'S ACTUALLY MISSING:
  The formula works for homonuclear because E_scale and sin(phase) are
  calibrated together through the SAME orbital. For heteronuclear, the
  two orbitals have different spatial extents and energies, and the simple
  sin((k1+k2)*R) doesn't capture the asymmetry correctly.

  The "ionic" term was adding: c * q^2 * 2*E_H/R
  This is proportional to: (energy difference)^2 / R
  It adds energy when the orbitals have DIFFERENT energies.

  From a one-force perspective, this means the overlap of two DIFFERENT
  standing waves produces more bonding than the simple sin formula predicts.
  Why? Because the electron density SHIFTS toward the more negative atom,
  creating an enhanced overlap in the region between the atoms.

  This is still wave resonance - just with an asymmetric amplitude profile.
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
def bohr_radius(orb): return get_n(orb)**2 / Z_eff[orb]


# =============================================================================
# APPROACH: Phase-weighted overlap with energy enhancement
# =============================================================================
#
# The key insight: sin((k1+k2)*R) = sin(k1*R)*cos(k2*R) + cos(k1*R)*sin(k2*R)
#
# Each term represents one atom "seeing" the other's wave.
# Term 1: wave 1 oscillation at R, modulated by wave 2's phase at R
# Term 2: wave 2 oscillation at R, modulated by wave 1's phase at R
#
# For asymmetric pairs, these two terms are NOT equal.
# The AMPLITUDE of each term is weighted by E_i:
#   Term 1 ~ E1 * sin(k1*R) * cos(k2*R) / E_ref
#   Term 2 ~ E2 * cos(k1*R) * sin(k2*R) / E_ref
#
# For homonuclear: E1 = E2 = E, terms add to E*sin(2kR)/E_ref = sin(2kR)
#   with E_ref = E. This gives E_scale = E, sin(2kR). Same as current.
#
# For heteronuclear:
#   D = C * [E1*sin(k1R)*cos(k2R) + E2*cos(k1R)*sin(k2R)] / E_ref
#
# What's E_ref? For homonuclear, E_ref = E gives sin(2kR).
# For the general case, E_ref = sqrt(E1*E2) preserves the geometric mean.
# Then: D = C * sqrt(E1*E2) * [sqrt(E1/E2)*sin(k1R)cos(k2R) + sqrt(E2/E1)*cos(k1R)sin(k2R)]
# Hmm, this doesn't simplify to sin((k1+k2)R) for homonuclear with E1=E2.
# Actually: sqrt(E1/E2) = 1 when E1=E2, so we get sin(k1R)cos(k2R) + cos(k1R)sin(k2R) = sin((k1+k2)R). YES!
#
# So the ENERGY-WEIGHTED overlap is:
#   O(R) = sqrt(E1/E2) * sin(k1*R) * cos(k2*R) + sqrt(E2/E1) * cos(k1*R) * sin(k2*R)
#
# For homonuclear: sin((k1+k2)*R) = sin(2kR). Same as current.
# For heteronuclear: the heavier atom's term gets enhanced.
#
# Let's test this!

print("=" * 100)
print("  ENERGY-WEIGHTED OVERLAP MODEL")
print("  O(R) = sqrt(E1/E2)*sin(k1R)*cos(k2R) + sqrt(E2/E1)*cos(k1R)*sin(k2R)")
print("  For homonuclear: reduces to sin((k1+k2)*R)")
print("  For heteronuclear: weights each atom's contribution by energy ratio")
print("=" * 100)
print()

# What energies to use? Options:
# A. E_H/n^a (current formula energy)
# B. E_orb = E_H*(Z/n)^2 (actual orbital energy)
# C. E_H/n^2 (unscreened hydrogen-like)

def compute_weighted(mol, E_type='formula', use_phase_fix=True):
    """Energy-weighted overlap model."""
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

    # Energy scale (geometric mean, as before)
    E_formula_1 = E_H / n1**a1
    E_formula_2 = E_H / n2**a2
    E_scale = np.sqrt(E_formula_1 * E_formula_2)

    k1 = 1.0 / n1**b1
    k2 = 1.0 / n2**b2

    # Energy ratio for weighting
    if E_type == 'formula':
        Er1, Er2 = E_formula_1, E_formula_2
    elif E_type == 'orbital':
        Er1, Er2 = orb_energy(orb1), orb_energy(orb2)
    elif E_type == 'hydrogen':
        Er1, Er2 = E_H / n1**2, E_H / n2**2
    else:
        Er1, Er2 = E_formula_1, E_formula_2

    # Weight factors
    if Er2 > 0.001 and Er1 > 0.001:
        w1 = np.sqrt(Er1 / Er2)  # weight for atom 1's sin term
        w2 = np.sqrt(Er2 / Er1)  # weight for atom 2's sin term
    else:
        w1, w2 = 1.0, 1.0

    npb = sum(c for bt, c in bonds if 'pi' in bt and 'anti' not in bt)
    npa = sum(c for bt, c in bonds if 'pi' in bt and 'anti' in bt)
    ifa = (npa >= npb) if npb > 0 else True

    D_cov = 0
    for bt, cnt in bonds:
        if 'sigma' in bt or bt in ('ss', 'sp'):
            ph1 = k1 * R
            ph2 = k2 * R
        else:
            ph1 = k1 * R * f_pi
            ph2 = k2 * R * f_pi

        # Energy-weighted overlap
        overlap = abs(w1 * np.sin(ph1) * np.cos(ph2) + w2 * np.cos(ph1) * np.sin(ph2))

        cont = C_bond * E_scale * overlap
        if 'anti' in bt:
            fa = 1.0 if ('sigma' in bt or ifa) else f_anti
            D_cov -= cnt * fa * cont
        else:
            D_cov += cnt * cont

    return D_cov, k1, k2


# Test with different energy types
for E_type in ['formula', 'orbital', 'hydrogen']:
    for pfix in [True, False]:
        pfix_label = 'fix' if pfix else 'old'
        errs = []
        for mol in molecules:
            D, _, _ = compute_weighted(mol, E_type=E_type, use_phase_fix=pfix)
            err = abs((D - mol[2]) / mol[2] * 100)
            errs.append(err)

        w5 = sum(1 for e in errs if e < 5)
        w10 = sum(1 for e in errs if e < 10)
        w20 = sum(1 for e in errs if e < 20)
        print(f"  E={E_type:>8s}, b={pfix_label}: avg={np.mean(errs):.1f}%, "
              f"med={np.median(errs):.1f}%, <5%:{w5}/24, <10%:{w10}/24, <20%:{w20}/24")


# Detail for the best one
print()
print("=" * 100)
print("  DETAIL: Energy-weighted overlap with orbital energies")
print("=" * 100)
print()

best_E = 'orbital'
best_pfix = True

print(f"{'Mol':<6} {'De_exp':>7} {'D_pred':>7} {'err%':>7} {'w1':>6} {'w2':>6} {'k1R':>6} {'k2R':>6}")
print("-" * 60)

for mol in molecules:
    name = mol[0]; De_exp = mol[2]
    D, k1, k2 = compute_weighted(mol, E_type=best_E, use_phase_fix=best_pfix)
    err = (D - De_exp) / De_exp * 100

    e1 = orb_energy(mol[4])
    e2 = orb_energy(mol[5])
    w1 = np.sqrt(e1/e2) if e2 > 0.001 else 1
    w2 = np.sqrt(e2/e1) if e1 > 0.001 else 1

    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
    print(f"{name:<6} {De_exp:7.3f} {D:7.3f} {err:+6.1f}% {w1:6.3f} {w2:6.3f} "
          f"{k1*mol[1]:6.3f} {k2*mol[1]:6.3f} {flag}")


# =============================================================================
# WHAT IF: The energy weighting goes into E_scale instead?
# =============================================================================
print()
print("=" * 100)
print("  ALTERNATIVE: E_scale uses orbital energies, phase unchanged")
print("  D = C * sqrt(E_orb1 * E_orb2) * |sin((k1+k2)*R)| / E_correction")
print("=" * 100)
print()

# For homonuclear 2p molecules:
# Current: sqrt(3.40 * 3.40) = 3.40
# Orbital: sqrt(E1 * E2)
#   N2: sqrt(50 * 50) = 50.0
# Need to rescale so homonuclear still works.
# Correction = E_orb / E_formula per atom
# For homonuclear: correction = (E_orb/E_form)^2 in the geometric mean
# So divide by E_orb/E_form: sqrt(E_orb1*E_orb2) / sqrt(E_form1*E_form2)
# which gives back E_form. Circular!

# Unless we use a DIFFERENT functional form with E_orb.
# What if: D = C * sqrt(E_orb1 * E_orb2) * |sin(phase)| / n_eff^p
# where n_eff is some quantum number combination and p adjusts the scale?

# For H2: E_orb = 13.6, D = pi/3 * 13.6 * sin(2.802) = 14.26 * 0.331 = 4.72
# For N2: E_orb = 50.0, without correction:
#   D_sigma = pi/3 * 50.0 * sin(2.074) = 52.36 * 0.879 = 46.0
#   D_pi = 2 * pi/3 * 50.0 * sin(0.9*2.074) = 2 * 52.36 * 0.958 = 100.3
#   Total = 146.3 vs 9.76 needed. Off by 15x!
# So E_orb is way too big for 2p orbitals. The formula E_H/n^a deliberately
# removes Z_eff to get the right scale. Using E_orb breaks everything.

# The key insight: E_H/n^a is the COUPLING ENERGY (how strongly the wave
# couples to external perturbations), not the ORBITAL ENERGY (total binding).
# The coupling goes as 1/n^a because higher n means more diffuse wave,
# weaker coupling. Z_eff doesn't enter because the coupling depends on
# the wave PATTERN (universally 1/n^a in d dimensions) not the depth
# of the potential well.

# So for energy weighting, we should use E_formula (not E_orb) as the
# overlap weight. But for homonuclear, E_formula_1 = E_formula_2, so
# w1 = w2 = 1, and the formula reduces to the current one.
# This means energy weighting with E_formula does NOTHING for heteronuclear
# 2p-2p pairs (which ALL have E = 3.40)!

print("INSIGHT: E_H/n^a gives the same energy to all 2p orbitals.")
print("This is BY DESIGN: the coupling energy depends on wave pattern,")
print("not the depth of the potential well.")
print()
print("For heteronuclear 2p-2p pairs (CO, NO, BF, CN):")
print("  E_scale = sqrt(3.40 * 3.40) = 3.40 for ALL of them")
print("  Phase also identical: k1 = k2 = 0.5 for all 2p")
print("  Only difference: bond list (which MOs are filled) and R")
print()
print("This means the formula CAN'T distinguish CO from BF through")
print("energy or wavevector - only through R and bond list.")
print("And it WORKS for these (BF: +7.7% with harm, CN: +21.1%)!")
print()
print("The REAL problem molecules involve DIFFERENT n values:")
print("  XH bonds (X is 2p/2s/3s, H is 1s)")
print("  LiF (2s + 2p)")
print("  NaCl (3s + 3p)")
print()

# =============================================================================
# FRESH APPROACH: What if k is the ONLY thing that needs fixing?
# =============================================================================
print("=" * 100)
print("  FRESH: What if we keep everything EXCEPT the wavevector definition?")
print("  For 2p orbitals (no nodes): k = 1/n = 0.5 [current, works]")
print("  For s orbitals with nodes: k = 1/n^1.95 [current, works for homonuclear]")
print("  For p orbitals with nodes: k = 1/n^0.05 -> BROKEN")
print("=" * 100)
print()

# The ONLY orbitals affected by the phase fix are:
# Cl_3p: n=3, l=1, h=1 -> b_old=0.05, b_fix=1.95
# Any future 3d, 4p, 4d, 4f etc would also be affected

# For our 24 molecules, only HCl and NaCl have Cl_3p.
# Cl2 also has Cl_3p but it's homonuclear (not affected by phase weighting).

# So the "phase fix" is really just about Cl_3p.
# It helps HCl (84% -> 2.3%) but slightly hurts NaCl (42% -> 47%).

# The OTHER failing molecules (LiH, CH, BH, NaH, LiF) have
# orbitals with b=1.0 (2p, no nodes) or b=1.95 (2s/3s, with nodes).
# These are NOT affected by the phase fix.

# For these, the problem is the COMBINATION of very different k's:
#   LiH: k(Li_2s) = 0.259, k(H_1s) = 1.000 -> ratio 3.86:1
#   BH:  k(B_2p) = 0.500, k(H_1s) = 1.000  -> ratio 2.00:1
#   CH:  k(C_2p) = 0.500, k(H_1s) = 1.000  -> ratio 2.00:1
#   NH:  k(N_2p) = 0.500, k(H_1s) = 1.000  -> ratio 2.00:1
#   NaH: k(Na_3s)= 0.117, k(H_1s) = 1.000  -> ratio 8.52:1

# BH, CH, NH all have k_ratio = 2. But their errors are very different:
#   BH: +2.3% (phase=1.112*pi, near node!)
#   CH: -40.4% (phase=1.010*pi, right at node!)
#   NH: -4.2% (phase=0.935*pi, OK)

# The pattern: phase near pi -> sin near 0 -> error blows up.
# BH is saved because phase=1.112*pi gives |sin|=0.345, still nonzero.
# CH at 1.010*pi gives |sin|=0.031, basically zero.

# So the REAL issue isn't k_ratio - it's phase proximity to pi.
# And phase = R*(k1+k2) depends on R (bond length).

# What determines which molecules have phase near pi?
print("Phase proximity to nearest multiple of pi:")
print()
print(f"{'Mol':<6} {'ph/pi':>6} {'dist_pi':>8} {'|sin|':>6} {'err%':>7} {'R':>5} {'k1+k2':>6}")
print("-" * 55)

for mol in molecules:
    name, R, De_exp, bonds, orb1, orb2 = mol
    n1 = get_n(orb1); h1 = has_nodes(orb1)
    n2 = get_n(orb2); h2 = has_nodes(orb2)
    b1 = 1 + beta_n * h1
    b2 = 1 + beta_n * h2
    k1 = 1.0 / n1**b1
    k2 = 1.0 / n2**b2
    phase = R * (k1 + k2)

    # Distance to nearest multiple of pi
    ph_mod = (phase / pi) % 1.0
    dist = min(ph_mod, 1 - ph_mod)

    sin_val = abs(np.sin(phase))

    D_cov = compute_weighted(mol, E_type='formula', use_phase_fix=True)[0]
    err = abs((D_cov - De_exp) / De_exp * 100)

    flag = ' <<<' if dist < 0.1 else ''
    print(f"{name:<6} {phase/pi:6.3f} {dist:8.3f} {sin_val:6.3f} {err:6.1f}% "
          f"{R:5.3f} {k1+k2:6.3f}{flag}")

print()
print("<<< = phase within 0.1*pi of a node (sin near zero)")
print()
print("CONCLUSION: The formula fails when R*(k1+k2) happens to land near")
print("a multiple of pi. This is a RESONANCE CONDITION issue, not a")
print("missing force. The sin function has nodes, and some molecules")
print("happen to sit near those nodes.")
print()
print("POSSIBLE FIX: Use a function that doesn't have nodes?")
print("  - |sin(phase)| has nodes at 0, pi, 2pi, ...")
print("  - sin^2(phase/2) has nodes at 0, 2pi, 4pi, ... (half as many)")
print("  - (1 - cos(phase))/2 = sin^2(phase/2), same thing")
print()

# Test sin^2(phase/2): this has period 2*pi instead of pi
# For H2: phase = 2.802, sin^2(1.401) = 0.9713
# Current: sin(2.802) = 0.3331
# Need to adjust C to match:
# C_new * sin^2(1.401) = C_old * sin(2.802)
# C_new = C_old * 0.3331 / 0.9713 = 0.343 * C_old
# That's NOT a clean d-fraction. Bad sign.

# Instead, try sin(phase) * |sin(phase)| = sin^2(phase) * sign(sin)
# Or: 2*sin(phase/2)*cos(phase/2) = sin(phase). Same thing. No help.

# What about taking the ENVELOPE of sin?
# The envelope of the beat pattern cos((k1-k2)R/2) * sin((k1+k2)R/2) is
# |cos((k1-k2)R/2)|.  For homonuclear: cos(0) = 1.
# For heteronuclear: cos((k1-k2)*R/2) is ALWAYS <= 1.
# So the envelope actually REDUCES the heteronuclear overlap. Wrong direction.

# Let me try something else: what if the overlap is
# sin((k1+k2)*R) but with the EFFECTIVE k corrected for the asymmetry?
#
# The beat pattern sin(k1R)*cos(k2R) + cos(k1R)*sin(k2R) = sin((k1+k2)R)
# But there's a SECOND solution: the overlap could also depend on |k1-k2|R
# via the beat envelope.
#
# Total overlap = sin((k1+k2)*R) * [1 + f(|k1-k2|/k_avg)]
# For homonuclear: f(0) = 0, so just sin(2kR). Good.
# For heteronuclear: f > 0 would enhance the overlap.
#
# What function f? The beat creates an additional resonance at the
# difference frequency. The energy associated with this is:
#   D_beat ~ E_scale * cos((k1-k2)*R/2) - but that's <= 1, not > 1.

# Actually, let me think about this differently.
# The CURRENT formula with ionic term gives:
#   D = C * E_scale * sin(phase) + c_ionic * q^2 * 2*E_H/R
#
# The ionic term is q^2/R proportional. For a given R,
# q depends on the ENERGY DIFFERENCE de = |eps1 - eps2|.
# de is essentially the electronegativity difference.
#
# What if this is really: the SECOND HARMONIC of the standing wave?
# A standing wave has modes at k, 2k, 3k, ...
# The fundamental (k) gives sin(kR) - the main resonance
# The second harmonic (2k) gives sin(2kR) - a correction
#
# For two different atoms, there might be CROSS-HARMONICS:
# sin(k1R) * sin(k2R) = [cos((k1-k2)R) - cos((k1+k2)R)] / 2
# This product has frequencies at k1+k2 and k1-k2.
# The k1+k2 term is what we already have (sin becomes cos, shifted).
# The k1-k2 term is NEW! It represents the BEAT frequency.
#
# The beat frequency contribution:
#   D_beat ~ C' * E_scale * cos((k1-k2)*R)
# For homonuclear: cos(0) = 1, so this adds a constant.
# But that constant would change the homonuclear results...
# Unless C' is proportional to |k1-k2|/(k1+k2) or similar,
# which vanishes for homonuclear.

print("=" * 100)
print("  BEAT FREQUENCY MODEL: D_beat ~ C * E * (dk/k_avg)^p * cos(dk*R)")
print("  dk = |k1-k2|, k_avg = (k1+k2)/2")
print("  Vanishes for homonuclear (dk=0). Adds energy for heteronuclear.")
print("=" * 100)
print()

def compute_beat(mol, C_beat, p_dk, use_phase_fix=True):
    """Current formula + beat frequency correction."""
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
    dk = abs(k1 - k2)
    k_avg = (k1 + k2) / 2

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

    # Beat correction (only for heteronuclear)
    dk_ratio = dk / k_avg if k_avg > 0.001 else 0
    D_beat = C_beat * E_scale * dk_ratio**p_dk * abs(np.cos(dk * R))

    return D_cov + D_beat, D_cov, D_beat


# Scan C_beat and p_dk
print("Scanning beat correction parameters:")
print(f"{'C_beat':>8} {'p_dk':>5} {'avg%':>7} {'med%':>7} {'<5%':>5} {'<10%':>5} {'<20%':>5}")
print("-" * 50)

best_beat = (999, 0, 0)
for C_b_frac, C_b_label in [
    (pi/d, 'pi/3'),
    (pi/(d+1), 'pi/4'),
    (1.0, '1'),
    (2.0, '2'),
    (pi, 'pi'),
    (2*pi/d, '2pi/3'),
    (pi/d**2, 'pi/9'),
    (1.0/d, '1/3'),
    (2.0/d, '2/3'),
]:
    for p in [0.5, 1.0, 1.5, 2.0]:
        errs = []
        for mol in molecules:
            D = compute_beat(mol, C_b_frac, p)[0]
            err = abs((D - mol[2]) / mol[2] * 100)
            errs.append(err)
        avg = np.mean(errs)
        w5 = sum(1 for e in errs if e < 5)
        w10 = sum(1 for e in errs if e < 10)
        w20 = sum(1 for e in errs if e < 20)
        if avg < best_beat[0]:
            best_beat = (avg, C_b_frac, p, C_b_label, w5, w10, w20)
        if w10 >= 15 or avg < 15:
            print(f"{C_b_label:>8} {p:5.1f} {avg:7.1f} {np.median(errs):7.1f} {w5:5d} {w10:5d} {w20:5d}")

print()
avg, C_b, p, C_b_label, w5, w10, w20 = best_beat
print(f"Best: C_beat={C_b_label}, p={p:.1f}, avg={avg:.1f}%, <5%:{w5}, <10%:{w10}, <20%:{w20}")

# Show detail
print()
print(f"{'Mol':<6} {'De_exp':>7} {'D_cov':>7} {'D_beat':>7} {'D_tot':>7} {'err%':>7} {'dk/k':>6}")
print("-" * 60)

for mol in molecules:
    D_tot, D_cov, D_beat = compute_beat(mol, C_b, p)
    err = (D_tot - mol[2]) / mol[2] * 100
    n1 = get_n(mol[4]); h1 = has_nodes(mol[4])
    n2 = get_n(mol[5]); h2 = has_nodes(mol[5])
    b1 = 1 + beta_n * h1
    b2 = 1 + beta_n * h2
    k1 = 1.0 / n1**b1
    k2 = 1.0 / n2**b2
    dk = abs(k1-k2); k_avg = (k1+k2)/2
    dk_ratio = dk/k_avg if k_avg > 0.001 else 0
    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
    print(f"{mol[0]:<6} {mol[2]:7.3f} {D_cov:7.3f} {D_beat:7.3f} {D_tot:7.3f} "
          f"{err:+6.1f}% {dk_ratio:6.3f} {flag}")
