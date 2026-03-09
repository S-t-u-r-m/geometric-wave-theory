"""
SECOND-ORDER PERTURBATION DERIVATION
=====================================
Can the "ionic" term be derived as a second-order correction
to the standing wave overlap?

QUANTUM MECHANICS REFRESHER (for non-specialists):
  When you have two quantum states with energies E1 and E2,
  and a coupling V between them, the EXACT solution is:

  E_bond = -V + (E1-E2)^2 / (4V)   [to second order]

  First order:  -V = the direct coupling = our sin(phase) term
  Second order: (dE)^2 / (4V) = correction for asymmetry

  The second-order term is ALWAYS positive (adds bonding energy)
  and grows with the energy difference between the orbitals.

  This is the SAME physics as "ionic bonding" - charge shifts
  toward the more electronegative atom, gaining Coulomb energy.
  But it's not a SEPARATE force - it's a correction to the
  wave resonance when the two waves are unequal.

IN GWT TERMS:
  First order: D_1 = C * E_scale * |sin(phase)|
  Second order: D_2 = C_2 * (delta_E)^2 / |D_1| * geometric_factor

  Where delta_E = |eps_1 - eps_2| = orbital energy difference
  And geometric_factor accounts for the 1/R Coulomb scaling.

  The 2-level model gives:
    q = delta_E / sqrt(delta_E^2 + 4*V^2)
  where V = D_cov (the covalent coupling).

  Then D_ionic = c * q^2 * 2*E_H/R

  Can we derive c = 1/(2d+1) = 1/7 from the perturbation expansion?

DERIVATION ATTEMPT:
  In d dimensions, the coupling tensor has d^2+1 modes.
  The first-order coupling uses 1 mode (sigma) or d^2 modes (pi).
  The second-order correction involves the REMAINING modes.

  For a sigma bond:
    First order uses 1 out of d^2+1 modes
    Second order uses the remaining d^2 modes
    Ratio: d^2 / (d^2+1) = 9/10 = f_pi. Interesting!

  But c_ionic = 1/7 = 1/(2d+1), not 9/10.

  Alternative: The Coulomb coupling c = 1/(2d+1) might come from:
    In d dimensions, the charge density couples through 2d+1 multipoles
    (monopole + 2d dipole directions). The monopole (q=0 mode) gives 1/(2d+1).

  Or: Consider the standing wave in d dimensions.
    The wave has 2d faces (for a d-cube) plus 1 interior = 2d+1 regions.
    The charge transfer corresponds to moving charge from one face to another.
    Each face has weight 1/(2d+1), so the Coulomb coupling is 1/(2d+1).

  Let's test this numerically.
"""
import numpy as np

pi = np.pi
E_H = 13.6057
d = 3

# Derived constants
C_bond = pi / d
f_pi = d**2 / (d**2 + 1)
alpha_n = 1 - f_pi / d
beta_n = (1 + f_pi) / 2
f_anti = 2*d / (2*d - 1)
c_ionic = 1.0 / (2*d + 1)  # = 1/7

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


# =============================================================================
# THE 2-LEVEL MODEL: Exact solution vs perturbation expansion
# =============================================================================
print("=" * 90)
print("  2-LEVEL QUANTUM SYSTEM: Exact vs Perturbation Expansion")
print("=" * 90)
print()
print("Two states |1> and |2> with energies eps_1, eps_2 and coupling V:")
print()
print("  H = | eps_1   V  |")
print("      |  V    eps_2 |")
print()
print("Exact eigenvalues:")
print("  E_+/- = (eps_1+eps_2)/2 +/- sqrt((eps_1-eps_2)^2/4 + V^2)")
print()
print("Bonding energy = E_- - (eps_1+eps_2)/2 = -sqrt(dE^2/4 + V^2)")
print()
print("Perturbation expansion (V >> dE, covalent limit):")
print("  E_bond = -V * sqrt(1 + (dE/2V)^2)")
print("         = -V - dE^2/(8V) + ...")
print()
print("  First order:  -V (our covalent term)")
print("  Second order: -dE^2/(8V) (the correction)")
print()
print("  The correction is: dE^2/(8V)")
print("  With dE = |eps_1 - eps_2| and V = D_cov")
print()

# Compare this perturbation form with our ionic term
print("Our ionic term: c_ionic * q^2 * 2*E_H/R")
print("  where q = dE/sqrt(dE^2 + 4*V^2)")
print()
print("Perturbation form: dE^2 / (8*V)")
print()
print("These are DIFFERENT functional forms!")
print("  Ionic: proportional to q^2/R ~ dE^2/(dE^2+4V^2) / R")
print("  Perturb: proportional to dE^2/V")
print()
print("Which one fits better? Let's test both.")
print()


def compute_cov(mol):
    """Compute covalent-only energy."""
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

    eps1 = orb_energy(orb1)
    eps2 = orb_energy(orb2)
    dE = abs(eps1 - eps2)

    return D_cov, dE, R, sigma_phase


# =============================================================================
# TEST: Perturbation correction vs ionic term
# =============================================================================
print("=" * 90)
print("  Model A: current ionic   D_ion = c * q^2 * 2*E_H/R")
print("  Model B: perturbation    D_pt  = c * dE^2 / (8*V)")
print("  Model C: exact 2-level   D_ex  = sqrt(V^2 + dE^2/4) - V")
print("=" * 90)
print()

# For each heteronuclear molecule, compute what c gives the right answer
print(f"{'Mol':<6} {'De_exp':>7} {'D_cov':>7} {'D_need':>7} {'dE':>7} "
      f"{'D_ion':>7} {'D_pert':>7} {'D_exact':>7}")
print("-" * 70)

for mol in molecules:
    name, R, De_exp = mol[0], mol[1], mol[2]
    orb1, orb2 = mol[4], mol[5]
    D_cov, dE, R_val, phase = compute_cov(mol)

    D_need = De_exp - D_cov
    V = max(abs(D_cov), 0.01)

    # Model A: current ionic formula
    q = dE / np.sqrt(dE**2 + (2*V)**2) if dE > 0 else 0
    D_ion = c_ionic * q**2 * 2 * E_H / R_val

    # Model B: perturbation correction dE^2 / (8V)
    # Need a coefficient to match. Try c_ionic for now.
    D_pert_raw = dE**2 / (8 * V) if V > 0.01 else 0
    # Scale: raw perturbation gives huge values for large dE
    # Need to normalize somehow

    # Model C: exact 2-level solution
    D_exact = np.sqrt(V**2 + dE**2/4) - V

    if orb1 != orb2:
        print(f"{name:<6} {De_exp:7.3f} {D_cov:7.3f} {D_need:7.3f} {dE:7.2f} "
              f"{D_ion:7.3f} {D_pert_raw:7.3f} {D_exact:7.3f}")


# =============================================================================
# THE EXACT 2-LEVEL SOLUTION
# =============================================================================
print()
print("=" * 90)
print("  EXACT 2-LEVEL: D_total = sqrt(V^2 + dE^2/4) where V = D_cov, dE = |eps1-eps2|")
print("  This is the FULL quantum mechanical solution, no approximation")
print("=" * 90)
print()

# But the exact 2-level uses V = D_cov and dE = orbital energy difference.
# The orbital energies use Z_eff (observed). Can we use E_formula instead?

# Test: D_total = sqrt(D_cov^2 + c * dE^2) for various c and dE definitions
print("Scanning 2-level formula: D = sqrt(D_cov^2 + c * dE^2)")
print()

for de_type, de_label in [('orbital', 'dE = |E_orb1 - E_orb2|'),
                           ('formula', 'dE = |E_H/n1^a - E_H/n2^a|')]:
    print(f"--- {de_label} ---")
    for c_val, c_label in [
        (1/4, '1/4 (exact 2-level)'),
        (1/(2*d+1), '1/7'),
        (c_ionic**2, '(1/7)^2'),
        (1/(d**2), '1/9'),
        (1/(d*(d+1)), '1/12'),
        (1/(2*d+1)**2, '1/49'),
        (1/(d**2+1), '1/10'),
    ]:
        errs = []
        for mol in molecules:
            name, R, De_exp = mol[0], mol[1], mol[2]
            orb1, orb2 = mol[4], mol[5]
            D_cov, dE_orb, R_val, phase = compute_cov(mol)

            if de_type == 'orbital':
                dE = dE_orb
            else:
                n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
                n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)
                a1 = 2 + (1 - 2*l1) * alpha_n * h1
                a2 = 2 + (1 - 2*l2) * alpha_n * h2
                dE = abs(E_H / n1**a1 - E_H / n2**a2)

            D_tot = np.sqrt(D_cov**2 + c_val * dE**2)
            # Keep the sign of D_cov
            if D_cov < 0:
                D_tot = -D_tot
            err = abs((D_tot - De_exp) / De_exp * 100)
            errs.append(err)

        avg = np.mean(errs)
        w5 = sum(1 for e in errs if e < 5)
        w10 = sum(1 for e in errs if e < 10)
        w20 = sum(1 for e in errs if e < 20)
        print(f"  c={c_val:.6f} ({c_label:>20s}): avg={avg:5.1f}%, "
              f"<5%:{w5:2d}/24, <10%:{w10:2d}/24, <20%:{w20:2d}/24")
    print()


# =============================================================================
# BEST MODEL: D = sqrt(D_cov^2 + c * dE_formula^2)
# =============================================================================
# dE_formula = |E_H/n1^a - E_H/n2^a| only depends on quantum numbers, not Z_eff
# This is fully self-consistent: NO observed values as input!

print("=" * 90)
print("  SELF-CONSISTENT MODEL: D = sqrt(D_cov^2 + c * dE_formula^2)")
print("  dE_formula = |E_H/n1^a - E_H/n2^a|  (no Z_eff!)")
print("=" * 90)
print()

# But wait: for all 2p-2p pairs, dE_formula = 0 (same n, same a)!
# So this can't help CO, NO, BF, CN.
# Those molecules NEED an energy difference to get the correction.
# Without Z_eff, all 2p orbitals look the same.

# What about using the PHASE difference instead?
# The mismatch between k1 and k2 creates a "detuning"
# that could play the role of dE in the 2-level model.
print("Alternative: Use wavevector mismatch as 'detuning'")
print("  delta_k = |k1 - k2| / (k1 + k2)")
print("  D = sqrt(D_cov^2 + c * (delta_k * E_scale)^2)")
print("  For homonuclear: delta_k = 0, D = |D_cov|. Same as current.")
print()

for c_val, c_label in [
    (1.0, '1'), (4.0, '4'), (d**2, 'd^2=9'),
    ((2*d+1), '2d+1=7'), (d**2+1, 'd^2+1=10'),
    (4*d, '4d=12'), (pi**2, 'pi^2'),
]:
    errs = []
    for mol in molecules:
        name, R, De_exp = mol[0], mol[1], mol[2]
        orb1, orb2 = mol[4], mol[5]
        D_cov, dE_orb, R_val, phase = compute_cov(mol)

        n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
        n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)
        b1 = 1 + beta_n * h1
        b2 = 1 + beta_n * h2
        k1 = 1.0 / n1**b1
        k2 = 1.0 / n2**b2
        dk = abs(k1 - k2) / (k1 + k2) if (k1+k2) > 0 else 0

        a1 = 2 + (1 - 2*l1) * alpha_n * h1
        a2 = 2 + (1 - 2*l2) * alpha_n * h2
        E_scale = np.sqrt(E_H / n1**a1 * E_H / n2**a2)

        detuning = dk * E_scale  # wavevector mismatch in energy units

        D_tot = np.sqrt(D_cov**2 + c_val * detuning**2)
        if D_cov < 0:
            D_tot = -D_tot
        err = abs((D_tot - De_exp) / De_exp * 100)
        errs.append(err)

    avg = np.mean(errs)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"  c={c_val:6.2f} ({c_label:>10s}): avg={avg:5.1f}%, "
          f"<5%:{w5:2d}/24, <10%:{w10:2d}/24, <20%:{w20:2d}/24")


# Detail for best
print()
print("Detail for c=4d=12:")
print(f"{'Mol':<6} {'De_exp':>7} {'D_cov':>7} {'detuning':>9} {'D_tot':>7} {'err%':>7}")
print("-" * 55)

for mol in molecules:
    name, R, De_exp = mol[0], mol[1], mol[2]
    orb1, orb2 = mol[4], mol[5]
    D_cov, dE_orb, R_val, phase = compute_cov(mol)

    n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
    n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)
    b1 = 1 + beta_n * h1
    b2 = 1 + beta_n * h2
    k1 = 1.0 / n1**b1
    k2 = 1.0 / n2**b2
    dk = abs(k1 - k2) / (k1 + k2) if (k1+k2) > 0 else 0

    a1 = 2 + (1 - 2*l1) * alpha_n * h1
    a2 = 2 + (1 - 2*l2) * alpha_n * h2
    E_scale = np.sqrt(E_H / n1**a1 * E_H / n2**a2)
    detuning = dk * E_scale

    c_best = 4*d
    D_tot = np.sqrt(D_cov**2 + c_best * detuning**2)
    if D_cov < 0:
        D_tot = -D_tot
    err = (D_tot - De_exp) / De_exp * 100
    flag = '***' if abs(err) < 2 else ' **' if abs(err) < 5 else '  *' if abs(err) < 10 else ''
    print(f"{name:<6} {De_exp:7.3f} {D_cov:7.3f} {detuning:9.4f} {D_tot:7.3f} {err:+6.1f}% {flag}")


# =============================================================================
# FINAL COMPARISON TABLE
# =============================================================================
print()
print("=" * 90)
print("  SUMMARY: One-force models vs current (with ionic)")
print("=" * 90)
print()
print(f"{'Model':<45} {'avg%':>6} {'<5%':>5} {'<10%':>5} {'<20%':>5}")
print("-" * 70)

models = {
    'Current + ionic (baseline)': 'baseline',
    '|sin| only (no ionic)': 'sin_only',
    'sqrt(D_cov^2 + 12*dk_E^2) (one force)': 'quadrature',
    'Harmonic phase, no ionic': 'harmonic',
}

for label, model in models.items():
    errs = []
    for mol in molecules:
        name, R, De_exp = mol[0], mol[1], mol[2]
        orb1, orb2 = mol[4], mol[5]
        D_cov, dE_orb, R_val, phase = compute_cov(mol)

        if model == 'baseline':
            V = max(abs(D_cov), 0.01)
            q = dE_orb / np.sqrt(dE_orb**2 + (2*V)**2) if dE_orb > 0 else 0
            D_ion = c_ionic * q**2 * 2 * E_H / R_val
            D_tot = D_cov + D_ion

        elif model == 'sin_only':
            D_tot = D_cov

        elif model == 'quadrature':
            n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
            n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)
            b1 = 1 + beta_n * h1; b2 = 1 + beta_n * h2
            k1 = 1.0/n1**b1; k2 = 1.0/n2**b2
            dk = abs(k1-k2)/(k1+k2) if (k1+k2)>0 else 0
            a1 = 2 + (1-2*l1)*alpha_n*h1; a2 = 2 + (1-2*l2)*alpha_n*h2
            E_s = np.sqrt(E_H/n1**a1 * E_H/n2**a2)
            detuning = dk * E_s
            D_tot = np.sqrt(D_cov**2 + 12*detuning**2)
            if D_cov < 0: D_tot = -D_tot

        elif model == 'harmonic':
            # Recompute with harmonic mean phase
            n1, l1, h1 = get_n(orb1), get_l(orb1), has_nodes(orb1)
            n2, l2, h2 = get_n(orb2), get_l(orb2), has_nodes(orb2)
            a1 = 2 + (1-2*l1)*alpha_n*h1; a2 = 2 + (1-2*l2)*alpha_n*h2
            b1 = 1 + beta_n*h1; b2 = 1 + beta_n*h2
            E_s = np.sqrt(E_H/n1**a1 * E_H/n2**a2)
            k1 = 1.0/n1**b1; k2 = 1.0/n2**b2
            sigma_phase = 4*R*k1*k2/(k1+k2) if (k1+k2)>0 else 0

            npb = sum(c for bt, c in mol[3] if 'pi' in bt and 'anti' not in bt)
            npa = sum(c for bt, c in mol[3] if 'pi' in bt and 'anti' in bt)
            ifa = (npa >= npb) if npb > 0 else True
            D_tot = 0
            for bt, cnt in mol[3]:
                ph = sigma_phase if ('sigma' in bt or bt in ('ss','sp')) else sigma_phase*f_pi
                cont = C_bond * E_s * abs(np.sin(ph))
                if 'anti' in bt:
                    fa = 1.0 if ('sigma' in bt or ifa) else f_anti
                    D_tot -= cnt * fa * cont
                else:
                    D_tot += cnt * cont

        err = abs((D_tot - De_exp) / De_exp * 100)
        errs.append(err)

    avg = np.mean(errs)
    w5 = sum(1 for e in errs if e < 5)
    w10 = sum(1 for e in errs if e < 10)
    w20 = sum(1 for e in errs if e < 20)
    print(f"{label:<45} {avg:6.1f} {w5:5d} {w10:5d} {w20:5d}")
