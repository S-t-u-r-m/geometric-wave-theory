"""
GWT Diamagnetic Susceptibility Calculator
==========================================
Predicts diamagnetic susceptibility (chi_dia) of atoms from Z_eff.

The key quantity is <r^2> — the mean square radius of each harmonic mode.
For a hydrogenic subshell with effective nuclear charge Z_eff:

  <r^2>_{n,l} = (a_0^2 / Z_eff^2) * n^2 * [5*n^2 + 1 - 3*l*(l+1)] / 2

The total diamagnetic susceptibility:

  chi_dia = -(e^2 / 6*m_e*c^2) * sum_modes(N_mode * <r^2>_mode)

In CGS units:  chi_dia = -(e^2 / 6*m_e*c^2) * <r^2>_total
             = -alpha^2 * a_0^2 / 6 * <r^2>_total  (in Bohr^2)

Convention: chi in units of 10^-6 cm^3/mol (CGS molar susceptibility).

This calculation feeds into bonding: <r^2> determines harmonic cloud extent,
which determines overlap between neighboring kinks = bond energy.

All inputs from GWT (d=3 Lagrangian). Zero free parameters.
"""

import sys, os, io
import numpy as np
from math import factorial

if not hasattr(sys.stdout, '_gwt_wrapped'):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stdout._gwt_wrapped = True

PI = np.pi
d = 3
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/PI**2 + np.log(2*d)))
a_0_cm = 5.29177e-9       # Bohr radius in cm
N_A = 6.02214e23           # Avogadro's number
e_cgs = 4.80321e-10        # electron charge in esu
m_e_cgs = 9.10938e-28      # electron mass in grams
c_cgs = 2.99792e10         # speed of light in cm/s

# Prefactor for chi_dia in CGS molar units (10^-6 cm^3/mol):
# chi = -N_A * e^2 / (6 * m_e * c^2) * <r^2>_total
# where <r^2> is in cm^2
CHI_PREFACTOR = -N_A * e_cgs**2 / (6 * m_e_cgs * c_cgs**2)
# This gives chi in cm^3/mol; multiply by 1e6 to get 10^-6 cm^3/mol

sys.path.insert(0, os.path.dirname(__file__))
from atomic_data import ATOMS
from z_eff_subshell import compute_all_levels, compute_subshell_energy

l_names = {0: 's', 1: 'p', 2: 'd', 3: 'f'}


def r_squared_hydrogenic(n, l, Z_eff):
    """
    <r^2> for a hydrogenic orbital with effective charge Z_eff.

    Formula (exact for hydrogen-like):
      <r^2> = (a_0/Z_eff)^2 * n^2 * [5*n^2 + 1 - 3*l*(l+1)] / 2

    Returns <r^2> in units of a_0^2 (Bohr radii squared).
    """
    return (n / Z_eff)**2 * (5 * n**2 + 1 - 3 * l * (l + 1)) / 2


def compute_diamagnetic_chi(symbol, use_quantum_defect=False):
    """
    Compute diamagnetic susceptibility for an atom.

    Args:
        symbol: element symbol (e.g., 'He', 'Ne', 'Ar')
        use_quantum_defect: if True, use quantum defect Z_eff instead of power-law

    Returns:
        chi_dia in 10^-6 cm^3/mol (CGS molar susceptibility)
        Also returns per-subshell breakdown.
    """
    atom = ATOMS[symbol]
    Z = atom['Z']
    config = atom['config']

    subshell_data = []
    r2_total = 0.0

    for n_sh, l_sh, count in config:
        if count == 0:
            continue

        # Get Z_eff for this subshell
        if use_quantum_defect:
            # Use quantum defect formula
            n_core = max((n2 for n2, l2, c in config if n2 < n_sh and c > 0), default=0)
            n_p_core = sum(1 for n2, l2, c in config if l2 == 1 and n2 < n_sh and c > 0)
            paired = count >= 2 and l_sh == 0
            from z_eff_subshell import z_eff_from_defect
            Z_eff = z_eff_from_defect(n_sh, l_sh, n_core, n_p_core, paired)
        else:
            # Use power-law model
            E, Z_eff, Z_net, alpha, S = compute_subshell_energy(Z, config, n_sh, l_sh)

        r2 = r_squared_hydrogenic(n_sh, l_sh, Z_eff)
        r2_contribution = count * r2
        r2_total += r2_contribution

        subshell_data.append({
            'n': n_sh, 'l': l_sh, 'count': count,
            'Z_eff': Z_eff, 'r2': r2,
            'r2_total': r2_contribution,
            'r_rms': np.sqrt(r2) * a_0_cm * 1e8  # in Angstroms
        })

    # chi in cm^3/mol
    chi_cm3 = CHI_PREFACTOR * r2_total * a_0_cm**2
    # Convert to 10^-6 cm^3/mol
    chi = chi_cm3 * 1e6

    return chi, r2_total, subshell_data


# ============================================================
# OBSERVED VALUES (10^-6 cm^3/mol, molar diamagnetic susceptibility)
# Source: CRC Handbook of Chemistry and Physics
# ============================================================
OBSERVED_CHI = {
    # Noble gases
    'He': -1.88,
    'Ne': -7.2,
    'Ar': -19.6,
    'Kr': -28.8,
    'Xe': -45.5,
    # Other closed-shell / diamagnetic atoms/ions
    'Be': -9.0,
    'Mg': -11.6,  # Note: these are for bulk metal, atomic may differ
    'Zn': -9.15,  # diamagnetic metal
}


if __name__ == '__main__':
    print("GWT Diamagnetic Susceptibility from Z_eff")
    print("=" * 75)
    print(f"  Prefactor: {CHI_PREFACTOR:.4e} cm^3/mol per a_0^2")
    print(f"  a_0 = {a_0_cm:.5e} cm")
    print()

    # Test on noble gases first (fully paired = purely diamagnetic)
    noble_gases = ['He', 'Ne', 'Ar', 'Kr', 'Xe']

    print(f"{'Atom':>4} {'Z':>3} | {'chi_GWT':>10} {'chi_obs':>10} {'error':>8} | {'<r2>':>8} | Subshell breakdown")
    print("-" * 95)

    for elem in noble_gases:
        if elem not in ATOMS:
            print(f"  {elem}: not in atomic_data.py")
            continue

        chi_gwt, r2_total, subshells = compute_diamagnetic_chi(elem)
        chi_obs = OBSERVED_CHI.get(elem, None)

        err_str = ""
        if chi_obs is not None:
            err = (chi_gwt - chi_obs) / abs(chi_obs) * 100
            err_str = f"{err:+.1f}%"

        # Subshell breakdown string
        sub_str = "  ".join(
            f"{s['n']}{l_names[s['l']]}({s['count']}): Z={s['Z_eff']:.2f}, r={s['r_rms']:.3f}A"
            for s in subshells
        )

        print(f"{elem:>4} {ATOMS[elem]['Z']:>3} | {chi_gwt:>10.2f} {chi_obs:>10.2f} {err_str:>8} | {r2_total:>8.3f} | {sub_str}")

    # Detailed breakdown for each noble gas
    print(f"\n{'='*75}")
    print("DETAILED SUBSHELL BREAKDOWN")
    print(f"{'='*75}")

    for elem in noble_gases:
        if elem not in ATOMS:
            continue

        chi_gwt, r2_total, subshells = compute_diamagnetic_chi(elem)
        chi_obs = OBSERVED_CHI.get(elem, None)

        print(f"\n{elem} (Z={ATOMS[elem]['Z']}): chi_GWT = {chi_gwt:.2f}, chi_obs = {chi_obs}")
        print(f"  {'Shell':>6} {'N_e':>4} {'Z_eff':>7} {'<r2>':>10} {'N*<r2>':>10} {'r_rms(A)':>9} {'% of total':>10}")
        print(f"  {'-'*62}")

        for s in subshells:
            pct = s['r2_total'] / r2_total * 100
            print(f"  {s['n']}{l_names[s['l']]:1s}    {s['count']:>4} {s['Z_eff']:>7.3f} "
                  f"{s['r2']:>10.4f} {s['r2_total']:>10.4f} {s['r_rms']:>9.4f} {pct:>9.1f}%")

        print(f"  {'Total':>6} {sum(s['count'] for s in subshells):>4} {'':>7} "
              f"{'':>10} {r2_total:>10.4f}")

    # Also try quantum defect Z_eff for comparison
    print(f"\n{'='*75}")
    print("COMPARISON: Power-law Z_eff vs Quantum Defect Z_eff")
    print(f"{'='*75}")
    print(f"{'Atom':>4} | {'chi_power':>10} {'chi_QD':>10} {'chi_obs':>10} | {'err_pow':>8} {'err_QD':>8}")
    print("-" * 70)

    for elem in noble_gases:
        if elem not in ATOMS:
            continue

        chi_pow, _, _ = compute_diamagnetic_chi(elem, use_quantum_defect=False)
        chi_qd, _, _ = compute_diamagnetic_chi(elem, use_quantum_defect=True)
        chi_obs = OBSERVED_CHI.get(elem, None)

        err_pow = (chi_pow - chi_obs) / abs(chi_obs) * 100 if chi_obs else 0
        err_qd = (chi_qd - chi_obs) / abs(chi_obs) * 100 if chi_obs else 0

        obs_str = f"{chi_obs:.2f}" if chi_obs else "N/A"
        print(f"{elem:>4} | {chi_pow:>10.2f} {chi_qd:>10.2f} {obs_str:>10} | {err_pow:>+7.1f}% {err_qd:>+7.1f}%")

    # ============================================================
    # INVERSE ANALYSIS: What Z_eff does observed chi require?
    # ============================================================
    print(f"\n{'='*75}")
    print("INVERSE ANALYSIS: What spatial Z_eff does observed chi require?")
    print(f"{'='*75}")
    print()
    print("If we assume ALL modes need the SAME correction factor f")
    print("such that Z_eff(spatial) = Z_eff(energy) * f:")
    print("Then <r2> scales as 1/f^2, so chi scales as 1/f^2.")
    print()

    for elem in noble_gases:
        if elem not in ATOMS:
            continue
        chi_gwt, r2_gwt, subshells = compute_diamagnetic_chi(elem)
        chi_obs = OBSERVED_CHI.get(elem, None)
        if chi_obs is None:
            continue

        # chi_gwt / chi_obs = f^(-2) => f = sqrt(chi_gwt / chi_obs)
        # But both are negative, so use magnitudes
        f = np.sqrt(abs(chi_gwt) / abs(chi_obs))
        Z = ATOMS[elem]['Z']
        print(f"  {elem:>2} (Z={Z:>2}): chi_gwt/chi_obs = {abs(chi_gwt)/abs(chi_obs):>7.1f}x"
              f"  => f = {f:.4f}  (Z_eff needs ×{f:.2f})")

    print()
    print("Pattern check: is f related to Z, d, or phase space?")
    print()
    factors = []
    for elem in noble_gases:
        if elem not in ATOMS or elem not in OBSERVED_CHI:
            continue
        chi_gwt, _, _ = compute_diamagnetic_chi(elem)
        chi_obs = OBSERVED_CHI[elem]
        f = np.sqrt(abs(chi_gwt) / abs(chi_obs))
        Z = ATOMS[elem]['Z']
        n_max = max(n for n, l, c in ATOMS[elem]['config'] if c > 0)
        factors.append((elem, Z, n_max, f))
        # Check various ratios
        print(f"  {elem:>2}: f={f:.3f}, Z={Z}, n_max={n_max}, "
              f"f/sqrt(Z)={f/np.sqrt(Z):.3f}, "
              f"f/Z^(1/3)={f/Z**(1/3):.3f}, "
              f"f/n_max={f/n_max:.3f}, "
              f"f/sqrt(n_max)={f/np.sqrt(n_max):.3f}")

    # Check if the valence shell alone could explain it
    print(f"\n{'='*75}")
    print("VALENCE SHELL DOMINANCE CHECK")
    print(f"{'='*75}")
    print("The outermost shell dominates <r2>. What Z_eff does it need?")
    print()

    for elem in noble_gases:
        if elem not in ATOMS or elem not in OBSERVED_CHI:
            continue

        chi_obs = OBSERVED_CHI[elem]
        chi_gwt, r2_gwt, subshells = compute_diamagnetic_chi(elem)
        Z = ATOMS[elem]['Z']

        # The required total <r2>
        r2_needed = abs(chi_obs) / abs(CHI_PREFACTOR) / a_0_cm**2 * 1e-6

        # The outermost shell contribution
        outer = subshells[-1]
        n_out, l_out = outer['n'], outer['l']

        # What Z_eff would give the right <r2> for JUST the outer shell?
        # Assume outer shell dominates and inner shells contribute ~10%
        # <r2>_outer * N_outer = r2_needed * 0.9 (rough)
        # <r2>_outer = (n/Z_eff)^2 * [5n^2+1-3l(l+1)]/2
        # So Z_eff_needed = n * sqrt([5n^2+1-3l(l+1)] * N_outer / (2 * r2_needed))
        r2_factor = (5 * n_out**2 + 1 - 3 * l_out * (l_out + 1)) / 2
        z_needed = n_out * np.sqrt(r2_factor * outer['count'] / r2_needed)

        print(f"  {elem:>2}: r2_needed = {r2_needed:.3f} a0^2, "
              f"outer = {n_out}{l_names[l_out]}({outer['count']}), "
              f"Z_eff(energy) = {outer['Z_eff']:.3f}, "
              f"Z_eff(spatial needed) ~ {z_needed:.3f}, "
              f"ratio = {z_needed/outer['Z_eff']:.3f}")
