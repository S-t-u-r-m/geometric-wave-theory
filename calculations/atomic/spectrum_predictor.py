"""
GWT Atomic Emission Spectrum Predictor
========================================
Predicts the emission spectrum of any atom from the Lagrangian.

Each breather (electron) occupies a mode in the kink well (Pöschl-Teller).
Energy levels: E(n,l) = -Z_eff(n,l)^2 * E_H / n^2
Selection rule: photon = T1u, so initial ⊗ T1u must contain final.
  This gives Δl = ±1 (derived from Oh tensor product, not postulated).
Wavelength: λ = hc / ΔE = 1240 / ΔE(eV) nm

All inputs from GWT: alpha (derived), Z_eff (from IE model), Oh irreps.
"""

import numpy as np
from math import factorial

PI = np.pi
d = 3

# GWT constants
alpha = np.exp(-(2/factorial(d)) * (2**(2*d+1)/PI**2 + np.log(2*d)))
E_H = alpha**2 * 0.51100e6 / 2  # 13.60 eV (Rydberg energy)
hc = 1240.0  # eV * nm


# ============================================================
# Z_eff MODEL (simplified from z_eff_v20.py)
# ============================================================
# For spectral predictions, we need Z_eff for each subshell (n, l).
# The screening comes from inner electrons reducing the nuclear charge.
# Full model is in calculations/atomic/z_eff_v20.py.
# Here we use the Slater screening rules adapted with GWT corrections.

def slater_z_eff(Z, n, l, config):
    """
    Compute effective nuclear charge for subshell (n, l).

    config: list of (n, l, count) for all occupied subshells
    Screening rules from GWT (Oh-derived):
      Same shell, same l:  s = 2/d = 0.667
      Same shell, diff l:  s = 4/(2d+1) = 0.571
      One shell below:     s = 2/(d+2) = 0.400 (closed shell)
      Two+ shells below:   s = 1.000 (full screening)
    """
    g_same = 2.0 / d           # 0.667
    g_diff = 4.0 / (2*d + 1)   # 0.571
    g_inner = 2.0 / (d + 2)    # 0.400

    S = 0  # total screening
    for ni, li, count in config:
        if ni == n and li == l:
            # Same subshell: other electrons screen partially
            S += (count - 1) * g_same  # -1 for the electron itself
        elif ni == n and li != l:
            # Same shell, different subshell
            S += count * g_diff
        elif ni == n - 1:
            # One shell below: partial screening
            S += count * g_inner
        elif ni < n - 1:
            # Deep inner shells: full screening
            S += count * 1.0

    return max(Z - S, 1.0)


# ============================================================
# ELECTRON CONFIGURATIONS
# ============================================================

# Standard electron configurations: (n, l, count)
# l: 0=s, 1=p, 2=d, 3=f
configs = {
    'H':  [(1,0,1)],
    'He': [(1,0,2)],
    'Li': [(1,0,2), (2,0,1)],
    'Be': [(1,0,2), (2,0,2)],
    'B':  [(1,0,2), (2,0,2), (2,1,1)],
    'C':  [(1,0,2), (2,0,2), (2,1,2)],
    'N':  [(1,0,2), (2,0,2), (2,1,3)],
    'O':  [(1,0,2), (2,0,2), (2,1,4)],
    'F':  [(1,0,2), (2,0,2), (2,1,5)],
    'Ne': [(1,0,2), (2,0,2), (2,1,6)],
    'Na': [(1,0,2), (2,0,2), (2,1,6), (3,0,1)],
    'Mg': [(1,0,2), (2,0,2), (2,1,6), (3,0,2)],
    'Al': [(1,0,2), (2,0,2), (2,1,6), (3,0,2), (3,1,1)],
    'Si': [(1,0,2), (2,0,2), (2,1,6), (3,0,2), (3,1,2)],
    'Fe': [(1,0,2), (2,0,2), (2,1,6), (3,0,2), (3,1,6), (3,2,6), (4,0,2)],
    'Ca': [(1,0,2), (2,0,2), (2,1,6), (3,0,2), (3,1,6), (4,0,2)],
}

# Oh irrep labels for each l
l_irreps = {0: 'A1g', 1: 'T1u', 2: 'T2g+Eg', 3: 'A2u+T1u+T2u'}
l_names = {0: 's', 1: 'p', 2: 'd', 3: 'f'}

# Selection rule: Δl = ±1 (from Oh: initial ⊗ T1u must contain final)
def is_allowed(l_initial, l_final):
    return abs(l_initial - l_final) == 1


# ============================================================
# SPECTRUM CALCULATOR
# ============================================================

def compute_spectrum(element, Z, config, max_n_excited=5):
    """
    Compute emission spectrum for an atom.

    For each occupied subshell, compute transitions to nearby empty/available
    subshells following the Oh selection rule Δl = ±1.

    Returns list of (wavelength_nm, transition_label, energy_eV)
    """
    lines = []

    # Compute energy levels for ground state and nearby excited states
    levels = {}

    # Ground state subshells
    for n, l, count in config:
        z_eff = slater_z_eff(Z, n, l, config)
        E = -z_eff**2 * E_H / n**2
        levels[(n, l)] = {'E': E, 'Z_eff': z_eff, 'occ': count, 'max': 2*(2*l+1)}

    # Excited state subshells (empty ones nearby)
    max_n = max(n for n, l, c in config)
    for n in range(1, max_n_excited + 1):
        for l in range(min(n, 4)):  # l < n
            if (n, l) not in levels:
                z_eff = slater_z_eff(Z, n, l, config)
                E = -z_eff**2 * E_H / n**2
                levels[(n, l)] = {'E': E, 'Z_eff': z_eff, 'occ': 0, 'max': 2*(2*l+1)}

    # Find transitions: occupied -> higher (absorption) or higher -> occupied (emission)
    # For emission: excited state falls to lower state
    sorted_levels = sorted(levels.items(), key=lambda x: x[1]['E'])

    for i, ((n_lo, l_lo), info_lo) in enumerate(sorted_levels):
        for j, ((n_hi, l_hi), info_hi) in enumerate(sorted_levels):
            if info_hi['E'] <= info_lo['E']:
                continue  # must be higher energy
            if not is_allowed(l_lo, l_hi):
                continue  # Oh selection rule

            dE = info_hi['E'] - info_lo['E']
            if dE <= 0.01:
                continue  # too small

            lam = hc / dE  # wavelength in nm
            if lam > 10000 or lam < 10:
                continue  # outside useful range

            label = f"{n_hi}{l_names[l_hi]}->{n_lo}{l_names[l_lo]}"
            lines.append((lam, label, dE,
                         info_hi['Z_eff'], info_lo['Z_eff']))

    # Sort by wavelength
    lines.sort(key=lambda x: x[0])
    return lines, levels


# ============================================================
# RUN: PREDICT SPECTRA
# ============================================================

print("=" * 70)
print("GWT ATOMIC EMISSION SPECTRUM PREDICTOR")
print("=" * 70)
print(f"All from alpha = 1/{1/alpha:.3f}, E_H = {E_H:.4f} eV")
print(f"Selection rule: Delta_l = +/-1 (from Oh: photon = T1u)")
print()

Z_values = {'H': 1, 'He': 2, 'Li': 3, 'C': 6, 'N': 7, 'O': 8,
            'Na': 11, 'Fe': 26, 'Ca': 20}

for element in ['H', 'He', 'Na', 'Fe']:
    Z = Z_values[element]
    config = configs[element]

    print(f"\n{'='*60}")
    print(f"{element} (Z={Z}) — Predicted Emission Spectrum")
    print(f"{'='*60}")

    lines, levels = compute_spectrum(element, Z, config)

    # Show energy levels
    print(f"\n  Energy levels:")
    for (n, l), info in sorted(levels.items(), key=lambda x: x[1]['E']):
        occ_str = f"[{info['occ']}/{info['max']}]" if info['occ'] > 0 else "     "
        print(f"    {n}{l_names[l]} {occ_str}: E = {info['E']:8.3f} eV, Z_eff = {info['Z_eff']:.2f}")

    # Show spectral lines
    if lines:
        print(f"\n  Emission lines ({len(lines)} predicted):")
        print(f"  {'Lambda':>8} {'Transition':>10} {'Energy':>8} {'Region':>10}")
        print(f"  {'-'*40}")
        for lam, label, dE, z_hi, z_lo in lines[:20]:  # show first 20
            region = ('UV' if lam < 380 else
                     ('Violet' if lam < 450 else
                     ('Blue' if lam < 495 else
                     ('Green' if lam < 570 else
                     ('Yellow' if lam < 590 else
                     ('Orange' if lam < 620 else
                     ('Red' if lam < 750 else 'IR')))))))
            print(f"  {lam:8.1f} nm {label:>10} {dE:8.3f} eV {region:>10}")
        if len(lines) > 20:
            print(f"  ... and {len(lines)-20} more lines")
    else:
        print("  No emission lines in 10-10000 nm range")

# Specific comparison: Hydrogen Balmer series
print(f"\n{'='*60}")
print("HYDROGEN BALMER SERIES — GWT vs Observed")
print(f"{'='*60}")
balmer_obs = {'H-alpha': 656.28, 'H-beta': 486.13, 'H-gamma': 434.05, 'H-delta': 410.17}
for n, (name, obs) in zip(range(3, 7), balmer_obs.items()):
    dE = E_H * (1/4 - 1/n**2)
    lam = hc / dE
    err = (lam - obs)/obs * 100
    print(f"  {name}: GWT = {lam:.2f} nm, Obs = {obs:.2f} nm, Err = {err:+.3f}%")

# Sodium D-lines (the classic yellow doublet)
print(f"\n{'='*60}")
print("SODIUM D-LINES — GWT Prediction")
print(f"{'='*60}")
lines_Na, _ = compute_spectrum('Na', 11, configs['Na'])
for lam, label, dE, z_hi, z_lo in lines_Na:
    if 580 < lam < 600:
        print(f"  {label}: {lam:.1f} nm ({dE:.4f} eV)")
print(f"  Observed: Na D1 = 589.6 nm, Na D2 = 589.0 nm")
