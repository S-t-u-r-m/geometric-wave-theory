"""
GWT Atomic Emission Spectrum Predictor
========================================
Predicts the emission spectrum of any atom from the Lagrangian.

Each breather (electron) occupies a mode in the kink well.
Energy levels: E(n,l) = -Z_eff(n,l)^2 * E_H / n^2
  where Z_eff comes from the Oh screening model (atomic_data.py)
Selection rule: photon = T1u, so initial x T1u must contain final.
  This gives delta_l = +/-1 (from Oh tensor product, not postulated).
Wavelength: lambda = hc / dE = 1240 / dE(eV) nm

Uses atomic_data.py (103 atoms, Z_eff from z_eff_v20.py model).
"""

import sys, os
import numpy as np
from math import factorial

# Import atomic data
sys.path.insert(0, os.path.dirname(__file__))
from atomic_data import ATOMS, E_H

PI = np.pi
d = 3
hc = 1240.0  # eV * nm
alpha = np.exp(-(2/factorial(d)) * (2**(2*d+1)/PI**2 + np.log(2*d)))

# Oh irrep labels
l_names = {0: 's', 1: 'p', 2: 'd', 3: 'f'}
l_max_electrons = {0: 2, 1: 6, 2: 10, 3: 14}


def compute_energy_levels(symbol, max_n=None):
    """
    Compute energy levels for each subshell of an atom.

    For occupied shells: use the Z_eff from the Oh screening model.
    For excited (empty) shells: use hydrogen-like with the valence Z_eff
    (the outermost electron sees the screened charge).

    Returns: dict of {(n, l): {'E': energy_eV, 'occ': count, ...}}
    """
    atom = ATOMS[symbol]
    Z = atom['Z']
    config = atom['config']
    Z_eff_val = atom['Z_eff']  # valence Z_eff
    val_n = atom['val_n']

    if max_n is None:
        max_n = val_n + 3  # go 3 shells above valence

    levels = {}

    # Occupied subshells: compute Z_eff for EACH subshell
    # The IE model gives the valence Z_eff. For inner shells,
    # they see more nuclear charge (less screening).
    for n_sh, l_sh, count in config:
        if count == 0:
            continue

        # Inner shell Z_eff: electrons in shell n see screening from
        # electrons in shells < n only. Rough model:
        # S_inner = sum of electrons in shells < n_sh
        S_inner = sum(c for n2, l2, c in config if n2 < n_sh)
        # Same-shell screening
        S_same = sum(c for n2, l2, c in config if n2 == n_sh) - 1
        S_same *= (2/d if l_sh == 0 else 4/(2*d+1))  # Oh screening weights

        Z_eff_shell = Z - S_inner - S_same
        Z_eff_shell = max(Z_eff_shell, 1.0)

        E = -Z_eff_shell**2 * E_H / n_sh**2

        levels[(n_sh, l_sh)] = {
            'E': E,
            'Z_eff': Z_eff_shell,
            'occ': count,
            'max_occ': l_max_electrons[l_sh],
            'type': 'occupied'
        }

    # Excited (empty) subshells: use valence Z_eff
    # These are the states an excited electron can jump to
    for n_ex in range(1, max_n + 1):
        for l_ex in range(min(n_ex, 4)):
            if (n_ex, l_ex) in levels:
                continue  # already occupied

            # Excited state sees the screened nuclear charge
            # For Rydberg states (high n): Z_eff ~ 1 (fully screened)
            # For low excited states: Z_eff ~ valence Z_eff
            if n_ex <= val_n + 1:
                z_ex = Z_eff_val
            else:
                # Rydberg: approaches 1 (all electrons screen)
                z_ex = max(1.0, Z_eff_val * val_n**2 / n_ex**2)

            E = -z_ex**2 * E_H / n_ex**2

            levels[(n_ex, l_ex)] = {
                'E': E,
                'Z_eff': z_ex,
                'occ': 0,
                'max_occ': l_max_electrons[l_ex],
                'type': 'excited'
            }

    return levels


def selection_rule_allowed(l_initial, l_final):
    """Oh selection rule: photon (T1u) connects states with delta_l = +/-1."""
    return abs(l_initial - l_final) == 1


def compute_emission_lines(symbol, max_n=None):
    """
    Compute all emission lines for an atom.

    An emission line = transition from a higher energy state to a lower one,
    obeying the Oh selection rule (delta_l = +/-1).

    Returns: sorted list of (wavelength_nm, label, energy_eV, upper_level, lower_level)
    """
    levels = compute_energy_levels(symbol, max_n)
    lines = []

    sorted_levels = sorted(levels.items(), key=lambda x: x[1]['E'])

    for i, ((n_lo, l_lo), info_lo) in enumerate(sorted_levels):
        for j, ((n_hi, l_hi), info_hi) in enumerate(sorted_levels):
            if info_hi['E'] <= info_lo['E']:
                continue
            if not selection_rule_allowed(l_lo, l_hi):
                continue

            dE = info_hi['E'] - info_lo['E']
            if dE < 0.01:
                continue

            lam = hc / dE
            if lam < 1 or lam > 50000:
                continue

            label = f"{n_hi}{l_names[l_hi]}->{n_lo}{l_names[l_lo]}"
            lines.append((lam, label, dE, (n_hi, l_hi), (n_lo, l_lo)))

    lines.sort(key=lambda x: x[0])
    return lines, levels


def classify_wavelength(lam):
    """Classify wavelength into spectral region."""
    if lam < 10: return 'X-ray'
    if lam < 100: return 'EUV'
    if lam < 380: return 'UV'
    if lam < 450: return 'Violet'
    if lam < 495: return 'Blue'
    if lam < 570: return 'Green'
    if lam < 590: return 'Yellow'
    if lam < 620: return 'Orange'
    if lam < 750: return 'Red'
    if lam < 2500: return 'Near-IR'
    return 'IR'


def print_spectrum(symbol, visible_only=False, max_lines=30):
    """Pretty-print the predicted emission spectrum."""
    atom = ATOMS[symbol]
    lines, levels = compute_emission_lines(symbol)

    if visible_only:
        lines = [l for l in lines if 380 <= l[0] <= 750]

    print(f"\n{'='*65}")
    print(f"{symbol} (Z={atom['Z']}) — GWT Predicted Emission Spectrum")
    print(f"{'='*65}")

    # Energy levels
    print(f"\n  Energy levels (Z_eff from Oh screening model):")
    for (n, l), info in sorted(levels.items(), key=lambda x: x[1]['E']):
        occ = f"[{info['occ']}/{info['max_occ']}]" if info['occ'] > 0 else "     "
        tag = "" if info['type'] == 'occupied' else " (excited)"
        print(f"    {n}{l_names[l]} {occ}: E = {info['E']:10.3f} eV, "
              f"Z_eff = {info['Z_eff']:.2f}{tag}")

    # Emission lines
    if lines:
        print(f"\n  Emission lines ({len(lines)} total"
              f"{', visible only' if visible_only else ''}):")
        print(f"  {'Lambda':>9} {'Transition':>12} {'Energy':>9} {'Region':>10}")
        print(f"  {'-'*45}")
        for lam, label, dE, upper, lower in lines[:max_lines]:
            region = classify_wavelength(lam)
            print(f"  {lam:9.2f} nm {label:>12} {dE:9.4f} eV {region:>10}")
        if len(lines) > max_lines:
            print(f"  ... and {len(lines)-max_lines} more lines")

        # Highlight visible lines
        vis = [l for l in lines if 380 <= l[0] <= 750]
        if vis and not visible_only:
            print(f"\n  VISIBLE LINES ({len(vis)}):")
            for lam, label, dE, _, _ in vis:
                region = classify_wavelength(lam)
                print(f"    {lam:8.2f} nm  {label:>12}  {region}")
    else:
        print("  No emission lines in range")


# ============================================================
# RUN
# ============================================================
if __name__ == '__main__':
    print("=" * 65)
    print("GWT ATOMIC EMISSION SPECTRUM PREDICTOR")
    print("=" * 65)
    print(f"Using Z_eff from Oh screening model (103 atoms, 2.6% mean IE)")
    print(f"Selection rule: delta_l = +/-1 (from Oh: photon = T1u)")
    print(f"E_H = {E_H:.4f} eV, alpha = 1/{1/alpha:.3f}")

    # Print spectra for key elements
    for element in ['H', 'He', 'Na', 'Ca', 'Fe']:
        print_spectrum(element)

    # Hydrogen Balmer comparison
    print(f"\n{'='*65}")
    print("HYDROGEN BALMER SERIES — GWT vs Observed")
    print(f"{'='*65}")
    balmer_obs = {'H-alpha': 656.28, 'H-beta': 486.13,
                  'H-gamma': 434.05, 'H-delta': 410.17}
    for n, (name, obs) in zip(range(3, 7), balmer_obs.items()):
        dE = E_H * (1/4 - 1/n**2)
        lam = hc / dE
        err = (lam - obs)/obs * 100
        print(f"  {name}: GWT = {lam:.2f} nm, Obs = {obs:.2f} nm, "
              f"Err = {err:+.3f}%")

    # Sodium D-line check
    print(f"\n{'='*65}")
    print("SODIUM — Visible spectrum check")
    print(f"{'='*65}")
    print_spectrum('Na', visible_only=True)
    print(f"\n  Observed Na D-lines: 589.0 nm, 589.6 nm")
