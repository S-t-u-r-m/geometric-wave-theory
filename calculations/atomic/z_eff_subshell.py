"""
GWT Per-Subshell Z_eff Calculator
===================================
Extends z_eff_v20.py to compute Z_eff for ANY subshell (n, l),
not just the valence shell. Needed for spectroscopy and improved
bond energies.

For a target subshell (n_t, l_t), the screening comes from all
OTHER electrons:
  - Electrons in lower shells: screen strongly
  - Electrons in the same shell: screen partially (Oh CG weights)
  - Electrons in the same subshell: screen by same-subshell weight

The Oh screening matrix is the same as z_eff_v20.py.
The alpha exponent is computed per-subshell.

Usage:
  levels = compute_all_levels('Na')
  for (n, l), info in levels.items():
      print(f"{n}{l_name[l]}: E = {info['E']:.3f} eV, Z_eff = {info['Z_eff']:.3f}")
"""

import sys, os, io
import numpy as np
from math import factorial, comb

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1) * PI - 2)
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/PI**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6

# Oh channel weights
w_sigma = 1.0
w_pi = np.cos(PI / d)         # 0.5
w_delta = np.cos(2 * PI / d)  # -0.5

l_names = {0: 's', 1: 'p', 2: 'd', 3: 'f'}
l_max_e = {0: 2, 1: 6, 2: 10, 3: 14}

# Import atom database
sys.path.insert(0, os.path.dirname(__file__))
from atomic_data import ATOMS


def screening_for_target(Z, config, n_target, l_target):
    """
    Compute screening of the nuclear charge seen by an electron
    in subshell (n_target, l_target).

    Screening comes from ALL other electrons:
      - Same subshell: (count - 1) * same_weight
      - Same shell, different l: count * diff_weight
      - Inner shells: count * inner_weight (Oh CG dependent)
    """
    S = 0.0

    for n_sh, l_sh, count in config:
        if count == 0:
            continue

        if n_sh == n_target and l_sh == l_target:
            # Same subshell: other electrons screen partially
            # Oh: same-mode coupling A1g(T1u^2)=1 for p, A1g(A1g^2)=1 for s
            # Screening weight = 2/d = 2/3 for same subshell
            same_weight = 2.0 / d  # 0.667
            S += (count - 1) * same_weight

        elif n_sh == n_target:
            # Same shell, different subshell
            # Oh: cross-subshell coupling
            # s-p: A1g x T1u = T1u (A1g = 0, weak coupling)
            # p-s: T1u x A1g = T1u (same)
            diff_weight = 4.0 / (2*d + 1)  # 0.571
            S += count * diff_weight

        elif n_sh == n_target - 1:
            # One shell below: partial screening
            # Depends on l_sh and l_target (Oh selection rules)
            if l_sh <= 1:
                # s,p inner shell: screens via T1u mediator
                S += count * w_pi
            elif l_sh == 2:
                # d inner shell: depends on target l
                if l_target == 1:
                    # d -> p: Oh ALLOWED
                    S += count * w_pi
                else:
                    # d -> s/d: Oh FORBIDDEN -> anti-screening
                    S += count * w_delta
            elif l_sh == 3:
                # f inner shell
                if l_target == 1:
                    S += min(count, 3) * w_pi + max(0, min(count, 7) - 3) * w_delta / d
                else:
                    S += min(count, 7) * w_delta / d

        elif n_sh < n_target - 1:
            # Deep inner shell: nearly full screening
            # Blend toward 1.0 (full screening) at depth
            depth = n_target - n_sh
            blend = 1.0 - 1.0 / (d * depth)  # approaches 1 for large depth
            S += count * blend

        # Outer shells (n_sh > n_target): no screening from outer electrons

    return S


def compute_alpha_for_target(config, n_target, l_target, S_core):
    """
    Compute the alpha exponent for a specific subshell.

    This determines how Z_net maps to Z_eff: Z_eff = Z_net^alpha
    """
    # Count electrons in the target shell
    count_in_shell = sum(c for n, l, c in config if n == n_target)
    count_in_subshell = sum(c for n, l, c in config if n == n_target and l == l_target)

    # Base alpha from the subshell type
    if l_target == 0:
        # s-orbital
        if count_in_subshell >= 2:
            C = 1.0   # paired: constructive
        else:
            C = -1.0  # single: destructive
        alpha = (d * n_target + C) / (d**2 * n_target)

    elif l_target == 1:
        # p-orbital: more complex coupling
        pa = min(count_in_subshell, d)  # unpaired
        pl = max(0, count_in_subshell - d)  # pairs beyond half-fill

        # N_eff from Oh coupling
        N_eff = pa
        if pl > 0:
            # First pair: Hund penalty
            w1 = (n_target**2 - d) / n_target**2
            N_eff += min(pl, 1) * w1
            # Remaining pairs
            if pl > 1:
                N_eff += (pl - 1) * (1 + w_pi)

        alpha = (d + w_pi * N_eff - 1) / d**2

    elif l_target == 2:
        # d-orbital
        if count_in_subshell <= d:
            C = w_delta  # underfilled: anti-screening character
        else:
            C = 1.0  # overfilled: constructive
        alpha = (d * n_target + C) / (d**2 * n_target)

    else:
        # f-orbital
        alpha = (d * n_target - 1) / (d**2 * n_target)

    # Clamp
    alpha = max(alpha, 0.05)
    alpha = min(alpha, 0.8)

    return alpha


def compute_subshell_energy(Z, config, n_target, l_target):
    """
    Compute the energy of a specific subshell (n, l).

    E(n, l) = -(Z_eff / n)^2 * E_H

    where Z_eff = (Z - S_core)^alpha
    """
    S = screening_for_target(Z, config, n_target, l_target)
    Z_net = max(Z - S, 0.5)
    alpha = compute_alpha_for_target(config, n_target, l_target, S)
    Z_eff = Z_net ** alpha
    E = -(Z_eff / n_target)**2 * E_H

    return E, Z_eff, Z_net, alpha, S


def compute_all_levels(symbol, max_n=None):
    """
    Compute energy levels for ALL subshells of an atom:
    occupied ones from the config, excited ones up to max_n.

    Returns: dict of {(n, l): {'E': eV, 'Z_eff': float, 'occ': int, ...}}
    """
    atom = ATOMS[symbol]
    Z = atom['Z']
    config = atom['config']
    val_n = atom['val_n']

    if max_n is None:
        max_n = val_n + 3

    levels = {}

    # Occupied subshells
    for n_sh, l_sh, count in config:
        E, Z_eff, Z_net, alpha, S = compute_subshell_energy(Z, config, n_sh, l_sh)
        levels[(n_sh, l_sh)] = {
            'E': E, 'Z_eff': Z_eff, 'Z_net': Z_net, 'alpha': alpha,
            'S': S, 'occ': count, 'max_occ': l_max_e[l_sh], 'type': 'occupied'
        }

    # Excited subshells (empty)
    # Key: the excited electron sees ALL ground-state electrons as screening.
    # It's NOT in the config (it's empty), so screening = sum of ALL electrons.
    # This is different from occupied shells where one electron is subtracted.
    #
    # Create a config that includes a virtual electron in the excited state
    # so the screening function correctly counts (count-1) for same-subshell.
    for n_ex in range(1, max_n + 1):
        for l_ex in range(min(n_ex, 4)):
            if (n_ex, l_ex) in levels:
                continue

            # Add a virtual electron to the excited subshell
            config_excited = list(config) + [(n_ex, l_ex, 1)]
            E, Z_eff, Z_net, alpha, S = compute_subshell_energy(Z, config_excited, n_ex, l_ex)
            levels[(n_ex, l_ex)] = {
                'E': E, 'Z_eff': Z_eff, 'Z_net': Z_net, 'alpha': alpha,
                'S': S, 'occ': 0, 'max_occ': l_max_e[l_ex], 'type': 'excited'
            }

    return levels


# ============================================================
# TEST
# ============================================================
if __name__ == '__main__':
    print("GWT Per-Subshell Z_eff Calculator")
    print("=" * 65)

    for element in ['H', 'He', 'Li', 'Na', 'Fe', 'Ca']:
        atom = ATOMS[element]
        levels = compute_all_levels(element)

        print(f"\n{element} (Z={atom['Z']}):")
        print(f"  {'n,l':>4} {'occ':>5} {'Z_eff':>6} {'Z_net':>6} {'alpha':>6} {'S':>6} {'E(eV)':>10}")
        print(f"  {'-'*50}")

        for (n, l), info in sorted(levels.items(), key=lambda x: x[1]['E']):
            occ = f"[{info['occ']}/{info['max_occ']}]" if info['occ'] > 0 else "     "
            tag = "" if info['type'] == 'occupied' else " *"
            print(f"  {n}{l_names[l]:1s} {occ:>5} {info['Z_eff']:6.2f} {info['Z_net']:6.2f} "
                  f"{info['alpha']:6.4f} {info['S']:6.2f} {info['E']:10.3f}{tag}")

    # Sodium D-line check
    print(f"\n{'='*65}")
    print("SODIUM D-LINE CHECK")
    print(f"{'='*65}")
    levels_Na = compute_all_levels('Na')
    E_3s = levels_Na[(3, 0)]['E']
    E_3p = levels_Na[(3, 1)]['E']
    dE = E_3p - E_3s
    lam = 1240.0 / dE if dE > 0 else float('inf')
    print(f"  3s: E = {E_3s:.4f} eV, Z_eff = {levels_Na[(3,0)]['Z_eff']:.4f}")
    print(f"  3p: E = {E_3p:.4f} eV, Z_eff = {levels_Na[(3,1)]['Z_eff']:.4f}")
    print(f"  dE(3p-3s) = {dE:.4f} eV")
    print(f"  Lambda = {lam:.1f} nm")
    print(f"  Observed Na D: 589.0 nm (2.104 eV)")
    if dE > 0:
        print(f"  Error: {(lam - 589.0)/589.0*100:+.1f}%")

    # Hydrogen Balmer check
    print(f"\n{'='*65}")
    print("HYDROGEN BALMER CHECK")
    print(f"{'='*65}")
    levels_H = compute_all_levels('H', max_n=7)
    for n_hi in range(3, 7):
        E_hi = levels_H[(n_hi, 1)]['E']  # p state
        E_lo = levels_H[(2, 0)]['E']     # 2s state
        dE = E_hi - E_lo
        lam = 1240.0 / dE if dE > 0 else float('inf')
        obs = {3: 656.28, 4: 486.13, 5: 434.05, 6: 410.17}[n_hi]
        err = (lam - obs)/obs*100 if dE > 0 else 999
        print(f"  {n_hi}p->2s: {lam:.2f} nm (obs: {obs:.2f}, err: {err:+.3f}%)")
