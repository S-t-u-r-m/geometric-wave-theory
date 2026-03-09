"""
GWT Aufbau Test: Does wave energy self-consistently find the correct bond config?
==================================================================================

For each molecule with p-orbital bonds, we:
1. Take ONLY: orbital types, R, and total ne_pp (p-channel wave quanta)
2. Enumerate ALL possible distributions of ne_pp across channels:
   - pp_sigma bonding (0, 1, or 2 quanta)
   - pi bonding (0, 1, 2, 3, or 4 quanta — two degenerate channels)
   - pp_sigma antibonding (0, 1, or 2 quanta)
   - pi antibonding (0, 1, 2, 3, or 4 quanta)
3. Compute D_e for each configuration using the V6 formula
4. Pick the configuration that MAXIMIZES bond energy (most stable harmonic)
5. Compare to the known MO theory bond list

If GWT self-consistently recovers the correct filling order,
that's strong evidence the wave picture is fundamental.
"""

import numpy as np
from itertools import product

pi = np.pi
E_H = 13.6057

d = 3
C_bond = pi / d
f_pi = d**2 / (d**2 + 1)
alpha = 1 - f_pi / d
beta = (1 + f_pi) / 2
f_anti = 2*d / (2*d - 1)
c_ionic = 1.0 / (2*d + 1)
overlap_floor = 1.0 / (d + 1)
ionic_threshold = 1.0 / d**3
c_ionic_enhanced = d / (2*d + 1)
phase_ext_power = d - 1
asymm_node_exp = 1.0 / d

Z_eff = {
    'H_1s': 1.0000, 'Li_2s': 1.2792, 'B_2p': 2.4214,
    'C_2p': 3.1358, 'N_2p': 3.8340, 'O_2p': 4.4532,
    'F_2p': 5.0998, 'Na_3s': 2.5074, 'Cl_3p': 6.1161,
}

def get_n(orb): return int(orb.split('_')[1][0])
def get_l(orb): return {'s': 0, 'p': 1}[orb.split('_')[1][1]]
def has_nodes(orb):
    n, l = get_n(orb), get_l(orb)
    return min(n - l - 1, 1)
def orbital_energy(orb):
    n, z = get_n(orb), Z_eff[orb]
    return E_H * (z / n)**2


def compute_De_from_config(R, orb1, orb2, ne, config):
    """Compute bond energy for a given channel configuration.

    config = dict with keys:
      'pp_sigma': (bonding_count, antibonding_count)  -- each 0, 1, or 2
      'pi':       (bonding_count, antibonding_count)  -- each 0-4

    bonding_count = number of wave quanta in bonding channel
    antibonding_count = number in antibonding channel

    For half-filled channels: count can be 0.5, 1, 1.5, 2 etc.
    But wave quanta are integers, so we use integer counts and let
    the formula handle the effective contribution.
    """
    n1, l1 = get_n(orb1), get_l(orb1)
    n2_, l2 = get_n(orb2), get_l(orb2)
    h1, h2 = has_nodes(orb1), has_nodes(orb2)

    a1 = 2 + (1 - 2*l1) * alpha * h1
    a2 = 2 + (1 - 2*l2) * alpha * h2
    b1 = 1 + beta * h1
    b2 = 1 + beta * h2

    E1 = E_H / n1**a1
    E2 = E_H / n2_**a2
    E_scale = np.sqrt(E1 * E2)
    sigma_phase = R / n1**b1 + R / n2_**b2

    z1, z2 = Z_eff[orb1], Z_eff[orb2]
    is_het = (orb1 != orb2)

    # Phase extension for heteronuclear pp bonds
    if l1 == 1 and l2 == 1 and n1 == n2_ and is_het:
        base = 2*np.sqrt(z1*z2)/(z1+z2)
        sigma_phase = sigma_phase / base**phase_ext_power

    sigma_S = abs(np.sin(sigma_phase))
    pi_phase = sigma_phase * f_pi
    pi_S = abs(np.sin(pi_phase))

    # Node corrections
    if (h1 + h2 > 0) and sigma_phase > pi:
        n_lobes = int(np.ceil(sigma_phase / pi))
        if h1 != h2:
            exp = 1.0 + asymm_node_exp
            sigma_S = sigma_S / (n_lobes ** exp)
        else:
            sigma_S = sigma_S / n_lobes
    if sigma_S < overlap_floor:
        sigma_S = overlap_floor

    if (h1 + h2 > 0) and pi_phase > pi:
        n_lobes_pi = int(np.ceil(pi_phase / pi))
        if h1 != h2:
            exp = 1.0 + asymm_node_exp
            pi_S = pi_S / (n_lobes_pi ** exp)
        else:
            pi_S = pi_S / n_lobes_pi
    if pi_S < overlap_floor:
        pi_S = overlap_floor

    # Build bond list from config
    sig_bond = config.get('pp_sigma_bond', 0)
    sig_anti = config.get('pp_sigma_anti', 0)
    pi_bond = config.get('pi_bond', 0)
    pi_anti = config.get('pi_anti', 0)

    # Compute D_cov
    # Each quantum in a bonding channel contributes positively
    # Each quantum in an antibonding channel contributes negatively

    # For sigma: each quantum contributes C_bond * E_scale * S
    # A channel holds max 2 quanta. If 1 quantum, it's "half-filled" = count 0.5?
    # No — in our framework, count is the number of CHANNELS, not quanta.
    # A sigma channel with 2 quanta = count 1 (one filled channel)
    # A sigma channel with 1 quantum = count 0.5 (half-filled channel)
    # A sigma channel with 0 quanta = count 0 (empty)

    # For pi: TWO degenerate channels, each holds 2 quanta (total 4)
    # 4 quanta = count 2 (two filled channels)
    # 3 quanta = count 1.5
    # 2 quanta = count 1
    # 1 quantum = count 0.5
    # 0 quanta = count 0

    sig_bond_count = sig_bond / 2.0    # 0, 0.5, or 1
    sig_anti_count = sig_anti / 2.0
    pi_bond_count = pi_bond / 2.0      # 0, 0.5, 1, 1.5, or 2
    pi_anti_count = pi_anti / 2.0

    # Determine if antibonding is "full" (affects f_anti)
    is_full_anti = (pi_anti_count >= pi_bond_count) if pi_bond_count > 0 else True

    D_cov = 0.0

    # Sigma bonding
    D_cov += sig_bond_count * C_bond * E_scale * sigma_S

    # Pi bonding
    D_cov += pi_bond_count * C_bond * E_scale * pi_S

    # Sigma antibonding
    D_cov -= sig_anti_count * C_bond * E_scale * sigma_S

    # Pi antibonding
    f_a = 1.0 if is_full_anti else f_anti
    D_cov -= pi_anti_count * f_a * C_bond * E_scale * pi_S

    # Ionic correction
    eps1, eps2 = orbital_energy(orb1), orbital_energy(orb2)
    delta_eps = abs(eps1 - eps2)
    V_cov = max(abs(D_cov), 0.01)
    q = delta_eps / np.sqrt(delta_eps**2 + (2*V_cov)**2) if delta_eps > 0 else 0
    ratio = abs(D_cov) / delta_eps if delta_eps > 0.01 else float('inf')
    c_ion = c_ionic_enhanced if ratio < ionic_threshold else c_ionic
    D_ion = c_ion * q**2 * 2 * E_H / R

    return D_cov + D_ion


def config_to_bondlist(config):
    """Convert config dict to readable bond list string."""
    parts = []
    sig_b = config.get('pp_sigma_bond', 0)
    sig_a = config.get('pp_sigma_anti', 0)
    pi_b = config.get('pi_bond', 0)
    pi_a = config.get('pi_anti', 0)

    if sig_b > 0:
        parts.append(f"sig({sig_b})")
    if pi_b > 0:
        parts.append(f"pi({pi_b})")
    if sig_a > 0:
        parts.append(f"sig*({sig_a})")
    if pi_a > 0:
        parts.append(f"pi*({pi_a})")

    bo = (sig_b + pi_b - sig_a - pi_a) / 2.0
    return ' + '.join(parts) if parts else 'empty', bo


def config_to_mo_string(config):
    """Convert to MO-theory-like notation."""
    sig_b = config.get('pp_sigma_bond', 0)
    pi_b = config.get('pi_bond', 0)
    sig_a = config.get('pp_sigma_anti', 0)
    pi_a = config.get('pi_anti', 0)

    parts = []
    if sig_b > 0: parts.append(f"sigma^{sig_b}")
    if pi_b > 0: parts.append(f"pi^{pi_b}")
    if sig_a > 0: parts.append(f"sigma*^{sig_a}")
    if pi_a > 0: parts.append(f"pi*^{pi_a}")
    return ' '.join(parts) if parts else '(empty)'


# =============================================================================
# PP-BOND MOLECULES TO TEST
# =============================================================================
# Only molecules where both orbitals are p-type (pp bonds)
pp_molecules = [
    ('B2',  3.005, 3.02,  'B_2p',  'B_2p',  10, [('pi', 2)]),
    ('C2',  2.348, 6.32,  'C_2p',  'C_2p',  12, [('pi', 2), ('sp_sigma', 1), ('sp_sigma_anti', 1)]),
    ('N2',  2.074, 9.759, 'N_2p',  'N_2p',  14, [('pp_sigma', 1), ('pi', 2)]),
    ('O2',  2.282, 5.213, 'O_2p',  'O_2p',  16, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)]),
    ('F2',  2.668, 1.660, 'F_2p',  'F_2p',  18, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)]),
    ('Cl2', 3.757, 2.514, 'Cl_3p', 'Cl_3p', 34, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 2)]),
    ('CO',  2.132, 11.225,'C_2p',  'O_2p',  14, [('pp_sigma', 1), ('pi', 2)]),
    ('NO',  2.175, 6.497, 'N_2p',  'O_2p',  15, [('pp_sigma', 1), ('pi', 2), ('pi_anti', 1)]),
    ('BF',  2.386, 7.81,  'B_2p',  'F_2p',  14, [('pp_sigma', 1), ('pi', 2)]),
    ('CN',  2.214, 7.738, 'C_2p',  'N_2p',  13, [('pp_sigma', 1), ('pi', 2)]),
]


# =============================================================================
# MAIN: For each molecule, find the most stable configuration
# =============================================================================
print("=" * 90)
print("  GWT AUFBAU TEST: Does wave energy find the correct bond configuration?")
print("  Input: orbital types + R + ne_pp (p-channel wave quanta)")
print("  Output: configuration that maximizes bond energy")
print("=" * 90)

for name, R, De_exp, orb1, orb2, ne, known_bonds in pp_molecules:
    l1, l2 = get_l(orb1), get_l(orb2)
    n1, n2 = get_n(orb1), get_n(orb2)

    # Determine core electrons
    if n1 == 2 and n2 == 2:
        core = 8
    elif n1 == 3 and n2 == 3:
        core = 24
    elif n1 == 2 and n2 == 3:
        core = 16  # approximate
    elif n1 == 3 and n2 == 2:
        core = 16
    else:
        core = 8

    ne_pp = ne - core

    # Generate all valid configurations
    # Filling order: quanta go into channels
    # Max per channel: sigma_bond=2, pi_bond=4, sigma_anti=2, pi_anti=4
    # Total must equal ne_pp
    # Constraint: antibonding can't exceed bonding in same type
    #   (can't have sigma* without sigma filled first)

    best_De = -1e10
    best_config = None
    all_configs = []

    for sig_b in range(0, min(ne_pp, 2) + 1):
        for pi_b in range(0, min(ne_pp - sig_b, 4) + 1):
            for sig_a in range(0, min(ne_pp - sig_b - pi_b, 2) + 1):
                pi_a = ne_pp - sig_b - pi_b - sig_a
                if pi_a < 0 or pi_a > 4:
                    continue
                # Physical constraint: can't have more antibonding than bonding
                # (antibonding channel only exists if bonding is occupied)
                if sig_a > sig_b:
                    continue
                if pi_a > pi_b:
                    continue

                config = {
                    'pp_sigma_bond': sig_b,
                    'pi_bond': pi_b,
                    'pp_sigma_anti': sig_a,
                    'pi_anti': pi_a,
                }

                De = compute_De_from_config(R, orb1, orb2, ne, config)
                all_configs.append((De, config))

                if De > best_De:
                    best_De = De
                    best_config = config

    # Also find 2nd and 3rd best
    all_configs.sort(key=lambda x: -x[0])

    # Convert known bonds to config for comparison
    known_config = {'pp_sigma_bond': 0, 'pi_bond': 0, 'pp_sigma_anti': 0, 'pi_anti': 0}
    for btype, count in known_bonds:
        if btype == 'pp_sigma':
            known_config['pp_sigma_bond'] = count * 2  # count=1 means 2 quanta
        elif btype == 'pi':
            known_config['pi_bond'] = count * 2  # count=2 means 4 quanta
        elif btype == 'pi_anti':
            known_config['pi_anti'] = count * 2
        elif btype == 'sp_sigma':
            known_config['pp_sigma_bond'] = count * 2
        elif btype == 'sp_sigma_anti':
            known_config['pp_sigma_anti'] = count * 2

    best_str, best_bo = config_to_bondlist(best_config)
    known_str, known_bo = config_to_bondlist(known_config)

    match = (best_config == known_config)

    print(f"\n  {name} (ne_pp={ne_pp}, De_exp={De_exp:.3f} eV):")
    print(f"    Known (MO theory): {config_to_mo_string(known_config):30s}  BO={known_bo:.1f}")
    print(f"    GWT best:          {config_to_mo_string(best_config):30s}  BO={best_bo:.1f}  "
          f"De={best_De:.3f} eV  {'MATCH' if match else 'MISMATCH'}")

    if not match:
        # Show what the known config gives
        known_De = compute_De_from_config(R, orb1, orb2, ne, known_config)
        print(f"    Known config De:   {known_De:.3f} eV")
        print(f"    Difference:        {best_De - known_De:.3f} eV "
              f"({'GWT prefers its own' if best_De > known_De else 'MO theory is more stable'})")

    # Show top 3 configs
    print(f"    Top 3 configurations:")
    for rank, (De_val, cfg) in enumerate(all_configs[:3]):
        cfg_str = config_to_mo_string(cfg)
        _, bo = config_to_bondlist(cfg)
        marker = " <-- KNOWN" if cfg == known_config else ""
        marker2 = " <-- GWT BEST" if cfg == best_config else ""
        print(f"      {rank+1}. {cfg_str:30s}  BO={bo:.1f}  De={De_val:.3f} eV{marker}{marker2}")


# =============================================================================
# SUMMARY
# =============================================================================
print(f"\n{'='*90}")
print(f"  SUMMARY")
print(f"{'='*90}")

matches = 0
mismatches = 0
for name, R, De_exp, orb1, orb2, ne, known_bonds in pp_molecules:
    l1, l2 = get_l(orb1), get_l(orb2)
    n1, n2 = get_n(orb1), get_n(orb2)
    if n1 == 2 and n2 == 2: core = 8
    elif n1 == 3 and n2 == 3: core = 24
    else: core = 16
    ne_pp = ne - core

    known_config = {'pp_sigma_bond': 0, 'pi_bond': 0, 'pp_sigma_anti': 0, 'pi_anti': 0}
    for btype, count in known_bonds:
        if btype == 'pp_sigma': known_config['pp_sigma_bond'] = count * 2
        elif btype == 'pi': known_config['pi_bond'] = count * 2
        elif btype == 'pi_anti': known_config['pi_anti'] = count * 2
        elif btype == 'sp_sigma': known_config['pp_sigma_bond'] = count * 2
        elif btype == 'sp_sigma_anti': known_config['pp_sigma_anti'] = count * 2

    best_De = -1e10
    best_config = None
    for sig_b in range(0, min(ne_pp, 2) + 1):
        for pi_b in range(0, min(ne_pp - sig_b, 4) + 1):
            for sig_a in range(0, min(ne_pp - sig_b - pi_b, 2) + 1):
                pi_a = ne_pp - sig_b - pi_b - sig_a
                if pi_a < 0 or pi_a > 4: continue
                if sig_a > sig_b: continue
                if pi_a > pi_b: continue
                config = {
                    'pp_sigma_bond': sig_b, 'pi_bond': pi_b,
                    'pp_sigma_anti': sig_a, 'pi_anti': pi_a,
                }
                De = compute_De_from_config(R, orb1, orb2, ne, config)
                if De > best_De:
                    best_De = De
                    best_config = config

    if best_config == known_config:
        matches += 1
        print(f"  {name:5s}: MATCH")
    else:
        mismatches += 1
        _, best_bo = config_to_bondlist(best_config)
        _, known_bo = config_to_bondlist(known_config)
        print(f"  {name:5s}: MISMATCH  (GWT: BO={best_bo:.1f}, MO: BO={known_bo:.1f})")

print(f"\n  Result: {matches}/{matches+mismatches} match MO theory filling order")
print(f"  {'ALL MATCH' if mismatches == 0 else f'{mismatches} mismatches'}")
