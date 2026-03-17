#!/usr/bin/env python3
"""
GWT Wave Force Simulator v0.3 — Full periodic table.

Hybrid: simulation finds Z_net, v19 formula computes alpha.
Self-consistent radial force balance with GWT coupling weights.
Deep-core only screening (same-period handled by alpha).
103 atoms, ~0ms per atom.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial, comb
import time

d = 3
gamma = np.pi / (2**(d+1)*np.pi - 2)
S_tunnel = 2**(2*d+1) / np.pi**2
BZ = np.log(2*d)
exponent = (2/factorial(d)) * (S_tunnel + BZ)
alpha_em = np.exp(-exponent)
E_H = (alpha_em**2 / 2) * 0.51100e6
w_pi = np.cos(np.pi / d)
w_delta = np.cos(2*np.pi / d)
J_hund = 2 / (d + 2)
n_ref = d + 1
W_F = {'s': w_delta / d, 'p': w_pi, 'd': w_delta / d}

# Import atom data and v19 calc from z_eff_v19.py
# We read the file and extract what we need
v19_code = open('calculations/z_eff_v19.py').read()
v19_setup = v19_code.split('# === RUN AND DISPLAY ===')[0]

# Execute v19 setup to get atoms list and calc_all function
# Save and restore stdout around exec (v19 wraps stdout and closes it)
_saved_stdout = sys.stdout
_saved_buffer = sys.stdout.buffer
exec_globals = {}
exec(v19_setup, exec_globals)
sys.stdout = io.TextIOWrapper(open(1, 'wb', closefd=False), encoding='utf-8', errors='replace')
atoms_v19 = exec_globals['atoms']
calc_all_v19 = exec_globals['calc_all']

# Run v19 formula for all atoms
v19_results = calc_all_v19()
v19_lookup = {r[0]: r for r in v19_results}  # Z -> (Z, sym, E_obs, E_pred, err, alpha)

# ============================================================
# SIMULATION: Self-consistent Z_net from wave forces
# ============================================================

def solve_Znet_sim(Z, config, val_block='s', val_l=0, core_ref_n=None, is_pd=False):
    """Self-consistent Z_net using v19 angular coupling + sim radial overlap.

    Uses the FULL v19 core_screening logic for angular weights (t2g/eg,
    valence-dependent f, deep d10 blending) but modulates each shell's
    contribution by the self-consistent radial overlap factor.
    """
    if not config:
        return float(Z), 0.0

    val_n = max(nn for nn, ll, c in config)
    if core_ref_n is None:
        core_ref_n = val_n

    # All shells for radial solve
    all_shells = [(nn, ll, c) for nn, ll, c in config if c > 0]
    if not all_shells:
        return float(Z), 0.0

    # Initialize radii (hydrogen-like)
    r = {}
    for nn, ll, c in all_shells:
        key = (nn, ll)
        r[key] = nn**2 / max(Z * 0.3, 1.0)

    # Self-consistent radial positions
    # Use simple w_pi screening for radial solve (fast convergence)
    for iteration in range(300):
        r_old = dict(r)
        for nn_i, ll_i, ne_i in all_shells:
            key_i = (nn_i, ll_i)
            S = 0.0
            for nn_j, ll_j, ne_j in all_shells:
                key_j = (nn_j, ll_j)
                if key_j == key_i or nn_j == nn_i:
                    continue
                n_ch = min(ne_j, 2*ll_j + 1)
                r_j, r_i = r[key_j], r[key_i]
                sigma = 0.15 * max(r_i, 0.01)
                overlap = 1.0 / (1.0 + np.exp((r_j - r_i) / max(sigma, 0.001)))
                S += n_ch * w_pi * overlap  # w_pi for radial solve
            r[key_i] = nn_i**2 / max(Z - S, 0.5)

        max_dr = max(abs(r[k] - r_old[k]) / (r_old[k] + 1e-10) for k in r)
        if max_dr < 1e-8:
            break

    # Now compute S_core using v19 angular logic + sim overlap factors
    # Get v19's per-shell screening contributions
    v19_S = exec_globals['core_screening'](config, core_ref_n, is_pd, val_l, val_block)
    if is_pd:
        dc_val = sum(c for nn, ll, c in config if nn == val_n - 1 and ll == 2)
        v19_S += exec_globals['d_screen_t2g'](dc_val - 1, 2)

    # Compute per-shell v19 contributions and their overlap factors
    # We compute the overlap correction: how much does the sim's spatial
    # arrangement differ from the formula's assumption (overlap = 1)?
    val_shells = [(nn, ll) for nn, ll, c in config if nn == core_ref_n or
                  (nn == val_n and nn >= core_ref_n)]
    val_key = max(((nn, ll) for nn, ll, c in all_shells if nn >= core_ref_n),
                  key=lambda x: (x[0], x[1]), default=(val_n, 0))

    if val_key in r:
        r_val = r[val_key]
    else:
        r_val = val_n**2 / max(Z * 0.3, 1.0)

    # Average overlap factor for core shells
    total_weight = 0.0
    weighted_overlap = 0.0
    for nn_j, ll_j, ne_j in all_shells:
        if nn_j >= core_ref_n:
            continue
        key_j = (nn_j, ll_j)
        if key_j not in r:
            continue
        r_j = r[key_j]
        sigma = 0.15 * max(r_val, 0.01)
        overlap = 1.0 / (1.0 + np.exp((r_j - r_val) / max(sigma, 0.001)))
        n_ch = min(ne_j, 2*ll_j + 1)
        total_weight += n_ch
        weighted_overlap += n_ch * overlap

    if total_weight > 0:
        avg_overlap = weighted_overlap / total_weight
    else:
        avg_overlap = 1.0

    # Apply overlap correction to v19's S_core
    S_core_sim = v19_S * avg_overlap
    Z_net_sim = Z - S_core_sim
    return Z_net_sim, S_core_sim

# ============================================================
# RUN ALL 103 ATOMS
# ============================================================
print(f"GWT Wave Force Simulator v0.3 — {len(atoms_v19)} atoms")
print(f"{'='*85}")

t0 = time.time()
results = []

for Z, sym, E_obs, config in atoms_v19:
    # V19 formula result
    v19 = v19_lookup[Z]
    E_v19, err_v19, alpha_v19 = v19[3], v19[4], v19[5]

    # Determine valence block (same logic as v19)
    val_n = max(nn for nn, ll, c in config)
    active = [ll for nn, ll, c in config if nn == val_n and c > 0]
    val_l = max(active) if active else 0
    s_count = sum(c for nn, ll, c in config if nn == val_n and ll == 0)
    p_count = sum(c for nn, ll, c in config if nn == val_n and ll == 1)
    dc_val = sum(c for nn, ll, c in config if nn == val_n - 1 and ll == 2)
    is_pd = (s_count == 0 and dc_val == 10)
    n_used_local = val_n - 1 if is_pd else val_n
    has_p = val_l >= 1 and p_count > 0 and not is_pd
    if has_p: val_block = 'p'
    elif dc_val > 0 or is_pd: val_block = 'd'
    else: val_block = 's'

    # Simulation Z_net with v19 angular coupling
    Znet_sim, Score_sim = solve_Znet_sim(Z, config, val_block, val_l, n_used_local, is_pd)

    # V19's Z_net for comparison
    Znet_v19 = Z - (Z - v19[3]**(1/v19[5]) * max(nn for nn,_,_ in config) / E_H**0.5)
    # Actually easier to get from the v19 internals... let me just use the formula
    # Z_net_v19 is implicit. Use: E_pred = (Z_net^alpha / n)^2 * E_H
    # Z_net^alpha = n * sqrt(E_pred / E_H)
    n_val = max(nn for nn, ll, c in config)
    # For Pd: n is val_n - 1
    s_count = sum(c for nn, ll, c in config if nn == n_val and ll == 0)
    dc_val = sum(c for nn, ll, c in config if nn == n_val - 1 and ll == 2)
    is_pd = (s_count == 0 and dc_val == 10)
    n_used = n_val - 1 if is_pd else n_val

    # Hybrid: sim Z_net + v19 alpha
    Z_eff_hyb = max(Znet_sim, 0.01) ** alpha_v19
    E_hybrid = (Z_eff_hyb / n_used)**2 * E_H
    err_hyb = (E_hybrid - E_obs) / E_obs * 100

    results.append((Z, sym, E_obs, E_v19, err_v19, E_hybrid, err_hyb, alpha_v19))

elapsed = time.time() - t0

# Display
print(f"\n  {'Z':>3} {'Sym':<3} {'E_obs':>7} {'v19':>7} {'v19%':>7} {'hybrid':>7} {'hyb%':>7} {'better':>6}")
print(f"  {'-'*55}")
n_better = 0
n_worse = 0
for r in results:
    Z, sym, obs, v19, ev, hyb, eh, alpha = r
    better = abs(eh) < abs(ev)
    if better: n_better += 1
    else: n_worse += 1
    mark = '  <' if better and abs(ev - eh) > 0.1 else ''
    if Z in [3,11,19,21,31,37,39,49,55,57,72,81,87,89]:
        print()
    flag = ' ***' if abs(eh) > 10 else (' **' if abs(eh) > 5 else '')
    print(f"  {Z:3d} {sym:<3} {obs:7.3f} {v19:7.3f} {ev:+7.1f} {hyb:7.3f} {eh:+7.1f}{flag}{mark}")

# Stats
v19_all = np.mean([abs(r[4]) for r in results])
hyb_all = np.mean([abs(r[6]) for r in results])
v19_p15 = np.mean([abs(r[4]) for r in results if r[0] <= 54])
hyb_p15 = np.mean([abs(r[6]) for r in results if r[0] <= 54])
v19_p6 = np.mean([abs(r[4]) for r in results if 55 <= r[0] <= 86])
hyb_p6 = np.mean([abs(r[6]) for r in results if 55 <= r[0] <= 86])
v19_p7 = np.mean([abs(r[4]) for r in results if r[0] >= 87])
hyb_p7 = np.mean([abs(r[6]) for r in results if r[0] >= 87])

print(f"\n{'='*85}")
print(f"  SUMMARY: {len(results)} atoms, {elapsed*1000:.1f}ms total ({elapsed/len(results)*1000:.2f}ms/atom)")
print(f"{'='*85}")
print(f"  {'':>15} {'v19':>8} {'hybrid':>8} {'winner':>8}")
print(f"  {'ALL':>15} {v19_all:8.2f}% {hyb_all:8.2f}% {'SIM' if hyb_all < v19_all else 'V19':>8}")
print(f"  {'P1-5':>15} {v19_p15:8.2f}% {hyb_p15:8.2f}% {'SIM' if hyb_p15 < v19_p15 else 'V19':>8}")
print(f"  {'P6':>15} {v19_p6:8.2f}% {hyb_p6:8.2f}% {'SIM' if hyb_p6 < v19_p6 else 'V19':>8}")
print(f"  {'P7':>15} {v19_p7:8.2f}% {hyb_p7:8.2f}% {'SIM' if hyb_p7 < v19_p7 else 'V19':>8}")
print(f"\n  Sim better: {n_better}/{len(results)} atoms")
print(f"  v19 outliers >10%: {sum(1 for r in results if abs(r[4]) > 10)}")
print(f"  Hybrid outliers >10%: {sum(1 for r in results if abs(r[6]) > 10)}")
