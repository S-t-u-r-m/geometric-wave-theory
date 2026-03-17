#!/usr/bin/env python3
"""
GWT Force Balance Simulator v2 — IE directly from equilibrium.

Uses FULL v19 angular coupling (including anti-screening from d/f shells).
Position-dependent Z_eff(r) sampled across radial width gives effective alpha.

NO explicit alpha formula. Energy from force balance.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial, comb
import time

d = 3
gamma_gwt = np.pi / (2**(d+1)*np.pi - 2)
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6
w_pi = np.cos(np.pi / d)
w_delta = np.cos(2*np.pi / d)
n_ref = d + 1
W_F = {'s': w_delta / d, 'p': w_pi, 'd': w_delta / d}

# Import v19 functions and atom data
code = open('calculations/z_eff_v19.py').read()
setup = code.split('# === RUN AND DISPLAY ===')[0]
g = {}
exec(setup, g)
sys.stdout = io.TextIOWrapper(open(1, 'wb', closefd=False), encoding='utf-8', errors='replace')

atoms = g['atoms']
core_screening = g['core_screening']
d_screen_t2g = g['d_screen_t2g']

v19_results = g['calc_all']()
v19_lookup = {r[0]: r for r in v19_results}

def get_atom_info(Z, config):
    """Extract valence info (same logic as v19)."""
    val_n = max(nn for nn, ll, c in config)
    active = [ll for nn, ll, c in config if nn == val_n and c > 0]
    val_l = max(active) if active else 0
    s_count = sum(c for nn, ll, c in config if nn == val_n and ll == 0)
    p_count = sum(c for nn, ll, c in config if nn == val_n and ll == 1)
    dc_val = sum(c for nn, ll, c in config if nn == val_n - 1 and ll == 2)
    is_pd = (s_count == 0 and dc_val == 10)
    n = val_n - 1 if is_pd else val_n
    has_p = val_l >= 1 and p_count > 0 and not is_pd
    if has_p: val_block = 'p'
    elif dc_val > 0 or is_pd: val_block = 'd'
    else: val_block = 's'
    return val_n, val_l, n, s_count, p_count, dc_val, is_pd, has_p, val_block

def solve_radii(Z, config):
    """Self-consistent radii with l-dependent centrifugal barrier."""
    shells = [(nn, ll, c) for nn, ll, c in config if c > 0]
    if not shells:
        return {}
    r = {}
    for nn, ll, c in shells:
        Z_g = max(Z * 0.3, 1.0)
        n_eff = nn + ll * (ll + 1) / (2 * nn * Z_g)
        r[(nn, ll)] = n_eff**2 / Z_g

    for _ in range(500):
        r_old = dict(r)
        for nn_i, ll_i, ne_i in shells:
            S = 0.0
            for nn_j, ll_j, ne_j in shells:
                if (nn_j, ll_j) == (nn_i, ll_i) or nn_j == nn_i:
                    continue
                n_ch = min(ne_j, 2*ll_j + 1)
                ri = r[(nn_i, ll_i)]
                rj = r[(nn_j, ll_j)]
                width = 0.15 * max(ri, 0.01)
                overlap = 1.0 / (1.0 + np.exp((rj - ri) / max(width, 0.001)))
                S += n_ch * w_pi * overlap
            Z_eff = max(Z - S, 0.5)
            n_eff = nn_i + ll_i * (ll_i + 1) / (2 * nn_i * Z_eff)
            r[(nn_i, ll_i)] = n_eff**2 / Z_eff
        if max(abs(r[k] - r_old[k]) / (r_old[k] + 1e-10) for k in r) < 1e-9:
            break
    return r

def v19_S_core_per_shell(config, core_ref_n, is_pd, val_l, val_block, dc_val):
    """Get per-shell contributions to S_core using v19 logic.
    Returns list of (nn, ll, count, S_contribution)."""
    w_f = W_F[val_block]
    contributions = []
    for nn, ll, count in config:
        S_shell = 0.0
        if nn < core_ref_n:
            if ll <= 1:
                S_shell = min(count, 2*ll+1) * w_pi
            elif ll == 2:
                dn = core_ref_n - nn
                if dn >= 2 and count == 10:
                    Sa = d_screen_t2g(count, val_l)
                    Sn = 5 * w_pi
                    S_shell = Sn + (Sa - Sn) / d
                else:
                    S_shell = d_screen_t2g(count, val_l)
            elif ll == 3:
                S_shell = min(count, 7) * w_f
        elif nn == core_ref_n and is_pd and ll != 2:
            S_shell = min(count, 2*ll+1) * w_pi
        if abs(S_shell) > 1e-10:
            contributions.append((nn, ll, count, S_shell))
    # Pd d-shell self-screening
    if is_pd:
        S_pd = d_screen_t2g(dc_val - 1, 2)
        contributions.append((core_ref_n, 2, dc_val, S_pd))
    return contributions

def Z_eff_at_position(Z, r_probe, shell_contributions, radii, r_val_default):
    """Position-dependent Z_eff using v19 angular weights + spatial overlap.

    Each core shell's contribution is weighted by how much of it
    is inside r_probe. Anti-screening shells (negative S) are also
    position-dependent — they push LESS when far away.
    """
    S = 0.0
    for nn, ll, count, S_full in shell_contributions:
        key = (nn, ll)
        if key in radii:
            r_shell = radii[key]
        else:
            r_shell = nn**2 / max(Z * 0.3, 1.0)

        # Overlap: fraction of this shell inside r_probe
        width = max(r_shell / (2 * max(nn, 1)), 0.01)
        overlap = 1.0 / (1.0 + np.exp((r_shell - r_probe) / max(width, 0.001)))

        S += S_full * overlap

    return Z - S

def solve_atom(Z, config):
    """Find IE from force balance with position-dependent Z_eff."""
    val_n, val_l, n, s_count, p_count, dc_val, is_pd, has_p, val_block = get_atom_info(Z, config)

    # Get per-shell screening contributions (v19 angular coupling)
    shell_S = v19_S_core_per_shell(config, n, is_pd, val_l, val_block, dc_val)

    # Self-consistent radii
    radii = solve_radii(Z, config)

    # Valence shell key
    if is_pd:
        val_key = (val_n - 1, 2)
    elif has_p:
        val_key = (val_n, 1)
    else:
        val_key = (val_n, 0)
    if val_key not in radii:
        val_key = max(radii.keys(), key=lambda k: (k[0], k[1]))

    r_eq = radii[val_key]

    # Position-dependent Z_eff at equilibrium and across radial width
    Z_eq = Z_eff_at_position(Z, r_eq, shell_S, radii, r_eq)

    # Sample Z_eff across the mode's radial extent
    # Width from quantum uncertainty: delta_r ~ r_eq / (2*n)
    width = r_eq / (2 * max(n, 1))
    n_samples = 5
    r_samples = np.linspace(max(r_eq - width, 0.01), r_eq + width, n_samples)

    # Weight by |psi|^2 ~ r^(2l) * exp(-2r/r_eq) (hydrogen-like envelope)
    Z_samples = []
    weights = []
    for r_s in r_samples:
        Z_s = Z_eff_at_position(Z, r_s, shell_S, radii, r_eq)
        Z_samples.append(Z_s)
        # Weight: probability density ~ r^2l * exp(-2*Z_eff*r/n)
        # Simplified: 1/r^2 weighting (inner regions have more weight)
        w = 1.0 / max(r_s, 0.01)**2
        weights.append(w)

    Z_samples = np.array(Z_samples)
    weights = np.array(weights)
    Z_eff_avg = np.sum(Z_samples * weights) / np.sum(weights)

    # IE from averaged Z_eff
    # E = (Z_eff_avg / n)^2 * E_H
    # But Z_eff_avg is the position-averaged effective charge,
    # which naturally accounts for the "alpha < 1" effect
    IE = (max(Z_eff_avg, 0.01) / n)**2 * E_H

    # Lattice shear
    N_core = sum(c for nn, ll, c in config if nn < n)
    IE *= (1 + gamma_gwt * N_core / (d**2 * n**2))

    return IE, r_eq, Z_eq, Z_eff_avg

# === RUN ===
print(f"GWT Force Balance v2 — Full angular coupling")
print(f"{'='*80}")

t0 = time.time()
results = []

for Z, sym, E_obs, config in atoms:
    IE_sim, r_eq, Z_eq, Z_avg = solve_atom(Z, config)
    v19 = v19_lookup[Z]
    err_sim = (IE_sim - E_obs) / E_obs * 100
    err_v19 = v19[4]

    # Implied alpha from Z_avg
    _, _, n, _, _, dc_val, is_pd, _, val_block = get_atom_info(Z, config)
    S_v19 = core_screening(config, n, is_pd,
              max([ll for nn, ll, c in config if nn == max(nn2 for nn2,_,_ in config) and c > 0], default=0),
              val_block)
    if is_pd: S_v19 += d_screen_t2g(dc_val - 1, 2)
    Z_net = Z - S_v19
    a_impl = np.log(max(Z_avg, 0.01)) / np.log(max(Z_net, 1.01)) if Z_net > 1 else 0

    results.append((Z, sym, E_obs, IE_sim, err_sim, err_v19, Z_avg, a_impl, v19[5]))

elapsed = time.time() - t0

# Display
print(f"\n  {'Z':>3} {'Sym':<3} {'obs':>7} {'sim':>7} {'sim%':>7} {'v19%':>7} {'Z_avg':>6} {'a_imp':>6} {'a_v19':>6}")
print(f"  {'-'*65}")
n_better = 0
for r in results:
    Z, sym, obs, sim, es, ev, za, ai, av = r
    if abs(es) < abs(ev): n_better += 1
    mark = ' <' if abs(es) < abs(ev) and abs(ev - es) > 0.5 else ''
    if Z in [3,11,19,21,31,37,39,49,55,57,72,81,87,89]:
        print()
    flag = ' ***' if abs(es) > 10 else (' **' if abs(es) > 5 else '')
    print(f"  {Z:3d} {sym:<3} {obs:7.3f} {sim:7.3f} {es:+7.1f}{flag} {ev:+7.1f} {za:6.2f} {ai:6.4f} {av:6.4f}{mark}")

# Stats
sm = np.mean([abs(r[4]) for r in results])
vm = np.mean([abs(r[5]) for r in results])
print(f"\n{'='*80}")
print(f"  Force balance: {sm:.2f}%")
print(f"  v19+shear:     {vm:.2f}%")
print(f"  Sim better: {n_better}/{len(results)} atoms")
print(f"  Sim <5%: {sum(1 for r in results if abs(r[4]) < 5)}/{len(results)}")
print(f"  Sim <10%: {sum(1 for r in results if abs(r[4]) < 10)}/{len(results)}")
print(f"  Time: {elapsed*1000:.0f}ms ({elapsed/len(results)*1000:.1f}ms/atom)")
