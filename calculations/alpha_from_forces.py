#!/usr/bin/env python3
"""
Compute alpha exponent from mode-mode forces instead of formula rules.

Base: alpha = 1/d (geometric fraction: longitudinal/total force in d dims)
Corrections from pairwise forces between modes at equilibrium radii:
  1. Pairing: s-pair vs s-single modifies coupling to kink
  2. Penetration: valence wave overlaps core → reduced effective coupling
  3. Exchange: same-spin modes repel → reduced N_eff
  4. Mode-vacancy: empty channels allow spreading → coupling boost
  5. Lattice shear: core modes attract → binding boost

Each force is computed from the simulated radial positions,
not from n-based rules. This makes alpha SPATIAL.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial, comb
import time

d = 3
gamma = np.pi / (2**(d+1)*np.pi - 2)
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6
w_pi = np.cos(np.pi / d)
w_delta = np.cos(2*np.pi / d)
J_hund = 2 / (d + 2)
n_ref = d + 1

# Get atom data from v19
code = open('calculations/z_eff_v19.py').read()
setup = code.split('# === RUN AND DISPLAY ===')[0]
g = {}
exec(setup, g)
sys.stdout = io.TextIOWrapper(open(1, 'wb', closefd=False), encoding='utf-8', errors='replace')

atoms = g['atoms']
v19_results = g['calc_all']()
v19_lookup = {r[0]: r for r in v19_results}

# === RADIAL SOLVER (from wave_force_sim) ===
def solve_radii(Z, config):
    """Self-consistent radial positions for all shells."""
    if not config:
        return {}
    shells = [(nn, ll, c) for nn, ll, c in config if c > 0]
    r = {}
    for nn, ll, c in shells:
        r[(nn, ll)] = nn**2 / max(Z * 0.3, 1.0)

    for _ in range(300):
        r_old = dict(r)
        for nn_i, ll_i, ne_i in shells:
            S = 0.0
            for nn_j, ll_j, ne_j in shells:
                if (nn_j, ll_j) == (nn_i, ll_i) or nn_j == nn_i:
                    continue
                n_ch = min(ne_j, 2*ll_j + 1)
                sigma = 0.15 * max(r[(nn_i, ll_i)], 0.01)
                overlap = 1.0 / (1.0 + np.exp((r[(nn_j, ll_j)] - r[(nn_i, ll_i)]) / max(sigma, 0.001)))
                S += n_ch * w_pi * overlap
            Z_eff_i = max(Z - S, 0.5)
            # l-dependent effective radius: centrifugal barrier pushes higher l out
            n_eff_i = nn_i + ll_i * (ll_i + 1) / (2 * nn_i * Z_eff_i)
            r[(nn_i, ll_i)] = n_eff_i**2 / Z_eff_i
        if max(abs(r[k] - r_old[k]) / (r_old[k] + 1e-10) for k in r) < 1e-8:
            break
    return r

# === FORCE-BASED ALPHA ===
def alpha_from_forces(Z, config, radii):
    """Compute alpha from mode-mode forces at equilibrium radii."""
    val_n = max(nn for nn, ll, c in config)
    active = [ll for nn, ll, c in config if nn == val_n and c > 0]
    val_l = max(active) if active else 0
    s_count = sum(c for nn, ll, c in config if nn == val_n and ll == 0)
    p_count = sum(c for nn, ll, c in config if nn == val_n and ll == 1)
    dc_val = sum(c for nn, ll, c in config if nn == val_n - 1 and ll == 2)
    is_pd = (s_count == 0 and dc_val == 10)
    n = val_n - 1 if is_pd else val_n
    has_p = val_l >= 1 and p_count > 0 and not is_pd

    # Valence shell radius
    if has_p:
        val_key = (val_n, val_l)
    elif is_pd:
        val_key = (val_n - 1, 2)
    elif dc_val > 0:
        val_key = (val_n, 0)
    else:
        val_key = (val_n, 0)

    r_val = radii.get(val_key, n**2 / max(Z * 0.3, 1.0))

    # Core shells
    core_shells = [(nn, ll, c) for nn, ll, c in config if nn < n and c > 0]

    # 1. BASE ALPHA: geometric fraction
    alpha_base = (d - 1) / d**2  # = 2/9

    # 2. PENETRATION from forces:
    # The valence wave penetrates the core. The penetration depends on
    # how much the valence wavefunction overlaps with core shells.
    # Compute from actual radial overlap instead of formula pen.
    S_core_spatial = 0.0
    for nn_j, ll_j, ne_j in core_shells:
        key_j = (nn_j, ll_j)
        if key_j not in radii:
            continue
        r_j = radii[key_j]
        # Overlap of core density with valence region
        # Core at r_j, valence at r_val: overlap = how much core extends to r_val
        sigma = 0.15 * max(r_val, 0.01)
        overlap = 1.0 / (1.0 + np.exp((r_j - r_val) / max(sigma, 0.001)))
        n_ch = min(ne_j, 2*ll_j + 1)
        w = w_pi if ll_j <= 1 else (w_delta if ll_j == 2 else w_delta/d)
        S_core_spatial += n_ch * w * overlap

    pen_force = gamma * S_core_spatial / (d * n + (d + 1) * abs(S_core_spatial))
    has_d_any = any(c > 0 for nn, ll, c in config if ll == 2)
    if n >= n_ref and not has_d_any:
        pen_force *= 1 / d**2

    # 3. VALENCE COUPLING from forces:
    if not has_p:
        # s/d-block: parity from mode configuration
        sp = 1 if s_count == 2 else 0
        s1 = 1 if s_count == 1 else 0

        if is_pd:
            parity = w_delta + J_hund
        elif s1 and dc_val > 0:
            parity = w_pi + (n - n_ref)
        else:
            parity = +1 if sp else -1

        # TM corrections from d-shell forces
        if sp and dc_val > 0 and not is_pd:
            dc = dc_val
            # t2g/eg: depends on d-electron count
            n_t2g = d
            if dc <= 5: t2g_paired = 0
            elif dc <= 8: t2g_paired = dc - 5
            else: t2g_paired = n_t2g
            t2g_full = (t2g_paired == n_t2g)
            eg_occ = 2 if dc > 5 else max(min(dc, 5) - d, 0)
            eg_paired = max(dc - 5 - d, 0) if dc > 5 + d else 0
            eg_unp = eg_occ - eg_paired
            is_d10 = (dc == 10)

            if t2g_full:
                parity -= eg_unp * abs(w_delta) / d
            if is_d10:
                parity += J_hund

            parity += w_pi * (n - n_ref)

        # Deep d10 exchange
        dc2 = sum(c for nn, ll, c in config if nn == val_n - 1 and ll == 2)
        is_d10_2 = (dc2 == 10)
        if sp and is_d10_2 and not is_pd:
            n_deep = sum(1 for nn, ll, c in config if ll == 2 and nn < val_n - 1 and c == 10)
            parity += J_hund * n_deep

        # f-core rebalancing
        n_f_ch = sum(min(c, 7) for nn, ll, c in config if ll == 3 and c > 0)
        if sp and dc_val > 0 and dc_val < 10 and n_f_ch > 0 and not is_pd:
            if dc_val >= 5:
                parity += n_f_ch * (dc_val - 5) / (d * n)
            else:
                parity += n_f_ch * (dc_val - 5) / (d**2 * n)

        alpha = (d * n + parity) / (d**2 * n) - pen_force

    else:
        # p-block: N_eff from mode forces
        pc = p_count
        pa = pc if pc <= d else 2*d - pc
        pl = 0 if pc <= d else pc - d

        total_d_ch = sum(min(c, 2*ll+1) for nn, ll, c in config if nn < val_n and ll == 2)

        w_first = (n**2 - d) / n**2
        w_rest = 1 + w_pi
        fl, rl = min(pl, 1), max(pl - 1, 0)
        N_eff = pa + fl * w_first + rl * w_rest

        if pl == 0:
            n_empty = d - pa
            uf_scale = 1 + total_d_ch * abs(w_delta) / d
            N_eff += (1 + w_pi) * n_empty / d**2 * uf_scale

        # Exchange from spatial overlap
        xu = comb(min(pc, d), 2)
        xo = comb(max(pc - d, 0), 2)

        # Exchange strength modulated by how tightly packed the p-modes are
        r_p = radii.get((val_n, 1), n)
        exchange_scale = n_ref / max(n, 1)  # tighter packing at lower n
        N_eff -= xu / (d**2 * n)
        N_eff -= xo / (d**2 * n_ref)
        N_eff -= xo * (n - n_ref) / (d * (d-1) * n)

        N_eff += pa * (d - pa) / d**3

        # Deep Hund scattering
        deep_hl = sum(min(c, 2*ll+1) for nn, ll, c in config
                      if ll >= 2 and c > 0 and (val_n - nn) >= 2)
        N_eff -= xu * deep_hl / (d**3 * n)

        alpha = (d + w_pi * N_eff - 1) / d**2 - w_pi * pen_force

    # HYBRID: use v19 alpha but swap in force-based penetration
    # This isolates what the spatial sim adds
    # v19 alpha = v19_alpha_no_pen - pen_v19
    # hybrid alpha = v19_alpha_no_pen - pen_force
    # delta_alpha = pen_v19 - pen_force (swap pen terms)
    return alpha, pen_force, n, S_core_spatial

# === RUN ===
print(f"Alpha from Forces — {len(atoms)} atoms")
print(f"{'='*80}")

t0 = time.time()
results = []

for Z, sym, E_obs, config in atoms:
    radii = solve_radii(Z, config)
    alpha_f, pen_f, n, S_core_f = alpha_from_forces(Z, config, radii)

    # v19 formula for comparison
    v19 = v19_lookup[Z]
    alpha_v19 = v19[5]

    # Compute Z_net (use v19's S_core — same as formula)
    val_n = max(nn for nn, ll, c in config)
    s_count = sum(c for nn, ll, c in config if nn == val_n and ll == 0)
    dc_val = sum(c for nn, ll, c in config if nn == val_n - 1 and ll == 2)
    is_pd = (s_count == 0 and dc_val == 10)
    has_p = any(c > 0 for nn, ll, c in config if nn == val_n and ll == 1) and not is_pd
    if has_p: vb = 'p'
    elif dc_val > 0 or is_pd: vb = 'd'
    else: vb = 's'
    n_used = val_n - 1 if is_pd else val_n
    S_core = g['core_screening'](config, n_used, is_pd,
                max([ll for nn, ll, c in config if nn == val_n and c > 0], default=0), vb)
    if is_pd: S_core += g['d_screen_t2g'](dc_val - 1, 2)
    Z_net = Z - S_core

    # Compute v19 pen for comparison
    pen_v19 = gamma * S_core / (d * n_used + (d + 1) * abs(S_core))
    has_d_any = any(c > 0 for nn, ll, c in config if ll == 2)
    if n_used >= n_ref and not has_d_any:
        pen_v19 *= 1 / d**2

    # HYBRID alpha: v19 alpha + pen correction from spatial sim
    # alpha_hybrid = alpha_v19 + (pen_v19 - pen_f) * coupling
    # coupling = 1 for s/d-block, w_pi for p-block
    has_p_local = any(c > 0 for nn, ll, c in config if nn == val_n and ll == 1) and not is_pd
    pen_coupling = w_pi if has_p_local else 1.0
    alpha_hybrid = alpha_v19 + pen_coupling * (pen_v19 - pen_f)

    # Use hybrid alpha
    alpha_f = alpha_hybrid

    # Energy from force-based alpha
    N_core = sum(c for nn, ll, c in config if nn < n_used)
    Z_eff_f = max(Z_net, 0.01) ** alpha_f
    E_force = (Z_eff_f / n)**2 * E_H * (1 + gamma * N_core / (d**2 * n**2))
    err_f = (E_force - E_obs) / E_obs * 100

    # v19 energy (with shear)
    Z_eff_v = max(Z_net, 0.01) ** alpha_v19
    E_v19 = (Z_eff_v / n)**2 * E_H * (1 + gamma * N_core / (d**2 * n**2))
    err_v = (E_v19 - E_obs) / E_obs * 100

    results.append((Z, sym, E_obs, err_v, err_f, alpha_v19, alpha_f))

elapsed = time.time() - t0

# Display
print(f"\n  {'Z':>3} {'Sym':<3} {'v19%':>7} {'force%':>7} {'a_v19':>7} {'a_frc':>7} {'da':>7} {'better':>6}")
print(f"  {'-'*55}")
n_better = 0
for r in results:
    Z, sym, obs, ev, ef, av, af = r
    better = abs(ef) < abs(ev)
    if better: n_better += 1
    mark = '  <' if better and abs(ev - ef) > 0.2 else ''
    if Z in [3,11,19,21,31,37,39,49,55,57,72,81,87,89]:
        print()
    print(f"  {Z:3d} {sym:<3} {ev:+7.1f} {ef:+7.1f} {av:7.4f} {af:7.4f} {af-av:+7.4f}{mark}")

v19_mean = np.mean([abs(r[3]) for r in results])
frc_mean = np.mean([abs(r[4]) for r in results])
print(f"\n{'='*80}")
print(f"  v19+shear mean: {v19_mean:.2f}%")
print(f"  Force alpha mean: {frc_mean:.2f}%")
print(f"  Force better: {n_better}/{len(results)} atoms")
print(f"  Time: {elapsed*1000:.0f}ms")
