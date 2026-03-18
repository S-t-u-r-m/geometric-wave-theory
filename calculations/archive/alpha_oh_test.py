#!/usr/bin/env python3
"""
Test: Can we match V19 using Oh A1g content for alpha corrections?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial, comb
import importlib.util, contextlib

d = 3
gamma = np.pi / (2**(d+1) * np.pi - 2)
alpha_em_c = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em_c**2 / 2) * 0.51100e6
w_pi = np.cos(np.pi/d)
w_delta = np.cos(2*np.pi/d)
J_hund = 2 / (d + 2)
n_ref = d + 1

# Oh closed forms
def a1g_T1u(p):
    if p <= 0 or p % 2 == 1: return 0
    return (3**p + 15) // 24

def a1g_T2g(n):
    if n <= 0: return 0
    return (3**n + 6 + 9*(-1)**n) // 24

def a1g_Eg(n):
    if n <= 0: return 0
    return (2**n + 2*(-1)**n) // 6

# Load V19
spec = importlib.util.spec_from_file_location("zeff", "calculations/z_eff_final.py")
mod = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(mod)


def alpha_oh(Z, config):
    """
    Compute alpha using V19 structure but with Oh A1g quantities
    replacing ad-hoc corrections where possible.

    Strategy: keep V19's framework (it works), but tag each correction
    with its Oh origin. This is DOCUMENTATION, not revolution.
    The formula IS the Oh tensor product — V19 just didn't know it.
    """
    val_n = max(nn for nn, ll, c in config)
    s_count = sum(c for nn, ll, c in config if nn == val_n and ll == 0)
    p_count = sum(c for nn, ll, c in config if nn == val_n and ll == 1)
    dc_val = sum(c for nn, ll, c in config if nn == val_n - 1 and ll == 2)
    is_pd = (s_count == 0 and dc_val == 10)
    n = val_n - 1 if is_pd else val_n
    val_l = max([ll for nn, ll, c in config if nn == val_n and c > 0], default=0)
    has_p = p_count > 0 and not is_pd

    # Screening (V19 — already Oh CG matrix)
    S_core = mod.screening(config, n, is_pd, val_l)

    # Penetration (V19)
    pen = gamma * S_core / (d*n + (d+1)*abs(S_core)) if S_core != 0 else 0
    has_d_any = any(c > 0 for nn, ll, c in config if ll == 2)
    if n >= n_ref and not has_d_any:
        pen *= 1 / d**2

    # Use EXACTLY V19's compute_alpha — it already IS the Oh tensor product
    # We just recognize each term's origin
    alpha = mod.compute_alpha(config, n, is_pd, val_l, has_p,
                              s_count, dc_val, p_count, S_core)

    # Now ADD the A1g-based N-body correction on top
    # This is NEW — V19 doesn't have explicit N-body correction

    # P-block N-body: A1g(T1u^p) / (dim * d)
    if has_p:
        a1g_p = a1g_T1u(p_count)
        if a1g_p > 0:
            dim_p = 3**p_count
            # The A1g correction modifies the effective coupling
            # Direction: higher A1g → more internal coupling → LOWER alpha
            # (energy tied up in internal correlations)
            nbody_p = a1g_p / (dim_p * d**2)
            alpha -= w_pi * nbody_p

    # D-shell N-body: A1g(T2g^t2g) + A1g(Eg^eg)
    if dc_val > 0 and not is_pd:
        if dc_val <= 3: t2g_e, eg_e = dc_val, 0
        elif dc_val <= 5: t2g_e, eg_e = 3, dc_val - 3
        elif dc_val <= 8: t2g_e, eg_e = dc_val - 2, 2
        else: t2g_e, eg_e = 6, dc_val - 6

        a_t2g = a1g_T2g(t2g_e)
        a_eg = a1g_Eg(eg_e)

        if a_t2g > 0:
            dim_t2g = 3**t2g_e
            alpha -= a_t2g / (dim_t2g * d**3 * n)

        if a_eg > 0:
            dim_eg = 2**eg_e
            alpha -= a_eg / (dim_eg * d**3 * n)

    return max(alpha, 0.05), S_core


# ============================================================
# TEST
# ============================================================
print("V19 + Oh N-body corrections: All 103 atoms")
print("=" * 70)

errs_v19 = []
errs_oh = []
changes = []

for Z, sym, IE_obs, config in mod.atoms:
    # V19
    IE_v19, a_v19, S_v19 = mod.ionization_energy(Z, config)
    ev = (IE_v19 - IE_obs) / IE_obs * 100
    errs_v19.append(abs(ev))

    # Oh-enhanced
    a_oh, S_oh = alpha_oh(Z, config)
    val_n = max(nn for nn, ll, c in config)
    s_c = sum(c for nn, ll, c in config if nn == val_n and ll == 0)
    dc = sum(c for nn, ll, c in config if nn == val_n - 1 and ll == 2)
    is_pd = (s_c == 0 and dc == 10)
    n = val_n - 1 if is_pd else val_n

    Z_net = Z - S_oh
    Z_eff = Z_net ** a_oh
    IE_oh = (Z_eff / n)**2 * E_H
    N_core = sum(c for nn, ll, c in config if nn < n)
    IE_oh *= 1 + gamma * N_core / (d**2 * n**2)
    eo = (IE_oh - IE_obs) / IE_obs * 100
    errs_oh.append(abs(eo))

    if abs(abs(eo) - abs(ev)) > 0.5:
        changes.append((Z, sym, ev, eo, a_v19, a_oh))

print(f"V19:        mean={np.mean(errs_v19):.2f}%  <3%={sum(1 for e in errs_v19 if e<3)}/103  <5%={sum(1 for e in errs_v19 if e<5)}/103  <10%={sum(1 for e in errs_v19 if e<10)}/103")
print(f"V19+Oh_Nb:  mean={np.mean(errs_oh):.2f}%  <3%={sum(1 for e in errs_oh if e<3)}/103  <5%={sum(1 for e in errs_oh if e<5)}/103  <10%={sum(1 for e in errs_oh if e<10)}/103")
print()

# Did N-body corrections help or hurt?
improved = sum(1 for i in range(len(errs_v19)) if errs_oh[i] < errs_v19[i] - 0.1)
worsened = sum(1 for i in range(len(errs_v19)) if errs_oh[i] > errs_v19[i] + 0.1)
print(f"Improved: {improved}  Worsened: {worsened}  Unchanged: {103-improved-worsened}")
print()

if changes:
    print(f"Significant changes (>{0.5}%):")
    print(f"{'Z':>3} {'Sym':>3} {'V19':>8} {'Oh+Nb':>8} {'change':>8} {'better?':>8}")
    for Z, sym, ev, eo, av, ao in sorted(changes, key=lambda x: abs(x[3]) - abs(x[2])):
        better = "YES" if abs(eo) < abs(ev) else "no"
        print(f"{Z:3d} {sym:>3} {ev:+8.2f}% {eo:+8.2f}% {abs(eo)-abs(ev):+8.2f}% {better:>8}")
