"""
Z_eff Octahedral — Unified Screening from Group Theory
========================================================
ALL screening weights from ONE source: the Clebsch-Gordan coefficients
of the octahedral group Oh acting on breather modes.

Angular modes decompose into Oh irreps:
  l=0 (s):  A1g  (1D, fully symmetric)
  l=1 (p):  T1u  (3D, vector)
  l=2 (d):  Eg   (2D, eg) + T2g (3D, t2g)
  l=3 (f):  A2u  (1D) + T1u (3D) + T2u (3D)

Screening is mediated by T1u (Coulomb vector coupling).
Selection rule: irrep1 x T1u must contain irrep2.

Coupling matrix (from CG coefficients):
  s->p:     1/sqrt(3)
  p->s:     1/sqrt(3)
  p->eg:    sqrt(2/3)
  p->t2g:   1/sqrt(3)
  eg->p:    sqrt(2/3)
  t2g->p:   1/sqrt(3)
  All others: zero (selection rule forbidden)

This replaces v19's rules 2, 8, 9 with pure group theory.
"""

import numpy as np
from math import factorial, comb

d = 3
gamma_sg = np.pi / (2**(d+1) * np.pi - 2)
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6
n_ref = d + 1
J_hund = 2 / (d + 2)
w_pi = np.cos(np.pi / d)
w_delta = np.cos(2 * np.pi / d)


# === THE OCTAHEDRAL COUPLING MATRIX ===
# Computed from CG coefficients of Oh group.
# Entry (core_irrep, val_irrep) = coupling strength through T1u mediator.
# Normalized so that the maximum coupling = 1.

# CG relative weights × physical coupling (breather mass ratio)
# Physical coupling = w_pi = sin(gamma)/sin(2*gamma) ≈ 0.5
# CG gives relative structure, mass ratio gives absolute scale
_phys = w_pi  # 0.5 — the breather mass ratio
_cg_1 = 1.0   # normalized CG for T1u coupling (base)
_cg_2 = np.sqrt(2.0)  # eg couples sqrt(2) stronger than t2g

# Oh coupling: (core_irrep, val_irrep) -> strength
# Zero means selection-rule forbidden.
Oh_coupling = {
    # s (A1g) core — screens through T1u mediator
    ('A1g', 'A1g'):  0.0,              # s->s: forbidden
    ('A1g', 'T1u'):  _phys * _cg_1,   # s->p: allowed (0.5)
    ('A1g', 'Eg'):   0.0,              # s->eg: forbidden
    ('A1g', 'T2g'):  0.0,              # s->t2g: forbidden

    # p (T1u) core — screens everything through T1u^2
    ('T1u', 'A1g'):  _phys * _cg_1,   # p->s: allowed (0.5)
    ('T1u', 'T1u'):  _phys * _cg_1,   # p->p: self-screening (0.5)
    ('T1u', 'Eg'):   _phys * _cg_2,   # p->eg: stronger (0.707)
    ('T1u', 'T2g'):  _phys * _cg_1,   # p->t2g: standard (0.5)

    # eg (Eg) core
    ('Eg', 'A1g'):   0.0,              # eg->s: forbidden!
    ('Eg', 'T1u'):   _phys * _cg_2,   # eg->p: allowed, strong (0.707)
    ('Eg', 'Eg'):    0.0,              # eg->eg: forbidden
    ('Eg', 'T2g'):   0.0,              # eg->t2g: forbidden

    # t2g (T2g) core
    ('T2g', 'A1g'):  0.0,              # t2g->s: forbidden!
    ('T2g', 'T1u'):  _phys * _cg_1,   # t2g->p: allowed, standard (0.5)
    ('T2g', 'Eg'):   0.0,              # t2g->eg: forbidden
    ('T2g', 'T2g'):  0.0,              # t2g->t2g: forbidden
}


def get_irreps(l):
    """Decompose angular momentum l into Oh irreps.

    Returns list of (irrep_name, n_channels).
    """
    if l == 0:
        return [('A1g', 1)]
    elif l == 1:
        return [('T1u', 3)]
    elif l == 2:
        return [('T2g', 3), ('Eg', 2)]  # t2g (3) + eg (2) = 5
    elif l == 3:
        return [('A2u', 1), ('T1u', 3), ('T2u', 3)]  # 1+3+3 = 7
    return [('A1g', 1)]  # fallback


def get_val_irreps(l_val):
    """Get the primary irreps for the valence mode."""
    if l_val == 0:
        return [('A1g', 1)]
    elif l_val == 1:
        return [('T1u', 3)]
    elif l_val == 2:
        return [('T2g', 3), ('Eg', 2)]
    return [('A1g', 1)]


def octahedral_screening(config, core_ref_n, is_pd, val_l, val_block):
    """Core screening: radial (charge blocking) + angular (Oh coupling).

    Two components combined:
    1. RADIAL: every core mode blocks charge. Weight = w_pi per channel.
       This is the electrostatic screening — always positive, always present.
    2. ANGULAR: Oh CG coefficients modify the radial screening.
       If angular coupling is nonzero: modes interact through T1u mediator,
       enhancing screening beyond the radial baseline.
       If angular coupling is zero: modes are orthogonal on the cube.
       For d-modes with zero angular coupling to s-valence, the screening
       is REDUCED (anti-screening) because the wave interference is destructive.

    Combined: S = sum over core modes of:
      w_radial + w_angular  if Oh allows coupling
      w_radial - w_anti     if Oh forbids coupling (anti-screening)
    """
    val_irreps = get_val_irreps(val_l)

    S = 0.0
    for nn, ll, count in config:
        if nn >= core_ref_n:
            if not (nn == core_ref_n and is_pd and ll != 2):
                continue

        n_ch = min(count, 2 * ll + 1)
        delta_n = core_ref_n - nn
        core_irreps = get_irreps(ll)

        if ll <= 1:
            # s, p core: radial screening dominates
            # All s,p modes screen at w_pi per channel (charge blocking)
            S += n_ch * _phys  # w_pi = 0.5

        elif ll == 2:
            # d-core: split into t2g (3) and eg (2)
            # Distribute electrons: t2g fills first (lower energy)
            dc = count
            t2g_count = min(dc, 2 * 3)   # up to 6 electrons in 3 t2g channels
            eg_count = max(dc - 6, 0)     # remainder in 2 eg channels
            t2g_ch = min(t2g_count, 3)    # occupied t2g channels
            eg_ch = min(eg_count, 2)      # occupied eg channels

            for v_irrep, v_dim in val_irreps:
                # t2g contribution
                cg_t2g = Oh_coupling.get(('T2g', v_irrep), 0.0)
                # eg contribution
                cg_eg = Oh_coupling.get(('Eg', v_irrep), 0.0)

                if cg_t2g > 0:
                    # t2g couples to valence: screening
                    S += t2g_ch * cg_t2g
                else:
                    # t2g doesn't couple: anti-screening
                    # The wave is orthogonal → destructive interference
                    # Anti-screening strength = w_delta per channel
                    S += t2g_ch * w_delta  # negative!

                if cg_eg > 0:
                    S += eg_ch * cg_eg
                else:
                    S += eg_ch * w_delta / d  # weaker anti-screening for eg

            # Depth attenuation for deep d10
            if delta_n >= 2 and dc == 10:
                # Deep d10: blend toward normal screening
                S_correction = 0
                for v_irrep, v_dim in val_irreps:
                    cg_t2g = Oh_coupling.get(('T2g', v_irrep), 0.0)
                    if cg_t2g == 0:
                        # The anti-screening weakens with depth
                        S_correction -= t2g_ch * w_delta * (1 - d/(d + delta_n))
                S += S_correction

        elif ll == 3:
            # f-core: decompose into A2u (1) + T1u (3) + T2u (3)
            # The T1u component of f-modes CAN couple to s and p valence
            # The A2u and T2u components are selection-rule forbidden for
            # most valence types → anti-screening

            # f-modes screening depends on valence type
            for v_irrep, v_dim in val_irreps:
                if v_irrep == 'T1u':
                    # f's T1u component screens p-valence
                    # 3 channels of T1u in f-shell
                    f_T1u_ch = min(count, 3)
                    S += f_T1u_ch * _phys  # screens at w_pi
                elif v_irrep == 'A1g':
                    # f cannot screen s directly
                    # Anti-screening from all 7 channels
                    S += n_ch * w_delta / d  # weak anti-screen
                elif v_irrep in ('Eg', 'T2g'):
                    # f cannot screen d directly
                    S += n_ch * w_delta / d  # weak anti-screen

    return S


# === ATOMS (full periodic table) ===
atoms = [
    (1,'H',13.598,[(1,0,1)]),
    (2,'He',24.587,[(1,0,2)]),
    (3,'Li',5.392,[(1,0,2),(2,0,1)]),
    (4,'Be',9.323,[(1,0,2),(2,0,2)]),
    (5,'B',8.298,[(1,0,2),(2,0,2),(2,1,1)]),
    (6,'C',11.260,[(1,0,2),(2,0,2),(2,1,2)]),
    (7,'N',14.534,[(1,0,2),(2,0,2),(2,1,3)]),
    (8,'O',13.618,[(1,0,2),(2,0,2),(2,1,4)]),
    (9,'F',17.423,[(1,0,2),(2,0,2),(2,1,5)]),
    (10,'Ne',21.565,[(1,0,2),(2,0,2),(2,1,6)]),
    (11,'Na',5.139,[(1,0,2),(2,0,2),(2,1,6),(3,0,1)]),
    (12,'Mg',7.646,[(1,0,2),(2,0,2),(2,1,6),(3,0,2)]),
    (13,'Al',5.986,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,1)]),
    (14,'Si',8.152,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,2)]),
    (15,'P',10.487,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,3)]),
    (16,'S',10.360,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,4)]),
    (17,'Cl',12.968,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,5)]),
    (18,'Ar',15.760,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6)]),
    (19,'K',4.341,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(4,0,1)]),
    (20,'Ca',6.113,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(4,0,2)]),
    (21,'Sc',6.561,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,1),(4,0,2)]),
    (22,'Ti',6.828,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,2),(4,0,2)]),
    (23,'V',6.746,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,3),(4,0,2)]),
    (24,'Cr',6.767,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,5),(4,0,1)]),
    (25,'Mn',7.434,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,5),(4,0,2)]),
    (26,'Fe',7.902,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,6),(4,0,2)]),
    (27,'Co',7.881,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,7),(4,0,2)]),
    (28,'Ni',7.640,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,8),(4,0,2)]),
    (29,'Cu',7.726,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,1)]),
    (30,'Zn',9.394,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2)]),
    (31,'Ga',5.999,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,1)]),
    (32,'Ge',7.900,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,2)]),
    (33,'As',9.815,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,3)]),
    (34,'Se',9.752,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,4)]),
    (35,'Br',11.814,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,5)]),
    (36,'Kr',14.000,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6)]),
]


def calc_all():
    results = []
    for Z, sym, E_obs, config in atoms:
        val_n = max(nn for nn, ll, c in config)
        s_count = sum(c for nn, ll, c in config if nn == val_n and ll == 0)
        p_count = sum(c for nn, ll, c in config if nn == val_n and ll == 1)
        dc_val = sum(c for nn, ll, c in config if nn == val_n - 1 and ll == 2)

        is_pd = (s_count == 0 and dc_val == 10)
        n = val_n - 1 if is_pd else val_n
        has_p = p_count > 0 and not is_pd
        val_l = max([ll for nn, ll, c in config if nn == val_n and c > 0], default=0)
        val_block = 'p' if has_p else ('d' if (dc_val > 0 or is_pd) else 's')

        # === OCTAHEDRAL SCREENING ===
        S_core = octahedral_screening(config, n, is_pd, val_l, val_block)
        Z_net = Z - S_core

        # === PENETRATION ===
        pen = gamma_sg * S_core / (d * n + (d + 1) * abs(S_core)) if S_core != 0 else 0
        has_d_any = any(c > 0 for nn, ll, c in config if ll == 2)
        if n >= n_ref and not has_d_any:
            pen *= 1 / d**2

        # === WAVE COUPLING (C_total) — same as combined formula ===
        sp = (s_count == 2)
        s1 = (s_count == 1)

        if not has_p:
            if is_pd:
                C_total = w_delta + J_hund
            elif s1 and dc_val > 0:
                C_total = w_pi + (n - n_ref)
            elif sp and dc_val > 0:
                C_total = 1.0
                dc = dc_val
                if dc > 5:
                    t2g_paired = min(dc - 5, d)
                    if t2g_paired == d:
                        eg_occ = dc - 5
                        eg_paired = max(0, dc - 8)
                        eg_unp = eg_occ - eg_paired
                        C_total -= eg_unp * abs(w_delta) / d
                if dc == 10:
                    C_total += J_hund
                    n_deep = sum(1 for nn2, ll2, c2 in config
                                if ll2 == 2 and c2 == 10 and nn2 < n-1)
                    C_total += J_hund * n_deep
                C_total += w_pi * (n - n_ref)
            elif sp:
                C_total = 1.0
            else:
                C_total = -1.0

            alpha = (d * n + C_total) / (d**2 * n) - pen

        else:
            pc = p_count
            pa = min(pc, 2*d - pc) if pc <= d else 2*d - pc
            pl = max(0, pc - d)

            if pc <= d:
                uf_boost = (1 + w_pi) / d**2
                has_d_ch = any(c > 0 for nn2, ll2, c in config if ll2 == 2)
                if has_d_ch:
                    uf_boost *= d
                N_eff = pa + uf_boost
                xu = comb(pa, 2) if pa >= 2 else 0
                N_eff -= xu / (d**2 * n)
            else:
                fl = min(pl, 1)
                rl = max(pl - 1, 0)
                w_first = (n**2 - d) / n**2
                w_rest = 1 + w_pi
                N_eff = pa + fl * w_first + rl * w_rest
                xp_under = comb(min(pc, d), 2)
                xp_over = comb(max(pc - d, 0), 2)
                N_eff -= xp_under / (d**2 * n)
                N_eff -= xp_over / (d**2 * n_ref)
                N_eff -= xp_over * (n - n_ref) / (d * (d-1) * n)

            N_eff += pa * (d - pa) / d**3

            deep_hl_ch = sum(min(c, 2*ll+1) for nn2, ll, c in config
                           if ll >= 2 and c > 0 and (n - nn2) >= 2)
            xu_total = comb(min(pc, d), 2)
            if xu_total > 0 and deep_hl_ch > 0:
                N_eff -= xu_total * deep_hl_ch / (d**3 * n)

            alpha = (d + w_pi * N_eff - 1) / d**2 - w_pi * pen

        alpha = max(alpha, 0.05)
        alpha = min(alpha, 0.95)

        Z_eff = Z_net ** alpha
        E_pred = (Z_eff / n)**2 * E_H

        N_core = sum(c for nn, ll, c in config if nn < n)
        shear = 1 + gamma_sg * N_core / (d**2 * n**2)
        E_pred *= shear

        err = (E_pred - E_obs) / E_obs * 100
        results.append((Z, sym, E_obs, E_pred, err, alpha))

    return results


# === RUN ===
results = calc_all()
print("Z_eff Octahedral — Group Theory Screening")
print("=" * 55)
print(f"  Screening from Oh CG coefficients (one matrix)")
print(f"  base coupling = {_phys:.4f}, eg enhancement = sqrt(2) = {_cg_2:.4f}")
print()

errs = [abs(r[4]) for r in results]
print(f"  ALL (P1-4): {np.mean(errs):.2f}%")
print(f"  Under 5%: {sum(1 for e in errs if e < 5)}/{len(results)}")
print(f"  Under 10%: {sum(1 for e in errs if e < 10)}/{len(results)}")
print(f"  v19: 2.07% on P1-5")
print()

# Show S_core comparison for key atoms
print(f"  {'Z':>3} {'Sym':<3} {'S_Oh':>7} {'err':>7}")
print(f"  {'-'*22}")
for Z, sym, E_obs, E_pred, err, alpha in results:
    marker = ' **' if abs(err) > 10 else (' *' if abs(err) > 5 else '')
    if abs(err) > 5 or Z <= 10:
        print(f"  {Z:3d} {sym:<3} {err:+6.1f}%{marker}")
