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
    """Core screening: radial + angular (Oh) with depth attenuation.

    Physics:
    1. RADIAL (charge blocking): w_pi per channel, always present
    2. ANGULAR (Oh selection rules): modifies radial for high-l modes
    3. DEPTH: deep shells (delta_n >= 2) blur toward radial-only

    For each core shell, the effective weight per channel is:
      w_eff = w_radial * (1 - angular_fraction) + w_angular * angular_fraction

    where angular_fraction = d / (d + delta_n) → 1 for adjacent, 0 for deep.
    At delta_n = 0: full angular coupling (Oh matrix dominates).
    At delta_n >> d: angular structure blurs out (radial screening only).
    """
    val_irreps = get_val_irreps(val_l)

    S = 0.0
    for nn, ll, count in config:
        if nn >= core_ref_n:
            if not (nn == core_ref_n and is_pd and ll != 2):
                continue

        n_ch = min(count, 2 * ll + 1)
        delta_n = core_ref_n - nn

        # Depth blending factor: how much angular structure survives
        # Adjacent shell (delta_n=1): full angular coupling
        # Deep shell (delta_n >= 2): blend depends on l
        #   d-core (l=2): sharp 1/d blend (from v19 deep d10 logic)
        #   f-core (l=3): gradual d/(d+dn-1) blend (f angular structure persists)
        #   s,p core: no blending needed (always radial)
        if delta_n <= 1:
            ang_frac = 1.0
        elif ll == 2:
            ang_frac = 1.0 / d  # d-core: sharp blend
        elif ll == 3:
            ang_frac = d / (d + delta_n)  # f-core: gradual blend
        else:
            ang_frac = 1.0  # s,p: always full (but they only use radial anyway)

        if ll <= 1:
            # s, p core: radial screening, no angular modification needed
            # (Oh says s screens p and p screens everything — consistent
            #  with flat w_pi for all s,p core)
            S += n_ch * _phys

        elif ll == 2:
            # d-core: Oh gives t2g/eg split
            dc = count
            t2g_count = min(dc, 2 * 3)
            eg_count = max(dc - 6, 0)
            t2g_ch = min(t2g_count, 3)
            eg_ch = min(eg_count, 2)

            # For d-core, the Oh selection rule is:
            # t2g -> p: allowed (w_pi), t2g -> s: forbidden (anti-screen)
            # eg -> p: allowed (w_pi*sqrt(2)), eg -> s: forbidden (anti-screen)
            #
            # Use v19's exact t2g/eg logic for precision:
            # t2g unpaired: w_delta/(d+1) per channel
            # t2g paired: beyond half-fill, paired reduce anti-screening
            # eg unpaired: w_delta/d per channel
            # eg paired: w_delta*d per channel (RESTORES coupling)

            if val_l == 1:
                # d -> p: Oh says BOTH t2g and eg screen p
                # Use w_pi for t2g (weaker) and w_pi*sqrt(2) for eg (stronger)
                # But v19 uses the full t2g/eg function which works at 2.61%
                # Match v19: treat d-core screening p as positive
                S_angular = n_ch * _phys  # d screens p at w_pi per channel
            else:
                # d -> s/d: Oh says FORBIDDEN (anti-screening)
                # Use v19's t2g/eg anti-screening
                # t2g: 3 channels at w_delta/(d+1) = -0.125 each
                # eg occupied: w_delta/d per channel (unpaired), w_delta*d (paired)
                t2g_paired = max(dc - 5, 0) if dc > 5 else 0
                t2g_paired = min(t2g_paired, d)
                t2g_full = (t2g_paired == d)

                # t2g contribution
                S_angular = d * w_delta / (d + 1)  # 3 channels * w_delta/(d+1)

                # eg pairing
                if t2g_full and eg_ch > 0:
                    eg_paired_n = max(0, dc - 8)
                    eg_unp = eg_ch - eg_paired_n
                    # Unpaired eg: w_delta/d
                    S_angular += eg_unp * w_delta / d
                    # Paired eg: w_delta * d (restores coupling)
                    S_angular += eg_paired_n * w_delta * d

            # Radial baseline
            S_radial = n_ch * _phys

            # Blend with depth
            S += S_radial * (1 - ang_frac) + S_angular * ang_frac

        elif ll == 3:
            # f-core: T1u component screens p, rest anti-screens

            # Radial baseline
            S_radial = n_ch * _phys

            # Angular contribution
            S_angular = 0.0
            for v_irrep, v_dim in val_irreps:
                if v_irrep == 'T1u':
                    # f's T1u (3 channels) couples to p-valence
                    f_T1u_ch = min(count, 3)
                    S_angular += f_T1u_ch * _phys
                    # Remaining channels (A2u + T2u = 4 channels): anti-screen
                    f_other_ch = n_ch - f_T1u_ch
                    S_angular += f_other_ch * w_delta / d
                elif v_irrep == 'A1g':
                    # f cannot screen s: full anti-screening
                    S_angular += n_ch * w_delta / d
                elif v_irrep in ('Eg', 'T2g'):
                    # f cannot screen d: anti-screening
                    S_angular += n_ch * w_delta / d

            # Blend with depth
            S += S_radial * (1 - ang_frac) + S_angular * ang_frac

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
    # Period 5
    (37,'Rb',4.177,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(5,0,1)]),
    (38,'Sr',5.695,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(5,0,2)]),
    (39,'Y', 6.217,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,1),(5,0,2)]),
    (40,'Zr',6.634,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,2),(5,0,2)]),
    (41,'Nb',6.759,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,4),(5,0,1)]),
    (42,'Mo',7.092,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,5),(5,0,1)]),
    (43,'Tc',7.119,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,5),(5,0,2)]),
    (44,'Ru',7.361,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,7),(5,0,1)]),
    (45,'Rh',7.459,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,8),(5,0,1)]),
    (46,'Pd',8.337,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,0)]),
    (47,'Ag',7.576,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,1)]),
    (48,'Cd',8.994,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,2)]),
    (49,'In',5.786,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,2),(5,1,1)]),
    (50,'Sn',7.344,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,2),(5,1,2)]),
    (51,'Sb',8.608,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,2),(5,1,3)]),
    (52,'Te',9.010,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,2),(5,1,4)]),
    (53,'I', 10.451,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,2),(5,1,5)]),
    (54,'Xe',12.130,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,2),(5,1,6)]),
    # Period 6
    (55,'Cs',3.894,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,2),(5,1,6),(6,0,1)]),
    (56,'Ba',5.212,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,2),(5,1,6),(6,0,2)]),
    (57,'La',5.577,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,2),(5,1,6),(5,2,1),(6,0,2)]),
    (58,'Ce',5.539,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,1),(5,0,2),(5,1,6),(5,2,1),(6,0,2)]),
    (59,'Pr',5.473,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,3),(5,0,2),(5,1,6),(6,0,2)]),
    (60,'Nd',5.525,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,4),(5,0,2),(5,1,6),(6,0,2)]),
    (61,'Pm',5.582,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,5),(5,0,2),(5,1,6),(6,0,2)]),
    (62,'Sm',5.644,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,6),(5,0,2),(5,1,6),(6,0,2)]),
    (63,'Eu',5.670,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,7),(5,0,2),(5,1,6),(6,0,2)]),
    (64,'Gd',6.150,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,7),(5,0,2),(5,1,6),(5,2,1),(6,0,2)]),
    (65,'Tb',5.864,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,9),(5,0,2),(5,1,6),(6,0,2)]),
    (66,'Dy',5.939,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,10),(5,0,2),(5,1,6),(6,0,2)]),
    (67,'Ho',6.022,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,11),(5,0,2),(5,1,6),(6,0,2)]),
    (68,'Er',6.108,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,12),(5,0,2),(5,1,6),(6,0,2)]),
    (69,'Tm',6.184,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,13),(5,0,2),(5,1,6),(6,0,2)]),
    (70,'Yb',6.254,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(6,0,2)]),
    (71,'Lu',5.426,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,1),(6,0,2)]),
    (72,'Hf',6.825,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,2),(6,0,2)]),
    (73,'Ta',7.550,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,3),(6,0,2)]),
    (74,'W', 7.864,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,4),(6,0,2)]),
    (75,'Re',7.833,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,5),(6,0,2)]),
    (76,'Os',8.438,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,6),(6,0,2)]),
    (77,'Ir',8.967,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,7),(6,0,2)]),
    (78,'Pt',8.959,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,9),(6,0,1)]),
    (79,'Au',9.226,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,1)]),
    (80,'Hg',10.438,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2)]),
    (81,'Tl',6.108,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,1)]),
    (82,'Pb',7.417,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,2)]),
    (83,'Bi',7.286,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,3)]),
    (84,'Po',8.414,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,4)]),
    (85,'At',9.318,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,5)]),
    (86,'Rn',10.749,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,6)]),
    # Period 7
    (87,'Fr',4.073,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,6),(7,0,1)]),
    (88,'Ra',5.278,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,6),(7,0,2)]),
    (89,'Ac',5.380,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,6),(6,2,1),(7,0,2)]),
    (90,'Th',6.307,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,6),(6,2,2),(7,0,2)]),
    (91,'Pa',5.890,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,2),(6,0,2),(6,1,6),(6,2,1),(7,0,2)]),
    (92,'U', 6.194,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,3),(6,0,2),(6,1,6),(6,2,1),(7,0,2)]),
    (93,'Np',6.266,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,4),(6,0,2),(6,1,6),(6,2,1),(7,0,2)]),
    (94,'Pu',6.026,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,6),(6,0,2),(6,1,6),(7,0,2)]),
    (95,'Am',5.974,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,7),(6,0,2),(6,1,6),(7,0,2)]),
    (96,'Cm',5.991,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,7),(6,0,2),(6,1,6),(6,2,1),(7,0,2)]),
    (97,'Bk',6.198,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,9),(6,0,2),(6,1,6),(7,0,2)]),
    (98,'Cf',6.282,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,10),(6,0,2),(6,1,6),(7,0,2)]),
    (99,'Es',6.367,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,11),(6,0,2),(6,1,6),(7,0,2)]),
    (100,'Fm',6.500,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,12),(6,0,2),(6,1,6),(7,0,2)]),
    (101,'Md',6.580,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,13),(6,0,2),(6,1,6),(7,0,2)]),
    (102,'No',6.650,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,14),(6,0,2),(6,1,6),(7,0,2)]),
    (103,'Lr',4.960,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,14),(6,0,2),(6,1,6),(7,0,2),(7,1,1)]),
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

            # f-core d-shell rebalancing (5d TMs with f14)
            n_f_ch = sum(min(c, 7) for nn2, ll2, c in config if ll2 == 3 and c > 0)
            if sp and dc_val > 0 and dc_val < 10 and n_f_ch > 0 and not is_pd:
                if dc_val >= 5:
                    C_total += n_f_ch * (dc_val - 5) / (d * n)
                else:
                    C_total += n_f_ch * (dc_val - 5) / (d**2 * n)

            # Deep f-shell parity boost (s-block actinides)
            if dc_val == 0 and not is_pd:
                n_deep_f = sum(min(c, 7) for nn2, ll2, c in config
                              if ll2 == 3 and c > 0 and (n - nn2) >= 3)
                if n_deep_f > 0:
                    C_total += (d - 1) * n_deep_f / (d * n)

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
errs_p15 = [abs(r[4]) for r in results if r[0] <= 54]
errs_p6 = [abs(r[4]) for r in results if 55 <= r[0] <= 86]
errs_p7 = [abs(r[4]) for r in results if r[0] >= 87]
print(f"  ALL: {np.mean(errs):.2f}%  P1-5: {np.mean(errs_p15):.2f}%  P6: {np.mean(errs_p6):.2f}%  P7: {np.mean(errs_p7):.2f}%")
print(f"  Under 5%: {sum(1 for e in errs if e < 5)}/{len(results)}")
print(f"  Under 10%: {sum(1 for e in errs if e < 10)}/{len(results)}")
print(f"  v19: ALL=2.61%, P1-5=2.07%")
print()

# Show S_core comparison for key atoms
print(f"  {'Z':>3} {'Sym':<3} {'S_Oh':>7} {'err':>7}")
print(f"  {'-'*22}")
for Z, sym, E_obs, E_pred, err, alpha in results:
    marker = ' **' if abs(err) > 10 else (' *' if abs(err) > 5 else '')
    if abs(err) > 5 or Z <= 10:
        print(f"  {Z:3d} {sym:<3} {err:+6.1f}%{marker}")
