#!/usr/bin/env python3
"""
v18: Z_eff from d=3 — clean rewrite with all rules explicit.

Rules (7 total, all from d=3):
  1. Channel weights: w_l = cos(l·π/d)
     w_s=1, w_pi=0.5, w_delta=-0.5
  2. Core screening:
     - s,p core: w_pi per channel
     - d core: t2g/eg cubic split (w_delta vs w_delta/d)
     - f core: VALENCE-DEPENDENT
       s-valence: w_f = w_delta/d = -1/6  (anti-screen)
       p-valence: w_f = w_pi = +1/2       (screen — same angular parity)
       d-valence: w_f = w_delta/d = -1/6  (anti-screen)
     - Deep d10 (2+ shells below): partial extra anti-screening
  3. Penetration: pen = γ·S_core / (d·n + (d+1)·|S_core|)
     - Diluted by 1/d² when d-channel available but empty
  4. s/d-block parity:
     - s-pair: +1, s-single: -1
     - s-single + d: w_pi + (n - n_ref)
     - s-pair + d: eg correction + J_hund(d10) + w_pi·(n - n_ref)
     - Pd-type: w_delta + J_hund
     - Deep d10 exchange: J_hund per deep d10 (when valence d10)
  5. p-block N_eff:
     - Underfill: pa active + empty-channel boost from d-shells
     - Overfill: first lone pair at (n²-d)/n², rest at 1+w_pi
     - Exchange: underfill pairs/(d²n), overfill pairs/(d²·n_ref) + n-scaling
     - Mode-vacancy coupling: pa·(d-pa)/d³
  6. Alpha: (d·n + parity)/(d²·n) - pen  [s/d-block]
           (d + w_pi·N_eff - 1)/d² - w_pi·pen  [p-block]
  7. E_pred = (Z_net^α / n)² × E_H

New in v18: valence-dependent f-screening (#2 above).
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial, comb

# === CONSTANTS FROM d=3 ===
d = 3
gamma = np.pi / (2**(d+1)*np.pi - 2)
S_tunnel = 2**(2*d+1) / np.pi**2
BZ = np.log(2*d)
exponent = (2/factorial(d)) * (S_tunnel + BZ)
alpha_em = np.exp(-exponent)
E_H = (alpha_em**2 / 2) * 0.51100e6
w_pi = np.cos(np.pi / d)       # 0.5
w_delta = np.cos(2*np.pi / d)  # -0.5
J_hund = 2 / (d + 2)           # 0.4
n_ref = d + 1                   # 4

# f-screening weights by valence type
W_F = {
    's': w_delta / d,   # -1/6  (anti-screen: different angular parity)
    'p': w_pi,          # +1/2  (screen: same angular parity)
    'd': w_delta / d,   # -1/6  (anti-screen: different angular parity)
}

atoms = [
    # Period 1
    (1,'H',13.598,[(1,0,1)]), (2,'He',24.587,[(1,0,2)]),
    # Period 2
    (3,'Li',5.392,[(1,0,2),(2,0,1)]), (4,'Be',9.323,[(1,0,2),(2,0,2)]),
    (5,'B',8.298,[(1,0,2),(2,0,2),(2,1,1)]), (6,'C',11.260,[(1,0,2),(2,0,2),(2,1,2)]),
    (7,'N',14.534,[(1,0,2),(2,0,2),(2,1,3)]), (8,'O',13.618,[(1,0,2),(2,0,2),(2,1,4)]),
    (9,'F',17.423,[(1,0,2),(2,0,2),(2,1,5)]), (10,'Ne',21.565,[(1,0,2),(2,0,2),(2,1,6)]),
    # Period 3
    (11,'Na',5.139,[(1,0,2),(2,0,2),(2,1,6),(3,0,1)]),
    (12,'Mg',7.646,[(1,0,2),(2,0,2),(2,1,6),(3,0,2)]),
    (13,'Al',5.986,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,1)]),
    (14,'Si',8.152,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,2)]),
    (15,'P',10.487,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,3)]),
    (16,'S',10.360,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,4)]),
    (17,'Cl',12.968,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,5)]),
    (18,'Ar',15.760,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6)]),
    # Period 4
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
    # Period 6 — s-block
    (55,'Cs',3.894,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,2),(5,1,6),(6,0,1)]),
    (56,'Ba',5.212,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,2),(5,1,6),(6,0,2)]),
    # Lanthanides
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
    # 5d TMs
    (72,'Hf',6.825,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,2),(6,0,2)]),
    (73,'Ta',7.550,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,3),(6,0,2)]),
    (74,'W', 7.864,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,4),(6,0,2)]),
    (75,'Re',7.833,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,5),(6,0,2)]),
    (76,'Os',8.438,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,6),(6,0,2)]),
    (77,'Ir',8.967,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,7),(6,0,2)]),
    (78,'Pt',8.959,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,9),(6,0,1)]),
    (79,'Au',9.226,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,1)]),
    (80,'Hg',10.438,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2)]),
    # 6p block
    (81,'Tl',6.108,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,1)]),
    (82,'Pb',7.417,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,2)]),
    (83,'Bi',7.286,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,3)]),
    (84,'Po',8.414,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,4)]),
    (85,'At',9.318,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,5)]),
    (86,'Rn',10.749,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,6)]),
]


# === SCREENING FUNCTIONS ===

def d_screen_t2g(dc, val_l):
    """d-shell screening with t2g/eg cubic split."""
    n_t2g, n_eg = d, 2
    w_t2g, w_eg = w_delta, w_delta / d
    if val_l != 0:
        return min(dc, 5) * w_pi
    if dc <= n_t2g:
        return dc * w_t2g
    elif dc <= 5:
        return n_t2g * w_t2g + (dc - n_t2g) * w_eg
    elif dc <= 5 + n_t2g:
        np_ = dc - 5
        return (n_t2g - np_) * w_t2g + np_ * w_t2g / (d + 1) + n_eg * w_eg
    else:
        nep = dc - 5 - n_t2g
        return n_t2g * w_t2g / (d + 1) + (n_eg - nep) * w_eg + nep * w_eg * d


def core_screening(config, core_ref_n, is_pd, val_l, val_block):
    """Core screening sum. val_block is 's', 'p', or 'd'."""
    w_f = W_F[val_block]
    S = 0.0
    for nn, ll, count in config:
        if nn < core_ref_n:
            if ll <= 1:
                S += min(count, 2 * ll + 1) * w_pi
            elif ll == 2:
                delta_n = core_ref_n - nn
                if delta_n >= 2 and count == 10:
                    # Deep d10: partial extra anti-screening
                    S_anti = d_screen_t2g(count, val_l)
                    S_normal = 5 * w_pi
                    S += S_normal + (S_anti - S_normal) / d
                else:
                    S += d_screen_t2g(count, val_l)
            elif ll == 3:
                S += min(count, 7) * w_f
        elif nn == core_ref_n and is_pd and ll != 2:
            S += min(count, 2 * ll + 1) * w_pi
    return S


# === D-SHELL INFO ===

def get_d_shell_info(config, val_n):
    dc = sum(c for nn, ll, c in config if nn == val_n - 1 and ll == 2)
    if dc == 0:
        return 0, 0, False, False
    n_t2g = d
    if dc <= 5:
        t2g_paired = 0
    elif dc <= 8:
        t2g_paired = dc - 5
    else:
        t2g_paired = n_t2g
    eg_occ = max(min(dc, 5) - n_t2g, 0)
    if dc > 5:
        eg_occ = 2
    eg_paired = max(dc - 5 - n_t2g, 0) if dc > 5 + n_t2g else 0
    eg_unp = eg_occ - eg_paired
    t2g_full = (t2g_paired == n_t2g)
    is_d10 = (dc == 10)
    return dc, eg_unp, t2g_full, is_d10


def count_deep_d10(config, val_n):
    val_d_n = val_n - 1
    return sum(1 for nn, ll, c in config if ll == 2 and nn < val_d_n and c == 10)


# === MAIN CALCULATION ===

def calc_all():
    results = []
    for Z, sym, E_ion, config in atoms:
        val_n = max(nn for nn, ll, c in config)
        val_shells = sorted(
            [(nn, ll, c) for nn, ll, c in config if nn == val_n],
            key=lambda x: x[1], reverse=True
        )
        val_l = val_shells[0][1]
        n = val_n

        s_count = sum(c for nn, ll, c in config if nn == val_n and ll == 0)
        p_count = sum(c for nn, ll, c in config if nn == val_n and ll == 1)
        dc_val = sum(c for nn, ll, c in config if nn == val_n - 1 and ll == 2)

        # Pd special case
        is_pd = (s_count == 0 and dc_val == 10)
        if is_pd:
            n = val_n - 1

        # Determine valence block for f-screening
        has_p = val_l >= 1 and p_count > 0 and not is_pd
        if has_p:
            val_block = 'p'
        elif dc_val > 0 or is_pd:
            val_block = 'd'
        else:
            val_block = 's'

        # Core screening
        S_core = core_screening(config, n, is_pd, val_l, val_block)
        if is_pd:
            S_core += d_screen_t2g(dc_val - 1, 2)

        Z_net = Z - S_core

        # Penetration
        pen = gamma * S_core / (d * n + (d + 1) * abs(S_core))
        has_d_any = any(c > 0 for nn, ll, c in config if ll == 2)
        if n >= n_ref and not has_d_any:
            pen *= 1 / d**2

        sp = 1 if s_count == 2 else 0
        s1 = 1 if s_count == 1 else 0

        if not has_p:
            # === s/d-block: parity ===
            if is_pd:
                parity = w_delta + J_hund
            elif s1 and dc_val > 0:
                parity = w_pi + (n - n_ref)
            else:
                parity = +1 if sp else -1

            # s-pair TM corrections
            if sp and dc_val > 0 and not is_pd:
                dc, eg_unp, t2g_full, is_d10 = get_d_shell_info(config, val_n)
                if dc > 0:
                    if t2g_full:
                        parity -= eg_unp * abs(w_delta) / d
                    if is_d10:
                        parity += J_hund
                parity += w_pi * (n - n_ref)

            # Deep d10 exchange (when valence d10 + s-pair)
            dc2, _, _, is_d10_2 = get_d_shell_info(config, val_n)
            if sp and is_d10_2 and not is_pd:
                n_deep = count_deep_d10(config, val_n)
                parity += J_hund * n_deep

            alpha = (d * n + parity) / (d**2 * n)
            alpha -= pen
        else:
            # === p-block: N_eff ===
            pc = p_count
            if pc <= d:
                pa, pl = pc, 0
            else:
                pa, pl = 2 * d - pc, pc - d

            total_d_ch = sum(
                min(c, 2 * ll + 1)
                for nn, ll, c in config if nn < val_n and ll == 2
            )

            w_first = (n**2 - d) / n**2
            w_rest = 1 + w_pi
            fl = min(pl, 1)
            rl = max(pl - 1, 0)
            N_eff = pa + fl * w_first + rl * w_rest

            # Underfill empty-channel boost
            if pl == 0:
                n_empty = d - pa
                uf_scale = 1 + total_d_ch * abs(w_delta) / d
                N_eff += (1 + w_pi) * n_empty / d**2 * uf_scale

            # Exchange
            xp_under = comb(min(pc, d), 2)
            xp_over = comb(max(pc - d, 0), 2)
            N_eff -= xp_under / (d**2 * n)
            N_eff -= xp_over / (d**2 * n_ref)
            N_eff -= xp_over * (n - n_ref) / (d * (d - 1) * n)

            # Mode-vacancy coupling
            N_eff += pa * (d - pa) / d**3

            alpha = (d + w_pi * N_eff - 1) / d**2
            alpha -= w_pi * pen

        Z_eff = max(Z_net, 0.01) ** alpha
        E_pred = (Z_eff / n)**2 * E_H
        err = (E_pred - E_ion) / E_ion * 100
        results.append((Z, sym, E_ion, E_pred, err, alpha))
    return results


# === RUN AND DISPLAY ===

results = calc_all()

print("=" * 85)
print("  v18: + valence-dependent f-screening (w_f_p = w_pi)")
print(f"  f-screen: s-val→w_delta/d={W_F['s']:.4f}  p-val→w_pi={W_F['p']:.4f}  "
      f"d-val→w_delta/d={W_F['d']:.4f}")
print("=" * 85)

print(f"\n  {'Z':>3} {'Sym':<3} {'E_obs':>8} {'E_pred':>8} {'err%':>7}")
for r in results:
    Z, sym, E_obs, E_pred, err, alpha = r
    flag = ' ***' if abs(err) > 10 else (' **' if abs(err) > 5 else '')
    if Z in [3, 11, 19, 21, 31, 37, 39, 49, 55, 57, 72, 81]:
        print()
    if Z == 37: print("  --- PERIOD 5 (PREDICTION) ---")
    if Z == 55: print("  --- PERIOD 6 ---")
    if Z == 57: print("  --- LANTHANIDES ---")
    if Z == 72: print("  --- 5d TMs ---")
    if Z == 81: print("  --- 6p ---")
    print(f"  {Z:3d} {sym:<3} {E_obs:8.3f} {E_pred:8.3f} {err:+7.1f}%{flag}")

# === STATS ===
def stats(label, rlist):
    errs = [abs(r[4]) for r in rlist]
    return (f"    {label:>18}: {np.mean(errs):5.2f}%  max={max(errs):5.1f}%  "
            f"<5%:{sum(1 for e in errs if e<5)}/{len(errs)}")

all_e = [abs(r[4]) for r in results]
print(f"\n{'='*85}")
print(f"  STATISTICS:")
print(stats("Periods 1-4", [r for r in results if r[0] <= 36]))
print(stats("Period 5", [r for r in results if 37 <= r[0] <= 54]))
print(stats("Period 6", [r for r in results if r[0] >= 55]))
print(stats("ALL (86 atoms)", results))
print()
print(stats("3d TM (Sc-Zn)", [r for r in results if 21 <= r[0] <= 30]))
print(stats("4d TM (Y-Cd)", [r for r in results if 39 <= r[0] <= 48]))
print(stats("5d TM (Hf-Hg)", [r for r in results if 72 <= r[0] <= 80]))
print(stats("Lanthanides", [r for r in results if 57 <= r[0] <= 71]))
print(stats("2p (B-Ne)", [r for r in results if 5 <= r[0] <= 10]))
print(stats("3p (Al-Ar)", [r for r in results if 13 <= r[0] <= 18]))
print(stats("4p (Ga-Kr)", [r for r in results if 31 <= r[0] <= 36]))
print(stats("5p (In-Xe)", [r for r in results if 49 <= r[0] <= 54]))
print(stats("6p (Tl-Rn)", [r for r in results if 81 <= r[0] <= 86]))

p15 = [r for r in results if r[0] <= 54]
p6 = [r for r in results if r[0] >= 55]
print(f"\n    P1-5: {np.mean([abs(r[4]) for r in p15]):.2f}%")
print(f"    P6:   {np.mean([abs(r[4]) for r in p6]):.2f}%")

print(f"\n  PROGRESSION:")
print(f"    v15: 2.92% (54 atoms)")
print(f"    v16: 1.96% (54 atoms), 3.89% (86 atoms)")
print(f"    v17: 2.08% (54 atoms), 4.09% (86 atoms)")
print(f"    v18: {np.mean([abs(r[4]) for r in p15]):.2f}% (54 atoms), "
      f"{np.mean(all_e):.2f}% (86 atoms)")
print(f"\n    <5%: {sum(1 for e in all_e if e<5)}/86  "
      f"<10%: {sum(1 for e in all_e if e<10)}/86  "
      f"max: {max(all_e):.1f}%")

# Remaining outliers
print(f"\n  REMAINING OUTLIERS (>10%):")
for r in results:
    if abs(r[4]) > 10:
        print(f"    {r[1]:3s} Z={r[0]:2d}: {r[4]:+.1f}%")

print(f"{'='*85}")
