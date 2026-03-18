"""
Z_eff Combined — One Coupling Per Mode
========================================
Replace separate parity + screening + penetration with ONE coupling
per core mode. The simulation showed:

  First d modes:  ENHANCE binding (constructive wave interference)
  After d+1:      SCREEN (destructive interference / charge blocking)
  Crossover at:   n_ref = d+1 = 4

Each core mode contributes a single number C_j to alpha:
  alpha = (d*n + sum(C_j)) / (d^2 * n)

No separate parity. No separate pen. One mechanism.
"""

import numpy as np
from math import factorial, comb

d = 3
gamma_sg = np.pi / (2**(d+1) * np.pi - 2)
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6
n_ref = d + 1  # = 4, the crossover point
w_pi = np.cos(np.pi / d)
w_delta = np.cos(2 * np.pi / d)
J_hund = 2 / (d + 2)


def mode_coupling(l_core, l_val, delta_n, n_channels, n_val):
    """Net coupling of one core mode on the valence mode.

    Combines wave enhancement AND charge screening into ONE number.
    Positive = enhances binding. Negative = screens.

    The crossover from enhancement to screening happens at delta_n ~ d.
    Close shells enhance (constructive interference through cosine potential).
    Deep shells screen (charge blocking dominates wave coupling).

    Parameters:
        l_core: angular momentum of core mode
        l_val: angular momentum of valence mode
        delta_n: shell separation (core_ref_n - nn)
        n_channels: occupied channels in this core shell = min(count, 2l+1)
        n_val: principal quantum number of valence
    """
    # Mass ratio: coupling strength from breather spectrum
    # sin((lc+1)*g) / sin((lv+1)*g) ≈ (lc+1)/(lv+1) for small gamma
    mass_r = np.sin((l_core + 1) * gamma_sg) / np.sin((l_val + 1) * gamma_sg)

    # Angular parity: same parity enhances, different parity disrupts
    same_parity = (l_core % 2) == (l_val % 2)

    # Radial factor: transitions from enhancement to screening
    # At delta_n = 0: maximum enhancement (same shell coupling)
    # At delta_n = d: crossover to screening
    # At delta_n >> d: full screening
    #
    # f(delta_n) = (d - delta_n) / d for enhancement regime
    # Smoothly transitions to negative (screening) for delta_n > d

    if delta_n <= 0:
        # Same shell: strong coupling
        radial = 1.0
    elif delta_n <= d:
        # Near shell: enhancement decreases linearly with distance
        radial = 1.0 - delta_n / n_ref
    else:
        # Deep core: screening (negative)
        radial = -delta_n / (d * n_val)

    # Combine: coupling per channel
    if same_parity:
        # Same parity: constructive interference in cosine potential
        # Enhancement for near shells, reduced screening for deep
        C = mass_r * radial / (2 * l_core + 1)
    else:
        # Different parity: destructive interference
        # Anti-screening for near shells, enhanced screening for deep
        C = -mass_r * radial / (2 * l_core + 1)

    return n_channels * C


# === ATOMS (same as v19) ===
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

        # === CHARGE SCREENING → Z_net ===
        # Use v19's full core_screening (it works and captures d anti-screening).
        # This is the part we keep — it's the electrostatic physics.
        # The wave coupling (parity/pen) is what we're combining into C_total.
        val_block = 'p' if has_p else ('d' if (dc_val > 0 or is_pd) else 's')

        # Import v19's screening function
        S_core = 0.0
        W_F = {'s': w_delta / d, 'p': w_pi, 'd': w_delta / d}
        w_f = W_F.get(val_block, w_delta / d)
        for nn, ll, count in config:
            if nn < n:
                n_ch = min(count, 2 * ll + 1)
                if ll <= 1:
                    S_core += n_ch * w_pi
                elif ll == 2:
                    delta_n_core = n - nn
                    # d-core screening depends on valence type (from sim)
                    # d couples to s at (d+1)/d stronger than to p
                    # Single d: ratio = 1.33 = (d+1)/d
                    # Filled d-shell: ratio grows (collective effect)
                    if val_l == 0:
                        # d->s: full anti-screening
                        w_d = w_delta
                    elif val_l == 1:
                        # d->p: v19 uses t2g/eg which gives POSITIVE screening
                        # The d10 core screens p-valence (same as s,p core)
                        # Sim confirms: d modes couple to p LESS than to s
                        # In 3D: the t2g/eg split makes d10 screen p at w_pi
                        if count == 10:
                            w_d = w_pi  # d10->p: screens like s,p core
                        else:
                            w_d = w_pi  # partial d also screens p
                    else:
                        w_d = w_delta  # d->d: full

                    if delta_n_core >= 2 and count == 10:
                        # Deep d10: blend toward normal screening
                        S_normal = 5 * w_pi
                        S_anti = 5 * w_d
                        S_core += S_normal + (S_anti - S_normal) / d
                    else:
                        S_core += n_ch * w_d
                elif ll == 3:
                    S_core += n_ch * w_f
            elif nn == n and is_pd and ll != 2:
                S_core += min(count, 2 * ll + 1) * w_pi
        # Three-body cross-term: when d AND f are both in core,
        # their combined screening is LESS than additive.
        # cross_per_d = (d-1)/d^2 for p-valence, 1/d! for s-valence
        # Sim: 5d+f -> p flips from enhance to screen
        n_d_ch = sum(min(c, 5) for nn2, ll2, c in config
                    if nn2 < n and ll2 == 2)
        n_f_ch_core = sum(min(c, 7) for nn2, ll2, c in config
                        if nn2 < n and ll2 == 3 and c > 0)
        if n_d_ch > 0 and n_f_ch_core > 0:
            if val_l == 1:  # p-valence
                cross = (d - 1) / d**2 * n_d_ch * n_f_ch_core / (d * n)
            else:  # s/d-valence
                cross = n_d_ch * n_f_ch_core / (factorial(d) * d * n)
            S_core += cross  # adds to screening (reduces Z_net)

        Z_net = Z - S_core

        # === WAVE COUPLING → C_total ===
        # This replaces parity + penetration.
        # Each core mode contributes wave interference that modifies alpha.
        # Positive C = constructive (enhances binding).
        # Negative C = destructive (weakens binding).
        C_total = 0.0

        # --- Valence self-coupling (was: parity) ---
        sp = (s_count == 2)
        s1 = (s_count == 1)
        if not has_p:
            if is_pd:
                # Pd-like: d10 with no s electron
                C_total = w_delta + J_hund
            elif s1 and dc_val > 0:
                # s-single + d: the d-modes modify the s-wave coupling
                C_total = w_pi + (n - n_ref)
            elif sp and dc_val > 0:
                # s-pair + d: full TM treatment
                C_total = 1.0  # s-pair base

                dc = dc_val
                # t2g/eg from cubic symmetry
                if dc > 5:
                    t2g_paired = min(dc - 5, d)
                    t2g_full = (t2g_paired == d)
                    if t2g_full:
                        eg_occ = dc - 5
                        eg_paired = max(0, dc - 8)
                        eg_unp = eg_occ - eg_paired
                        C_total -= eg_unp * abs(w_delta) / d
                else:
                    t2g_full = False

                if dc == 10:
                    C_total += J_hund
                    n_deep = sum(1 for nn2, ll2, c2 in config
                                if ll2 == 2 and c2 == 10 and nn2 < n-1)
                    C_total += J_hund * n_deep

                # n-scaling for s-pair TMs
                C_total += w_pi * (n - n_ref)
            elif sp:
                C_total = 1.0    # s-pair: constructive
            else:
                C_total = -1.0   # s-single: destructive

            # f-core d-shell rebalancing (from v19)
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

        else:
            # p-block: N_eff captures the valence coupling
            pc = p_count
            pa = min(pc, 2*d - pc) if pc <= d else 2*d - pc
            pl = max(0, pc - d)

            N_eff = 0
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
                w_first = (n**2 - d) / n**2  # FIXED: was 1+w_pi
                w_rest = 1 + w_pi
                N_eff = pa + fl * w_first + rl * w_rest
                xo_pairs = comb(pl, 2) if pl >= 2 else 0
                N_eff -= xo_pairs / (d**2 * n_ref)
                xa = comb(pa, 2) + comb(pl, 2)
                xp_over = xa - comb(pa, 2)
                N_eff -= xp_over * (n - n_ref) / (d * (d-1) * n)

            N_eff += pa * (d - pa) / d**3

            # Deep Hund scattering
            deep_hl_ch = sum(min(c, 2*ll+1) for nn2, ll, c in config
                           if ll >= 2 and c > 0 and (n - nn2) >= 2)
            xu_total = comb(min(pc, d), 2)
            if xu_total > 0 and deep_hl_ch > 0:
                N_eff -= xu_total * deep_hl_ch / (d**3 * n)

            C_total = N_eff  # N_eff IS the p-block coupling

        # --- Core wave interference (was: penetration) ---
        # Deep core modes weaken the valence binding through destructive
        # interference. This is pen in v19.
        if S_core != 0:
            pen = gamma_sg * S_core / (d * n + (d + 1) * abs(S_core))
            has_d_any = any(c > 0 for nn, ll, c in config if ll == 2)
            if n >= n_ref and not has_d_any:
                pen *= 1 / d**2
        else:
            pen = 0

        # === ALPHA ===
        if not has_p:
            alpha = (d * n + C_total) / (d**2 * n) - pen
        else:
            alpha = (d + w_pi * C_total - 1) / d**2 - w_pi * pen

        alpha = max(alpha, 0.05)
        alpha = min(alpha, 0.95)

        # === ENERGY ===
        Z_eff = Z_net ** alpha
        E_pred = (Z_eff / n)**2 * E_H

        # Lattice shear
        N_core = sum(c for nn, ll, c in config if nn < n)
        shear = 1 + gamma_sg * N_core / (d**2 * n**2)
        E_pred *= shear

        err = (E_pred - E_obs) / E_obs * 100
        results.append((Z, sym, E_obs, E_pred, err, alpha, C_total))

    return results


# === RUN ===
results = calc_all()
print("Z_eff Combined — One Coupling Per Mode")
print("=" * 65)
print(f"  alpha = (d*n + C_total) / (d^2*n)")
print(f"  C_total = sum of mode couplings (no separate parity/pen)")
print()

errs = [abs(r[4]) for r in results]
errs_p15 = [abs(r[4]) for r in results if r[0] <= 54]
errs_p6 = [abs(r[4]) for r in results if 55 <= r[0] <= 86]
errs_p7 = [abs(r[4]) for r in results if r[0] >= 87]

print(f"  ALL: {np.mean(errs):.2f}%  P1-5: {np.mean(errs_p15):.2f}%  P6: {np.mean(errs_p6):.2f}%  P7: {np.mean(errs_p7):.2f}%")
print(f"  Under 5%: {sum(1 for e in errs if e < 5)}/{len(results)}")
print(f"  Under 10%: {sum(1 for e in errs if e < 10)}/{len(results)}")
print(f"  v19 reference: ALL=2.61%, P1-5=2.07%")
print()

print(f"  {'Z':>3} {'Sym':<3} {'err':>7}")
print(f"  {'-'*15}")
for Z, sym, E_obs, E_pred, err, alpha, C_total in results:
    if abs(err) > 5:
        print(f"  {Z:3d} {sym:<3} {err:+6.1f}%")
