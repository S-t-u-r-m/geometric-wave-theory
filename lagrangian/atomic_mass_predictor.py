#!/usr/bin/env python3
"""
Aether Framework — Atomic Mass Predictor
=========================================
Atoms are single unified standing waves in the aether.
When nucleons merge into one coherent waveform, energy is released
(binding energy). Mass = wave energy / c².

M(A,Z) = Z·m_H + N·m_n − B(A,Z)/c²

Binding energy from the framework's standing wave formula (SEMF):
  B = a_V·A − a_S·A^(2/3) − a_C·Z(Z−1)/A^(1/3) − a_A·(A−2Z)²/A + δ

Framework coefficients (from Aether_Equations_Quick_Reference.md):
  a_V = 15.56 MeV   (volume:     bulk wave energy per nucleon)
  a_S = 17.23 MeV   (surface:    bubble boundary cost)
  a_C = 0.697 MeV   (Coulomb:    same-charge mode repulsion)
  a_A = 23.29 MeV   (asymmetry:  yin-yang mode balance penalty)
  δ   = ±12/√A      (pairing:    paired modes resonate better)
"""

import math

# ── Physical constants ──────────────────────────────────────────────────────
m_H       = 938.783    # Hydrogen atom mass (MeV/c²) = proton + electron
m_n       = 939.565    # Neutron mass (MeV/c²)
MeV_per_u = 931.494    # MeV per atomic mass unit (u)

# ── Framework SEMF coefficients ─────────────────────────────────────────────
a_V = 15.56
a_S = 17.23
a_C = 0.697
a_A = 23.29

def pairing(A, Z):
    """Pairing term δ = ±12/√A  (paired standing wave modes resonate better)"""
    N = A - Z
    if A % 2 == 1:                          # odd-A: no pairing effect
        return 0.0
    elif Z % 2 == 0 and N % 2 == 0:        # even-even: favorable
        return 12.0 / math.sqrt(A)
    else:                                    # odd-odd: unfavorable
        return -12.0 / math.sqrt(A)

def binding_energy(A, Z):
    """Standing wave binding energy B(A,Z) in MeV."""
    if A < 2 or Z < 1 or Z >= A:
        return 0.0
    N = A - Z
    B = (a_V * A
         - a_S * A**(2/3)
         - a_C * Z * (Z - 1) / A**(1/3)
         - a_A * (A - 2*Z)**2 / A
         + pairing(A, Z))
    return max(B, 0.0)

def predict_mass(A, Z):
    """Predicted atomic mass in u.  M = (Z·m_H + N·m_n − B) / MeV_per_u"""
    N = A - Z
    return (Z * m_H + N * m_n - binding_energy(A, Z)) / MeV_per_u

# ── Element data ────────────────────────────────────────────────────────────
# (symbol, name, Z, A_ref, measured_mass_u)
# A_ref  = most abundant stable isotope (or most stable for synthetic elements)
# mass_u = measured atomic mass from AME2020/NIST  (u)
# Note: for Z > 100 masses are approximate / estimated from systematics
elements = [
    ("H",  "Hydrogen",       1,   1,   1.00782503),
    ("He", "Helium",         2,   4,   4.00260325),
    ("Li", "Lithium",        3,   7,   7.01600344),
    ("Be", "Beryllium",      4,   9,   9.01218307),
    ("B",  "Boron",          5,  11,  11.00930536),
    ("C",  "Carbon",         6,  12,  12.00000000),
    ("N",  "Nitrogen",       7,  14,  14.00307401),
    ("O",  "Oxygen",         8,  16,  15.99491462),
    ("F",  "Fluorine",       9,  19,  18.99840316),
    ("Ne", "Neon",          10,  20,  19.99244018),
    ("Na", "Sodium",        11,  23,  22.98976928),
    ("Mg", "Magnesium",     12,  24,  23.98504170),
    ("Al", "Aluminium",     13,  27,  26.98153853),
    ("Si", "Silicon",       14,  28,  27.97692653),
    ("P",  "Phosphorus",    15,  31,  30.97376200),
    ("S",  "Sulfur",        16,  32,  31.97207100),
    ("Cl", "Chlorine",      17,  35,  34.96885268),
    ("Ar", "Argon",         18,  40,  39.96238312),
    ("K",  "Potassium",     19,  39,  38.96370668),
    ("Ca", "Calcium",       20,  40,  39.96259086),
    ("Sc", "Scandium",      21,  45,  44.95590828),
    ("Ti", "Titanium",      22,  48,  47.94794198),
    ("V",  "Vanadium",      23,  51,  50.94395778),
    ("Cr", "Chromium",      24,  52,  51.94050623),
    ("Mn", "Manganese",     25,  55,  54.93804391),
    ("Fe", "Iron",          26,  56,  55.93493633),
    ("Co", "Cobalt",        27,  59,  58.93319429),
    ("Ni", "Nickel",        28,  58,  57.93534241),
    ("Cu", "Copper",        29,  63,  62.92959747),
    ("Zn", "Zinc",          30,  64,  63.92914201),
    ("Ga", "Gallium",       31,  69,  68.92557354),
    ("Ge", "Germanium",     32,  74,  73.92117776),
    ("As", "Arsenic",       33,  75,  74.92159457),
    ("Se", "Selenium",      34,  80,  79.91652176),
    ("Br", "Bromine",       35,  79,  78.91833710),
    ("Kr", "Krypton",       36,  84,  83.91149773),
    ("Rb", "Rubidium",      37,  85,  84.91178974),
    ("Sr", "Strontium",     38,  88,  87.90561226),
    ("Y",  "Yttrium",       39,  89,  88.90584033),
    ("Zr", "Zirconium",     40,  90,  89.90469018),
    ("Nb", "Niobium",       41,  93,  92.90637309),
    ("Mo", "Molybdenum",    42,  98,  97.90540490),
    ("Tc", "Technetium",    43,  98,  97.90721200),
    ("Ru", "Ruthenium",     44, 102, 101.90434930),
    ("Rh", "Rhodium",       45, 103, 102.90550400),
    ("Pd", "Palladium",     46, 106, 105.90348300),
    ("Ag", "Silver",        47, 107, 106.90509700),
    ("Cd", "Cadmium",       48, 114, 113.90336510),
    ("In", "Indium",        49, 115, 114.90387780),
    ("Sn", "Tin",           50, 120, 119.90220163),
    ("Sb", "Antimony",      51, 121, 120.90381200),
    ("Te", "Tellurium",     52, 130, 129.90622748),
    ("I",  "Iodine",        53, 127, 126.90447277),
    ("Xe", "Xenon",         54, 132, 131.90415400),
    ("Cs", "Cesium",        55, 133, 132.90545196),
    ("Ba", "Barium",        56, 138, 137.90524700),
    ("La", "Lanthanum",     57, 139, 138.90635227),
    ("Ce", "Cerium",        58, 140, 139.90543910),
    ("Pr", "Praseodymium",  59, 141, 140.90764800),
    ("Nd", "Neodymium",     60, 142, 141.90772910),
    ("Pm", "Promethium",    61, 145, 144.91274880),
    ("Sm", "Samarium",      62, 152, 151.91973270),
    ("Eu", "Europium",      63, 153, 152.92123800),
    ("Gd", "Gadolinium",    64, 158, 157.92410130),
    ("Tb", "Terbium",       65, 159, 158.92534680),
    ("Dy", "Dysprosium",    66, 164, 163.92917510),
    ("Ho", "Holmium",       67, 165, 164.93032800),
    ("Er", "Erbium",        68, 166, 165.93029500),
    ("Tm", "Thulium",       69, 169, 168.93421800),
    ("Yb", "Ytterbium",     70, 174, 173.93886210),
    ("Lu", "Lutetium",      71, 175, 174.94077700),
    ("Hf", "Hafnium",       72, 180, 179.94655010),
    ("Ta", "Tantalum",      73, 181, 180.94799510),
    ("W",  "Tungsten",      74, 184, 183.95093310),
    ("Re", "Rhenium",       75, 187, 186.95575020),
    ("Os", "Osmium",        76, 192, 191.96148110),
    ("Ir", "Iridium",       77, 193, 192.96292600),
    ("Pt", "Platinum",      78, 195, 194.96479410),
    ("Au", "Gold",          79, 197, 196.96656870),
    ("Hg", "Mercury",       80, 202, 201.97064300),
    ("Tl", "Thallium",      81, 205, 204.97442760),
    ("Pb", "Lead",          82, 208, 207.97665220),
    ("Bi", "Bismuth",       83, 209, 208.98039910),
    ("Po", "Polonium",      84, 209, 208.98243080),
    ("At", "Astatine",      85, 210, 209.98714800),
    ("Rn", "Radon",         86, 222, 222.01757820),
    ("Fr", "Francium",      87, 223, 223.01973600),
    ("Ra", "Radium",        88, 226, 226.02541030),
    ("Ac", "Actinium",      89, 227, 227.02775200),
    ("Th", "Thorium",       90, 232, 232.03805600),
    ("Pa", "Protactinium",  91, 231, 231.03588430),
    ("U",  "Uranium",       92, 238, 238.05078826),
    ("Np", "Neptunium",     93, 237, 237.04817340),
    ("Pu", "Plutonium",     94, 244, 244.06420430),
    ("Am", "Americium",     95, 243, 243.06138130),
    ("Cm", "Curium",        96, 247, 247.07035400),
    ("Bk", "Berkelium",     97, 247, 247.07030700),
    ("Cf", "Californium",   98, 251, 251.07958700),
    ("Es", "Einsteinium",   99, 252, 252.08298000),
    ("Fm", "Fermium",      100, 257, 257.09510500),
    ("Md", "Mendelevium",  101, 258, 258.09843100),
    ("No", "Nobelium",     102, 259, 259.10103000),
    ("Lr", "Lawrencium",   103, 262, 262.10963000),
    ("Rf", "Rutherfordium",104, 267, 267.12179000),
    ("Db", "Dubnium",      105, 268, 268.12544000),
    ("Sg", "Seaborgium",   106, 271, 271.13347000),
    ("Bh", "Bohrium",      107, 272, 272.13803000),
    ("Hs", "Hassium",      108, 269, 269.13361000),
    ("Mt", "Meitnerium",   109, 278, 278.15631000),
    ("Ds", "Darmstadtium", 110, 281, 281.16451000),
    ("Rg", "Roentgenium",  111, 282, 282.16912000),
    ("Cn", "Copernicium",  112, 285, 285.17712000),
    ("Nh", "Nihonium",     113, 286, 286.18221000),
    ("Fl", "Flerovium",    114, 289, 289.19042000),
    ("Mc", "Moscovium",    115, 290, 290.19659000),
    ("Lv", "Livermorium",  116, 293, 293.20449000),
    ("Ts", "Tennessine",   117, 294, 294.21073000),
    ("Og", "Oganesson",    118, 294, 294.21398000),
]

# ── Run predictions ─────────────────────────────────────────────────────────
print("AETHER FRAMEWORK — ATOMIC MASS PREDICTIONS")
print("Atoms as unified standing waves in the aether")
print("M(A,Z) = Z·m_H + N·m_n − B(A,Z)/c²")
print()
print(f"{'Z':>3} {'Sym':<4} {'Name':<16} {'A':>3}  "
      f"{'Predicted (u)':>14} {'Measured (u)':>14} {'Δ (MeV)':>9} {'Δ (%)':>7}  "
      f"{'B/A pred':>9} {'B/A meas':>9}")
print("─" * 100)

errors_light  = []   # Z = 1–9
errors_medium = []   # Z = 10–83
errors_heavy  = []   # Z = 84–118

for sym, name, Z, A, mass_meas in elements:
    mass_pred   = predict_mass(A, Z)
    delta_u     = mass_pred - mass_meas
    delta_MeV   = delta_u * MeV_per_u
    delta_pct   = (delta_u / mass_meas) * 100

    B_pred = binding_energy(A, Z)
    B_meas = (Z * m_H + (A - Z) * m_n) - mass_meas * MeV_per_u
    BA_pred = B_pred / A if A > 0 else 0
    BA_meas = B_meas / A if A > 0 else 0

    flag = ""
    if abs(delta_pct) > 0.1:
        flag = " ◄"   # highlight large deviations

    print(f"{Z:>3} {sym:<4} {name:<16} {A:>3}  "
          f"{mass_pred:>14.8f} {mass_meas:>14.8f} "
          f"{delta_MeV:>+9.3f} {delta_pct:>+7.4f}%  "
          f"{BA_pred:>9.4f} {BA_meas:>9.4f}{flag}")

    if Z <= 9:
        errors_light.append(abs(delta_pct))
    elif Z <= 83:
        errors_medium.append(abs(delta_pct))
    else:
        errors_heavy.append(abs(delta_pct))

# ── Summary statistics ───────────────────────────────────────────────────────
print()
print("─" * 100)
print("SUMMARY")
print(f"  Z =   1–9   (light nuclei, SEMF less accurate):  avg |Δ| = {sum(errors_light)/len(errors_light):.4f}%")
print(f"  Z =  10–83  (main periodic table):               avg |Δ| = {sum(errors_medium)/len(errors_medium):.4f}%")
print(f"  Z =  84–118 (heavy / superheavy):                avg |Δ| = {sum(errors_heavy)/len(errors_heavy):.4f}%")
all_errors = errors_light + errors_medium + errors_heavy
print(f"  All 118 elements:                                 avg |Δ| = {sum(all_errors)/len(all_errors):.4f}%")
print()
print("Note: Δ is error in TOTAL atomic mass (u). Binding energy errors are larger.")
print("Framework claim: < 1% for A > 10  ← check B/A columns above.")
