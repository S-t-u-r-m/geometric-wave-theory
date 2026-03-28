"""
N-Body Correction from Oh Tensor Products
==========================================
For N breather modes on the d=3 cubic lattice:
1. Each mode transforms as an Oh irrep
2. The N-fold tensor product decomposes into Oh irreps
3. The A1g (symmetric) content = the "coherent" coupling fraction
4. This fraction IS the N-body correction — analogous to Wyler's alpha

Character table of Oh (order 48):
Classes: E(1), C3(8), C2(6), C4(6), C2'(3), i(1), S6(8), σd(6), S4(6), σh(3)
"""

import numpy as np
from math import factorial

d = 3

# Oh character table
# Rows: irreps [A1g, A2g, Eg, T1g, T2g, A1u, A2u, Eu, T1u, T2u]
# Cols: classes [E, C3, C2, C4, C2', i, S6, σd, S4, σh]
class_sizes = np.array([1, 8, 6, 6, 3, 1, 8, 6, 6, 3])
Oh_order = 48

chars = {
    'A1g': np.array([1,  1,  1,  1,  1,  1,  1,  1,  1,  1]),
    'A2g': np.array([1,  1, -1, -1,  1,  1,  1, -1, -1,  1]),
    'Eg':  np.array([2, -1,  0,  0,  2,  2, -1,  0,  0,  2]),
    'T1g': np.array([3,  0, -1,  1, -1,  3,  0, -1,  1, -1]),
    'T2g': np.array([3,  0,  1, -1, -1,  3,  0,  1, -1, -1]),
    'A1u': np.array([1,  1,  1,  1,  1, -1, -1, -1, -1, -1]),
    'A2u': np.array([1,  1, -1, -1,  1, -1, -1,  1,  1, -1]),
    'Eu':  np.array([2, -1,  0,  0,  2, -2,  1,  0,  0, -2]),
    'T1u': np.array([3,  0, -1,  1, -1, -3,  0,  1, -1,  1]),
    'T2u': np.array([3,  0,  1, -1, -1, -3,  0, -1,  1,  1]),
}

dims = {k: v[0] for k, v in chars.items()}


def tensor_product_char(irrep_list):
    """Character of the tensor product of multiple Oh irreps.
    χ_product(g) = product of χ_i(g) for each irrep.
    """
    result = np.ones(10)
    for irrep in irrep_list:
        result *= chars[irrep]
    return result


def decompose(product_char):
    """Decompose a character into Oh irreps.
    mult(Γ) = (1/|Oh|) × Σ_classes size(class) × χ_Γ*(class) × χ_product(class)
    """
    decomp = {}
    for irrep_name, irrep_char in chars.items():
        mult = np.sum(class_sizes * irrep_char * product_char) / Oh_order
        mult = round(mult)
        if mult > 0:
            decomp[irrep_name] = mult
    return decomp


def a1g_fraction(irrep_list):
    """Compute the A1g content of the tensor product as a fraction.

    A1g fraction = dim(A1g component) / total_dim
    = mult(A1g) / product(dims)
    """
    prod_char = tensor_product_char(irrep_list)
    decomp = decompose(prod_char)

    total_dim = 1
    for irrep in irrep_list:
        total_dim *= dims[irrep]

    a1g_mult = decomp.get('A1g', 0)

    return a1g_mult, total_dim, a1g_mult / total_dim if total_dim > 0 else 0, decomp


def atom_irreps(config):
    """Map atomic configuration to Oh irreps.

    s-modes (l=0): A1g (1 per s-electron pair or single)
    p-modes (l=1): T1u (1 per p-direction: x, y, z)
    d-modes (l=2): T2g (xy, xz, yz) + Eg (z², x²-y²)
    f-modes (l=3): A2u + T1u + T2u

    Returns list of irreps, one per occupied channel.
    """
    val_n = max(nn for nn, ll, c in config)
    irreps = []

    for nn, ll, count in config:
        if nn != val_n and not (nn == val_n - 1 and ll == 2):
            continue  # only valence modes

        if ll == 0:
            # s-channel: A1g
            for i in range(min(count, 1)):
                irreps.append('A1g')
            if count >= 2:
                irreps.append('A1g')  # paired
        elif ll == 1:
            # p-channels: T1u for each occupied direction
            for i in range(min(count, 3)):
                irreps.append('T1u')
            # Paired (overfill): same irrep again
            for i in range(max(count - 3, 0)):
                irreps.append('T1u')
        elif ll == 2:
            # d-channels: T2g (first 3) + Eg (next 2)
            t2g_count = min(count, 6)
            eg_count = max(count - 6, 0)
            for i in range(min(t2g_count, 3)):
                irreps.append('T2g')
            if t2g_count > 3:
                for i in range(t2g_count - 3):
                    irreps.append('T2g')  # paired t2g
            for i in range(min(eg_count, 2)):
                irreps.append('Eg')
            if eg_count > 2:
                for i in range(eg_count - 2):
                    irreps.append('Eg')  # paired eg

    return irreps


# === TEST ===
print("N-Body Correction from Oh Tensor Products")
print("=" * 55)
print()

# Test simple atoms
test_atoms = [
    ('H',  1, 13.598, [(1,0,1)]),
    ('He', 2, 24.587, [(1,0,2)]),
    ('Li', 3, 5.392,  [(1,0,2),(2,0,1)]),
    ('Be', 4, 9.323,  [(1,0,2),(2,0,2)]),
    ('B',  5, 8.298,  [(1,0,2),(2,0,2),(2,1,1)]),
    ('C',  6, 11.260, [(1,0,2),(2,0,2),(2,1,2)]),
    ('N',  7, 14.534, [(1,0,2),(2,0,2),(2,1,3)]),
    ('O',  8, 13.618, [(1,0,2),(2,0,2),(2,1,4)]),
    ('F',  9, 17.423, [(1,0,2),(2,0,2),(2,1,5)]),
    ('Ne', 10,21.565, [(1,0,2),(2,0,2),(2,1,6)]),
]

# Known errors from our formula
errors = {'H':0.0,'He':2.5,'Li':4.5,'Be':-4.0,'B':-1.8,'C':-2.8,
          'N':-0.1,'O':-1.3,'F':-2.2,'Ne':2.4}

print(f"{'Sym':>3} {'irreps':>30} {'A1g':>5} {'total':>7} {'frac':>8} {'err%':>6}")
print("-" * 65)

for sym, Z, IE_obs, config in test_atoms:
    irreps = atom_irreps(config)

    if len(irreps) == 0:
        print(f"  {sym:>3} (no valence irreps)")
        continue

    a1g, total, frac, decomp = a1g_fraction(irreps)
    err = errors.get(sym, 0)

    irrep_str = '+'.join(irreps) if len(irreps) <= 6 else f"{len(irreps)} modes"
    print(f"  {sym:>3} {irrep_str:>30} {a1g:5d} {total:7d} {frac:8.4f} {err:+5.1f}%")

print()
print("A1g fraction = how much of the coupling space is symmetric.")
print("Higher fraction = more coherent coupling = more accurate formula?")
print()

# Does the fraction correlate with error?
fracs = []
errs_list = []
for sym, Z, IE_obs, config in test_atoms:
    irreps = atom_irreps(config)
    if len(irreps) == 0: continue
    a1g, total, frac, _ = a1g_fraction(irreps)
    fracs.append(frac)
    errs_list.append(abs(errors.get(sym, 0)))

if len(fracs) > 2:
    corr = np.corrcoef(fracs, errs_list)[0,1]
    print(f"Correlation(A1g_fraction, |error|): r = {corr:+.4f}")
    print(f"  Negative = higher fraction -> lower error (expected)")
