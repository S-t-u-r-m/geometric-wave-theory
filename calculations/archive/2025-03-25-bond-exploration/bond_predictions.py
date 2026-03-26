"""Quick bond prediction comparison: Hessian multimode vs V8."""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

PI = np.pi; d = 3; gamma = PI / (2**(d+1)*PI - 2)
N = 512; x = np.arange(N, dtype=np.float64); center = N//2; kw = 3

def kink(x, pos):
    return (4.0/PI)*(np.arctan(np.exp(x-pos+kw/2))-np.arctan(np.exp(x-pos-kw/2)))

def hess_ev(phi, k=10):
    diag = 2.0 + np.cos(PI*phi)
    off = -np.ones(N-1)
    H = sparse.diags([off,diag,off],[-1,0,1],shape=(N,N),format='lil')
    H[0,N-1]=-1.0; H[N-1,0]=-1.0
    ev,_ = eigsh(H.tocsr(), k=k, which='SM')
    return np.sort(ev)

phi_s = kink(x, center)
ev_s = hess_ev(phi_s, k=10)

# Compute V_n(R) for each mode
R_range = list(range(4, 50))
V_modes = {0: {}, 1: {}, 2: {}}
for R in R_range:
    pA = center - R//2; pB = center + R//2
    phi_d = kink(x, pA) + kink(x, pB)
    ev_d = hess_ev(phi_d, k=8)
    for i in range(3):
        V_modes[i][R] = ev_d[2*i] - ev_s[i]

# V_total at key R values
print("V_total at key R values:")
print(f"{'R':>3} {'V(1occ)':>10} {'V(2occ)':>10} {'V(3occ)':>10}")
for R in range(4, 12):
    v1 = V_modes[0][R]
    v2 = V_modes[0][R] + V_modes[1][R]
    v3 = V_modes[0][R] + V_modes[1][R] + V_modes[2][R]
    print(f"{R:3d} {v1:+10.6f} {v2:+10.6f} {v3:+10.6f}")

print()

# R_eq for each occupation
for n_occ in [1, 2, 3]:
    V_total = {R: sum(V_modes[i][R] for i in range(n_occ)) for R in R_range}
    best_R = min(V_total, key=V_total.get)
    De = -V_total[best_R]
    print(f"{n_occ}-mode: R_eq={best_R}, D_e={De:.6f}")

print()

# Effective couplings
V0_ref = -V_modes[0][7]
c1 = -V_modes[0][7] / V0_ref
c2 = -(V_modes[0][7] + V_modes[1][7]) / V0_ref
c3 = -(V_modes[0][6] + V_modes[1][6] + V_modes[2][6]) / V0_ref

print(f"Effective couplings:")
print(f"  bo=1 (R=7): {c1:.4f}  (V8: 1.0)")
print(f"  bo=2 (R=7): {c2:.4f}  (V8: 1.5)")
print(f"  bo=3 (R=6): {c3:.4f}  (V8: 2.0)")
print()

C_BOND = PI/d**2
cpl_map = {1: c1, 2: c2, 3: c3}

mols = [
    ("H2",   13.598, 1, 4.478),
    ("N2",   14.534, 3, 9.759),
    ("O2",   13.618, 2, 5.116),
    ("F2",   17.423, 1, 1.602),
    ("HF",   15.271, 1, 5.869),
    ("HCl",  13.277, 1, 4.434),
    ("C-C",  11.260, 1, 3.600),
    ("C=C",  11.260, 2, 6.360),
    ("CC3",  11.260, 3, 8.700),
    ("CO",   12.328, 3, 11.09),
    ("C=O",  12.328, 2, 7.710),
    ("Cl2",  12.968, 1, 2.514),
    ("Li2",   5.392, 1, 1.046),
    ("NaCl",  7.361, 1, 4.230),
    ("LiH",   7.722, 1, 2.429),
]

print(f"{'mol':>6} bo  Eharm   cpl    Dpred  Dobs   err%    V8cpl  V8pred V8err%")
print("-"*75)

errs_new = []; errs_v8 = []
for name, Eh, bo, Dobs in mols:
    cpl = cpl_map[bo]
    Dpred = C_BOND * Eh * cpl
    v8c = 1.0 + (bo-1)*0.5
    v8p = C_BOND * Eh * v8c
    err = (Dpred-Dobs)/Dobs*100
    v8e = (v8p-Dobs)/Dobs*100
    errs_new.append(abs(err))
    errs_v8.append(abs(v8e))
    print(f"{name:>6} {bo}  {Eh:6.2f} {cpl:5.3f} {Dpred:7.3f} {Dobs:6.3f} {err:+6.1f}%  {v8c:5.3f} {v8p:7.3f} {v8e:+6.1f}%")

print()
print(f"HESSIAN: mean={np.mean(errs_new):.1f}%, median={np.median(errs_new):.1f}%, max={np.max(errs_new):.1f}%")
print(f"V8base: mean={np.mean(errs_v8):.1f}%, median={np.median(errs_v8):.1f}%, max={np.max(errs_v8):.1f}%")
print(f"  <10%: Hessian {sum(1 for e in errs_new if e<10)}/15, V8 {sum(1 for e in errs_v8 if e<10)}/15")
print(f"  <20%: Hessian {sum(1 for e in errs_new if e<20)}/15, V8 {sum(1 for e in errs_v8 if e<20)}/15")
