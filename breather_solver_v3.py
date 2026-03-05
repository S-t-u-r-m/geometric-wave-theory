"""
GWT 3D Breather Solver — v3: Targeted Kink Analysis
=====================================================
Key finding from v2: Width=3 kink-antikink is most stable in 1D.
Now test 3D configurations systematically.

Questions to answer:
1. Does a 3D kink (displaced cube) stay localized?
2. What is the minimum stable excitation in 3D?
3. Does d=3 appear as a magic number?
"""

import numpy as np
from scipy.integrate import solve_ivp
import time as timer

PI = np.pi
K_NN = 1.0
K_2NN = 0.5
M_KINK = 8.0 / PI**2

# ============================================================
# LATTICE FUNCTIONS (from v2)
# ============================================================
def forces_1d(u, k_nn):
    N = len(u)
    F = np.zeros(N)
    du_right = np.diff(u)
    f_right = (k_nn / PI) * np.sin(PI * du_right)
    F[:-1] += f_right
    F[1:] -= f_right
    return F

def energy_1d(u, p, k_nn):
    du = np.diff(u)
    return 0.5 * np.sum(p**2) + np.sum((k_nn/PI**2) * (1 - np.cos(PI*du)))

def eom_1d(t, y, N, k_nn):
    u, p = y[:N], y[N:]
    return np.concatenate([p, forces_1d(u, k_nn)])

def forces_3d(u, shape, k_nn, k_2nn):
    Nx, Ny, Nz = shape
    F = np.zeros_like(u)
    for axis in range(3):
        for sign in [1, -1]:
            u_shifted = np.roll(u, -sign, axis=axis)
            delta_u = u_shifted - u
            dhat = np.zeros(3)
            dhat[axis] = 1.0
            proj = np.einsum('...i,i->...', delta_u, dhat)
            f_mag = (k_nn / PI) * np.sin(PI * proj)
            for comp in range(3):
                F[..., comp] += f_mag * dhat[comp]
            if sign == 1:
                if axis == 0: F[-1, :, :, :] = 0
                elif axis == 1: F[:, -1, :, :] = 0
                else: F[:, :, -1, :] = 0
            else:
                if axis == 0: F[0, :, :, :] = 0
                elif axis == 1: F[:, 0, :, :] = 0
                else: F[:, :, 0, :] = 0
    for a1 in range(3):
        for a2 in range(a1+1, 3):
            for s1 in [1, -1]:
                for s2 in [1, -1]:
                    u_shifted = np.roll(np.roll(u, -s1, axis=a1), -s2, axis=a2)
                    delta_u = u_shifted - u
                    dhat = np.zeros(3)
                    dhat[a1] = s1 / np.sqrt(2)
                    dhat[a2] = s2 / np.sqrt(2)
                    proj = np.einsum('...i,i->...', delta_u, dhat)
                    f_mag = (k_2nn / (PI * np.sqrt(2))) * np.sin(PI * proj / np.sqrt(2))
                    for comp in range(3):
                        F[..., comp] += f_mag * dhat[comp]
    return F

def energy_3d(u, p, shape, k_nn, k_2nn):
    Nx, Ny, Nz = shape
    E = 0.5 * np.sum(p**2)
    for axis in range(3):
        u_shifted = np.roll(u, -1, axis=axis)
        delta_u = u_shifted - u
        dhat = np.zeros(3)
        dhat[axis] = 1.0
        proj = np.einsum('...i,i->...', delta_u, dhat)
        mask = np.ones((Nx, Ny, Nz), dtype=bool)
        if axis == 0: mask[-1, :, :] = False
        elif axis == 1: mask[:, -1, :] = False
        else: mask[:, :, -1] = False
        E += np.sum((k_nn / PI**2) * (1 - np.cos(PI * proj)) * mask)
    for a1 in range(3):
        for a2 in range(a1+1, 3):
            for s2 in [-1, 1]:
                u_shifted = np.roll(np.roll(u, -1, axis=a1), -s2, axis=a2)
                delta_u = u_shifted - u
                dhat = np.zeros(3)
                dhat[a1] = 1.0 / np.sqrt(2)
                dhat[a2] = s2 / np.sqrt(2)
                proj = np.einsum('...i,i->...', delta_u, dhat)
                mask = np.ones((Nx, Ny, Nz), dtype=bool)
                if a1 == 0: mask[-1, :, :] = False
                elif a1 == 1: mask[:, -1, :] = False
                else: mask[:, :, -1] = False
                if s2 == 1:
                    if a2 == 0: mask[-1, :, :] = False
                    elif a2 == 1: mask[:, -1, :] = False
                    else: mask[:, :, -1] = False
                else:
                    if a2 == 0: mask[0, :, :] = False
                    elif a2 == 1: mask[:, 0, :] = False
                    else: mask[:, :, 0] = False
                E += np.sum((k_2nn / PI**2) * (1 - np.cos(PI * proj / np.sqrt(2))) * mask)
    return E

def eom_3d_flat(t, y_flat, shape, k_nn, k_2nn):
    N_total = np.prod(shape) * 3
    u = y_flat[:N_total].reshape(*shape, 3)
    p = y_flat[N_total:].reshape(*shape, 3)
    F = forces_3d(u, shape, k_nn, k_2nn)
    return np.concatenate([p.ravel(), F.ravel()])

def local_energy_3d(u, p, shape, k_nn, center, hw):
    """Energy in a cube of half-width hw around center."""
    c = center
    Nx, Ny, Nz = shape
    E_loc = 0.0
    for ix in range(max(0, c[0]-hw), min(Nx, c[0]+hw+1)):
        for iy in range(max(0, c[1]-hw), min(Ny, c[1]+hw+1)):
            for iz in range(max(0, c[2]-hw), min(Nz, c[2]+hw+1)):
                E_loc += 0.5 * np.sum(p[ix,iy,iz,:]**2)
                for axis in range(3):
                    jx, jy, jz = ix, iy, iz
                    if axis == 0: jx += 1
                    elif axis == 1: jy += 1
                    else: jz += 1
                    if jx < Nx and jy < Ny and jz < Nz:
                        du = u[jx,jy,jz,axis] - u[ix,iy,iz,axis]
                        E_loc += (k_nn/PI**2) * (1 - np.cos(PI*du))
    return E_loc

# ============================================================
# TEST 1: Width=3 in 1D — Why is it special?
# ============================================================
print("=" * 70)
print("TEST 1: Why Width=3 is Special in 1D")
print("=" * 70)
print()

N = 128
center = N // 2

# Width=3 kink-antikink: u = [0,...,0, 1, 1, 1, 0,...,0]
# This has TWO kinks (0->1 and 1->0), each with energy 2/pi^2
# But width=3 means the TOPOLOGICAL charge spans d=3 sites

print("Energy vs width for 1D kink-antikink pairs:")
print(f"{'Width':>6} {'E_total':>10} {'E/M_kink':>10} {'Kinks':>6} {'E_per_kink':>12}")
print("-" * 50)
for w in range(1, 8):
    u = np.zeros(N)
    for i in range(w):
        u[center + i] = 1.0
    p = np.zeros(N)
    E = energy_1d(u, p, K_NN)
    print(f"{w:6d} {E:10.6f} {E/M_KINK:10.4f} {2:6d} {E/2:12.6f}")

print()
print("NOTE: Energy is ALWAYS 2/pi^2 * 2 = 4/pi^2 = 0.4053 regardless of width!")
print("(Each kink boundary contributes exactly 2/pi^2 = M_kink/4)")
print()

# What about AMPLITUDE variation?
print("Energy vs displacement amplitude (width=3):")
print(f"{'u_max':>6} {'E_total':>10} {'E/M_kink':>10}")
print("-" * 30)
for u_max in [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]:
    u = np.zeros(N)
    for i in range(3):
        u[center + i] = u_max
    p = np.zeros(N)
    E = energy_1d(u, p, K_NN)
    print(f"{u_max:6.2f} {E:10.6f} {E/M_KINK:10.4f}")

print()

# ============================================================
# TEST 2: 1D Stability vs Width (long time)
# ============================================================
print("=" * 70)
print("TEST 2: 1D Localization vs Width (t=200)")
print("=" * 70)
print()

print(f"{'Width':>6} {'E_init':>10} {'frac(t=50)':>12} {'frac(t=100)':>12} {'frac(t=200)':>12}")
print("-" * 60)
for w in range(1, 7):
    u_init = np.zeros(N)
    for i in range(w):
        u_init[center + i] = 1.0
    p_init = np.zeros(N)
    E0 = energy_1d(u_init, p_init, K_NN)

    y0 = np.concatenate([u_init, p_init])
    sol = solve_ivp(eom_1d, [0, 200], y0, args=(N, K_NN),
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    t_eval=[50, 100, 200])

    fracs = []
    for idx_t in range(sol.y.shape[1]):
        us = sol.y[:N, idx_t]
        ps = sol.y[N:, idx_t]
        Es = energy_1d(us, ps, K_NN)

        hw = w + 4
        E_loc = 0.5 * np.sum(ps[center-hw:center+w+hw]**2)
        for i in range(center-hw, center+w+hw):
            if 0 <= i < N-1:
                du = us[i+1] - us[i]
                E_loc += (K_NN/PI**2) * (1 - np.cos(PI*du))
        fracs.append(E_loc / Es)

    print(f"{w:6d} {E0:10.6f} {fracs[0]:12.3f} {fracs[1]:12.3f} {fracs[2]:12.3f}")

print()

# ============================================================
# TEST 3: 3D Cubic Kink (displaced cube)
# ============================================================
print("=" * 70)
print("TEST 3: 3D Cubic Kink Configurations")
print("=" * 70)
print()

N3 = 16
shape3 = (N3, N3, N3)
c3 = N3 // 2

# A "3D kink" = displacing a cube of sites
# If we displace a w x w x w cube by u=1 along x:
# - Creates 2*w^2 kink faces (front and back in x-direction)
# - Each face has w^2 bonds, each contributing 2/pi^2 energy

print("3D Cubic kink: displace w x w x w cube by u_x = 1")
print(f"{'w':>4} {'E_total':>10} {'E/(w^2*M_kink)':>16} {'E/M_kink':>10}")
print("-" * 45)

for w in [1, 2, 3, 4, 5]:
    u = np.zeros((N3, N3, N3, 3))
    for ix in range(w):
        for iy in range(w):
            for iz in range(w):
                u[c3+ix, c3+iy, c3+iz, 0] = 1.0
    p = np.zeros_like(u)
    E = energy_3d(u, p, shape3, K_NN, K_2NN)
    print(f"{w:4d} {E:10.4f} {E/(w**2 * M_KINK):16.4f} {E/M_KINK:10.4f}")

print()

# ============================================================
# TEST 4: 3D Kink Stability (w=1,2,3 cubes)
# ============================================================
print("=" * 70)
print("TEST 4: 3D Kink Stability (time evolution)")
print("=" * 70)
print()

# Use smaller lattice for speed
N3s = 12
shape3s = (N3s, N3s, N3s)
c3s = N3s // 2

for w in [1, 2, 3]:
    u = np.zeros((N3s, N3s, N3s, 3))
    for ix in range(w):
        for iy in range(w):
            for iz in range(w):
                u[c3s+ix, c3s+iy, c3s+iz, 0] = 1.0
    p = np.zeros_like(u)
    E0 = energy_3d(u, p, shape3s, K_NN, K_2NN)

    y0 = np.concatenate([u.ravel(), p.ravel()])
    t_checks = [5, 10, 20, 30]
    t0 = timer.time()
    sol = solve_ivp(eom_3d_flat, [0, 30], y0, args=(shape3s, K_NN, K_2NN),
                    method='DOP853', rtol=1e-9, atol=1e-11,
                    max_step=0.3, t_eval=t_checks)
    elapsed = timer.time() - t0

    print(f"3D Cube w={w}: E0={E0:.4f} ({elapsed:.1f}s)")
    for idx_t in range(sol.y.shape[1]):
        ys = sol.y[:, idx_t]
        us = ys[:N3s**3*3].reshape(N3s, N3s, N3s, 3)
        ps = ys[N3s**3*3:].reshape(N3s, N3s, N3s, 3)
        Es = energy_3d(us, ps, shape3s, K_NN, K_2NN)
        E_loc = local_energy_3d(us, ps, shape3s, K_NN, (c3s, c3s, c3s), w+2)
        frac = E_loc / Es if Es > 1e-15 else 0
        print(f"  t={t_checks[idx_t]:3d}: E={Es:.4f}, loc={frac:.3f}, max|u|={np.max(np.abs(us)):.4f}")
    print()

# ============================================================
# TEST 5: The d=3 Connection
# ============================================================
print("=" * 70)
print("TEST 5: The d=3 Connection")
print("=" * 70)
print()

# Key ratios
print("Critical ratios involving d=3:")
print(f"  sqrt(2d)/2sqrt(d) = sqrt(6)/(2*sqrt(3)) = 1/sqrt(2) = {1/np.sqrt(2):.6f}")
print(f"  This is the Koide amplitude |b| = 1/sqrt(2) = {1/np.sqrt(2):.6f}")
print()
print(f"  Effective well: omega_eff = sqrt(2d) = sqrt(6) = {np.sqrt(6):.6f}")
print(f"  Phonon band top: 2*sqrt(d) = 2*sqrt(3) = {2*np.sqrt(3):.6f}")
print(f"  Ratio: {np.sqrt(6)/(2*np.sqrt(3)):.6f} = 1/sqrt(2)")
print()

# The number of breather states is floor(pi/(2*gamma)) = floor((16*pi - 2)/2) = 24
# But in 3D with d=3, we need to count ORIENTED states
# Each 1D breather can propagate in 3 directions (x,y,z)
# and have 2 polarizations (longitudinal/transverse per axis)
# For cubic symmetry: representations of S_3

n_1d = 24  # 1D breather count
print(f"  1D breather states: {n_1d}")
print(f"  3D oriented states: {n_1d} x d = {n_1d * 3} = 72")
print(f"  But cubic symmetry S_3 identifies orientations...")
print(f"  Distinct representations: {n_1d} x {3} / 3 = {n_1d} (if fully symmetric)")
print(f"  OR with polarization: {n_1d} x (d-1) transverse = {n_1d * 2}")
print()

# Width=3 kink-antikink mass in various dimensions
print("Kink-antikink (width=d) energy in d dimensions:")
print("  (using 1D chain, width = d sites)")
for d in [1, 2, 3, 4, 5, 6]:
    u = np.zeros(N)
    for i in range(d):
        u[center + i] = 1.0
    E = energy_1d(u, np.zeros(N), K_NN)
    print(f"  d={d}: E = {E:.6f} = {E/M_KINK:.4f} * M_kink")

print()

# ============================================================
# TEST 6: Normal mode analysis of kink-antikink
# ============================================================
print("=" * 70)
print("TEST 6: Kink-Antikink Internal Mode Frequencies")
print("=" * 70)
print()

# For a width-w kink-antikink, linearize around the static config
# and find the internal oscillation frequencies
# H_ij = d^2V/du_i du_j

for w in [1, 2, 3, 4, 5]:
    u_eq = np.zeros(N)
    for i in range(w):
        u_eq[center + i] = 1.0

    # Hessian matrix (only compute for sites near the kink)
    hw = w + 6
    sites = list(range(center - hw, center + w + hw))
    n_sites = len(sites)
    H = np.zeros((n_sites, n_sites))

    for a, i in enumerate(sites):
        for b, j in enumerate(sites):
            if abs(i - j) > 1:
                continue
            if i == j:
                # Diagonal: sum of bond stiffnesses
                if i > 0:
                    du = u_eq[i] - u_eq[i-1]
                    H[a, a] += (K_NN) * np.cos(PI * du)
                if i < N-1:
                    du = u_eq[i+1] - u_eq[i]
                    H[a, a] += (K_NN) * np.cos(PI * du)
            elif j == i + 1:
                du = u_eq[j] - u_eq[i]
                H[a, b] = -(K_NN) * np.cos(PI * du)
            elif j == i - 1:
                du = u_eq[i] - u_eq[j]
                H[a, b] = -(K_NN) * np.cos(PI * du)

    evals = np.linalg.eigvalsh(H)
    # Filter out near-zero (translation mode) and negative eigenvalues
    pos_evals = evals[evals > 0.01]
    freqs = np.sqrt(pos_evals)

    # Internal modes are those below the phonon band edge (omega = 2 in 1D)
    internal = freqs[freqs < 2.0]
    band_edge = 2.0 * np.sqrt(K_NN)

    print(f"Width={w}: {len(internal)} internal modes below band edge ({band_edge:.3f})")
    if len(internal) > 0:
        print(f"  Frequencies: {', '.join(f'{f:.4f}' for f in sorted(internal))}")
        print(f"  Lowest freq: {min(internal):.6f}")

    # Count zero modes
    zero_modes = np.sum(np.abs(evals) < 0.01)
    print(f"  Zero modes (translation): {zero_modes}")
    neg_modes = np.sum(evals < -0.01)
    if neg_modes > 0:
        print(f"  WARNING: {neg_modes} NEGATIVE eigenvalues (unstable!)")
        neg_evals = evals[evals < -0.01]
        print(f"  Negative eigenvalues: {', '.join(f'{e:.4f}' for e in neg_evals)}")
    print()

print()
print("=" * 70)
print("FINAL SUMMARY")
print("=" * 70)
print()
print("The 1D kink-antikink of width w has CONSTANT energy = 4/pi^2 = M_kink/2")
print("regardless of width (each boundary contributes 2/pi^2).")
print()
print("In 3D, a cubic kink of width w has energy proportional to w^2")
print("(surface area of the displaced cube).")
print()
print("The ratio omega_eff/omega_band = sqrt(6)/(2*sqrt(3)) = 1/sqrt(2)")
print("is EXACTLY the Koide amplitude |b| = sqrt(kappa/k).")
print()
print("d=3 appears as:")
print("  - Optimal kink width for 1D stability")
print("  - Koide phase delta = 2/d^2 = 2/9")
print("  - Number of spatial dimensions")
print("  - Cubic symmetry S_3")
print("  - Q = (d-1)/d = 2/3")
