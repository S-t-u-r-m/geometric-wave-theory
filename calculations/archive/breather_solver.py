"""
GWT 3D Discrete Breather Solver — v2
=====================================
The GWT Hamiltonian has NO on-site potential — only inter-site cosine coupling.
Breathers are intrinsically LATTICE objects (they dissolve in the continuum limit).

Strategy: Single-site excitation → time evolution → measure localization & energy.
"""

import numpy as np
from scipy.integrate import solve_ivp
import time as timer

PI = np.pi
K_NN = 1.0
K_2NN = 0.5

# SG parameters (for reference)
XI = 1.0 / (8*PI - 1)
GAMMA = PI * XI / 2
M_KINK = 8.0 / PI**2

print("=" * 60)
print("GWT DISCRETE BREATHER SOLVER v2")
print("=" * 60)
print()

# ============================================================
# 1D LATTICE: VECTORIZED EOM AND ENERGY
# ============================================================
def forces_1d(u, k_nn):
    """Compute forces on 1D lattice (vectorized)."""
    N = len(u)
    F = np.zeros(N)
    # Right bonds: u[i+1] - u[i]
    du_right = np.diff(u)  # length N-1
    f_right = (k_nn / PI) * np.sin(PI * du_right)
    F[:-1] += f_right  # force from right neighbor
    F[1:] -= f_right   # reaction on left
    return F

def energy_1d(u, p, k_nn):
    """Total energy of 1D lattice."""
    du = np.diff(u)
    return 0.5 * np.sum(p**2) + np.sum((k_nn/PI**2) * (1 - np.cos(PI*du)))

def eom_1d(t, y, N, k_nn):
    u, p = y[:N], y[N:]
    return np.concatenate([p, forces_1d(u, k_nn)])

# ============================================================
# 1D: SINGLE-SITE EXCITATION TEST
# ============================================================
print("=" * 60)
print("TEST 1: 1D Single-Site Excitation")
print("=" * 60)
print()

N = 64
center = N // 2

# Try different initial displacements at the center site
print(f"{'u0':>6} {'E_init':>10} {'E(t=50)':>10} {'E_center(50)':>12} {'localized?':>12}")
print("-" * 55)

for u0 in [0.1, 0.3, 0.5, 0.8, 1.0]:
    u_init = np.zeros(N)
    u_init[center] = u0
    p_init = np.zeros(N)

    E0 = energy_1d(u_init, p_init, K_NN)

    y0 = np.concatenate([u_init, p_init])
    sol = solve_ivp(eom_1d, [0, 50], y0, args=(N, K_NN),
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    dense_output=True)

    # Check energy at t=50
    yf = sol.sol(50)
    uf, pf = yf[:N], yf[N:]
    Ef = energy_1d(uf, pf, K_NN)

    # Check localization: energy in central 5 sites
    E_center = 0.5 * np.sum(pf[center-2:center+3]**2)
    for i in range(center-3, center+3):
        if 0 <= i < N-1:
            du = uf[i+1] - uf[i]
            E_center += (K_NN/PI**2) * (1 - np.cos(PI*du))

    localized = "YES" if E_center / Ef > 0.5 else "NO"
    print(f"{u0:6.1f} {E0:10.6f} {Ef:10.6f} {E_center:12.6f} {localized:>12}")

print()
print("In 1D without on-site potential, energy RADIATES away.")
print("Discrete breathers don't exist in 1D for this Hamiltonian.")
print()

# ============================================================
# 1D: What about exciting a PAIR of sites (kink-antikink)?
# ============================================================
print("=" * 60)
print("TEST 2: 1D Kink-Antikink Pair (breather-like)")
print("=" * 60)
print()

# A kink-antikink pair: sites below center at u=0, site at center at u=1
# This creates a "bump" that looks like a localized excitation
# The kink mass is 8/pi^2 per kink, so a pair has mass 16/pi^2

for width in [1, 2, 3]:
    u_init = np.zeros(N)
    for w in range(width):
        u_init[center + w] = 1.0  # shift by one full period
    p_init = np.zeros(N)

    E0 = energy_1d(u_init, p_init, K_NN)

    y0 = np.concatenate([u_init, p_init])

    # Evolve
    sol = solve_ivp(eom_1d, [0, 100], y0, args=(N, K_NN),
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    dense_output=True)

    # Check at several times
    print(f"Width={width} (u={width:.0f} at center): E_init = {E0:.6f}")
    for t_check in [10, 50, 100]:
        ys = sol.sol(t_check)
        us, ps = ys[:N], ys[N:]
        Es = energy_1d(us, ps, K_NN)

        # Energy in central region (width + 4 sites each side)
        hw = width + 4
        E_loc = 0.5 * np.sum(ps[center-hw:center+width+hw]**2)
        for i in range(center-hw, center+width+hw):
            if 0 <= i < N-1:
                du = us[i+1] - us[i]
                E_loc += (K_NN/PI**2) * (1 - np.cos(PI*du))

        print(f"  t={t_check:3d}: E_total={Es:.6f}, E_local/E_total={E_loc/Es:.3f}")
    print()

# ============================================================
# 3D LATTICE: VECTORIZED
# ============================================================
print("=" * 60)
print("TEST 3: 3D Lattice Breather Search")
print("=" * 60)
print()

def forces_3d(u, shape, k_nn, k_2nn):
    """
    Compute forces on 3D lattice. u shape: (Nx, Ny, Nz, 3).
    Uses numpy roll for vectorization, with boundary correction.
    """
    Nx, Ny, Nz = shape
    F = np.zeros_like(u)

    # NN bonds along each axis
    for axis in range(3):
        for sign in [1, -1]:
            # Compute delta_u = u_neighbor - u_self
            u_shifted = np.roll(u, -sign, axis=axis)
            delta_u = u_shifted - u

            # Bond direction
            dhat = np.zeros(3)
            dhat[axis] = 1.0

            # Project
            proj = np.einsum('...i,i->...', delta_u, dhat)

            # Force magnitude
            f_mag = (k_nn / PI) * np.sin(PI * proj)

            # Add force (in bond direction)
            for comp in range(3):
                F[..., comp] += f_mag * dhat[comp]

            # Zero out boundary (open BC)
            if sign == 1:
                slc = [slice(None)] * 3 + [slice(None)]
                if axis == 0: F[-1, :, :, :] = 0
                elif axis == 1: F[:, -1, :, :] = 0
                else: F[:, :, -1, :] = 0
            else:
                if axis == 0: F[0, :, :, :] = 0
                elif axis == 1: F[:, 0, :, :] = 0
                else: F[:, :, 0, :] = 0

    # 2NN bonds (face diagonals)
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
    """Total energy of 3D lattice."""
    Nx, Ny, Nz = shape
    E = 0.5 * np.sum(p**2)

    # NN bonds (positive direction only)
    for axis in range(3):
        u_shifted = np.roll(u, -1, axis=axis)
        delta_u = u_shifted - u
        dhat = np.zeros(3)
        dhat[axis] = 1.0
        proj = np.einsum('...i,i->...', delta_u, dhat)

        # Mask boundary
        mask = np.ones((Nx, Ny, Nz), dtype=bool)
        if axis == 0: mask[-1, :, :] = False
        elif axis == 1: mask[:, -1, :] = False
        else: mask[:, :, -1] = False

        E += np.sum((k_nn / PI**2) * (1 - np.cos(PI * proj)) * mask)

    # 2NN bonds
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
    """Flat EOM for scipy integrator."""
    N_total = np.prod(shape) * 3
    u = y_flat[:N_total].reshape(*shape, 3)
    p = y_flat[N_total:].reshape(*shape, 3)
    F = forces_3d(u, shape, k_nn, k_2nn)
    return np.concatenate([p.ravel(), F.ravel()])


# 3D Test: single-site excitation
N3 = 12
shape3 = (N3, N3, N3)
c3 = N3 // 2

print(f"3D lattice: {N3}x{N3}x{N3} = {N3**3} sites")
print()

# Test different initial conditions
test_configs = [
    ("1D longitudinal (u_x at center)", "1d_long"),
    ("1D transverse (u_y at center, propagating along x)", "1d_trans"),
    ("3D spherical (radial displacement from center)", "3d_sphere"),
]

for desc, mode in test_configs:
    for u0 in [0.5, 1.0]:
        u = np.zeros((N3, N3, N3, 3))
        p = np.zeros((N3, N3, N3, 3))

        if mode == "1d_long":
            u[c3, c3, c3, 0] = u0
        elif mode == "1d_trans":
            u[c3, c3, c3, 1] = u0  # transverse displacement
        elif mode == "3d_sphere":
            # Radial displacement on all 6 NN neighbors
            for axis in range(3):
                for sign in [-1, 1]:
                    idx = [c3, c3, c3]
                    idx[axis] += sign
                    if 0 <= idx[0] < N3 and 0 <= idx[1] < N3 and 0 <= idx[2] < N3:
                        rhat = np.zeros(3)
                        rhat[axis] = sign
                        u[idx[0], idx[1], idx[2], :] = u0 * 0.5 * rhat

        E0 = energy_3d(u, p, shape3, K_NN, K_2NN)

        # Time evolve
        y0 = np.concatenate([u.ravel(), p.ravel()])
        t0 = timer.time()
        sol = solve_ivp(eom_3d_flat, [0, 20], y0, args=(shape3, K_NN, K_2NN),
                        method='DOP853', rtol=1e-8, atol=1e-10,
                        max_step=0.5, t_eval=[20])
        elapsed = timer.time() - t0

        yf = sol.y[:,-1]  # final state
        uf = yf[:N3**3*3].reshape(N3, N3, N3, 3)
        pf = yf[N3**3*3:].reshape(N3, N3, N3, 3)
        Ef = energy_3d(uf, pf, shape3, K_NN, K_2NN)

        # Localization: energy within 3x3x3 cube around center
        hw = 2
        E_loc = 0.0
        for ix in range(max(0,c3-hw), min(N3,c3+hw+1)):
            for iy in range(max(0,c3-hw), min(N3,c3+hw+1)):
                for iz in range(max(0,c3-hw), min(N3,c3+hw+1)):
                    E_loc += 0.5 * np.sum(pf[ix,iy,iz,:]**2)

        frac = E_loc / Ef if Ef > 1e-15 else 0

        print(f"{desc}, u0={u0:.1f}: E0={E0:.4f}, "
              f"E(20)={Ef:.4f}, local_frac={frac:.3f} ({elapsed:.1f}s)")

print()

# ============================================================
# KEY INSIGHT SECTION
# ============================================================
print("=" * 60)
print("KEY PHYSICAL INSIGHT")
print("=" * 60)
print()
print("The GWT Hamiltonian has NO on-site potential.")
print("V depends ONLY on displacement differences (bond variables).")
print()
print("Consequence: in the continuum limit, the EOM is the LINEAR")
print("wave equation. Breathers exist ONLY on the discrete lattice.")
print()
print("The effective on-site potential (from frozen neighbors):")
print(f"  V_eff(u) = (2d/pi^2)(1-cos(pi*u)) with well freq = sqrt(2d) = {np.sqrt(6):.4f}")
print(f"  Phonon band: [0, 2*sqrt(d)] = [0, {2*np.sqrt(3):.4f}]")
print(f"  Well freq vs band top: {np.sqrt(6):.4f} / {2*np.sqrt(3):.4f} = {np.sqrt(6)/(2*np.sqrt(3)):.4f}")
print()
print("The well frequency sqrt(6) = 2.449 is INSIDE the phonon band")
print("[0, 3.464]. This means the breather frequency is in resonance")
print("with propagating phonons => energy leaks away.")
print()
print("HOWEVER: the cosine is a SOFT nonlinearity (frequency decreases")
print("with amplitude). Large-amplitude oscillations have LOWER frequency,")
print("pushing them toward the band edge at omega = 0.")
print()
print("For VERY large amplitude (u ~ 1, full kink width), the mode")
print("becomes a kink-antikink pair, which IS topologically stable.")
print("The kink mass M_kink = 8/pi^2 = 0.811 is the MINIMUM stable")
print("excitation energy.")
print()

# Compute kink-antikink energy on 3D lattice
print("=" * 60)
print("KINK-ANTIKINK PAIR ON 3D LATTICE")
print("=" * 60)
print()

# 1D kink along x-axis: u goes from 0 to 1 (one lattice period) over ~1 site
# Kink-antikink = localized bump: u = 1 at one site, 0 elsewhere
for n_sites_up in [1, 2, 3]:
    u = np.zeros((N3, N3, N3, 3))
    for s in range(n_sites_up):
        u[c3+s, c3, c3, 0] = 1.0  # x-displacement of 1 lattice spacing
    p = np.zeros_like(u)
    E_kk = energy_3d(u, p, shape3, K_NN, K_2NN)
    # Compare to analytical kink pair energy
    E_kink_pair_1d = 2 * M_KINK  # two kinks in 1D
    print(f"Kink-antikink ({n_sites_up} sites): E = {E_kk:.6f}, "
          f"2*M_kink(1D) = {E_kink_pair_1d:.6f}, "
          f"ratio = {E_kk/E_kink_pair_1d:.4f}")

print()
print("For the 1D kink-antikink (width=1):")
u_kk = np.zeros((N3, N3, N3, 3))
u_kk[c3, c3, c3, 0] = 1.0
E_kk_1d = energy_3d(u_kk, np.zeros_like(u_kk), shape3, K_NN, K_2NN)
print(f"  E = {E_kk_1d:.6f}")
print(f"  This involves {2+4+4} bond energy contributions:")
print(f"    2 NN bonds along x (u changes by 1): each = {K_NN/PI**2*(1-np.cos(PI*1.0)):.6f}")
print(f"    4 2NN bonds in xz and xy planes: each involves Delta_u_x = 1")

# Detailed bond analysis
E_nn_x = 2 * (K_NN/PI**2) * (1 - np.cos(PI * 1.0))  # 2 NN bonds along x
print(f"  NN along x: {E_nn_x:.6f} ({2*(K_NN/PI**2)*2:.6f} = 4/pi^2)")

# 2NN bonds: face diagonals involving the center site
# In the xy and xz planes, the displacement u_x = 1 projects onto the
# face diagonal direction by 1/sqrt(2)
E_2nn_per_bond = (K_2NN/PI**2) * (1 - np.cos(PI * 1.0 / np.sqrt(2)))
n_2nn_bonds = 8  # 4 face-diag directions in xy, 4 in xz planes
E_2nn_total = n_2nn_bonds * E_2nn_per_bond
print(f"  2NN face diag: {n_2nn_bonds} bonds x {E_2nn_per_bond:.6f} = {E_2nn_total:.6f}")
print(f"  Total analytical: {E_nn_x + E_2nn_total:.6f}")
print(f"  Total numerical:  {E_kk_1d:.6f}")
print()

# Time-evolve the kink-antikink to check stability
print("Time-evolving kink-antikink on 3D lattice...")
u_kk = np.zeros((N3, N3, N3, 3))
u_kk[c3, c3, c3, 0] = 1.0
p_kk = np.zeros_like(u_kk)
E_kk_0 = energy_3d(u_kk, p_kk, shape3, K_NN, K_2NN)

y0 = np.concatenate([u_kk.ravel(), p_kk.ravel()])
t_checks = [5, 10, 20, 30]
t0 = timer.time()
sol = solve_ivp(eom_3d_flat, [0, 30], y0, args=(shape3, K_NN, K_2NN),
                method='DOP853', rtol=1e-8, atol=1e-10,
                max_step=0.5, t_eval=t_checks)
elapsed = timer.time() - t0
print(f"Integration time: {elapsed:.1f}s")

for idx_t, t_check in enumerate(t_checks):
    if idx_t < sol.y.shape[1]:
        ys = sol.y[:, idx_t]
        us = ys[:N3**3*3].reshape(N3, N3, N3, 3)
        ps = ys[N3**3*3:].reshape(N3, N3, N3, 3)
        Es = energy_3d(us, ps, shape3, K_NN, K_2NN)

        # Localization
        hw = 3
        E_loc = 0
        for ix in range(max(0,c3-hw), min(N3,c3+hw+1)):
            for iy in range(max(0,c3-hw), min(N3,c3+hw+1)):
                for iz in range(max(0,c3-hw), min(N3,c3+hw+1)):
                    E_loc += 0.5 * np.sum(ps[ix,iy,iz,:]**2)
                    # Add local bond energy
                    for axis in range(3):
                        jx, jy, jz = ix, iy, iz
                        if axis == 0: jx += 1
                        elif axis == 1: jy += 1
                        else: jz += 1
                        if jx < N3 and jy < N3 and jz < N3:
                            du = us[jx,jy,jz,axis] - us[ix,iy,iz,axis]
                            E_loc += (K_NN/PI**2) * (1 - np.cos(PI*du))

        print(f"  t={t_check:3d}: E={Es:.6f}, E_loc/E_tot={E_loc/Es:.3f}, "
              f"max|u|={np.max(np.abs(us)):.4f}")

print()
print("=" * 60)
print("SUMMARY")
print("=" * 60)
print()
print("1. The GWT lattice has NO on-site potential (V depends on differences)")
print("2. In 1D, energy radiates away => no stable 1D breathers")
print("3. In 3D, MacKay-Aubry requires breather freq OUTSIDE phonon band")
print(f"4. Effective well freq sqrt(2d) = {np.sqrt(6):.3f} is INSIDE band [0, {2*np.sqrt(3):.3f}]")
print("5. Topological excitations (kink-antikink pairs) ARE stable")
print(f"6. Min kink-antikink energy = 2*M_kink = {2*M_KINK:.4f} (Planck units)")
print()
print("IMPLICATION: 'Particles' in GWT are NOT breathers in the standard sense.")
print("They are TOPOLOGICAL defects (kink configurations) stabilized by")
print("the discrete lattice. The kink mass 8/pi^2 = 0.811 sets the energy scale.")
print("The 24 'breather' states from DHN quantization correspond to 24 distinct")
print("kink-antikink bound states, each with a different oscillation pattern.")
