#!/usr/bin/env python3
"""
Collective Field Solver — v4 (3D damped-EOM time dynamics)
============================================================

The bond_3d_gpu.py archive script already does this; v4 is a focused
extension that runs H2/Cl2/F2 in parallel and properly tests the LP
emergence hypothesis.

Why this might work where SCF failed
-------------------------------------
Static mean-field SCF (v1-v3) failed because the cumulant expansion
<cos(pi*delta_phi)> ~= 1 - pi^2 <delta_phi^2>/2 is invalid when
<delta_phi^2> ~ O(1), which it always is for occupied bound states.

Time dynamics sidesteps the cumulant entirely: evolve the FULL nonlinear
EOM, no expansion, no approximation. cos(pi*phi) is computed exactly at
each timestep. The breather modes (which sit at non-trivial amplitude)
are handled correctly because they're part of the field, not perturbations
on top of it.

The damping term slowly drains kinetic energy until the field settles
into a stationary attractor — which is exactly the consolidation note's
framing of what a "stable configuration" actually is in real systems.

Why bond_3d_gpu.py's first H2 result is encouraging
----------------------------------------------------
First run with H2 vs F2 in atom_modes table:
  H atom: 1 mode (l=0, A_half)
  F atom: 3 modes (l=1, m=0 sigma; m=+1 LP full; m=-1 LP full)

V_int ratio F2/H2 = 0.13/0.34 = 0.39
Observed D_e ratio F2/H2 = 1.60/4.75 = 0.34

Match to ~15% with ZERO LP correction — the LP physics emerges from the
nonlinear field dynamics. If we calibrate H2 V_int -> 4.75 eV, F2 lands
at ~1.78 eV vs observed 1.60 eV (+11%).

This v4 extends to Cl2 and does cleaner R scans to pin down R_eq and
extract D_e systematically.
"""
import sys
import time
import numpy as np
from math import factorial

try:
    import cupy as cp
    GPU = True
    xp = cp
except ImportError:
    cp = np
    xp = np
    GPU = False

# ============================================================
# CONSTANTS
# ============================================================
PI = np.pi
d = 3
V_0 = 1.0 / PI**2
gamma_sg = PI / (2**(d+1) * PI - 2)
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/PI**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6  # 13.604 eV

omega = [float(np.cos(n * gamma_sg)) for n in range(1, 5)]
A_HALF = 0.10
A_FULL = 0.20

# ============================================================
# 3D GRID
# ============================================================
N = 48          # 48^3 = 110,592 (smaller than bond_3d_gpu's 64^3, faster)
BOX = 8.0       # box ±BOX Bohr
dx = 2 * BOX / N
dt = 0.25 * dx

print(f"GPU: {'CuPy + CUDA' if GPU else 'NumPy (CPU)'}")
print(f"Grid: {N}^3 = {N**3:,}, box=±{BOX} Bohr, dx={dx:.3f}")
print()

x1d = np.linspace(-BOX, BOX, N, endpoint=False)

def make_grid(center_z=0.0):
    X, Y, Z = np.meshgrid(x1d, x1d, x1d - center_z, indexing='ij')
    R = np.sqrt(X**2 + Y**2 + Z**2) + 1e-10
    return X, Y, Z, R

# ============================================================
# FIELD OPS
# ============================================================
def laplacian_3d(phi):
    return (
        xp.roll(phi, 1, 0) + xp.roll(phi, -1, 0) +
        xp.roll(phi, 1, 1) + xp.roll(phi, -1, 1) +
        xp.roll(phi, 1, 2) + xp.roll(phi, -1, 2) -
        6 * phi
    ) / dx**2

def total_energy(phi):
    GE = 0.5 * (
        (xp.roll(phi, 1, 0) - phi)**2 +
        (xp.roll(phi, 1, 1) - phi)**2 +
        (xp.roll(phi, 1, 2) - phi)**2
    )
    PE = V_0 * (1.0 - xp.cos(PI * phi))
    return float(xp.sum(GE + PE) * dx**3)

def evolve(phi_init, n_steps=3000, damping=0.02):
    phi = xp.asarray(phi_init.copy(), dtype=np.float64)
    phi_old = phi.copy()
    energies = []
    for step in range(n_steps):
        lap = laplacian_3d(phi)
        force = (1.0 / PI) * xp.sin(PI * phi)
        phi_new = (2 - damping*dt) * phi - (1 - damping*dt) * phi_old + dt**2 * (lap - force)
        phi_old = phi.copy()
        phi = phi_new
        if step > n_steps * 2 // 3 and step % 50 == 0:
            energies.append(total_energy(phi))
    return phi, np.mean(energies) if energies else total_energy(phi)

# ============================================================
# ATOM BUILDERS (from bond_3d_gpu.py)
# ============================================================
def kink_field(R_arr, Z_eff=1.0):
    R_g = xp.asarray(R_arr)
    return (4.0 / PI) * xp.arctan(xp.exp(-float(np.sqrt(Z_eff)) * R_g))

def breather_mode(X, Y, Z, R, l, m, amp, omega_n):
    X_g, Y_g, Z_g, R_g = xp.asarray(X), xp.asarray(Y), xp.asarray(Z), xp.asarray(R)
    eps = float(np.sqrt(max(1.0 - omega_n**2, 1e-12)))
    if l == 0:
        ang = xp.ones_like(R_g)
    elif l == 1:
        if m == 0:   ang = Z_g / R_g
        elif m == 1: ang = X_g / R_g
        elif m == -1: ang = Y_g / R_g
        else: ang = xp.ones_like(R_g)
    else:
        ang = xp.ones_like(R_g)
    radial = amp / (omega_n * xp.cosh(eps * R_g) + 1e-10)
    return (4.0 / PI) * xp.arctan(eps * ang * radial)

ATOM_MODES = {
    'H':  [(0, 0, A_HALF, omega[0])],
    'F':  [(1, 0, A_HALF, omega[1]),
           (1, 1, A_FULL, omega[1]),
           (1, -1, A_FULL, omega[1])],
    'Cl': [(1, 0, A_HALF, omega[1]),
           (1, 1, A_FULL, omega[1]),
           (1, -1, A_FULL, omega[1])],
    'O':  [(1, 0, A_FULL, omega[1]),
           (1, 1, A_HALF, omega[1]),
           (1, -1, A_HALF, omega[1])],
    'N':  [(1, 0, A_HALF, omega[1]),
           (1, 1, A_HALF, omega[1]),
           (1, -1, A_HALF, omega[1])],
    'C':  [(1, 0, A_HALF, omega[1]),
           (1, 1, A_HALF, omega[1])],
}

# Z_eff for sqrt(Z) kink scaling
Z_EFF = {'H': 1.0, 'C': 3.14, 'N': 3.83, 'O': 4.45,
         'F': 5.13, 'Cl': 6.12}

OBSERVED_DE = {
    ('H', 'H'): (4.478, 1.401),  # (D_e in eV, R_eq in Bohr)
    ('F', 'F'): (1.602, 1.412),
    ('Cl', 'Cl'): (2.514, 1.988),
    ('O', 'O'): (5.116, 1.207),
    ('N', 'N'): (9.759, 1.098),
}

# ============================================================
# BOND SIMULATION
# ============================================================
def build_atom(sym, center_z, Z_eff, sign_flips=None):
    """Build atom field. sign_flips: optional dict {(l,m): sign} to flip
    the sign of specific orbital modes (for MO construction)."""
    X, Y, Z, R = make_grid(center_z)
    phi = kink_field(R, Z_eff)
    sign_flips = sign_flips or {}
    for (l, m, amp, om) in ATOM_MODES.get(sym, []):
        sign = sign_flips.get((l, m), +1.0)
        phi = phi + sign * breather_mode(X, Y, Z, R, l, m, amp, om)
    return phi


def bond_mode_signs(sym_a, sym_b):
    """For each atom (A, B), return the sign flips needed to build proper
    bonding MOs from the orbital superposition.

    Standard chemistry MO rules:
      sigma (p_z along bond axis): bonding = p_z(A) - p_z(B). Flip B's p_z.
      pi (p_x, p_y perpendicular to bond axis): bonding = p_x(A) + p_x(B),
        same for y. No flip needed.
      s (l=0): no flip (always bonding when added).

    We assume the bond axis is z and the atoms are bonded with as many
    BOND pairs as their unpaired electrons allow. LP electrons (full pairs
    on one atom that don't participate) use naive superposition (no flip).
    """
    flips_A = {}
    flips_B = {}

    modes_a = ATOM_MODES.get(sym_a, [])
    modes_b = ATOM_MODES.get(sym_b, [])

    # Count unpaired electrons by orbital between the two atoms.
    # If both atoms have an unpaired electron in the same (l,m), they pair
    # into a bonding MO. If one atom has it doubly occupied (LP), it stays
    # localized.
    occupied_a = {(l, m): amp for (l, m, amp, om) in modes_a}
    occupied_b = {(l, m): amp for (l, m, amp, om) in modes_b}

    for (l, m), amp_a in occupied_a.items():
        amp_b = occupied_b.get((l, m), 0)
        # Both have it AND at least one is half-occupied (bond-forming)?
        if amp_b > 0 and (amp_a <= A_HALF + 1e-6 or amp_b <= A_HALF + 1e-6):
            # Bond pair. Apply MO sign convention.
            if l == 1 and m == 0:
                # Sigma from p_z: flip B
                flips_B[(l, m)] = -1.0
        # Otherwise: LP or non-bonding, no flip needed
    return flips_A, flips_B


def simulate_bond(sym_a, sym_b, R_bond, n_steps=2500, use_MO=True):
    """Returns (V_int_raw, time_seconds)."""
    Z_a = Z_EFF.get(sym_a, 1.0)
    Z_b = Z_EFF.get(sym_b, 1.0)
    if use_MO:
        flips_A, flips_B = bond_mode_signs(sym_a, sym_b)
    else:
        flips_A, flips_B = {}, {}
    t0 = time.time()
    # Isolated atoms (no flips — sign convention only matters for the dimer)
    phi_a = build_atom(sym_a, -R_bond/2, Z_a)
    _, E_a = evolve(phi_a, n_steps=n_steps // 2, damping=0.02)
    phi_b = build_atom(sym_b, +R_bond/2, Z_b)
    _, E_b = evolve(phi_b, n_steps=n_steps // 2, damping=0.02)
    # Combined with MO signs applied
    phi_ab = (build_atom(sym_a, -R_bond/2, Z_a, sign_flips=flips_A) +
              build_atom(sym_b, +R_bond/2, Z_b, sign_flips=flips_B))
    _, E_ab = evolve(phi_ab, n_steps=n_steps, damping=0.02)
    return E_ab - E_a - E_b, time.time() - t0

# ============================================================
# RUN: scan H2, Cl2, F2, and (for size variety) O2 and N2
# ============================================================
print("=" * 78)
print("COLLECTIVE FIELD SOLVER v4 — 3D damped-EOM time dynamics")
print("=" * 78)
print()

R_scan = [0.8, 1.0, 1.1, 1.2, 1.4, 1.8, 2.5, 3.5]   # extended down for N2/O2

molecules = [('H', 'H'), ('F', 'F'), ('Cl', 'Cl'), ('O', 'O'), ('N', 'N')]
results = {m: {} for m in molecules}

for mol in molecules:
    sym_a, sym_b = mol
    print(f"---- {sym_a}{sym_b} scan ----")
    for R in R_scan:
        V_int, elapsed = simulate_bond(sym_a, sym_b, R, n_steps=2000)
        results[mol][R] = V_int
        print(f"  R={R:.2f} Bohr  V_int={V_int:+.5f}  ({elapsed:.1f}s)")
    print()

# Pick most-negative V_int for each molecule
print("=" * 78)
print("SUMMARY (most-negative V_int across scan)")
print("=" * 78)
print(f"{'Mol':>6} {'V_int_min':>10} {'R_min':>7}   {'obs D_e':>10} {'obs R_eq':>10}")
V_int_H2 = None
for mol in molecules:
    V_min = min(results[mol].values())
    R_min = min(results[mol], key=lambda r: results[mol][r])
    obs_De, obs_R = OBSERVED_DE.get(mol, (None, None))
    if mol == ('H', 'H'):
        V_int_H2 = V_min
    obs_str = f"{obs_De:.3f} eV" if obs_De else "n/a"
    obsR_str = f"{obs_R:.3f} Bohr" if obs_R else "n/a"
    print(f"  {mol[0]}{mol[1]:>2} {V_min:>+10.5f} {R_min:>5.2f}   "
          f"{obs_str:>10} {obsR_str:>10}")
print()

# Calibrate to H2 and report
print("=" * 78)
print("CALIBRATED (V_int -> eV via H2 = 4.478 eV observed)")
print("=" * 78)
if V_int_H2 is not None:
    conv = 4.478 / abs(V_int_H2)
    print(f"  Conversion: D_e(eV) = |V_int| * {conv:.4f}")
    print(f"{'Mol':>6} {'V_int_min':>10} {'D_e calc':>10} {'D_e obs':>10} {'err':>8}")
    for mol in molecules:
        V_min = min(results[mol].values())
        De_calc = abs(V_min) * conv
        obs_De, _ = OBSERVED_DE.get(mol, (None, None))
        if obs_De:
            err = (De_calc - obs_De) / obs_De * 100
            print(f"  {mol[0]}{mol[1]:>2} {V_min:>+10.5f} "
                  f"{De_calc:>10.3f} {obs_De:>10.3f} {err:>+7.1f}%")
print()
print("Pass: errors comparable to V8's 7.8% mean. The damped-EOM approach")
print("would then BE the from-scratch derivation of V8's correction stack.")
