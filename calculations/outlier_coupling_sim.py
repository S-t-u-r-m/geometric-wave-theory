"""
Outlier Coupling Simulation
============================
Targeted simulations for the combined formula's worst outliers:
  As (+24%): half-fill p with d10 core
  Bi (+34%): half-fill p with d10 + f14 core
  Tl (-24%): underfill p with d10 + f14 core
  Lu (+28%): d1 on f14 core

Measure: how do f-modes and d-modes couple to p-valence and s-valence
when BOTH are present? This is the missing multi-mode coupling.
"""

import numpy as np
from math import factorial
import time

d = 3
gamma_sg = np.pi / (2**(d+1) * np.pi - 2)
V_0 = 1.0 / np.pi**2

Nx = 1500
x_max = 25.0
dx = 2 * x_max / Nx
x = np.linspace(-x_max, x_max, Nx)
dt = 0.4 * dx
N_steps = 12000

_Z_well = 1.0

def set_Z(Z):
    global _Z_well
    _Z_well = Z

def kink(x_grid, Z):
    return (4.0 / np.pi) * np.arctan(np.exp(np.sqrt(Z) * x_grid))

def add_breather(phi_bg, x_grid, omega, x0=0.0):
    eps = np.sqrt(max(1.0 - omega**2, 1e-12))
    num = eps
    den = omega * np.cosh(eps * (x_grid - x0))
    return phi_bg + (4.0 / np.pi) * np.arctan(num / (den + 1e-30))

def evolve(phi_init):
    phi = phi_init.copy()
    phi_old = phi.copy()
    energies = []
    for step in range(N_steps):
        lap = np.zeros_like(phi)
        lap[1:-1] = (phi[2:] - 2*phi[1:-1] + phi[:-2]) / dx**2
        force = _Z_well * (1.0 / np.pi) * np.sin(np.pi * phi)
        phi_new = 2*phi - phi_old + dt**2 * (lap - force)
        phi_new[0] = phi_init[0]
        phi_new[-1] = phi_init[-1]
        phi_old = phi.copy()
        phi = phi_new
        if step > N_steps // 2 and step % 50 == 0:
            KE = 0.5 * np.sum(((phi - phi_old)/dt)**2) * dx
            GE = 0.5 * np.sum(((phi[1:] - phi[:-1])/dx)**2) * dx
            PE = np.sum(_Z_well * V_0 * (1 - np.cos(np.pi * phi))) * dx
            energies.append(KE + GE + PE)
    return np.mean(energies) if energies else 0.0


omega_s = np.cos(1 * gamma_sg)  # s-mode (l=0)
omega_p = np.cos(2 * gamma_sg)  # p-mode (l=1)
omega_d = np.cos(3 * gamma_sg)  # d-mode (l=2)
omega_f = np.cos(4 * gamma_sg)  # f-mode (l=3)

Z = 10.0  # Use Z=10 for clear signal
set_Z(Z)

print("Outlier Coupling Simulation")
print("=" * 60)
print(f"Z={Z}, modes: s={omega_s:.4f} p={omega_p:.4f} d={omega_d:.4f} f={omega_f:.4f}")
print()

phi_k = kink(x, Z)
E_k = evolve(phi_k)

# Bare single-mode bindings
E_s = evolve(add_breather(phi_k, x, omega_s)) - E_k
E_p = evolve(add_breather(phi_k, x, omega_p)) - E_k
E_d = evolve(add_breather(phi_k, x, omega_d)) - E_k
E_f = evolve(add_breather(phi_k, x, omega_f)) - E_k
print(f"Single modes: E_s={E_s:.4f} E_p={E_p:.4f} E_d={E_d:.4f} E_f={E_f:.4f}")
print(flush=True)

# === TEST 1: f + d + p (Bi/Tl-like) ===
# How does having BOTH f and d in the core affect p-valence?
# Is the combined effect additive or is there a three-body term?

print("\n--- TEST 1: Three-body coupling (f + d + p) ---")
print("Does f+d together screen p differently than f alone + d alone?")

# d alone -> p
phi_d1 = add_breather(phi_k, x, omega_d)
E_d1 = evolve(phi_d1)
E_p_with_d = evolve(add_breather(phi_d1, x, omega_p)) - E_d1

# f alone -> p
phi_f1 = add_breather(phi_k, x, omega_f)
E_f1 = evolve(phi_f1)
E_p_with_f = evolve(add_breather(phi_f1, x, omega_p)) - E_f1

# d + f together -> p (three-body!)
phi_df = add_breather(phi_k, x, omega_d)
phi_df = add_breather(phi_df, x, omega_f)
E_df = evolve(phi_df)
E_p_with_df = evolve(add_breather(phi_df, x, omega_p)) - E_df

# Additive prediction: if no three-body, E_p_with_df = E_p + (d->p shift) + (f->p shift)
C_d_p = E_p_with_d / E_p
C_f_p = E_p_with_f / E_p
C_df_p = E_p_with_df / E_p
C_additive = (C_d_p - 1) + (C_f_p - 1) + 1  # sum of individual shifts

print(f"  d alone -> p:  C = {C_d_p:.4f}")
print(f"  f alone -> p:  C = {C_f_p:.4f}")
print(f"  d+f -> p:      C = {C_df_p:.4f}")
print(f"  Additive pred: C = {C_additive:.4f}")
print(f"  THREE-BODY:    {C_df_p - C_additive:+.4f} ({(C_df_p-C_additive)/C_additive*100:+.1f}%)")
print(flush=True)

# === TEST 2: Multiple d-modes + f + p (heavy p-block like Bi) ===
print("\n--- TEST 2: 5d + f -> p (Bi-like, heavy p-block) ---")

# Build 5 d-modes + 1 f-mode core
phi_core = phi_k.copy()
for i in range(5):
    phi_core = add_breather(phi_core, x, omega_d, x0=0.2*i-0.4)
phi_core = add_breather(phi_core, x, omega_f)
E_core = evolve(phi_core)

# Add p-mode
E_p_heavy = evolve(add_breather(phi_core, x, omega_p)) - E_core
C_heavy_p = E_p_heavy / E_p

# Compare: 5d alone -> p
phi_5d = phi_k.copy()
for i in range(5):
    phi_5d = add_breather(phi_5d, x, omega_d, x0=0.2*i-0.4)
E_5d = evolve(phi_5d)
E_p_5d = evolve(add_breather(phi_5d, x, omega_p)) - E_5d
C_5d_p = E_p_5d / E_p

print(f"  5d alone -> p: C = {C_5d_p:.4f}")
print(f"  5d+f -> p:     C = {C_heavy_p:.4f}")
print(f"  f contribution: {C_heavy_p - C_5d_p:+.4f}")
print(flush=True)

# === TEST 3: Same but for s-valence (Lu-like: f + d -> s) ===
print("\n--- TEST 3: f + d -> s (Lu-like) ---")

# f alone -> s
E_s_with_f = evolve(add_breather(phi_f1, x, omega_s)) - E_f1
C_f_s = E_s_with_f / E_s

# d alone -> s
E_s_with_d = evolve(add_breather(phi_d1, x, omega_s)) - E_d1
C_d_s = E_s_with_d / E_s

# d + f -> s
E_s_with_df = evolve(add_breather(phi_df, x, omega_s)) - E_df
C_df_s = E_s_with_df / E_s
C_additive_s = (C_d_s - 1) + (C_f_s - 1) + 1

print(f"  d alone -> s:  C = {C_d_s:.4f}")
print(f"  f alone -> s:  C = {C_f_s:.4f}")
print(f"  d+f -> s:      C = {C_df_s:.4f}")
print(f"  Additive pred: C = {C_additive_s:.4f}")
print(f"  THREE-BODY:    {C_df_s - C_additive_s:+.4f} ({(C_df_s-C_additive_s)/C_additive_s*100:+.1f}%)")
print(flush=True)

# === TEST 4: Half-fill effect ===
# Does adding 3 p-modes (half-fill) behave differently than 1+1+1?
print("\n--- TEST 4: Half-fill coupling (3 p-modes in d10 core) ---")

# Build d10-like core (5 d-modes)
phi_d10 = phi_k.copy()
for i in range(5):
    phi_d10 = add_breather(phi_d10, x, omega_d, x0=0.2*i-0.4)
E_d10 = evolve(phi_d10)

# Add p-modes one at a time
prev_E = E_d10
p_deltas = []
for ip in range(1, 4):
    phi_prev = phi_d10.copy()
    for j in range(ip):
        phi_prev = add_breather(phi_prev, x, omega_p, x0=0.3*j-0.3)
    E_new = evolve(phi_prev)
    E_p_incremental = E_new - prev_E
    C_inc = E_p_incremental / E_p
    p_deltas.append(C_inc)
    print(f"  p-mode {ip} (in d10 core): C = {C_inc:.4f}")
    prev_E = E_new

print(f"\n  Ratio p2/p1: {p_deltas[1]/p_deltas[0]:.4f}")
print(f"  Ratio p3/p2: {p_deltas[2]/p_deltas[1]:.4f}")
print(f"  At half-fill (p3): coupling = {p_deltas[2]:.4f}")
print(f"  If ratio = d: expect {p_deltas[0]*d**2:.4f} for p3")
print(flush=True)

print("\n" + "=" * 60)
print("SUMMARY")
print("  Three-body term tells us if f+d coupling is additive.")
print("  If NOT additive, the formula needs a cross-term for heavy atoms.")
print("  Half-fill tells us how exchange modifies in d10 core.")
