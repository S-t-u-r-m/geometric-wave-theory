"""
d->s vs d->p Coupling Simulation
================================
Measure how a d-mode (3rd harmonic) modifies:
  1. An s-mode (1st harmonic) — for TMs like Zn
  2. A p-mode (2nd harmonic) — for 4p block like Ga

This should reveal the t2g/eg split from pure wave dynamics.
Write results incrementally.
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


omega_s = np.cos(1 * gamma_sg)  # s-mode (n=1)
omega_p = np.cos(2 * gamma_sg)  # p-mode (n=2)
omega_d = np.cos(3 * gamma_sg)  # d-mode (n=3)
omega_f = np.cos(4 * gamma_sg)  # f-mode (n=4)

print("d->s vs d->p Coupling from Sine-Gordon Dynamics")
print("=" * 60)
print(f"  omega_s={omega_s:.6f}  omega_p={omega_p:.6f}  omega_d={omega_d:.6f}")
print(flush=True)

# Sweep Z values
for Z in [5.0, 10.0, 20.0]:
    set_Z(Z)
    print(f"\n--- Z = {Z:.0f} ---")
    t0 = time.time()

    phi_k = kink(x, Z)
    E_k = evolve(phi_k)

    # Single mode binding energies
    E_s_bare = evolve(add_breather(phi_k, x, omega_s)) - E_k
    E_p_bare = evolve(add_breather(phi_k, x, omega_p)) - E_k
    E_d_bare = evolve(add_breather(phi_k, x, omega_d)) - E_k

    print(f"  Singles: E_s={E_s_bare:.4f}  E_p={E_p_bare:.4f}  E_d={E_d_bare:.4f}")

    # d-mode present, then add s-mode (d->s screening)
    phi_kd = add_breather(phi_k, x, omega_d)
    E_kd = evolve(phi_kd)

    phi_kds = add_breather(phi_kd, x, omega_s)
    E_kds = evolve(phi_kds)
    E_s_with_d = E_kds - E_kd
    C_d_to_s = E_s_with_d / E_s_bare

    # d-mode present, then add p-mode (d->p screening)
    phi_kdp = add_breather(phi_kd, x, omega_p)
    E_kdp = evolve(phi_kdp)
    E_p_with_d = E_kdp - E_kd
    C_d_to_p = E_p_with_d / E_p_bare

    # s-mode present, then add d-mode (s->d screening)
    phi_ks = add_breather(phi_k, x, omega_s)
    E_ks = evolve(phi_ks)

    phi_ksd = add_breather(phi_ks, x, omega_d)
    E_ksd = evolve(phi_ksd)
    E_d_with_s = E_ksd - E_ks
    C_s_to_d = E_d_with_s / E_d_bare

    # p-mode present, then add d-mode (p->d screening)
    phi_kp = add_breather(phi_k, x, omega_p)
    E_kp = evolve(phi_kp)

    phi_kpd = add_breather(phi_kp, x, omega_d)
    E_kpd = evolve(phi_kpd)
    E_d_with_p = E_kpd - E_kp
    C_p_to_d = E_d_with_p / E_d_bare

    print(f"\n  Coupling ratios (C = E_with_other / E_bare):")
    print(f"    d present -> s binds at: C(d->s) = {C_d_to_s:.4f}")
    print(f"    d present -> p binds at: C(d->p) = {C_d_to_p:.4f}")
    print(f"    s present -> d binds at: C(s->d) = {C_s_to_d:.4f}")
    print(f"    p present -> d binds at: C(p->d) = {C_p_to_d:.4f}")

    print(f"\n  Net coupling (W = C - 1):")
    W_d_s = C_d_to_s - 1
    W_d_p = C_d_to_p - 1
    W_s_d = C_s_to_d - 1
    W_p_d = C_p_to_d - 1
    print(f"    W(d->s) = {W_d_s:+.4f}  W(s->d) = {W_s_d:+.4f}  ratio = {W_d_s/W_s_d:.4f}")
    print(f"    W(d->p) = {W_d_p:+.4f}  W(p->d) = {W_p_d:+.4f}  ratio = {W_d_p/W_p_d:.4f}")

    print(f"\n  KEY: d screens s vs p differently?")
    print(f"    W(d->s) / W(d->p) = {W_d_s/W_d_p:.4f}")
    print(f"    This ratio IS the t2g/eg split from wave dynamics!")
    print(f"    v19: d->s uses w_delta=-0.5, d->p uses different (t2g/eg)")

    # Also test: multiple d-modes (simulating d10 shell)
    print(f"\n  Multiple d-modes (simulating filled d-shell):")
    phi_d5 = phi_k.copy()
    for i in range(5):
        phi_d5 = add_breather(phi_d5, x, omega_d, x0=0.2*i - 0.4)
    E_d5 = evolve(phi_d5)

    phi_d5s = add_breather(phi_d5, x, omega_s)
    E_d5s = evolve(phi_d5s)
    E_s_with_d5 = E_d5s - E_d5
    C_d5_s = E_s_with_d5 / E_s_bare

    phi_d5p = add_breather(phi_d5, x, omega_p)
    E_d5p = evolve(phi_d5p)
    E_p_with_d5 = E_d5p - E_d5
    C_d5_p = E_p_with_d5 / E_p_bare

    print(f"    5 d-modes -> s: C = {C_d5_s:.4f}  W = {C_d5_s-1:+.4f}")
    print(f"    5 d-modes -> p: C = {C_d5_p:.4f}  W = {C_d5_p-1:+.4f}")
    print(f"    Ratio W(5d->s)/W(5d->p) = {(C_d5_s-1)/(C_d5_p-1):.4f}")

    # f->s, f->p, f->d for comparison
    print(f"\n  f-mode coupling:")
    phi_kf = add_breather(phi_k, x, omega_f)
    E_kf = evolve(phi_kf)

    for val_label, omega_val in [('s', omega_s), ('p', omega_p), ('d', omega_d)]:
        phi_kfv = add_breather(phi_kf, x, omega_val)
        E_kfv = evolve(phi_kfv)
        E_val_with_f = E_kfv - E_kf
        E_val_bare_lookup = {'s': E_s_bare, 'p': E_p_bare, 'd': E_d_bare}[val_label]
        C_f_val = E_val_with_f / E_val_bare_lookup
        print(f"    f->{val_label}: C = {C_f_val:.4f}  W = {C_f_val-1:+.4f}")

    print(f"\n  Time: {time.time()-t0:.1f}s", flush=True)

print("\n" + "=" * 60)
print("SUMMARY")
print("If W(d->s)/W(d->p) != 1, the d-mode couples differently to")
print("s vs p valence. This IS the t2g/eg effect from wave dynamics.")
