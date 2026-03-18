"""
Kink-Breather Simulation v2 — Targeted Physics Extraction
===========================================================

Sim1 showed: pairing required, screening strong, alpha~0.14 for 1D.
Now dig deeper into specific questions:

A) BINDING vs MODE NUMBER — what's the functional form?
   Does E_bind follow sin(n*gamma), cos(n*gamma), or something else?
   This maps to: why does alpha depend on mode structure?

B) PAIRED vs UNPAIRED binding to kink — quantify the pairing enhancement
   Ratio of paired/unpaired binding = the "pairing penalty" in Z_eff formula
   Our formula says: first lock weight = w_pi, rest = 1+w_pi. Does sim agree?

C) TWO DIFFERENT modes at kink — cross-mode screening
   When mode_inner is present, how much does mode_outer's binding change?
   This is the core screening S_core = channels * w_pi. Does sim give w_pi?

D) BINDING ENERGY RATIOS — the alpha measurement
   For a paired mode at Z kinks, measure E_bind(Z) / E_bind(1).
   Fit alpha from this. Compare for different modes.

E) BREATHER OFFSET FROM KINK — radial structure
   Place breather at distance r from kink. Binding vs r.
   Does it decay as exp(-eta*r)? What's the equilibrium distance?
"""

import numpy as np
import time as timer

pi = np.pi
d = 3

V_0 = 1.0 / pi**2
N_br = 24
gamma = pi / (N_br + 1)

# Lattice
N = 4000
L = 160.0
dx = L / N
dt = 0.08 * dx
x = np.linspace(-L/2, L/2, N)

n_evolve = 8000
measure_interval = 200


def kink_phi(x_arr, center=0.0):
    xi = x_arr - center
    return (4.0 / pi) * np.arctan(np.exp(np.clip(xi, -50, 50)))

def breather_phi(x_arr, n_mode, center=0.0, t=0.0, sign=1.0):
    w = np.cos(n_mode * gamma)
    eta = np.sin(n_mode * gamma)
    xi = eta * (x_arr - center)
    cosh_xi = np.cosh(np.clip(xi, -50, 50))
    arg = sign * (eta / w) * np.sin(w * t) / cosh_xi
    return (4.0 / pi) * np.arctan(arg)

def breather_dphi(x_arr, n_mode, center=0.0, t=0.0, sign=1.0):
    w = np.cos(n_mode * gamma)
    eta = np.sin(n_mode * gamma)
    xi = eta * (x_arr - center)
    cosh_xi = np.cosh(np.clip(xi, -50, 50))
    f = sign * (eta / w) * np.sin(w * t) / cosh_xi
    dfdt = sign * eta * np.cos(w * t) / cosh_xi
    return (4.0 / pi) * dfdt / (1 + f**2)

def total_energy(phi, dphi_dt):
    E_kin = 0.5 * np.sum(dphi_dt**2) * dx
    grad = np.diff(phi) / dx
    E_grad = 0.5 * np.sum(grad**2) * dx
    E_pot = V_0 * np.sum(1.0 - np.cos(pi * phi)) * dx
    return E_kin + E_grad + E_pot

def leapfrog_step(phi, dphi_dt, dt_step):
    acc = np.zeros_like(phi)
    acc[1:-1] = (phi[:-2] - 2*phi[1:-1] + phi[2:]) / dx**2
    acc -= V_0 * pi * np.sin(pi * phi)
    dphi_half = dphi_dt + 0.5 * dt_step * acc
    phi_new = phi + dt_step * dphi_half
    acc_new = np.zeros_like(phi_new)
    acc_new[1:-1] = (phi_new[:-2] - 2*phi_new[1:-1] + phi_new[2:]) / dx**2
    acc_new -= V_0 * pi * np.sin(pi * phi_new)
    dphi_new = dphi_half + 0.5 * dt_step * acc_new
    return phi_new, dphi_new

def evolve_and_measure(phi0, dphi0, n_steps=8000):
    phi = phi0.copy()
    dphi = dphi0.copy()
    phi_left, phi_right = phi[0], phi[-1]
    energies = []
    for step in range(n_steps):
        phi, dphi = leapfrog_step(phi, dphi, dt)
        phi[0] = phi_left; phi[-1] = phi_right
        dphi[0] = 0.0; dphi[-1] = 0.0
        if step % measure_interval == 0 and step > n_steps // 3:
            energies.append(total_energy(phi, dphi))
    return np.mean(energies) if energies else total_energy(phi0, dphi0)

t_init = 0.25  # fraction of half-period


# =============================================================================
# BASELINES
# =============================================================================
print("=" * 80)
print("  KINK-BREATHER SIM v2 — PHYSICS EXTRACTION")
print("=" * 80)

# Kink energy
phi_k = kink_phi(x)
E_kink = evolve_and_measure(phi_k, np.zeros_like(x))
print(f"\n  E_kink = {E_kink:.6f} (theory: {8/pi**2:.6f})")

# Single breather energies
print(f"\n  Single breather energies:")
single_E = {}
for n in range(1, 11):
    w_n = np.cos(n * gamma)
    eta_n = np.sin(n * gamma)
    t0 = t_init * pi / w_n
    phi_br = breather_phi(x, n, t=t0)
    dphi_br = breather_dphi(x, n, t=t0)
    E = evolve_and_measure(phi_br, dphi_br)
    M_th = (16/pi**2) * np.sin(n * gamma)
    single_E[n] = E
    print(f"    n={n:2d}: E={E:.6f}, M_theory={M_th:.6f}, err={abs(E-M_th)/M_th*100:.2f}%")


# =============================================================================
# A) BINDING vs MODE NUMBER — functional form
# =============================================================================
print(f"\n{'='*80}")
print("  A) SINGLE-MODE BINDING TO KINK vs MODE NUMBER")
print("  E_bind(n) = E(kink+br_n) - E(kink) - E(br_n)")
print("  Looking for: does E_bind scale with sin, cos, eta, w?")
print(f"{'='*80}\n")

print(f"  {'n':>3} {'w_n':>7} {'eta_n':>7} {'E_bind':>10} "
      f"{'|Eb|/M':>8} {'|Eb|/sin':>8} {'|Eb|/cos':>8} {'|Eb|/eta^2':>10}")

bind_1 = {}
for n in range(1, 11):
    w_n = np.cos(n * gamma)
    eta_n = np.sin(n * gamma)
    t0 = t_init * pi / w_n

    phi0 = kink_phi(x) + breather_phi(x, n, t=t0)
    dphi0 = np.zeros_like(x) + breather_dphi(x, n, t=t0)
    E_kb = evolve_and_measure(phi0, dphi0)
    E_bind = E_kb - E_kink - single_E[n]
    bind_1[n] = E_bind

    M_n = (16/pi**2) * np.sin(n * gamma)
    r_M = abs(E_bind) / M_n
    r_sin = abs(E_bind) / np.sin(n * gamma) if np.sin(n * gamma) > 0.01 else 0
    r_cos = abs(E_bind) / np.cos(n * gamma) if np.cos(n * gamma) > 0.01 else 0
    r_eta2 = abs(E_bind) / eta_n**2

    print(f"  {n:3d} {w_n:7.4f} {eta_n:7.4f} {E_bind:+10.6f} "
          f"{r_M:8.4f} {r_sin:8.4f} {r_eta2:10.4f}")

# Check power law: E_bind ~ eta^p or sin^p
print(f"\n  Power law fit: E_bind = A * eta^p")
log_eta = np.array([np.log(np.sin(n * gamma)) for n in range(1, 9)])
log_Eb = np.array([np.log(abs(bind_1[n])) for n in range(1, 9)])
# Linear fit
p, logA = np.polyfit(log_eta, log_Eb, 1)
A = np.exp(logA)
print(f"  p = {p:.4f}, A = {A:.6f}")
print(f"  E_bind ~ {A:.4f} * eta^{p:.2f}")

# What GWT values is p close to?
print(f"\n  p = {p:.4f} candidates:")
for label, val in [('1', 1.0), ('2', 2.0), ('d-1=2', 2.0), ('d=3', 3.0),
                    ('1/Gamma=9', 9.0), ('pi', pi), ('2*pi', 2*pi)]:
    print(f"    {label:>12}: {val:.4f}  (delta={abs(p-val):.3f})")


# =============================================================================
# B) PAIRING ENHANCEMENT — quantify the ratio
# =============================================================================
print(f"\n{'='*80}")
print("  B) PAIRING ENHANCEMENT")
print("  E_pair = E(kink + br+ + br-) - E(kink) - 2*E(br)")
print("  E_unpr = E(kink + br+ + br+) - E(kink) - 2*E(br)")
print("  Enhancement = E_pair / (2 * E_bind_single)")
print(f"{'='*80}\n")

print(f"  {'n':>3} {'E_pair':>10} {'E_unpr':>10} {'2*E_1':>10} "
      f"{'pair/2E1':>9} {'unpr/2E1':>9} {'pair/unpr':>10}")

for n in range(1, 9):
    w_n = np.cos(n * gamma)
    t0 = t_init * pi / w_n

    # Paired (opposite sign)
    phi_p = kink_phi(x) + breather_phi(x, n, t=t0, sign=+1) + \
            breather_phi(x, n, t=t0*1.03, sign=-1)
    dphi_p = breather_dphi(x, n, t=t0, sign=+1) + \
             breather_dphi(x, n, t=t0*1.03, sign=-1)
    E_pair = evolve_and_measure(phi_p, dphi_p) - E_kink - 2*single_E[n]

    # Unpaired (same sign)
    phi_u = kink_phi(x) + breather_phi(x, n, t=t0, sign=+1) + \
            breather_phi(x, n, t=t0*1.03, sign=+1)
    dphi_u = breather_dphi(x, n, t=t0, sign=+1) + \
             breather_dphi(x, n, t=t0*1.03, sign=+1)
    E_unpr = evolve_and_measure(phi_u, dphi_u) - E_kink - 2*single_E[n]

    two_E1 = 2 * bind_1[n]
    r_pair = E_pair / two_E1 if abs(two_E1) > 1e-10 else float('nan')
    r_unpr = E_unpr / two_E1 if abs(two_E1) > 1e-10 else float('nan')
    r_pu = E_pair / E_unpr if abs(E_unpr) > 1e-10 else float('nan')

    print(f"  {n:3d} {E_pair:+10.6f} {E_unpr:+10.6f} {two_E1:+10.6f} "
          f"{r_pair:9.4f} {r_unpr:9.4f} {r_pu:10.4f}")


# =============================================================================
# C) CROSS-MODE SCREENING — inner mode reduces outer binding
# =============================================================================
print(f"\n{'='*80}")
print("  C) CROSS-MODE SCREENING")
print("  S(in,out) = E_bind(out|in present) / E_bind(out|alone)")
print("  This is the 'screening fraction' — what fraction of binding survives")
print("  GWT predicts: core channels screen at w_pi = 0.5")
print(f"{'='*80}\n")

# Kink + inner mode baselines
kink_mode_E = {}
for n in range(1, 9):
    w_n = np.cos(n * gamma)
    t0 = t_init * pi / w_n
    phi0 = kink_phi(x) + breather_phi(x, n, t=t0)
    dphi0 = breather_dphi(x, n, t=t0)
    kink_mode_E[n] = evolve_and_measure(phi0, dphi0)

print(f"  {'in':>3} {'out':>3} {'dE_bare':>10} {'dE_screen':>10} "
      f"{'S_frac':>8} {'1-S':>8} {'w_pi?':>6}")

# Test: each mode as inner, each higher mode as outer
for n_in in [1, 2, 3, 4, 5]:
    for n_out in [n_in+1, n_in+2, n_in+3]:
        if n_out > 8:
            continue
        w_in = np.cos(n_in * gamma)
        w_out = np.cos(n_out * gamma)
        t0_in = t_init * pi / w_in
        t0_out = t_init * pi / w_out

        # Kink + inner + outer
        phi0 = kink_phi(x) + breather_phi(x, n_in, t=t0_in) + \
               breather_phi(x, n_out, t=t0_out)
        dphi0 = breather_dphi(x, n_in, t=t0_in) + \
                breather_dphi(x, n_out, t=t0_out)
        E_both = evolve_and_measure(phi0, dphi0)

        # Binding of outer with inner present
        dE_screen = E_both - kink_mode_E[n_in] - single_E[n_out]
        # Bare binding
        dE_bare = bind_1[n_out]
        # Screening fraction
        S_frac = dE_screen / dE_bare if abs(dE_bare) > 1e-10 else float('nan')
        screen_amount = 1 - S_frac
        w_pi_match = "~w_pi" if abs(screen_amount - 0.5) < 0.15 else ""

        print(f"  {n_in:3d} {n_out:3d} {dE_bare:+10.6f} {dE_screen:+10.6f} "
              f"{S_frac:8.4f} {screen_amount:+8.4f} {w_pi_match:>6}")


# =============================================================================
# D) Z-SCALING with paired mode — better alpha measurement
# =============================================================================
print(f"\n{'='*80}")
print("  D) Z-SCALING: PAIRED MODE BINDING vs KINK COUNT")
print("  E_bind(Z) for a paired (+-) breather at Z kinks")
print("  Fit: E_bind ~ Z^p => alpha = p/2")
print(f"{'='*80}\n")

kink_spacing = 0.3

for test_n in [1, 2, 3, 5]:
    if test_n not in single_E:
        continue
    w_n = np.cos(test_n * gamma)
    t0 = t_init * pi / w_n

    print(f"\n  --- Mode n={test_n} (w={w_n:.4f}) ---")
    print(f"  {'Z':>4} {'E_bind_pair':>12} {'ratio':>8}")

    Z_vals = []
    Eb_vals = []

    for Z in [1, 2, 3, 4, 5, 6, 8]:
        # Z kinks
        phi_Z = np.zeros_like(x)
        for k in range(Z):
            c = (k - (Z-1)/2) * kink_spacing
            phi_Z += kink_phi(x, center=c)

        # Z kinks alone
        E_Z = evolve_and_measure(phi_Z, np.zeros_like(x))

        # Z kinks + paired breather
        phi_Zp = phi_Z + breather_phi(x, test_n, t=t0, sign=+1) + \
                 breather_phi(x, test_n, t=t0*1.03, sign=-1)
        dphi_Zp = breather_dphi(x, test_n, t=t0, sign=+1) + \
                  breather_dphi(x, test_n, t=t0*1.03, sign=-1)
        E_Zp = evolve_and_measure(phi_Zp, dphi_Zp)

        E_bind = E_Zp - E_Z - 2*single_E[test_n]
        Z_vals.append(Z)
        Eb_vals.append(abs(E_bind))

        ratio = E_bind / Eb_vals[0] if Eb_vals[0] > 1e-10 else 0
        print(f"  {Z:4d} {E_bind:+12.6f} {ratio:8.4f}")

    # Fit power law
    if len(Z_vals) > 2 and all(e > 1e-10 for e in Eb_vals):
        logZ = np.log(np.array(Z_vals, dtype=float))
        logE = np.log(np.array(Eb_vals))
        p_fit, _ = np.polyfit(logZ, logE, 1)
        alpha_fit = p_fit / 2
        print(f"  => p = {p_fit:.4f}, alpha = p/2 = {alpha_fit:.4f}")
        print(f"     GWT candidates: 1/d={1/d:.4f}, Gamma={1/d**2:.4f}, "
              f"w_pi/d={0.5/d:.4f}, (d-1)/d^2={2/9:.4f}")


# =============================================================================
# E) BREATHER OFFSET FROM KINK — equilibrium distance
# =============================================================================
print(f"\n{'='*80}")
print("  E) BINDING vs OFFSET — where does the breather sit?")
print("  Place breather at distance r from kink, measure binding")
print(f"{'='*80}\n")

for test_n in [1, 3, 5]:
    if test_n not in single_E:
        continue
    w_n = np.cos(test_n * gamma)
    eta_n = np.sin(test_n * gamma)
    t0 = t_init * pi / w_n

    print(f"\n  --- Mode n={test_n} (eta={eta_n:.4f}, width={1/eta_n:.2f}) ---")
    print(f"  {'r':>6} {'E_bind':>10} {'|Eb|/|Eb0|':>12} {'exp(-eta*r)':>12}")

    E_bind_0 = None
    for r in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0]:
        phi0 = kink_phi(x) + breather_phi(x, test_n, center=r, t=t0)
        dphi0 = breather_dphi(x, test_n, center=r, t=t0)
        E_kb = evolve_and_measure(phi0, dphi0)
        E_bind = E_kb - E_kink - single_E[test_n]

        if E_bind_0 is None:
            E_bind_0 = E_bind

        ratio = abs(E_bind / E_bind_0) if abs(E_bind_0) > 1e-10 else 0
        exp_decay = np.exp(-eta_n * r)

        print(f"  {r:6.1f} {E_bind:+10.6f} {ratio:12.4f} {exp_decay:12.4f}")


# =============================================================================
# F) TRIPLE MODE: inner pair + outer single (the atom!)
# =============================================================================
print(f"\n{'='*80}")
print("  F) ATOM MODEL: Kink + inner PAIR + outer SINGLE")
print("  Like He-core + valence: (1s^2) + outer mode")
print("  How much does the 1s pair screen the outer mode?")
print(f"{'='*80}\n")

# Inner pair: mode n_in, opposite signs
# Outer single: mode n_out
n_in = 1  # "1s pair"

print(f"  Inner pair: mode {n_in} (paired)")
print(f"  {'n_out':>5} {'dE_bare':>10} {'dE_w_pair':>10} {'screen':>8} {'1-screen':>8}")

# Kink + inner pair baseline
w_in = np.cos(n_in * gamma)
t0_in = t_init * pi / w_in
phi_kp = kink_phi(x) + breather_phi(x, n_in, t=t0_in, sign=+1) + \
         breather_phi(x, n_in, t=t0_in*1.03, sign=-1)
dphi_kp = breather_dphi(x, n_in, t=t0_in, sign=+1) + \
          breather_dphi(x, n_in, t=t0_in*1.03, sign=-1)
E_kink_pair = evolve_and_measure(phi_kp, dphi_kp)
print(f"  E(kink + 1s pair) = {E_kink_pair:.6f}")
print()

for n_out in range(2, 9):
    w_out = np.cos(n_out * gamma)
    t0_out = t_init * pi / w_out

    # Kink + pair + outer
    phi0 = kink_phi(x) + breather_phi(x, n_in, t=t0_in, sign=+1) + \
           breather_phi(x, n_in, t=t0_in*1.03, sign=-1) + \
           breather_phi(x, n_out, t=t0_out)
    dphi0 = breather_dphi(x, n_in, t=t0_in, sign=+1) + \
            breather_dphi(x, n_in, t=t0_in*1.03, sign=-1) + \
            breather_dphi(x, n_out, t=t0_out)
    E_all = evolve_and_measure(phi0, dphi0)

    dE_screen = E_all - E_kink_pair - single_E[n_out]
    dE_bare = bind_1[n_out]
    S = dE_screen / dE_bare if abs(dE_bare) > 1e-10 else float('nan')

    print(f"  {n_out:5d} {dE_bare:+10.6f} {dE_screen:+10.6f} "
          f"{S:8.4f} {1-S:+8.4f}")


# =============================================================================
# SUMMARY TABLE
# =============================================================================
print(f"\n{'='*80}")
print("  SUMMARY: KEY NUMBERS FOR Z_eff FORMULA")
print(f"{'='*80}")
print(f"\n  w_pi = cos(pi/d) = {np.cos(pi/d):.4f}")
print(f"  Gamma = 1/d^2 = {1/d**2:.4f}")
print(f"  1/d = {1/d:.4f}")
print()
