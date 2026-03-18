#!/usr/bin/env python3
"""
Radial Curve Pattern Search — Can 40 curves become 1?
======================================================
Compute radial curves for key combinations and look for:
1. Does V(R; L) = f(R/L) * g(L)?  → scaling collapse
2. Does sigma/pi ratio stay constant?  → one curve + Oh weight
3. Does LP/bond ratio have a d=3 form? → amplitude scaling

If patterns exist, the entire force engine reduces to ONE master
curve + algebraic scaling. Same Wyler simplification.
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from math import factorial
import time

try:
    import cupy as cp; xp = cp; GPU = True; print("GPU active")
except: xp = np; GPU = False; print("CPU mode")

d = 3; V_0 = 1.0/np.pi**2
gamma_sg = np.pi/(2**(d+1)*np.pi - 2)
omega_p = float(np.cos(2*gamma_sg))
eps_p = float(np.sqrt(max(1-omega_p**2, 1e-12)))

# ============================================================
# 3D LATTICE SETUP
# ============================================================
N = 80  # balance speed vs accuracy
BOX = 8.0
dx = 2*BOX/N; dt = 0.2*dx

x1d = xp.linspace(-BOX, BOX, N, endpoint=False, dtype=np.float64)
X, Y, Z = xp.meshgrid(x1d, x1d, x1d, indexing='ij')

def make_mode(center_z, channel, L, amp):
    """Breather mode with width L, at center_z, in given channel."""
    Zs = Z - center_z
    R = xp.sqrt(X**2 + Y**2 + Zs**2) + 1e-10
    # Scale eps by L: wider breather = smaller eps
    eps_scaled = eps_p / L  # breather width scales with L
    omega_scaled = float(np.sqrt(max(1 - (eps_scaled*L)**2 / L**2, 0.01)))

    if channel == 'sigma':
        ang = Zs / R
    else:  # pi
        ang = X / R

    rad = amp / (omega_p * xp.cosh(eps_p * R / L) + 1e-10)
    return (4.0/np.pi) * xp.arctan(eps_p * ang * rad)

def total_energy(phi):
    GE = 0.5*((xp.roll(phi,1,0)-phi)**2+(xp.roll(phi,1,1)-phi)**2+(xp.roll(phi,1,2)-phi)**2)
    PE = V_0*(1-xp.cos(np.pi*phi))
    return float(xp.sum(GE+PE)*dx**3)

def evolve(phi0, ns=2500, damp=0.02):
    phi = phi0.copy(); po = phi.copy(); es = []
    for s in range(ns):
        lap = (xp.roll(phi,1,0)+xp.roll(phi,-1,0)+xp.roll(phi,1,1)+
               xp.roll(phi,-1,1)+xp.roll(phi,1,2)+xp.roll(phi,-1,2)-6*phi)/dx**2
        f = (1.0/np.pi)*xp.sin(np.pi*phi)
        pn = (2-damp*dt)*phi-(1-damp*dt)*po+dt**2*(lap-f)
        po = phi.copy(); phi = pn
        if s > ns*3//4 and s%50==0: es.append(total_energy(phi))
    return np.mean(es) if es else total_energy(phi)

def compute_V(R, L_a, L_b, amp_a, amp_b, channel):
    """Interaction energy at separation R."""
    phi_a = make_mode(-R/2, channel, L_a, amp_a)
    phi_b = make_mode(+R/2, channel, L_b, amp_b)
    E_a = evolve(phi_a, ns=1500)
    E_b = evolve(phi_b, ns=1500)
    E_ab = evolve(phi_a + phi_b, ns=2500)
    return E_ab - E_a - E_b

# ============================================================
# PATTERN 1: L-scaling — does V(R; L) = f(R/L)?
# ============================================================
print("="*70)
print("  PATTERN SEARCH: Can 40 radial curves become 1?")
print("="*70)
print()

AMP = 1.0  # fixed amplitude for comparison

# Compute sigma bond curve for different L values
print("PATTERN 1: L-scaling (sigma, half+half, amp=1.0)")
print("  Does V(R; L) collapse when plotted as V vs R/L?")
print()

R_points = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0])

curves = {}
for L in [0.5, 1.0, 1.5]:
    print(f"  L = {L}:")
    V_arr = []
    for R in R_points:
        t0 = time.time()
        V = compute_V(R, L, L, AMP, AMP, 'sigma')
        V_arr.append(V)
        print(f"    R={R:.1f} R/L={R/L:.2f} V={V:+.6f} ({time.time()-t0:.0f}s)")
    curves[f'sig_L{L}'] = np.array(V_arr)
    print()

# Check scaling: normalize each curve by its value at R/L = 2
print("  SCALING TEST: V / V(R/L=2) for each L")
print(f"  {'R/L':>6}", end='')
for L in [0.5, 1.0, 1.5]:
    print(f"  {'L='+str(L):>10}", end='')
print()

for i, R in enumerate(R_points):
    for L in [0.5, 1.0, 1.5]:
        rl = R/L
        # Find value at R/L closest to 2 for normalization
    # Just print the raw values at matching R/L
    pass

# Better: print V at same R/L values
print()
print("  R/L-matched comparison (V values at same R/L):")
# For L=0.5: R=1.0 → R/L=2.0
# For L=1.0: R=2.0 → R/L=2.0
# For L=1.5: R=3.0 → R/L=2.0
rl_targets = [2.0, 3.0, 4.0]
print(f"  {'R/L':>6}", end='')
for L in [0.5, 1.0, 1.5]:
    print(f"  L={L:>4}", end='')
print("  ratio(L1/L0.5)")
for rl in rl_targets:
    print(f"  {rl:>6.1f}", end='')
    vals = []
    for L in [0.5, 1.0, 1.5]:
        R_needed = rl * L
        # Find closest R in our grid
        idx = np.argmin(np.abs(R_points - R_needed))
        if abs(R_points[idx] - R_needed) < 0.1:
            V = curves[f'sig_L{L}'][idx]
            vals.append(V)
            print(f"  {V:+.5f}", end='')
        else:
            vals.append(None)
            print(f"  {'N/A':>8}", end='')
    if vals[0] and vals[1] and vals[0] != 0:
        print(f"  {vals[1]/vals[0]:>6.3f}", end='')
    print()


# ============================================================
# PATTERN 2: sigma vs pi ratio
# ============================================================
print()
print("PATTERN 2: sigma/pi ratio (L=1.0, half+half, amp=1.0)")
print("  Oh predicts: pi/sigma = cos(pi/3) = 0.5")
print()

print(f"  {'R':>5} {'V_sigma':>10} {'V_pi':>10} {'pi/sig':>8}")
for R in [1.5, 2.0, 2.5, 3.0, 4.0]:
    t0 = time.time()
    V_s = compute_V(R, 1.0, 1.0, AMP, AMP, 'sigma')
    V_p = compute_V(R, 1.0, 1.0, AMP, AMP, 'pi')
    ratio = V_p/V_s if abs(V_s) > 1e-10 else 0
    print(f"  {R:5.1f} {V_s:+10.6f} {V_p:+10.6f} {ratio:+8.3f}  ({time.time()-t0:.0f}s)")


# ============================================================
# PATTERN 3: bond vs LP (amplitude scaling)
# ============================================================
print()
print("PATTERN 3: bond vs LP (L=1.0, sigma channel)")
print("  bond = half+half (amp=1), LP = full+full (amp=2)")
print("  If linear: LP/bond = (2/1)^2 = 4")
print()

print(f"  {'R':>5} {'V_bond':>10} {'V_LP':>10} {'LP/bond':>8} {'LP/4bond':>8}")
for R in [1.5, 2.0, 2.5, 3.0, 4.0]:
    t0 = time.time()
    V_b = compute_V(R, 1.0, 1.0, AMP, AMP, 'sigma')
    V_lp = compute_V(R, 1.0, 1.0, AMP*2, AMP*2, 'sigma')
    ratio = V_lp/V_b if abs(V_b) > 1e-10 else 0
    quarter = V_lp/(4*V_b) if abs(V_b) > 1e-10 else 0
    print(f"  {R:5.1f} {V_b:+10.6f} {V_lp:+10.6f} {ratio:+8.3f} {quarter:+8.3f}  ({time.time()-t0:.0f}s)")

print()
print("PATTERN 3b: LP in pi channel (perpendicular — the repulsive case)")
print(f"  {'R':>5} {'V_bond_pi':>10} {'V_LP_pi':>10} {'LP/bond':>8}")
for R in [1.5, 2.0, 3.0]:
    t0 = time.time()
    V_bp = compute_V(R, 1.0, 1.0, AMP, AMP, 'pi')
    V_lpp = compute_V(R, 1.0, 1.0, AMP*2, AMP*2, 'pi')
    ratio = V_lpp/V_bp if abs(V_bp) > 1e-10 else 0
    print(f"  {R:5.1f} {V_bp:+10.6f} {V_lpp:+10.6f} {ratio:+8.3f}  ({time.time()-t0:.0f}s)")
