#!/usr/bin/env python3
"""
Lu simulation v2: f14 core coupling to d1 valence.

Key fix: f-shell must be INSIDE d-shell (4f inside 5d).
Use different effective Z for inner (f) vs outer (d) shells,
so the f-electrons are buried deep and the d-electron sits outside.

Measure: how does the f14 charge density modify the d-electron
binding energy? Compare screening/anti-screening with formula.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np

try:
    import cupy as cp
    GPU = True
    print("GPU: CuPy available")
except ImportError:
    cp = np
    GPU = False
    print("GPU: using CPU")

# === PARAMETERS ===
N = 1024
dr = 0.02
dt = 0.00005
n_relax = 80000

r = cp.arange(1, N+1, dtype=cp.float64) * dr
d = 3

def normalize(psi, r, dr):
    norm = cp.sqrt(cp.sum(cp.abs(psi)**2 * r**2 * dr))
    return psi / (norm + 1e-30)

def kinetic(psi, dr):
    d2 = cp.zeros_like(psi)
    d2[1:-1] = (psi[2:] - 2*psi[1:-1] + psi[:-2]) / dr**2
    d2[0] = (psi[1] - 2*psi[0]) / dr**2
    d2[-1] = (-2*psi[-1] + psi[-2]) / dr**2
    return -0.5 * d2

def evolve(psi, V, r, dr, dt):
    Hpsi = kinetic(psi, dr) + V * psi
    return normalize(psi - dt * Hpsi, r, dr)

def energy(psi, V, r, dr):
    Kpsi = kinetic(psi, dr)
    return float(cp.sum(psi * (Kpsi + V * psi) * r**2 * dr))

def ground_state(V, r, dr, dt, n_steps, l=0, label=""):
    psi = r**(l+1) * cp.exp(-r * 0.5)
    psi = normalize(psi, r, dr)
    for step in range(n_steps):
        psi = evolve(psi, V, r, dr, dt)
    E = energy(psi, V, r, dr)
    rho = cp.abs(psi)**2
    r_mean = float(cp.sum(rho * r**3 * dr) / cp.sum(rho * r**2 * dr))
    if label:
        print(f"  {label}: E = {E:.6f}, <r> = {r_mean:.3f}")
    return psi, E, r_mean

def screening_pot(rho, r, dr, n_e):
    Q_in = cp.cumsum(rho * r**2 * dr) * n_e
    V_in = Q_in / r
    rho_r = rho * r * dr * n_e
    V_out = cp.flip(cp.cumsum(cp.flip(rho_r)))
    return V_in + V_out

print(f"\nGrid: {N} pts, dr={dr}, rmax={N*dr}")
print(f"Relaxation: {n_relax} steps, dt={dt}")

# === REALISTIC SETUP FOR Lu ===
# In Lu: 4f14 is deeply buried, 5d1 is outer valence
# The 4f sees a large effective charge (Z_eff ~ 10-15 for 4f in Lu)
# The 5d sees a much smaller charge (Z_eff ~ 2-4 after screening)
#
# We simulate three scenarios:
# (a) Different core charges to position f inside d
# (b) Measure d-electron energy with and without f14 screening

# Effective charges (approximate, from Slater's rules as guide)
Z_f = 12.0   # 4f effective charge (deeply bound)
Z_d = 3.5    # 5d effective charge (valence, after core screening)
l_f = 3
l_d = 2

print(f"\nEffective charges: Z_f={Z_f} (inner), Z_d={Z_d} (outer)")

# Step 1: Find f-shell ground state (deeply bound)
print(f"\n{'='*60}")
print("Step 1: f-shell (l=3, inner shell)")
V_f = -Z_f / r + l_f * (l_f + 1) / (2 * r**2)
psi_f, E_f, r_f = ground_state(V_f, r, dr, dt, n_relax, l=l_f, label="f-shell")

# Step 2: d-electron ALONE (no f-shell)
print(f"\n{'='*60}")
print("Step 2: d-electron alone (l=2, no f-screening)")
V_d_bare = -Z_d / r + l_d * (l_d + 1) / (2 * r**2)
psi_d_bare, E_d_bare, r_d = ground_state(V_d_bare, r, dr, dt, n_relax, l=l_d, label="d (bare)")

print(f"\n  f inside d? r_f={r_f:.3f} < r_d={r_d:.3f}: {r_f < r_d}")

# Step 3: d-electron with f14 screening
print(f"\n{'='*60}")
print("Step 3: d-electron + f14 screening")
rho_f = cp.abs(psi_f)**2
V_f14 = screening_pot(rho_f, r, dr, n_e=14)
V_d_f14 = V_d_bare + V_f14
psi_d_f14, E_d_f14, r_d_f14 = ground_state(V_d_f14, r, dr, dt, n_relax, l=l_d, label="d+f14")

# Step 4: d-electron with f7 (half-fill)
print(f"\n{'='*60}")
print("Step 4: d-electron + f7 screening (half-fill)")
V_f7 = screening_pot(rho_f, r, dr, n_e=7)
V_d_f7 = V_d_bare + V_f7
psi_d_f7, E_d_f7, r_d_f7 = ground_state(V_d_f7, r, dr, dt, n_relax, l=l_d, label="d+f7")

# Step 5: d-electron with d10 screening (same l, inner shell)
print(f"\n{'='*60}")
print("Step 5: d-electron + d10 screening (inner d-shell)")
# Inner d-shell at higher Z (like 3d10 inside 5d)
Z_d_inner = 8.0
V_d_inner = -Z_d_inner / r + l_d * (l_d + 1) / (2 * r**2)
psi_d_inner, E_d_inner, r_d_inner = ground_state(V_d_inner, r, dr, dt, n_relax, l=l_d, label="d-inner")
rho_d_inner = cp.abs(psi_d_inner)**2
V_d10 = screening_pot(rho_d_inner, r, dr, n_e=10)
V_d_d10 = V_d_bare + V_d10
psi_d_d10, E_d_d10, r_d_d10 = ground_state(V_d_d10, r, dr, dt, n_relax, l=l_d, label="d+d10")

# Step 6: d-electron with p6 screening (inner p-shell)
print(f"\n{'='*60}")
print("Step 6: d-electron + p6 screening (inner p-shell)")
Z_p_inner = 10.0
V_p_inner = -Z_p_inner / r + 1 * (1 + 1) / (2 * r**2)
psi_p_inner, E_p_inner, r_p_inner = ground_state(V_p_inner, r, dr, dt, n_relax, l=1, label="p-inner")
rho_p_inner = cp.abs(psi_p_inner)**2
V_p6 = screening_pot(rho_p_inner, r, dr, n_e=6)
V_d_p6 = V_d_bare + V_p6
psi_d_p6, E_d_p6, r_d_p6 = ground_state(V_d_p6, r, dr, dt, n_relax, l=l_d, label="d+p6")

# === ANALYSIS ===
print(f"\n{'='*60}")
print("ANALYSIS")
print("=" * 60)

print(f"\n  Shell positions:")
print(f"    p-inner: <r> = {r_p_inner:.3f}")
print(f"    f-shell: <r> = {r_f:.3f}")
print(f"    d-inner: <r> = {r_d_inner:.3f}")
print(f"    d-outer: <r> = {r_d:.3f}")

dE_f14 = E_d_f14 - E_d_bare
dE_f7 = E_d_f7 - E_d_bare
dE_d10 = E_d_d10 - E_d_bare
dE_p6 = E_d_p6 - E_d_bare

print(f"\n  Energy shifts (positive = screening, negative = anti-screening):")
print(f"    f14: dE = {dE_f14:+.6f} ({dE_f14/abs(E_d_bare)*100:+.1f}%)")
print(f"    f7:  dE = {dE_f7:+.6f} ({dE_f7/abs(E_d_bare)*100:+.1f}%)")
print(f"    d10: dE = {dE_d10:+.6f} ({dE_d10/abs(E_d_bare)*100:+.1f}%)")
print(f"    p6:  dE = {dE_p6:+.6f} ({dE_p6/abs(E_d_bare)*100:+.1f}%)")

# Per-electron screening weight (normalized)
w_f14 = dE_f14 / (14 * abs(E_d_bare))
w_f7 = dE_f7 / (7 * abs(E_d_bare))
w_d10 = dE_d10 / (10 * abs(E_d_bare))
w_p6 = dE_p6 / (6 * abs(E_d_bare))

print(f"\n  Per-electron screening weight (dE / n_e / |E_bare|):")
print(f"    p6:  {w_p6:+.6f}")
print(f"    d10: {w_d10:+.6f}")
print(f"    f7:  {w_f7:+.6f}")
print(f"    f14: {w_f14:+.6f}")

print(f"\n  Screening ratios (relative to p6 per electron):")
if abs(w_p6) > 1e-10:
    print(f"    d10/p6: {w_d10/w_p6:.4f}")
    print(f"    f7/p6:  {w_f7/w_p6:.4f}")
    print(f"    f14/p6: {w_f14/w_p6:.4f}")

print(f"\n  GWT model predictions:")
print(f"    p screening: w_pi = +0.5 per channel")
print(f"    d screening: w_delta = -0.5 per channel (anti-screen)")
print(f"    f screening on d-valence: w_f_d = -1/6 = -0.167 per channel")
print(f"    Ratio d/p (model): -0.5/0.5 = -1.0")
print(f"    Ratio f/p (model): -0.167/0.5 = -0.333")

# CRITICAL: is f14 screening or anti-screening the d-electron?
print(f"\n  CRITICAL QUESTION: Does f14 SCREEN or ANTI-SCREEN the d-electron?")
if dE_f14 > 0:
    print(f"    ANSWER: f14 SCREENS (shields) the d-electron -> less bound")
    print(f"    This means w_f_d should be POSITIVE (screening)")
    print(f"    But our formula uses w_f_d = -1/6 (anti-screening)!")
    print(f"    THE FORMULA HAS THE WRONG SIGN FOR f->d COUPLING")
elif dE_f14 < 0:
    print(f"    ANSWER: f14 ANTI-SCREENS the d-electron -> more bound")
    print(f"    Consistent with w_f_d = -1/6 (anti-screening)")
else:
    print(f"    ANSWER: negligible effect")

# Pairing test
print(f"\n  PAIRING TEST: f14 vs 2*f7")
print(f"    dE(f14) = {dE_f14:.6f}")
print(f"    2*dE(f7) = {2*dE_f7:.6f}")
print(f"    Ratio = {dE_f14/(2*dE_f7):.4f}  (1.0 = linear, <1 = pairing reduces screening)")

# Spatial overlap
print(f"\n  OVERLAP: f-shell position relative to d-shell")
print(f"    r_f / r_d = {r_f / r_d:.4f}")
if r_f < r_d:
    print(f"    f is INSIDE d (correct for 4f inside 5d)")
    frac_inside = float(cp.sum(rho_f * (r < r_d).astype(float) * r**2 * dr))
    print(f"    Fraction of f-density inside d peak: {frac_inside:.3f}")
else:
    print(f"    WARNING: f is OUTSIDE d (unphysical)")
