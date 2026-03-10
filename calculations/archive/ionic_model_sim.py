"""
GWT Ionic Model Simulation — GPU
==================================

Tests the TWO-STATE charge transfer model used in V6:
    q = delta_eps / sqrt(delta_eps^2 + (2*V_cov)^2)
    D_ion = c_ion * q^2 * 2*E_H / R

Key questions:
1. Does E_int scale linearly with |sin(phase)| (V6 assumption)?
   Or are there nonlinear corrections at large amplitude?
2. In an asymmetric potential (modeling electronegativity),
   does amplitude transfer match the two-state q formula?
3. Does the transferred amplitude depend on wave type (s vs p)?

If the formula overestimates charge transfer, D_ion is too large,
explaining the OH +9.5% overshoot.
"""

import cupy as cp
import numpy as np
import time as timer

pi = float(np.pi)
k_long = 1.0
k_trans = 0.5


def make_field(N, L, wave, center_z):
    """Create a 3-component wave field."""
    dx = L / N
    coords = cp.linspace(-L/2, L/2, N, dtype=cp.float32)
    X, Y, Z = cp.meshgrid(coords, coords, coords, indexing='ij')

    field = cp.zeros((*X.shape, 3), dtype=cp.float32)
    rz = Z - center_z
    r = cp.sqrt(X**2 + Y**2 + rz**2)
    r = cp.maximum(r, dx/2)

    n = wave['n']
    l = wave['l']
    width = wave['width']
    amp = wave['amplitude']

    envelope = amp * cp.exp(-r**2 / (2 * width**2))
    n_nodes = n - l - 1
    if n_nodes > 0:
        envelope *= cp.cos(n_nodes * pi * r / (3 * width))

    channels = wave.get('channels', {0: 1.0})
    ch_to_comp = {0: 2, 1: 0, 2: 1}

    for ch_idx, ch_amp in channels.items():
        comp = ch_to_comp[ch_idx]
        if l == 0:
            field[:,:,:,comp] += ch_amp * envelope
        elif l == 1:
            if ch_idx == 0:
                field[:,:,:,comp] += ch_amp * envelope * (rz / r)
            elif ch_idx == 1:
                field[:,:,:,comp] += ch_amp * envelope * (X / r)
            elif ch_idx == 2:
                field[:,:,:,comp] += ch_amp * envelope * (Y / r)
    return field, X, Y, Z, dx


def compute_energy(phi, dx):
    E = cp.float32(0.0)
    V_0 = 1.0 / pi**2
    for axis in range(3):
        dphi = cp.diff(phi, axis=axis)
        E += 0.5 * k_long * cp.sum(dphi[:,:,:,axis]**2) * dx
        for comp in range(3):
            if comp != axis:
                E += 0.5 * k_trans * cp.sum(dphi[:,:,:,comp]**2) * dx
    for comp in range(3):
        E += V_0 * cp.sum(1 - cp.cos(pi * phi[:,:,:,comp])) * dx**3
    return float(E)


def compute_force(phi, dx):
    force = cp.zeros_like(phi)
    V_0 = 1.0 / pi**2
    for axis in range(3):
        pp = cp.roll(phi, -1, axis=axis)
        pm = cp.roll(phi, 1, axis=axis)
        lap = (pp + pm - 2*phi) / dx**2
        force[:,:,:,axis] += k_long * lap[:,:,:,axis]
        for comp in range(3):
            if comp != axis:
                force[:,:,:,comp] += k_trans * lap[:,:,:,comp]
    for comp in range(3):
        force[:,:,:,comp] -= V_0 * pi * cp.sin(pi * phi[:,:,:,comp])
    return force


def evolve_and_measure(phi_init, dx, n_steps, dt_factor=0.12):
    dt = dt_factor * dx
    phi = phi_init.copy()
    dphi = cp.zeros_like(phi)
    force = compute_force(phi, dx)
    dphi += 0.5 * dt * force
    E_samples = []
    for step in range(n_steps):
        phi += dt * dphi
        force = compute_force(phi, dx)
        dphi += dt * force
        if step % 50 == 0 and step > n_steps // 4:
            dphi_half = dphi - 0.5 * dt * force
            E_kin = 0.5 * float(cp.sum(dphi_half**2)) * dx**3
            E_pot = compute_energy(phi, dx)
            E_samples.append(E_kin + E_pot)
    return np.mean(E_samples) if E_samples else compute_energy(phi_init, dx)


def run_pair(N, L, R, wave1, wave2, n_steps):
    """Compute interaction energy for a pair of waves."""
    dx = L / N
    phi1, X, Y, Z, _ = make_field(N, L, wave1, -R/2)
    E1 = evolve_and_measure(phi1, dx, n_steps)

    phi2, _, _, _, _ = make_field(N, L, wave2, +R/2)
    E2 = evolve_and_measure(phi2, dx, n_steps)

    phi_both = phi1 + phi2
    E_both = evolve_and_measure(phi_both, dx, n_steps)

    E_int = E_both - E1 - E2
    cp.cuda.Stream.null.synchronize()
    return E_int


# =============================================================================
print("=" * 80)
print("  IONIC MODEL SIMULATION")
print("=" * 80)

_ = cp.zeros(10)
cp.cuda.Stream.null.synchronize()

N = 48
L = 16.0
n_steps = 2500


# =============================================================================
# TEST 1: Amplitude scaling — is E_int proportional to amp^2?
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 1: Amplitude scaling for s+p bond at R=4")
print(f"  V6 assumes D_cov ~ amp^2 (through E_scale)")
print(f"  Nonlinear potential: deviations at large amp?")
print(f"{'='*80}")

print(f"\n  {'amp':>5}  {'E_int':>10}  {'E/amp^2':>10}  {'ratio_to_ref':>12}")

E_ref_per_amp2 = None
for amp_val in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]:
    w_s = {'amplitude': amp_val, 'width': 2.0, 'n': 1, 'l': 0, 'channels': {0: 1.0}}
    w_p = {'amplitude': amp_val, 'width': 2.0, 'n': 2, 'l': 1, 'channels': {0: 1.0}}
    t0 = timer.time()
    E_int = run_pair(N, L, 4.0, w_s, w_p, n_steps)
    dt = timer.time() - t0

    E_per_amp2 = E_int / amp_val**2 if amp_val > 0 else 0
    if E_ref_per_amp2 is None:
        E_ref_per_amp2 = E_per_amp2

    ratio = E_per_amp2 / E_ref_per_amp2 if abs(E_ref_per_amp2) > 1e-10 else 0
    print(f"  {amp_val:5.2f}  {E_int:10.4f}  {E_per_amp2:10.4f}  {ratio:12.4f}  ({dt:.1f}s)",
          flush=True)


# =============================================================================
# TEST 2: Asymmetric amplitudes — model charge transfer
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 2: Asymmetric amplitudes at R=4")
print(f"  Higher amplitude on one side models electronegativity")
print(f"  Does asymmetry change E_int beyond simple scaling?")
print(f"{'='*80}")

print(f"\n  {'amp1':>5}  {'amp2':>5}  {'E_int':>10}  {'E/(a1*a2)':>10}  {'sym_ref':>10}  {'ratio':>8}")

# Reference: symmetric
w_s_sym = {'amplitude': 0.5, 'width': 2.0, 'n': 1, 'l': 0, 'channels': {0: 1.0}}
w_p_sym = {'amplitude': 0.5, 'width': 2.0, 'n': 2, 'l': 1, 'channels': {0: 1.0}}
E_sym = run_pair(N, L, 4.0, w_s_sym, w_p_sym, n_steps)
E_sym_norm = E_sym / (0.5 * 0.5)

for a1, a2 in [(0.5, 0.5), (0.3, 0.7), (0.2, 0.8), (0.1, 0.9),
               (0.7, 0.3), (0.8, 0.2), (0.4, 0.6), (0.6, 0.4)]:
    w_s = {'amplitude': a1, 'width': 2.0, 'n': 1, 'l': 0, 'channels': {0: 1.0}}
    w_p = {'amplitude': a2, 'width': 2.0, 'n': 2, 'l': 1, 'channels': {0: 1.0}}
    t0 = timer.time()
    E_int = run_pair(N, L, 4.0, w_s, w_p, n_steps)
    dt = timer.time() - t0

    E_norm = E_int / (a1 * a2) if a1*a2 > 0 else 0
    ratio = E_norm / E_sym_norm if abs(E_sym_norm) > 1e-10 else 0
    print(f"  {a1:5.2f}  {a2:5.2f}  {E_int:10.4f}  {E_norm:10.4f}  {E_sym_norm:10.4f}  {ratio:8.4f}  ({dt:.1f}s)",
          flush=True)


# =============================================================================
# TEST 3: Amplitude transfer in asymmetric potential
# Model charge transfer: add a site-dependent potential V_asym
# that makes one side deeper (more electronegative)
# Measure how much wave amplitude shifts from one side to other
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 3: Charge transfer test — evolve with asymmetric potential")
print(f"  Start with symmetric waves, evolve with deeper potential on one side")
print(f"  Measure how much amplitude transfers")
print(f"{'='*80}")

dx = L / N
coords = cp.linspace(-L/2, L/2, N, dtype=cp.float32)
X, Y, Z = cp.meshgrid(coords, coords, coords, indexing='ij')

R_test = 4.0
w_base = 2.0
amp_base = 0.5

# Create initial symmetric configuration
w1 = {'amplitude': amp_base, 'width': w_base, 'n': 1, 'l': 0, 'channels': {0: 1.0}}
w2 = {'amplitude': amp_base, 'width': w_base, 'n': 2, 'l': 1, 'channels': {0: 1.0}}
phi1, _, _, _, _ = make_field(N, L, w1, -R_test/2)
phi2, _, _, _, _ = make_field(N, L, w2, +R_test/2)
phi_init = phi1 + phi2

# Measure initial amplitude on each side
left_mask = (Z < 0).astype(cp.float32)
right_mask = (Z >= 0).astype(cp.float32)

def measure_sides(phi):
    """Measure total |phi|^2 on each side of z=0."""
    phi_sq = cp.sum(phi**2, axis=3)  # sum over components
    left = float(cp.sum(phi_sq * left_mask)) * dx**3
    right = float(cp.sum(phi_sq * right_mask)) * dx**3
    return left, right

L0, R0 = measure_sides(phi_init)
print(f"\n  Initial: left={L0:.4f}, right={R0:.4f}, ratio L/R={L0/R0:.4f}")

# Now evolve with modified force that includes asymmetric potential
# V_asym = -epsilon * phi^2 for z > 0 (deeper well on right = more electronegative)
# This is like adding a "mass" term that makes the right side prefer larger phi

print(f"\n  {'epsilon':>8}  {'L_final':>8}  {'R_final':>8}  {'L/R':>8}  {'transfer%':>10}")

dt = 0.12 * dx
n_steps_ct = 3000

for epsilon in [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50]:
    phi = phi_init.copy()
    dphi = cp.zeros_like(phi)

    # Modified force: standard + asymmetric potential
    def compute_force_asym(phi, eps):
        force = compute_force(phi, dx)
        # Add asymmetric potential: F_asym = 2*epsilon*phi on right side
        for comp in range(3):
            force[:,:,:,comp] += 2 * eps * phi[:,:,:,comp] * right_mask
        return force

    force = compute_force_asym(phi, epsilon)
    dphi += 0.5 * dt * force

    for step in range(n_steps_ct):
        phi += dt * dphi
        force = compute_force_asym(phi, epsilon)
        dphi += dt * force
        # Damping to reach equilibrium
        if step < n_steps_ct // 2:
            dphi *= 0.999

    L_f, R_f = measure_sides(phi)
    transfer = (R_f/(L_f+R_f) - R0/(L0+R0)) * 100  # % of total that shifted right
    lr_ratio = L_f / R_f if R_f > 0 else float('inf')
    print(f"  {epsilon:8.3f}  {L_f:8.4f}  {R_f:8.4f}  {lr_ratio:8.4f}  {transfer:+9.2f}%",
          flush=True)

cp.cuda.Stream.null.synchronize()


# =============================================================================
# TEST 4: Does amplitude nonlinearity depend on wave type?
# Compare s+s vs s+p at same R but varying amplitude
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 4: Nonlinearity comparison: s+s vs s+p")
print(f"  If nonlinear corrections differ by wave type,")
print(f"  the V6 formula's uniform treatment is wrong")
print(f"{'='*80}")

print(f"\n  {'amp':>5}  {'E_ss':>10}  {'E_ss/a2':>10}  {'E_sp':>10}  {'E_sp/a2':>10}  {'sp_dev%':>8}  {'ss_dev%':>8}")

E_ss_ref = None
E_sp_ref = None

for amp_val in [0.1, 0.2, 0.3, 0.5, 0.7, 1.0]:
    w_s1 = {'amplitude': amp_val, 'width': 2.0, 'n': 1, 'l': 0, 'channels': {0: 1.0}}
    w_s2 = {'amplitude': amp_val, 'width': 2.0, 'n': 1, 'l': 0, 'channels': {0: 1.0}}
    w_p = {'amplitude': amp_val, 'width': 2.0, 'n': 2, 'l': 1, 'channels': {0: 1.0}}

    t0 = timer.time()
    E_ss = run_pair(N, L, 4.0, w_s1, w_s2, n_steps)
    E_sp = run_pair(N, L, 4.0, w_s1, w_p, n_steps)
    dt = timer.time() - t0

    E_ss_n = E_ss / amp_val**2
    E_sp_n = E_sp / amp_val**2

    if E_ss_ref is None:
        E_ss_ref = E_ss_n
        E_sp_ref = E_sp_n

    ss_dev = (E_ss_n / E_ss_ref - 1) * 100
    sp_dev = (E_sp_n / E_sp_ref - 1) * 100

    print(f"  {amp_val:5.2f}  {E_ss:10.4f}  {E_ss_n:10.4f}  {E_sp:10.4f}  {E_sp_n:10.4f}  {sp_dev:+7.1f}%  {ss_dev:+7.1f}%  ({dt:.1f}s)",
          flush=True)


print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
