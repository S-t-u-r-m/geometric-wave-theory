"""
GWT sp Bond — Asymmetric Width Test
=====================================

Real sp bonds have asymmetric wave sizes:
  H(1s): width ~ a0/Z_H = a0 (large, diffuse)
  X(2p): width ~ 2*a0/Z_X (compact for O,F; wider for B,C)

The V6 formula ignores this size mismatch.
Does it affect the interaction energy?

If yes: this could explain why OH (+9.5%) overshoots more than HF (+0.5%),
since O(2p) is moderately compact while F(2p) is very compact.
"""

import cupy as cp
import numpy as np
import time as timer

pi = float(np.pi)
k_long = 1.0
k_trans = 0.5

def run_bond_sim(N_grid, L, R_sep, wave1, wave2, n_steps=2000, dt_factor=0.12):
    """Run 3-component 3D bond simulation on GPU."""
    dx = L / N_grid
    dt = dt_factor * dx

    coords = cp.linspace(-L/2, L/2, N_grid, dtype=cp.float32)
    X, Y, Z = cp.meshgrid(coords, coords, coords, indexing='ij')

    def make_wave(X, Y, Z, wave, center_z):
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
        return field

    def compute_energy(phi):
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

    def compute_force(phi):
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

    def evolve_and_measure(phi_init, n_steps):
        phi = phi_init.copy()
        dphi = cp.zeros_like(phi)
        force = compute_force(phi)
        dphi += 0.5 * dt * force
        E_samples = []
        for step in range(n_steps):
            phi += dt * dphi
            force = compute_force(phi)
            dphi += dt * force
            if step % 50 == 0 and step > n_steps // 4:
                dphi_half = dphi - 0.5 * dt * force
                E_kin = 0.5 * float(cp.sum(dphi_half**2)) * dx**3
                E_pot = compute_energy(phi)
                E_samples.append(E_kin + E_pot)
        return np.mean(E_samples) if E_samples else compute_energy(phi_init)

    phi1 = make_wave(X, Y, Z, wave1, -R_sep/2)
    E1 = evolve_and_measure(phi1, n_steps)

    phi2 = make_wave(X, Y, Z, wave2, +R_sep/2)
    E2 = evolve_and_measure(phi2, n_steps)

    phi_both = make_wave(X, Y, Z, wave1, -R_sep/2) + make_wave(X, Y, Z, wave2, +R_sep/2)
    E_both = evolve_and_measure(phi_both, n_steps)

    E_int = E_both - E1 - E2
    cp.cuda.Stream.null.synchronize()
    return E1, E2, E_both, E_int


# =============================================================================
print("=" * 80)
print("  SP BOND — ASYMMETRIC WIDTH TEST")
print("  Does wave size mismatch affect interaction energy?")
print("=" * 80)

_ = cp.zeros(10)
cp.cuda.Stream.null.synchronize()

N = 48
L = 16.0
n_steps = 2500
amp = 0.5

# =============================================================================
# TEST 1: Symmetric width reference
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 1: Reference — symmetric widths (both = 2.0)")
print(f"{'='*80}")

R_test = 4.0
w_ref = 2.0

w_s = {'amplitude': amp, 'width': w_ref, 'n': 1, 'l': 0, 'channels': {0: 1.0}}
w_p = {'amplitude': amp, 'width': w_ref, 'n': 2, 'l': 1, 'channels': {0: 1.0}}

t0 = timer.time()
_, _, _, E_ref = run_bond_sim(N, L, R_test, w_s, w_p, n_steps)
dt = timer.time() - t0
print(f"  E_int(symmetric, w=2.0) = {E_ref:.4f}  ({dt:.1f}s)")


# =============================================================================
# TEST 2: Vary p-wave width while keeping s-wave fixed
# Models different Z_eff for the p-atom
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 2: Fixed s-wave (w=2.0), varying p-wave width")
print(f"  Smaller p-width = higher Z_eff = more compact orbital")
print(f"  R = {R_test}")
print(f"{'='*80}")

print(f"\n  {'w_s':>4}  {'w_p':>4}  {'ratio':>6}  {'E_int':>10}  {'vs_ref':>8}")

for w_p_val in [1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5]:
    w_s_t = {'amplitude': amp, 'width': 2.0, 'n': 1, 'l': 0, 'channels': {0: 1.0}}
    w_p_t = {'amplitude': amp, 'width': w_p_val, 'n': 2, 'l': 1, 'channels': {0: 1.0}}
    t0 = timer.time()
    _, _, _, E_int = run_bond_sim(N, L, R_test, w_s_t, w_p_t, n_steps)
    dt = timer.time() - t0

    ratio = w_p_val / 2.0
    pct = (E_int - E_ref) / abs(E_ref) * 100 if abs(E_ref) > 1e-10 else 0
    print(f"  {2.0:4.1f}  {w_p_val:4.1f}  {ratio:6.2f}  {E_int:10.4f}  {pct:+7.1f}%  ({dt:.1f}s)",
          flush=True)


# =============================================================================
# TEST 3: Model real H-X molecules
# Map Z_eff to wave width: width ~ n / Z_eff (in some scale)
# Use w_H = scale/1.0, w_X = 2*scale/Z_eff
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 3: Model real H-X bonds with Z_eff-scaled widths")
print(f"  w_H = scale * 1/Z_H = scale (Z_H=1)")
print(f"  w_X = scale * n/Z_X")
print(f"{'='*80}")

Z_effs = {
    'B': (2.4214, 2.329, 3.42),
    'C': (3.1358, 2.116, 3.47),
    'N': (3.834,  1.958, 3.57),
    'O': (4.4532, 1.834, 4.392),
    'F': (5.0998, 1.733, 5.869),
}

# Try different scale factors to map Z_eff to sim width
for scale in [2.0, 3.0, 4.0]:
    print(f"\n  Scale = {scale}")
    print(f"  {'Atom':>4}  {'Z_eff':>5}  {'w_H':>5}  {'w_X':>5}  {'w_X/w_H':>7}  {'R':>5}  {'E_int':>10}")

    for atom, (Z, R_real, De_exp) in Z_effs.items():
        w_H = scale * 1.0    # H: n=1, Z=1
        w_X = scale * 2.0 / Z  # X(2p): n=2

        # Scale R proportionally
        R_sim = R_real * scale / 2.0  # map real R to sim R

        w_s_t = {'amplitude': amp, 'width': w_H, 'n': 1, 'l': 0, 'channels': {0: 1.0}}
        w_p_t = {'amplitude': amp, 'width': w_X, 'n': 2, 'l': 1, 'channels': {0: 1.0}}

        t0 = timer.time()
        _, _, _, E_int = run_bond_sim(N, L, R_sim, w_s_t, w_p_t, n_steps)
        dt = timer.time() - t0

        print(f"  {atom:>4}  {Z:5.2f}  {w_H:5.2f}  {w_X:5.2f}  {w_X/w_H:7.3f}  {R_sim:5.2f}  {E_int:10.4f}  ({dt:.1f}s)",
              flush=True)


# =============================================================================
# TEST 4: Direct comparison — same widths vs asymmetric widths
# Use a fixed R and compare symmetric vs asymmetric for each molecule
# =============================================================================
print(f"\n{'='*80}")
print(f"  TEST 4: Symmetric vs asymmetric width at R=4")
print(f"  Does width ratio affect bonding strength?")
print(f"{'='*80}")

print(f"\n  {'w_X/w_H':>7}  {'E_sym':>10}  {'E_asym':>10}  {'ratio':>8}")

w_H_base = 2.5  # H-like: broader

for w_ratio in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5]:
    w_X_val = w_H_base * w_ratio

    # Symmetric (both at geometric mean)
    w_geom = np.sqrt(w_H_base * w_X_val)
    w_sym_s = {'amplitude': amp, 'width': w_geom, 'n': 1, 'l': 0, 'channels': {0: 1.0}}
    w_sym_p = {'amplitude': amp, 'width': w_geom, 'n': 2, 'l': 1, 'channels': {0: 1.0}}

    # Asymmetric (actual widths)
    w_asym_s = {'amplitude': amp, 'width': w_H_base, 'n': 1, 'l': 0, 'channels': {0: 1.0}}
    w_asym_p = {'amplitude': amp, 'width': w_X_val, 'n': 2, 'l': 1, 'channels': {0: 1.0}}

    t0 = timer.time()
    _, _, _, E_sym = run_bond_sim(N, L, 4.0, w_sym_s, w_sym_p, n_steps)
    _, _, _, E_asym = run_bond_sim(N, L, 4.0, w_asym_s, w_asym_p, n_steps)
    dt = timer.time() - t0

    ratio = E_asym / E_sym if abs(E_sym) > 1e-6 else 0
    print(f"  {w_ratio:7.2f}  {E_sym:10.4f}  {E_asym:10.4f}  {ratio:8.4f}  ({dt:.1f}s)",
          flush=True)


print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
