"""
GWT Wave Simulator v2 — Mode Coupling from Cosine Potential
============================================================
Simulate the GWT Lagrangian:
  L = (1/2)(dphi)^2 + (1/pi^2)(1 - cos(pi*phi))

An atom = kink (charge Z) + breather modes (bound wave oscillations).

The simulation evolves the sine-Gordon FIELD on a 1D radial grid.
Mode-mode coupling happens automatically through the cosine potential.
No alpha formula, no parity rules — those should EMERGE.

Step 1: Solve for a single breather in a kink potential of charge Z_net.
         Measure binding energy vs Z_net. Extract the alpha relationship.
Step 2: Add multiple modes and let them self-consistently couple.
Step 3: IE = energy with N modes - energy with N-1 modes.
"""

import numpy as np
from math import factorial

# === GWT CONSTANTS (all from d=3 Lagrangian) ===
d = 3
gamma_sg = np.pi / (2**(d+1) * np.pi - 2)       # sine-Gordon coupling
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
V_0 = 1.0 / np.pi**2                              # potential depth
M_kink = 8.0 / np.pi**2                           # kink mass (BPS bound)
E_H_eV = (alpha_em**2 / 2) * 0.51100e6            # 13.604 eV

# Channel weights
w_pi = np.cos(np.pi / d)         # 0.5
w_delta = np.cos(2 * np.pi / d)  # -0.5

print("GWT Wave Simulator v2 — Sine-Gordon Field Dynamics")
print("=" * 65)
print(f"  L = (1/2)(dphi)^2 + (1/pi^2)(1 - cos(pi*phi))")
print(f"  gamma = {gamma_sg:.6f}, M_kink = {M_kink:.4f}")
print(f"  V_0 = 1/pi^2 = {V_0:.6f}")
print()

# === STEP 1: Single breather in kink potential ===
#
# The sine-Gordon equation on a 1D grid:
#   d^2(phi)/dt^2 = d^2(phi)/dx^2 - (1/pi) * sin(pi * phi)
#
# The kink solution: phi_kink(x) = (4/pi) * arctan(exp(x))
# The kink creates a potential well. A breather is a bound oscillation:
#   phi_breather(x,t) = (4/pi) * arctan(sin(omega*t) / (cosh(x) * epsilon))
# where omega = cos(gamma_sg) and epsilon = sin(gamma_sg)
#
# For an atom: the kink has "charge" Z, meaning the potential well is Z times deeper.
# The effective kink: phi_kink(x) = (4/pi) * arctan(exp(Z_eff * x))
# The breather binding energy scales with the well depth.

# Grid setup
Nx = 2000
x_max = 30.0
dx = 2 * x_max / Nx
x = np.linspace(-x_max, x_max, Nx)

# Time step (CFL condition: dt < dx)
dt = 0.5 * dx
N_steps = 20000  # enough for breather to complete several oscillations


def sine_gordon_kink(x, Z_eff=1.0):
    """Kink solution in potential well of charge Z_eff.
    Z protons = Z stacked kinks = potential well Z times deeper.
    Kink width shrinks as 1/sqrt(Z) (deeper well = steeper walls).
    """
    return (4.0 / np.pi) * np.arctan(np.exp(np.sqrt(Z_eff) * x))


# Global Z for potential scaling — set before each simulation
_Z_well = 1.0

def _update_Z(Z):
    global _Z_well
    _Z_well = Z

def sine_gordon_potential(phi):
    """V(phi) = Z * (1/pi^2)(1 - cos(pi*phi))
    Z kink-charges create a well Z times deeper.
    """
    return _Z_well * V_0 * (1.0 - np.cos(np.pi * phi))


def sine_gordon_force(phi):
    """dV/dphi = Z * (1/pi) * sin(pi*phi)"""
    return _Z_well * (1.0 / np.pi) * np.sin(np.pi * phi)


def add_breather(phi_bg, x, omega, x0=0.0, t0=0.0):
    """Add a breather perturbation to background field.

    omega = breather frequency (0 < omega < 1)
    Binding energy increases as omega decreases from 1.
    """
    epsilon = np.sqrt(1.0 - omega**2)
    # Breather solution (at t=t0):
    numerator = epsilon * np.cos(omega * t0)
    denominator = omega * np.cosh(epsilon * (x - x0))
    phi_br = (4.0 / np.pi) * np.arctan(numerator / (denominator + 1e-30))
    return phi_bg + phi_br


def evolve_sine_gordon(phi_init, N_steps, dt, dx):
    """Evolve sine-Gordon field using leapfrog (Verlet) integration.

    Returns time-averaged energy.
    """
    phi = phi_init.copy()
    phi_old = phi.copy()  # start with zero velocity

    energies = []

    for step in range(N_steps):
        # Laplacian (second derivative in space)
        laplacian = np.zeros_like(phi)
        laplacian[1:-1] = (phi[2:] - 2*phi[1:-1] + phi[:-2]) / dx**2

        # Force from cosine potential
        force = sine_gordon_force(phi)

        # Leapfrog update: phi_new = 2*phi - phi_old + dt^2 * (laplacian - force)
        phi_new = 2*phi - phi_old + dt**2 * (laplacian - force)

        # Boundary: fixed at kink asymptotes
        phi_new[0] = phi_init[0]
        phi_new[-1] = phi_init[-1]

        phi_old = phi.copy()
        phi = phi_new

        # Measure energy periodically
        if step > N_steps // 2 and step % 100 == 0:
            KE = 0.5 * np.sum(((phi - phi_old) / dt)**2) * dx
            GE = 0.5 * np.sum(((phi[1:] - phi[:-1]) / dx)**2) * dx
            PE = np.sum(sine_gordon_potential(phi)) * dx
            E_total = KE + GE + PE
            energies.append(E_total)

    return np.mean(energies) if energies else 0.0


def measure_binding_energy(Z_eff, omega=None):
    """Measure binding energy of a breather in a kink of charge Z_eff.

    Z_eff protons means:
    - Potential well is Z_eff times deeper
    - Kink walls are sqrt(Z_eff) times steeper
    - Breather binding energy scales with the well depth

    1. Compute energy of kink alone
    2. Compute energy of kink + breather
    3. Binding energy = difference
    """
    global _Z_well
    _Z_well = Z_eff  # Set potential depth

    if omega is None:
        # Breather frequency in scaled potential
        # Deeper well → lower frequency (more tightly bound)
        omega = np.cos(gamma_sg)

    # Kink alone
    phi_kink = sine_gordon_kink(x, Z_eff)
    E_kink = evolve_sine_gordon(phi_kink, N_steps, dt, dx)

    # Kink + breather
    phi_kb = add_breather(phi_kink, x, omega, x0=0.0)
    E_kink_breather = evolve_sine_gordon(phi_kb, N_steps, dt, dx)

    # Binding energy (should be negative — breather is bound)
    E_bind = E_kink_breather - E_kink

    return E_bind, E_kink, E_kink_breather


# === STEP 1: Single breather binding vs Z ===
print("Step 1: Single breather binding energy vs kink charge Z")
print("-" * 60)

z_values = [1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]
bind_1 = []

print(f"  {'Z':>5} {'E_bind':>10} {'E_bind/Z':>10} {'log fit':>10}")
print(f"  {'-'*45}")

for Z in z_values:
    E_bind, E_k, E_kb = measure_binding_energy(Z)
    bind_1.append(E_bind)
    print(f"  {Z:5.1f} {E_bind:10.4f} {E_bind/Z:10.4f}", flush=True)

# Fit power law
z_arr = np.array(z_values)
e_arr = np.array(bind_1)
log_z = np.log(z_arr)
log_e = np.log(np.abs(e_arr))
fit = np.polyfit(log_z, log_e, 1)
two_alpha_1D = fit[0]
alpha_1D = two_alpha_1D / 2

print(f"\n  1D fit: E ~ Z^{two_alpha_1D:.3f}  =>  alpha_1D = {alpha_1D:.4f}")
print(f"  3D correction: alpha_3D = alpha_1D * (d-1)/d = {alpha_1D * (d-1)/d:.4f}")
print(f"  v19 base: (d-1)/d^2 = {(d-1)/d**2:.4f}")


# === STEP 2: Two breathers — does pairing emerge? ===
print("\n" + "=" * 60)
print("Step 2: Two breathers in same kink — mode coupling")
print("-" * 60)
print("Adding a second breather. Same mode (s-pair) vs different mode.")
print()


def measure_two_breather_binding(Z_eff, omega1=None, omega2=None, x1=0.0, x2=0.0):
    """Measure energy of kink + 2 breathers.

    Compare:
    - E(kink + 2 breathers) vs E(kink + 1 breather) vs E(kink alone)
    - The SECOND breather's binding = E(k+2b) - E(k+1b)
    - If it differs from the first breather's binding, that's mode coupling.
    """
    global _Z_well
    _Z_well = Z_eff

    if omega1 is None:
        omega1 = np.cos(gamma_sg)
    if omega2 is None:
        omega2 = np.cos(2 * gamma_sg)  # second harmonic

    # Kink alone
    phi_k = sine_gordon_kink(x, Z_eff)
    E_k = evolve_sine_gordon(phi_k, N_steps, dt, dx)

    # Kink + first breather
    phi_k1 = add_breather(phi_k, x, omega1, x0=x1)
    E_k1 = evolve_sine_gordon(phi_k1, N_steps, dt, dx)

    # Kink + two breathers
    phi_k2 = add_breather(phi_k1, x, omega2, x0=x2)
    E_k2 = evolve_sine_gordon(phi_k2, N_steps, dt, dx)

    E_bind_1st = E_k1 - E_k    # first breather binding
    E_bind_2nd = E_k2 - E_k1   # second breather binding (includes coupling)
    E_coupling = E_bind_2nd - E_bind_1st  # pure coupling energy

    return E_bind_1st, E_bind_2nd, E_coupling, E_k, E_k1, E_k2


# Test: He-like system (Z=2, two breathers in same mode)
print("He-like: Z=2, two breathers in same kink well")
print("  Same frequency (s-pair):  both at omega = cos(gamma)")
E1, E2, E_coup, _, _, _ = measure_two_breather_binding(
    2.0, omega1=np.cos(gamma_sg), omega2=np.cos(gamma_sg), x1=0.0, x2=0.0)
print(f"    1st breather binding: {E1:.4f}")
print(f"    2nd breather binding: {E2:.4f}")
print(f"    Coupling energy:      {E_coup:.4f}")
print(f"    Ratio E2/E1:          {E2/E1:.4f}  (1.0 = no coupling)")
print()

# Test: same mode but offset (modeling different spin)
print("  Offset breathers (spin separation):  x1=0, x2=2")
E1, E2, E_coup, _, _, _ = measure_two_breather_binding(
    2.0, omega1=np.cos(gamma_sg), omega2=np.cos(gamma_sg), x1=-1.0, x2=1.0)
print(f"    1st breather binding: {E1:.4f}")
print(f"    2nd breather binding: {E2:.4f}")
print(f"    Coupling energy:      {E_coup:.4f}")
print(f"    Ratio E2/E1:          {E2/E1:.4f}")
print()

# Test: different harmonics (s + p like)
print("  Different harmonics (s + p like):  omega1=cos(g), omega2=cos(2g)")
E1, E2, E_coup, _, _, _ = measure_two_breather_binding(
    2.0, omega1=np.cos(gamma_sg), omega2=np.cos(2*gamma_sg), x1=0.0, x2=0.0)
print(f"    1st breather binding: {E1:.4f}")
print(f"    2nd breather binding: {E2:.4f}")
print(f"    Coupling energy:      {E_coup:.4f}")
print(f"    Ratio E2/E1:          {E2/E1:.4f}")
print()

# Sweep Z with two breathers — does the coupling scale geometrically?
print("Coupling vs Z (two same-frequency breathers):")
print(f"  {'Z':>5} {'E_1st':>8} {'E_2nd':>8} {'coupling':>10} {'ratio':>7}")
print(f"  {'-'*42}")

for Z in [1.0, 2.0, 3.0, 5.0, 10.0]:
    E1, E2, E_coup, _, _, _ = measure_two_breather_binding(
        Z, omega1=np.cos(gamma_sg), omega2=np.cos(gamma_sg))
    ratio = E2 / E1 if abs(E1) > 1e-10 else 0
    print(f"  {Z:5.1f} {E1:8.4f} {E2:8.4f} {E_coup:10.4f} {ratio:7.4f}", flush=True)

print()
print("If ratio != 1.0, the second breather is affected by the first.")
print("This IS mode coupling — the origin of parity corrections.")
print("Ratio > 1 means repulsion (anti-screening).")
print("Ratio < 1 means attraction (screening).")

# === STEP 3: Sweep breather harmonics — extract channel weights ===
print("\n" + "=" * 60)
print("Step 3: Mode coupling matrix — do channel weights emerge?")
print("-" * 60)
print("Testing all combinations of breather harmonics n1, n2.")
print("n=1 ~ s-channel, n=2 ~ p-channel, n=3 ~ d-channel")
print()

Z_test = 5.0  # use Z=5 for clear signal

# First, measure single-breather binding for each harmonic
print("Single breather binding by harmonic:")
single_bind = {}
for n_mode in [1, 2, 3, 4, 5]:
    omega = np.cos(n_mode * gamma_sg)
    if omega <= 0 or omega >= 1:
        break
    E_bind, _, _ = measure_binding_energy(Z_test, omega=omega)
    single_bind[n_mode] = E_bind
    print(f"  n={n_mode}: omega={omega:.6f}  E_bind={E_bind:.4f}", flush=True)

print()
print("Two-breather coupling matrix (ratio = E_2nd / E_1st_alone):")
print(f"  {'n1':>3} {'n2':>3} {'omega1':>8} {'omega2':>8} {'E_1st':>8} {'E_2nd':>8} {'ratio':>7} {'coupling':>10}")
print(f"  {'-'*65}")

for n1 in [1, 2, 3]:
    omega1 = np.cos(n1 * gamma_sg)
    if omega1 <= 0 or omega1 >= 1:
        continue
    for n2 in [1, 2, 3]:
        omega2 = np.cos(n2 * gamma_sg)
        if omega2 <= 0 or omega2 >= 1:
            continue

        _update_Z(Z_test)

        # Kink alone
        phi_k = sine_gordon_kink(x, Z_test)
        E_k = evolve_sine_gordon(phi_k, N_steps, dt, dx)

        # Kink + first breather (n1)
        phi_k1 = add_breather(phi_k, x, omega1, x0=0.0)
        E_k1 = evolve_sine_gordon(phi_k1, N_steps, dt, dx)

        # Kink + both breathers (n1 + n2)
        phi_k2 = add_breather(phi_k1, x, omega2, x0=0.0)
        E_k2 = evolve_sine_gordon(phi_k2, N_steps, dt, dx)

        E_1st = E_k1 - E_k
        E_2nd = E_k2 - E_k1
        ratio = E_2nd / E_1st if abs(E_1st) > 1e-10 else 0
        coupling = E_2nd - single_bind.get(n2, E_2nd)

        label = ""
        if n1 == n2:
            label = f" (same channel, expect ~d={d})"
        else:
            label = f" (cross channel)"

        print(f"  {n1:3d} {n2:3d} {omega1:8.5f} {omega2:8.5f} {E_1st:8.4f} {E_2nd:8.4f} {ratio:7.3f} {coupling:10.4f}{label}", flush=True)

# Extract the coupling matrix
print()
print("Summary — does the coupling matrix match GWT channel weights?")
print("  Same channel (n1=n2):  ratio ~ d = 3")
print("  Cross channel (n1!=n2): ratio ~ d+1 = 4")
print("  Difference: 1.0 = one extra coupling unit for cross-channel")
print("  w_pi = 0.5 would show as ratio = d + 0.5 for s+p...")
