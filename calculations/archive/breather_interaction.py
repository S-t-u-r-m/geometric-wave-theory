"""
Two-Breather Interaction on 1D Sine-Gordon Lattice
===================================================
GWT says atoms are breathers (localized oscillations) on the Planck lattice.
Bond energy = interaction energy of two breathers at separation R.

The sine-Gordon lattice Lagrangian:
  L = sum [ (1/2)(dphi/dt)^2 - (1/2)(phi_i - phi_{i+1})^2 - V_0*(1-cos(pi*phi)) ]

where V_0 = 1/pi^2 (from GWT).

A breather is a localized oscillating solution. Two breathers separated by R
interact through the lattice coupling. The interaction energy E(R) is:
  E(R) = E_total(two breathers at R) - 2*E_single

This naturally includes decay, nodes, phase, everything.

We compute this numerically by:
1. Setting up a breather profile (static approximation)
2. Placing two breathers at separation R
3. Computing the total energy
4. Comparing to 2x isolated breather energy
"""

import numpy as np
from scipy.optimize import minimize_scalar

pi = np.pi

# =============================================================================
# LATTICE PARAMETERS (from GWT Lagrangian, d=3)
# =============================================================================
V_0 = 1.0 / pi**2          # potential depth
M_kink = 8.0 / pi**2       # kink mass
gamma_breather = pi / (2**3 * pi)  # gamma = pi/(8pi) = 1/8 ... actually
# gamma = pi * V_0 / sqrt(1) -- need to check

# Breather parameter: for sine-Gordon, breather frequency depends on mode n:
# omega_n = sin(n * pi / (N_breathers + 1))
# where N_breathers = floor(8*pi - 1) = 24

N_br = 24  # number of breather modes


# =============================================================================
# SINE-GORDON BREATHER PROFILE (continuum limit)
# =============================================================================
# The static breather of the sine-Gordon equation is:
#   phi(x,t) = (4/pi) * arctan[ (sqrt(1-w^2)/w) * sin(w*t) / cosh(sqrt(1-w^2)*x) ]
# At t=0 (maximum displacement):
#   phi(x,0) = 0  (for the oscillating breather)
# At t = pi/(2w) (maximum):
#   phi(x, t_max) = (4/pi) * arctan[ (sqrt(1-w^2)/w) / cosh(sqrt(1-w^2)*x) ]
#
# The breather mass (energy) is:
#   M_n = (16/pi^2) * sin(n * gamma)
# where gamma = pi / (2N+2) for N+1 total modes... actually in GWT:
#   gamma = pi * sqrt(V_0) / 1 = 1/pi (for lattice spacing = 1)
#   M_n = 2 * M_kink * sin(n * gamma)

# Let's use gamma from the breather spectrum
gamma = pi / (N_br + 1)  # = pi/25

def breather_mass(n):
    """Mass of breather mode n (in Planck units)."""
    return 2 * M_kink * np.sin(n * gamma)


def breather_profile(x, n, amplitude=1.0):
    """Static breather profile for mode n.

    This is the maximum-displacement snapshot of a breather oscillation.
    The profile is a localized bump with width ~ 1/sqrt(1-w^2).

    For a sine-Gordon breather at maximum displacement:
      phi(x) = (4/pi) * arctan[ eta / cosh(eta * x) ]
    where eta = sqrt(1 - w_n^2), w_n = cos(n*gamma).
    """
    w_n = np.cos(n * gamma)
    eta = np.sin(n * gamma)  # sqrt(1 - w_n^2)

    # Prevent division by zero
    eta = max(eta, 1e-10)

    return amplitude * (4.0 / pi) * np.arctan(eta / np.cosh(eta * x))


def breather_width(n):
    """Characteristic width of breather mode n."""
    eta = np.sin(n * gamma)
    return 1.0 / max(eta, 0.01)


# =============================================================================
# LATTICE ENERGY COMPUTATION
# =============================================================================

def lattice_energy(phi, dx=1.0):
    """Compute total energy of field configuration phi on lattice.

    E = sum [ (1/2)(phi_i - phi_{i+1})^2 / dx + V_0*(1-cos(pi*phi_i))*dx ]

    where dx is the lattice spacing (= 1 in Planck units).
    """
    # Gradient energy (nearest-neighbor coupling)
    grad = np.diff(phi)
    E_grad = 0.5 * np.sum(grad**2) / dx

    # Potential energy
    E_pot = V_0 * np.sum(1 - np.cos(pi * phi)) * dx

    return E_grad + E_pot


def two_breather_field(x, n1, n2, R, a1=1.0, a2=1.0):
    """Create field with two breathers at positions -R/2 and +R/2."""
    phi1 = breather_profile(x + R/2, n1, a1)
    phi2 = breather_profile(x - R/2, n2, a2)
    return phi1 + phi2


def single_breather_energy(n, N_lattice=2000, L=100):
    """Energy of a single breather on the lattice."""
    x = np.linspace(-L/2, L/2, N_lattice)
    dx = x[1] - x[0]
    phi = np.array([breather_profile(xi, n) for xi in x])
    return lattice_energy(phi, dx)


def interaction_energy(n1, n2, R, N_lattice=2000, L=100):
    """Interaction energy of two breathers at separation R.

    E_int = E(two breathers at R) - E(breather 1 alone) - E(breather 2 alone)
    """
    x = np.linspace(-L/2, L/2, N_lattice)
    dx = x[1] - x[0]

    # Two breathers
    phi_both = np.array([two_breather_field(xi, n1, n2, R) for xi in x])
    E_both = lattice_energy(phi_both, dx)

    # Individual breathers
    phi1 = np.array([breather_profile(xi, n1) for xi in x])
    E1 = lattice_energy(phi1, dx)

    phi2 = np.array([breather_profile(xi, n2) for xi in x])
    E2 = lattice_energy(phi2, dx)

    return E_both - E1 - E2


# =============================================================================
# MAIN: Compute interaction energy vs separation
# =============================================================================

print("=" * 80)
print("  TWO-BREATHER INTERACTION ON SINE-GORDON LATTICE")
print("=" * 80)
print(f"\n  Lattice: V_0 = 1/pi^2 = {V_0:.6f}")
print(f"  Kink mass M_k = 8/pi^2 = {M_kink:.6f}")
print(f"  Breather gamma = pi/{N_br+1} = {gamma:.6f}")
print(f"  N_breathers = {N_br}")

# Show breather spectrum
print(f"\n  Breather masses (first 10 modes):")
for n in range(1, 11):
    w = breather_width(n)
    m = breather_mass(n)
    E_single = single_breather_energy(n)
    print(f"    n={n:2d}: M = {m:.6f}, width = {w:.2f}, E_lattice = {E_single:.6f}")


# =============================================================================
# INTERACTION ENERGY vs R for several mode pairs
# =============================================================================
print(f"\n{'='*80}")
print(f"  INTERACTION ENERGY vs SEPARATION")
print(f"{'='*80}")

mode_pairs = [
    (1, 1, "1s-1s (like H-H)"),
    (1, 2, "1s-2s (like H-Li)"),
    (2, 2, "2s-2s (like Li-Li)"),
    (1, 3, "1s-2p (like H-B/C/N/O/F)"),
    (3, 3, "2p-2p (like N-N, O-O)"),
    (3, 5, "2p-2p' (different Z, like B-F)"),
    (1, 4, "1s-3s (like H-Na)"),
]

R_values = np.arange(0.5, 15.1, 0.5)

for n1, n2, label in mode_pairs:
    print(f"\n  --- Mode ({n1},{n2}): {label} ---")
    print(f"  {'R':>6} {'E_int':>12} {'|E_int|':>10} {'sin(R)':>8} {'|sin|':>8}")

    E_vals = []
    for R in R_values:
        E_int = interaction_energy(n1, n2, R)
        E_vals.append(E_int)
        sin_R = np.sin(R)
        print(f"  {R:6.1f} {E_int:12.6f} {abs(E_int):10.6f} {sin_R:8.4f} {abs(sin_R):8.4f}")

    # Find the functional form: is E_int ~ exp(-R) * sin(kR)?
    E_arr = np.array(E_vals)
    R_arr = R_values

    # Find peaks and zeros
    peaks = []
    zeros = []
    for i in range(1, len(E_arr)-1):
        if abs(E_arr[i]) > abs(E_arr[i-1]) and abs(E_arr[i]) > abs(E_arr[i+1]):
            peaks.append((R_arr[i], E_arr[i]))
        if E_arr[i] * E_arr[i+1] < 0:
            zeros.append(R_arr[i])

    if peaks:
        print(f"  Peaks at: {', '.join(f'R={p[0]:.1f}(E={p[1]:.4f})' for p in peaks[:5])}")
    if zeros:
        print(f"  Zeros at: {', '.join(f'R={z:.1f}' for z in zeros[:5])}")

    # Fit envelope: |E_peak| vs R_peak
    if len(peaks) >= 2:
        R_peaks = np.array([p[0] for p in peaks])
        E_peaks = np.array([abs(p[1]) for p in peaks])
        if E_peaks[0] > 0 and E_peaks[-1] > 0:
            decay = np.log(E_peaks[0] / E_peaks[-1]) / (R_peaks[-1] - R_peaks[0])
            print(f"  Decay rate: {decay:.4f} (envelope ~ exp(-{decay:.3f}*R))")


# =============================================================================
# MAP BREATHER MODES TO ATOMS
# =============================================================================
print(f"\n{'='*80}")
print(f"  MAPPING: Breather modes to atomic orbitals")
print(f"{'='*80}")

print("""
  In GWT, each orbital type maps to a breather mode via quantum numbers:
    n_mode = n_principal (or some function of n, l)

  The question: which breather mode maps to which orbital?

  Atomic orbital  ->  Breather mode n  ->  Width (lattice units)
""")

orbital_map = [
    ('H_1s',  1),
    ('Li_2s', 2),
    ('B_2p',  3),
    ('C_2p',  3),
    ('N_2p',  3),
    ('O_2p',  3),
    ('F_2p',  3),
    ('Na_3s', 4),
    ('Cl_3p', 5),
]

for orb, n_mode in orbital_map:
    w = breather_width(n_mode)
    m = breather_mass(n_mode)
    print(f"  {orb:8s} -> mode {n_mode}: width = {w:.2f}, mass = {m:.6f}")


# =============================================================================
# COMPARE: Breather interaction to actual bond energies
# =============================================================================
print(f"\n{'='*80}")
print(f"  COMPARE: Breather interaction vs observed bond energies")
print(f"{'='*80}")

# Key molecules with their breather mode mapping
# R is in Bohr, but breather R is in lattice units (Planck lengths)
# Need a scale factor: R_lattice = R_bohr * (a_0 / l_P)
# a_0 = 0.529e-10 m, l_P = 1.616e-35 m
# R_lattice = R_bohr * 3.27e25 -- way too large!
#
# BUT: in GWT, the effective lattice for atomic physics isn't the Planck lattice.
# The breather IS the atom, and its width is in lattice units.
# The separation R should be measured in units of breather width.
#
# So: R_eff = R_bohr / (characteristic width)
# or simply: R_eff = R_bohr (treating bohr as the natural unit)

print(f"\n  Using R in Bohr as the natural length scale...")
print(f"  (The breather width sets the orbital size)")

molecules_test = [
    ('H2',   1.401,  4.745, 1, 1),
    ('LiH',  3.015,  2.515, 2, 1),
    ('BH',   2.329,  3.42,  3, 1),
    ('CH',   2.116,  3.47,  3, 1),
    ('NH',   1.958,  3.57,  3, 1),
    ('N2',   2.074,  9.759, 3, 3),
    ('HF',   1.733,  5.869, 3, 1),
    ('NaH',  3.566,  1.97,  4, 1),
    ('Li2',  5.051,  1.056, 2, 2),
    ('LiF',  2.955,  5.939, 2, 3),
    ('BF',   2.386,  7.81,  3, 3),
    ('NaCl', 4.461,  4.23,  4, 5),
    ('OH',   1.834,  4.392, 3, 1),
]

print(f"\n{'Mol':<7} {'R':>6} {'De_exp':>7} {'E_int':>10} {'|E_int|':>10} {'ratio':>8}")
print("-" * 55)

for name, R, De_exp, n1, n2 in molecules_test:
    E_int = interaction_energy(n1, n2, R)
    ratio = De_exp / abs(E_int) if abs(E_int) > 1e-8 else float('inf')
    print(f"{name:<7} {R:6.3f} {De_exp:7.3f} {E_int:10.6f} {abs(E_int):10.6f} {ratio:8.1f}")


# =============================================================================
# SCAN: What if R is in units of breather width?
# =============================================================================
print(f"\n{'='*80}")
print(f"  SCAN: R in units of breather width")
print(f"{'='*80}")

print(f"\n  R_scaled = R_bohr / avg_breather_width")

print(f"\n{'Mol':<7} {'R':>6} {'w1':>6} {'w2':>6} {'R/w_avg':>8} {'De_exp':>7} {'E_int':>10}")
print("-" * 60)

for name, R, De_exp, n1, n2 in molecules_test:
    w1 = breather_width(n1)
    w2 = breather_width(n2)
    w_avg = (w1 + w2) / 2
    R_scaled = R / w_avg
    E_int = interaction_energy(n1, n2, R_scaled)
    print(f"{name:<7} {R:6.3f} {w1:6.2f} {w2:6.2f} {R_scaled:8.3f} {De_exp:7.3f} {E_int:10.6f}")
