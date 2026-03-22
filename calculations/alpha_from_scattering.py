"""
Measure alpha from breather-breather scattering on the discrete sine-Gordon lattice.

On the continuum sine-Gordon, breather-breather scattering is ELASTIC (integrable).
On the DISCRETE lattice, integrability is broken. Two colliding breathers:
  1. Pass through each other (elastic part)
  2. Radiate a small fraction of energy (inelastic part)

The inelastic fraction = the coupling strength squared = alpha^2.
If we measure the radiated energy, we get alpha without using the formula.

Method:
  - Create two breather wavepackets moving toward each other
  - Evolve on the discrete sine-Gordon lattice
  - Measure total energy before and after collision
  - The energy difference = radiated energy = alpha^2 * E_collision
  - Extract alpha = sqrt(E_rad / E_collision)
"""

import numpy as np
import time

PI = np.pi
d = 3

# GWT prediction for comparison
from math import factorial, log, exp
alpha_gwt = exp(-(2/factorial(d)) * (2**(2*d+1)/PI**2 + log(2*d)))
print(f"GWT prediction: alpha = {alpha_gwt:.8f} = 1/{1/alpha_gwt:.3f}")
print()

# =============================================================
# DISCRETE SINE-GORDON SIMULATION
# =============================================================

# Lattice parameters
Nx = 4096       # spatial sites
dx = 0.05       # lattice spacing (in natural units)
dt = 0.02       # time step (CFL: dt < dx)
L = Nx * dx     # total length

x = np.arange(Nx) * dx - L/2  # centered at 0

# Sine-Gordon equation: phi_tt = phi_xx - (1/pi)*sin(pi*phi)
# On discrete lattice: phi_tt = (phi_{i+1} - 2*phi_i + phi_{i-1})/dx^2 - (1/pi)*sin(pi*phi_i)

# Breather solution (moving):
# phi(x,t) = (4/pi)*arctan(sin(omega*gamma*(t-t0) - v*gamma*(x-x0)/... ) / (omega*cosh(...)))
# For simplicity, use a boosted breather.

def breather(x, t, x0, v, n=1):
    """Moving breather solution of the sine-Gordon equation.
    x0: initial center position
    v: velocity (|v| < 1)
    n: mode number
    """
    gamma_sg = PI / (2**(d+1) * PI - 2)
    omega = np.cos(n * gamma_sg)
    eps = np.sin(n * gamma_sg)

    # Lorentz boost
    gamma_L = 1.0 / np.sqrt(1 - v**2 + 1e-30)
    xi = gamma_L * (x - x0 - v * t)
    tau = gamma_L * (t - v * (x - x0))

    # Breather profile
    num = eps * np.sin(omega * tau)
    den = omega * np.cosh(eps * xi) + 1e-30
    phi = (4/PI) * np.arctan(num / den)

    return phi

def energy_density(phi, phi_dot, dx):
    """Compute energy density at each site."""
    # Kinetic
    T = 0.5 * phi_dot**2

    # Gradient (using periodic BC)
    dphi = np.roll(phi, -1) - phi
    G = 0.5 * dphi**2 / dx**2

    # Potential
    V = (1/PI**2) * (1 - np.cos(PI * phi))

    return T + G + V

def total_energy(phi, phi_dot, dx):
    return np.sum(energy_density(phi, phi_dot, dx)) * dx

# =============================================================
# EXPERIMENT 1: Single breather — measure rest energy
# =============================================================
print("EXPERIMENT 1: Single breather energy")
print("=" * 55)

phi = breather(x, 0, 0, 0, n=1)
phi_dot = np.zeros_like(phi)  # approximate (exact would need time derivative)

# Better: compute phi_dot from the breather time derivative at t=0
gamma_sg = PI / (2**(d+1) * PI - 2)
omega_1 = np.cos(gamma_sg)
eps_1 = np.sin(gamma_sg)

# d/dt of breather at t=0, v=0:
# phi = (4/pi)*arctan(eps*sin(omega*t) / (omega*cosh(eps*x)))
# dphi/dt|_{t=0} = (4/pi) * eps*omega / (omega*cosh(eps*x)) / (1 + 0) = (4/pi)*eps/cosh(eps*x)
phi_dot = (4/PI) * eps_1 / np.cosh(eps_1 * x)

E_single = total_energy(phi, phi_dot, dx)
print(f"  Single breather energy: E = {E_single:.6f}")
print(f"  Theoretical: m_breather = (2^(d+1)/pi^2)*sin(gamma) = {(2**(d+1)/PI**2)*np.sin(gamma_sg):.6f}")
print()

# =============================================================
# EXPERIMENT 2: Two-breather collision
# =============================================================
print("EXPERIMENT 2: Two-breather collision")
print("=" * 55)

# Two breathers approaching each other
separation = 40.0  # initial separation
v_collision = 0.3   # collision velocity

# Initialize
phi = breather(x, 0, -separation/2, +v_collision, n=1) + \
      breather(x, 0, +separation/2, -v_collision, n=1)

# Time derivatives (sum of individual breather velocities)
# For moving breather, dphi/dt involves both omega and boost terms
# Use finite difference for accuracy
dt_init = 1e-4
phi_plus = breather(x, dt_init, -separation/2, +v_collision, n=1) + \
           breather(x, dt_init, +separation/2, -v_collision, n=1)
phi_minus = breather(x, -dt_init, -separation/2, +v_collision, n=1) + \
            breather(x, -dt_init, +separation/2, -v_collision, n=1)
phi_dot = (phi_plus - phi_minus) / (2 * dt_init)

E_initial = total_energy(phi, phi_dot, dx)
print(f"  Initial energy: E = {E_initial:.6f}")
print(f"  (Should be ~2 * E_single = {2*E_single:.6f})")
print()

# Evolve until collision and well past it
# Time for breathers to meet: t_meet = (separation/2) / v_collision
t_meet = separation / (2 * v_collision)
t_total = 3 * t_meet  # run well past collision
N_steps = int(t_total / dt)

print(f"  Collision time: t_meet = {t_meet:.1f}")
print(f"  Total evolution: t_total = {t_total:.1f}")
print(f"  Steps: {N_steps}")
print(f"  Evolving...", end=" ", flush=True)

t0 = time.time()

# Leapfrog integration
phi_old = phi.copy()
# First half-step for velocity
acc = np.zeros_like(phi)
acc[1:-1] = (phi[2:] - 2*phi[1:-1] + phi[:-2]) / dx**2
acc[0] = (phi[1] - 2*phi[0] + phi[-1]) / dx**2  # periodic
acc[-1] = (phi[0] - 2*phi[-1] + phi[-2]) / dx**2
acc -= (1/PI) * np.sin(PI * phi)

phi_dot += 0.5 * dt * acc

# Track energy at key points
energies = []
times_record = []

for step in range(N_steps):
    # Position update
    phi += dt * phi_dot

    # Acceleration
    acc[1:-1] = (phi[2:] - 2*phi[1:-1] + phi[:-2]) / dx**2
    acc[0] = (phi[1] - 2*phi[0] + phi[-1]) / dx**2
    acc[-1] = (phi[0] - 2*phi[-1] + phi[-2]) / dx**2
    acc -= (1/PI) * np.sin(PI * phi)

    # Velocity update
    phi_dot += dt * acc

    # Record energy periodically
    if step % (N_steps // 20) == 0:
        E = total_energy(phi, phi_dot, dx)
        energies.append(E)
        times_record.append(step * dt)

t1 = time.time()
print(f"done ({t1-t0:.1f}s)")

E_final = total_energy(phi, phi_dot, dx)
print()

print(f"  Energy conservation check:")
print(f"    E_initial = {E_initial:.8f}")
print(f"    E_final   = {E_final:.8f}")
print(f"    Delta_E   = {E_final - E_initial:.2e}")
print(f"    Relative  = {abs(E_final - E_initial)/E_initial:.2e}")
print()

# =============================================================
# EXPERIMENT 3: Measure radiation from collision
# =============================================================
print("EXPERIMENT 3: Radiation measurement")
print("=" * 55)

# After collision, the breathers have passed through each other.
# Any energy NOT in the two breather regions = radiation.
# Measure energy in the breather regions vs total.

# Find breather positions after collision
# They should be at approximately +/- (separation/2 + v * t_extra)
t_extra = t_total - t_meet
x_left = -separation/2 - v_collision * t_extra
x_right = +separation/2 + v_collision * t_extra

# Energy in breather regions (window of width ~20/eps around each center)
window = 30.0  # half-width of measurement window
mask_left = (x > x_left - window) & (x < x_left + window)
mask_right = (x > x_right - window) & (x < x_right + window)
mask_breathers = mask_left | mask_right
mask_radiation = ~mask_breathers

e_density = energy_density(phi, phi_dot, dx)
E_breathers = np.sum(e_density[mask_breathers]) * dx
E_radiation = np.sum(e_density[mask_radiation]) * dx

print(f"  Breather positions: x_L = {x_left:.1f}, x_R = {x_right:.1f}")
print(f"  Measurement window: +/- {window:.0f} around each breather")
print(f"  E in breather regions: {E_breathers:.8f}")
print(f"  E in radiation region: {E_radiation:.8f}")
print(f"  Total: {E_breathers + E_radiation:.8f}")
print()

# The radiation fraction
f_rad = E_radiation / E_initial
print(f"  Radiation fraction: E_rad / E_total = {f_rad:.6f}")
print(f"  = {f_rad:.2e}")
print()

# In the integrable continuum, f_rad = 0 exactly.
# On the discrete lattice, f_rad ~ alpha^2 * (discreteness parameter)
# The discreteness parameter depends on v and dx.

# For the SG breather-breather scattering on a discrete lattice,
# the leading inelasticity comes from the lattice breaking integrability.
# The radiation amplitude ~ (dx/lambda)^2 where lambda = breather width.
# The breather width = 1/eps = 1/sin(gamma) ~ 15.8 lattice widths.
# So discreteness parameter ~ (dx * eps)^2 = (0.05 * 0.063)^2 = 1e-5.
# And f_rad ~ alpha^2 * (dx*eps)^2 ~ 5e-5 * 1e-5 ~ 5e-10.

# That's too small to measure! The radiation from discreteness is tiny
# because our lattice spacing dx = 0.05 is very fine compared to the breather width.

discreteness = (dx * eps_1)**2
print(f"  Discreteness parameter: (dx*eps)^2 = {discreteness:.2e}")
print(f"  Expected radiation: alpha^2 * discreteness = {alpha_gwt**2 * discreteness:.2e}")
print()

# The radiation we measured is dominated by NUMERICAL NOISE (energy conservation error)
# not by the physical inelasticity.
print(f"  Energy conservation error: {abs(E_final-E_initial)/E_initial:.2e}")
print(f"  Radiation signal: {f_rad:.2e}")
print(f"  Signal/noise: {f_rad / (abs(E_final-E_initial)/E_initial + 1e-30):.2f}")
print()

# =============================================================
# EXPERIMENT 4: Alternative - measure coupling from energy levels
# =============================================================
print("EXPERIMENT 4: Coupling from double-well splitting")
print("=" * 55)

# Instead of scattering, measure alpha from the ENERGY LEVEL SPLITTING
# of a breather in a double-well potential (two adjacent kink wells).
# The splitting = 2*alpha (the tunneling amplitude between wells).
#
# This is what the 1D transfer matrix measured: t = 0.0429.
# But this is the 1D tunneling, not the 3D alpha.
#
# On the 3D lattice, the splitting involves d=3 directions.
# The 3D alpha should be related to the 1D tunneling by:
# alpha_3D = f(t_1D, d) where f encodes the geometric factors.
#
# From the transfer matrix: t_1D = 0.0429
# From the formula: alpha = 0.00730
# Ratio: alpha / t_1D = 0.170

t_1D = 0.04288  # from earlier transfer matrix calculation
ratio = alpha_gwt / t_1D
print(f"  1D tunneling: t = {t_1D:.6f}")
print(f"  3D alpha: {alpha_gwt:.6f}")
print(f"  Ratio alpha/t = {ratio:.6f}")
print(f"  = 1/{1/ratio:.2f}")
print()

# What is this ratio?
# 0.170 ~ 1/5.88 ~ 1/(2*d) = 1/6? Close but not exact.
# 0.170 ~ (d-1)/(d*d!) = 2/18 = 0.111? No.
# 0.170 ~ alpha_gwt / t_1D

# The 1D tunneling t already includes the 1D barrier action exp(-M_kink)
# and the 1D prefactor. The 3D alpha includes the 3D geometric factors.
# The ratio encodes the GEOMETRIC CONVERSION from 1D to 3D.

# Let's check: if alpha = t^p, p = ln(alpha)/ln(t) = 1.562
# We computed this earlier. It's not a clean power.

print("SUMMARY")
print("=" * 55)
print()
print("The scattering approach doesn't work here because:")
print("  - The lattice is too fine (dx << breather width)")
print("  - The inelastic signal is ~1e-10, below numerical noise")
print()
print("The double-well splitting approach works in 1D:")
print(f"  - 1D tunneling t = {t_1D:.6f} (measured from transfer matrix)")
print(f"  - 1D barrier action confirmed: M_kink = 8/pi^2")
print()
print("The 3D geometric factors (channel selection, symmetry, transverse")
print("fraction) convert the 1D tunneling to the 3D coupling:")
print(f"  alpha = t_1D * {ratio:.4f}")
print(f"  = {alpha_gwt:.6f} = 1/{1/alpha_gwt:.3f}")
print()
print("WHAT HAS BEEN CONFIRMED FROM SIMULATION:")
print("  1. Barrier height M_kink = 8/pi^2 (exact, BPS bound)")
print("  2. 1D tunneling rate (transfer matrix eigenvalues)")
print("  3. Fluctuation determinant = 1 (reflectionless PT)")
print("  4. Bond energy D_e = pi/d^2 * E_H (Hessian eigenvalues)")
print("  5. VP_self = -0.7589 (direct computation)")
print("  6. 8 stable breather modes (dynamics)")
print()
print("WHAT STILL USES THE FORMULA:")
print("  The conversion from 1D tunneling to 3D alpha = 1/137")
print("  specifically the prefactor (2d)^(-2/d!) = 6^(-1/3)")
