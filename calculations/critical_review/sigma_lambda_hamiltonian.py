#!/usr/bin/env python3
"""
Sigma-Lambda from a coupled-oscillator Hamiltonian.

Pre-registered approach (NOT fitting to 77 MeV):
  - Model the uds baryon as 3 coupled oscillators on a triangular geometry
  - 2 "light" modes (frequency w_l)
  - 1 "heavy" mode for strange (frequency w_h > w_l)
  - Symmetric coupling J between all pairs

The two stable configurations:
  - Sigma: light pair oscillating IN PHASE with each other (symmetric)
  - Lambda: light pair oscillating OUT OF PHASE (antisymmetric)

These are the two eigenstates of the same underlying torus structure.

Parameters chosen from framework primitives (NOT from target 77 MeV):
  - omega_0 = 294 MeV (light angular quantum, framework-derived)
  - w_l = omega_0 / d = 98 MeV (light mode frequency, per axis)
  - w_h = w_l * (m_strange/m_light scaling factor)
  - J = alpha_s * omega_0 (strong-coupling-mediated coupling)

We then DIAGONALIZE the Hamiltonian to find the two normal mode
configurations and their energies. The gap between symmetric and
antisymmetric eigenmodes = our prediction for Sigma-Lambda.

This is the proper version of last night's exploration. We compute,
then compare. No reverse-engineering.
"""
import numpy as np

PI = np.pi
d = 3
alpha_s = 0.11794  # framework strong coupling
m_p = 938.272
omega_0 = m_p * (2*d*PI) / (2*d*PI**2 + 1)  # light angular quantum

# Pre-registered parameter choices (framework-motivated, not fit to target)
w_l = omega_0 / d                 # 98 MeV per axis for light mode
                                  # (this is the one of d sub-circulations)
mass_ratio = 1.5                  # heavy/light frequency ratio
                                  # (justification: strange has higher
                                  #  effective mass; 1.5x is a placeholder
                                  #  not fit to target. Real value should
                                  #  derive from breather mass formulas.)
w_h = w_l * mass_ratio            # 147 MeV for "strange" mode
J   = alpha_s * omega_0           # 34.7 MeV coupling

print('=' * 70)
print('SIGMA-LAMBDA HAMILTONIAN ATTEMPT')
print('=' * 70)
print()
print(f'Pre-registered framework parameters:')
print(f'  omega_0 = {omega_0:.2f} MeV  (light angular quantum)')
print(f'  w_l     = omega_0/d = {w_l:.2f} MeV  (light mode)')
print(f'  w_h     = w_l * {mass_ratio} = {w_h:.2f} MeV  (heavy mode)')
print(f'  J       = alpha_s * omega_0 = {J:.2f} MeV  (coupling)')
print()
print('NOTE: mass_ratio = 1.5 is a placeholder. A real derivation would')
print('use the framework breather mass formulas for u/d vs s modes.')
print('We commit to this number BEFORE looking at the result.')
print()

# Build the Hamiltonian for 3 coupled harmonic oscillators
# H = sum_i (1/2)(p_i^2 + w_i^2 x_i^2) + J * sum_{i<j} x_i x_j
#
# In the basis (x_1, x_2, x_3) where 1,2 are light and 3 is heavy,
# the potential matrix is:
#   V = diag(w_l^2, w_l^2, w_h^2) + J * (off-diag couplings)
#
# Coupling matrix (off-diagonal): symmetric, J everywhere
# So V_ij = w_i^2 * delta_ij + J * (1 - delta_ij)

V = np.array([
    [w_l**2, J,      J     ],
    [J,      w_l**2, J     ],
    [J,      J,      w_h**2],
])

print('Potential matrix V (in MeV^2):')
for row in V:
    print('  ', '  '.join(f'{x:>10.2f}' for x in row))
print()

# Diagonalize: eigenfrequencies are sqrt of eigenvalues of V
# Eigenvectors tell us the normal mode structure
eigvals, eigvecs = np.linalg.eigh(V)
eigfreqs = np.sqrt(eigvals)

print('Normal mode analysis:')
print(f'  {"Mode":<6} {"Frequency (MeV)":>18} {"Eigenvector (x1,x2,x3)":>30}')
for i in range(3):
    vec_str = '(' + ', '.join(f'{v:+.3f}' for v in eigvecs[:, i]) + ')'
    print(f'  {i+1:<6} {eigfreqs[i]:>18.2f} {vec_str:>30}')
print()

# Identify which mode is "Lambda-like" (antisymmetric in light pair)
# and which is "Sigma-like" (symmetric in light pair)
# By inspection:
#   - Antisymmetric mode: light components have opposite signs, heavy ~ 0
#   - Symmetric modes: light components same sign

antisym_idx = None
sym_with_heavy = None
sym_against_heavy = None

for i in range(3):
    v = eigvecs[:, i]
    # antisymmetric: x1 = -x2, |x3| small
    if abs(v[0] + v[1]) < 0.01 and abs(v[2]) < 0.1:
        antisym_idx = i
    # symmetric: x1 = +x2
    elif abs(v[0] - v[1]) < 0.01:
        if v[0] * v[2] > 0:
            sym_with_heavy = i
        else:
            sym_against_heavy = i

print('Mode identification:')
if antisym_idx is not None:
    print(f'  Antisymmetric (Lambda candidate): mode {antisym_idx+1}, '
          f'frequency = {eigfreqs[antisym_idx]:.2f} MeV')
if sym_with_heavy is not None:
    print(f'  Symmetric (with heavy):           mode {sym_with_heavy+1}, '
          f'frequency = {eigfreqs[sym_with_heavy]:.2f} MeV')
if sym_against_heavy is not None:
    print(f'  Symmetric (against heavy):        mode {sym_against_heavy+1}, '
          f'frequency = {eigfreqs[sym_against_heavy]:.2f} MeV')
print()

# Compute the Sigma-Lambda gap
# In a quantum harmonic oscillator interpretation:
# - Each mode contributes (n + 1/2) * hbar * omega to total energy
# - Total mass ~ sum of mode energies in lowest state
# - Energy gap between two configurations = difference in zero-point energy
#   if we use symmetric vs antisymmetric for the LIGHT pair

# Simplest interpretation:
# - Lambda config = antisymmetric light pair: uses the antisymmetric mode
# - Sigma config = symmetric light pair: uses the symmetric mode (with heavy)
# - Gap = difference in mode frequencies (times some occupation factor)

if antisym_idx is not None and sym_with_heavy is not None:
    f_lambda = eigfreqs[antisym_idx]
    f_sigma  = eigfreqs[sym_with_heavy]
    gap = f_sigma - f_lambda

    print('='*70)
    print('PREDICTION')
    print('='*70)
    print()
    print(f'  Lambda config (antisymmetric light pair): omega = {f_lambda:.2f} MeV')
    print(f'  Sigma config (symmetric light pair):       omega = {f_sigma:.2f} MeV')
    print(f'  Energy gap (Sigma - Lambda):               {gap:.2f} MeV')
    print()
    print(f'  Observed:                                  77.0 MeV')
    print(f'  Predicted:                                 {gap:.2f} MeV')
    err = (gap - 77.0) / 77.0 * 100
    print(f'  Residual:                                  {err:+.1f}%')
    print()

print('='*70)
print('HONEST ASSESSMENT')
print('='*70)
print()
print('What this calculation IS:')
print('  - A toy 3-coupled-oscillator model with framework parameters')
print('  - Honest diagonalization (no parameter tuning)')
print('  - First-cut estimate of energy gap structure')
print()
print('What this calculation IS NOT:')
print('  - A real lattice sine-Gordon calculation')
print('  - A derivation that w_l = omega_0/d (chosen, not derived)')
print('  - A derivation that mass_ratio = 1.5 (placeholder)')
print('  - A derivation that J = alpha_s * omega_0 (motivated, not proven)')
print()
print('To be a real derivation, we would need:')
print('  - Lattice sine-Gordon Hamiltonian for 3 coupled breather modes')
print('  - Mass parameters derived from framework breather formulas')
print('  - Coupling derived from torus geometry, not assumed')
print('  - All from FIRST PRINCIPLES, not framework-natural-but-chosen')
print()
print('Result (whatever it is) tells us if the STRUCTURE of the calc')
print('produces reasonable scales, not whether the formula is right.')
