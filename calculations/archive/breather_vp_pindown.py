"""
Pin Down the VP Geometric Fraction
====================================
The VP/E ratio scales as alpha^2 × G, where G is the geometric fraction.
From the first run: G ≈ 0.25-0.30. GWT predicts G = 1/2^(d/2) = 0.354.

To pin it down:
  1. Vary lattice size N (16, 32, 64, 128, 256, 512) → extrapolate N→∞
  2. Exact sum over all k-modes (not average approximation)
  3. Second-order perturbation theory for the eigenvalue shift
  4. Compare 1D vs 3D to isolate the d-dependence
  5. Separate the geometric fraction from the numerical prefactor

The key formula:
  VP_shift = sum_k [delta_H(k)] / [2 * omega_k]
  where delta_H(k) = sum_i |psi_k(i)|^2 * [sin(pi*phi_x(i))/(pi*phi_x(i)) - 1]
  and psi_k(i) = (1/sqrt(N)) * exp(i*k*x_i) are the phonon eigenmodes

For a quasi-1D breather on a 3D lattice:
  The kx modes see the full perturbation.
  The ky,kz modes see a uniform perturbation (delta_H independent of y,z).
  So the sum factorizes.
"""
import sys, io, os, time
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)
alpha = np.exp(-(2/6) * (2**7/PI**2 + np.log(6)))

outfile = os.path.join(os.path.dirname(__file__), "breather_vp_pindown_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("PIN DOWN THE VP GEOMETRIC FRACTION")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"alpha = {alpha:.6f}, alpha^2 = {alpha**2:.6e}")
report(f"1/2^(d/2) = {1/2**(d/2):.6f} (GWT confined prediction)")
report(f"8/9 = {8/9:.6f} (GWT free prediction)")
report("")

# ============================================================
# EXACT 1D VP CALCULATION
# ============================================================
report("PART 1: 1D EXACT VP — N-DEPENDENCE")
report("-" * 60)
report("Exact sum over all k-modes, no approximations.")
report("")

def compute_vp_1d(N, n_mode):
    """Exact VP correction on 1D discrete lattice of N sites."""
    omega_n = np.cos(n_mode * gamma)
    eps_n = np.sin(n_mode * gamma)
    center = N // 2

    # Breather peak profile
    x = np.arange(N, dtype=np.float64) - center
    phi = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * x + 1e-30))

    # Hessian perturbation at each site
    abs_phi = np.abs(phi) + 1e-30
    delta_H = np.sin(PI * abs_phi) / (PI * abs_phi) - 1.0

    # Breather energy
    E_br = np.sum(phi**2) * omega_n / 2

    # Exact VP shift: sum over all k-modes
    # For each k: delta_omega_k = (1/N) * sum_i delta_H(i) * exp(0) / (2*omega_k)
    # Since |psi_k(i)|^2 = 1/N for plane waves and delta_H is real:
    # <k|delta_H|k> = (1/N) * sum_i delta_H(i)  [k-independent for uniform weights!]
    #
    # Wait — this assumes plane waves see the perturbation uniformly.
    # Actually: <k|delta_H|k> = (1/N) sum_i delta_H(i) ONLY if we ignore
    # the phase exp(ikx). The correct expression:
    # <k|delta_H|k> = (1/N) sum_i delta_H(i)
    # This IS correct because the perturbation is diagonal in position space
    # and the diagonal matrix elements of a position-space perturbation
    # in a plane-wave basis are just the spatial average (for each k separately).

    delta_H_avg = np.sum(delta_H) / N  # same for all k

    # Sum over k: delta_E = (1/2) * sum_k delta_H_avg / (2*omega_k) * N_modes
    # Wait, let me be more careful.
    #
    # The Hamiltonian for the y-component is H = T + V where
    # T = -Laplacian, V = scalar potential at each site
    # In the scalar model: V_i = 1 (for all i)
    # In the vector model: V_i = sin(pi*|phi_x_i|)/(pi*|phi_x_i|)
    # Perturbation: delta_V_i = V_vector_i - V_scalar_i = delta_H(i)
    #
    # The eigenvalues of H_scalar: omega_k^2 = 1 + 4*sin^2(k/2)
    # The first-order eigenvalue shift: delta(omega_k^2) = <k|delta_V|k> = (1/N)*sum delta_H
    # So: delta(omega_k) = delta(omega_k^2) / (2*omega_k) = delta_H_avg / (2*omega_k)
    #
    # ZPE shift = (1/2) * sum_k delta(omega_k) = (1/2) * delta_H_avg * sum_k 1/(2*omega_k)
    #           = (delta_H_avg / 4) * sum_k 1/omega_k

    k_arr = 2*PI * np.arange(N) / N
    omega_k = np.sqrt(1.0 + 4*np.sin(k_arr/2)**2)

    sum_inv_omega = np.sum(1.0 / omega_k)

    # VP shift for ONE transverse component
    delta_E_1comp = (delta_H_avg / 4.0) * sum_inv_omega

    # Total: 2 transverse components (y and z in d=3)
    delta_E_VP = (d - 1) * delta_E_1comp

    vp_over_E = delta_E_VP / E_br if E_br > 0 else 0
    G = vp_over_E / alpha**2

    return {
        'N': N, 'n': n_mode, 'E_br': E_br,
        'delta_E': delta_E_VP, 'vp_E': vp_over_E, 'G': G,
        'sum_dH': np.sum(delta_H), 'avg_dH': delta_H_avg,
        'sum_inv_omega': sum_inv_omega
    }

report(f"{'N':>6} {'n':>3} {'sum(dH)':>12} {'sum(1/w)':>12} "
       f"{'delta_E_VP':>14} {'VP/E':>12} {'G=VP/(E*a^2)':>14}")
report("-" * 85)

for n_mode in [1]:
    for N in [16, 32, 64, 128, 256, 512, 1024, 2048]:
        r = compute_vp_1d(N, n_mode)
        report(f"{r['N']:6d} {r['n']:3d} {r['sum_dH']:12.6f} {r['sum_inv_omega']:12.4f} "
               f"{r['delta_E']:+14.8f} {r['vp_E']:12.4e} {r['G']:14.4f}")

report("")

# Now check n-dependence at fixed large N
report("N-DEPENDENCE OF G (at N=2048):")
report(f"{'n':>3} {'G':>12} {'eps':>10} {'width':>8}")
report("-" * 36)
for n_mode in range(1, 9):
    r = compute_vp_1d(2048, n_mode)
    eps = np.sin(n_mode * gamma)
    width = 1.0 / eps if eps > 0.01 else 999
    report(f"{n_mode:3d} {r['G']:12.4f} {eps:10.6f} {width:8.1f}")

report("")

# ============================================================
# PART 2: PROPER 3D CALCULATION
# ============================================================
report("PART 2: 3D QUASI-1D VP — N-DEPENDENCE")
report("-" * 60)
report("")

def compute_vp_3d(N3, n_mode):
    """Exact VP on 3D lattice for quasi-1D breather."""
    omega_n = np.cos(n_mode * gamma)
    eps_n = np.sin(n_mode * gamma)
    center = N3 // 2

    x = np.arange(N3, dtype=np.float64) - center
    phi_1d = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * x + 1e-30))

    abs_phi = np.abs(phi_1d) + 1e-30
    delta_H_1d = np.sin(PI * abs_phi) / (PI * abs_phi) - 1.0

    # For quasi-1D: breather uniform in y,z
    # delta_H(i,j,k) = delta_H_1d(i)
    # <kx,ky,kz|delta_H|kx,ky,kz> = (1/N) sum_i delta_H_1d(i) × delta(ky,ky') × delta(kz,kz')
    # Wait — since delta_H is independent of y,z:
    # <k|delta_H|k> = (1/Nx) sum_i delta_H_1d(i)  [independent of ky, kz]

    delta_H_avg = np.sum(delta_H_1d) / N3  # average over x-direction only

    # 3D dispersion: omega_k^2 = 1 + 4*(sin^2(kx/2) + sin^2(ky/2) + sin^2(kz/2))
    k_arr = 2*PI * np.arange(N3) / N3

    # Sum over all 3D k-points: sum_{kx,ky,kz} 1/omega_k
    # Factorization fails because omega_k depends on all three k's.
    # Compute directly (N3^3 terms).
    sin2_arr = np.sin(k_arr/2)**2
    sum_inv_omega_3d = 0.0
    for ikx in range(N3):
        for iky in range(N3):
            omega_kz = np.sqrt(1.0 + 4*(sin2_arr[ikx] + sin2_arr[iky] + sin2_arr))
            sum_inv_omega_3d += np.sum(1.0 / omega_kz)

    # VP shift for one transverse component:
    # delta_E = (delta_H_avg / 4) * sum_{all 3D k} 1/omega_k / N3^2
    # The N3^2 comes from normalization: 3D eigenfunction norm = 1/N3^3,
    # but delta_H is uniform in y,z so ky,kz integrals give N3^2.
    # Net: delta_E = (1/4) * (sum_i delta_H / N3) * (sum_k 1/omega_k) / N3^2
    delta_E_1comp = (delta_H_avg / 4.0) * sum_inv_omega_3d / N3**2

    # 2 transverse components
    delta_E_VP = (d - 1) * delta_E_1comp

    # Breather energy (quasi-1D on 3D lattice: multiply by N3^2)
    E_br = N3**2 * np.sum(phi_1d**2) * omega_n / 2

    vp_over_E = delta_E_VP / E_br if E_br > 0 else 0
    G = vp_over_E / alpha**2

    return {'N3': N3, 'n': n_mode, 'delta_E': delta_E_VP, 'E_br': E_br,
            'vp_E': vp_over_E, 'G': G}

report(f"{'N':>4} {'n':>3} {'delta_E_VP':>14} {'E_br':>12} {'VP/E':>12} {'G':>12}")
report("-" * 64)

for n_mode in [1]:
    for N3 in [8, 12, 16, 24, 32, 48, 64]:
        r = compute_vp_3d(N3, n_mode)
        report(f"{r['N3']:4d} {r['n']:3d} {r['delta_E']:+14.8f} {r['E_br']:12.4f} "
               f"{r['vp_E']:12.4e} {r['G']:12.4f}")

report("")

# n-dependence at N=48
report("N-DEPENDENCE OF G IN 3D (at N=48):")
report(f"{'n':>3} {'G_3d':>12} {'G_1d':>12} {'ratio':>10}")
report("-" * 40)
for n_mode in range(1, 9):
    r3 = compute_vp_3d(48, n_mode)
    r1 = compute_vp_1d(2048, n_mode)
    ratio = r3['G'] / r1['G'] if r1['G'] != 0 else 0
    report(f"{n_mode:3d} {r3['G']:12.4f} {r1['G']:12.4f} {ratio:10.4f}")

report("")

# ============================================================
# PART 3: INTERPRETATION
# ============================================================
report("PART 3: INTERPRETATION")
report("-" * 60)
report("")

# The geometric fraction G should relate to Oh channel decomposition
r1_conv = compute_vp_1d(2048, 1)
r3_conv = compute_vp_3d(48, 1)

report(f"Converged values (n=1):")
report(f"  G_1d (N=2048) = {r1_conv['G']:.6f}")
report(f"  G_3d (N=48)   = {r3_conv['G']:.6f}")
report("")
report(f"GWT predictions:")
report(f"  1/2^(d/2) = {1/2**(d/2):.6f}  (confined VP, mass ratio)")
report(f"  8/9       = {8/9:.6f}  (free VP, alpha dressing)")
report(f"  8/3       = {8/3:.6f}  (gluon VP, alpha_s dressing)")
report(f"  (d-1)/6   = {(d-1)/6:.6f}  (from -(pi*phi)^2/6, 2 components)")
report(f"  1/d       = {1/d:.6f}")
report(f"  1/d^2     = {1/d**2:.6f}")
report("")

# The leading-order analytic result:
# delta_H ≈ -(pi*phi_x)^2/6
# sum delta_H ≈ -(pi^2/6) * sum phi_x^2
# E_br = (omega/2) * sum phi_x^2
# VP/E = -(pi^2/6) * (d-1) * <1/omega> / (4 * omega * N_eff)
# where N_eff depends on lattice size

# Direct calculation of the leading-order coefficient:
x = np.arange(2048, dtype=np.float64) - 1024
eps_1 = np.sin(gamma)
omega_1 = np.cos(gamma)
phi = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_1 * x + 1e-30))

sum_phi2 = np.sum(phi**2)
sum_phi4 = np.sum(phi**4)
sum_dH = np.sum(np.sin(PI*np.abs(phi)+1e-30)/(PI*np.abs(phi)+1e-30) - 1)

# Ratio of leading-order term to exact:
leading_order = -PI**2/6 * sum_phi2
report(f"Leading-order check:")
report(f"  sum(delta_H) exact    = {sum_dH:.6f}")
report(f"  sum(delta_H) leading  = {leading_order:.6f}  (-pi^2/6 * sum(phi^2))")
report(f"  Ratio: {sum_dH/leading_order:.6f}  (1.0 = perfect leading-order)")
report("")

# The phi^4 correction:
next_order = PI**4/120 * sum_phi4
report(f"  Next-order correction: {next_order:.6f}  (+pi^4/120 * sum(phi^4))")
report(f"  sum(delta_H) ≈ leading + next = {leading_order + next_order:.6f}")
report(f"  Ratio with correction: {sum_dH/(leading_order+next_order):.6f}")
report("")

report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))

log.close()
print(f"\nResults saved to: {outfile}")
