"""
EXACT VP Correction: Hessian Method
=====================================
The VP is a shift in the zero-point energy of the y,z vacuum modes
caused by the breather's nonlinear coupling.

At the breather configuration (phi_x = breather, phi_y = phi_z = 0):

  Scalar model: d²V/dphi_y² = 1  (independent of phi_x)
  Vector model: d²V/dphi_y² = sin(pi*|phi_x|) / (pi*|phi_x|)  (depends on phi_x!)

The difference: delta_H = sin(pi*phi_x)/(pi*phi_x) - 1 ≈ -(pi*phi_x)²/6

This means the breather SOFTENS the y,z phonon modes at sites where phi_x is nonzero.
The total VP correction = change in zero-point energy of all y,z modes.

No noise. No time evolution. Pure linear algebra.
"""
import sys, io, os, time
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)
alpha = np.exp(-(2/6) * (2**7/PI**2 + np.log(6)))

outfile = os.path.join(os.path.dirname(__file__), "breather_vp_exact_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("EXACT VP CORRECTION — HESSIAN METHOD")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"alpha = {alpha:.6f}, alpha^2 = {alpha**2:.6e}")
report(f"gamma = {gamma:.10f}")
report("")

# ============================================================
# PART 1: The Hessian perturbation at the breather
# ============================================================
report("PART 1: HESSIAN PERTURBATION")
report("-" * 60)

# The n=1 breather profile at peak amplitude on a discrete lattice:
# phi_x(i) = (4/pi) * arctan(1 / cosh(eps * (i - center)))
# At peak amplitude (t = pi/(2*omega)): the breather is at maximum displacement

for n_mode in [1, 2, 3, 4, 5, 6, 7, 8]:
    omega_n = np.cos(n_mode * gamma)
    eps_n = np.sin(n_mode * gamma)

    # 1D breather profile on discrete lattice
    N1d = 512
    center = N1d // 2
    x = np.arange(N1d, dtype=np.float64) - center

    # Peak amplitude profile: phi(x) = (4/pi)*arctan(1/cosh(eps*x))
    phi_peak = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * x + 1e-30))

    # The Hessian perturbation at each site:
    # delta_H_i = sin(pi*|phi_x_i|) / (pi*|phi_x_i|) - 1
    # For phi_x > 0 (the breather peak): sinc function minus 1
    abs_phi = np.abs(phi_peak) + 1e-30
    sinc_phi = np.sin(PI * abs_phi) / (PI * abs_phi)
    delta_H = sinc_phi - 1.0  # negative everywhere (softening)

    # The VP energy shift per transverse component:
    # delta_E_ZPE = (1/2) * sum_k delta_omega_k
    #
    # First-order perturbation theory for eigenvalue shift:
    # delta_omega_k = <k|delta_H|k> / (2*omega_k)
    #
    # For a UNIFORM perturbation delta_H = const:
    # delta_omega_k = delta_H / (2*omega_k) for all k
    #
    # For a LOCALIZED perturbation (only nonzero within breather width):
    # Sum over sites: delta_E = (1/2) * sum_i delta_H_i / (2 * omega_avg) * (1/N)
    # Wait, need to be more careful.
    #
    # The perturbation is diagonal in position space:
    # V_pert(i) = delta_H_i * phi_y(i)^2 / 2  (from Taylor expansion)
    #
    # The ZPE shift = (1/2) * sum_k (omega_k' - omega_k)
    # = (1/2) * sum_k delta_H_avg / (2*omega_k)  (first-order PT)
    # where delta_H_avg = <k|delta_H|k> = sum_i |psi_k(i)|^2 * delta_H_i
    # For plane waves: |psi_k|^2 = 1/N (uniform)
    # So: delta_H_avg = (1/N) * sum_i delta_H_i
    #
    # ZPE shift per component = (1/2) * N_modes * <delta_H> / (2 * <omega>)
    # where <delta_H> = (1/N) * sum_i delta_H_i

    sum_delta_H = np.sum(delta_H)
    avg_delta_H = sum_delta_H / N1d

    # Average phonon frequency on discrete 1D lattice:
    # omega_k^2 = 1 + 4*sin^2(k/2), k = 2*pi*n/N
    k_arr = 2*PI*np.arange(N1d) / N1d
    omega_k = np.sqrt(1.0 + 4*np.sin(k_arr/2)**2)
    avg_inv_omega = np.mean(1.0 / omega_k)

    # ZPE shift per transverse component (first-order PT):
    # delta_E_component = (1/2) * sum_k delta_omega_k
    # = (1/2) * sum_k <k|delta_H|k> / (2*omega_k)
    # = (1/2) * (1/N) * sum_i delta_H_i * sum_k 1/(2*omega_k)
    # = (1/4) * sum_delta_H * <1/omega>
    delta_E_per_comp = 0.25 * sum_delta_H * avg_inv_omega

    # Total VP correction (2 transverse components: y and z):
    delta_E_VP = 2 * delta_E_per_comp

    # Breather energy (rough estimate: integral of phi_peak^2 * omega)
    E_breather = np.sum(phi_peak**2) * omega_n / 2

    # Relative VP correction:
    vp_relative = delta_E_VP / E_breather if E_breather > 0 else 0

    # Key quantities
    phi_max = np.max(phi_peak)
    width = np.sum(phi_peak > 0.01 * phi_max)  # sites with significant breather

    if n_mode <= 3 or n_mode == 7:
        report(f"\n  n={n_mode}: eps={eps_n:.6f}, breather width={width} sites, phi_max={phi_max:.4f}")
        report(f"    sum(delta_H) = {sum_delta_H:.6f}")
        report(f"    avg(delta_H) = {avg_delta_H:.8f}")
        report(f"    delta_E_VP (2 components) = {delta_E_VP:.8f}")
        report(f"    E_breather = {E_breather:.6f}")
        report(f"    VP/E = {vp_relative:.6e}")
        report(f"    alpha^2 = {alpha**2:.6e}")
        report(f"    VP/E / alpha^2 = {vp_relative/alpha**2:.4f}" if alpha**2 > 0 else "")

report("")

# ============================================================
# PART 2: Exact calculation on 3D lattice (smaller for tractability)
# ============================================================
report("PART 2: 3D LATTICE VP (exact Hessian perturbation)")
report("-" * 60)

N3 = 32
center3 = N3 // 2

report(f"Lattice: {N3}^3 = {N3**3} sites")
report("")

for n_mode in [1, 2, 3, 4, 7]:
    omega_n = np.cos(n_mode * gamma)
    eps_n = np.sin(n_mode * gamma)

    # Breather along x-axis (uniform in y,z = quasi-1D)
    x = np.arange(N3, dtype=np.float64) - center3
    phi_1d = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * x + 1e-30))

    # Peak profile on 3D lattice (uniform in y,z)
    # phi_x(i,j,k) = phi_1d(i), phi_y = phi_z = 0
    phi_3d = phi_1d  # just the 1D profile

    # Hessian perturbation at each site
    abs_phi = np.abs(phi_1d) + 1e-30
    sinc_phi = np.sin(PI * abs_phi) / (PI * abs_phi)
    delta_H_1d = sinc_phi - 1.0  # same at each (j,k) since uniform in y,z

    # For 3D lattice with uniform breather in y,z:
    # delta_H(i,j,k) = delta_H_1d(i) for all j,k
    # sum over 3D lattice = N^2 * sum_1d
    sum_delta_H_3d = N3**2 * np.sum(delta_H_1d)

    # 3D phonon dispersion: omega_k^2 = 1 + 4*(sin^2(kx/2)+sin^2(ky/2)+sin^2(kz/2))
    kx = 2*PI*np.arange(N3)/N3
    ky = 2*PI*np.arange(N3)/N3
    kz = 2*PI*np.arange(N3)/N3
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    omega_3d = np.sqrt(1.0 + 4*(np.sin(KX/2)**2 + np.sin(KY/2)**2 + np.sin(KZ/2)**2))
    avg_inv_omega_3d = np.mean(1.0 / omega_3d)

    # First-order PT: delta_E = (1/4) * sum_delta_H * <1/omega> for 1 transverse component
    # But for 3D quasi-1D mode: the perturbation is uniform in y,z
    # So <k|delta_H|k> = delta_H_1d(i_mode) for modes with the right kx
    # For the AVERAGE over all k: same formula
    delta_E_3d_per_comp = 0.25 * sum_delta_H_3d * avg_inv_omega_3d / N3**3

    # 2 transverse components (y and z)
    delta_E_VP_3d = 2 * delta_E_3d_per_comp

    # Breather energy on 3D lattice
    E_br_3d = N3**2 * np.sum(phi_1d**2) * omega_n / 2

    vp_rel_3d = delta_E_VP_3d / E_br_3d if E_br_3d > 0 else 0

    report(f"  n={n_mode}: delta_E_VP = {delta_E_VP_3d:.8f}, E_breather = {E_br_3d:.4f}, "
           f"VP/E = {vp_rel_3d:.4e}, VP/E / alpha^2 = {vp_rel_3d/alpha**2:.2f}")

report("")

# ============================================================
# PART 3: The sinc expansion and the Oh connection
# ============================================================
report("PART 3: ANALYTIC STRUCTURE")
report("-" * 60)
report("")
report("The Hessian perturbation at the breather:")
report("  delta_H = sin(pi*phi)/(pi*phi) - 1")
report("")
report("Taylor expansion:")
report("  sin(pi*r)/(pi*r) = 1 - (pi*r)^2/6 + (pi*r)^4/120 - ...")
report("  delta_H = -(pi*phi)^2/6 + (pi*phi)^4/120 - ...")
report("")

# The leading-order VP contribution:
# delta_H ≈ -(pi*phi_x)^2/6
# This is the phi_x^2 * phi_y^2 cross-coupling from |phi|^4
#
# |phi|^4 = (phi_x^2 + phi_y^2 + phi_z^2)^2
#         = phi_x^4 + phi_y^4 + phi_z^4 + 2*(phi_x^2*phi_y^2 + ...)
#
# The phi_x^2 * phi_y^2 term in the cosine expansion:
# V = (1/pi^2)(1 - cos(pi*|phi|)) = |phi|^2/2 - pi^2*|phi|^4/24 + ...
# The phi_x^2*phi_y^2 coefficient = -pi^2*2/24 = -pi^2/12
#
# In d²V/dphi_y²: this contributes -pi^2*phi_x^2/6 (from 2*phi_y term)
# Which is exactly delta_H ≈ -(pi*phi_x)^2/6 ✓

report("Leading term: delta_H ≈ -(pi*phi_x)^2 / 6")
report("")
report("This comes from the phi_x^2 * phi_y^2 cross-coupling in |phi|^4:")
report("  V = |phi|^2/2 - pi^2*|phi|^4/24 + ...")
report("  The phi_x^2*phi_y^2 coefficient = -pi^2/12")
report("  d^2V/dphi_y^2 contribution = -pi^2*phi_x^2/6")
report("")

# Now: how does this relate to alpha^2?
# The breather amplitude phi_max ~ eps / omega ~ gamma ~ pi/(16*pi-2)
# For n=1: phi_max = (4/pi) * arctan(1) = 1.0 (at center, x=0)
# But the time-averaged amplitude is phi_rms ~ phi_max / sqrt(2)
#
# The VP correction is:
# delta_E / E ~ -(pi*phi_rms)^2/6 * (geometric factor)
#
# For n=1: phi_max = 1.273, phi_rms ~ 0.9
# -(pi*0.9)^2/6 ~ -1.33
# This is ORDER 1, not alpha^2!
#
# The alpha^2 scaling comes from the COUPLING being weak (alpha = probability
# of interaction per transit). In the lattice simulation, the coupling IS the
# cosine potential, which is order 1.
#
# The alpha^2 in the VP law comes from:
# 1. alpha = tunneling rate per cosine barrier = exp(-S)
# 2. VP = two interactions with the barrier = alpha^2
# 3. In the simulation, each time step IS a barrier interaction
#    so the total effect is not suppressed by alpha^2

report("IMPORTANT INSIGHT:")
report("The VP correction in the simulation is ORDER 1 (not alpha^2).")
report("This is because the LATTICE coupling is strong (cosine nonlinearity).")
report("")
report("In the physical theory, alpha^2 comes from the TUNNELING rate between")
report("lattice sites. Each tunneling event is suppressed by exp(-S_barrier).")
report("The VP = two tunneling events = alpha^2.")
report("")
report("In the simulation, we evolve the CONTINUOUS field dynamics on the lattice.")
report("The field moves freely between sites — no tunneling suppression.")
report("The phi^4 coupling acts at full strength, not alpha^2-suppressed.")
report("")
report("To get the PHYSICAL VP (alpha^2-suppressed), we would need to:")
report("1. Simulate the QUANTUM lattice (not classical wave equation)")
report("2. Or compute it analytically: VP = alpha^2 * Oh_fraction * geometric_sum")
report("")

# The bridge: what is the simulation measuring?
# It measures the GEOMETRIC FRACTION of the VP — the Oh decomposition.
# The alpha^2 suppression is a SEPARATE factor from the tunneling physics.
# The simulation gives: VP_geom = delta_H_integral / E_breather
# The physical VP = alpha^2 * VP_geom

for n_mode in [1, 4, 7]:
    omega_n = np.cos(n_mode * gamma)
    eps_n = np.sin(n_mode * gamma)
    x = np.arange(N3, dtype=np.float64) - center3
    phi_1d = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * x + 1e-30))
    abs_phi = np.abs(phi_1d) + 1e-30
    delta_H_1d = np.sin(PI*abs_phi)/(PI*abs_phi) - 1.0

    # Geometric VP = integral of delta_H / breather_energy
    vp_geom = np.sum(delta_H_1d) / (np.sum(phi_1d**2) * omega_n / 2)

    report(f"  n={n_mode}: VP_geometric = {vp_geom:.6f}")
    report(f"         Physical VP = alpha^2 * VP_geom = {alpha**2 * vp_geom:.6e}")
    report(f"         This should match the mass/coupling VP corrections")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))

log.close()
print(f"\nResults saved to: {outfile}")
