"""
VP Geometric Fraction — Clean Derivation
==========================================
The VP correction softens the transverse modes by:
  delta_H(i) = sin(pi*phi_x(i)) / (pi*phi_x(i)) - 1 ≈ -(pi*phi_x)^2/6

The PHYSICAL VP = how much this softening affects the breather ITSELF.
Not the total lattice ZPE (extensive), but the breather self-energy (intensive).

The breather self-energy correction = overlap of the perturbation with
the breather's own wavefunction:

  VP_self = <breather|delta_H|breather> = sum_i phi(i)^2 * delta_H(i) / sum_i phi(i)^2

This is a PURE NUMBER — independent of lattice size, alpha, everything.
It's the geometric fraction from the Oh decomposition.
"""
import sys, io, os
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
gamma = PI / (2**(d+1)*PI - 2)
alpha = np.exp(-(2/6) * (2**7/PI**2 + np.log(6)))

outfile = os.path.join(os.path.dirname(__file__), "breather_vp_pindown2_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("VP GEOMETRIC FRACTION — CLEAN DERIVATION")
report("=" * 70)
report(f"alpha = {alpha:.6f}, alpha^2 = {alpha**2:.6e}")
report(f"GWT predictions:")
report(f"  1/2^(d/2) = {1/2**(d/2):.6f}  (confined, mass ratio)")
report(f"  (d^2-1)/d^2 = 8/9 = {8/9:.6f}  (free, alpha)")
report(f"  (d^2-1)/d = 8/3 = {8/3:.6f}  (gluon, alpha_s)")
report("")

# ============================================================
# THE BREATHER SELF-ENERGY: WEIGHTED AVERAGE
# ============================================================
report("THE BREATHER SELF-ENERGY OVERLAP")
report("-" * 60)
report("")
report("VP_self = <phi^2 * delta_H> / <phi^2>")
report("       = sum[phi(i)^2 * delta_H(i)] / sum[phi(i)^2]")
report("")
report("This is the fraction of the on-site potential that's modified by")
report("the cross-component coupling, weighted by where the breather lives.")
report("")

N = 4096  # large lattice for convergence
center = N // 2

report(f"{'n':>3} {'VP_self':>12} {'VP_lead':>12} {'phi_max':>8} {'width':>6} "
       f"{'sum_phi2':>10} {'sum_phi2dH':>12}")
report("-" * 70)

vp_selfs = []

for n_mode in range(1, 25):
    omega_n = np.cos(n_mode * gamma)
    eps_n = np.sin(n_mode * gamma)

    if omega_n < 0.005:
        continue

    x = np.arange(N, dtype=np.float64) - center
    phi = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * x + 1e-30))

    abs_phi = np.abs(phi) + 1e-30
    delta_H = np.sin(PI * abs_phi) / (PI * abs_phi) - 1.0

    sum_phi2 = np.sum(phi**2)
    sum_phi2_dH = np.sum(phi**2 * delta_H)

    VP_self = sum_phi2_dH / sum_phi2  # weighted average

    # Leading-order: delta_H ≈ -(pi*phi)^2/6
    # So VP_self_lead ≈ -(pi^2/6) * <phi^4> / <phi^2>
    sum_phi4 = np.sum(phi**4)
    VP_lead = -(PI**2/6) * sum_phi4 / sum_phi2

    phi_max = np.max(phi)
    width = np.sum(phi > 0.01 * phi_max)

    vp_selfs.append((n_mode, VP_self, VP_lead, omega_n, eps_n))

    report(f"{n_mode:>3} {VP_self:12.6f} {VP_lead:12.6f} {phi_max:8.4f} {width:6d} "
           f"{sum_phi2:10.4f} {sum_phi2_dH:12.6f}")

report("")

# ============================================================
# ANALYSIS: What determines VP_self?
# ============================================================
report("ANALYSIS: WHAT DETERMINES VP_SELF?")
report("-" * 60)
report("")

# VP_self ≈ -(pi^2/6) * <phi^4>/<phi^2>
# For the sine-Gordon breather profile phi = (4/pi)*arctan(1/cosh(eps*x)):
# <phi^4>/<phi^2> depends on the SHAPE of the breather, not its size.
# Let's compute the kurtosis: <phi^4> / <phi^2>^2 * N_eff
# and also <phi^4>/<phi^2> directly.

report(f"{'n':>3} {'<phi^4>/<phi^2>':>16} {'phi_max^2':>10} {'ratio':>10} "
       f"{'kurtosis':>10} {'VP_self':>10}")
report("-" * 62)

for n_mode, VP_self, VP_lead, omega_n, eps_n in vp_selfs:
    x = np.arange(N, dtype=np.float64) - center
    phi = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * x + 1e-30))

    phi2 = np.sum(phi**2)
    phi4 = np.sum(phi**4)
    ratio_42 = phi4 / phi2  # <phi^4>/<phi^2>
    phi_max = np.max(phi)
    kurtosis = phi4 * np.sum(np.ones_like(phi)) / phi2**2  # excess kurtosis proxy

    # The ratio <phi^4>/<phi^2> for a sech^2 profile:
    # phi = A * sech(eps*x), then:
    # <phi^4>/<phi^2> = A^2 * integral(sech^4) / integral(sech^2) = A^2 * (2/3)
    # For our profile: A = 4/pi * arctan(1) at x=0... but it's not exactly sech.

    report(f"{n_mode:>3} {ratio_42:16.6f} {phi_max**2:10.4f} "
           f"{ratio_42/phi_max**2:10.4f} {kurtosis:10.2f} {VP_self:10.6f}")

report("")

# ============================================================
# THE KEY RATIO: VP_self normalized by the breather peak
# ============================================================
report("KEY RATIO: VP_self / phi_max^2")
report("-" * 60)
report("")
report("If VP_self = -(pi^2/6) * C * phi_max^2, then C = VP_self / (-(pi^2/6) * phi_max^2)")
report("C encodes the SHAPE of the breather (kurtosis-like).")
report("")

report(f"{'n':>3} {'C_shape':>12} {'VP_self':>12} {'-(pi^2/6)*phimax^2':>20}")
report("-" * 50)

C_values = []
for n_mode, VP_self, VP_lead, omega_n, eps_n in vp_selfs:
    x = np.arange(N, dtype=np.float64) - center
    phi = (4.0/PI) * np.arctan(1.0 / np.cosh(eps_n * x + 1e-30))
    phi_max = np.max(phi)

    leading = -(PI**2/6) * phi_max**2
    C = VP_self / leading if abs(leading) > 1e-30 else 0

    C_values.append((n_mode, C))
    report(f"{n_mode:>3} {C:12.6f} {VP_self:12.6f} {leading:20.6f}")

report("")

# ============================================================
# THE PHYSICAL VP: putting it all together
# ============================================================
report("PHYSICAL VP CORRECTION")
report("-" * 60)
report("")
report("The VP correction to the breather frequency/mass:")
report("  delta_m/m = alpha^2 * (d-1) * VP_self")
report("")
report("(d-1) = 2 transverse components (y and z)")
report("alpha^2 comes from the tunneling physics (not in the Hessian)")
report("VP_self is the geometric overlap (from the Hessian)")
report("")

report(f"{'n':>3} {'VP_self':>12} {'(d-1)*VP_self':>14} {'alpha^2*..':>14} {'GWT pred':>12}")
report("-" * 60)

for n_mode, VP_self, VP_lead, omega_n, eps_n in vp_selfs[:10]:
    correction = (d-1) * VP_self
    physical = alpha**2 * correction
    # GWT prediction for the mass ratio: alpha^2 / 2^(d/2)
    gwt_pred = alpha**2 / 2**(d/2)

    report(f"{n_mode:>3} {VP_self:12.6f} {correction:14.6f} {physical:14.6e} "
           f"{gwt_pred:12.6e}")

report("")

# What if we include the (d-1) factor in the geometric fraction?
report("COMPARISON:")
for n_mode, VP_self, VP_lead, omega_n, eps_n in vp_selfs[:1]:
    report(f"  n={n_mode}: (d-1) * VP_self = {(d-1)*VP_self:.6f}")
    report(f"  1/2^(d/2)           = {1/2**(d/2):.6f}")
    report(f"  Ratio               = {(d-1)*VP_self / (1/2**(d/2)):.4f}")
    report(f"  So VP_self          = {VP_self:.6f}")
    report(f"  1/[2*(d-1)*2^(d/2)] = {1/(2*(d-1)*2**(d/2)):.6f}")
    report("")

# What's VP_self for n=1?
n1_vp = vp_selfs[0][1]
report(f"VP_self(n=1) = {n1_vp:.8f}")
report("")

# Is it a clean fraction of pi?
report("Searching for clean expression:")
candidates = {
    "-1/6": -1/6,
    "-pi^2/6^2": -PI**2/36,
    "-2/3*pi^2/6": -2*PI**2/18,
    "-(4/pi)^2/6": -(4/PI)**2/6,
    "-(4/pi)^2 * 2/(3*pi^2)": -(4/PI)**2 * 2/(3*PI**2),
    "-8/(3*pi^2)": -8/(3*PI**2),
    "-2/(3*pi)": -2/(3*PI),
    "-1/(2*sqrt(2))": -1/(2*np.sqrt(2)),
    "-1/pi^2": -1/PI**2,
    "-2/pi^2": -2/PI**2,
    "-4/(3*pi^2)": -4/(3*PI**2),
    "-8/(3*pi^4)": -8/(3*PI**4),
    "-(4/pi)^2 * int(sech^4)/int(sech^2)": -(4/PI)**2 * (2/3),
    "-16/(3*pi^2)": -16/(3*PI**2),
    "-32/(9*pi^2)": -32/(9*PI**2),
    "-(16/pi^2)*(2/3)/(pi^2/6)": -(16/PI**2)*(2/3)/(PI**2/6),
}

for name, val in sorted(candidates.items(), key=lambda x: abs(x[1] - n1_vp)):
    ratio = n1_vp / val if abs(val) > 1e-30 else 999
    if 0.9 < abs(ratio) < 1.1:
        report(f"  {name:>35} = {val:12.8f}, ratio = {ratio:.6f} {'<-- MATCH!' if abs(ratio-1) < 0.02 else ''}")
    elif abs(ratio - 1) < 0.2:
        report(f"  {name:>35} = {val:12.8f}, ratio = {ratio:.6f}")

report("")

# Also compute the integral analytically for comparison
# For phi = (4/pi)*arctan(1/cosh(eps*x)):
# At continuous limit (eps -> 0, width -> inf):
# <phi^4>/<phi^2> -> (4/pi)^2 * integral(arctan(1/cosh(t))^4 dt) / integral(arctan(1/cosh(t))^2 dt)
#
# Let's compute this numerically with very fine grid
t = np.linspace(-100, 100, 1000000)
f = np.arctan(1.0 / np.cosh(t))
int_f2 = np.trapezoid(f**2, t)
int_f4 = np.trapezoid(f**4, t)
ratio_cont = int_f4 / int_f2

report(f"Continuum limit integrals:")
report(f"  int(arctan(1/cosh(t))^2 dt) = {int_f2:.8f}")
report(f"  int(arctan(1/cosh(t))^4 dt) = {int_f4:.8f}")
report(f"  Ratio int(f^4)/int(f^2)     = {ratio_cont:.8f}")
report(f"  (4/pi)^2 * ratio            = {(4/PI)**2 * ratio_cont:.8f}")
report(f"  VP_self_continuum = -(pi^2/6) * (4/pi)^2 * ratio = {-(PI**2/6)*(4/PI)**2*ratio_cont:.8f}")
report("")

# With d-1 transverse components:
vp_cont = -(PI**2/6) * (4/PI)**2 * ratio_cont
report(f"  (d-1) * VP_self_continuum = {(d-1)*vp_cont:.8f}")
report(f"  1/2^(d/2) = {1/2**(d/2):.8f}")
report(f"  Ratio: {(d-1)*vp_cont / (1/2**(d/2)):.6f}")

report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))

log.close()
print(f"\nResults saved to: {outfile}")
