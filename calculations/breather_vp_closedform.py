"""
Find the Exact Form of VP_self
================================
VP_self = integral[f^2 * (sinc(pi*f) - 1) dt] / integral[f^2 dt]

where f(t) = (4/pi) * arctan(1/cosh(t))  [sine-Gordon breather peak profile]

This is a pure number determined by the Lagrangian. Compute it to high
precision and search for closed-form expressions.

IMPORTANT: Only accept candidates with a clear derivation.
"""
import sys, io, os
import numpy as np
from scipy import integrate
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3

outfile = os.path.join(os.path.dirname(__file__), "breather_vp_closedform_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("EXACT FORM OF VP_SELF")
report("=" * 70)
report("")

# ============================================================
# STEP 1: High-precision computation of VP_self
# ============================================================
report("STEP 1: HIGH-PRECISION COMPUTATION")
report("-" * 60)

def f_breather(t):
    """Sine-Gordon breather peak profile: (4/pi)*arctan(sech(t))"""
    return (4.0/PI) * np.arctan(1.0 / np.cosh(t))

def integrand_num(t):
    """f^2 * (sinc(pi*f) - 1)"""
    f = f_breather(t)
    f = np.where(np.abs(f) < 1e-15, 1e-15, f)
    sinc_val = np.sin(PI * f) / (PI * f)
    return f**2 * (sinc_val - 1)

def integrand_den(t):
    """f^2"""
    return f_breather(t)**2

# Compute integrals using scipy.integrate.quad (high precision)
I_num, err_num = integrate.quad(integrand_num, 0, np.inf, limit=200, epsabs=1e-15, epsrel=1e-15)
I_den, err_den = integrate.quad(integrand_den, 0, np.inf, limit=200, epsabs=1e-15, epsrel=1e-15)

# Factor of 2 for symmetric integral (-inf to inf) / same for both => cancels
VP_self = I_num / I_den

report(f"int[f^2 * (sinc(pi*f)-1) dt, 0..inf] = {I_num:.15f} +/- {err_num:.2e}")
report(f"int[f^2 dt, 0..inf]                   = {I_den:.15f} +/- {err_den:.2e}")
report(f"VP_self = {VP_self:.15f}")
report("")

# Also compute the sub-integrals for the Taylor expansion
def integrand_f4(t):
    return f_breather(t)**4

def integrand_f6(t):
    return f_breather(t)**6

I_f2 = I_den  # same as denominator
I_f4, _ = integrate.quad(integrand_f4, 0, np.inf, limit=200, epsabs=1e-15)
I_f6, _ = integrate.quad(lambda t: f_breather(t)**6, 0, np.inf, limit=200, epsabs=1e-15)
I_f8, _ = integrate.quad(lambda t: f_breather(t)**8, 0, np.inf, limit=200, epsabs=1e-15)

report("Sub-integrals (all from 0 to inf):")
report(f"  I_f2 = int[f^2] = {I_f2:.15f}")
report(f"  I_f4 = int[f^4] = {I_f4:.15f}")
report(f"  I_f6 = int[f^6] = {I_f6:.15f}")
report(f"  I_f8 = int[f^8] = {I_f8:.15f}")
report("")

# Ratios (these are the key shape parameters)
R42 = I_f4 / I_f2
R62 = I_f6 / I_f2
R82 = I_f8 / I_f2
report("Shape ratios:")
report(f"  <f^4>/<f^2> = {R42:.15f}")
report(f"  <f^6>/<f^2> = {R62:.15f}")
report(f"  <f^8>/<f^2> = {R82:.15f}")
report("")

# Verify Taylor expansion:
# sinc(pi*f) - 1 = -(pi*f)^2/6 + (pi*f)^4/120 - (pi*f)^6/5040 + ...
VP_taylor2 = -(PI**2/6) * R42
VP_taylor4 = VP_taylor2 + (PI**4/120) * R62
VP_taylor6 = VP_taylor4 - (PI**6/5040) * R82

report("Taylor expansion check:")
report(f"  Leading -(pi^2/6)*<f^4>/<f^2>:           {VP_taylor2:.15f}")
report(f"  + (pi^4/120)*<f^6>/<f^2>:                {VP_taylor4:.15f}")
report(f"  + -(pi^6/5040)*<f^8>/<f^2>:              {VP_taylor6:.15f}")
report(f"  Exact VP_self:                            {VP_self:.15f}")
report(f"  6th-order residual:                       {VP_self - VP_taylor6:.6e}")
report("")

# ============================================================
# STEP 2: Identify the sub-integrals
# ============================================================
report("STEP 2: IDENTIFY SUB-INTEGRALS")
report("-" * 60)
report("")

# f(t) = (4/pi)*arctan(sech(t))
# Let u = sech(t), du = -sech(t)*tanh(t)*dt = -u*sqrt(1-u^2)*dt
# At t=0: u=1. At t=inf: u=0.
# dt = -du / (u*sqrt(1-u^2))
# f = (4/pi)*arctan(u)
#
# int_0^inf f^2 dt = int_0^1 [(4/pi)*arctan(u)]^2 * du/(u*sqrt(1-u^2))
#                  = (16/pi^2) * int_0^1 arctan(u)^2 / (u*sqrt(1-u^2)) du

# Let's verify:
def integrand_u_f2(u):
    if u < 1e-15 or u > 1-1e-15:
        return 0
    return (16/PI**2) * np.arctan(u)**2 / (u * np.sqrt(1-u**2))

I_f2_check, _ = integrate.quad(integrand_u_f2, 0, 1, limit=200, epsabs=1e-15)
report(f"I_f2 via u-substitution: {I_f2_check:.15f} (should match {I_f2:.15f})")

# Try known integrals involving arctan(u)/sqrt(1-u^2):
# int_0^1 arctan(u)^2 / (u*sqrt(1-u^2)) du = ?
#
# Related Catalan-type integrals:
# int_0^1 arctan(u)/sqrt(1-u^2) du = pi/2 * ln(1+sqrt(2))  [known]
# int_0^1 arctan(u)^2/sqrt(1-u^2) du = ?

# Compute these:
I_at_sqrt, _ = integrate.quad(lambda u: np.arctan(u)/np.sqrt(1-u**2+1e-30), 0, 1-1e-10, limit=200)
known_val = PI/2 * np.log(1 + np.sqrt(2))
report(f"\nKnown integral: int_0^1 arctan(u)/sqrt(1-u^2) du")
report(f"  Computed: {I_at_sqrt:.15f}")
report(f"  pi/2*ln(1+sqrt(2)): {known_val:.15f}")
report(f"  Match: {abs(I_at_sqrt - known_val) < 1e-10}")
report("")

# ============================================================
# STEP 3: Search for closed form of VP_self
# ============================================================
report("STEP 3: SYSTEMATIC SEARCH FOR CLOSED FORM")
report("-" * 60)
report(f"Target: VP_self = {VP_self:.15f}")
report("")

# Build candidates from the constants that appear in the theory:
# pi, sqrt(2), sqrt(3), ln(2), Catalan G, gamma_SG, d=3

G_catalan = 0.915965594177219  # Catalan's constant
ln2 = np.log(2)
sqrt2 = np.sqrt(2)
sqrt3 = np.sqrt(3)
gamma_sg = PI / (2**(d+1)*PI - 2)

# Systematic search: a*X + b*Y where a,b are small rationals and X,Y are constants
constants = {
    '1': 1.0,
    'pi': PI,
    'pi^2': PI**2,
    '1/pi': 1/PI,
    '1/pi^2': 1/PI**2,
    'sqrt(2)': sqrt2,
    'sqrt(3)': sqrt3,
    'ln(2)': ln2,
    'G_cat': G_catalan,
    'pi*ln2': PI*ln2,
    'pi^2/6': PI**2/6,
    'ln(1+sqrt2)': np.log(1+sqrt2),
}

# Try VP_self = -p/q * const for small p,q
report("Single-constant search: VP_self = -(p/q) * C")
report(f"{'expression':>40} {'value':>18} {'error':>12}")
report("-" * 72)

matches = []
target = -VP_self  # work with positive target

for name, val in constants.items():
    for p in range(1, 20):
        for q in range(1, 20):
            candidate = p * val / q
            if abs(candidate) < 0.01:
                continue
            err = abs(candidate - target) / target
            if err < 0.001:  # within 0.1%
                expr = f"({p}/{q}) * {name}"
                matches.append((err, expr, candidate))

# Also try two-constant combinations
for n1, v1 in constants.items():
    for n2, v2 in constants.items():
        if n1 >= n2:
            continue
        for p1 in range(-5, 6):
            for q1 in range(1, 8):
                for p2 in range(-5, 6):
                    for q2 in range(1, 8):
                        if p1 == 0 and p2 == 0:
                            continue
                        candidate = abs(p1*v1/q1 + p2*v2/q2)
                        if candidate < 0.01:
                            continue
                        err = abs(candidate - target) / target
                        if err < 0.0001:  # within 0.01%
                            expr = f"({p1}/{q1})*{n1} + ({p2}/{q2})*{n2}"
                            matches.append((err, expr, candidate))

matches.sort()
for err, expr, val in matches[:20]:
    report(f"  {expr:>40} = {-val:18.12f}  err={err:.2e}")

report("")

# ============================================================
# STEP 4: Check specific physics-motivated candidates
# ============================================================
report("STEP 4: PHYSICS-MOTIVATED CANDIDATES")
report("-" * 60)
report("")

candidates = {
    # From the sinc expansion and breather profile
    "-(pi^2/6) * <f^4>/<f^2> (leading order)": VP_taylor2,
    # Clean fractions
    "-3/4": -3/4,
    "-pi/4": -PI/4,
    # From Oh group theory
    "-(d^2-1)/d^2 = -8/9": -8/9,
    "-1/2^(d/2) = -1/(2*sqrt(2))": -1/(2*sqrt2),
    # From breather physics
    "-(4/pi)^2 * 2/3 / pi": -(4/PI)**2 * 2/(3*PI),
    # Involving Catalan's constant
    "-G_catalan * pi / (d+1)": -G_catalan * PI / 4,
    "-2*G_catalan/pi": -2*G_catalan/PI,
    # From the known integral pi/2*ln(1+sqrt2)
    "-pi/2 * ln(1+sqrt2) / pi": -0.5 * np.log(1+sqrt2),
    "-ln(1+sqrt2)^2 / 2": -np.log(1+sqrt2)**2 / 2,
    "-2*ln(1+sqrt2)/pi": -2*np.log(1+sqrt2)/PI,
    # Composite
    "-pi^2/12 - 1/4": -PI**2/12 - 1/4,
    "-(pi^2-6)/(2*pi)": -(PI**2-6)/(2*PI),
    "-pi^2/4/pi = -pi/4": -PI/4,
    # Kink mass related
    "-8/pi^2 * something": -8/PI**2 * PI * np.sqrt(2)/3,
    # Try: VP_self = -(2/3) * (Catalan/1)
    "-2*G_cat/3": -2*G_catalan/3,
    # sin and cos of specific angles
    "-sin(pi/4) = -1/sqrt(2)": -1/sqrt2,
    "-cos(pi/d) = -cos(60) = -1/2": -0.5,
    "-cos(pi/4) = -1/sqrt(2)": -1/sqrt2,
    # From lattice: 2/pi^2 is the barrier height
    "-(2/pi^2)*pi^2*..": -(2/PI**2) * PI**3 / 8,
}

report(f"{'candidate':>45} {'value':>16} {'error':>12}")
report("-" * 75)

results = []
for name, val in candidates.items():
    err = abs(val - VP_self) / abs(VP_self)
    results.append((err, name, val))

results.sort()
for err, name, val in results:
    marker = " <== MATCH" if err < 0.005 else ""
    report(f"  {name:>45} {val:16.12f} {err:12.6e}{marker}")

report("")

# ============================================================
# STEP 5: Direct analytical approach
# ============================================================
report("STEP 5: DIRECT ANALYTICAL APPROACH")
report("-" * 60)
report("")
report("VP_self = integral[f^2 * (sinc(pi*f)-1)] / integral[f^2]")
report("       = integral[f^2 * (sin(4*h)/(4*h) - 1)] / integral[f^2]")
report("where h = arctan(sech(t)), f = (4/pi)*h")
report("")

# Let's try computing VP_self as a series in terms of known integrals
# sinc(pi*f) = sinc(4*h) = sum_{n=1}^inf (-1)^n * (4h)^{2n} / (2n+1)!
# f^2 * (sinc(pi*f) - 1) = (4/pi)^2 * h^2 * sum_{n=1}^inf (-1)^n * (4h)^{2n} / (2n+1)!
# = (4/pi)^2 * sum_{n=1}^inf (-1)^n * 4^{2n} / (2n+1)! * h^{2n+2}
# = (16/pi^2) * sum_{n=1}^inf (-1)^n * 16^n / (2n+1)! * h^{2n+2}

# So VP_self = sum_{n=1}^inf (-1)^n * 16^n / (2n+1)! * J_{2n+2} / J_2
# where J_k = int_0^inf [arctan(sech(t))]^k dt

# Compute J_k integrals:
J = {}
for k in range(2, 18, 2):
    val, _ = integrate.quad(lambda t: np.arctan(1.0/np.cosh(t))**k, 0, np.inf,
                            limit=200, epsabs=1e-15)
    J[k] = val

report("J_k = int_0^inf [arctan(sech(t))]^k dt:")
for k, val in J.items():
    # Try to identify J_k
    report(f"  J_{k} = {val:.15f}")

report("")

# VP_self from the series:
VP_series = 0
report("Series: VP_self = sum (-1)^n * 16^n / (2n+1)! * J_{2n+2} / J_2")
for n in range(1, 8):
    term = (-1)**n * 16**n / np.math.factorial(2*n+1) * J.get(2*n+2, 0) / J[2]
    VP_series += term
    report(f"  n={n}: term = {term:+.10f}, cumulative = {VP_series:.10f}")

report(f"\n  Exact VP_self = {VP_self:.10f}")
report(f"  Series (7 terms) = {VP_series:.10f}")
report(f"  Residual = {VP_self - VP_series:.2e}")
report("")

# Try to identify J_2
report("Identifying J_2:")
J2 = J[2]
report(f"  J_2 = {J2:.15f}")
report(f"  pi^2/16 = {PI**2/16:.15f}")
report(f"  Ratio J_2/(pi^2/16) = {J2/(PI**2/16):.10f}")
report(f"  pi*ln(1+sqrt2)/4 = {PI*np.log(1+sqrt2)/4:.15f}")

# Actually, int_0^inf arctan(sech(t))^2 dt should be a known integral
# Let me check: this is related to the integral of the Gudermannian squared
# gd(t) = 2*arctan(tanh(t/2)), so arctan(sech(t)) is related but not identical

# Check: is J_2 = pi*G_catalan/4 ?
report(f"  pi*G_catalan/4 = {PI*G_catalan/4:.15f}")
report(f"  G_catalan/2 = {G_catalan/2:.15f}")
# Check: J_2 = pi/2 * ln(1+sqrt2) - pi^2/16 ?
report(f"  pi/2*ln(1+sqrt2) - pi^2/16 = {PI/2*np.log(1+sqrt2) - PI**2/16:.15f}")

report("")
report("Completed.")
log.close()
print(f"\nResults saved to: {outfile}")
