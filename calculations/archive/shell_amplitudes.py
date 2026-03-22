"""
Shell Breather Amplitudes from Kink Hessian
=============================================
For each atom (Z), build the kink well on the discrete lattice and solve
for the bound state wavefunctions. The AMPLITUDE of each wavefunction at
each shell's characteristic radius gives the sinc factor for screening.

The kink well: cos(pi * phi_kink(x)) where phi_kink = (4/pi)*arctan(1/cosh(beta*x))
beta = sqrt(Z) controls the well depth.

The bound states are the Hessian eigenvectors with eigenvalue < 1 (below mass gap).
The eigenVECTOR |psi_n(x)|^2 tells us where each bound electron "lives."

The breather amplitude at position x IS the kink profile phi_kink(x).
So sinc_factor(shell_n) = |sinc(pi * phi_kink(r_shell_n))| where r_shell_n
is the radius of the n-th shell.

Key insight: we don't need to solve for the eigenvectors at all!
The sinc factor depends on the KINK PROFILE at each shell's radius,
not on the electron wavefunctions. The kink profile is known analytically:
  phi_kink(r) = (4/pi) * arctan(1/cosh(sqrt(Z) * r))

So sinc_factor(r) = |sinc(pi * (4/pi) * arctan(1/cosh(sqrt(Z) * r)))|
                   = |sin(4*arctan(sech(sqrt(Z)*r))) / (4*arctan(sech(sqrt(Z)*r)))|

At r=0 (nucleus): phi=1, sinc(pi)=0  (maximum softening)
At r=inf:         phi=0, sinc(0)=1   (no softening)
At r = r_n:       interpolation       (shell-dependent)
"""
import sys, io, os
import numpy as np
from math import factorial
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/PI**2 + np.log(2*d)))

outfile = os.path.join(os.path.dirname(__file__), "shell_amplitudes_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("SHELL BREATHER AMPLITUDES FROM KINK PROFILE")
report("=" * 70)
report("")

# ============================================================
# THE KINK PROFILE AND SINC FACTOR
# ============================================================
def kink_profile(r, Z):
    """Kink field at distance r from nucleus with charge Z.
    phi(r) = (4/pi) * arctan(1/cosh(sqrt(Z) * r))
    """
    beta = np.sqrt(Z)
    return (4.0/PI) * np.arctan(1.0 / np.cosh(beta * r + 1e-30))

def sinc_factor(r, Z):
    """The VP sinc softening at distance r from a Z-nucleus.
    Returns |sinc(pi * phi_kink(r))| = |sin(4*arctan(sech(beta*r))) / (4*arctan(sech(beta*r)))|

    At r=0: phi=1, sinc(pi)=0 (full softening)
    At r=inf: phi=0, sinc(0)=1 (no softening)
    """
    phi = kink_profile(r, Z)
    if np.isscalar(phi):
        if abs(phi) < 1e-15:
            return 1.0
        return abs(np.sin(PI * phi) / (PI * phi))
    else:
        result = np.ones_like(phi)
        mask = np.abs(phi) > 1e-15
        result[mask] = np.abs(np.sin(PI * phi[mask]) / (PI * phi[mask]))
        return result

# ============================================================
# SHELL RADII
# ============================================================
# The characteristic radius of shell n in a Z-atom:
# In atomic units: r_n ~ n^2 / Z (Bohr model)
# On the lattice (Planck units): need conversion
#
# But the kink profile has its own scale: beta = sqrt(Z)
# The kink width is ~ 1/beta = 1/sqrt(Z)
# The n-th shell in kink units: r_n ~ n / sqrt(Z)
# (because the Bohr radius scales as 1/Z, and the kink as 1/sqrt(Z),
#  giving r_n_kink = n * a_0 * sqrt(Z) ~ n/sqrt(Z) * Z ... )
#
# Actually, let me think about this more carefully.
# The kink profile on the discrete lattice has scale 1/sqrt(Z).
# The hydrogen Bohr radius a_0 = 1/(alpha * m_e) in natural units.
# On the lattice, a_0 maps to ~ 1/alpha in lattice units (very large).
# The n-th shell is at r_n ~ n^2 * a_0 / Z ~ n^2/(alpha*Z) in lattice units.
#
# But the kink width is 1/sqrt(Z). The ratio:
# r_n / kink_width = n^2 / (alpha * Z) * sqrt(Z) = n^2 / (alpha * sqrt(Z))
#
# For Z=1, n=1: r/w = 1/alpha ~ 137. The electron is 137 kink widths away!
# The sinc factor at r=137/sqrt(Z) is essentially sinc(0) = 1 (no softening).
#
# This means: for the VALENCE electron, the sinc correction is negligible.
# The softening only matters for INNER shells near the nucleus.
#
# But wait — the SCREENING is from INNER shells affecting the VALENCE electron.
# The question is: at what radius does the INNER shell sit relative to the kink?
#
# Inner shell at n=1, Z=80: r ~ 1/(alpha*80) in kink widths... still very far.
#
# Hmm. Maybe the conversion is different. The kink IS the nucleus. The electron
# breathers orbit at their Bohr radii. But the kink profile decays as
# 1/cosh(sqrt(Z)*r). At the Bohr radius, the kink field is essentially zero.
#
# So the sinc factor is ~1 for ALL shells? That would mean no correction at all!

# Let me check numerically:
report("KINK PROFILE AT SHELL RADII")
report("-" * 55)
report("")

for Z in [1, 10, 30, 80]:
    report(f"Z = {Z}:")
    beta = np.sqrt(Z)
    report(f"  Kink width = 1/sqrt(Z) = {1/beta:.4f} lattice units")

    for n_shell in [1, 2, 3, 4, 5, 6]:
        # Shell radius in lattice units
        # The Bohr model: r_n = n^2 * a_0 / Z_eff
        # In Planck units: a_0 = 1/(alpha * m_e_planck)
        # But on our lattice, we set a=1 (lattice spacing = 1 Planck length)
        # The Bohr radius in lattice units: a_0 = 1/alpha ~ 137 sites
        # Shell n with effective Z_eff: r_n = n^2 * a_0 / Z_eff = n^2 * 137 / Z_eff

        # But Z_eff for inner shells is approximately Z (unscreened)
        # r_1 = 137/Z for the 1s shell
        r_bohr = n_shell**2 / (alpha_em * Z)  # in lattice units
        phi_at_r = kink_profile(r_bohr, Z)
        sinc_at_r = sinc_factor(r_bohr, Z)

        report(f"  n={n_shell}: r_Bohr={r_bohr:.1f}, phi={phi_at_r:.6f}, sinc={sinc_at_r:.6f}")
    report("")

# The radii are HUGE (r_1 ~ 137/Z ~ 2-137 lattice units)
# and the kink decays in ~1/sqrt(Z) ~ 0.1-1 units.
# So phi is essentially 0 and sinc is essentially 1 at all shell radii.

report("CONCLUSION: The shell radii (in Bohr) are MUCH larger than the")
report("kink width (in Planck). The kink profile is essentially zero at")
report("all electron shell positions. The sinc factor is ~1 everywhere.")
report("")
report("This means the VP sinc correction to screening is NEGLIGIBLE")
report("when using physical shell radii. The correction would only matter")
report("if shells were within ~1/sqrt(Z) of the nucleus, which they aren't.")
report("")

# ============================================================
# ALTERNATIVE: RELATIVE AMPLITUDE
# ============================================================
report("ALTERNATIVE INTERPRETATION")
report("-" * 55)
report("")
report("Maybe the sinc factor doesn't use the SPATIAL position of the shell")
report("but the FIELD AMPLITUDE of the breather mode itself.")
report("")
report("Each electron in a filled shell contributes to the total breather")
report("field phi at the nucleus. The more electrons, the larger the total")
report("phi, and the stronger the sinc softening.")
report("")
report("For a shell with occupation f (0 to 1):")
report("  phi_total = f × phi_max = f × 1.0")
report("  sinc_factor = |sinc(pi × f)|")
report("")
report("But this is what we tried before (too aggressive).")
report("The issue: phi_max = 1 assumes ALL electrons contribute at the")
report("nucleus. In reality, higher shells (larger n) have smaller")
report("amplitude at the nucleus (the wavefunction is spread out).")
report("")

# In hydrogen, |psi_ns(0)|^2 = (Z/n)^3 / (pi * a_0^3)
# The amplitude at the nucleus scales as (Z/n)^{3/2}
# For the sinc argument: phi_effective = sum_shells (Z/n)^{3/2} × count / max
# This is the NUCLEAR CONTACT DENSITY

report("NUCLEAR CONTACT DENSITY approach:")
report("  phi_eff(shell n, l) = (Z_eff/n)^(3/2) for s-electrons (l=0)")
report("  phi_eff = 0 for l > 0 (p, d, f don't reach nucleus)")
report("")
report("Only s-electrons contribute to phi at the nucleus!")
report("The sinc softening is from the s-electron density.")
report("")

# For an atom with s-electrons at shells n=1,2,...:
# phi_total ~ sum_n count_s(n) × (Z_eff_n/n)^{3/2}
# Normalized to phi_max: phi_total / max_possible

# For the SCREENING of one shell by another:
# The sinc factor should use the amplitude of the SCREENING shell
# at the location of the SCREENED shell (not at the nucleus).

# Shell n screens shell n': the screening field has amplitude
# proportional to the overlap integral of the two wavefunctions.
# For s-shells: the overlap ~ (Z_eff)^3 / n^3 × (Z_eff')^3 / n'^3 × ...
# This is getting complicated.

# Let me try the simplest physically motivated version:
# The sinc factor for a d-shell screening an s-shell depends on the
# d-shell's amplitude AT the s-shell's radius.
# Since d-electrons have l=2, their amplitude at the nucleus is 0.
# But their amplitude at r ~ n^2/Z is significant.
# The overlap between d(n) and s(n') depends on n, n', Z.

# Actually, the simplest insight might be:
# The FRACTION of the potential that's nonlinear (and thus affected by sinc)
# is the fraction where the kink field is significant: r < 1/sqrt(Z).
# All electron physics happens at r >> 1/sqrt(Z).
# So the sinc correction is negligible for ALL atoms.

# UNLESS the "sinc" operates in FIELD SPACE, not physical space.
# That is: the coupling between two breather modes goes through the
# nonlinear potential V = (1/pi^2)(1-cos(pi*phi)).
# The effective coupling is d^2V/dphi^2 = cos(pi*phi).
# At phi=0: cos(0) = 1 (linear coupling).
# At phi=1: cos(pi) = -1 (anti-screening!).
# The transition from screening to anti-screening happens at phi=0.5:
# cos(pi*0.5) = 0.

# So the SIGN of the coupling (screening vs anti-screening) depends on
# the local kink field phi. At phi < 0.5: screening. At phi > 0.5: anti.
# This IS what the Oh model captures through w_pi (+0.5) and w_delta (-0.5).

report("FIELD-SPACE INTERPRETATION:")
report("  cos(pi*phi) at the kink profile:")
for phi_val in [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]:
    report(f"    phi = {phi_val:.2f}: cos(pi*phi) = {np.cos(PI*phi_val):+.4f}"
           f"  {'screening' if np.cos(PI*phi_val) > 0 else 'ANTI-screening'}")

report("")
report("The kink profile goes from phi=0 (far) to phi~1-2 (center).")
report("At phi=0.5: coupling switches from screening to anti-screening.")
report("The Oh weights (w_pi=+0.5, w_delta=-0.5) encode this transition.")
report("")
report("The VP_self sinc correction is about the SECOND-ORDER effect:")
report("how the phi^4 nonlinearity modifies the coupling STRENGTH, not sign.")
report("This is an order alpha^2 ~ 5e-5 correction — too small for 5-20% errors.")

report("")
report("FINAL ASSESSMENT:")
report("The remaining V21 outliers (Lu, Pd, La, Ce) are NOT from missing VP.")
report("They are from structural issues in the Oh CHANNEL COUNTING:")
report("  - Lu: f14->d1 three-body weight (unknown Oh fraction)")
report("  - Pd: d10 ionizing d-electron (unique case)")
report("  - La/Ce: f0/f1+d1 transition coupling")
report("These need the multi-electron Hessian (exact eigenvalues) to resolve.")

log.close()
print(f'\nResults saved to: {outfile}')
