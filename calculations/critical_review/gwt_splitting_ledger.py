#!/usr/bin/env python3
"""
GWT mass-splitting ledger + dimensional-rule unification test.

A splitting is a difference of two sectors of the SAME standing wave, so the
common-mode bare scale (and its derivation noise) cancels, leaving the
correction structure exposed against a clean observable.

Hypothesis (the "one geometry rule"): for each splitting, the correction's
COUPLING TYPE, SIGN, and COEFFICIENT are all fixed by which axis / winding
the transition activates -- derived from geometry, not fit to the gap.

This script decomposes each GWT splitting as  base x (1 +/- c*g)  and asks
whether one axis-indexed scheme accounts for (g, sign, c) across all of them.
"""
import numpy as np

PI = np.pi
d  = 3
alpha   = 1/137.036
alpha_s = 0.11794
m_e = 0.510999
m_p = 938.2720

omega0 = m_p * (2*d*PI) / (2*d*PI**2 + 1)          # fundamental angular quantum
strange_base_frac = 3*PI / (6*PI**2 + 1)            # 3pi/(6pi^2+1)

# ---- observed isospin-averaged splittings (PDG), MeV ----
obs = {
    'n - p            (charge change, light)': 939.5654 - 938.2720,
    'Delta - p        (spin flip, light)'    : 1232.0   - 938.2720,
    'Sigma* - Sigma   (spin flip, 1 strange)': 1384.6   - 1193.2,
    'Xi* - Xi         (spin flip, 2 strange)': 1533.4   - 1318.3,
    'Sigma - Lambda   (isospin recouple)'    : 1192.64  - 1115.68,  # out-of-sample
}

# ---- GWT predictions + decomposition (base, coupling g, sign, coeff c) ----
# each entry: (predicted, base_value, base_expr, g_name, g_val, sign, c_name, c_val)
rows = []

# n-p : EM correction, charge axis (U(1) phase mode)
base = m_e * 8/3
rows.append(('n - p            (charge change, light)',
             base*(1 - 7*alpha), base, 'm_e * 8/3',
             'alpha (EM, U(1))', alpha, '-', '2d+1', 7))

# Delta-p : pure geometric angular quantum, rotational axis (SU(2)); no coupling corr
rows.append(('Delta - p        (spin flip, light)',
             omega0, omega0, 'omega_0',
             '(none, geometric)', 0.0, '', '-', 0))

# Sigma*-Sigma : strong correction, strange-color axis (SU(3))
base = 1193.2 * strange_base_frac
rows.append(('Sigma* - Sigma   (spin flip, 1 strange)',
             base*(1 + alpha_s/(d+1)), base, 'm_Sigma * 3pi/(6pi^2+1)',
             'alpha_s (strong, SU(3))', alpha_s, '+', 'd+1', 4))

# Xi*-Xi : same form, heavier baryon
base = 1318.3 * strange_base_frac
rows.append(('Xi* - Xi         (spin flip, 2 strange)',
             base*(1 + alpha_s/(d+1)), base, 'm_Xi * 3pi/(6pi^2+1)',
             'alpha_s (strong, SU(3))', alpha_s, '+', 'd+1', 4))

# Sigma-Lambda : OUT OF SAMPLE -- GWT has no published formula. Do NOT fit.
rows.append(('Sigma - Lambda   (isospin recouple)',
             np.nan, np.nan, '(no GWT formula)',
             '?', np.nan, '?', '?', np.nan))

print('='*92)
print('GWT SPLITTING LEDGER  --  base x (1 +/- c*g)  vs observed')
print('='*92)
hdr = f'{"splitting":40} {"obs":>8} {"GWT":>9} {"err%":>7}  {"coupling g":>22} {"sign":>4} {"coeff c":>6}'
print(hdr); print('-'*92)
for name, pred, base, bexpr, gname, gval, sign, cname, cval in rows:
    o = obs[name]
    if np.isnan(pred):
        print(f'{name:40} {o:8.2f} {"--":>9} {"--":>7}  {gname:>22} {str(sign):>4} {str(cname):>6}')
    else:
        err = 100*(pred-o)/o
        print(f'{name:40} {o:8.2f} {pred:9.3f} {err:7.3f}  {gname:>22} {sign:>4} {cname:>6}')
print('-'*92)

print()
print('UNIFICATION TEST -- does ONE axis-indexed rule set (coupling, sign, coeff)?')
print('-'*92)
print("""
  transition         activated axis        ->  coupling     sign    coeff
  charge change      U(1)  phase mode           alpha          -     2d+1 = 7
  light spin flip    SU(2) rotational           (geometric)    n/a   --
  strange spin flip  SU(3) color (heavy sub)    alpha_s        +     d+1  = 4

  VERDICT:
   [WIN]     Coupling TYPE follows the activated axis cleanly:
             charge->U(1)->alpha,  strange-spin->SU(3)->alpha_s,
             light-spin->SU(2)->pure geometric (omega_0, no coupling term).
             This is a real structural unification and it is exactly what the
             'sectors of one wave' ontology predicts -- the gauge axis that the
             transition touches dictates which coupling renormalizes it.

   [WIN]     SIGN follows the vacuum sign rule: charge change is particle-like
             (loses energy -> minus); mass-raising spin flip gains -> plus.

   [OPEN]    The BASE factor is NOT unified: m_e*8/3, omega_0, m_B*3pi/(6pi^2+1)
             are three different geometric expressions, one per transition class.

   [OPEN]    The COEFFICIENTS (2d+1=7 vs d+1=4) are different integer counts with
             individually-plausible labels ('exchange paths' vs 'spacetime divisor'),
             NOT yet output by one axis-counting scheme. This is the same gap as
             vacuum_corrections.md Open Q#5, now localized: the *selection principle*
             that maps an activated axis to its integer coefficient is what's missing.

  STATUS: principle-level unification YES, quantitative one-formula rule NOT YET.
""")

print('-'*92)
print('THE TEST THAT WOULD CARRY IT (out of sample):  Sigma - Lambda = %.2f MeV'
      % obs['Sigma - Lambda   (isospin recouple)'])
print("""  Same uds content, no charge change, no strangeness change -- a pure light-quark
  isospin/spin recoupling. The axis-rule says it must be an SU(2)/SU(3) light-sector
  effect with a sign and coefficient PREDICTABLE from the same scheme. GWT has no
  formula for it yet. Deriving it from the axis-counting WITHOUT fitting -- and
  landing near 77 MeV -- is the clean, above-noise, real-observable win the whole
  dimensional idea needs. If the scheme can't produce it, that's the honest limit.""")
print('='*92)
