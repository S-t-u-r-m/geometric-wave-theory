"""
GWT Decay Rates from Harmonic Instability
==========================================
In GWT, decay = transition from an unstable standing wave to a
more stable set of harmonics. The decay rate should depend on:
  1. Energy available (mass difference = instability depth)
  2. Coupling strength (how easily modes can exchange energy)
  3. Phase space (how many final states are available)
  4. Harmonic mismatch (how far the parent is from a stable mode)

Can we find a universal pattern relating decay rates to wave properties?
"""

import numpy as np
import math

d = 3
gamma_sg = np.pi / (16*np.pi - 2)

def m_fermion(n, p):
    return (16.0/np.pi**2) * np.sin(n*gamma_sg) * np.exp(-16*p/np.pi**2) * 1.2209e22

alpha = d**2 / (2**(d+1) * math.factorial(d+2)**(1.0/(d+1)) * np.pi**((d**2+d-1)/(d+1)))
m_e = m_fermion(16, 32)
m_p = 6 * np.pi**5 * m_e

hbar_MeV_s = 6.582119569e-22  # hbar in MeV*s

# ==============================================================
# TABLE OF UNSTABLE PARTICLES AND THEIR DECAYS
# ==============================================================
# For each: mass, lifetime, dominant decay channel, coupling type
#
# KEY QUESTION: What determines the HUGE range of lifetimes?
# tau ranges from 10^-25 s (strong decay) to 880 s (neutron)
# That's 27 orders of magnitude! Can GWT explain this?

print("=" * 90)
print("DECAY LIFETIMES: SEARCHING FOR HARMONIC PATTERNS")
print("=" * 90)

# Organize by decay type (coupling that mediates the instability)
particles = [
    # WEAK DECAYS (W boson mediates the mode transition)
    # (name, mass_MeV, lifetime_s, Q_MeV, decay_channel, coupling)
    ("neutron",    939.565, 878.4,    1.293,  "p e nu",         "weak"),
    ("muon",       105.658, 2.197e-6, 105.1,  "e nu nu",        "weak"),
    ("tau",        1776.86, 2.903e-13, 1776.3, "hadrons/leptons","weak"),
    ("pi+/-",      139.570, 2.603e-8, 139.0,  "mu nu",          "weak"),
    ("K+/-",       493.677, 1.238e-8, 493.2,  "mu nu / pi pi",  "weak"),
    ("Lambda",    1115.683, 2.631e-10, 176.1,  "p pi-",          "weak"),
    ("Sigma+",    1189.37,  0.802e-10, 189.4,  "p pi0 / n pi+",  "weak"),
    ("Xi-",       1321.71,  1.639e-10, 205.0,  "Lambda pi-",     "weak"),
    ("Omega-",    1672.45,  0.821e-10, 211.0,  "Lambda K-",      "weak"),
    ("D+/-",      1869.66,  1.040e-12, 1869.0, "K* + ...",       "weak"),
    ("B+/-",      5279.34,  1.638e-12, 5279.0, "D* + ...",       "weak"),

    # EM DECAYS (photon mediates)
    ("pi0",        134.977, 8.43e-17,  134.977,"gamma gamma",    "EM"),
    ("eta",        547.862, 5.02e-19,  547.0,  "gamma gamma/3pi","EM"),
    ("Sigma0",    1192.642, 7.4e-20,   76.96,  "Lambda gamma",   "EM"),

    # STRONG DECAYS (gluon/pion mediates)
    ("Delta",     1232.0,   5.63e-24,  294.0,  "N pi",           "strong"),
    ("rho",        775.26,  4.45e-24,  630.0,  "pi pi",          "strong"),
    ("omega",      782.65,  7.75e-23,  643.0,  "3pi",            "strong"),  # isospin violating
    ("phi",       1019.46,  1.55e-22,  632.0,  "K Kbar",         "strong"),  # OZI suppressed
    ("Sigma*",    1385.0,   1.7e-23,   193.0,  "Lambda pi",      "strong"),
]

# Compute log(lifetime) and look for patterns with mass, Q-value, coupling
print(f"\n{'Particle':>10s}  {'Mass':>8s}  {'tau (s)':>12s}  {'Q (MeV)':>8s}  {'log10(tau)':>10s}  {'Coupling':>8s}")
print("-" * 90)

for name, mass, tau, Q, channel, coupling in particles:
    print(f"{name:>10s}  {mass:8.1f}  {tau:12.3e}  {Q:8.1f}  {np.log10(tau):10.2f}  {coupling:>8s}")


# ==============================================================
# PATTERN SEARCH: What determines the lifetime?
# ==============================================================
print("\n\n" + "=" * 90)
print("PATTERN ANALYSIS")
print("=" * 90)

# HYPOTHESIS 1: Lifetime scales with coupling strength
# Weak: G_F ~ 10^-5 GeV^-2 → tau ~ 1/(G_F^2 * Q^5) (Fermi theory)
# EM: alpha ~ 1/137 → tau ~ 1/(alpha^2 * Q)
# Strong: alpha_s ~ 1 → tau ~ 1/Q (immediate)
#
# In GWT: all couplings come from d=3. The coupling hierarchy IS the
# harmonic hierarchy — weak transitions require tunneling through MORE
# barriers (higher p in the m(n,p) formula), EM through fewer, strong = direct.
#
# The key GWT insight: the TUNNELING FACTOR T^2 = exp(-16/pi^2) = 0.1977
# sets the coupling hierarchy.
#
# Weak coupling ~ T^4 (two extra barriers to cross)
# EM coupling ~ T^2 (one barrier)
# Strong coupling ~ T^0 = 1 (no barrier)

T2 = np.exp(-16/np.pi**2)  # = 0.1977
print(f"\n  GWT tunneling factor: T^2 = exp(-16/pi^2) = {T2:.4f}")
print(f"  T^4 = {T2**2:.6f}")
print(f"  T^8 = {T2**4:.10f}")

# HYPOTHESIS 2: Universal decay rate formula
# In GWT, a standing wave decays when it can tunnel to a lower-energy
# configuration. The rate depends on:
#   Gamma = (coupling)^2 * Q^(2L+1) * (phase_space) / mass^(2L)
#
# where L = angular momentum change (harmonic selection rule)
# Q = energy release (instability depth)
# coupling = alpha (EM), G_F*Q (weak), alpha_s (strong)
#
# Can we unify all decays as:
#   Gamma = (T^2)^(2*n_barriers) * Q^(2L+1) / (hbar * mass^(2L))
# where n_barriers = number of tunneling barriers?

print(f"\n  --- Testing: log(Gamma) vs barrier count ---")
print(f"\n{'Particle':>10s}  {'Gamma(MeV)':>12s}  {'Q/M':>8s}  {'Barriers':>8s}  {'T^(2*barr)*Q':>14s}")
print("-" * 90)

for name, mass, tau, Q, channel, coupling in particles:
    Gamma = hbar_MeV_s / tau  # decay width in MeV
    Q_over_M = Q / mass

    # Assign barrier count from coupling type:
    # Strong = 0 barriers (direct harmonic transition)
    # EM = 2 barriers (one for each photon vertex: alpha = T^2 * geometric)
    # Weak = 4 barriers (W boson exchange crosses 4 barriers)
    if coupling == "strong":
        n_barr = 0
    elif coupling == "EM":
        n_barr = 2
    else:  # weak
        n_barr = 4

    predicted_scale = T2**(n_barr) * Q  # crude prediction for Gamma
    ratio = Gamma / predicted_scale if predicted_scale > 0 else 0

    print(f"{name:>10s}  {Gamma:12.3e}  {Q_over_M:8.4f}  {n_barr:8d}  {predicted_scale:14.4e}  ratio={ratio:.3e}")


# ==============================================================
# DEEPER PATTERN: Fermi theory and phase space
# ==============================================================
print("\n\n" + "=" * 90)
print("WEAK DECAYS: FERMI THEORY PATTERN")
print("=" * 90)

# For weak decays, the standard formula is:
# Gamma ~ G_F^2 * Q^5 / (192*pi^3)  (for 3-body leptonic)
# Gamma ~ G_F^2 * Q^5 * |V_CKM|^2 * (1 + corrections)
#
# In GWT: G_F = 1/(sqrt(2)*v^2) where v = Higgs VEV = m(3,23)
# So G_F is already derived from d=3.
#
# The Q^5 comes from the 3-body phase space (dimensional analysis).
# In GWT: this is the density of final-state modes available = Q^(2d-1) = Q^5

G_F = 1.0 / (np.sqrt(2) * (m_fermion(3,23))**2)  # MeV^-2

print(f"\n  G_F (GWT) = {G_F:.4e} MeV^-2")
print(f"\n  Fermi scaling: Gamma ~ G_F^2 * Q^5 / (192*pi^3)")
print(f"\n{'Particle':>10s}  {'tau obs':>12s}  {'Q^5':>14s}  {'tau_Fermi':>12s}  {'Ratio':>8s}")
print("-" * 90)

for name, mass, tau, Q, channel, coupling in particles:
    if coupling != "weak":
        continue
    # Crude Fermi estimate: Gamma = G_F^2 * Q^5 / (192*pi^3)
    Gamma_Fermi = G_F**2 * Q**5 / (192 * np.pi**3)
    tau_Fermi = hbar_MeV_s / Gamma_Fermi if Gamma_Fermi > 0 else np.inf
    ratio = tau / tau_Fermi

    print(f"{name:>10s}  {tau:12.3e}  {Q**5:14.3e}  {tau_Fermi:12.3e}  {ratio:8.2f}")


# ==============================================================
# GWT HARMONIC ANALYSIS: What makes a mode unstable?
# ==============================================================
print("\n\n" + "=" * 90)
print("HARMONIC STABILITY ANALYSIS")
print("=" * 90)

# In GWT, the breather spectrum is m(n,p) = (16/pi^2)*sin(n*gamma)*exp(-16p/pi^2)*m_P
# A particle is STABLE if there's no combination of lighter breathers
# that has the same quantum numbers and lower total energy.
#
# The neutron is unstable because:
# n(udd) -> p(uud) + e + nu
# This is a d-quark -> u-quark transition: breather n=5 -> breather n=13
# The harmonic mismatch: Delta_n = 13 - 5 = 8 = 2^d
# And Delta_p = 31 - 30 = 1
#
# The muon is unstable because:
# mu -> e + nu + nu
# breather n=4,p=28 -> breather n=16,p=32
# Delta_n = 16 - 4 = 12 = d*2^(d-1) = N/2
# Delta_p = 32 - 28 = 4 = 2^(d-1)
#
# The pion is unstable because:
# pi -> mu + nu
# The pion is a qq-bar pair, not a single breather.
# Its mass ~ sqrt(m_u + m_d) * condensate scale
#
# PATTERN: How many "harmonic steps" does each decay require?

print(f"\n  Key question: Is the decay rate related to the harmonic")
print(f"  distance (Delta_n, Delta_p) between initial and final modes?")

print(f"\n  Neutron -> Proton:")
print(f"    d-quark (n=5, p=30) -> u-quark (n=13, p=31)")
print(f"    Delta_n = 8 = 2^d")
print(f"    Delta_p = 1")
print(f"    Tunneling cost: T^(2*Delta_p) = T^2 = {T2:.4f}")

print(f"\n  Muon -> Electron:")
print(f"    muon (n=4, p=28) -> electron (n=16, p=32)")
print(f"    Delta_n = 12 = N/2 = d*2^(d-1)")
print(f"    Delta_p = 4 = 2^(d-1)")
print(f"    Tunneling cost: T^(2*Delta_p) = T^8 = {T2**4:.6e}")

print(f"\n  Tau -> (hadrons/electron):")
print(f"    tau (n=18, p=27) -> e/mu + hadrons")
print(f"    Delta_p (tau->e) = 5")
print(f"    Tunneling cost: T^(2*5) = T^10 = {T2**5:.6e}")

# Now test: does tau ~ 1/(T^(2*Delta_p) * Q^5 * G_F^2)?
print(f"\n\n  --- Testing: tau vs T^(2*Delta_p) ---")
print(f"\n{'Decay':>15s}  {'Delta_p':>8s}  {'T^(2dp)':>12s}  {'tau_pred':>12s}  {'tau_obs':>12s}  {'Ratio':>8s}")
print("-" * 90)

# Neutron: d->u, Delta_p = 1
Q_n = 1.293
dp_n = 1
Gamma_n = G_F**2 * Q_n**5 / (192*np.pi**3) * T2**(2*dp_n)
tau_n = hbar_MeV_s / Gamma_n
print(f"{'n -> p e nu':>15s}  {dp_n:8d}  {T2**(2*dp_n):12.4e}  {tau_n:12.3e}  {878.4:12.3e}  {tau_n/878.4:8.2f}")

# Muon: mu->e, Delta_p = 4
Q_mu = 105.1
dp_mu = 4
Gamma_mu = G_F**2 * Q_mu**5 / (192*np.pi**3) * T2**(2*dp_mu)
tau_mu = hbar_MeV_s / Gamma_mu
print(f"{'mu -> e nu nu':>15s}  {dp_mu:8d}  {T2**(2*dp_mu):12.4e}  {tau_mu:12.3e}  {2.197e-6:12.3e}  {tau_mu/2.197e-6:8.2f}")

# Tau: tau->X, Delta_p = 5 (to electron)
Q_tau = 1776.3
dp_tau = 5  # tau p=27 to electron p=32
Gamma_tau = G_F**2 * Q_tau**5 / (192*np.pi**3) * T2**(2*dp_tau)
tau_tau = hbar_MeV_s / Gamma_tau
print(f"{'tau -> X':>15s}  {dp_tau:8d}  {T2**(2*dp_tau):12.4e}  {tau_tau:12.3e}  {2.903e-13:12.3e}  {tau_tau/2.903e-13:8.2f}")

# Pion: pi -> mu nu. The pion is not a single breather.
# But its decay involves a quark annihilation: u + dbar -> W -> mu + nu
# The quark transition: same as creating a mu (p=28) from quarks (p~30)
# Delta_p ~ 2?
Q_pi = 139.0
dp_pi = 2
Gamma_pi = G_F**2 * Q_pi**5 / (192*np.pi**3) * T2**(2*dp_pi)
tau_pi = hbar_MeV_s / Gamma_pi
print(f"{'pi -> mu nu':>15s}  {dp_pi:8d}  {T2**(2*dp_pi):12.4e}  {tau_pi:12.3e}  {2.603e-8:12.3e}  {tau_pi/2.603e-8:8.2f}")

# Kaon: K -> mu nu. Strange quark s(n=4,p=28) -> u(n=13,p=31), Delta_p = 3
Q_K = 493.2
dp_K = 3  # s(p=28) to u(p=31)
V_us = 0.2243  # CKM suppression
Gamma_K = G_F**2 * Q_K**5 / (192*np.pi**3) * T2**(2*dp_K) * V_us**2
tau_K = hbar_MeV_s / Gamma_K
print(f"{'K -> mu nu':>15s}  {dp_K:8d}  {T2**(2*dp_K):12.4e}  {tau_K:12.3e}  {1.238e-8:12.3e}  {tau_K/1.238e-8:8.2f}")


# ==============================================================
# DOES THE TUNNELING FACTOR REPLACE G_F?
# ==============================================================
print("\n\n" + "=" * 90)
print("KEY QUESTION: Is G_F^2 itself a tunneling factor?")
print("=" * 90)

# In GWT: G_F = 1/(sqrt(2)*v^2) where v = Higgs VEV
# v = m(3,23) = (16/pi^2)*sin(3*gamma)*exp(-16*23/pi^2)*m_Planck
# The exp(-16*23/pi^2) factor IS a tunneling factor with p=23.
#
# G_F ~ 1/v^2 ~ exp(+2*16*23/pi^2) / m_Planck^2
# = exp(+32*23/pi^2) / m_Planck^2
# = T^(-2*23) / m_Planck^2
#
# So G_F^2 ~ T^(-4*23) / m_Planck^4 = T^{-92} / m_Planck^4
#
# And the decay rate involves:
# Gamma ~ G_F^2 * Q^5 * T^(2*Delta_p)
#       = T^{-92} * Q^5 * T^(2*Delta_p) / m_Planck^4
#       = T^{2*Delta_p - 92} * Q^5 / m_Planck^4
#
# The TOTAL tunneling exponent is 2*Delta_p - 92.
# For neutron: 2*1 - 92 = -90 (deeper tunneling = slower decay)
# For muon: 2*4 - 92 = -84 (shallower = faster)
# For tau: 2*5 - 92 = -82 (even shallower = even faster)
#
# Wait — this means ALL weak decays have the SAME overall tunneling structure,
# and the difference is just how many barriers the decay mode must cross
# relative to the G_F "background" of 92 barriers.

v_MeV = m_fermion(3, 23)
p_v = 23  # Higgs VEV tunneling depth
print(f"\n  Higgs VEV: v = m(3, {p_v}) = {v_MeV:.0f} MeV")
print(f"  G_F = 1/(sqrt(2)*v^2) contains tunneling T^(-2*{p_v}) = T^(-{2*p_v})")
print(f"  G_F^2 contains T^(-{4*p_v})")
print(f"\n  Total tunneling exponent for decay = 2*Delta_p - {4*p_v}")
print(f"\n  Neutron (dp=1): total = {2*1 - 4*p_v}")
print(f"  Muon (dp=4):    total = {2*4 - 4*p_v}")
print(f"  Tau (dp=5):     total = {2*5 - 4*p_v}")
print(f"  Pion (dp=2):    total = {2*2 - 4*p_v}")
print(f"\n  All have similar total tunneling ~ -90. The Q^5 phase space")
print(f"  then determines the RELATIVE rates between different decays.")

# Let's test: strip out G_F and just use T^(2*p_total) * Q^5
# where p_total = p(initial) - p(final) for the quark transition
# Neutron: d(p=30) -> u(p=31): p_total = 30 (initial quark depth)
# Actually, the DECAY rate should involve the DIFFERENCE between
# initial and final configurations in the tunneling landscape.

print("\n\n  --- PURE TUNNELING MODEL ---")
print(f"\n  Rate = C * T^(2*p_initial) * Q^(2d-1) / m_Planck^4")
print(f"  where C is a geometric factor from d=3")
print(f"  and p_initial = tunneling depth of the decaying quark")
print(f"\n  Testing against weak decay lifetimes:")
print(f"\n{'Decay':>15s}  {'p_init':>6s}  {'Q^5 (MeV^5)':>14s}  {'T^(2p)*Q^5':>14s}")
print("-" * 70)

# The idea: the decay rate is set by how deep the decaying quark sits.
# Deeper quarks = more barriers = slower decay.
decays = [
    ("n->p e nu",    30, 1.293,   878.4),      # d-quark at p=30
    ("mu->e nu nu",  28, 105.1,   2.197e-6),    # muon at p=28
    ("tau->X",       27, 1776.3,  2.903e-13),   # tau at p=27
    ("pi->mu nu",    30, 139.0,   2.603e-8),    # d-quark (in pion) at p=30
    ("K->mu nu",     28, 493.2,   1.238e-8),    # s-quark at p=28
    ("D->K+",        27, 1869.0,  1.040e-12),   # c-quark at p=27
    ("B->D+",        26, 5279.0,  1.638e-12),   # b-quark at p=26
]

for name, p_init, Q, tau_obs in decays:
    rate_scale = T2**(2*p_init) * Q**5
    tau_scale = hbar_MeV_s / rate_scale
    print(f"{name:>15s}  {p_init:6d}  {Q**5:14.3e}  {rate_scale:14.3e}  tau_scale={tau_scale:.3e}  obs={tau_obs:.3e}")

print(f"\n  Note: tau_scale is NOT the actual lifetime (missing geometric prefactors).")
print(f"  But the RATIOS should show if the pattern works.")

# Compare ratios of lifetimes
print(f"\n\n  --- LIFETIME RATIOS (model vs observed) ---")
ref_name, ref_p, ref_Q, ref_tau = decays[0]  # neutron as reference
ref_scale = T2**(2*ref_p) * ref_Q**5

for name, p_init, Q, tau_obs in decays:
    scale = T2**(2*p_init) * Q**5
    ratio_model = ref_scale / scale  # tau_ref / tau_i model
    ratio_obs = ref_tau / tau_obs  # tau_ref / tau_i observed
    print(f"  tau_n / tau_{name:>12s}: model = {ratio_model:12.3e}  obs = {ratio_obs:12.3e}  (model/obs = {ratio_model/ratio_obs:.2f})")


# ==============================================================
# NEUTRON LIFETIME FROM PURE GWT
# ==============================================================
print("\n\n" + "=" * 90)
print("NEUTRON LIFETIME: PURE GWT APPROACH")
print("=" * 90)

# If we can express the neutron lifetime in terms of the muon lifetime
# (which is well-measured and well-predicted by Fermi theory), we can
# bypass g_A entirely.
#
# Both are weak decays. The ratio should be:
# tau_n / tau_mu = (Q_mu / Q_n)^5 * T^(2*(p_d - p_mu)) * (CKM & spin factors)
#
# p_d = 30 (d-quark), p_mu = 28 (muon)
# Q_n = 1.293 MeV, Q_mu = 105.1 MeV
#
# tau_n / tau_mu = (105.1/1.293)^5 * T^(2*(30-28)) * (1/V_ud^2) * (1/(1+3*g_A^2))
#
# Hmm, g_A appears again because the neutron is a composite (spin matters).
# The muon is elementary (spin is simple).
#
# But in GWT, the spin factor IS a geometric factor from the standing wave.
# For an elementary fermion (muon): the spin factor = 1
# For a composite (neutron): the spin factor = (1 + 3*g_A^2)
# In the quark model: g_A = 5/3 -> 1 + 3*(5/3)^2 = 28/3 = 9.33
# This is a pure SU(6) number from d=3.
# In GWT: the proton is a j_0 mode with 3 quarks. The (1+3*g_A^2) factor
# counts the number of helicity states that contribute.
#
# Actually, (1 + 3*g_A^2) has a simple meaning:
# 1 = Fermi (vector) contribution
# 3*g_A^2 = Gamow-Teller (axial) contribution
# The factor of 3 comes from d=3 spin orientations.
# g_A^2 is the axial form factor squared.
#
# For the bare SU(6): 1 + 3*(5/3)^2 = 1 + 25/3 = 28/3
# For cavity (g_A=5/9): 1 + 3*(5/9)^2 = 1 + 25/27 = 52/27 = 1.926

# Let me compute the neutron lifetime as a RATIO to the muon lifetime.
# This way, many common factors (G_F, etc.) cancel.

tau_mu_obs = 2.1970e-6  # muon lifetime (seconds)
Q_mu = 105.658  # muon Q-value (mass, since m_e << m_mu)
Q_n = 1.293  # neutron Q-value (observed mn-mp)

# tau_n / tau_mu = (Q_mu/Q_n)^5 * (m_mu/m_e)^5 correction * (spin factors)
# Actually, more carefully:
# Gamma_mu = G_F^2 * m_mu^5 / (192*pi^3)  (exactly)
# Gamma_n = G_F^2 * V_ud^2 * m_e^5 * f * (1+3*g_A^2) / (2*pi^3)
# where f is the phase space integral (depends on Q_n/m_e)
#
# Ratio: tau_n/tau_mu = Gamma_mu/Gamma_n
# = [m_mu^5 / (192*pi^3)] / [V_ud^2 * m_e^5 * f * (1+3*g_A^2) / (2*pi^3)]
# = m_mu^5 / (96 * V_ud^2 * m_e^5 * f * (1+3*g_A^2))
#
# Everything here is GWT-derivable EXCEPT g_A.
# Let's see what g_A would NEED to be to give the right tau_n.

from scipy import integrate

def fermi_function(Z, E):
    p = np.sqrt(E**2 - 1)
    eta = 2 * np.pi * alpha * Z * E / p
    return eta / (1 - np.exp(-eta))

def phase_space_integrand(E, E0):
    if E >= E0 or E <= 1:
        return 0
    p = np.sqrt(E**2 - 1)
    return fermi_function(1, E) * p * E * (E0 - E)**2

m_e_obs = 0.51100
E0_obs = Q_n / m_e_obs
f_ps, _ = integrate.quad(phase_space_integrand, 1.0, E0_obs, args=(E0_obs,))

# What g_A gives tau_n = 878.4 s?
# tau_n = tau_mu * m_mu^5 / (96 * V_ud^2 * m_e^5 * f * (1+3*g_A^2))
# (1+3*g_A^2) = tau_mu * m_mu^5 / (96 * V_ud^2 * m_e^5 * f * tau_n)

V_ud_obs = 0.97373
spin_factor_needed = tau_mu_obs * Q_mu**5 / (96 * V_ud_obs**2 * m_e_obs**5 * f_ps * 878.4)
g_A_needed = np.sqrt((spin_factor_needed - 1) / 3)

print(f"\n  Muon-to-neutron ratio approach:")
print(f"  tau_mu = {tau_mu_obs:.4e} s")
print(f"  Q_mu = {Q_mu:.1f} MeV, Q_n = {Q_n:.3f} MeV")
print(f"  Phase space f = {f_ps:.4f}")
print(f"  (1+3*g_A^2) needed for tau_n = 878.4 s: {spin_factor_needed:.4f}")
print(f"  g_A needed: {g_A_needed:.4f}")
print(f"  Observed g_A: 1.2756")

# What if we use g_A = 5(d-1)/(6d) = 5/9 (bare cavity)?
spin_factor_cavity = 1 + 3*(5.0/9)**2
tau_n_cavity = tau_mu_obs * Q_mu**5 / (96 * V_ud_obs**2 * m_e_obs**5 * f_ps * spin_factor_cavity)
print(f"\n  With cavity g_A = 5/9:")
print(f"  (1+3*g_A^2) = {spin_factor_cavity:.4f}")
print(f"  tau_n = {tau_n_cavity:.1f} s  (obs: 878.4)")
print(f"  Error: {(tau_n_cavity-878.4)/878.4*100:+.1f}%")

# With SU(6) g_A = 5/3:
spin_factor_su6 = 1 + 3*(5.0/3)**2
tau_n_su6 = tau_mu_obs * Q_mu**5 / (96 * V_ud_obs**2 * m_e_obs**5 * f_ps * spin_factor_su6)
print(f"\n  With SU(6) g_A = 5/3:")
print(f"  (1+3*g_A^2) = {spin_factor_su6:.4f}")
print(f"  tau_n = {tau_n_su6:.1f} s  (obs: 878.4)")
print(f"  Error: {(tau_n_su6-878.4)/878.4*100:+.1f}%")

# Interesting: the spin factor matters A LOT.
# With g_A=5/9: 1+3*(25/81) = 1 + 75/81 = 1.926 → tau ~ 2300 s
# With g_A=5/3: 1+3*(25/9) = 1 + 75/9 = 9.333 → tau ~ 476 s
# With g_A=1.276: 1+3*1.628 = 5.885 → tau ~ 878 s
# So the answer is BETWEEN the two geometric values 5/9 and 5/3.
# In fact: 1+3*g_A^2 = 5.885 ≈ 6 ≈ 2d
# Is this a coincidence?

print(f"\n  Note: 1+3*g_A_obs^2 = {1+3*1.2756**2:.3f} ~ {2*d} = 2d")
print(f"  If (1+3*g_A^2) = 2d exactly: g_A = sqrt((2d-1)/3) = sqrt(5/3) = {np.sqrt(5/3):.4f}")
print(f"  Error vs observed: {(np.sqrt(5/3)-1.2756)/1.2756*100:+.1f}%")

g_A_geometric = np.sqrt((2*d-1)/3)
spin_factor_geo = 1 + 3*g_A_geometric**2  # = 2d = 6 exactly
tau_n_geo = tau_mu_obs * Q_mu**5 / (96 * V_ud_obs**2 * m_e_obs**5 * f_ps * spin_factor_geo)
print(f"\n  With g_A = sqrt((2d-1)/3) = sqrt(5/3):")
print(f"  (1+3*g_A^2) = {spin_factor_geo:.4f} = 2d = {2*d}")
print(f"  tau_n = {tau_n_geo:.1f} s  (obs: 878.4)")
print(f"  Error: {(tau_n_geo-878.4)/878.4*100:+.1f}%")

# WOW. If g_A = sqrt(5/3) = sqrt((2d-1)/3), we get:
# 1+3*g_A^2 = 1 + (2d-1) = 2d = 6 (exactly at d=3)
# This is the number of faces of a d-cube!
# Physical interpretation: the total spin response of the neutron
# is 2d = 6 (one per face of the unit cell in the lattice).
# The Fermi (vector) part is 1 (the scalar), and the GT (axial) part
# is 2d-1 = 5 (the (2d-1) transverse directions on the d-cube faces).
