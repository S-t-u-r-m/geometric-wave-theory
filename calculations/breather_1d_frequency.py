"""
1D Exact Breather Frequency Analysis (CPU/NumPy)
=================================================
Measures oscillation frequency of sine-Gordon breathers on a 1D discrete lattice.
1D breathers are EXACT solutions (unlike 3D which radiate into phonons).
The ratio omega_discrete / omega_continuous gives the lattice correction factor.
"""

import numpy as np
import time

gamma_gwt = np.pi / (16 * np.pi - 2)
m_Planck_MeV = 1.2209e22

n_values = [4, 5, 7, 11, 12, 13, 16, 18]
particle_names = {
    4: "muon/strange", 5: "down", 7: "bottom", 11: "charm",
    12: "top", 13: "up", 16: "electron", 18: "tau"
}

N_1d = 4096
a = 1.0
dt = 0.001
N_steps = 200000
SKIP_TRANSIENT = 5000
RECORD_EVERY = 2

print(f"1D Lattice: {N_1d} sites, dt={dt}, {N_steps} steps")
n_recorded = (N_steps - SKIP_TRANSIENT) // RECORD_EVERY
dt_sample = dt * RECORD_EVERY
print(f"Samples: {n_recorded}, freq resolution: {1.0/(n_recorded*dt_sample):.6f}")
print(f"Nyquist: {1.0/(2*dt_sample):.1f}")
print()


def make_breather_1d(N_grid, n):
    center = N_grid // 2
    omega = np.cos(n * gamma_gwt)
    eta = np.sin(n * gamma_gwt)
    x = np.arange(N_grid, dtype=np.float64) - center
    phi = 4.0 * np.arctan(eta / (omega * np.cosh(np.minimum(eta * np.abs(x), 50))))
    phi_t = np.zeros(N_grid, dtype=np.float64)
    return phi, phi_t, center


def accel_1d(phi, a_sp):
    return (np.roll(phi, 1) + np.roll(phi, -1) - 2*phi) / a_sp**2 - np.sin(phi)


def evolve_and_record_1d(phi, phi_t, dt_step, a_sp, n_steps, center_idx,
                          skip=5000, rec_every=2):
    series = []
    times = []
    acc = accel_1d(phi, a_sp)
    for step in range(n_steps):
        phi_t += 0.5 * dt_step * acc
        phi += dt_step * phi_t
        acc = accel_1d(phi, a_sp)
        phi_t += 0.5 * dt_step * acc
        if step >= skip and step % rec_every == 0:
            series.append(phi[center_idx])
            times.append(step * dt_step)
    return np.array(series), np.array(times)


def extract_frequency(signal, times):
    signal = signal - np.mean(signal)
    window = np.hanning(len(signal))
    windowed = signal * window
    fft_vals = np.fft.rfft(windowed)
    fft_power = np.abs(fft_vals)**2
    dt_s = times[1] - times[0]
    freqs = np.fft.rfftfreq(len(windowed), d=dt_s)
    peak_idx = np.argmax(fft_power[1:]) + 1
    peak_freq = freqs[peak_idx]
    if 1 < peak_idx < len(fft_power) - 1:
        al = np.log(fft_power[peak_idx - 1] + 1e-30)
        be = np.log(fft_power[peak_idx] + 1e-30)
        gp = np.log(fft_power[peak_idx + 1] + 1e-30)
        delta_p = 0.5 * (al - gp) / (al - 2*be + gp)
        peak_freq = freqs[peak_idx] + delta_p * (freqs[1] - freqs[0])
    return 2 * np.pi * peak_freq


print("=" * 100)
print("1D EXACT BREATHER FREQUENCY ANALYSIS")
print("=" * 100)
print()
hdr = '  n        Particle        w_cont     w_1D_disc       Ratio    Correction    Time' 
print(hdr)
print("=" * 100)

results_1d = []
for n_val in n_values:
    name = particle_names[n_val]
    omega_cont = np.cos(n_val * gamma_gwt)
    t0 = time.time()
    phi, phi_t, center = make_breather_1d(N_1d, n_val)
    series, times = evolve_and_record_1d(phi, phi_t, dt, a, N_steps, center, skip=SKIP_TRANSIENT, rec_every=RECORD_EVERY)
    omega_disc = extract_frequency(series, times)
    ratio = omega_disc / omega_cont
    correction = (ratio - 1.0) * 100
    elapsed = time.time() - t0
    print(f"{n_val:3d}  {name:>14s}  {omega_cont:12.6f}  {omega_disc:12.6f}  {ratio:10.6f}  {correction:+11.4f}%  {elapsed:5.1f}s")
    results_1d.append({"n": n_val, "name": name, "omega_cont": omega_cont, "omega_disc": omega_disc, "ratio": ratio, "correction": correction})

gamma_sg = np.pi / (16 * np.pi - 2)
fermion_data = [
    (16, 32, "electron", 0.511),
    (13, 31, "up",       2.16),
    (5,  30, "down",     4.67),
    (4,  28, "muon",     105.66),
    (4,  28, "strange",  93.4),
    (11, 27, "charm",    1271),
    (18, 27, "tau",      1776.86),
    (7,  26, "bottom",   4183),
    (12, 24, "top",      172760),
]

print()
print("=" * 100)
print("CORRECTED PHYSICAL MASSES WITH 1D DISCRETENESS CORRECTIONS")
print("=" * 100)
print()
print(f"{'Particle':>10s}  {'n':>3s} {'p':>3s}  {'m_cont(MeV)':>14s}  {'m_1D(MeV)':>14s}  {'Observed':>12s}  {'Err_cont':>10s}  {'Err_1D':>10s}")
print("-" * 100)

total_err_cont = 0
total_err_1d = 0
count = 0
for n_val, p_val, name, obs_MeV in fermion_data:
    m_cont = (16.0/np.pi**2) * np.sin(n_val * gamma_sg) * np.exp(-16*p_val/np.pi**2) * m_Planck_MeV
    ratio = 1.0
    for r in results_1d:
        if r["n"] == n_val:
            ratio = r["ratio"]
            break
    m_1d = m_cont * ratio
    err_cont = (m_cont - obs_MeV) / obs_MeV * 100
    err_1d = (m_1d - obs_MeV) / obs_MeV * 100
    print(f"{name:>10s}  {n_val:3d} {p_val:3d}  {m_cont:14.4f}  {m_1d:14.4f}  {obs_MeV:12.4f}  {err_cont:+9.2f}%  {err_1d:+9.2f}%")
    total_err_cont += abs(err_cont)
    total_err_1d += abs(err_1d)
    count += 1
print(f"  Mean |error|:  continuous = {total_err_cont/count:.2f}%,  1D-corrected = {total_err_1d/count:.2f}%")

# Neutrino implications
print()
print("=" * 100)
print("ELECTRON CORRECTION -> NEUTRINO IMPLICATIONS")
print("=" * 100)
e_ratio = 1.0
for r in results_1d:
    if r["n"] == 16:
        e_ratio = r["ratio"]
        e_corr = r["correction"]
        break
d = 3
m_e_cont = (16.0/np.pi**2) * np.sin(16 * gamma_sg) * np.exp(-16*32/np.pi**2) * m_Planck_MeV
m_e_1d = m_e_cont * e_ratio
m_p_MeV = 938.272
print(f"  Electron mass (continuous): {m_e_cont:.4f} MeV (error: {(m_e_cont-0.511)/0.511*100:+.2f}%)")
print(f"  Electron mass (1D disc):    {m_e_1d:.4f} MeV (error: {(m_e_1d-0.511)/0.511*100:+.2f}%)")
print(f"  Frequency correction: {e_corr:+.4f}%")
M_nu_cont = m_e_cont**3 / (d * m_p_MeV**2)
M_nu_1d = m_e_1d**3 / (d * m_p_MeV**2)
wyler = 1 + 1/(d * 2 * np.pi**2)
M_eff_cont_meV = M_nu_cont * 1e3 * wyler
M_eff_1d_meV = M_nu_1d * 1e3 * wyler
N_eff = 25 * (1 + 1/(2*np.pi**2))
dm31_cont = (1 - 1/N_eff) * (M_eff_cont_meV * 1e-3)**2
dm31_1d = (1 - 1/N_eff) * (M_eff_1d_meV * 1e-3)**2
dm21_cont = (d/(4*N_eff)) * (M_eff_cont_meV * 1e-3)**2
dm21_1d = (d/(4*N_eff)) * (M_eff_1d_meV * 1e-3)**2
obs_dm31 = 2.534e-3
obs_dm21 = 7.53e-5
print(f"  Neutrino M_eff (continuous): {M_eff_cont_meV:.4f} meV")
print(f"  Neutrino M_eff (1D disc):    {M_eff_1d_meV:.4f} meV")
print(f"  Dm31 (continuous): {dm31_cont:.4e} eV^2  (error: {(dm31_cont-obs_dm31)/obs_dm31*100:+.2f}%)")
print(f"  Dm31 (1D disc):    {dm31_1d:.4e} eV^2  (error: {(dm31_1d-obs_dm31)/obs_dm31*100:+.2f}%)")
print(f"  Dm21 (continuous): {dm21_cont:.4e} eV^2  (error: {(dm21_cont-obs_dm21)/obs_dm21*100:+.2f}%)")
print(f"  Dm21 (1D disc):    {dm21_1d:.4e} eV^2  (error: {(dm21_1d-obs_dm21)/obs_dm21*100:+.2f}%)")
print()
print("Done.")
