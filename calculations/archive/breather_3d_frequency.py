"""
3D Breather Frequency Analysis
================================
Measures the oscillation frequency of each breather mode on the 3D cubic lattice
by tracking phi at the center site and computing the FFT.

The ratio omega_discrete / omega_continuous gives the mass correction:
  m_corrected = m_continuous * (omega_discrete / omega_continuous)
"""

import numpy as np
import time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = 3
gamma_gwt = np.pi / (16 * np.pi - 2)
m_Planck_MeV = 1.2209e22

n_values = [4, 5, 7, 11, 12, 13, 16, 18]
particle_names = {
    4: "muon/strange", 5: "down", 7: "bottom", 11: "charm",
    12: "top", 13: "up", 16: "electron", 18: "tau"
}

N = 48
a = 1.0
dt = 0.002
N_steps = 50000
RECORD_EVERY = 1
SKIP_TRANSIENT = 2000

print(f"3D Lattice: {N}^3 = {N**3} sites")
print(f"Time step: {dt}, total time: {N_steps * dt}")
print(f"Frequency resolution: df = {1.0 / ((N_steps - SKIP_TRANSIENT) * dt):.4f}")
print(f"Nyquist frequency: {1.0 / (2 * dt):.1f}")
print()


def make_breather_3d(N_grid, n, a_spacing):
    center = np.array([N_grid // 2, N_grid // 2, N_grid // 2])
    omega = np.cos(n * gamma_gwt)
    eta = np.sin(n * gamma_gwt)
    ix, iy, iz = np.meshgrid(
        np.arange(N_grid), np.arange(N_grid), np.arange(N_grid), indexing='ij')
    r = a_spacing * np.sqrt(
        (ix - center[0])**2 + (iy - center[1])**2 + (iz - center[2])**2)
    r = np.maximum(r, 1e-10)
    phi = 4.0 * np.arctan(eta / (omega * np.cosh(np.minimum(eta * r, 50))))
    phi_t = np.zeros_like(phi)
    return phi, phi_t, center


def acceleration_3d_free(phi, a_spacing):
    acc = -6 * phi / a_spacing**2 - np.sin(phi)
    for axis in range(3):
        acc += (np.roll(phi, 1, axis=axis) + np.roll(phi, -1, axis=axis)) / a_spacing**2
    return acc


def evolve_and_record(phi, phi_t, dt_step, a_spacing, n_steps, center,
                      skip_transient=2000, record_every=1):
    center_series = []
    times = []
    acc = acceleration_3d_free(phi, a_spacing)
    cx, cy, cz = center
    for step in range(n_steps):
        phi_t += 0.5 * dt_step * acc
        phi += dt_step * phi_t
        acc = acceleration_3d_free(phi, a_spacing)
        phi_t += 0.5 * dt_step * acc
        if step >= skip_transient and step % record_every == 0:
            center_series.append(phi[cx, cy, cz])
            times.append(step * dt_step)
    return np.array(center_series), np.array(times)


def extract_frequency(signal, times):
    signal = signal - np.mean(signal)
    window = np.hanning(len(signal))
    windowed = signal * window
    fft_vals = np.fft.rfft(windowed)
    fft_power = np.abs(fft_vals)**2
    dt_sample = times[1] - times[0]
    freqs = np.fft.rfftfreq(len(windowed), d=dt_sample)
    peak_idx = np.argmax(fft_power[1:]) + 1
    peak_freq = freqs[peak_idx]
    if 1 < peak_idx < len(fft_power) - 1:
        al = np.log(fft_power[peak_idx - 1] + 1e-30)
        be = np.log(fft_power[peak_idx] + 1e-30)
        gp = np.log(fft_power[peak_idx + 1] + 1e-30)
        delta_p = 0.5 * (al - gp) / (al - 2*be + gp)
        peak_freq = freqs[peak_idx] + delta_p * (freqs[1] - freqs[0])
    omega_measured = 2 * np.pi * peak_freq
    return omega_measured, freqs, fft_power


print('=' * 100)
print('3D BREATHER FREQUENCY ANALYSIS')
print('=' * 100)

hdr = f"{'n':>3s}  {'Particle':>14s}  {'omega_cont':>12s}  {'omega_3D':>12s}  {'Ratio':>10s}  {'Correction':>12s}  {'Time':>6s}"
print()
print(hdr)
print('=' * 100)

results = []

for n_val in n_values:
    name = particle_names[n_val]
    omega_continuous = np.cos(n_val * gamma_gwt)
    t0 = time.time()
    phi, phi_t, center = make_breather_3d(N, n_val, a)
    center_series, times = evolve_and_record(
        phi, phi_t, dt, a, N_steps, center,
        skip_transient=SKIP_TRANSIENT, record_every=RECORD_EVERY)
    omega_3d, freqs, fft_power = extract_frequency(center_series, times)
    ratio = omega_3d / omega_continuous
    correction_pct = (ratio - 1.0) * 100
    elapsed = time.time() - t0
    print(f"{n_val:3d}  {name:>14s}  {omega_continuous:12.6f}  {omega_3d:12.6f}  "
          f"{ratio:10.6f}  {correction_pct:+11.4f}%  {elapsed:5.1f}s")
    results.append({
        'n': n_val, 'name': name,
        'omega_continuous': omega_continuous, 'omega_3d': omega_3d,
        'ratio': ratio, 'correction_pct': correction_pct,
        'center_series': center_series, 'times': times,
        'freqs': freqs, 'fft_power': fft_power,
    })


# CORRECTED PHYSICAL MASSES
print()
print('=' * 100)
print('CORRECTED PHYSICAL MASSES WITH 3D FREQUENCY SHIFTS')
print('=' * 100)

gamma_sg = np.pi / (16 * np.pi - 2)
fermion_data = [
    (16, 32, 'electron', 0.511,   'free'),
    (13, 31, 'up',       2.16,    'free'),
    (5,  30, 'down',     4.67,    'free'),
    (4,  28, 'muon',     105.66,  'free'),
    (4,  28, 'strange',  93.4,    'free'),
    (11, 27, 'charm',    1271,    'free'),
    (18, 27, 'tau',      1776.86, 'free'),
    (7,  26, 'bottom',   4183,    'free'),
    (12, 24, 'top',      172760,  'free'),
]

mhdr = f"{'Particle':>10s}  {'n':>3s} {'p':>3s}  {'m_cont(MeV)':>14s}  {'m_3D(MeV)':>14s}  {'Observed':>12s}  {'Err_cont':>10s}  {'Err_3D':>10s}"
print()
print(mhdr)
print('-' * 100)

total_err_cont = 0
total_err_3d = 0
count = 0

for n_val, p_val, name, obs_MeV, ptype in fermion_data:
    m_cont = (16.0/np.pi**2) * np.sin(n_val * gamma_sg) * np.exp(-16*p_val/np.pi**2) * m_Planck_MeV
    for r in results:
        if r['n'] == n_val:
            ratio = r['ratio']
            break
    m_3d = m_cont * ratio
    if name == 'muon':
        m_3d = m_3d * np.sqrt(1.126)
    elif name == 'strange':
        m_3d = m_3d / np.sqrt(1.126)
    err_cont = (m_cont - obs_MeV) / obs_MeV * 100
    err_3d = (m_3d - obs_MeV) / obs_MeV * 100
    print(f"{name:>10s}  {n_val:3d} {p_val:3d}  {m_cont:14.4f}  "
          f"{m_3d:14.4f}  {obs_MeV:12.4f}  {err_cont:+9.2f}%  {err_3d:+9.2f}%")
    total_err_cont += abs(err_cont)
    total_err_3d += abs(err_3d)
    count += 1

print(f"  Mean |error|:  continuous = {total_err_cont/count:.2f}%,  3D-corrected = {total_err_3d/count:.2f}%")


# ELECTRON -> NEUTRINO IMPLICATIONS
print()
print('=' * 100)
print('ELECTRON CORRECTION -> NEUTRINO IMPLICATIONS')
print('=' * 100)

for r in results:
    if r['n'] == 16:
        e_ratio = r['ratio']
        e_corr = r['correction_pct']
        break

m_e_cont = (16.0/np.pi**2) * np.sin(16 * gamma_sg) * np.exp(-16*32/np.pi**2) * m_Planck_MeV
m_e_3d = m_e_cont * e_ratio
m_p_MeV = 938.272

print(f"  Electron mass (continuous): {m_e_cont:.4f} MeV (error: {(m_e_cont-0.511)/0.511*100:+.2f}%)")
print(f"  Electron mass (3D lattice): {m_e_3d:.4f} MeV (error: {(m_e_3d-0.511)/0.511*100:+.2f}%)")
print(f"  Frequency correction: {e_corr:+.4f}%")

M_nu_cont = m_e_cont**3 / (d * m_p_MeV**2)
M_nu_3d = m_e_3d**3 / (d * m_p_MeV**2)
M_nu_cont_meV = M_nu_cont * 1e3
M_nu_3d_meV = M_nu_3d * 1e3

wyler = 1 + 1/(d * 2 * np.pi**2)
M_eff_cont = M_nu_cont_meV * wyler
M_eff_3d = M_nu_3d_meV * wyler

N_eff = 25 * (1 + 1/(2*np.pi**2))
dm31_cont = (1 - 1/N_eff) * (M_eff_cont * 1e-3)**2
dm31_3d = (1 - 1/N_eff) * (M_eff_3d * 1e-3)**2
dm21_cont = (d/(4*N_eff)) * (M_eff_cont * 1e-3)**2
dm21_3d = (d/(4*N_eff)) * (M_eff_3d * 1e-3)**2

obs_dm31 = 2.534e-3
obs_dm21 = 7.53e-5

print(f"  Neutrino scale (continuous): M_eff = {M_eff_cont:.2f} meV")
print(f"  Neutrino scale (3D):         M_eff = {M_eff_3d:.2f} meV")
print(f"  Dm31 (continuous): {dm31_cont:.4e} eV^2  (error: {(dm31_cont-obs_dm31)/obs_dm31*100:+.2f}%)")
print(f"  Dm31 (3D):         {dm31_3d:.4e} eV^2  (error: {(dm31_3d-obs_dm31)/obs_dm31*100:+.2f}%)")
print(f"  Dm21 (continuous): {dm21_cont:.4e} eV^2  (error: {(dm21_cont-obs_dm21)/obs_dm21*100:+.2f}%)")
print(f"  Dm21 (3D):         {dm21_3d:.4e} eV^2  (error: {(dm21_3d-obs_dm21)/obs_dm21*100:+.2f}%)")


# PLOTS
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

ax1 = axes[0, 0]
for r in results:
    if r['n'] == 16:
        ax1.plot(r['times'][:500], r['center_series'][:500], 'b-', linewidth=0.5)
        break
ax1.set_xlabel('Time (Planck units)')
ax1.set_ylabel('phi(center)')
ax1.set_title('Electron (n=16): Center-site oscillation')

ax2 = axes[0, 1]
for r in results:
    if r['n'] == 16:
        fmask = r['freqs'] > 0
        ax2.semilogy(2*np.pi*r['freqs'][fmask], r['fft_power'][fmask], 'b-')
        ax2.axvline(x=r['omega_continuous'], color='r', linestyle='--',
                    label=f"omega_cont = {r['omega_continuous']:.4f}")
        ax2.axvline(x=r['omega_3d'], color='g', linestyle=':',
                    label=f"omega_3D = {r['omega_3d']:.4f}")
        break
ax2.set_xlabel('Angular frequency omega')
ax2.set_ylabel('FFT Power')
ax2.set_title('Electron (n=16): Frequency spectrum')
ax2.legend()
ax2.set_xlim(0, 2)

ax3 = axes[1, 0]
ns = [r['n'] for r in results]
corrs = [r['correction_pct'] for r in results]
colors = ['blue' if n in [4, 16, 18] else 'red' for n in ns]
ax3.bar(range(len(ns)), corrs, color=colors, alpha=0.7)
ax3.set_xticks(range(len(ns)))
ax3.set_xticklabels([f"n={n}" for n in ns], fontsize=7)
ax3.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax3.set_ylabel('Frequency correction (%)')
ax3.set_title('3D Discreteness: omega_3D / omega_cont - 1')

ax4 = axes[1, 1]
names_plot = []
errs_cont_plot = []
errs_3d_plot = []
for n_val, p_val, name, obs_MeV, ptype in fermion_data:
    m_cont = (16.0/np.pi**2) * np.sin(n_val * gamma_sg) * np.exp(-16*p_val/np.pi**2) * m_Planck_MeV
    for r in results:
        if r['n'] == n_val:
            ratio = r['ratio']
            break
    m_3d = m_cont * ratio
    if name == 'muon':
        m_3d *= np.sqrt(1.126)
    elif name == 'strange':
        m_3d /= np.sqrt(1.126)
    names_plot.append(name)
    errs_cont_plot.append((m_cont - obs_MeV) / obs_MeV * 100)
    errs_3d_plot.append((m_3d - obs_MeV) / obs_MeV * 100)

x = np.arange(len(names_plot))
ax4.bar(x - 0.15, [abs(e) for e in errs_cont_plot], 0.3, label='Continuous', alpha=0.7, color='gray')
ax4.bar(x + 0.15, [abs(e) for e in errs_3d_plot], 0.3, label='3D corrected', alpha=0.7, color='green')
ax4.set_xticks(x)
ax4.set_xticklabels(names_plot, rotation=45, fontsize=7)
ax4.set_ylabel('|Error| (%)')
ax4.set_title('Mass predictions: Continuous vs 3D-corrected')
ax4.legend()

plt.tight_layout()
plt.savefig('calculations/breather_3d_frequency.png', dpi=150)
print(f"Plot saved to calculations/breather_3d_frequency.png")
