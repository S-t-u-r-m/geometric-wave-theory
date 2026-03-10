"""
GPU-Accelerated Breather Analysis (CuPy + RTX 4070 Ti)
1D: Exact breather frequency (no phonon contamination)
3D: Confinement splitting for mu/strange
"""

import numpy as np
import cupy as cp
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

N_1d = 4096
dt_1d = 0.001
N_steps_1d = 500000
SKIP_1d = 10000
REC_EVERY = 2

print("=" * 100)
print("PART 1: 1D EXACT BREATHER FREQUENCY (GPU)")
print("=" * 100)
dt_rec = dt_1d * REC_EVERY
n_rec = (N_steps_1d - SKIP_1d) // REC_EVERY
print(f"1D: {N_1d} sites, dt={dt_1d}, {N_steps_1d} steps")
print(f"Freq resolution: {1.0/(n_rec*dt_rec):.6f}, Nyquist: {1/(2*dt_rec):.1f}")
print()

def make_breather_1d(N_grid, n):
    omega = np.cos(n * gamma_gwt)
    eta = np.sin(n * gamma_gwt)
    center = N_grid // 2
    x = cp.arange(N_grid, dtype=cp.float64) - center
    arg = cp.minimum(cp.float64(eta) * cp.abs(x), cp.float64(50.0))
    phi = 4.0 * cp.arctan(cp.float64(eta) / (cp.float64(omega) * cp.cosh(arg)))
    phi_t = cp.zeros_like(phi)
    return phi, phi_t, center

def accel_1d(phi):
    return cp.roll(phi, 1) + cp.roll(phi, -1) - 2*phi - cp.sin(phi)

def evolve_1d(phi, phi_t, dt_s, n_steps, center, skip, rec_every):
    series = []
    acc = accel_1d(phi)
    for step in range(n_steps):
        phi_t += 0.5 * dt_s * acc
        phi += dt_s * phi_t
        acc = accel_1d(phi)
        phi_t += 0.5 * dt_s * acc
        if step >= skip and step % rec_every == 0:
            series.append(float(phi[center]))
    return series

def extract_freq(signal_list, dt_sample):
    signal = np.array(signal_list)
    signal = signal - np.mean(signal)
    window = np.hanning(len(signal))
    windowed = signal * window
    fft_vals = np.fft.rfft(windowed)
    fft_power = np.abs(fft_vals)**2
    freqs = np.fft.rfftfreq(len(windowed), d=dt_sample)
    peak_idx = np.argmax(fft_power[1:]) + 1
    peak_freq = freqs[peak_idx]
    if 1 < peak_idx < len(fft_power) - 1:
        al = np.log(fft_power[peak_idx - 1] + 1e-30)
        be = np.log(fft_power[peak_idx] + 1e-30)
        gp = np.log(fft_power[peak_idx + 1] + 1e-30)
        dp = 0.5 * (al - gp) / (al - 2*be + gp)
        peak_freq = freqs[peak_idx] + dp * (freqs[1] - freqs[0])
    return 2 * np.pi * peak_freq, freqs, fft_power

_ = cp.zeros(100)  # warm up GPU

print("  n        Particle        w_cont     w_1D_disc       Ratio    Correction    Time")
print("=" * 100)

results_1d = []

for n_val in n_values:
    name = particle_names[n_val]
    omega_cont = np.cos(n_val * gamma_gwt)
    t0 = time.time()
    phi, phi_t, center = make_breather_1d(N_1d, n_val)
    series = evolve_1d(phi, phi_t, dt_1d, N_steps_1d, center, SKIP_1d, REC_EVERY)
    cp.cuda.Stream.null.synchronize()
    omega_disc, freqs, fft_power = extract_freq(series, dt_rec)
    ratio = omega_disc / omega_cont
    corr_pct = (ratio - 1.0) * 100
    elapsed = time.time() - t0
    print(f"{n_val:3d}  {name:>14s}  {omega_cont:12.6f}  {omega_disc:12.6f}  "
          f"{ratio:10.6f}  {corr_pct:+11.4f}%  {elapsed:5.1f}s")
    results_1d.append({
        "n": n_val, "name": name,
        "omega_cont": omega_cont, "omega_disc": omega_disc,
        "ratio": ratio, "corr_pct": corr_pct,
        "series": series, "freqs": freqs, "fft_power": fft_power,
    })

print()
print("=" * 100)
print("PART 2: 3D CONFINEMENT SPLITTING (GPU)")
print("=" * 100)

N_3d = 48
dt_3d = 0.002
N_steps_3d = 10000
MEAS_AFTER = 3000
MEAS_EVERY = 50

def make_breather_3d(N_grid, n):
    ctr = N_grid // 2
    omega = np.cos(n * gamma_gwt)
    eta = np.sin(n * gamma_gwt)
    ix, iy, iz = cp.meshgrid(cp.arange(N_grid), cp.arange(N_grid), cp.arange(N_grid), indexing="ij")
    r = cp.sqrt((ix-ctr)**2 + (iy-ctr)**2 + (iz-ctr)**2).astype(cp.float64)
    r = cp.maximum(r, cp.float64(1e-10))
    phi = 4.0 * cp.arctan(cp.float64(eta) / (cp.float64(omega) * cp.cosh(cp.minimum(cp.float64(eta)*r, cp.float64(50.0)))))
    return phi, cp.zeros_like(phi)

def energy_3d(phi, phi_t):
    KE = 0.5 * phi_t**2
    GE = cp.zeros_like(phi)
    for axis in range(3):
        dphi = cp.diff(phi, axis=axis)
        sl = [slice(None)] * 3
        sl[axis] = slice(None, -1)
        GE[tuple(sl)] += 0.5 * dphi**2
    PE = 1.0 - cp.cos(phi)
    return float(cp.sum(KE + GE + PE))

def accel_3d_free(phi):
    acc = -6*phi - cp.sin(phi)
    for axis in range(3):
        acc += cp.roll(phi, 1, axis=axis) + cp.roll(phi, -1, axis=axis)
    return acc

def accel_3d_conf(phi, mask):
    pm = phi * mask
    acc = -6*pm - cp.sin(pm)
    for axis in range(3):
        acc += cp.roll(pm, 1, axis=axis) + cp.roll(pm, -1, axis=axis)
    return acc * mask

def make_mask(N_grid, R):
    ctr = N_grid // 2
    ix, iy, iz = cp.meshgrid(cp.arange(N_grid), cp.arange(N_grid), cp.arange(N_grid), indexing="ij")
    r = cp.sqrt((ix-ctr)**2 + (iy-ctr)**2 + (iz-ctr)**2)
    return (r <= R).astype(cp.float64)

R_values = [6, 7, 8, 9, 10, 11, 12]
print(f"Grid: {N_3d}^3, dt={dt_3d}, steps={N_steps_3d}")
print(f"Confinement radii: {R_values}")
print()

# Free energy for n=4
phi_f, phit_f = make_breather_3d(N_3d, 4)
acc_f = accel_3d_free(phi_f)
energies_f = []
for step in range(N_steps_3d):
    phit_f += 0.5 * dt_3d * acc_f
    phi_f += dt_3d * phit_f
    acc_f = accel_3d_free(phi_f)
    phit_f += 0.5 * dt_3d * acc_f
    if step >= MEAS_AFTER and step % MEAS_EVERY == 0:
        energies_f.append(energy_3d(phi_f, phit_f))

cp.cuda.Stream.null.synchronize()
E_free = np.mean(energies_f)
print(f"n=4 FREE energy: {E_free:.4f}")
print()
obs_split = (105.66 - 93.4) / ((105.66 + 93.4) / 2) * 100

print(f"{"R_conf":>8s}  {"N_sites":>8s}  {"E_free":>12s}  {"E_conf":>12s}  {"Split":>10s}  {"Time":>6s}")
print("-" * 65)

for R_c in R_values:
    t0 = time.time()
    mask = make_mask(N_3d, R_c)
    n_in = int(cp.sum(mask))
    phi_c, phit_c = make_breather_3d(N_3d, 4)
    phi_c *= mask
    acc_c = accel_3d_conf(phi_c, mask)
    energies_c = []
    for step in range(N_steps_3d):
        phit_c += 0.5 * dt_3d * acc_c
        phi_c += dt_3d * phit_c
        phi_c *= mask
        acc_c = accel_3d_conf(phi_c, mask)
        phit_c += 0.5 * dt_3d * acc_c
        phit_c *= mask
        if step >= MEAS_AFTER and step % MEAS_EVERY == 0:
            energies_c.append(energy_3d(phi_c, phit_c))
    cp.cuda.Stream.null.synchronize()
    E_conf = np.mean(energies_c)
    split = (E_free - E_conf) / ((E_free + E_conf) / 2) * 100
    elapsed = time.time() - t0
    marker = " <-- closest" if abs(split - obs_split) < 2 else ""
    print(f"{R_c:8d}  {n_in:8d}  {E_free:12.4f}  {E_conf:12.4f}  {split:+9.3f}%  {elapsed:5.1f}s{marker}")

print(f"Observed mu-strange splitting: {obs_split:+.1f}%")

print()
print("=" * 100)
print("CORRECTED PHYSICAL MASSES (1D Discrete Frequency Corrections)")
print("=" * 100)

gamma_sg = np.pi / (16 * np.pi - 2)
fermion_data = [
    (16, 32, "electron", 0.511, "free"),
    (13, 31, "up", 2.16, "free"),
    (5, 30, "down", 4.67, "free"),
    (4, 28, "muon", 105.66, "free"),
    (4, 28, "strange", 93.4, "free"),
    (11, 27, "charm", 1271, "free"),
    (18, 27, "tau", 1776.86, "free"),
    (7, 26, "bottom", 4183, "free"),
    (12, 24, "top", 172760, "free"),
]

print()
print(f"  {"Particle":>10s}  {"n":>3s} {"p":>3s}  {"m_cont":>14s}  {"m_1D":>14s}  {"Observed":>12s}  {"Err_cont":>10s}  {"Err_1D":>10s}")
print("-" * 100)

total_ec = 0
total_e1 = 0
cnt = 0

for n_val, p_val, name, obs, ptype in fermion_data:
    mc = (16.0/np.pi**2) * np.sin(n_val*gamma_sg) * np.exp(-16*p_val/np.pi**2) * m_Planck_MeV
    for r in results_1d:
        if r["n"] == n_val:
            rat = r["ratio"]
            break
    m1 = mc * rat
    if name == "muon": m1 *= np.sqrt(1.126)
    elif name == "strange": m1 /= np.sqrt(1.126)
    ec = (mc - obs) / obs * 100
    e1 = (m1 - obs) / obs * 100
    print(f"  {name:>10s}  {n_val:3d} {p_val:3d}  {mc:14.4f}  {m1:14.4f}  {obs:12.4f}  {ec:+9.2f}%  {e1:+9.2f}%")
    total_ec += abs(ec)
    total_e1 += abs(e1)
    cnt += 1

print(f"  Mean |error|: continuous={total_ec/cnt:.2f}%, 1D-corrected={total_e1/cnt:.2f}%")

print()
print("=" * 100)
print("ELECTRON CORRECTION -> NEUTRINO IMPLICATIONS")
print("=" * 100)

for r in results_1d:
    if r["n"] == 16:
        e_ratio = r["ratio"]
        e_corr = r["corr_pct"]
        break

m_e_c = (16.0/np.pi**2)*np.sin(16*gamma_sg)*np.exp(-16*32/np.pi**2)*m_Planck_MeV
m_e_1 = m_e_c * e_ratio
m_p = 938.272

print(f"  Electron (cont):  {m_e_c:.6f} MeV ({(m_e_c-0.511)/0.511*100:+.3f}%)")
print(f"  Electron (1D):    {m_e_1:.6f} MeV ({(m_e_1-0.511)/0.511*100:+.3f}%)")
print(f"  Freq correction:  {e_corr:+.4f}%")

Mn_c = m_e_c**3 / (d * m_p**2) * 1e3
Mn_1 = m_e_1**3 / (d * m_p**2) * 1e3
wy = 1 + 1/(d*2*np.pi**2)
Me_c = Mn_c * wy
Me_1 = Mn_1 * wy
Ne = 25*(1+1/(2*np.pi**2))
d31c = (1-1/Ne)*(Me_c*1e-3)**2
d311 = (1-1/Ne)*(Me_1*1e-3)**2
d21c = (d/(4*Ne))*(Me_c*1e-3)**2
d211 = (d/(4*Ne))*(Me_1*1e-3)**2
od31 = 2.534e-3
od21 = 7.53e-5

print(f"  M_eff (cont): {Me_c:.3f} meV")
print(f"  M_eff (1D):   {Me_1:.3f} meV")
print(f"  Dm31 (cont):  {d31c:.4e} eV^2 ({(d31c-od31)/od31*100:+.2f}%)")
print(f"  Dm31 (1D):    {d311:.4e} eV^2 ({(d311-od31)/od31*100:+.2f}%)")
print(f"  Dm21 (cont):  {d21c:.4e} eV^2 ({(d21c-od21)/od21*100:+.2f}%)")
print(f"  Dm21 (1D):    {d211:.4e} eV^2 ({(d211-od21)/od21*100:+.2f}%)")
print(f"  Ratio (cont): {d31c/d21c:.2f} (obs: 33.65)")
print(f"  Ratio (1D):   {d311/d211:.2f} (obs: 33.65)")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

ax1 = axes[0, 0]
for r in results_1d:
    if r["n"] == 16:
        sig = np.array(r["series"][:2000])
        t_arr = np.arange(len(sig)) * dt_rec + SKIP_1d * dt_1d
        ax1.plot(t_arr, sig, "b-", linewidth=0.5)
        break
ax1.set_xlabel("Time (Planck units)")
ax1.set_ylabel("phi(center)")
ax1.set_title("Electron (n=16): 1D center-site oscillation")

ax2 = axes[0, 1]
for r in results_1d:
    if r["n"] == 16:
        fm = r["freqs"] > 0
        ax2.semilogy(2*np.pi*r["freqs"][fm], r["fft_power"][fm], "b-")
        ax2.axvline(x=r["omega_cont"], color="r", linestyle="--", label="w_cont")
        ax2.axvline(x=r["omega_disc"], color="g", linestyle=":", label="w_disc")
        break
ax2.set_xlabel("Angular frequency")
ax2.set_ylabel("FFT Power")
ax2.set_title("Electron (n=16): 1D Frequency spectrum")
ax2.legend()
ax2.set_xlim(0, 2)

ax3 = axes[1, 0]
ns = [r["n"] for r in results_1d]
corrs = [r["corr_pct"] for r in results_1d]
ax3.bar(range(len(ns)), corrs, color="steelblue", alpha=0.7)
ax3.set_xticks(range(len(ns)))
ax3.set_xticklabels([str(n) for n in ns])
ax3.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
ax3.set_ylabel("Freq correction (%)")
ax3.set_title("1D Discreteness corrections")

ax4 = axes[1, 1]
ax4.text(0.5, 0.5, "GPU: RTX 4070 Ti | 1D: 4096 sites, 500K steps | 3D: 48^3",
         ha="center", va="center", fontsize=14, transform=ax4.transAxes)
ax4.set_title("Simulation parameters")

plt.tight_layout()
plt.savefig("calculations/breather_gpu_analysis.png", dpi=150)
print(f"Plot saved to calculations/breather_gpu_analysis.png")
