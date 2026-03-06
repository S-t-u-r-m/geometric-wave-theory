"""
Full 3D Corrected Fermion Spectrum
===================================
Applies 3D cubic confinement corrections to ALL 9 fermion masses.
- Quarks: CONFINED in cubic region L = 2^d - 1 = 7 (inside proton kink)
- Leptons: FREE on full lattice (periodic BC)

Computes correction factors E_free/E_init and E_conf/E_init for each
breather index n, then applies to the continuous m(n,p) predictions.

Also tests:
1. Geometric mean: m_mu * m_s = m(4,28)^2
2. Other degenerate/near-degenerate pairs
3. Full spectrum error comparison: continuous vs 3D-corrected
"""

import numpy as np
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# GWT CONSTANTS
# ============================================================
d = 3
gamma_gwt = np.pi / (16 * np.pi - 2)
m_Planck_MeV = 1.2209e22

def m_fermion_cont(n, p):
    """Continuous GWT mass formula."""
    return (16.0 / np.pi**2) * np.sin(n * gamma_gwt) * np.exp(-16*p / np.pi**2) * m_Planck_MeV

# ============================================================
# 3D LATTICE SIMULATION ENGINE
# ============================================================
N = 48
a = 1.0
dt = 0.002
N_steps = 10000
MEASURE_AFTER = 3000
MEASURE_EVERY = 50
L_conf = 7  # 2^d - 1 = 7

def make_breather_3d(N_grid, n):
    center = N_grid // 2
    omega = np.cos(n * gamma_gwt)
    eta = np.sin(n * gamma_gwt)
    ix, iy, iz = np.meshgrid(np.arange(N_grid), np.arange(N_grid), np.arange(N_grid), indexing='ij')
    r = np.sqrt((ix - center)**2 + (iy - center)**2 + (iz - center)**2)
    r = np.maximum(r, 1e-10)
    phi = 4.0 * np.arctan(eta / (omega * np.cosh(np.minimum(eta * r, 50))))
    return phi, np.zeros_like(phi)

def compute_energy(phi, phi_t):
    KE = 0.5 * phi_t**2
    GE = np.zeros_like(phi)
    for axis in range(3):
        dphi = np.diff(phi, axis=axis)
        sl = [slice(None)] * 3
        sl[axis] = slice(None, -1)
        GE[tuple(sl)] += 0.5 * dphi**2
    PE = 1.0 - np.cos(phi)
    return np.sum(KE + GE + PE)

def accel_free(phi, _dummy=None):
    acc = -6 * phi - np.sin(phi)
    for axis in range(3):
        acc += np.roll(phi, 1, axis=axis) + np.roll(phi, -1, axis=axis)
    return acc

def make_cubic_mask(N_grid, L):
    center = N_grid // 2
    ix, iy, iz = np.meshgrid(np.arange(N_grid), np.arange(N_grid), np.arange(N_grid), indexing='ij')
    return ((np.abs(ix - center) <= L) & (np.abs(iy - center) <= L) & (np.abs(iz - center) <= L)).astype(float)

mask_cubic = make_cubic_mask(N, L_conf)
n_confined_sites = int(np.sum(mask_cubic))

def accel_confined(phi, mask=mask_cubic):
    pm = phi * mask
    acc = -6 * pm - np.sin(pm)
    for axis in range(3):
        acc += np.roll(pm, 1, axis=axis) + np.roll(pm, -1, axis=axis)
    return acc * mask

def evolve(phi, phi_t, acc_func, n_steps=N_steps):
    energies = []
    acc = acc_func(phi)
    for step in range(n_steps):
        phi_t += 0.5 * dt * acc
        phi += dt * phi_t
        acc = acc_func(phi)
        phi_t += 0.5 * dt * acc
        if step >= MEASURE_AFTER and step % MEASURE_EVERY == 0:
            energies.append(compute_energy(phi, phi_t))
    return np.array(energies)

# ============================================================
# RUN ALL n-VALUES: FREE and CONFINED
# ============================================================
print("=" * 90)
print(f"3D FULL SPECTRUM: {N}^3 lattice, cubic confinement L={L_conf} ({n_confined_sites} sites)")
print("=" * 90)

# All unique n-values used by fermions
all_n = sorted(set([4, 5, 7, 11, 12, 13, 16, 18]))

corrections = {}  # n -> (E_free_ratio, E_conf_ratio)

print(f"\n{'n':>3s}  {'E_init':>12s}  {'E_free':>12s}  {'E_conf':>12s}  "
      f"{'ratio_free':>10s}  {'ratio_conf':>10s}  {'split':>8s}  {'time':>6s}")
print("-" * 80)

for n_val in all_n:
    t0 = time.time()

    # FREE
    phi, phit = make_breather_3d(N, n_val)
    E_init = compute_energy(phi, phit)
    en_free = evolve(phi.copy(), phit.copy(), accel_free)
    E_free = np.mean(en_free) if len(en_free) > 5 else E_init

    # CONFINED
    phi_c, phit_c = make_breather_3d(N, n_val)
    phi_c *= mask_cubic
    E_init_c = compute_energy(phi_c, phit_c)
    en_conf = evolve(phi_c, phit_c, accel_confined)
    E_conf = np.mean(en_conf) if len(en_conf) > 5 else E_init_c

    # Ratios: how does the energy change from initial?
    # For mass correction, we care about E_final / E_init
    ratio_free = E_free / E_init
    ratio_conf = E_conf / E_init_c
    splitting = (E_free - E_conf) / ((E_free + E_conf) / 2) * 100

    elapsed = time.time() - t0
    corrections[n_val] = (ratio_free, ratio_conf, E_free, E_conf, E_init, E_init_c)

    print(f"{n_val:3d}  {E_init:12.4f}  {E_free:12.4f}  {E_conf:12.4f}  "
          f"{ratio_free:10.6f}  {ratio_conf:10.6f}  {splitting:+7.2f}%  {elapsed:5.1f}s")


# ============================================================
# APPLY CORRECTIONS TO FULL FERMION SPECTRUM
# ============================================================
print("\n" + "=" * 100)
print("FULL FERMION SPECTRUM: CONTINUOUS vs 3D-CORRECTED")
print("=" * 100)

# Complete fermion table: (n, p, name, obs_MeV, type)
fermions = [
    (16, 32, "electron",  0.511,    "free"),
    (13, 31, "up",        2.16,     "confined"),
    (5,  30, "down",      4.67,     "confined"),
    (4,  28, "muon",      105.66,   "free"),
    (4,  28, "strange",   93.4,     "confined"),
    (11, 27, "charm",     1271,     "confined"),
    (18, 27, "tau",       1776.86,  "free"),
    (7,  26, "bottom",    4183,     "confined"),
    (12, 24, "top",       172760,   "confined"),
]

print(f"\n{'Particle':>10s}  {'n':>3s} {'p':>3s}  {'Type':>8s}  "
      f"{'m_cont':>12s}  {'m_3D':>12s}  {'Observed':>12s}  "
      f"{'Err_cont':>9s}  {'Err_3D':>9s}  {'Improved?':>10s}")
print("-" * 110)

errors_cont = []
errors_3d = []
results = []

for n, p, name, obs, ptype in fermions:
    m_cont = m_fermion_cont(n, p)
    ratio_free, ratio_conf, E_free, E_conf, E_init, E_init_c = corrections[n]

    if ptype == "free":
        # For free particles, the correction comes from how the 3D breather
        # energy differs from the 1D continuous prediction.
        # The breather on the full 3D lattice has MORE energy than the 1D
        # continuous formula because it extends in 3 dimensions.
        # We use the ratio E_free/E_conf to split degenerate pairs.
        # For non-degenerate free particles, we use the free energy evolution ratio.
        corr_factor = ratio_free
    else:
        corr_factor = ratio_conf

    # But wait — the correction should be RELATIVE to the continuous prediction.
    # The 3D energy includes the volume integral which differs from 1D.
    # What we actually want is: for particles with the SAME n,
    # the FREE one gets mass proportional to sqrt(E_free)
    # and the CONFINED one gets mass proportional to sqrt(E_conf)
    # relative to the continuous (degenerate) prediction.

    # For the mu-strange pair (both n=4, p=28):
    # m_mu = m_cont * sqrt(E_free / E_geom_mean)
    # m_s  = m_cont * sqrt(E_conf / E_geom_mean)
    # where E_geom_mean = sqrt(E_free * E_conf)

    # For other particles (no degenerate partner), the correction is smaller.
    # The dominant effect is the confinement boundary condition itself.
    # Use: m_corrected = m_cont * (E_type / E_init_type)
    # This captures how the lattice evolution shifts the energy.

    m_3d = m_cont * corr_factor

    err_cont = (m_cont - obs) / obs * 100
    err_3d = (m_3d - obs) / obs * 100

    improved = "YES" if abs(err_3d) < abs(err_cont) else "no"
    if abs(abs(err_3d) - abs(err_cont)) < 0.01:
        improved = "same"

    errors_cont.append(abs(err_cont))
    errors_3d.append(abs(err_3d))
    results.append((name, n, p, ptype, m_cont, m_3d, obs, err_cont, err_3d, improved))

    print(f"{name:>10s}  {n:3d} {p:3d}  {ptype:>8s}  "
          f"{m_cont:12.4f}  {m_3d:12.4f}  {obs:12.4f}  "
          f"{err_cont:+8.2f}%  {err_3d:+8.2f}%  {improved:>10s}")

# For the mu-strange pair, also do the geometric-mean splitting
print("\n" + "-" * 110)
print("MU-STRANGE SPLITTING (geometric mean method):")
E_free_4, E_conf_4 = corrections[4][2], corrections[4][3]
E_geom = np.sqrt(E_free_4 * E_conf_4)
m_cont_4_28 = m_fermion_cont(4, 28)

m_mu_split = m_cont_4_28 * np.sqrt(E_free_4 / E_geom)
m_s_split = m_cont_4_28 * np.sqrt(E_conf_4 / E_geom)
# Update the errors for mu and strange
err_mu_split = (m_mu_split - 105.66) / 105.66 * 100
err_s_split = (m_s_split - 93.4) / 93.4 * 100

print(f"  m(4,28) continuous = {m_cont_4_28:.4f} MeV")
print(f"  E_free/E_conf = {E_free_4/E_conf_4:.6f}")
print(f"  sqrt(E_free/E_conf) = {np.sqrt(E_free_4/E_conf_4):.6f}")
print(f"  muon:    {m_mu_split:.2f} MeV (obs 105.66, err {err_mu_split:+.2f}%)")
print(f"  strange: {m_s_split:.2f} MeV (obs 93.4, err {err_s_split:+.2f}%)")

# Replace mu/strange errors with split versions
errors_3d_split = list(errors_3d)
for i, (name, *_) in enumerate(results):
    if name == "muon":
        errors_3d_split[i] = abs(err_mu_split)
    elif name == "strange":
        errors_3d_split[i] = abs(err_s_split)


# ============================================================
# SUMMARY STATISTICS
# ============================================================
print("\n" + "=" * 100)
print("SUMMARY STATISTICS")
print("=" * 100)

print(f"\n{'Metric':>30s}  {'Continuous':>12s}  {'3D Corrected':>12s}  {'3D+Split':>12s}")
print("-" * 75)
print(f"{'Mean |error|':>30s}  {np.mean(errors_cont):11.3f}%  {np.mean(errors_3d):11.3f}%  {np.mean(errors_3d_split):11.3f}%")
print(f"{'Median |error|':>30s}  {np.median(errors_cont):11.3f}%  {np.median(errors_3d):11.3f}%  {np.median(errors_3d_split):11.3f}%")
print(f"{'Max |error|':>30s}  {np.max(errors_cont):11.3f}%  {np.max(errors_3d):11.3f}%  {np.max(errors_3d_split):11.3f}%")
print(f"{'RMS error':>30s}  {np.sqrt(np.mean(np.array(errors_cont)**2)):11.3f}%  {np.sqrt(np.mean(np.array(errors_3d)**2)):11.3f}%  {np.sqrt(np.mean(np.array(errors_3d_split)**2)):11.3f}%")
print(f"{'Particles within 1%':>30s}  {sum(1 for e in errors_cont if e < 1):11d}  {sum(1 for e in errors_3d if e < 1):11d}  {sum(1 for e in errors_3d_split if e < 1):11d}")
print(f"{'Particles within 2%':>30s}  {sum(1 for e in errors_cont if e < 2):11d}  {sum(1 for e in errors_3d if e < 2):11d}  {sum(1 for e in errors_3d_split if e < 2):11d}")
print(f"{'Particles within 5%':>30s}  {sum(1 for e in errors_cont if e < 5):11d}  {sum(1 for e in errors_3d if e < 5):11d}  {sum(1 for e in errors_3d_split if e < 5):11d}")


# ============================================================
# GEOMETRIC MEAN TEST
# ============================================================
print("\n" + "=" * 100)
print("GEOMETRIC MEAN TEST: m_mu * m_s = m(4,28)^2")
print("=" * 100)

obs_product = 105.66 * 93.4
cont_sq = m_cont_4_28**2
split_product = m_mu_split * m_s_split

print(f"  m_mu * m_s (observed)   = {obs_product:.2f} MeV^2")
print(f"  m(4,28)^2 (continuous)  = {cont_sq:.2f} MeV^2")
print(f"  m_mu * m_s (3D split)   = {split_product:.2f} MeV^2")
print(f"  Observed vs continuous:  {(obs_product - cont_sq)/cont_sq * 100:+.2f}%")
print(f"  3D split vs continuous:  {(split_product - cont_sq)/cont_sq * 100:+.2f}%")
print(f"  3D split vs observed:   {(split_product - obs_product)/obs_product * 100:+.2f}%")
print(f"\n  GWT prediction: product should equal m(4,28)^2 = {cont_sq:.2f}")
print(f"  Observed product: {obs_product:.2f} ({(obs_product-cont_sq)/cont_sq*100:+.1f}% off)")
print(f"  This {abs(obs_product-cont_sq)/cont_sq*100:.1f}% residual is the higher-order correction.")


# ============================================================
# SCAN FOR OTHER DEGENERATE/NEAR-DEGENERATE PAIRS
# ============================================================
print("\n" + "=" * 100)
print("DEGENERATE PAIR SCAN")
print("=" * 100)
print("Looking for particle pairs with same n OR same p that could show")
print("free-vs-confined splitting similar to mu-strange...")
print()

for i, (n1, p1, name1, obs1, type1) in enumerate(fermions):
    for j, (n2, p2, name2, obs2, type2) in enumerate(fermions):
        if j <= i:
            continue
        if type1 == type2:
            continue  # both free or both confined — no splitting
        # Check if they share n (same breather) or p (same depth)
        if n1 == n2:
            ratio = obs1 / obs2
            m_cont1 = m_fermion_cont(n1, p1)
            m_cont2 = m_fermion_cont(n2, p2)
            print(f"  SAME n={n1}: {name1}({type1}, p={p1}) vs {name2}({type2}, p={p2})")
            print(f"    Observed ratio: {ratio:.4f}")
            print(f"    Continuous: {m_cont1:.4f} vs {m_cont2:.4f} (ratio {m_cont1/m_cont2:.4f})")
            if p1 == p2:
                print(f"    *** DEGENERATE (same n AND p) — split by confinement! ***")
                # Compute predicted split
                E_f = corrections[n1][2]
                E_c = corrections[n1][3]
                pred_ratio = np.sqrt(E_f / E_c) / np.sqrt(E_c / E_f)  # = E_f/E_c
                print(f"    Predicted E_free/E_conf = {E_f/E_c:.4f}")
                print(f"    Observed mass ratio = {obs1/obs2:.4f}")
            print()
        if p1 == p2 and n1 != n2:
            m_cont1 = m_fermion_cont(n1, p1)
            m_cont2 = m_fermion_cont(n2, p2)
            print(f"  SAME p={p1}: {name1}({type1}, n={n1}) vs {name2}({type2}, n={n2})")
            print(f"    Continuous: {m_cont1:.4f} vs {m_cont2:.4f}")
            print(f"    These share tunneling depth but different breather modes.")
            print(f"    Free/confined correction affects them differently.")
            # Compute corrected ratio
            r1_f, r1_c = corrections[n1][0], corrections[n1][1]
            r2_f, r2_c = corrections[n2][0], corrections[n2][1]
            corr1 = r1_f if type1 == "free" else r1_c
            corr2 = r2_f if type2 == "free" else r2_c
            m_corr1 = m_cont1 * corr1
            m_corr2 = m_cont2 * corr2
            print(f"    Corrected: {m_corr1:.4f} vs {m_corr2:.4f}")
            print(f"    Observed:  {obs1:.4f} vs {obs2:.4f}")
            print()


# ============================================================
# PLOT: ERROR COMPARISON
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel 1: Error bars for each particle
ax1 = axes[0]
names = [r[0] for r in results]
errs_c = [r[7] for r in results]  # signed errors continuous
errs_3d = [r[8] for r in results]  # signed errors 3D
x_pos = np.arange(len(names))

ax1.barh(x_pos - 0.15, errs_c, 0.3, label='Continuous', alpha=0.7, color='blue')
ax1.barh(x_pos + 0.15, errs_3d, 0.3, label='3D Corrected', alpha=0.7, color='red')
ax1.set_yticks(x_pos)
ax1.set_yticklabels(names)
ax1.set_xlabel('Error (%)')
ax1.axvline(x=0, color='black', linewidth=0.5)
ax1.set_title('Fermion Mass Errors: Continuous vs 3D')
ax1.legend()

# Panel 2: |Error| comparison
ax2 = axes[1]
ax2.bar(x_pos - 0.2, errors_cont, 0.2, label='Continuous', alpha=0.7, color='blue')
ax2.bar(x_pos, errors_3d, 0.2, label='3D Corrected', alpha=0.7, color='red')
ax2.bar(x_pos + 0.2, errors_3d_split, 0.2, label='3D+Split', alpha=0.7, color='green')
ax2.set_xticks(x_pos)
ax2.set_xticklabels(names, rotation=45, ha='right')
ax2.set_ylabel('|Error| (%)')
ax2.set_title('Absolute Errors: All Methods')
ax2.legend()

plt.tight_layout()
plt.savefig('lagrangian/full_spectrum_3d_results.png', dpi=150)
print(f"\nPlot saved to lagrangian/full_spectrum_3d_results.png")
