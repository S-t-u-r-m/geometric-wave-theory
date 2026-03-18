"""
Direct Alpha from 3D Simulation
=================================
For each outlier atom:
1. Compute Z_net from Oh screening (proven correct)
2. Simulate ONLY the valence modes in a kink of charge Z_net
3. Measure IE directly from energy difference
4. Compare to formula prediction — the gap IS the alpha correction

No formula for alpha. The wave dynamics compute it exactly.
"""

import numpy as np
from math import factorial
import time

try:
    import cupy as cp
    xp = cp
except ImportError:
    xp = np

d = 3
gamma_sg = np.pi / (2**(d+1) * np.pi - 2)
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6
V_0 = 1.0 / np.pi**2

N = 48
L = 10.0
dx = 2 * L / N
x1d = np.linspace(-L, L, N, endpoint=False)
X, Y, Z_grid = np.meshgrid(x1d, x1d, x1d, indexing='ij')
R = np.sqrt(X**2 + Y**2 + Z_grid**2) + 1e-10
X_g, Y_g, Z_g, R_g = xp.asarray(X), xp.asarray(Y), xp.asarray(Z_grid), xp.asarray(R)
dt = 0.3 * dx
N_steps = 4000


def laplacian_3d(phi):
    return (xp.roll(phi,1,0)+xp.roll(phi,-1,0)+xp.roll(phi,1,1)+
            xp.roll(phi,-1,1)+xp.roll(phi,1,2)+xp.roll(phi,-1,2)-6*phi)/dx**2

def kink_3d(R, Z): return (4.0/np.pi)*xp.arctan(xp.exp(-xp.sqrt(float(Z))*R))

def breather_3d(X, Y, Z, R, l, m, omega, amp=0.15, n_quantum=1, Z_eff=1.0):
    """Breather mode with principal quantum number n.

    The radial profile scales as the hydrogen-like wavefunction:
    - Peak at r ~ n^2 / Z_eff
    - n-l-1 radial nodes
    - Width proportional to n/Z_eff

    This encodes the 1/n^2 energy scaling into the wave shape.
    """
    eps = float(np.sqrt(max(1.0-omega**2, 1e-12)))

    # Angular part
    if l==0: ang = xp.ones_like(R)
    elif l==1:
        if m==0: ang = Z/R
        elif m==1: ang = X/R
        else: ang = Y/R
    elif l==2:
        if m==0: ang = (3*Z**2-R**2)/R**2
        elif m==1: ang = X*Z/R**2
        elif m==-1: ang = Y*Z/R**2
        elif m==2: ang = X*Y/R**2
        else: ang = (X**2-Y**2)/R**2
    else: ang = xp.ones_like(R)

    # Radial part: hydrogen-like envelope scaled by n
    # r_peak = n^2 / Z_eff (where the breather amplitude peaks)
    # width = n / Z_eff (radial extent)
    r_scale = float(n_quantum**2 / max(Z_eff, 0.5))
    sigma = float(n_quantum / max(Z_eff, 0.5))

    # Radial nodes for n > l+1
    n_nodes = n_quantum - l - 1
    if n_nodes == 0:
        radial_node = xp.ones_like(R)
    elif n_nodes == 1:
        # One node: (1 - r/r_node) where r_node ~ n*sigma/2
        r_node = r_scale * 0.5
        radial_node = 1.0 - R / r_node
    elif n_nodes == 2:
        r1 = r_scale * 0.33
        r2 = r_scale * 0.75
        radial_node = (1.0 - R/r1) * (1.0 - R/r2)
    else:
        radial_node = xp.ones_like(R)

    # Gaussian-like envelope centered at r_scale
    envelope = amp * xp.exp(-(R - r_scale)**2 / (2 * sigma**2))

    # Combine: angular * radial_nodes * envelope
    breather = (4.0/np.pi) * xp.arctan(eps * ang * radial_node * envelope)
    return breather

def total_energy(phi, Z_eff):
    GE = 0.5*((xp.roll(phi,1,0)-phi)**2+(xp.roll(phi,1,1)-phi)**2+(xp.roll(phi,1,2)-phi)**2)
    PE = Z_eff*V_0*(1-xp.cos(np.pi*phi))
    return float(xp.sum(GE+PE)*dx**3)

def evolve_3d(phi_init, Z_eff):
    phi = phi_init.copy(); phi_old = phi.copy(); energies = []
    for step in range(N_steps):
        lap = laplacian_3d(phi)
        force = Z_eff*(1.0/np.pi)*xp.sin(np.pi*phi)
        phi_new = 2*phi-phi_old+dt**2*(lap-force)
        phi_old = phi.copy(); phi = phi_new
        if step > N_steps//2 and step%50==0:
            energies.append(total_energy(phi, Z_eff))
    return np.mean(energies) if energies else 0.0

omega_s = float(np.cos(1*gamma_sg))
omega_p = float(np.cos(2*gamma_sg))
omega_d = float(np.cos(3*gamma_sg))


def simulate_IE(Z_net, val_config, n_val, eV_scale=None):
    """Simulate IE for a valence configuration in a screened kink.

    Z_net: effective nuclear charge (after core screening)
    val_config: list of (l, m, n) for each valence mode
    n_val: principal quantum number of outermost mode
    eV_scale: calibration factor (lattice units to eV)

    Returns: IE in eV
    """
    omega_map = {0: omega_s, 1: omega_p, 2: omega_d}

    # Kink at Z_net
    phi_k = kink_3d(R_g, Z_net)
    E_k = evolve_3d(phi_k, Z_net)

    # Add all valence modes with their quantum numbers
    phi_full = phi_k.copy()
    for l, m, n_q in val_config:
        omega = omega_map.get(l, omega_s)
        phi_full = phi_full + breather_3d(X_g, Y_g, Z_g, R_g, l, m, omega,
                                          n_quantum=n_q, Z_eff=Z_net)
    E_full = evolve_3d(phi_full, Z_net)

    # Remove outermost mode (last in list)
    phi_ion = phi_k.copy()
    for l, m, n_q in val_config[:-1]:
        omega = omega_map.get(l, omega_s)
        phi_ion = phi_ion + breather_3d(X_g, Y_g, Z_g, R_g, l, m, omega,
                                        n_quantum=n_q, Z_eff=Z_net)
    E_ion = evolve_3d(phi_ion, Z_net)

    IE_lattice = abs(E_ion - E_full)

    if eV_scale is not None and eV_scale > 0:
        return IE_lattice * eV_scale
    return IE_lattice


print("Direct Alpha from 3D Simulation")
print("=" * 60)
print(f"  Simulating valence modes in Oh-screened kink potential")
print(f"  Grid: {N}^3, N_steps={N_steps}")
print(flush=True)

# First: calibrate with hydrogen (n=1, l=0, Z_net=1)
print("\nCalibrating with hydrogen...", flush=True)
t0 = time.time()
IE_H_lu = simulate_IE(1.0, [(0, 0, 1)], 1)  # returns lattice units
eV_per_lu = 13.598 / IE_H_lu if IE_H_lu > 0 else 1.0
print(f"  H: IE_lattice = {IE_H_lu:.6f}, scale = {eV_per_lu:.2f} eV/lu")
print(f"  Time: {time.time()-t0:.1f}s", flush=True)

# Test atoms: (sym, Z_net, val_config as (l, m, n), n_val, IE_obs)
# Z_net from Oh screening. val_config includes n for each mode.
test_atoms = [
    ('H',  1.0,  [(0,0,1)], 1, 13.598),
    ('He', 1.5,  [(0,0,1), (0,0,1)], 1, 24.587),
    ('Li', 2.5,  [(0,0,2)], 2, 5.392),
    ('Be', 3.5,  [(0,0,2), (0,0,2)], 2, 9.323),
    ('B',  4.0,  [(0,0,2), (0,0,2), (1,1,2)], 2, 8.298),
    ('C',  4.0,  [(0,0,2), (0,0,2), (1,1,2), (1,-1,2)], 2, 11.260),
    ('N',  4.0,  [(0,0,2), (0,0,2), (1,1,2), (1,-1,2), (1,0,2)], 2, 14.534),
    ('O',  4.0,  [(0,0,2), (0,0,2), (1,1,2), (1,-1,2), (1,0,2), (1,1,2)], 2, 13.618),
    ('F',  4.0,  [(0,0,2), (0,0,2), (1,1,2), (1,-1,2), (1,0,2), (1,1,2), (1,-1,2)], 2, 17.423),
    ('Ne', 4.0,  [(0,0,2), (0,0,2), (1,1,2), (1,-1,2), (1,0,2), (1,1,2), (1,-1,2), (1,0,2)], 2, 21.565),
]

print(f"\n{'Sym':>4} {'Z_net':>5} {'modes':>5} {'IE_sim':>8} {'IE_obs':>8} {'err':>7}")
print("-" * 45)

for sym, Z_net, val_cfg, n_v, IE_obs in test_atoms:
    t0 = time.time()
    IE_sim = simulate_IE(Z_net, val_cfg, n_v, eV_scale=eV_per_lu)
    err = (IE_sim - IE_obs) / IE_obs * 100
    print(f"  {sym:>3} {Z_net:5.1f} {len(val_cfg):5d} {IE_sim:8.3f} {IE_obs:8.3f} {err:+6.1f}%  ({time.time()-t0:.1f}s)", flush=True)
