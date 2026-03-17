"""
1D Sine-Gordon Coupling Matrix — Extended
==========================================
Measure the full coupling matrix between breather harmonics n=1..7
at multiple Z values. Write results incrementally to file.

Goal: does the coupling pattern reveal a closed-form expression
that unifies the v19 screening rules?
"""

import numpy as np
from math import factorial
import time, os

d = 3
gamma_sg = np.pi / (2**(d+1) * np.pi - 2)
V_0 = 1.0 / np.pi**2

# Grid
Nx = 2000
x_max = 30.0
dx = 2 * x_max / Nx
x = np.linspace(-x_max, x_max, Nx)
dt = 0.5 * dx
N_steps = 15000  # slightly fewer for speed

_Z_well = 1.0

def set_Z(Z):
    global _Z_well
    _Z_well = Z

def kink(x, Z):
    return (4.0 / np.pi) * np.arctan(np.exp(np.sqrt(Z) * x))

def add_breather(phi_bg, x, omega, x0=0.0):
    eps = np.sqrt(max(1.0 - omega**2, 1e-12))
    num = eps
    den = omega * np.cosh(eps * (x - x0))
    return phi_bg + (4.0 / np.pi) * np.arctan(num / (den + 1e-30))

def evolve(phi_init, N_steps, dt, dx):
    phi = phi_init.copy()
    phi_old = phi.copy()
    energies = []
    for step in range(N_steps):
        lap = np.zeros_like(phi)
        lap[1:-1] = (phi[2:] - 2*phi[1:-1] + phi[:-2]) / dx**2
        force = _Z_well * (1.0 / np.pi) * np.sin(np.pi * phi)
        phi_new = 2*phi - phi_old + dt**2 * (lap - force)
        phi_new[0] = phi_init[0]
        phi_new[-1] = phi_init[-1]
        phi_old = phi.copy()
        phi = phi_new
        if step > N_steps // 2 and step % 100 == 0:
            KE = 0.5 * np.sum(((phi - phi_old)/dt)**2) * dx
            GE = 0.5 * np.sum(((phi[1:] - phi[:-1])/dx)**2) * dx
            PE = np.sum(_Z_well * V_0 * (1 - np.cos(np.pi * phi))) * dx
            energies.append(KE + GE + PE)
    return np.mean(energies) if energies else 0.0

# Output file
outfile = os.path.join(os.path.dirname(__file__), 'coupling_matrix_results.txt')

with open(outfile, 'w') as f:
    f.write("1D Sine-Gordon Coupling Matrix\n")
    f.write(f"gamma = {gamma_sg:.6f}, Nx = {Nx}, N_steps = {N_steps}\n")
    f.write("=" * 70 + "\n\n")
    f.flush()

    # Test Z values
    Z_values = [3.0, 5.0, 10.0]
    n_max = 7

    for Z in Z_values:
        set_Z(Z)
        f.write(f"\n--- Z = {Z:.1f} ---\n")
        f.write(f"{'n':>3} {'omega':>10} {'E_single':>10}\n")
        f.flush()
        print(f"\n=== Z = {Z:.1f} ===")

        # Single breather energies
        phi_k = kink(x, Z)
        E_k = evolve(phi_k, N_steps, dt, dx)

        E_single = {}
        for n in range(1, n_max + 1):
            omega = np.cos(n * gamma_sg)
            if omega <= 0.01:
                break
            phi_kb = add_breather(phi_k, x, omega)
            E_kb = evolve(phi_kb, N_steps, dt, dx)
            E_single[n] = E_kb - E_k
            f.write(f"{n:3d} {omega:10.6f} {E_single[n]:10.4f}\n")
            f.flush()
            print(f"  n={n}: E_single = {E_single[n]:.4f}")

        # Coupling matrix
        f.write(f"\nCoupling matrix C(n1,n2) = E_2nd(n2|n1) / E_single(n2):\n")
        header = "n1\\n2" + "".join(f"  {n:>7}" for n in E_single.keys())
        f.write(header + "\n")
        f.write("-" * len(header) + "\n")
        f.flush()

        print(f"\n  Coupling matrix:")
        coupling = {}

        for n1 in E_single:
            omega1 = np.cos(n1 * gamma_sg)
            phi_k1 = add_breather(phi_k, x, omega1)
            E_k1 = evolve(phi_k1, N_steps, dt, dx)

            row = f"  {n1:2d}  "
            for n2 in E_single:
                omega2 = np.cos(n2 * gamma_sg)
                phi_k2 = add_breather(phi_k1, x, omega2)
                E_k2 = evolve(phi_k2, N_steps, dt, dx)

                E_2nd = E_k2 - E_k1
                C = E_2nd / E_single[n2] if abs(E_single[n2]) > 1e-10 else 0
                coupling[(n1, n2)] = C
                row += f"  {C:7.3f}"

            f.write(row + "\n")
            f.flush()
            print(row)

        # Asymmetry ratios: C(n1,n2) / C(n2,n1)
        f.write(f"\nAsymmetry ratios C(n1->n2) / C(n2->n1):\n")
        print(f"\n  Asymmetry ratios:")
        for n1 in E_single:
            for n2 in E_single:
                if n2 > n1:
                    c12 = coupling.get((n1,n2), 0)
                    c21 = coupling.get((n2,n1), 0)
                    ratio = (c12 - 1) / (c21 - 1) if abs(c21 - 1) > 0.01 else 0
                    line = f"  W({n1}->{n2})/W({n2}->{n1}) = {ratio:.4f}"
                    # Compare to cos(dl*pi/d)
                    dl = n2 - n1
                    predicted = np.cos(dl * np.pi / d)
                    line += f"  cos({dl}*pi/{d}) = {predicted:.4f}"
                    f.write(line + "\n")
                    print(line)

        # Diagonal trend
        f.write(f"\nDiagonal (same-channel) coupling:\n")
        print(f"\n  Diagonal coupling:")
        for n in E_single:
            C_nn = coupling.get((n,n), 0)
            line = f"  C({n},{n}) = {C_nn:.4f}  W = {C_nn-1:.4f}"
            f.write(line + "\n")
            print(line)

        f.write("\n")
        f.flush()

    # Final analysis
    f.write("\n" + "=" * 70 + "\n")
    f.write("SUMMARY\n")
    f.write("=" * 70 + "\n")
    f.flush()

print(f"\nResults written to: {outfile}")
print("Done!")
