"""
Two-Center Standing Wave Overlap Simulation
=============================================
Instead of approximating overlap with sin(kR), compute it directly
from the actual wave profiles.

Each atom creates a standing wave: psi(r) ~ r^l * exp(-Z*r/n) * sin(kr)
(hydrogen-like radial wavefunction, simplified).

The bond coupling = integral of psi_1(r) * psi_2(R-r) dr
This naturally includes:
  - Exponential decay (no need for separate envelope)
  - Node structure (radial nodes reduce overlap properly)
  - Phase behavior (no sin(pi)=0 coincidence)

Compare the exact overlap to our |sin(phase)| approximation.
"""

import numpy as np
from scipy import integrate

pi = np.pi
E_H = 13.6057
d = 3

C_bond  = pi / d
f_pi    = d**2 / (d**2 + 1)
alpha   = 1 - f_pi / d
beta    = (1 + f_pi) / 2
f_anti  = 2*d / (2*d - 1)
c_ionic = 1.0 / (2*d + 1)

Z_eff = {
    'H_1s': 1.0000, 'Li_2s': 1.2792, 'B_2p': 2.4214,
    'C_2p': 3.1358, 'N_2p': 3.8340, 'O_2p': 4.4532,
    'F_2p': 5.0998, 'Na_3s': 2.5074, 'Cl_3p': 6.1161,
}

def get_n(orb): return int(orb.split('_')[1][0])
def get_l(orb): return {'s': 0, 'p': 1}[orb.split('_')[1][1]]
def has_nodes(orb): return min(get_n(orb) - get_l(orb) - 1, 1)
def orbital_energy(orb): return E_H * (Z_eff[orb] / get_n(orb))**2


# =============================================================================
# HYDROGEN-LIKE RADIAL WAVEFUNCTIONS
# =============================================================================

def hydrogen_radial(r, n, l, Z):
    """Simplified hydrogen-like radial wavefunction R_nl(r).

    Uses the exact forms for low n,l:
      1s: R = 2*Z^(3/2) * exp(-Zr)
      2s: R = (Z/2)^(3/2) * (1/(2*sqrt(2))) * (2 - Zr) * exp(-Zr/2)
      2p: R = (Z/2)^(3/2) * (1/(2*sqrt(6))) * Zr * exp(-Zr/2)
      3s: R = (Z/3)^(3/2) * (2/(9*sqrt(3))) * (3 - 2Zr/3 + 2(Zr)^2/27) * exp(-Zr/3)
      3p: R = (Z/3)^(3/2) * (4/(9*sqrt(6))) * Zr/3 * (1 - Zr/6) * exp(-Zr/3)

    Returns R_nl(r) * r (the radial probability amplitude, unnormalized).
    We normalize numerically afterwards.
    """
    rho = Z * r / n  # scaled variable

    if n == 1 and l == 0:  # 1s
        return np.exp(-rho) * r
    elif n == 2 and l == 0:  # 2s
        return (2 - 2*rho) * np.exp(-rho) * r
    elif n == 2 and l == 1:  # 2p
        return rho * np.exp(-rho) * r
    elif n == 3 and l == 0:  # 3s
        return (3 - 6*rho + 2*rho**2) * np.exp(-rho) * r
    elif n == 3 and l == 1:  # 3p
        return rho * (2 - rho) * np.exp(-rho) * r
    else:
        raise ValueError(f"n={n}, l={l} not implemented")


def make_wavefunction(orb, center=0, direction=1):
    """Create a 1D wavefunction centered at 'center'.

    Returns a function psi(x) that gives the wavefunction value at position x.
    'direction' = +1 or -1 (which side the bonding lobe faces).
    """
    n = get_n(orb)
    l = get_l(orb)
    Z = Z_eff[orb]

    def psi(x):
        r = np.abs(x - center)
        r = np.maximum(r, 1e-10)  # avoid r=0
        val = hydrogen_radial(r, n, l, Z)

        # For p orbitals, include angular sign (bonding lobe direction)
        if l == 1:
            sign = np.sign(x - center) * direction
            val = val * sign

        return val

    return psi


def compute_overlap(orb1, orb2, R, bond_type='sigma'):
    """Compute the overlap integral of two wavefunctions at separation R.

    S = integral psi_1(x) * psi_2(x) dx

    where psi_1 is centered at x=0 and psi_2 at x=R.

    For sigma bonds: both wavefunctions point along bond axis
    For pi bonds: both perpendicular (overlap is reduced by geometry)
    """
    # Create wavefunctions
    # Atom 1 at origin, bonding lobe toward atom 2 (+x direction)
    psi1 = make_wavefunction(orb1, center=0, direction=+1)
    # Atom 2 at x=R, bonding lobe toward atom 1 (-x direction)
    psi2 = make_wavefunction(orb2, center=R, direction=-1)

    # Normalize each wavefunction
    def norm_integrand1(x):
        return psi1(x)**2
    def norm_integrand2(x):
        return psi2(x)**2

    # Integration range: well beyond both atoms
    x_max = R + 20 * max(get_n(orb1)**2/Z_eff[orb1], get_n(orb2)**2/Z_eff[orb2])
    x_min = -20 * get_n(orb1)**2/Z_eff[orb1]

    N1_sq, _ = integrate.quad(norm_integrand1, x_min, x_max, limit=200)
    N2_sq, _ = integrate.quad(norm_integrand2, x_min, x_max, limit=200)

    N1 = np.sqrt(N1_sq) if N1_sq > 0 else 1
    N2 = np.sqrt(N2_sq) if N2_sq > 0 else 1

    # Overlap integral
    def overlap_integrand(x):
        return psi1(x) * psi2(x) / (N1 * N2)

    S, err = integrate.quad(overlap_integrand, x_min, x_max, limit=200)

    return S


def compute_overlap_energy(orb1, orb2, R, bond_type='sigma'):
    """Compute the resonance/coupling integral (energy-weighted overlap).

    H_12 = integral psi_1(x) * V(x) * psi_2(x) dx

    where V(x) is the potential from the other nucleus.
    For atom 2's potential acting on atom 1's electron:
    V(x) = -Z_2 / |x - R|  (Coulomb, in atomic units -> multiply by E_H)
    """
    psi1 = make_wavefunction(orb1, center=0, direction=+1)
    psi2 = make_wavefunction(orb2, center=R, direction=-1)

    x_max = R + 20 * max(get_n(orb1)**2/Z_eff[orb1], get_n(orb2)**2/Z_eff[orb2])
    x_min = -20 * get_n(orb1)**2/Z_eff[orb1]

    N1_sq, _ = integrate.quad(lambda x: psi1(x)**2, x_min, x_max, limit=200)
    N2_sq, _ = integrate.quad(lambda x: psi2(x)**2, x_min, x_max, limit=200)
    N1 = np.sqrt(max(N1_sq, 1e-30))
    N2 = np.sqrt(max(N2_sq, 1e-30))

    Z2 = Z_eff[orb2]

    def coupling_integrand(x):
        r2 = abs(x - R)
        if r2 < 0.01:
            r2 = 0.01
        V = -Z2 / r2  # Coulomb potential from atom 2
        return psi1(x) * V * psi2(x) / (N1 * N2)

    H12, err = integrate.quad(coupling_integrand, x_min, x_max, limit=200)

    return H12 * E_H  # convert to eV


# =============================================================================
# MOLECULE DATA
# =============================================================================

molecules = [
    ('H2',   1.401,  4.745, [('ss',1)],                                         'H_1s',  'H_1s'),
    ('Li2',  5.051,  1.056, [('ss',1)],                                         'Li_2s', 'Li_2s'),
    ('B2',   3.005,  3.02,  [('pi',2)],                                         'B_2p',  'B_2p'),
    ('C2',   2.348,  6.32,  [('pi',2),('sp_sigma',1),('sp_sigma_anti',1)],      'C_2p',  'C_2p'),
    ('N2',   2.074,  9.759, [('pp_sigma',1),('pi',2)],                          'N_2p',  'N_2p'),
    ('O2',   2.282,  5.213, [('pp_sigma',1),('pi',2),('pi_anti',1)],            'O_2p',  'O_2p'),
    ('F2',   2.668,  1.660, [('pp_sigma',1),('pi',2),('pi_anti',2)],            'F_2p',  'F_2p'),
    ('Na2',  5.818,  0.746, [('ss',1)],                                         'Na_3s', 'Na_3s'),
    ('Cl2',  3.757,  2.514, [('pp_sigma',1),('pi',2),('pi_anti',2)],            'Cl_3p', 'Cl_3p'),
    ('HF',   1.733,  5.869, [('sp',1)],                                         'H_1s',  'F_2p'),
    ('CO',   2.132, 11.225, [('pp_sigma',1),('pi',2)],                          'C_2p',  'O_2p'),
    ('NO',   2.175,  6.497, [('pp_sigma',1),('pi',2),('pi_anti',1)],            'N_2p',  'O_2p'),
    ('OH',   1.834,  4.392, [('sp',1)],                                         'O_2p',  'H_1s'),
    ('HCl',  2.409,  4.434, [('sp',1)],                                         'H_1s',  'Cl_3p'),
    ('LiH',  3.015,  2.515, [('ss',1)],                                         'Li_2s', 'H_1s'),
    ('LiF',  2.955,  5.939, [('sp',1)],                                         'Li_2s', 'F_2p'),
    ('BH',   2.329,  3.42,  [('sp',1)],                                         'B_2p',  'H_1s'),
    ('CH',   2.116,  3.47,  [('sp',1)],                                         'C_2p',  'H_1s'),
    ('NH',   1.958,  3.57,  [('sp',1)],                                         'N_2p',  'H_1s'),
    ('BF',   2.386,  7.81,  [('pp_sigma',1),('pi',2)],                          'B_2p',  'F_2p'),
    ('CN',   2.214,  7.72,  [('pp_sigma',0.5),('pi',2)],                        'C_2p',  'N_2p'),
    ('NaH',  3.566,  1.97,  [('ss',1)],                                         'Na_3s', 'H_1s'),
    ('NaCl', 4.461,  4.23,  [('sp',1)],                                         'Na_3s', 'Cl_3p'),
    ('H2O',  1.809,  5.117, [('sp',1)],                                         'O_2p',  'H_1s'),
]


# =============================================================================
# COMPUTE OVERLAPS AND COMPARE TO sin(phase)
# =============================================================================

print("=" * 90)
print("  TWO-CENTER WAVE OVERLAP SIMULATION")
print("  Compare exact overlap integral to |sin(phase)| approximation")
print("=" * 90)

# First: compute overlaps for all molecules
print(f"\n{'Mol':<7} {'R':>6} {'S_exact':>9} {'|sin(ph)|':>10} {'ratio':>7} "
      f"{'H12_eV':>8} {'C*Es*|sin|':>10} {'De_exp':>7}")
print("-" * 80)

for mol in molecules:
    name, R, De_exp, bonds, orb1, orb2 = mol

    # Exact overlap
    S = compute_overlap(orb1, orb2, R)

    # Our formula's |sin(phase)|
    n1, n2 = get_n(orb1), get_n(orb2)
    h1, h2 = has_nodes(orb1), has_nodes(orb2)
    b1 = 1 + beta * h1
    b2 = 1 + beta * h2
    phase = R / n1**b1 + R / n2**b2
    sin_val = abs(np.sin(phase))

    # E_scale
    l1, l2 = get_l(orb1), get_l(orb2)
    a1 = 2 + (1-2*l1)*alpha*h1
    a2 = 2 + (1-2*l2)*alpha*h2
    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)

    # Our covalent formula (single bond)
    cov_formula = C_bond * E_scale * sin_val

    # Coupling integral
    H12 = compute_overlap_energy(orb1, orb2, R)

    ratio = abs(S) / sin_val if sin_val > 0.001 else float('inf')

    print(f"{name:<7} {R:6.3f} {S:9.4f} {sin_val:10.4f} {ratio:7.3f} "
          f"{H12:8.3f} {cov_formula:10.3f} {De_exp:7.3f}")


# =============================================================================
# SCAN: Overlap vs distance for key outliers
# =============================================================================
print(f"\n{'='*90}")
print(f"  OVERLAP vs DISTANCE for outlier molecules")
print(f"{'='*90}")

outlier_pairs = [
    ('CH',  'C_2p', 'H_1s', 2.116),
    ('LiH', 'Li_2s', 'H_1s', 3.015),
    ('BH',  'B_2p', 'H_1s', 2.329),
    ('NaH', 'Na_3s', 'H_1s', 3.566),
    ('LiF', 'Li_2s', 'F_2p', 2.955),
    ('BF',  'B_2p', 'F_2p', 2.386),
]

for mol_name, o1, o2, R_eq in outlier_pairs:
    n1, n2 = get_n(o1), get_n(o2)
    h1, h2 = has_nodes(o1), has_nodes(o2)
    b1 = 1 + beta * h1
    b2 = 1 + beta * h2

    print(f"\n  {mol_name} ({o1} + {o2}, R_eq = {R_eq:.3f}):")
    print(f"  {'R':>6} {'S_exact':>9} {'|sin(ph)|':>10} {'phase/pi':>9}")

    for R in np.arange(1.0, 5.5, 0.5):
        S = compute_overlap(o1, o2, R)
        phase = R / n1**b1 + R / n2**b2
        sin_val = abs(np.sin(phase))
        print(f"  {R:6.2f} {S:9.4f} {sin_val:10.4f} {phase/pi:9.3f}")


# =============================================================================
# KEY TEST: What C_bond would make the exact overlap match De_exp?
# =============================================================================
print(f"\n{'='*90}")
print(f"  CALIBRATION: What prefactor makes exact overlap predict De_exp?")
print(f"{'='*90}")

print(f"\n  If D_cov = C_eff * E_scale * |S_exact|, what C_eff matches experiment?")
print(f"  (Compare to our C = pi/3 = {C_bond:.4f})")
print(f"\n{'Mol':<7} {'|S|':>8} {'E_scale':>8} {'De_exp':>7} {'D_ion':>6} {'D_cov_need':>10} {'C_eff':>7}")
print("-" * 60)

for mol in molecules:
    name, R, De_exp, bonds, orb1, orb2 = mol

    # Only single bonds for clean comparison
    if len(bonds) > 1:
        continue

    S = abs(compute_overlap(orb1, orb2, R))

    n1, l1 = get_n(orb1), get_l(orb1)
    n2, l2 = get_n(orb2), get_l(orb2)
    h1, h2 = has_nodes(orb1), has_nodes(orb2)
    a1 = 2 + (1-2*l1)*alpha*h1; a2 = 2 + (1-2*l2)*alpha*h2
    E_scale = np.sqrt(E_H/n1**a1 * E_H/n2**a2)

    # Compute ionic part
    eps1, eps2 = orbital_energy(orb1), orbital_energy(orb2)
    delta = abs(eps1 - eps2)

    # Need D_cov = De_exp - D_ion
    # But D_ion depends on D_cov... use iterative approach
    D_cov_need = De_exp * 0.7  # initial guess
    for _ in range(20):
        V = max(abs(D_cov_need), 0.01)
        q = delta / np.sqrt(delta**2 + (2*V)**2) if delta > 0 else 0
        D_ion = c_ionic * q**2 * 2 * E_H / R
        D_cov_need = De_exp - D_ion

    if S > 0.001 and E_scale > 0.01:
        C_eff = D_cov_need / (E_scale * S)
    else:
        C_eff = float('inf')

    print(f"{name:<7} {S:8.4f} {E_scale:8.3f} {De_exp:7.3f} {D_ion:6.3f} {D_cov_need:10.3f} {C_eff:7.3f}")
