"""
Ionization Energy from Hessian Eigenvalues
============================================
Replace the algebraic Z_eff^alpha formula with EXACT bound state energies
from the kink well Hessian on the discrete lattice.

The atom = kink well of depth proportional to Z.
Inner electrons = breathers that modify the well potential.
Ionization energy = energy of the outermost bound state.

The Poeschl-Teller well from the sine-Gordon Lagrangian:
  U(x) = -Z * (2/pi^2) / cosh^2(x)

On a discrete lattice, the Hessian is:
  H_ij = (2 + cos(pi*phi_kink_i)) * delta_ij - delta_{i,j+1} - delta_{i,j-1}

where phi_kink encodes the nuclear charge (well depth) and screening
from inner electrons.

Key insight from today: the kink well parameters (s, V_0) are derived
from the Lagrangian, and the bound state energies emerge from the
Hessian eigenvalues — same method that derived the bond energy.
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from math import factorial
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3

# GWT constants
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/PI**2 + np.log(2*d)))
E_H_eV = alpha_em**2 / 2 * 0.511e6  # 13.604 eV

# Poeschl-Teller
s = (-1 + np.sqrt(1 + 8/PI**2)) / 2  # 0.17279
V_0 = 2/PI**2  # well depth per unit Z

outfile = os.path.join(os.path.dirname(__file__), "ie_hessian_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("IONIZATION ENERGY FROM HESSIAN EIGENVALUES")
report("=" * 70)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report(f"s = {s:.8f}, V_0 = 2/pi^2 = {V_0:.8f}")
report(f"E_H = {E_H_eV:.4f} eV, alpha = {alpha_em:.8f}")
report("")

N = 512  # lattice sites
center = N // 2
x = np.arange(N, dtype=np.float64) - center

# ============================================================
# THE ATOMIC POTENTIAL
# ============================================================
def atomic_potential(Z, screening_profile=None):
    """Build the on-site potential for an atom with nuclear charge Z.

    The kink (nucleus) creates a Poeschl-Teller well:
      cos(pi*phi_kink) where phi_kink = (4/pi)*arctan(exp(-|x|*beta))

    For a bare nucleus: beta = sqrt(Z) * s_eff, depth ~ Z
    Inner electrons screen by reducing the effective well depth.

    Returns: the diagonal of the Hessian (2 + V(x)) at each site.
    """
    # Bare nuclear well: kink profile with Z-dependent depth
    # The kink amplitude scales with Z: larger Z = deeper well
    # Use: phi(x) = (4/pi) * arctan(1/cosh(sqrt(Z) * x))
    # This gives a well of depth ~ Z * V_0 = Z * 2/pi^2
    beta = np.sqrt(Z)
    phi_nuc = (4.0/PI) * np.arctan(1.0 / np.cosh(beta * x + 1e-30))

    # The on-site Hessian term: cos(pi*phi)
    # At phi=0 (far from nucleus): cos(0) = 1 (mass gap)
    # At phi=1 (kink peak): cos(pi) = -1 (deep well)
    V_site = np.cos(PI * phi_nuc)

    if screening_profile is not None:
        # Inner electrons modify the potential:
        # Each inner electron acts as a breather that partially fills the well.
        # The effective potential becomes: cos(pi*(phi_nuc - phi_screen))
        # where phi_screen is the screening field from inner electrons.
        V_site = np.cos(PI * (phi_nuc - screening_profile))

    return V_site

def build_atom_hessian(Z, screening_profile=None):
    """Build the full Hessian for an atom."""
    V = atomic_potential(Z, screening_profile)

    # H_ij = (2 + V_i) * delta_ij - delta_{i,j+1} - delta_{i,j-1}
    diag = 2.0 + V
    off = -np.ones(N - 1)
    H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
    H[0, N-1] = -1.0
    H[N-1, 0] = -1.0
    return H.tocsr()

def screening_from_electrons(Z, electron_config):
    """Build the screening profile from inner electrons.

    Each occupied shell contributes a breather-like screening field.
    The screening field has the shape of the shell's wavefunction,
    scaled by the number of electrons and the Oh coupling weight.

    electron_config: list of (n, l, count) tuples
    The OUTERMOST shell is what we're ionizing — all others screen.
    """
    if len(electron_config) <= 1:
        return None  # H or single-electron: no screening

    # The screening profile: inner electrons partially cancel the nuclear potential
    # Each inner electron at shell (n, l) with Z_eff contributes:
    #   phi_screen(x) ~ count * (4/pi) * arctan(1/cosh(sqrt(Z_eff_n) * x / n))
    # where Z_eff_n is the effective nuclear charge seen by that shell.

    # Simple model: each shell sees Z minus screening from deeper shells.
    # Build up from inside out.
    phi_screen = np.zeros(N)
    Z_remaining = float(Z)

    # All shells except the last (which is being ionized)
    inner_shells = electron_config[:-1]

    for n, l, count in inner_shells:
        # Effective Z for this shell: the nuclear charge minus screening from deeper shells
        # For simplicity, use Z_remaining (decremented after each shell)
        Z_eff = max(Z_remaining - count/2, 1.0)  # rough: each electron screens half

        # The screening field from this shell:
        # Shape: breather-like, width ~ n/sqrt(Z_eff), amplitude ~ count/Z
        beta_n = np.sqrt(Z_eff) / n  # spatial scale
        amplitude = count / Z  # fraction of nuclear charge screened

        phi_screen += amplitude * (4.0/PI) * np.arctan(1.0 / np.cosh(beta_n * x + 1e-30))

        Z_remaining -= count * 0.85  # rough screening: each electron screens ~0.85

    return phi_screen

def get_ionization_energy(Z, electron_config):
    """Compute ionization energy from Hessian eigenvalues.

    Returns IE in lattice units (omega^2 of the outermost bound state).
    """
    screening = screening_from_electrons(Z, electron_config)
    H = build_atom_hessian(Z, screening)

    # Find lowest eigenvalues
    n_eig = min(20, N - 2)
    try:
        evals, evecs = eigsh(H, k=n_eig, which='SM')
        evals = np.sort(evals)
    except:
        return None, 0

    # Count bound states (below mass gap omega^2 = 1)
    bound = evals[evals < 1.0]
    n_bound = len(bound)

    if n_bound == 0:
        return None, 0

    # The ionization energy = difference between outermost bound state and the continuum
    # Continuum starts at omega^2 = 1 (mass gap)
    outermost = bound[-1]  # highest bound state = outermost electron
    IE_lattice = 1.0 - outermost  # energy to reach the continuum

    return IE_lattice, n_bound

# ============================================================
# TEST ATOMS
# ============================================================
report("PERIOD 1-2 ATOMS")
report("-" * 65)

# Electron configs: (n, l, count)
atoms = [
    (1,  'H',  13.598, [(1,0,1)]),
    (2,  'He', 24.587, [(1,0,2)]),
    (3,  'Li',  5.392, [(1,0,2), (2,0,1)]),
    (4,  'Be',  9.323, [(1,0,2), (2,0,2)]),
    (5,  'B',   8.298, [(1,0,2), (2,0,2), (2,1,1)]),
    (6,  'C',  11.260, [(1,0,2), (2,0,2), (2,1,2)]),
    (7,  'N',  14.534, [(1,0,2), (2,0,2), (2,1,3)]),
    (8,  'O',  13.618, [(1,0,2), (2,0,2), (2,1,4)]),
    (9,  'F',  17.423, [(1,0,2), (2,0,2), (2,1,5)]),
    (10, 'Ne', 21.565, [(1,0,2), (2,0,2), (2,1,6)]),
    (11, 'Na',  5.139, [(1,0,2), (2,0,2), (2,1,6), (3,0,1)]),
    (12, 'Mg',  7.646, [(1,0,2), (2,0,2), (2,1,6), (3,0,2)]),
    (13, 'Al',  5.986, [(1,0,2), (2,0,2), (2,1,6), (3,0,2), (3,1,1)]),
    (18, 'Ar', 15.760, [(1,0,2), (2,0,2), (2,1,6), (3,0,2), (3,1,6)]),
]

# Also test outliers
outlier_atoms = [
    (30, 'Zn', 9.394, [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2)]),
    (46, 'Pd', 8.337, [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10)]),
    (48, 'Cd', 8.994, [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(5,0,2)]),
]

report(f"{'Z':>4} {'Sym':>3} {'IE_obs':>8} {'IE_hess':>10} {'err':>8} {'n_bound':>8}")
report("-" * 50)

# First, calibrate: for H (Z=1, no screening), the Hessian should give IE
# Let's see what the raw eigenvalue gives
for Z, sym, IE_obs, config in atoms + outlier_atoms:
    IE_lattice, n_bound = get_ionization_energy(Z, config)

    if IE_lattice is not None and IE_lattice > 0:
        # Convert to eV: we need a scale factor
        # For hydrogen: IE_lattice should map to 13.598 eV
        # We'll calibrate with H and use the same scale for all
        if Z == 1:
            scale_H = IE_obs / IE_lattice
            report(f"{Z:4d} {sym:>3} {IE_obs:8.3f} {IE_lattice:10.6f} {'(cal)':>8} {n_bound:8d}")
            report(f"     Scale: {scale_H:.4f} eV per lattice unit")
        else:
            IE_pred = IE_lattice * scale_H
            err = (IE_pred - IE_obs) / IE_obs * 100
            report(f"{Z:4d} {sym:>3} {IE_obs:8.3f} {IE_pred:10.3f} {err:+7.1f}% {n_bound:8d}")
    else:
        report(f"{Z:4d} {sym:>3} {IE_obs:8.3f} {'FAIL':>10} {'---':>8} {n_bound:8d}")

report("")

# ============================================================
# ANALYSIS
# ============================================================
report("ANALYSIS")
report("-" * 65)
report("")
report("The Hessian method computes bound state energies EXACTLY for a given potential.")
report("The key physics is in HOW the screening potential is built:")
report("  - Nuclear well: kink with depth ~ Z")
report("  - Inner electrons: breathers that partially fill the well")
report("  - The screening fraction per electron determines the IE")
report("")
report("The current V20 model uses algebraic alpha exponents for Z_eff.")
report("The Hessian method replaces this with exact eigenvalues.")
report("The screening MODEL determines the quality of the result.")
report("")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))

log.close()
print(f"\nResults saved to: {outfile}")
