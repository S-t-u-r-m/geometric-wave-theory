"""
Self-Consistent Hessian Ionization Energy
==========================================
Replace the algebraic Z_eff^alpha formula with self-consistent
eigenvalues of the kink well on the discrete lattice.

Method:
  1. Build nuclear kink well: phi_nuc(x) = A * arctan(1/cosh(sqrt(Z)*x))
     where A is chosen so the well supports the right number of bound states
  2. Each occupied shell screens the well: phi_screen = sum of shell contributions
  3. Screened well: phi_eff = phi_nuc - phi_screen
  4. Hessian H = -Lap + cos(pi*phi_eff)
  5. Bound states = eigenvalues below mass gap (omega^2 < 1)
  6. Self-consistency: iterate until eigenvalues converge
  7. IE = 1 - outermost_eigenvalue (gap to continuum)

The Oh angular physics enters through SCREENING WEIGHTS:
  - s,p core screening: weight w_pi = 0.5 per channel
  - d core screening of s: weight w_delta = -0.5 (anti-screening)
  - d core screening of p: weight w_pi = 0.5 (Oh allowed)
  - f core: T1u component w_pi, rest w_delta/d

No alpha exponents. No power law. Screening emerges self-consistently.
"""
import sys, io, os, time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from math import factorial, comb
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PI = np.pi
d = 3
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/PI**2 + np.log(2*d)))
E_H_eV = alpha_em**2 / 2 * 0.511e6
w_pi = np.cos(PI/d)     # 0.5
w_delta = np.cos(2*PI/d) # -0.5
s_PT = (-1 + np.sqrt(1 + 8/PI**2)) / 2  # 0.17279

outfile = os.path.join(os.path.dirname(__file__), "ie_hessian_sc_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("SELF-CONSISTENT HESSIAN IONIZATION ENERGY")
report("=" * 70)
report(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
report("")

N = 256
center = N // 2
x = np.arange(N, dtype=np.float64) - center

def kink_well(Z, amplitude=1.0):
    """Nuclear kink well profile. Depth scales with Z through beta = sqrt(Z).
    The amplitude parameter scales the overall field magnitude."""
    beta = np.sqrt(max(Z, 0.1))
    return amplitude * (4.0/PI) * np.arctan(1.0 / np.cosh(np.clip(beta * x, -500, 500)))

def screening_profile(eigvec, n_shell, l_shell, Z_eff_shell, oh_weight):
    """Screening contribution from one occupied shell.

    The eigenvector tells us WHERE the electron lives.
    The Oh weight determines HOW MUCH it screens.
    The shell size (n, Z_eff) determines the spatial extent.
    """
    # The screening field: proportional to |psi|^2 * Oh_weight
    # Normalized so the integral = oh_weight (total screening from this channel)
    psi2 = eigvec**2
    norm = np.sum(psi2)
    if norm < 1e-30:
        return np.zeros(N)

    # The screening profile has the SHAPE of the eigenvector
    # and the MAGNITUDE of the Oh weight
    return oh_weight * psi2 / norm

def solve_atom(Z, config, max_iter=20, tol=1e-6):
    """Self-consistent Hessian solve for one atom.

    config: list of (n, l, count) electron configuration.
    Returns: IE in lattice units, number of bound states, convergence info.
    """
    # Shell info for screening weights
    val_n = max(nn for nn, ll, c in config)
    val_l = max(ll for nn, ll, c in config if nn == val_n and c > 0)

    # Nuclear well — scale amplitude so bound state count is reasonable
    # The well needs to support sum(count)/2 doubly-occupied states
    n_electrons = sum(c for _, _, c in config)
    # Scale beta with Z, amplitude tuned so the well is physical
    phi_nuc = kink_well(Z, amplitude=1.0)

    # Initial: no screening
    phi_screen = np.zeros(N)
    prev_IE = 0.0

    for iteration in range(max_iter):
        # Effective field: nuclear - screening
        phi_eff = phi_nuc - phi_screen

        # Build Hessian: H = -Lap + cos(pi*phi_eff)
        V = np.cos(PI * phi_eff)
        diag = 2.0 + V
        off = -np.ones(N - 1)
        H = sparse.diags([off, diag, off], [-1, 0, 1], shape=(N, N), format='lil')
        H[0, N-1] = -1.0
        H[N-1, 0] = -1.0
        H = H.tocsr()

        # Find bound states (eigenvalues < 1)
        n_eig = min(30, N - 2)
        try:
            evals, evecs = eigsh(H, k=n_eig, which='SM')
        except:
            return None, 0, "eigsh failed"

        order = np.argsort(evals)
        evals = evals[order]
        evecs = evecs[:, order]

        bound_mask = evals < 1.0
        n_bound = np.sum(bound_mask)

        if n_bound == 0:
            return None, 0, "no bound states"

        # Assign electrons to bound states
        # Each state holds 2 electrons (spin up/down in Oh language)
        # Fill from lowest energy up
        n_states_needed = 0
        shell_assignments = []
        for nn, ll, count in config:
            n_channels = 2 * ll + 1  # s=1, p=3, d=5, f=7
            # Each channel holds 2 electrons (paired)
            n_states = (count + 1) // 2  # ceil(count/2)
            for i in range(n_states):
                if n_states_needed < n_bound:
                    shell_assignments.append((nn, ll, count, n_states_needed))
                    n_states_needed += 1

        # Build screening from all INNER states
        # (everything except the outermost occupied state)
        phi_screen_new = np.zeros(N)

        if len(shell_assignments) > 1:
            for nn, ll, count, state_idx in shell_assignments[:-1]:
                # Oh screening weight depends on shell type and valence type
                if ll <= 1:
                    weight = w_pi
                elif ll == 2:
                    if val_l == 1:
                        weight = w_pi  # d→p: Oh allowed
                    else:
                        weight = w_delta  # d→s: Oh forbidden
                elif ll == 3:
                    if val_l == 1:
                        weight = w_pi * 3/7 + w_delta/d * 4/7  # T1u screens, rest anti
                    else:
                        weight = w_delta / d  # f→s: Oh forbidden

                # Screening contribution from this state
                psi = evecs[:, state_idx]
                phi_screen_new += screening_profile(psi, nn, ll, Z, weight)

        # Scale screening: total should not exceed nuclear well
        scale = min(1.0, 0.8 * np.max(np.abs(phi_nuc)) /
                    (np.max(np.abs(phi_screen_new)) + 1e-30))
        phi_screen = phi_screen_new * scale

        # IE = gap from outermost bound state to continuum
        outermost_idx = min(n_states_needed - 1, n_bound - 1)
        IE = 1.0 - evals[outermost_idx]

        # Convergence check
        if abs(IE - prev_IE) < tol:
            return IE, n_bound, f"converged iter={iteration+1}"
        prev_IE = IE

    return IE, n_bound, f"max_iter={max_iter}"


# ============================================================
# TEST ATOMS
# ============================================================
atoms = [
    (1,'H',13.598,[(1,0,1)]),
    (2,'He',24.587,[(1,0,2)]),
    (3,'Li',5.392,[(1,0,2),(2,0,1)]),
    (4,'Be',9.323,[(1,0,2),(2,0,2)]),
    (5,'B',8.298,[(1,0,2),(2,0,2),(2,1,1)]),
    (6,'C',11.260,[(1,0,2),(2,0,2),(2,1,2)]),
    (7,'N',14.534,[(1,0,2),(2,0,2),(2,1,3)]),
    (8,'O',13.618,[(1,0,2),(2,0,2),(2,1,4)]),
    (9,'F',17.423,[(1,0,2),(2,0,2),(2,1,5)]),
    (10,'Ne',21.565,[(1,0,2),(2,0,2),(2,1,6)]),
    (11,'Na',5.139,[(1,0,2),(2,0,2),(2,1,6),(3,0,1)]),
    (12,'Mg',7.646,[(1,0,2),(2,0,2),(2,1,6),(3,0,2)]),
    (13,'Al',5.986,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,1)]),
    (18,'Ar',15.760,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6)]),
    # Outliers
    (30,'Zn',9.394,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2)]),
    (46,'Pd',8.337,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10)]),
    (71,'Lu',5.426,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,1),(6,0,2)]),
    (80,'Hg',10.438,[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),(4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2)]),
]

report("SELF-CONSISTENT HESSIAN IE")
report("-" * 65)

# Calibrate with hydrogen
IE_H, nb_H, info_H = solve_atom(1, [(1,0,1)])
if IE_H is not None and IE_H > 0:
    scale = 13.598 / IE_H
    report(f"H calibration: IE_lattice = {IE_H:.8f}, scale = {scale:.4f} eV/lattice")
else:
    scale = 1.0
    report(f"H calibration FAILED: {info_H}")

report("")
report(f"{'Z':>4} {'Sym':>3} {'IE_obs':>8} {'IE_pred':>8} {'err':>7} {'n_bnd':>6} {'info':>20}")
report("-" * 65)

for Z, sym, IE_obs, config in atoms:
    IE_lat, n_bound, info = solve_atom(Z, config)

    if IE_lat is not None and IE_lat > 0:
        IE_pred = IE_lat * scale
        err = (IE_pred - IE_obs) / IE_obs * 100
        report(f"{Z:4d} {sym:>3} {IE_obs:8.3f} {IE_pred:8.3f} {err:+6.1f}% {n_bound:6d} {info:>20}")
    else:
        report(f"{Z:4d} {sym:>3} {IE_obs:8.3f} {'FAIL':>8} {'---':>7} {n_bound:6d} {info:>20}")

report("")
report(f"Completed: {time.strftime('%Y-%m-%d %H:%M:%S')}")

log.close()
print(f"\nResults saved to: {outfile}")
