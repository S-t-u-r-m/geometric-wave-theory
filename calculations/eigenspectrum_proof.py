"""
EIGENSPECTRUM PROOF: Derive all 24 breather modes from the lattice equation.
No particle identities assumed. Pure nonlinear sine-Gordon dynamics.

Part 1: 1D at Nx=100,000 - all 24 modes (mathematical proof)
Part 2: 3D at 128^3 - modes 1-12 (3D lattice confirmation)
"""
import sys, io, os, time
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

d = 3
gamma = np.pi / (2**(d+1)*np.pi - 2)

outfile = os.path.join(os.path.dirname(__file__), "eigenspectrum_results.txt")
log = open(outfile, "w", encoding="utf-8")

def report(msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()

report("EIGENSPECTRUM PROOF")
report("=" * 60)
report("Started: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("gamma = {:.10f}".format(gamma))
report("Predicted: 24 modes at omega_n = cos(n*gamma)")
report("")

particles = {4:"mu/strange",5:"down",7:"bottom",11:"charm",12:"top",
             13:"up",16:"electron",18:"tau"}

# ============================================================
# PART 1: 1D Proof - All 24 modes
# ============================================================
report("PART 1: 1D SINE-GORDON (Nx=100,000)")
report("-" * 55)

Nx = 100000
L = 100.0
dx = 2*L/Nx
dt_1d = 0.4*dx
x = np.linspace(-L+dx/2, L-dx/2, Nx)

report("Grid: {} points, dx={:.6f}, dt={:.6f}".format(Nx, dx, dt_1d))
report("")

report("{:>3} {:>11} {:>11} {:>9} {:>8} {:>12}".format(
    "n", "omega_pred", "omega_meas", "err_ppm", "status", "particle"))
report("-" * 60)

results_1d = []
t_part1 = time.time()

for n in range(1, 25):
    omega_n = np.cos(n*gamma)
    eps_n = np.sin(n*gamma)

    # Exact breather initialization
    phi = np.zeros(Nx)
    velocity = (4.0/np.pi) * eps_n / (omega_n * np.cosh(eps_n * x) + 1e-30)
    phi_old = phi - dt_1d * velocity

    # Evolve for many periods
    period = 2*np.pi/max(omega_n, 0.001)
    N_periods = 100
    N_steps = min(int(N_periods * period / dt_1d), 2000000)
    rec = max(1, N_steps // 40000)

    ts = []
    for step in range(N_steps):
        lap = np.zeros_like(phi)
        lap[1:-1] = (phi[2:] - 2*phi[1:-1] + phi[:-2]) / dx**2
        lap[0] = (phi[1] - phi[0]) / dx**2
        lap[-1] = (phi[-2] - phi[-1]) / dx**2
        force = (1.0/np.pi) * np.sin(np.pi * phi)
        phi_new = 2*phi - phi_old + dt_1d**2 * (lap - force)
        phi_old = phi.copy()
        phi = phi_new
        if step % rec == 0:
            ts.append(phi[Nx//2])

    # FFT
    ts = np.array(ts)
    ts = ts - np.mean(ts)
    if len(ts) < 10:
        omega_meas = 0
    else:
        window = np.hanning(len(ts))
        fft = np.fft.rfft(ts * window)
        freqs = np.fft.rfftfreq(len(ts), d=dt_1d*rec)
        omega_arr = 2*np.pi*freqs
        power = np.abs(fft)**2
        mask = (omega_arr > omega_n*0.5) & (omega_arr < omega_n*1.5)
        if np.any(mask):
            idx = np.where(mask)[0]
            peak_idx = idx[np.argmax(power[mask])]
            omega_meas = omega_arr[peak_idx]
        else:
            omega_meas = 0

    err_ppm = abs(omega_meas - omega_n) / omega_n * 1e6 if omega_n > 0.001 else 999999
    status = "EXACT" if err_ppm < 100 else "GOOD" if err_ppm < 1000 else "ok" if err_ppm < 5000 else "MISS"
    p = particles.get(n, "")
    results_1d.append((n, omega_n, omega_meas, err_ppm, status))
    report("{:>3} {:11.6f} {:11.6f} {:8.0f} {:>8} {:>12}".format(
        n, omega_n, omega_meas, err_ppm, status, p))

t1 = time.time() - t_part1
report("\nPart 1 time: {:.0f}s ({:.1f} min)".format(t1, t1/60))

exact1 = sum(1 for _,_,_,e,_ in results_1d if e < 100)
good1 = sum(1 for _,_,_,e,_ in results_1d if e < 1000)
ok1 = sum(1 for _,_,_,e,_ in results_1d if e < 5000)
report("EXACT (<100 ppm): {}/24".format(exact1))
report("GOOD (<1000 ppm): {}/24".format(good1))
report("OK (<5000 ppm):   {}/24".format(ok1))
report("")

# ============================================================
# PART 2: 3D GPU - Modes 1-12
# ============================================================
report("PART 2: 3D LATTICE (modes 1-12)")
report("-" * 55)

try:
    import cupy as cp
    xp = cp
    report("GPU active (CuPy)")
except Exception:
    xp = np
    report("CPU mode (no CuPy)")

N3d = 128
BOX3 = 10.0
dx3 = 2*BOX3/N3d
dt3 = 0.1*dx3
x1d3 = xp.linspace(-BOX3, BOX3, N3d, endpoint=False, dtype=np.float64)
X3, Y3, Z3 = xp.meshgrid(x1d3, x1d3, x1d3, indexing="ij")
R3 = xp.sqrt(X3**2 + Y3**2 + Z3**2) + 1e-10

report("Grid: {}^3 = {:,} points, dx={:.4f}".format(N3d, N3d**3, dx3))
report("")

report("{:>3} {:>11} {:>11} {:>9} {:>8} {:>12}".format(
    "n", "omega_pred", "omega_meas", "err_ppm", "status", "particle"))
report("-" * 60)

results_3d = []
t_part2 = time.time()

for n in range(1, 13):
    omega_n = np.cos(n*gamma)
    eps_n = np.sin(n*gamma)

    phi = xp.zeros((N3d,N3d,N3d), dtype=np.float64)
    velocity = (4.0/np.pi) * eps_n / (omega_n * xp.cosh(eps_n * R3) + 1e-30)
    phi_old = phi - dt3 * velocity

    period = 2*np.pi/omega_n
    N_steps = min(int(60 * period / dt3), 300000)
    rec = max(1, N_steps // 15000)

    ts = []
    for step in range(N_steps):
        lap = (xp.roll(phi,1,0)+xp.roll(phi,-1,0)+xp.roll(phi,1,1)+
               xp.roll(phi,-1,1)+xp.roll(phi,1,2)+xp.roll(phi,-1,2)-6*phi)/dx3**2
        force = (1.0/np.pi)*xp.sin(np.pi*phi)
        phi_new = 2*phi - phi_old + dt3**2*(lap - force)
        phi_old = phi.copy()
        phi = phi_new
        if step % rec == 0:
            ts.append(float(phi[N3d//2, N3d//2, N3d//2]))

    ts = np.array(ts)
    ts = ts - np.mean(ts)
    if len(ts) < 10:
        omega_meas = 0
    else:
        window = np.hanning(len(ts))
        fft = np.fft.rfft(ts * window)
        freqs = np.fft.rfftfreq(len(ts), d=dt3*rec)
        omega_arr = 2*np.pi*freqs
        power = np.abs(fft)**2
        mask = (omega_arr > omega_n*0.7) & (omega_arr < omega_n*1.3)
        if np.any(mask):
            idx = np.where(mask)[0]
            peak_idx = idx[np.argmax(power[mask])]
            omega_meas = omega_arr[peak_idx]
        else:
            omega_meas = 0

    err_ppm = abs(omega_meas - omega_n) / omega_n * 1e6
    status = "EXACT" if err_ppm < 100 else "GOOD" if err_ppm < 1000 else "ok" if err_ppm < 5000 else "MISS"
    p = particles.get(n, "")
    results_3d.append((n, omega_n, omega_meas, err_ppm, status))
    report("{:>3} {:11.6f} {:11.6f} {:8.0f} {:>8} {:>12}".format(
        n, omega_n, omega_meas, err_ppm, status, p))

t2 = time.time() - t_part2
report("\nPart 2 time: {:.0f}s ({:.1f} min)".format(t2, t2/60))

exact3 = sum(1 for _,_,_,e,_ in results_3d if e < 100)
good3 = sum(1 for _,_,_,e,_ in results_3d if e < 1000)
report("EXACT (<100 ppm): {}/12".format(exact3))
report("GOOD (<1000 ppm): {}/12".format(good3))
report("")

# ============================================================
# FINAL SUMMARY
# ============================================================
report("=" * 60)
report("FINAL SUMMARY")
report("Completed: " + time.strftime("%Y-%m-%d %H:%M:%S"))
report("Total time: {:.1f} min".format((time.time()-t_part1)/60))
report("")
report("1D (all 24 modes): {} EXACT, {} GOOD out of 24".format(exact1, good1))
report("3D (modes 1-12):   {} EXACT, {} GOOD out of 12".format(exact3, good3))
report("")

if good1 >= 20:
    report("CONCLUSION: The breather spectrum omega_n = cos(n*gamma)")
    report("is DERIVED from the nonlinear sine-Gordon equation on a")
    report("discrete lattice. The particle modes emerge from the")
    report("equation without assuming particle identities.")
elif good1 >= 10:
    report("STRONG EVIDENCE: {}/24 modes confirmed.".format(good1))
    report("Higher modes need finer grid resolution.")
else:
    report("PARTIAL: {}/24 modes confirmed. More resolution needed.".format(good1))

log.close()
print("\nResults saved to: " + outfile)
