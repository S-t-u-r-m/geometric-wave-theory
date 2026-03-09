"""
WAVE MODE COUPLING v2: Testing harmonic decay laws
====================================================
The coupling between wave modes at different n should follow
a harmonic decay law. For period 2, (n_i/n_v)^2 works perfectly.

But maybe the coupling depends on MORE than just the frequency ratio.
In wave physics, the overlap integral between two modes depends on:
1. Frequency ratio (n_i/n_v)
2. Number of nodes between them (|n_v - n_i|)
3. Whether modes share angular patterns

KEY QUESTION: Is the coupling (n_i/n_v)^2, or does it have a
node-dependent factor?

Also test: does the coupling depend on the TOTAL number of
inner modes (many-body effect)?
"""
import numpy as np

pi = np.pi
E_H = 13.6057
dd = 3

atoms = {
    'H':  {'Z': 1,  'n': 1, 'l': 0, 'Z_eff': 1.0000},
    'Li': {'Z': 3,  'n': 2, 'l': 0, 'Z_eff': 1.2792},
    'B':  {'Z': 5,  'n': 2, 'l': 1, 'Z_eff': 2.4214},
    'C':  {'Z': 6,  'n': 2, 'l': 1, 'Z_eff': 3.1358},
    'N':  {'Z': 7,  'n': 2, 'l': 1, 'Z_eff': 3.8340},
    'O':  {'Z': 8,  'n': 2, 'l': 1, 'Z_eff': 4.4532},
    'F':  {'Z': 9,  'n': 2, 'l': 1, 'Z_eff': 5.0998},
    'Na': {'Z': 11, 'n': 3, 'l': 0, 'Z_eff': 2.5074},
    'Cl': {'Z': 17, 'n': 3, 'l': 1, 'Z_eff': 6.1161},
}

# Inner mode configs: (n_i, l_i, count)
configs = {
    'H':  [],
    'Li': [(1, 0, 2)],
    'B':  [(1, 0, 2), (2, 0, 2)],
    'C':  [(1, 0, 2), (2, 0, 2), (2, 1, 1)],
    'N':  [(1, 0, 2), (2, 0, 2), (2, 1, 2)],
    'O':  [(1, 0, 2), (2, 0, 2), (2, 1, 3)],
    'F':  [(1, 0, 2), (2, 0, 2), (2, 1, 4)],
    'Na': [(1, 0, 2), (2, 0, 2), (2, 1, 6)],
    'Cl': [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 4)],
}

# =============================================================================
# EXACT SCREENING VALUES NEEDED
# =============================================================================
print("=" * 80)
print("  EXACT COUPLING VALUES IMPLIED BY DATA")
print("=" * 80)
print()
print("  Z_eff = Z - sum(s_i)")
print("  s_i = 1 - coupling_i")
print("  coupling_i = g * (n_i/n_v)^p")
print()

# From each atom, extract what total coupling is needed
print(f"{'Atom':>4} {'Z':>3} {'n_v':>3} {'N_in':>4} {'Z_eff':>7} {'S':>7} "
      f"{'coupling_total':>14} {'coupling/mode':>13}")
print("-" * 70)

for name in atoms:
    info = atoms[name]
    Z = info['Z']; n_v = info['n']; Ze = info['Z_eff']
    N_inner = sum(cnt for _, _, cnt in configs[name])
    S = Z - Ze
    coupling_total = N_inner - S  # sum of all g_i * (n_i/n_v)^2
    coupling_per = coupling_total / N_inner if N_inner > 0 else 0
    print(f"{name:>4} {Z:3d} {n_v:3d} {N_inner:4d} {Ze:7.4f} {S:7.4f} "
          f"{coupling_total:14.4f} {coupling_per:13.4f}")


# =============================================================================
# WHAT (n_i/n_v)^p RATIOS ARE AVAILABLE?
# =============================================================================
print()
print("=" * 80)
print("  HARMONIC RATIOS at each shell gap")
print("=" * 80)

# For period 2 atoms (n_v=2): only ratio (1/2)^p
# For period 3 atoms (n_v=3): ratios (1/3)^p and (2/3)^p and (3/3)^p=1
# The model has to work with these specific ratios

for p in [1, 2, 3, 4]:
    print(f"\n  p={p}:")
    for n_v in [2, 3]:
        for n_i in range(1, n_v+1):
            ratio = (n_i/n_v)**p
            print(f"    (n_i={n_i}, n_v={n_v}): ({n_i}/{n_v})^{p} = {ratio:.6f}")


# =============================================================================
# TEST: Single g, variable power p
# =============================================================================
print()
print("=" * 80)
print("  SCAN: s = 1 - g * (n_i/n_v)^p  (single g)")
print("=" * 80)

# For same subshell: g -> g_same; for different: g -> g_diff
# But let's first see if a SINGLE g with different p works

best = (999, 0, 0, 0)
for p in np.arange(1.0, 6.1, 0.1):
    for gs in np.arange(0.2, 1.0, 0.01):
        for gd in np.arange(0.2, 1.0, 0.01):
            errs = []
            for name in atoms:
                info = atoms[name]
                Z = info['Z']; n_v = info['n']; l_v = info['l']
                S = 0
                for (n_i, l_i, cnt) in configs[name]:
                    if n_i == n_v and l_i == l_v:
                        g = gs
                    else:
                        g = gd
                    s = 1 - g * (n_i/n_v)**p
                    S += cnt * s
                Ze_pred = Z - S
                errs.append((Ze_pred - info['Z_eff'])**2)
            rms = np.sqrt(np.mean(errs))
            if rms < best[0]:
                best = (rms, p, gs, gd)

print(f"\n  Best: p={best[1]:.1f}, g_same={best[2]:.4f}, g_diff={best[3]:.4f}, RMS={best[0]:.4f}")

# Check GWT constant matches
p_b = best[1]; gs_b = best[2]; gd_b = best[3]
print(f"\n  g_same = {gs_b:.4f}")
gwt_vals = {
    '1/d': 1/dd, '2/d': 2/dd, '1/(d-1)': 1/(dd-1), '(d-1)/d': (dd-1)/dd,
    '2/(2d+1)': 2/(2*dd+1), '4/(2d+1)': 4/(2*dd+1), '1/(d+1)': 1/(dd+1),
    '2/(d+1)': 2/(dd+1), 'd/(d+1)': dd/(dd+1), '1/pi': 1/pi, '2/pi': 2/pi,
    '(d-1)/(d+1)': (dd-1)/(dd+1), 'pi/d^2': pi/dd**2,
}
for label, val in gwt_vals.items():
    if abs(val - gs_b) < 0.03:
        print(f"    ~ {label} = {val:.4f}")

print(f"\n  g_diff = {gd_b:.4f}")
for label, val in gwt_vals.items():
    if abs(val - gd_b) < 0.03:
        print(f"    ~ {label} = {val:.4f}")

# Predictions
print(f"\n  {'Atom':>4} {'Z':>3} {'Z_pred':>7} {'Z_real':>7} {'err':>7}")
print("  " + "-" * 35)
for name in atoms:
    info = atoms[name]
    Z = info['Z']; n_v = info['n']; l_v = info['l']
    S = 0
    for (n_i, l_i, cnt) in configs[name]:
        if n_i == n_v and l_i == l_v:
            g = gs_b
        else:
            g = gd_b
        s = 1 - g * (n_i/n_v)**p_b
        S += cnt * s
    Ze_pred = Z - S
    print(f"  {name:>4} {Z:3d} {Ze_pred:7.4f} {info['Z_eff']:7.4f} {Ze_pred-info['Z_eff']:+7.4f}")


# =============================================================================
# TEST: Node-dependent coupling
# =============================================================================
print()
print("=" * 80)
print("  NODE-DEPENDENT COUPLING: g * (n_i/n_v)^(2*|n_v-n_i|)")
print("=" * 80)

# In wave physics, modes separated by more nodes have weaker overlap
# The number of radial nodes between modes n_i and n_v is |n_v - n_i|
# Maybe the coupling power depends on this

print("\n  Coupling = g * (n_i/n_v)^(2*delta_n) where delta_n = n_v - n_i")
print("  For same shell (delta_n=0): coupling = g * 1 = g")
print("  For adjacent (delta_n=1): coupling = g * (n_i/n_v)^2")
print("  For gap of 2 (delta_n=2): coupling = g * (n_i/n_v)^4")

best_node = (999, 0, 0)
for gs in np.arange(0.2, 1.0, 0.01):
    for gd in np.arange(0.2, 1.0, 0.01):
        errs = []
        for name in atoms:
            info = atoms[name]
            Z = info['Z']; n_v = info['n']; l_v = info['l']
            S = 0
            for (n_i, l_i, cnt) in configs[name]:
                delta_n = n_v - n_i
                if n_i == n_v and l_i == l_v:
                    g = gs
                else:
                    g = gd
                p = 2 * max(delta_n, 1)  # at least p=2 for same-n diff-l
                coupling = g * (n_i/n_v)**p
                S += cnt * (1 - coupling)
            Ze_pred = Z - S
            errs.append((Ze_pred - info['Z_eff'])**2)
        rms = np.sqrt(np.mean(errs))
        if rms < best_node[0]:
            best_node = (rms, gs, gd)

gs_n = best_node[1]; gd_n = best_node[2]
print(f"\n  Best: g_same={gs_n:.4f}, g_diff={gd_n:.4f}, RMS={best_node[0]:.4f}")

print(f"\n  g_same = {gs_n:.4f}")
for label, val in gwt_vals.items():
    if abs(val - gs_n) < 0.03:
        print(f"    ~ {label} = {val:.4f}")

print(f"\n  g_diff = {gd_n:.4f}")
for label, val in gwt_vals.items():
    if abs(val - gd_n) < 0.03:
        print(f"    ~ {label} = {val:.4f}")

print(f"\n  {'Atom':>4} {'Z':>3} {'Z_pred':>7} {'Z_real':>7} {'err':>7}")
print("  " + "-" * 35)
for name in atoms:
    info = atoms[name]
    Z = info['Z']; n_v = info['n']; l_v = info['l']
    S = 0
    for (n_i, l_i, cnt) in configs[name]:
        delta_n = n_v - n_i
        if n_i == n_v and l_i == l_v:
            g = gs_n
        else:
            g = gd_n
        p = 2 * max(delta_n, 1)
        coupling = g * (n_i/n_v)**p
        S += cnt * (1 - coupling)
    Ze_pred = Z - S
    print(f"  {name:>4} {Z:3d} {Ze_pred:7.4f} {info['Z_eff']:7.4f} {Ze_pred-info['Z_eff']:+7.4f}")


# =============================================================================
# TEST: Coupling as sin function (truly harmonic)
# =============================================================================
print()
print("=" * 80)
print("  HARMONIC WAVE COUPLING: g * |sin(pi * n_i / n_v)|^p")
print("=" * 80)

# In wave physics, the overlap integral between sin(n_i*x) and sin(n_v*x)
# involves sin(pi * n_i / n_v) when modes have different boundary conditions
# This is the TRUE harmonic coupling

print("\n  sin(pi * n_i / n_v) values:")
for n_v in [2, 3]:
    for n_i in range(1, n_v+1):
        val = abs(np.sin(pi * n_i / n_v))
        print(f"    n_i={n_i}, n_v={n_v}: |sin(pi*{n_i}/{n_v})| = {val:.6f}")

# sin(pi*1/2) = 1.000
# sin(pi*1/3) = 0.866
# sin(pi*2/3) = 0.866
# sin(pi*3/3) = 0.000 -- same frequency = zero coupling!?

# That's wrong for same-shell. sin(pi*n/n) = sin(pi) = 0 means zero coupling.
# But same-shell electrons DO screen.

# Try: cos instead? cos(pi*n_i/n_v)?
print("\n  cos(pi * n_i / (n_v+1)) values:")
for n_v in [1, 2, 3]:
    for n_i in range(1, n_v+1):
        val = abs(np.cos(pi * n_i / (n_v+1)))
        print(f"    n_i={n_i}, n_v={n_v}: |cos(pi*{n_i}/{n_v+1})| = {val:.6f}")

# cos(pi*1/2) = 0  (1s screening 1s = 0, not useful for H)
# cos(pi*1/3) = 0.500 = 1/2
# cos(pi*2/3) = 0.500 = 1/2 (same!)
# cos(pi*1/4) = 0.707 = 1/sqrt(2)
# cos(pi*2/4) = 0.000
# cos(pi*3/4) = 0.707

# What about sin(pi * n_i / (2*n_v))?
print("\n  sin(pi * n_i / (2*n_v)) values:")
for n_v in [1, 2, 3]:
    for n_i in range(1, n_v+1):
        val = np.sin(pi * n_i / (2*n_v))
        print(f"    n_i={n_i}, n_v={n_v}: sin(pi*{n_i}/{2*n_v}) = {val:.6f}")

# sin(pi*1/2) = 1.000 (1s: full coupling to itself)
# sin(pi*1/4) = 0.707 = 1/sqrt(2) (1s coupling to n=2)
# sin(pi*2/4) = 1.000 (n=2 coupling to n=2 = same shell)
# sin(pi*1/6) = 0.500 (1s coupling to n=3)
# sin(pi*2/6) = 0.866 = sqrt(3)/2 (n=2 coupling to n=3)
# sin(pi*3/6) = 1.000 (n=3 coupling to n=3 = same shell)

# This is beautiful! sin(pi*n_i/(2*n_v)) gives:
# - 1.0 for same shell (maximum coupling)
# - Decreasing for deeper shells
# - The rate depends on the ratio

print("\n  Testing: coupling = g * sin^2(pi * n_i / (2*n_v))")
best_sin = (999, 0, 0)
for gs in np.arange(0.2, 1.0, 0.01):
    for gd in np.arange(0.2, 1.0, 0.01):
        errs = []
        for name in atoms:
            info = atoms[name]
            Z = info['Z']; n_v = info['n']; l_v = info['l']
            S = 0
            for (n_i, l_i, cnt) in configs[name]:
                harm = np.sin(pi * n_i / (2*n_v))**2
                if n_i == n_v and l_i == l_v:
                    g = gs
                else:
                    g = gd
                coupling = g * harm
                S += cnt * (1 - coupling)
            Ze_pred = Z - S
            errs.append((Ze_pred - info['Z_eff'])**2)
        rms = np.sqrt(np.mean(errs))
        if rms < best_sin[0]:
            best_sin = (rms, gs, gd)

gs_s = best_sin[1]; gd_s = best_sin[2]
print(f"\n  Best: g_same={gs_s:.4f}, g_diff={gd_s:.4f}, RMS={best_sin[0]:.4f}")

# Note: sin^2(pi*n_i/(2*n_v)) for n_i=n_v = sin^2(pi/2) = 1
# For n_i=1, n_v=2: sin^2(pi/4) = 0.5 = 1/2
# For n_i=1, n_v=3: sin^2(pi/6) = 0.25 = 1/4
# For n_i=2, n_v=3: sin^2(pi/3) = 0.75 = 3/4

# Compare: (n_i/n_v)^2 gives:
# n_i=1, n_v=2: 1/4
# n_i=1, n_v=3: 1/9
# n_i=2, n_v=3: 4/9

# sin^2(pi*n/(2N)) decays SLOWER than (n/N)^2
# This might be what Na needs: more coupling back from inner shells

print(f"\n  {'Atom':>4} {'Z':>3} {'Z_pred':>7} {'Z_real':>7} {'err':>7}")
print("  " + "-" * 35)
for name in atoms:
    info = atoms[name]
    Z = info['Z']; n_v = info['n']; l_v = info['l']
    S = 0
    for (n_i, l_i, cnt) in configs[name]:
        harm = np.sin(pi * n_i / (2*n_v))**2
        if n_i == n_v and l_i == l_v:
            g = gs_s
        else:
            g = gd_s
        coupling = g * harm
        S += cnt * (1 - coupling)
    Ze_pred = Z - S
    print(f"  {name:>4} {Z:3d} {Ze_pred:7.4f} {info['Z_eff']:7.4f} {Ze_pred-info['Z_eff']:+7.4f}")

print(f"\n  g_same = {gs_s:.4f}")
for label, val in gwt_vals.items():
    if abs(val - gs_s) < 0.03:
        print(f"    ~ {label} = {val:.4f}")
print(f"\n  g_diff = {gd_s:.4f}")
for label, val in gwt_vals.items():
    if abs(val - gd_s) < 0.03:
        print(f"    ~ {label} = {val:.4f}")


# =============================================================================
# COMPARISON TABLE: All models
# =============================================================================
print()
print("=" * 80)
print("  COMPARISON: All coupling models")
print("=" * 80)

models = {
    '(n_i/n_v)^2, 2g': lambda ni, nv, same: (2/dd if same else 4/(2*dd+1)) * (ni/nv)**2,
    'sin^2(pi*n/(2N)), 2g': lambda ni, nv, same: (gs_s if same else gd_s) * np.sin(pi*ni/(2*nv))**2,
    '(n_i/n_v)^p_best, 2g': lambda ni, nv, same: (gs_b if same else gd_b) * (ni/nv)**p_b,
    'node-dep 2*dn, 2g': lambda ni, nv, same: (gs_n if same else gd_n) * (ni/nv)**(2*max(nv-ni,1)),
}

print(f"\n{'Model':>30} {'RMS_Ze':>8}  H    Li    B     C     N     O     F     Na    Cl")
print("-" * 110)

for mname, coupling_func in models.items():
    preds = {}
    for name in atoms:
        info = atoms[name]
        Z = info['Z']; n_v = info['n']; l_v = info['l']
        S = 0
        for (n_i, l_i, cnt) in configs[name]:
            same = (n_i == n_v and l_i == l_v)
            c = coupling_func(n_i, n_v, same)
            S += cnt * (1 - c)
        preds[name] = Z - S

    errs = [(preds[n] - atoms[n]['Z_eff'])**2 for n in atoms]
    rms = np.sqrt(np.mean(errs))
    vals = " ".join(f"{preds[n]-atoms[n]['Z_eff']:+5.2f}" for n in atoms)
    print(f"{mname:>30} {rms:8.4f}  {vals}")
