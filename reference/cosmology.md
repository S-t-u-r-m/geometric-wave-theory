# Cosmological Parameters

*Part of the [GWT Complete Reference](../gwt_complete_reference.md). See also: [Forces](forces.md), [Coupling Constants](coupling_constants.md).*

| Parameter | Formula | Predicted | Observed | Error |
|-----------|---------|-----------|----------|-------|
| Dark energy fraction Omega_Lambda | (d-1)/d | 0.667 | 0.685 | -2.7% |
| Hubble constant H_0 (bare) | (c/l_P) * exp(-1/alpha) / d^3 | 64.5 km/s/Mpc | — | DERIVED |
| Hubble constant H_0 (CMB) | H_0(bare) * d^3/(d^3-1) | 67.0 km/s/Mpc | 67.4 | -0.6% [SPECULATIVE] |
| Hubble constant H_0 (local) | H_0(bare) * (d/(d-1))^(1/d) | 73.8 km/s/Mpc | 73.0 | +1.1% [SPECULATIVE] |
| Hubble tension ratio | H_0(local)/H_0(CMB) | 1.102 | 1.083 | +1.8% [SPECULATIVE] |
| Cosmic age t_0 | Friedmann + Omega_Lambda=2/3 | 13.66 Gyr | 13.8 Gyr | -1.0% |
| Cosmological constant Lambda | 2*H_0^2/c^2 | 1.048×10^-52 m^-2 | 1.088×10^-52 | -3.6% |
| Dark energy density u_DE | k*a/(8*R_H^2) | 5.05×10^-10 J/m^3 | 5.26×10^-10 | -4.0% |
| Deceleration parameter q_0 | -1/(d-1) | -0.500 | -0.55 | -9.1% |
| Dark energy EOS w | -1 exactly | -1.00 | -1 ± 0.1 | exact |
| Dark energy EOS w_a | 0 exactly | 0 | 0 ± 0.3 | exact |
| MOND acceleration a_0 | c*H_0/(pi*sqrt(d)) | 1.196×10^-10 m/s^2 | 1.2×10^-10 | -0.3% |
| G_eff in halos | (Omega_m/Omega_b)*G_N | 6.8 G_N | consistent | — |
| Baryon asymmetry eta_B | J × alpha^2 × d/2^d | 5.86×10^-10 | 6.1×10^-10 | -4.0% |

---

### Hubble constant — three values from one formula

The Hubble constant H_0 is the expansion rate of the lattice. The base formula is:
```
H_0(bare) = (c/l_P) * exp(-1/alpha) / d^3

  c/l_P = Planck frequency = 1.855 × 10^43 Hz
  exp(-1/alpha) = exp(-137.042) bridges ~60 orders of magnitude
  d^3 = 27: total orientations on the d-cube

  H_0(bare) = 64.5 km/s/Mpc
```

This is the pure lattice expansion rate. No observer sees this directly because
all observers are made of matter (kinks = tori) that interacts with the lattice.

**TWO independent corrections** give the two measured values:

#### H_0(CMB) = 67.0 km/s/Mpc — the torus correction [SPECULATIVE, -0.6%]
```
H_0(CMB) = H_0(bare) * d^3 / (d^3 - 1) = 64.5 * 27/26 = 67.0 km/s/Mpc

  Observed (Planck 2018): 67.4 ± 0.5 km/s/Mpc. Error: -0.6%.
```
The CMB measures expansion as seen by matter. Matter IS kinks (toroidal vortices).
A kink wraps one direction of the d-cube, locking that direction out of the
expansion. Of d^3 = 27 total orientations, only d^3 - 1 = 26 are free to expand.

This is the SAME 26/27 = (d^3-1)/d^3 torus correction as the proton charge radius:
```
  r_p = (d+1) * (1 - alpha*(d^3-1)/(d^3*pi^2)) * hbar*c/m_p
```
The proton locks one direction (kink wrapping). The expansion locks one direction
(matter frame). Same geometry, different physics.

#### H_0(local) = 73.8 km/s/Mpc — the MOND bias [SPECULATIVE, +1.1%]
```
H_0(local) = H_0(bare) * (d/(d-1))^(1/d) = 64.5 * (3/2)^(1/3) = 73.8 km/s/Mpc

  Observed (SH0ES 2022): 73.0 ± 1.0 km/s/Mpc. Error: +1.1%.
```
The local distance ladder (Cepheids + Type Ia SNe) measures expansion through
nearby galaxies. These galaxies have outer regions in the MOND regime (a < a_0),
where the effective gravitational coupling changes by the generation factor.

The bias factor is the generation factor d/(d-1) = 3/2 taken to the 1/d = 1/3
power — the bulk (3D) projection. This is the SAME d/(d-1) that gives the
muon-electron mass ratio, the Koide formula, the strange-light quark ratio,
and the kaon mass. Taken to the 1/d power because the measurement traverses
all d spatial dimensions (unlike the 1/(d-1) surface power used for CKM).

**Physical mechanism:** The Cepheid period-luminosity calibration passes through
the a_0 transition zone. At a < a_0, the effective gravity strengthens
(MOND regime), modifying the stellar oscillation period. The distance ladder
then propagates this calibration bias into the inferred H_0.

#### The Hubble tension — SPECULATIVE RESOLUTION
```
H_0(local) / H_0(CMB) = [(d/(d-1))^(1/d)] / [d^3/(d^3-1)]
                       = 1.1447 / 1.0385
                       = 1.102

Observed: 73.0 / 67.4 = 1.083. Error: +1.8%.
```
There is no tension. There are two CORRECT values of H_0, each with a different
geometric correction:
  - CMB measures through matter geometry (torus correction, 27/26)
  - Local measures through MOND transition (generation factor, (3/2)^(1/3))
Both corrections are derived from d=3. The 8% gap between them is not a
discrepancy — it is a prediction.

**Falsifiable:** If future measurements converge to a single H_0 (eliminating
the tension through systematic error reduction), GWT's two-value prediction
would be falsified. The prediction is that the tension is REAL and has the
specific ratio (d/(d-1))^(1/d) × (d^3-1)/d^3 = 1.102.

---

### Cosmic age
```
t_0 = (2/(3*H_0)) * arcsinh(sqrt(Omega_Lambda/(1-Omega_Lambda))) / sqrt(Omega_Lambda)

With Omega_Lambda = 2/3 and H_0(CMB) = 67.0:
  Omega_Lambda/(1-Omega_Lambda) = 2, arcsinh(sqrt(2)) = ln(sqrt(2)+sqrt(3)) = 1.1462
  sqrt(Omega_Lambda) = sqrt(2/3) = 0.8165

t_0 = (2/(3*H_0)) * 1.1462/0.8165 = (2/(3*H_0)) * 1.4037
t_0 = 13.66 Gyr   (obs: 13.8 Gyr, -1.0%)
```

### Cosmological constant
```
Lambda = 2*H_0^2 / c^2 = 1.048 × 10^-52 m^-2   (obs: 1.088e-52, -3.6%)
```
Follows directly from Omega_Lambda = 2/3 and the Friedmann equation.

### Dark energy equation of state
```
w = -1 exactly
w_a = 0 exactly
```
The lattice's L3 wave period (the largest standing wave) is far longer than the age of the universe. On cosmological timescales, this acts as a constant boundary pressure — giving w = -1 exactly, not approximately. Dark energy is the lattice's transverse restoring force, not a dynamical field. No time evolution → w_a = 0.

### Dark energy density
```
u_DE = k*a / (8*R_H^2)

where k*a = (2/pi)*l_P * (c^4/(pi*G)) and R_H = c/H_0.

Algebraic proof that Omega_Lambda = 2/3 exactly:
  u_DE = c^2*H_0^2 / (4*pi*G)
  u_crit = 3*c^2*H_0^2 / (8*pi*G)
  Omega_Lambda = u_DE/u_crit = 8/(4*3) = 2/3

H_0 cancels completely. The result is purely geometric.
```

### MOND acceleration
```
a_0 = c*H_0 / (pi*sqrt(d))
    = 1.196 × 10^-10 m/s^2   (obs: 1.2e-10, -0.3%)
```
Crossover between local gravitational wave gradient and cosmic carrier wave gradient. At a < a_0, the two gradients interfere, producing MOND behavior and flat rotation curves via v^4 = G_N*M*a_0 (Baryonic Tully-Fisher).

### Baryon asymmetry derivation
```
eta_B = J × alpha^2 × (d/2^d)
      = 2.93e-5 × 5.32e-5 × 0.375
      = 5.86e-10   (obs: 6.1e-10, -4.0%)
```

**Step-by-step chain:**

1. **Kink tunneling = baryon number violation.**
   The sine-Gordon kink phi(x) = (4/pi)arctan(exp(x)) connects adjacent
   potential minima. Topological charge Q = 1 = baryon number.
   Tunneling changes winding number -> changes baryon number.

2. **CP asymmetry from interference.**
   Tree-level tunneling is CP-symmetric (V(-phi) = V(phi)).
   Asymmetry requires interference between tree and loop amplitudes:
   |A(B)|^2 - |A(B~)|^2 = 4 Re(A_tree) Im(A_loop)

3. **J = Jarlskog invariant** (from CKM, fully derived from bare mass ratios).
   Measures the "area" of the CP-violation triangle. All CKM angles from
   quark mass ratios on the (d-1)-dimensional proton surface.

4. **alpha^2 = wave overlap probability.**
   - alpha^1: CP-violating loop amplitude (one loop with CKM phase)
   - alpha^1: coupling of that loop to the kink tunneling process
   - Total: alpha^2 = time-averaged probability of both waves overlapping
   - Cross-check: alpha^1 overshoots by 132x, alpha^3 undershoots by 143x.
     alpha^2 is the UNIQUE power that works.

5. **d/2^d = 3/8 = lattice projection factor.**
   Kink tunneling is 1D along one axis (d choices) acting on d-cube
   (2^d = 8 cells per vertex). Fraction participating = d/2^d.

6. **Coefficient = 1 from lattice quantization.**
   The lattice discretizes what QFT integrates: one tunneling event
   per cell per Hubble time. No continuous rate integral needed.

**Error budget:** The -4% comes from J being ~5% low (V_ub = sqrt(m_u/m_t)
is sensitive to GWT quark mass errors). Using PDG J gives eta_B = 6.15e-10 (+0.8%).

See: calculations/cosmology/kink_phase_baryogenesis.py for full derivation with cross-checks.
