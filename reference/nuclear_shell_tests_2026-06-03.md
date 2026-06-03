# Nuclear shell-structure tests of "shells = breathers" hypothesis

*Tests run 2026-06-03 (morning) to evaluate the conceptual principle from
session_2026-06-03 that nuclear shells emerge from breather configurations
on torus structures.*

**Status framework**: derived / forced-but-unproven / fitted / open / candidate

---

## Test 1: Neutron-excess curve N(Z) for beta-stable nuclei

### Hypothesis
If nuclei are torus-like structures filled by breather configurations,
larger nuclei need more uncharged "filler" breathers (neutrons). The
neutron-excess curve should fall out naturally.

### Method
- Compiled 36 beta-stable nuclei (Z from 1 to 92)
- Computed (N-Z)/A vs A^(2/3) (standard nuclear physics framing)
- Fitted linear relationship; compared coefficient to standard SEMF

### Result
```
(N-Z)/A = -0.027 + 0.00694 * A^(2/3)
```

| Source | b coefficient |
|--------|---------------|
| GWT framework (this fit) | 0.00694 |
| Standard SEMF | 0.00772 |
| Agreement | 10% off |

### Verdict: **CONSISTENT, not derived**

- ✓ Framework is consistent with observed neutron-excess curve
- ✓ The "torus needs neutral fillers" picture matches qualitative trend
- ✗ We have NOT derived the coefficient (0.00694 or 0.00772) from GWT primitives
- ✗ Standard SEMF already explains this via Coulomb + asymmetry balance

For this to be a real WIN over standard physics, GWT would need to derive:
- a_C (Coulomb constant) = 0.71 MeV from EM coupling between charged tori
- a_A (asymmetry constant) = 23 MeV from breather-pair imbalance energy
- Their ratio → b = 0.00772

This is real nuclear physics computation, not session work.

---

## Test 2: Magic numbers from harmonic-oscillator + spin-orbit

### Observed magic numbers
```
N or Z = 2, 8, 20, 28, 50, 82, 126
```

### Shell capacities (differences)
```
2, 6, 12, 8, 22, 32, 44
```

### 3D harmonic oscillator (no spin-orbit)
```
Cumulative (HO magic): 2, 8, 20, 40, 70, 112, 168
Observed:               2, 8, 20, 28, 50, 82, 126
```

### Verdict: **first 3 match naturally, rest require spin-orbit**

- ✓ Magic numbers **2, 8, 20** emerge from 3D harmonic oscillator structure
  (which GWT has via cube symmetry — d=3 oscillator naturally gives these)
- ✗ Higher magic numbers (28, 50, 82, 126) require spin-orbit-like coupling
  that shifts levels down
- ⚠ GWT has torus-twist which could provide spin-orbit, but the explicit
  derivation showing the shifts give exactly [28, 50, 82, 126] is NOT YET DONE

For this to be a real WIN over standard physics, GWT would need to:
1. Derive 3D harmonic-oscillator-like spectrum from torus modes
2. Include torus twist coupling → spin-orbit analog
3. Calculate level orderings → produce exactly [2, 8, 20, 28, 50, 82, 126]
4. Predict the next magic number (164 in standard physics; GWT might give same or different)

---

## Honest summary

The "shells = breathers" conceptual principle from session_2026-06-03 is:
- **Internally consistent** with both neutron-excess curve and magic numbers
- **Not yet a derivation** that produces these from GWT primitives
- **Compatible with** standard nuclear physics explanations
- **Not yet OUTPERFORMING** standard physics

To turn this conceptual principle into a real framework win:
- 1-2 weeks of nuclear shell model calculation in torus/breather formalism
- Explicit derivation of Coulomb and asymmetry constants from GWT couplings
- Calculation of magic number positions including spin-orbit analog
- Comparison to standard SEMF and shell model coefficients

This is a real research program. The conceptual framing is right; the
mathematics is the work.

---

## Predictive opportunity

If GWT can derive nuclear shell structure from torus breather modes, it
would:
- Unify atomic + nuclear + baryonic structure under one mechanism
- Potentially predict the next magic number beyond 126
- Predict properties of superheavy nuclei
- Be major paper material

But this requires the actual work, not just consistency checks.

---

*Status: tests show framework consistency, not derivation. The shells=
breathers hypothesis remains OPEN as a research target.*
