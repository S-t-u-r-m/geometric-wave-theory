# Slider applied to chemistry — failure analysis 2026-06-04

## Context

Jon's intuition: "applying this to bonding will give us a glimpse of the smooth transition"

The baryon slider (delta = (1/d) * max(0, (m_h - m_J/psi)/m_h)) closes
Σ-Λ type splittings at sub-1% across s, c, b. But the transition shape
is undetermined from 3 data points.

If the slider mechanism is real and universal, it should also apply
to chemistry V10 bonds. Testing this reveals:

## What we found

### V10 errors are BIDIRECTIONAL in chemistry

Unlike baryons (all over-predict before slider), chemistry V10 errors
go both ways:

**Over-predictions** (slider direction):
- H2: +6%, Li2: +80%, N2: +4%, O2: +5%, Cl2: +21%, S2: +10%
- LiH: +59%, PH: +13%, NH: +5%, SH: +3%

**Under-predictions** (slider can't fix):
- F2: -1.6%, NaCl: -13%, CO: -17%, NO: -11%, CH: -9%

The slider mechanism (which only DECREASES bond energy by rotating
w_pi toward less-binding angle) cannot fix under-predictions.

### Applying baryon slider directly: catastrophic over-correction

With delta = (1/d) * max(0, (mu_lh - m_J/psi)/mu_lh):

| Bond | V10 err | After slider | 
|------|---------|--------------|
| N2   | +4%     | -29%         |
| O2   | +5%     | -32%         |
| Cl2  | +21%    | -40%         |
| S2   | +10%    | -43%         |
| NaCl | -13%    | -54%         |

The (1/d) coefficient is far too aggressive for chemistry.

### Required slope is ~10x gentler

For Cl2 to close, needed delta = 5.5° at mu = 16511 MeV.
Implied coefficient X ≈ 1/8.5 (vs 1/3 for baryons).

For S2: X ≈ 1/14
For N2: X ≈ 1/16

Chemistry "slider coefficient" varies bond-by-bond. Not a clean
single-parameter mechanism.

## What this reveals

The baryon transition LOOKED sharp because:
- Only 3 data points (s, c, b)
- All errors monotonic in one direction
- Single mechanism dominates

Chemistry reveals:
- MULTIPLE competing mechanisms give bidirectional errors
- The slider, if real, has a MUCH gentler slope in chemistry
- The "smooth transition" Jon predicted is real — emerges from
  mechanism competition, not a single smooth function

## Implications

1. The baryon slider with sharp threshold at m_J/psi is a special
   case where one mechanism dominates
2. Chemistry needs at least TWO independent correction mechanisms
3. The "smooth transition" perception is correct - but it comes from
   superposed effects, not a gradual single-knob function
4. The chemistry slider, if isolated, would activate gently across
   all bonds (no clear threshold) rather than turning on sharply

## Status

**Open**: identifying the second chemistry mechanism that causes
under-predictions for F2, NaCl, CO, NO, CH. Possible candidates:
- Polar/ionic mixing not captured by current c_ionic
- Multiple bond order effects (CO, NO are partial multi-bonded)
- Specific anti-bonding orbital occupation (F2 has lone pair repulsion)

The slider extension to chemistry is INFORMATIVE about the transition
shape (gentle, multi-mechanism) but does not directly close the
chemistry residuals.
