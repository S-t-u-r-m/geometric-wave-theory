"""
Experimental De (Equilibrium Dissociation Energy) Reference
=============================================================

CRITICAL: De is the well depth (from bottom of potential to dissociation limit).
D0 is the spectroscopic dissociation energy (from v=0 to dissociation limit).
De = D0 + ZPE, where ZPE ≈ we/2 - wexe/4 for a Morse potential.

The V6 formula computes De (well depth), so ALL comparisons must use De.

Sources:
  [1] NIST Chemistry WebBook: https://webbook.nist.gov/
      Huber & Herzberg compilation (1979)
  [2] NIST CCCBDB: https://cccbdb.nist.gov/expdiatomicsx.asp
      Spectroscopic constants (we, wexe, Be, etc.)
  [3] Active Thermochemical Tables (ATcT): https://atct.anl.gov/
  [4] Luo, Yu-Ran. CRC Bond Dissociation Energies (2010)
  [5] PMC precision spectroscopy papers

Convention:
  - we, wexe in cm-1 (from CCCBDB [2] / NIST WebBook [1])
  - De in eV (= D0 + ZPE, or directly from source)
  - R_e in Angstrom (from NIST WebBook [1])
  - ZPE = we/2 - wexe/4 (Morse approximation)
  - 1 eV = 8065.544 cm-1 (CODATA)
"""

import numpy as np

eV_per_cm = 1.0 / 8065.544

# ============================================================================
# Complete De reference for all 23 V6 molecules
# ============================================================================
# Spectroscopic constants from CCCBDB [2] (consistent with Huber & Herzberg)
# D0 values from NIST WebBook [1] or Huber & Herzberg [1]
# De = D0 + ZPE computed here
#
# Format: name -> {we, wexe, Re, D0_cm, D0_source, De_eV, De_source, quality, notes}
# quality: 'high' (De or D0 from precision spectroscopy, ±0.01 eV)
#          'good' (D0 from H&H, ±0.03 eV)
#          'fair' (D0 uncertain, ±0.1 eV or more)

spectroscopic_data = {
    # =========================================================================
    # HOMONUCLEAR DIATOMICS
    # =========================================================================
    'H2': {
        'we': 4401.213, 'wexe': 121.336, 'Re': 0.7414,
        'D0_cm': 36118.3,
        'D0_source': 'NIST WebBook: very precise from dissociation limit',
        'De_eV': 4.747,  # D0 + ZPE = 36118 + 2170 = 38288 cm-1
        'quality': 'high',
        'notes': 'ZPE=2170 cm-1. De=38288 cm-1=4.747 eV.',
    },
    'Li2': {
        'we': 351.407, 'wexe': 2.583, 'Re': 2.6729,
        'D0_cm': 8517,
        'D0_source': 'Huber & Herzberg [1]',
        'De_eV': 1.078,  # D0 + ZPE = 8517 + 175 = 8692 cm-1
        'quality': 'good',
        'notes': 'ZPE=175 cm-1 (small). De=8692 cm-1=1.078 eV.',
    },
    'B2': {
        'we': 1059.68, 'wexe': 15.66, 'Re': 1.590,
        'D0_cm': 24370,
        'D0_source': 'Huber & Herzberg [1], ±900 cm-1',
        'De_eV': 3.086,  # D0 + ZPE = 24370 + 526 = 24896 cm-1
        'quality': 'fair',  # ±0.11 eV uncertainty on D0
        'notes': 'B2 3Sigma_g ground state. D0 uncertain ±0.11 eV. De=24896 cm-1=3.086 eV.',
    },
    'C2': {
        'we': 1855.066, 'wexe': 13.601, 'Re': 1.2425,
        'D0_cm': 50090,
        'D0_source': 'Huber & Herzberg [1]',
        'De_eV': 6.324,  # D0 + ZPE = 50090 + 924 = 51014 cm-1
        'quality': 'good',
        'notes': 'X 1Sigma+_g. De=51014 cm-1=6.324 eV.',
    },
    'N2': {
        'we': 2358.57, 'wexe': 14.324, 'Re': 1.0977,
        'D0_cm': 78688,
        'De_cm': 79864,
        'D0_source': 'PMC11955043 (2024): precision velocity-map imaging, ±3 cm-1',
        'De_eV': 9.901,  # De=79864 cm-1 directly from paper
        'quality': 'high',
        'notes': 'De=79864±3 cm-1=9.901 eV. ZPE=1176 cm-1. Highest precision.',
    },
    'O2': {
        'we': 1580.161, 'wexe': 11.951, 'Re': 1.2075,
        'D0_cm': 41268,
        'D0_source': 'NIST WebBook [1]',
        'De_eV': 5.214,  # D0 + ZPE = 41268 + 787 = 42055 cm-1
        'quality': 'high',
        'notes': 'ZPE=787 cm-1. De=42055 cm-1=5.214 eV.',
    },
    'F2': {
        'we': 916.929, 'wexe': 11.322, 'Re': 1.4119,
        'D0_cm': 12919,
        'D0_source': 'Huber & Herzberg [1]',
        'De_eV': 1.658,  # D0 + ZPE = 12919 + 456 = 13375 cm-1
        'quality': 'good',
        'notes': 'ZPE=456 cm-1. De=13375 cm-1=1.658 eV.',
    },
    'Na2': {
        'we': 159.086, 'wexe': 0.709, 'Re': 3.0789,
        'D0_cm': 5943,
        'D0_source': 'Huber & Herzberg [1]',
        'De_eV': 0.747,  # D0 + ZPE = 5943 + 79 = 6022 cm-1
        'quality': 'good',
        'notes': 'ZPE tiny (79 cm-1). De=6022 cm-1=0.747 eV.',
    },
    'Cl2': {
        'we': 559.751, 'wexe': 2.694, 'Re': 1.9879,
        'D0_cm': 20007,
        'D0_source': 'Huber & Herzberg [1]',
        'De_eV': 2.515,  # D0 + ZPE = 20007 + 279 = 20286 cm-1
        'quality': 'good',
        'notes': 'ZPE=279 cm-1. De=20286 cm-1=2.515 eV.',
    },

    # =========================================================================
    # HETERONUCLEAR (non-hydride)
    # =========================================================================
    'CO': {
        'we': 2169.756, 'wexe': 13.288, 'Re': 1.1283,
        'D0_cm': 89460,
        'D0_source': 'Huber & Herzberg [1], well-established',
        'De_eV': 11.226,  # D0 + ZPE = 89460 + 1082 = 90542 cm-1
        'quality': 'high',
        'notes': 'ZPE=1082 cm-1. De=90542 cm-1=11.226 eV.',
    },
    'NO': {
        'we': 1904.135, 'wexe': 14.088, 'Re': 1.1508,
        'D0_cm': 52405,
        'D0_source': 'NIST WebBook [1]: https://webbook.nist.gov/cgi/cbook.cgi?ID=C10102439&Mask=1000',
        'De_eV': 6.615,  # D0 + ZPE = 52405 + 949 = 53354 cm-1
        'quality': 'high',
        'notes': 'D0=6.497 eV from NIST confirmed. ZPE=949 cm-1. De=53354 cm-1=6.615 eV.',
    },
    'BF': {
        'we': 1402.159, 'wexe': 11.821, 'Re': 1.2626,
        'D0_cm': 62990,
        'D0_source': 'Huber & Herzberg [1]',
        'De_eV': 7.897,  # D0 + ZPE = 62990 + 698 = 63688 cm-1
        'quality': 'good',
        'notes': 'ZPE=698 cm-1. De=63688 cm-1=7.897 eV.',
    },
    'CN': {
        'we': 2068.648, 'wexe': 13.097, 'Re': 1.1718,
        'D0_cm': 62441,
        'D0_source': 'Huber & Herzberg [1], well-established',
        'De_eV': 7.866,  # D0 + ZPE = 62441 + 1031 = 63472 cm-1
        'quality': 'good',
        'notes': 'ZPE=1031 cm-1. De=63472 cm-1=7.866 eV.',
    },
    'NaCl': {
        'we': 364.684, 'wexe': 1.776, 'Re': 2.3609,
        'D0_cm': 34060,
        'D0_source': 'Huber & Herzberg [1]',
        'De_eV': 4.245,  # D0 + ZPE = 34060 + 182 = 34242 cm-1
        'quality': 'good',
        'notes': 'ZPE=182 cm-1 (small). De=34242 cm-1=4.245 eV.',
    },
    'LiF': {
        'we': 910.573, 'wexe': 8.208, 'Re': 1.5639,
        'D0_cm': 47399,
        'D0_source': 'Huber & Herzberg [1], moderate uncertainty',
        'De_eV': 5.932,  # D0 + ZPE = 47399 + 453 = 47852 cm-1
        'quality': 'good',
        'notes': 'ZPE=453 cm-1. De=47852 cm-1=5.932 eV. D0 has some spread in literature.',
    },

    # =========================================================================
    # HYDRIDES (CRITICAL: large ZPE, De >> D0)
    # =========================================================================
    'HF': {
        'we': 4138.385, 'wexe': 89.943, 'Re': 0.9168,
        'De_cm': 47333,
        'De_eV': 5.869,
        'De_source': 'NIST WebBook [1]: https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664393&Mask=1000',
        'quality': 'high',
        'notes': 'De=47333±60 cm-1=5.869 eV CONFIRMED from NIST. ZPE=2047 cm-1.',
    },
    'OH': {
        'we': 3737.761, 'wexe': 84.881, 'Re': 0.9697,
        'De_cm': 37280,
        'De_eV': 4.621,
        'De_source': 'NIST WebBook [1]: https://webbook.nist.gov/cgi/cbook.cgi?ID=C3352576&Mask=1000',
        'quality': 'high',
        'notes': 'De=37280 cm-1=4.621 eV CONFIRMED from NIST. ZPE=1848 cm-1.',
    },
    'HCl': {
        'we': 2990.925, 'wexe': 52.800, 'Re': 1.2746,
        'D0_cm': 35760,
        'D0_source': 'Huber & Herzberg [1]',
        'De_eV': 4.617,  # D0 + ZPE = 35760 + 1482 = 37242 cm-1
        'quality': 'good',
        'notes': 'ZPE=1482 cm-1 (large!). De=37242 cm-1=4.617 eV.',
    },
    'LiH': {
        'we': 1405.498, 'wexe': 23.168, 'Re': 1.5957,
        'D0_cm': 19589,
        'D0_source': 'Huber & Herzberg [1]',
        'De_eV': 2.515,  # D0 + ZPE = 19589 + 697 = 20286 cm-1
        'quality': 'good',
        'notes': 'ZPE=697 cm-1. De=20286 cm-1=2.515 eV. V6 already had De!',
    },
    'BH': {
        'we': 2366.729, 'wexe': 49.340, 'Re': 1.2324,
        'D0_cm': 27480,
        'D0_source': 'Huber & Herzberg [1], ±2000 cm-1 uncertainty',
        'De_eV': 3.565,  # D0 + ZPE = 27480 + 1171 = 28651 cm-1
        'quality': 'fair',  # D0 uncertain ±0.25 eV
        'notes': 'D0 very uncertain (±2000 cm-1). De≈3.565 eV. Spectroscopic De=3.60±0.05 (1937).',
    },
    'CH': {
        'we': 2860.751, 'wexe': 64.439, 'Re': 1.1199,
        'D0_cm': 27950,
        'D0_source': 'Huber & Herzberg [1], ±250 cm-1',
        'De_eV': 3.644,  # D0 + ZPE = 27950 + 1414 = 29364 cm-1
        'quality': 'good',
        'notes': 'ZPE=1414 cm-1. De=29364 cm-1=3.640 eV. ATcT D0≈3.469 eV consistent.',
    },
    'NH': {
        'we': 3282.721, 'wexe': 79.042, 'Re': 1.0362,
        'D0_cm': 27990,
        'D0_source': 'Huber & Herzberg [1], ±200 cm-1',
        'De_eV': 3.671,  # D0 + ZPE = 27990 + 1622 = 29612 cm-1
        'quality': 'good',
        'notes': 'ZPE=1622 cm-1 (large!). De=29612 cm-1=3.671 eV. V6 had 3.57 (WRONG).',
    },
    'NaH': {
        'we': 1171.968, 'wexe': 19.703, 'Re': 1.8874,
        'D0_cm': 15870,
        'D0_source': 'Huber & Herzberg [1]',
        'De_eV': 2.039,  # D0 + ZPE = 15870 + 581 = 16451 cm-1
        'quality': 'good',
        'notes': 'ZPE=581 cm-1. De=16451 cm-1=2.039 eV.',
    },
}


def compute_De(name):
    """Compute De from spectroscopic data."""
    d = spectroscopic_data[name]
    we = d['we']
    wexe = d['wexe']
    ZPE_cm = we/2.0 - wexe/4.0
    ZPE_eV = ZPE_cm * eV_per_cm

    if 'De_eV' in d:
        return d['De_eV'], ZPE_eV, 'direct'
    elif 'D0_cm' in d:
        De_cm = d['D0_cm'] + ZPE_cm
        return De_cm * eV_per_cm, ZPE_eV, 'D0+ZPE'
    elif 'De_cm' in d:
        return d['De_cm'] * eV_per_cm, ZPE_eV, 'De_cm'
    else:
        return None, ZPE_eV, 'unknown'


def print_full_table():
    """Print complete De reference table with quality flags."""
    print("=" * 110)
    print("  EXPERIMENTAL De REFERENCE — Complete sourced table")
    print("  De = D0 + ZPE (well depth), ZPE = we/2 - wexe/4")
    print("=" * 110)
    print(f"\n  {'Mol':>5} {'we':>8} {'wexe':>7} {'ZPE_cm':>7} {'ZPE_eV':>7} "
          f"{'D0_cm':>7} {'De_cm':>7} {'De_eV':>7} {'qual':>5}")
    print("-" * 110)

    for name, d in spectroscopic_data.items():
        we = d['we']
        wexe = d['wexe']
        ZPE_cm = we/2.0 - wexe/4.0
        ZPE_eV = ZPE_cm * eV_per_cm
        De_eV = d.get('De_eV', 0)
        D0_cm = d.get('D0_cm', 0)
        De_cm = d.get('De_cm', D0_cm + ZPE_cm if D0_cm else 0)
        qual = d.get('quality', '?')

        print(f"  {name:>5} {we:8.2f} {wexe:7.3f} {ZPE_cm:7.1f} {ZPE_eV:7.4f} "
              f"{D0_cm:7.0f} {De_cm:7.0f} {De_eV:7.3f} {qual:>5}")


if __name__ == '__main__':
    print_full_table()

    # Also compare with V6 values if available
    try:
        from v6_complete import molecules
        print(f"\n\n{'='*110}")
        print("  COMPARISON: V6 De values vs experimental De reference")
        print(f"{'='*110}")
        print(f"\n  {'Mol':>5} {'V6_De':>7} {'Ref_De':>7} {'diff':>7} {'diff%':>7} {'quality':>7}")
        print("-" * 70)

        for mol in molecules:
            name = mol[0]
            De_v6 = mol[2]

            if name in spectroscopic_data:
                De_ref = spectroscopic_data[name]['De_eV']
                qual = spectroscopic_data[name].get('quality', '?')
                diff = De_v6 - De_ref
                pct = diff / De_ref * 100
                flag = ' !!!' if abs(pct) > 1.0 else ''
                print(f"  {name:>5} {De_v6:7.3f} {De_ref:7.3f} {diff:+7.3f} {pct:+6.2f}% {qual:>7}{flag}")
            else:
                print(f"  {name:>5} {De_v6:7.3f}     ---                  missing")
    except ImportError:
        pass
