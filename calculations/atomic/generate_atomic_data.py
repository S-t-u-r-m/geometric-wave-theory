"""Generate the atomic_data.py lookup table from z_eff_v20.py"""
import sys, io, os, numpy as np
from math import factorial, comb

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

d = 3
gamma = np.pi / (2**(d+1) * np.pi - 2)
alpha_em = np.exp(-(2/factorial(d)) * (2**(2*d+1)/np.pi**2 + np.log(2*d)))
E_H = (alpha_em**2 / 2) * 0.51100e6

# Read the source and extract functions + atoms list
src_path = os.path.join(os.path.dirname(__file__), 'z_eff_v20.py')
with open(src_path, 'r', encoding='utf-8') as f:
    source = f.read()

# Remove the output/logging parts that cause I/O errors
# Replace the log file opening and report function
source = source.replace(
    'sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding=\'utf-8\', errors=\'replace\')',
    '# stdout redirect removed'
)

# Find and cut at the RUN section
for marker in ['print("GWT V20', 'print("  ALL:']:
    idx = source.find(marker)
    if idx > 0:
        source = source[:idx]
        break

# Execute to get functions and data
exec_env = {}
exec(source, exec_env)

atoms = exec_env['atoms']
ionization_energy = exec_env['ionization_energy']

print(f'Loaded {len(atoms)} atoms from z_eff_v20.py')

# Generate the data file
out_path = os.path.join(os.path.dirname(__file__), 'atomic_data.py')
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('"""\n')
    f.write('GWT Atomic Data -- Z_eff, IE, and configurations for 103 atoms\n')
    f.write('Generated from z_eff_v20.py (Oh tensor product model, 2.6% mean)\n')
    f.write('Usage: from atomic_data import ATOMS, E_H\n')
    f.write('"""\n\n')
    f.write(f'E_H = {E_H:.4f}  # alpha^2 * m_e / 2 in eV\n\n')
    f.write('ATOMS = {\n')

    for Z, sym, IE_obs, config in atoms:
        IE_pred, alpha_val, S_core = ionization_energy(Z, config)
        val_n = max(nn for nn, ll, c in config)
        s_count = sum(c for nn, ll, c in config if nn == val_n and ll == 0)
        dc_val = sum(c for nn, ll, c in config if nn == val_n - 1 and ll == 2)
        is_pd = (s_count == 0 and dc_val == 10)
        n = val_n - 1 if is_pd else val_n
        Z_net = Z - S_core
        Z_eff = Z_net ** alpha_val
        err = (IE_pred - IE_obs) / IE_obs * 100

        f.write(f'    "{sym}": {{"Z": {Z}, "IE_obs": {IE_obs}, "IE_pred": {IE_pred:.4f}, ')
        f.write(f'"err_pct": {err:.2f}, "val_n": {n}, "Z_eff": {Z_eff:.4f}, ')
        f.write(f'"Z_net": {Z_net:.4f}, "alpha": {alpha_val:.6f}, ')
        f.write(f'"config": {config}}},\n')

    f.write('}\n')

print(f'Saved {len(atoms)} atoms to {out_path}')

# Verify
exec(open(out_path, encoding='utf-8').read())
print(f'\nVerification:')
print(f'  H:  Z_eff={ATOMS["H"]["Z_eff"]:.4f}, IE={ATOMS["H"]["IE_pred"]:.3f} eV ({ATOMS["H"]["err_pct"]:+.1f}%)')
print(f'  Na: Z_eff={ATOMS["Na"]["Z_eff"]:.4f}, IE={ATOMS["Na"]["IE_pred"]:.3f} eV ({ATOMS["Na"]["err_pct"]:+.1f}%)')
print(f'  Fe: Z_eff={ATOMS["Fe"]["Z_eff"]:.4f}, IE={ATOMS["Fe"]["IE_pred"]:.3f} eV ({ATOMS["Fe"]["err_pct"]:+.1f}%)')
print(f'  Ca: Z_eff={ATOMS["Ca"]["Z_eff"]:.4f}, IE={ATOMS["Ca"]["IE_pred"]:.3f} eV ({ATOMS["Ca"]["err_pct"]:+.1f}%)')
