#!/usr/bin/env python3
# Convert GWT Markdown papers to LaTeX using pypandoc.
# Install: pip install pypandoc-binary
# Usage: python md_to_tex.py papers/gwt_lagrangian.md
import pypandoc, sys, os

if len(sys.argv) < 2:
    print('Usage: python md_to_tex.py <input.md> [output.tex]')
    sys.exit(1)

inp = sys.argv[1]
outp = sys.argv[2] if len(sys.argv) >= 3 else os.path.splitext(inp)[0] + '.tex'

pypandoc.convert_file(inp, 'latex', outputfile=outp,
    extra_args=['--standalone', '-V', 'geometry:margin=1in',
                '-V', 'author:Jonathan D. Wollenberg',
                '-V', 'date:March 2026'])

sz = os.path.getsize(outp)
print(f'Converted {inp} -> {outp} ({sz} bytes)')
