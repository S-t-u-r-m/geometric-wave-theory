#!/usr/bin/env python3
# Convert GWT Markdown papers to LaTeX and PDF using pypandoc.
# Install: pip install pypandoc-binary
# Usage: python md_to_tex.py papers/gwt_lagrangian.md
#   Generates both .tex and .pdf in the same directory.
import pypandoc, sys, os

if len(sys.argv) < 2:
    print("Usage: python md_to_tex.py <input.md> [output.tex]")
    sys.exit(1)

inp = sys.argv[1]
base = os.path.splitext(inp)[0]
outp_tex = sys.argv[2] if len(sys.argv) >= 3 else base + ".tex"
outp_pdf = os.path.splitext(outp_tex)[0] + ".pdf"

# LaTeX header for table styling and unicode
header = r"""
\usepackage{booktabs}
\usepackage{array}
\usepackage{textcomp}
\DeclareUnicodeCharacter{2264}{$\leq$}
\DeclareUnicodeCharacter{2265}{$\geq$}
\DeclareUnicodeCharacter{2212}{$-$}
\DeclareUnicodeCharacter{2013}{--}
\DeclareUnicodeCharacter{2014}{---}
\renewcommand{\arraystretch}{1.4}
\setlength{\tabcolsep}{8pt}
"""

# Write header to temp file
header_file = os.path.join(os.path.dirname(inp) or ".", "_header.tex")
with open(header_file, "w") as f:
    f.write(header)

common_args = ["--standalone", "-V", "geometry:margin=1in",
               "-V", "fontfamily:charter", "-V", "mathfont:charter",
               "-H", header_file]

# Generate .tex
pypandoc.convert_file(inp, "latex", outputfile=outp_tex, extra_args=common_args)
sz = os.path.getsize(outp_tex)
print(f"TeX:  {outp_tex} ({sz} bytes)")

# Generate .pdf
try:
    pypandoc.convert_file(inp, "pdf", outputfile=outp_pdf,
        extra_args=common_args + ["--pdf-engine=pdflatex"])
    sz = os.path.getsize(outp_pdf)
    print(f"PDF:  {outp_pdf} ({sz} bytes)")
except Exception as e:
    print(f"PDF generation failed: {e}")
    print("Install a TeX distribution (MiKTeX) to enable PDF output.")
