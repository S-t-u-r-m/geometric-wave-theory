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

common_args = ["--standalone", "-V", "geometry:margin=1in"]

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
