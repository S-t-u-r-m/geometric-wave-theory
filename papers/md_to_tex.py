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
\DeclareUnicodeCharacter{03B1}{$\alpha$}
\DeclareUnicodeCharacter{03B2}{$\beta$}
\DeclareUnicodeCharacter{03B3}{$\gamma$}
\DeclareUnicodeCharacter{03B4}{$\delta$}
\DeclareUnicodeCharacter{03B5}{$\epsilon$}
\DeclareUnicodeCharacter{03B7}{$\eta$}
\DeclareUnicodeCharacter{03B8}{$\theta$}
\DeclareUnicodeCharacter{03BB}{$\lambda$}
\DeclareUnicodeCharacter{03BC}{$\mu$}
\DeclareUnicodeCharacter{03BD}{$\nu$}
\DeclareUnicodeCharacter{03C0}{$\pi$}
\DeclareUnicodeCharacter{03C1}{$\rho$}
\DeclareUnicodeCharacter{03C3}{$\sigma$}
\DeclareUnicodeCharacter{03C4}{$\tau$}
\DeclareUnicodeCharacter{03C6}{$\phi$}
\DeclareUnicodeCharacter{03C7}{$\chi$}
\DeclareUnicodeCharacter{03C8}{$\psi$}
\DeclareUnicodeCharacter{03C9}{$\omega$}
\DeclareUnicodeCharacter{0393}{$\Gamma$}
\DeclareUnicodeCharacter{0394}{$\Delta$}
\DeclareUnicodeCharacter{039B}{$\Lambda$}
\DeclareUnicodeCharacter{03A3}{$\Sigma$}
\DeclareUnicodeCharacter{03A9}{$\Omega$}
\DeclareUnicodeCharacter{2202}{$\partial$}
\DeclareUnicodeCharacter{210F}{$\hbar$}
\DeclareUnicodeCharacter{2248}{$\approx$}
\DeclareUnicodeCharacter{2260}{$\neq$}
\DeclareUnicodeCharacter{00B1}{$\pm$}
\DeclareUnicodeCharacter{00D7}{$\times$}
\DeclareUnicodeCharacter{2192}{$\rightarrow$}
\DeclareUnicodeCharacter{221A}{$\sqrt{}$}
\DeclareUnicodeCharacter{221E}{$\infty$}
\DeclareUnicodeCharacter{2272}{$\lesssim$}
\DeclareUnicodeCharacter{2273}{$\gtrsim$}
\DeclareUnicodeCharacter{2261}{$\equiv$}
\DeclareUnicodeCharacter{207A}{$^{+}$}
\DeclareUnicodeCharacter{207B}{$^{-}$}
\DeclareUnicodeCharacter{2070}{$^{0}$}
\DeclareUnicodeCharacter{00B2}{$^{2}$}
\DeclareUnicodeCharacter{00B3}{$^{3}$}
\DeclareUnicodeCharacter{2074}{$^{4}$}
\DeclareUnicodeCharacter{2075}{$^{5}$}
\DeclareUnicodeCharacter{2076}{$^{6}$}
\DeclareUnicodeCharacter{00B9}{$^{1}$}
\DeclareUnicodeCharacter{2077}{$^{7}$}
\DeclareUnicodeCharacter{2078}{$^{8}$}
\DeclareUnicodeCharacter{2079}{$^{9}$}
\DeclareUnicodeCharacter{2080}{$_{0}$}
\DeclareUnicodeCharacter{2081}{$_{1}$}
\DeclareUnicodeCharacter{2082}{$_{2}$}
\DeclareUnicodeCharacter{2083}{$_{3}$}
\DeclareUnicodeCharacter{2084}{$_{4}$}
\DeclareUnicodeCharacter{2089}{$_{9}$}
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
