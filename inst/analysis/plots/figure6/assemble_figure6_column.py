"""Lay Figure 6 out as a column of panels beside the phase heatmap.

A and B share the top of the left column with the concordance scale under them,
the membership panel takes the width below, and the phase heatmap stands the
full height of the page on the right. An alternative to the two-by-two layout
of assemble_figure6_nature.py, built from the same panel directory.

Usage: python3 assemble_figure6_column.py [panel_dir] [out.pdf]

panel_dir defaults to $BAYRC_FIGURE_DIR/figure6 and out.pdf to
Figure_6_column.pdf in that directory. Requires PyMuPDF (pip install pymupdf).
"""
from pathlib import Path
import os
import sys
try:
    import fitz
except ImportError:
    raise SystemExit("PyMuPDF is required: pip install pymupdf")

if len(sys.argv) > 1:
    here = Path(sys.argv[1])
elif os.environ.get("BAYRC_FIGURE_DIR"):
    here = Path(os.environ["BAYRC_FIGURE_DIR"]) / "figure6"
else:
    raise SystemExit("give the panel directory, or set BAYRC_FIGURE_DIR")
out = Path(sys.argv[2]) if len(sys.argv) > 2 else here / "Figure_6_column.pdf"

LEFT, RIGHT, MARGIN, GAP, LETTER = 675, 446, 10, 14, 19
src = {n: fitz.open(here / f) for n, f in [
    ("A", "Fig6A_genomewide.pdf"), ("B", "Fig6B_circadian.pdf"),
    ("key", "Fig6_concordance_legend.pdf"),
    ("C", "Fig6C_circadian_membership.pdf"),
    ("D", "option8_balanced_global_long.pdf")]}
box = {n: d[0].rect for n, d in src.items()}

half = (LEFT - GAP) / 2
top = half * box["A"].height / box["A"].width
key_w = 400
key_h = key_w * box["key"].height / box["key"].width
c_h = LEFT * box["C"].height / box["C"].width
left_h = top + GAP + key_h + GAP + c_h
d_scale = RIGHT / box["D"].width
height = MARGIN * 2 + max(left_h, box["D"].height * d_scale)
width = MARGIN * 2 + LEFT + GAP + RIGHT

doc = fitz.open()
page = doc.new_page(width=width, height=height)


def place(name, x, y, w):
    b, s = box[name], w / box[name].width
    page.show_pdf_page(fitz.Rect(x, y, x + w, y + b.height * s), src[name], 0)


def letter(ch, x, y):
    page.insert_text((x, y + LETTER * 0.8), ch, fontname="hebo", fontsize=LETTER)


y = MARGIN
place("A", MARGIN, y, half)
place("B", MARGIN + half + GAP, y, half)
letter("A", MARGIN, y)
letter("B", MARGIN + half + GAP, y)
y += top + GAP
place("key", MARGIN + (LEFT - key_w) / 2, y, key_w)
y += key_h + GAP
place("C", MARGIN, y, LEFT)
letter("C", MARGIN, y)
place("D", MARGIN + LEFT + GAP, MARGIN, RIGHT)
letter("D", MARGIN + LEFT + GAP, MARGIN)

doc.save(out, garbage=4, deflate=True)
page.get_pixmap(matrix=fitz.Matrix(1.6, 1.6)).save(out.with_suffix(".png"))
print(f"{out.name}  {width:.0f}x{height:.0f} pt, type scales by {488.5 / width:.3f} at \\textwidth")
