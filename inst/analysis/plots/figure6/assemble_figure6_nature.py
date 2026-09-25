"""Lay Figure 6 out as a two-by-two block of panels.

Printed type is set by each panel's native width, so this script only places
the panels.

Usage: python3 assemble_figure6_nature.py [panel_dir] [out.pdf]

panel_dir defaults to $BAYRC_FIGURE_DIR/figure6, where replot_figure6.R,
figure6_panelC.R and explore_phase_groups.R write the panels; out.pdf defaults
to Figure_6.pdf in that directory. Requires PyMuPDF (pip install pymupdf).
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
out = Path(sys.argv[2]) if len(sys.argv) > 2 else here / "Figure_6.pdf"

PANELS = ["Fig6A_genomewide_nature.pdf", "Fig6B_circadian_nature.pdf",
          "Fig6C_circadian_membership_nature.pdf",
          "option8_balanced_global_nature.pdf"]
# slot width per row; the upper row is drawn narrower
ROW_SLOT = [480, 520]
SLOT, MARGIN, GUTTER, GAP_X, GAP_Y, LETTER = max(ROW_SLOT), 5, 15, 15, 62, 18
# the top row's shared concordance scale, centred in the band between the rows
SHARED, SHARED_SCALE = "Fig6_concordance_legend.pdf", 1.15
# panel letters use an embedded font file where one is found, else Helvetica Bold
LETTER_FONT = os.environ.get("BAYRC_LETTER_FONT",
                             "/System/Library/Fonts/Supplemental/Arial Bold.ttf")
letter_face = dict(fontname="panel", fontfile=LETTER_FONT) \
    if Path(LETTER_FONT).exists() else dict(fontname="hebo")

src = [fitz.open(here / p) for p in PANELS]
box = [d[0].rect for d in src]
scale = [ROW_SLOT[i // 2] / b.width for i, b in enumerate(box)]
rows = [max(box[i].height * scale[i] for i in (0, 1)),
        max(box[i].height * scale[i] for i in (2, 3))]

width = MARGIN * 2 + (GUTTER + SLOT) * 2 + GAP_X
height = MARGIN * 2 + rows[0] + rows[1] + GAP_Y
doc = fitz.open()
page = doc.new_page(width=width, height=height)

for i, (d, b, s) in enumerate(zip(src, box, scale)):
    x = MARGIN + GUTTER + (GUTTER + SLOT + GAP_X) * (i % 2)
    # panels centred in their slot and row; letters stay on the column edge
    inset = (SLOT - b.width * s) / 2
    y = MARGIN + (rows[0] + GAP_Y) * (i // 2)
    y += (rows[i // 2] - b.height * s) / 2
    page.show_pdf_page(fitz.Rect(x + inset, y, x + inset + b.width * s,
                                 y + b.height * s), d, 0)
    page.insert_text((x - GUTTER + 1, y + LETTER * 0.8), "ABCD"[i],
                     fontsize=LETTER, **letter_face)

key = fitz.open(here / SHARED)
kb = key[0].rect
kw, kh = kb.width * SHARED_SCALE, kb.height * SHARED_SCALE
kx = (width - kw) / 2
ky = MARGIN + rows[0] + (GAP_Y - kh) / 2
page.show_pdf_page(fitz.Rect(kx, ky, kx + kw, ky + kh), key, 0)

doc.save(out, garbage=4, deflate=True)
page.get_pixmap(matrix=fitz.Matrix(2, 2)).save(out.with_suffix(".png"))
print(f"{out.name}  {width:.0f}x{height:.0f} pt, type scales by {488.5 / width:.3f} at \\textwidth")
