"""Assemble Figure 5 from the two scatters and the three-species heatmap.

Printed type follows the page width, so the two arrangements are written side
by side and the type each one reaches is reported with them. The paper carries
Figure_5_column.pdf.

Usage: python3 assemble_figure5.py [panel_dir]

panel_dir defaults to $BAYRC_FIGURE_DIR/figure5, where make_figure5.R writes
the panels. Requires PyMuPDF (pip install pymupdf).
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
    here = Path(os.environ["BAYRC_FIGURE_DIR"]) / "figure5"
else:
    raise SystemExit("give the panel directory, or set BAYRC_FIGURE_DIR")
PANELS = ["Fig5A_human_baboon.pdf", "Fig5B_human_mouse.pdf", "Fig5C_heatmap.pdf"]
GAP, MARGIN, GUTTER, LETTER = 16, 6, 18, 24
TEXTWIDTH = 488.5

src = [fitz.open(here / p) for p in PANELS]
box = [d[0].rect for d in src]


def place(page, d, b, x, y, w):
    s = w / b.width
    page.show_pdf_page(fitz.Rect(x, y, x + w, y + b.height * s), d, 0)
    return b.height * s


def letter(page, ch, x, y):
    page.insert_text((x, y + LETTER * 0.8), ch, fontname="hebo", fontsize=LETTER)


def stacked(out):
    """Scatters across the top, heatmap at full width beneath them."""
    w = box[2].width
    half = (w - GAP - GUTTER) / 2
    doc = fitz.open()
    top = box[0].height * half / box[0].width
    page = doc.new_page(width=MARGIN * 2 + GUTTER + w,
                        height=MARGIN * 2 + top + GAP + box[2].height)
    x0 = MARGIN + GUTTER
    place(page, src[0], box[0], x0, MARGIN, half)
    letter(page, "A", MARGIN, MARGIN)
    place(page, src[1], box[1], x0 + half + GAP, MARGIN, half)
    letter(page, "B", x0 + half + GAP - GUTTER, MARGIN)
    place(page, src[2], box[2], x0, MARGIN + top + GAP, w)
    letter(page, "C", MARGIN, MARGIN + top + GAP)
    doc.save(out, garbage=4, deflate=True)
    return page.rect.width


def column(out):
    """Scatters in a left column, heatmap filling the height beside them."""
    col_h = box[0].height + GAP + box[1].height
    c_w = box[2].width * col_h / box[2].height
    doc = fitz.open()
    page = doc.new_page(width=MARGIN * 2 + GUTTER * 2 + box[0].width + GAP + c_w,
                        height=MARGIN * 2 + col_h)
    x0 = MARGIN + GUTTER
    place(page, src[0], box[0], x0, MARGIN, box[0].width)
    letter(page, "A", MARGIN, MARGIN)
    place(page, src[1], box[1], x0, MARGIN + box[0].height + GAP, box[1].width)
    letter(page, "B", MARGIN, MARGIN + box[0].height + GAP)
    xc = x0 + box[0].width + GAP + GUTTER
    place(page, src[2], box[2], xc, MARGIN, c_w)
    letter(page, "C", xc - GUTTER, MARGIN)
    doc.save(out, garbage=4, deflate=True)
    return page.rect.width


for name, fn in [("Figure_5_stacked.pdf", stacked), ("Figure_5_column.pdf", column)]:
    w = fn(here / name)
    p = fitz.open(here / name)[0]
    sz = [s["size"] * TEXTWIDTH / w for b in p.get_text("dict")["blocks"]
          for l in b.get("lines", []) for s in l["spans"] if s["text"].strip()]
    p.get_pixmap(matrix=fitz.Matrix(1.1, 1.1)).save(str(here / name).replace(".pdf", ".png"))
    print(f"  {name:24s} {w:5.0f} x {p.rect.height:4.0f} pt   min type {min(sz):.2f} pt   "
          f"under 5 pt: {sum(1 for x in sz if x < 5)} of {len(sz)}")
