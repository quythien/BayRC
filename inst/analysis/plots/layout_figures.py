"""Assemble the Figure 3 demo and reflow the Figure 5 legend from vector PDF
sources.

Crops are located from the text on each page.

Usage: python3 layout_figures.py [paper_dir]

paper_dir holds figures/ with the assembler's output and demos/figure3/ with
summary_panels.pdf. Requires PyMuPDF (pip install pymupdf).
"""
from pathlib import Path
import shutil
import sys
try:
    import fitz
except ImportError:
    raise SystemExit("PyMuPDF is required: pip install pymupdf")

# paper_dir defaults to paper/ at the repository root
paper = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    Path(__file__).resolve().parents[3] / "paper"
figures = paper / "figures"
archive = paper / "archive" / "before_legend_layout"
archive.mkdir(parents=True, exist_ok=True)

# the titles packLegend draws above each key, in the order they are packed
LEGEND_TITLES = ["Rhythmicity Status", "Phase Status", "Pr(", "Delta Peak"]


def place(page, src, clip, x, y, scale=1):
    clip = fitz.Rect(clip)
    dest = fitz.Rect(x, y, x + clip.width * scale, y + clip.height * scale)
    page.show_pdf_page(dest, src, 0, clip=clip)


def original(name):
    """The panel build as it left R, kept so a reflow is repeatable."""
    backup = archive / name
    if not backup.exists():
        shutil.copy2(figures / name, backup)
    return fitz.open(backup)


def legend_groups(page):
    """One clip per key: its title, the swatches under it and their labels."""
    titles = []
    for t in LEGEND_TITLES:
        # the lowest match on the page is the legend's own title
        hits = page.search_for(t)
        if hits:
            titles.append(max(hits, key=lambda r: r.y0))
    blocks = [fitz.Rect(b[:4]) for b in page.get_text("blocks")]
    if not titles:
        raise SystemExit(
            "no legend titles found on the page. This happens when figures/ "
            "holds an already reflowed file: empty archive/before_legend_layout "
            "and copy fresh assemble_figures.R output into figures/ first.")
    groups = []
    for t in titles:
        # a key runs down to the next key in its own column, or to the page
        below = [o.y0 for o in titles if o.y0 > t.y0 + 4 and abs(o.x0 - t.x0) < 30]
        ymax = min(below) - 6 if below else page.rect.y1
        own = [b for b in blocks
               if t.x0 - 24 <= b.x0 <= t.x0 + 60 and t.y0 - 6 <= b.y0 <= ymax
               and not any(o != t and o.intersects(b) for o in titles)]
        r = own[0]
        for b in own[1:]:
            r |= b
        groups.append(fitz.Rect(r.x0 - 6, r.y0 - 6, r.x1 + 6, r.y1 + 8))
    return groups


def content_bottom(page, groups):
    """Where the panels stop, which is the first legend title less a gap."""
    return min(g.y0 for g in groups) - 12


# Figure 3: the scatter row and its phase-class key, then the summary panels.
src = original("Figure_3.pdf")
page0 = src[0]
W = page0.rect.width
key = page0.search_for("Phase class")[0]
axis = max(r.y1 for r in page0.search_for("Peak Hour"))
scatters = (0, 0, W, axis + 12)          # scatter row, trimmed to its axis titles
legend = (0, key.y0 - 8, W, key.y1 + 8)  # the shared key, trimmed to itself

summary = fitz.open(paper / "demos" / "figure3" / "summary_panels.pdf")
scale = W / summary[0].rect.width
gap = 10                                  # space between the row and its key
demo = fitz.open()
h = (scatters[3] - scatters[1]) + gap + (legend[3] - legend[1]) + gap \
    + summary[0].rect.height * scale
page = demo.new_page(width=W, height=h)
y = 0
place(page, src, scatters, 0, y);  y += scatters[3] - scatters[1] + gap
place(page, src, legend, 0, y);    y += legend[3] - legend[1] + gap
place(page, summary, summary[0].rect, 0, y, scale)
demo.save(paper / "demos/figure3/Figure_3_demo.pdf", garbage=4, deflate=True)
demo.save(figures / "Figure_3.pdf", garbage=4, deflate=True)
page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5)).save(paper / "demos/figure3/Figure_3_demo.png")

# Figure 5 seats the shared legend below A, within B's vertical extent.
src = original("Figure_5_row.pdf")
page0 = src[0]
groups = legend_groups(page0)
strip = fitz.Rect(groups[0])
for g in groups[1:]:
    strip |= g
bottom = content_bottom(page0, groups)
# panel A is the left half, so its last row is where the legend can start
a_bottom = max(fitz.Rect(b[:4]).y1 for b in page0.get_text("blocks")
               if fitz.Rect(b[:4]).x1 < page0.rect.width / 2 and fitz.Rect(b[:4]).y1 < bottom)
scale = .86
out = fitz.open()
page = out.new_page(width=page0.rect.width, height=bottom + 9)
place(page, src, (0, 0, page0.rect.width, bottom), 0, 0)
place(page, src, strip, 27, a_bottom - 2, scale)
out.save(figures / "Figure_5_row.pdf", garbage=4, deflate=True)
page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5)).save(paper / "demos/Figure_5_row_preview.png")

# Figure 6 is laid out by plots/figure6/assemble_figure6.py and is not touched.
print("Wrote Figure 3 demo and compact Figure 5 layout.")
