"""Assemble the Figure 3 demo and reflow legends from vector PDF sources.

Every crop is read off the page rather than written down, so a panel that moves
or a legend that gains a row still lands correctly.
"""
from pathlib import Path
import shutil
import sys
import fitz

# the assembled figures live beside the analysis output, not in the package
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
        # the search ignores case, and a panel annotation can carry the same
        # words, so the strip's own copy is the lowest one on the page
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


# Figure 3 keeps the existing scatter panels and adds the summary panels.
# The scatter row ends at its lowest axis title and the phase-class key is
# found by its own text.
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

# Figure 6 spreads the keys across one bottom strip and drops the figure title.
src = original("Figure_6.pdf")
page0 = src[0]
groups = legend_groups(page0)
bottom = content_bottom(page0, groups)
title = page0.search_for("Cross-species lung")
if title:
    # cut above whatever the figure keeps rather than a fixed drop below the
    # title, so the panel letters sitting just under it keep their headroom
    keep = [fitz.Rect(b[:4]) for b in page0.get_text("blocks")
            if not title[0].intersects(fitz.Rect(b[:4]))]
    top = max(title[0].y1 + 2, min(r.y0 for r in keep) - 6)
else:
    top = 0
scale, gap = .82, 20
width = sum(g.width for g in groups) * scale + gap * (len(groups) - 1)
out = fitz.open()
page = out.new_page(width=page0.rect.width,
                    height=bottom - top + 14 + max(g.height for g in groups) * scale)
place(page, src, (0, top, page0.rect.width, bottom), 0, 0)
x, y = (page0.rect.width - width) / 2, bottom - top + 8
for g in groups:
    place(page, src, g, x, y, scale)
    x += g.width * scale + gap
out.save(figures / "Figure_6.pdf", garbage=4, deflate=True)
page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5)).save(paper / "demos/Figure_6_preview.png")
print("Wrote Figure 3 demo and compact Figure 5/6 layouts.")
