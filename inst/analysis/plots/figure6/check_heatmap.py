"""Check a three-condition heatmap before it is used.

Reads the drawn page rather than the script that made it, so a block that is
built from a different rule than its neighbour is caught.
Usage: python3 check_heatmap.py <pdf> [expected_group ...]
Requires PyMuPDF (pip install pymupdf).
"""
import sys
import collections
try:
    import fitz
except ImportError:
    raise SystemExit("PyMuPDF is required: pip install pymupdf")

GREY = (0.88, 0.88, 0.88)


def spans(page):
    return [(fitz.Rect(s["bbox"]), s["text"].strip())
            for b in page.get_text("dict")["blocks"] for l in b.get("lines", [])
            for s in l["spans"] if s["text"].strip()]


def main(path, expect):
    page = fitz.open(path)[0]
    txt = spans(page)
    fails = []

    offsets = sorted([(r, t) for r, t in txt if "(h)" in t and "-" in t.replace("−", "-")],
                     key=lambda rt: rt[0].x0)
    if len(offsets) < 2:
        fails.append(f"expected two offset blocks, found {len(offsets)}")

    counts = []
    for r, t in offsets:
        n = sum(1 for d in page.get_drawings()
                if d.get("fill") and tuple(round(c, 2) for c in d["fill"]) == GREY
                and r.x0 - 55 < d["rect"][0] < r.x1 + 55 and d["rect"][3] < r.y0 - 5
                and (d["rect"][3] - d["rect"][1]) < 40)
        counts.append((t, n))
    # the blocks share one rule, so neither may be empty while the other is full
    if len(counts) == 2:
        a, b = counts[0][1], counts[1][1]
        if max(a, b) and min(a, b) < 0.2 * max(a, b):
            fails.append(f"offset blocks drawn from different rules: {counts}")

    for g in expect:
        if not any(g in t for _, t in txt):
            fails.append(f"missing label: {g}")

    rows = collections.defaultdict(list)
    for r, t in txt:
        rows[round(r.y0)].append((r, t))
    for y, items in rows.items():
        items.sort(key=lambda rt: rt[0].x0)
        for (r1, t1), (r2, t2) in zip(items, items[1:]):
            if r1.x1 > r2.x0 + 1:
                fails.append(f"overlap at y={y}: {t1!r} into {t2!r}")

    print(f"  {path}")
    for t, n in counts:
        print(f"    offset block {t!r}: {n} unclassified bars")
    if fails:
        print("  FAIL")
        for f in fails[:8]:
            print("   -", f)
        return 1
    print("  PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2:]))
