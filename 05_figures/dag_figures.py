"""Figures drawn on the hand-tuned DAG layout (fig2_layout.json, saved from the layout editor).

Figure 3: where structure recovery failed, one panel per data set, plus a key.
Figure 5: where each attribution method put its credit, full-set, one panel per truth or method.

Both write SVG with live text; PDFs are made from the SVG by a browser so the label font is embedded.
"""
import argparse
import json
import math
import os
from pathlib import Path
from xml.sax.saxutils import escape

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("PIPELINE_ROOT", HERE.parent))
OUT = ROOT / "outputs" / "05_figures"
STRUCTURE = ROOT / "outputs" / "03_structure_validation"
ATTR = ROOT / "outputs" / "04_attribution_validation"
DATA = ROOT / "02_data"

W, H = 1830, 1000
PAD, INNER = 70, 0.92  # panel inset so labels at the left edge stay inside
FONT = "Barlow Semi Condensed, TeX Gyre Heros Cn, Arial Narrow, Arial, sans-serif"
INK, NODE, FAINT = "#0b0d12", "#00a6d6", "#c4ccd6"
DATASETS = ["full_set", "reference_subsample", "astronaut_set"]
DS_LABEL = {"full_set": "Full-set", "reference_subsample": "Reference-set", "astronaut_set": "Astronaut-set"}
OUTCOME = "nephrolithiasis"

# recovery categories: colour, dash pattern, arrow, legend label
CATS = {
    "absent": ("#1b1e24", "14 9", False, "Missed"),
    "reversed": ("#d55e00", None, True, "Recovered, direction reversed"),
    "overoriented": ("#e69f00", "18 6 4 6", True, "Direction imposed on an undirected edge"),
    "unoriented": ("#8a5cc2", None, False, "Left undirected"),
    "added": ("#cc79a7", "3 7", False, "Added (not in the true graph)"),
}

METHODS = [("truth", "True effect"), ("standard_shap", "Standard SHAP"), ("ordering_only_shap", "Ordering-only SHAP"),
           ("ng_style_shap", "Causal predictive SHAP"), ("interventional_shap", "Causal interventional SHAP")]


def load_layout(path):
    d = json.loads(Path(path).read_text())
    d = d.get("data", d)
    return d["nodes"], d.get("edges", {}), d["style"]


def label_lines(label):
    return [t.strip() for t in label.split("|")]


class Canvas:
    """One SVG panel in layout units, placed inside a larger figure by a transform."""

    def __init__(self, nodes, bends, style, fs=None):
        self.nodes, self.bends, self.s = nodes, bends, dict(style)
        if fs:
            self.s["fs"] = fs
        self.parts = []

    def edge(self, a, b, color, width, dash=None, arrow=True, reverse=False, opacity=1.0, r=None):
        if reverse:
            a, b = b, a
        A, B = self.nodes[a], self.nodes[b]
        r = r or self.s["nr"]
        bend = self.bends.get(f"{a}__{b}") or self.bends.get(f"{b}__{a}")
        mx, my = (A["x"] + B["x"]) / 2, (A["y"] + B["y"]) / 2
        cx, cy = (mx + bend["dx"], my + bend["dy"]) if bend else (mx, my)
        ux, uy = unit(cx - A["x"], cy - A["y"])
        vx, vy = unit(B["x"] - cx, B["y"] - cy)
        sx, sy = A["x"] + ux * (r + 2), A["y"] + uy * (r + 2)
        ex, ey = B["x"] - vx * (r + (4 if arrow else 2)), B["y"] - vy * (r + (4 if arrow else 2))
        d = f"M{sx:.1f},{sy:.1f} Q{cx:.1f},{cy:.1f} {ex:.1f},{ey:.1f}" if bend else f"M{sx:.1f},{sy:.1f} L{ex:.1f},{ey:.1f}"
        attrs = f'd="{d}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linecap="round"'
        if dash:
            attrs += f' stroke-dasharray="{dash}"'
        if opacity < 1:
            attrs += f' stroke-opacity="{opacity}"'
        if arrow:
            attrs += f' marker-end="url(#arr{color.strip("#")})"'
        self.parts.append(f"<path {attrs}/>")
        return color if arrow else None

    def node(self, nid, fill, stroke, r, ring=False, width=1.4):
        n = self.nodes[nid]
        self.parts.append(f'<circle cx="{n["x"]}" cy="{n["y"]}" r="{r:.1f}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>')
        if ring:
            self.parts.append(f'<circle cx="{n["x"]}" cy="{n["y"]}" r="{r + 3.5:.1f}" fill="none" stroke="{stroke}" stroke-width="{width * 0.7}"/>')

    def label(self, nid, r, color=INK, weight=700):
        n, s = self.nodes[nid], self.s
        lines = label_lines(n["label"])
        lh = s["fs"] * s["ld"]
        # the editor's offsets were set at the editor's type size; scale them with the type
        k = s["fs"] / s.get("fs0", s["fs"])
        # offset of the first baseline as the editor set it, at the editor's dot size
        base = -s["nr"] - s["gap"] * k + n.get("ly", 0) * k
        grow = max(r - s["nr"], 0)
        y0 = n["y"] + base + (grow if base > 0 else -grow) - (len(lines) - 1) * lh
        x = n["x"] + n.get("lx", 0) * k
        spans = "".join(f'<tspan x="{x:.1f}" dy="{0 if i == 0 else lh:.1f}">{escape(t)}</tspan>' for i, t in enumerate(lines))
        self.parts.append(f'<text x="{x:.1f}" y="{y0:.1f}" text-anchor="middle" font-family="{FONT}" font-size="{s["fs"]}" '
                          f'font-weight="{weight}" fill="{color}">{spans}</text>')

    def svg(self):
        return "\n".join(self.parts)


def unit(x, y):
    l = math.hypot(x, y) or 1
    return x / l, y / l


def markers(colors):
    m = []
    for c in sorted(set(colors)):
        m.append(f'<marker id="arr{c.strip("#")}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" '
                 f'markerUnits="strokeWidth" orient="auto"><path d="M0,1 L10,5 L0,9 z" fill="{c}"/></marker>')
    return "<defs>" + "".join(m) + "</defs>"


def text(x, y, s, size, weight=700, anchor="start", color=INK):
    return (f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" font-weight="{weight}" '
            f'text-anchor="{anchor}" fill="{color}">{escape(s)}</text>')


def document(width_mm, vb_w, vb_h, body, colors):
    h_mm = width_mm * vb_h / vb_w
    return (f'<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {vb_w} {vb_h}" '
            f'width="{width_mm}mm" height="{h_mm:.1f}mm">{markers(colors)}<rect width="{vb_w}" height="{vb_h}" fill="#ffffff"/>'
            f"{body}</svg>")


# ---------- Figure 3: recovery failures per data set ----------

def recovery(added_threshold):
    d = pd.read_csv(STRUCTURE / "difference_dag_edges.csv")
    true = d[d["true_edge"].astype(str).str.upper() == "TRUE"].copy()
    rec = pd.read_csv(STRUCTURE / "edge_recovery.csv")
    rec = rec[(rec["variant"] == "MGM-PC-Stable") & (rec["true_edge"].astype(str).str.upper() != "TRUE")]
    added = rec[rec["adjacency_rate"] >= added_threshold].copy()
    added[["a", "b"]] = added["edge"].str.split(" ", expand=True)
    return true, added


def fig3(layout, added_threshold, fs):
    nodes, bends, style = layout
    true, added = recovery(added_threshold)
    pw, ph, gap_x, gap_y, title_h = W, H, 90, 60, 70
    vb_w, vb_h = 2 * pw + gap_x, 2 * (ph + title_h) + gap_y
    body, colors = [], set()
    counts = {}
    slots = [(0, 0), (1, 0), (0, 1)]
    for (cx, cy), ds in zip(slots, DATASETS):
        c = Canvas(nodes, bends, {**style, "fs0": style["fs"]}, fs=fs)
        t = true[true["dataset"] == ds]
        hit = set()
        for e in t.itertuples():
            if e.modal_category == "correct":
                c.edge(e.parent, e.child, FAINT, 2.2, dash="2 7", arrow=False)
        cnt = {}
        for e in t.itertuples():
            cat = e.modal_category
            if cat == "correct":
                continue
            col, dash, arrow, _ = CATS[cat]
            arrow = arrow and e.true_mark != "undir" or cat == "overoriented"
            m = c.edge(e.parent, e.child, col, 3.4, dash=dash, arrow=arrow, reverse=(cat == "reversed"))
            if m:
                colors.add(m)
            hit |= {e.parent, e.child}
            cnt[cat] = cnt.get(cat, 0) + 1
        for e in added[added["dataset"] == ds].itertuples():
            col, dash, _, _ = CATS["added"]
            c.edge(e.a, e.b, col, 4.5, dash=dash, arrow=False)
            hit |= {e.a, e.b}
            cnt["added"] = cnt.get("added", 0) + 1
        counts[ds] = cnt
        r = style["nr"] * 1.15
        for nid in nodes:
            on = nid in hit
            if nid == OUTCOME:
                c.node(nid, NODE, NODE, r * 1.45, ring=True)
            else:
                c.node(nid, NODE if on else "#ffffff", NODE if on else FAINT, r, width=2)
        for nid in hit | {OUTCOME}:
            c.label(nid, r * (1.45 if nid == OUTCOME else 1))
        x0, y0 = cx * (pw + gap_x), cy * (ph + title_h + gap_y)
        letter = "abc"[DATASETS.index(ds)]
        body.append(text(x0 + 6, y0 + 46, f"{letter}", 48))
        body.append(text(x0 + 56, y0 + 46, DS_LABEL[ds], 40, weight=600))
        body.append(f'<g transform="translate({x0 + PAD},{y0 + title_h}) scale({INNER})">{c.svg()}</g>')
    # key panel
    x0, y0 = pw + gap_x, ph + title_h + gap_y
    body.append(text(x0 + 6, y0 + 46, "Key", 40, weight=600))
    yy = y0 + 150
    for cat in ["absent", "reversed", "overoriented", "unoriented", "added"]:
        col, dash, arrow, lab = CATS[cat]
        da = f' stroke-dasharray="{dash}"' if dash else ""
        mk = f' marker-end="url(#arr{col.strip("#")})"' if arrow else ""
        if arrow:
            colors.add(col)
        body.append(f'<path d="M{x0 + 20},{yy} L{x0 + 200},{yy}" stroke="{col}" stroke-width="5" fill="none"{da}{mk}/>')
        body.append(text(x0 + 230, yy + 13, lab, 38, weight=500))
        yy += 78
    body.append(f'<path d="M{x0 + 20},{yy} L{x0 + 200},{yy}" stroke="{FAINT}" stroke-width="3" stroke-dasharray="2 7" fill="none"/>')
    body.append(text(x0 + 230, yy + 13, "Recovered correctly", 38, weight=500))
    yy += 110
    head = ["", "Full", "Ref.", "Astro."]
    cols = [x0 + 20, x0 + 640, x0 + 780, x0 + 920]
    for x, h in zip(cols, head):
        body.append(text(x, yy, h, 34, weight=700, anchor="start" if x == cols[0] else "end"))
    yy += 52
    for cat in ["absent", "reversed", "overoriented", "unoriented", "added"]:
        body.append(text(cols[0], yy, CATS[cat][3], 34, weight=500))
        for x, ds in zip(cols[1:], DATASETS):
            body.append(text(x, yy, str(counts[ds].get(cat, 0)), 34, weight=500, anchor="end"))
        yy += 48
    body.append(text(cols[0], yy + 20, "Edges by most frequent outcome over 50 replicates (of 63).", 30, weight=500, color="#4a5260"))
    svg = document(183, vb_w, vb_h, "\n".join(body), colors | {NODE})
    return svg, counts


# ---------- Figure 5: where credit went ----------

def credit_shares(dataset):
    imp = pd.read_csv(ATTR / "importance.csv")
    imp = imp[(imp["dataset"] == dataset) & (imp["feature_set"] == "all")]
    imp = imp.assign(a=imp["importance"].abs())
    imp["share"] = imp["a"] / imp.groupby(["method", "replicate"])["a"].transform("sum")
    sh = imp.groupby(["method", "feature"])["share"].mean().unstack(0)
    truth = pd.read_csv(DATA / "truth_values.csv")
    pop = "selected" if dataset == "astronaut_set" else "source"
    t = truth[(truth["truth_type"] == "per_unit") & (truth["population"] == pop)].set_index("feature")["value"].abs()
    sh["truth"] = (t / t.sum()).reindex(sh.index).fillna(0)
    return sh


def colour(v, vmax):
    # white to deep cyan-blue, by square root so small shares stay visible
    f = min(1.0, math.sqrt(max(v, 0) / vmax)) if vmax else 0
    lo, hi = (255, 255, 255), (6, 54, 122)
    mid = (0, 166, 214)
    if f < 0.5:
        a, b, g = lo, mid, f / 0.5
    else:
        a, b, g = mid, hi, (f - 0.5) / 0.5
    return "#%02x%02x%02x" % tuple(round(a[i] + (b[i] - a[i]) * g) for i in range(3))


def fig5(layout, dataset, fs, label_min):
    nodes, bends, style = layout
    sh = credit_shares(dataset)
    vmax = 0.30
    pw, ph, gap_x, gap_y, title_h = W, H, 90, 60, 70
    vb_w, vb_h = 2 * pw + gap_x, 3 * (ph + title_h) + 2 * gap_y
    body = []
    edges = list(pd.read_csv(STRUCTURE / "analysis_graph_edges.csv")[["parent", "child"]].itertuples(index=False))
    for i, (m, title) in enumerate(METHODS):
        cx, cy = i % 2, i // 2
        c = Canvas(nodes, bends, {**style, "fs0": style["fs"]}, fs=fs)
        for a, b in edges:
            c.edge(a, b, FAINT, 1.6, arrow=False)
        for nid in nodes:
            v = float(sh.loc[nid, m]) if nid in sh.index else 0.0
            r = 7 + 34 * math.sqrt(v / vmax) if v > 0 else 5
            r = min(r, 44)
            stroke = INK if nid == OUTCOME else "#6b7686"
            c.node(nid, colour(v, vmax) if nid != OUTCOME else "#ffffff", stroke, r, ring=(nid == OUTCOME), width=1.6)
            if v >= label_min or nid == OUTCOME:
                c.label(nid, r)
        x0, y0 = cx * (pw + gap_x), cy * (ph + title_h + gap_y)
        body.append(text(x0 + 6, y0 + 46, "abcde"[i], 48))
        body.append(text(x0 + 56, y0 + 46, title, 40, weight=600))
        body.append(f'<g transform="translate({x0 + PAD},{y0 + title_h}) scale({INNER})">{c.svg()}</g>')
    # colour key
    x0, y0 = pw + gap_x, 2 * (ph + title_h + gap_y)
    body.append(text(x0 + 6, y0 + 46, "Share of credit", 40, weight=600))
    for j, v in enumerate([0.01, 0.05, 0.10, 0.20, 0.30]):
        r = 7 + 34 * math.sqrt(v / vmax)
        x = x0 + 80 + j * 190
        body.append(f'<circle cx="{x}" cy="{y0 + 220}" r="{r:.1f}" fill="{colour(v, vmax)}" stroke="#6b7686" stroke-width="1.6"/>')
        body.append(text(x, y0 + 310, f"{v:.0%}", 34, weight=500, anchor="middle"))
    body.append(text(x0 + 20, y0 + 400, f"{DS_LABEL[dataset]}, all-features runs, mean over 50 replicates.", 30, weight=500, color="#4a5260"))
    body.append(text(x0 + 20, y0 + 445, f"Labels shown for nodes with at least {label_min:.0%} of credit.", 30, weight=500, color="#4a5260"))
    body.append(text(x0 + 20, y0 + 490, "Ringed node: the outcome.", 30, weight=500, color="#4a5260"))
    return document(183, vb_w, vb_h, "\n".join(body), {FAINT})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--layout", default=str(HERE / "fig2_layout.json"))
    ap.add_argument("--added-threshold", type=float, default=0.10)
    ap.add_argument("--fs", type=float, default=36, help="label size in layout units for the small panels")
    ap.add_argument("--dataset", default="full_set", help="data set for Figure 5")
    ap.add_argument("--label-min", type=float, default=0.03)
    a = ap.parse_args()
    layout = load_layout(a.layout)
    OUT.mkdir(parents=True, exist_ok=True)
    svg3, counts = fig3(layout, a.added_threshold, a.fs)
    (OUT / "fig3_recovery_grid.svg").write_text(svg3)
    (OUT / f"fig5_credit_dag_{a.dataset}.svg").write_text(fig5(layout, a.dataset, a.fs, a.label_min))
    print(json.dumps(counts, indent=1))


if __name__ == "__main__":
    main()
