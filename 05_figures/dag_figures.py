"""Figures drawn on the portrait DAG layout (fig2_portrait.json, in points).

Figure 2: the true analysis graph, 4 in wide.
Figure 3: where structure recovery failed in the full-set, with a key; the reference-set and astronaut-set
          go in a supplementary figure.
Figure 6 (and Supplementary Figs 1-2): where each attribution method put its credit, one panel per truth or method.

All write SVG with live text; PDFs are made from the SVG by a browser so the label font is embedded.
"""
import argparse
import json
import math
import os
from pathlib import Path
from xml.sax.saxutils import escape

import pandas as pd

import dag_labels as dl

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("PIPELINE_ROOT", HERE.parent))
OUT = ROOT / "outputs" / "05_figures"
STRUCTURE = ROOT / "outputs" / "03_structure_validation"
ATTR = ROOT / "outputs" / "04_attribution_validation"
DATA = ROOT / "02_data"

FONT = "Inter, Helvetica Neue, Arial, sans-serif"  # OFL, Google Fonts
INK, GREY, BLUE, FAINT = "#0b0d12", "#4a5260", "#4cc0e8", "#c4ccd6"
BLUE_EDGE = "#a6e0f4"  # node blue 50% lighter, for edges recovered correctly
DATASETS = ["full_set", "reference_subsample", "astronaut_set"]
DS_LABEL = {"full_set": "Full-set", "reference_subsample": "Reference-set", "astronaut_set": "Astronaut-set"}
DS_N = {"full_set": "n = 10,000", "reference_subsample": "n = 900", "astronaut_set": "n = 900, selected"}
OUTCOME = "nephrolithiasis"

TEXT_W = 468.0  # printed width of Figures 3 and 6, pt (6.5 in)
PANEL_GAP, TITLE_H = 12.0, 16.0
PANEL_PAD = 14.0  # layout units added left and right of each panel for labels on the outer lanes
LABEL_PT = 5.0  # printed size of node labels in the panels

# recovery categories: colour, dash pattern (panel units), arrow, legend label
CATS = {
    "absent": ("#f29a92", "0.1 2.6", False, "Missed"),
    "reversed": ("#6fc2a7", None, True, "Recovered, direction reversed"),
    "overoriented": ("#ecb740", "5 1.8 1.2 1.8", True, "Direction imposed on an undirected edge"),
    "unoriented": ("#b9a3e3", None, False, "Left undirected"),
    "added": ("#cc79a7", "1 2.2", False, "Added (not in the true graph)"),
}

METHODS = [("truth", "True effect"), ("standard_shap", "Standard SHAP"), ("ordering_only_shap", "Ordering-only SHAP"),
           ("ng_style_shap", "Causal predictive SHAP"), ("interventional_shap", "Causal interventional SHAP")]


def unit(x, y):
    l = math.hypot(x, y) or 1
    return x / l, y / l


def marker(c, size=7):
    return (f'<marker id="arr{c.strip("#")}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="{size}" '
            f'markerHeight="{size}" markerUnits="strokeWidth" orient="auto"><path d="M0,1 L10,5 L0,9 z" fill="{c}"/></marker>')


def text(x, y, s, size, weight=700, anchor="start", color=INK):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" font-weight="{weight}" '
            f'text-anchor="{anchor}" fill="{color}">{escape(s)}</text>')


def label_svg(nid, x, y, lines, anchor, fs, color=INK):
    lh = fs * dl.LEAD
    spans = "".join(f'<tspan x="{x:.2f}" dy="{0 if i == 0 else lh:.2f}">{escape(t)}</tspan>' for i, t in enumerate(lines))
    return (f'<text id="label_{nid}" x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" font-family="{FONT}" '
            f'font-size="{fs:.2f}" font-weight="{dl.WEIGHT}" letter-spacing="{dl.TRACK}em" fill="{color}">{spans}</text>')


def svg_doc(w, h, body, colors, marker_size=7):
    defs = "<defs>" + "".join(marker(c, marker_size) for c in sorted(colors)) + "</defs>"
    return (f'<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w:.1f} {h:.1f}" '
            f'width="{w / 72:.3f}in" height="{h / 72:.3f}in">{defs}<rect width="{w:.1f}" height="{h:.1f}" fill="#ffffff"/>'
            + "\n".join(body) + "</svg>")


def load_portrait(path):
    return json.loads(Path(path).read_text())


def read_edges():
    return list(pd.read_csv(STRUCTURE / "analysis_graph_edges.csv")[["parent", "child"]].itertuples(index=False, name=None))


def ancestors_of(target, edges):
    par = {}
    for a, b in edges:
        par.setdefault(b, set()).add(a)
    seen, stack = set(), [target]
    while stack:
        for p in par.get(stack.pop(), ()):
            if p not in seen:
                seen.add(p)
                stack.append(p)
    return seen


class Panel:
    """One DAG panel in the portrait layout's units, placed in a figure by translate and scale."""

    def __init__(self, lay, scale, label_pt=LABEL_PT):
        self.lay, self.nodes, self.s, self.label_pt = lay, lay["nodes"], scale, label_pt
        self.P = {k: (n["x"] + PANEL_PAD, n["y"]) for k, n in self.nodes.items()}
        self.w = lay["width"] + 2 * PANEL_PAD
        self.rad, self.segs, self.parts, self.labels = {}, [], [], {}

    def edge(self, a, b, color, width, dash=None, arrow=True, reverse=False):
        if reverse:
            a, b = b, a
        (ax, ay), (bx, by) = self.P[a], self.P[b]
        ux, uy = unit(bx - ax, by - ay)
        ra, rb = self.rad.get(a, 2.6), self.rad.get(b, 2.6)
        sx, sy = ax + ux * (ra + 0.8), ay + uy * (ra + 0.8)
        ex, ey = bx - ux * (rb + (1.4 if arrow else 0.8)), by - uy * (rb + (1.4 if arrow else 0.8))
        attrs = f'd="M{sx:.2f},{sy:.2f} L{ex:.2f},{ey:.2f}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linecap="round"'
        if dash:
            attrs += f' stroke-dasharray="{dash}"'
        if arrow:
            attrs += f' marker-end="url(#arr{color.strip("#")})"'
        self.parts.append(f"<path {attrs}/>")
        self.segs.append(((ax, ay), (bx, by)))
        return color if arrow else None

    def node(self, nid, fill, stroke, r, ring=False, width=0.9):
        x, y = self.P[nid]
        self.parts.append(f'<circle cx="{x}" cy="{y}" r="{r:.2f}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>')
        if ring:
            self.parts.append(f'<circle cx="{x}" cy="{y}" r="{r + 1.8:.2f}" fill="none" stroke="{stroke}" stroke-width="{width * 0.6:.2f}"/>')

    def place_labels(self, labelled):
        fs = self.label_pt / self.s
        opts = {k: dl.variants(self.nodes[k]["label"], self.nodes[k]["label"].split("|")) for k in labelled}
        pref = {k: self.nodes[k].get("place") for k in labelled}
        lab = dl.place(self.P, self.rad, labelled, opts, self.segs, self.w, self.lay["height"], fs,
                       iters=10, preferred=pref, clear=0.45)
        for k, a in lab.items():
            x, y = self.P[k]
            self.parts.append(label_svg(k, x + a["lx"], y + a["ly"], a["lines"], a["anchor"], fs))

    def svg(self, x0, y0):
        return f'<g transform="translate({x0:.2f},{y0:.2f}) scale({self.s:.4f})">' + "".join(self.parts) + "</g>"


def panel_title(x0, y0, letter, title, sub=None):
    tail = f'<tspan font-size="6" font-weight="500" fill="{GREY}">  {escape(sub)}</tspan>' if sub else ""
    return (text(x0, y0 + 9, letter, 8.5) + f'<text x="{x0 + 10:.1f}" y="{y0 + 9:.1f}" font-family="{FONT}" font-size="7" '
            f'font-weight="700" fill="{INK}" xml:space="preserve">{escape(title)}{tail}</text>')


# ---------- Figure 2: the true analysis graph ----------

def fig2(lay, edges):
    nodes, f = lay["nodes"], lay["font"]
    anc = ancestors_of(OUTCOME, edges)
    body = []
    for a, b in edges:
        A, B = nodes[a], nodes[b]
        ux, uy = unit(B["x"] - A["x"], B["y"] - A["y"])
        ra, rb = (lay["outcome_r"] if n == OUTCOME else lay["node_r"] for n in (a, b))
        sx, sy = A["x"] + ux * (ra + 0.8), A["y"] + uy * (ra + 0.8)
        ex, ey = B["x"] - ux * (rb + 1.2), B["y"] - uy * (rb + 1.2)
        body.append(f'<path id="edge_{a}__{b}" d="M{sx:.2f},{sy:.2f} L{ex:.2f},{ey:.2f}" stroke="{BLUE}" '
                    f'stroke-width="0.55" fill="none" marker-end="url(#arr{BLUE.strip("#")})"/>')
    for nid, n in nodes.items():
        r = lay["outcome_r"] if nid == OUTCOME else lay["node_r"]
        filled = nid == OUTCOME or nid in anc
        g = [f'<circle cx="{n["x"]}" cy="{n["y"]}" r="{r}" fill="{BLUE if filled else "#ffffff"}" stroke="{BLUE}" stroke-width="0.8"/>']
        if nid == OUTCOME:
            g.append(f'<circle cx="{n["x"]}" cy="{n["y"]}" r="{r + 1.6}" fill="none" stroke="{BLUE}" stroke-width="0.5"/>')
        g.append(label_svg(nid, n["x"] + n["lx"], n["y"] + n["ly"], n["label"].split("|"), n["anchor"], f["size"]))
        body.append(f'<g id="node_{nid}">{"".join(g)}</g>')
    return svg_doc(lay["width"], lay["height"], body, {BLUE})


# ---------- Figure 3: recovery failures per data set ----------

def recovery(added_threshold):
    d = pd.read_csv(STRUCTURE / "difference_dag_edges.csv")
    true = d[d["true_edge"].astype(str).str.upper() == "TRUE"].copy()
    rec = pd.read_csv(STRUCTURE / "edge_recovery.csv")
    rec = rec[(rec["variant"] == "MGM-PC-Stable") & (rec["true_edge"].astype(str).str.upper() != "TRUE")]
    added = rec[rec["adjacency_rate"] >= added_threshold].copy()
    added[["a", "b"]] = added["edge"].str.split(" ", expand=True)
    return true, added


KEY_ROWS = [("absent", "Missed"), ("reversed", "Reversed"), ("overoriented", "Direction imposed"),
            ("unoriented", "Left undirected"), ("added", "Added")]


def recovery_panel(lay, s, ds, true, added, label_pt):
    """One data set's recovery panel; returns the panel, arrow colours used and category counts."""
    p = Panel(lay, s, label_pt=label_pt)
    nr = 2.6 + 1.4 * (1 - min(s, 1))
    for nid in p.nodes:
        p.rad[nid] = nr * 1.4 if nid == OUTCOME else nr
    t = true[true["dataset"] == ds]
    hit, cnt, colors = set(), {}, set()
    for e in t.itertuples():
        if e.modal_category == "correct":
            p.edge(e.parent, e.child, BLUE_EDGE, 0.9, arrow=False)
    for e in t.itertuples():
        cat = e.modal_category
        if cat == "correct":
            continue
        col, dash, arrow, _ = CATS[cat]
        arrow = arrow and e.true_mark != "undir" or cat == "overoriented"
        m = p.edge(e.parent, e.child, col, 1.1 if cat == "absent" else 1.3, dash=dash, arrow=arrow, reverse=(cat == "reversed"))
        if m:
            colors.add(m)
        hit |= {e.parent, e.child}
        cnt[cat] = cnt.get(cat, 0) + 1
    for e in added[added["dataset"] == ds].itertuples():
        col, dash, _, _ = CATS["added"]
        p.edge(e.a, e.b, col, 1.7, dash=dash, arrow=False)
        hit |= {e.a, e.b}
        cnt["added"] = cnt.get("added", 0) + 1
    for nid in p.nodes:
        if nid == OUTCOME:
            p.node(nid, BLUE, BLUE, p.rad[nid], ring=True)
        else:
            on = nid in hit
            p.node(nid, BLUE if on else "#ffffff", BLUE if on else FAINT, nr, width=0.9)
    p.place_labels(sorted(hit | {OUTCOME}))
    return p, colors, cnt


def recovery_key(x0, y0, datasets, counts, wide=False):
    """Line key plus a table of edge counts; wide puts the table to the right of the lines."""
    body, colors = [text(x0, y0 + 9, "Key", 7)], set()
    yy = y0 + 26
    for cat, _ in KEY_ROWS:
        col, dash, arrow, lab = CATS[cat]
        da = f' stroke-dasharray="{" ".join(str(round(float(v) * 0.6, 2)) for v in dash.split())}"' if dash else ""
        mk = f' marker-end="url(#arr{col.strip("#")})"' if arrow else ""
        if arrow:
            colors.add(col)
        body.append(f'<path d="M{x0},{yy - 2.2} L{x0 + 22},{yy - 2.2}" stroke="{col}" stroke-width="0.9" '
                    f'stroke-linecap="round" fill="none"{da}{mk}/>')
        body.append(text(x0 + 28, yy, lab, 6, weight=500))
        yy += 11
    body.append(f'<path d="M{x0},{yy - 2.2} L{x0 + 22},{yy - 2.2}" stroke="{BLUE_EDGE}" stroke-width="0.6" fill="none"/>')
    body.append(text(x0 + 28, yy, "Recovered correctly", 6, weight=500))
    lines_bottom = yy
    tx, ty = (x0 + 200, y0 + 26 - 11) if wide else (x0, yy + 26)
    short = {"full_set": "Full", "reference_subsample": "Ref.", "astronaut_set": "Astro."}
    cx = [tx] + [tx + 98 + 36 * j for j in range(len(datasets))]
    for x, h in zip(cx, ["Edges"] + [short[d] for d in datasets]):
        body.append(text(x, ty, h, 6, anchor="start" if x == cx[0] else "end"))
    for cat, lab in KEY_ROWS:
        ty += 10
        body.append(text(cx[0], ty, lab, 6, weight=500))
        for x, ds in zip(cx[1:], datasets):
            body.append(text(x, ty, str(counts[ds].get(cat, 0)), 6, weight=500, anchor="end"))
    body.append(text(cx[0], ty + 14, "By most frequent outcome over", 5.5, weight=500, color=GREY))
    body.append(text(cx[0], ty + 21.5, "50 replicates, of 63 true edges.", 5.5, weight=500, color=GREY))
    return body, colors, max(lines_bottom, ty + 21.5) - y0


def fig3(lay, added_threshold):
    """Main text: the full-set panel with the key to its right, 6.5 in wide."""
    true, added = recovery(added_threshold)
    pw = 312.0
    s = pw / (lay["width"] + 2 * PANEL_PAD)
    ds = "full_set"
    p, colors, cnt = recovery_panel(lay, s, ds, true, added, label_pt=6.0)
    ph = lay["height"] * s
    body = [f'<text x="0" y="9" font-family="{FONT}" font-size="7" font-weight="700" fill="{INK}" xml:space="preserve">'
            f'{DS_LABEL[ds]}<tspan font-size="6" font-weight="500" fill="{GREY}">  {DS_N[ds]}</tspan></text>', p.svg(0, TITLE_H)]
    kb, kc, _ = recovery_key(pw + 2, TITLE_H + 6, [ds], {ds: cnt})
    body += kb
    return svg_doc(TEXT_W, TITLE_H + ph, body, colors | kc | {FAINT}, marker_size=5), {ds: cnt}


def fig3_supp(lay, added_threshold):
    """Supplement: reference-set and astronaut-set panels side by side, 6.5 in wide, key beneath."""
    true, added = recovery(added_threshold)
    pw = (TEXT_W - PANEL_GAP) / 2
    s = pw / (lay["width"] + 2 * PANEL_PAD)
    ph = lay["height"] * s
    body, colors, counts = [], set(), {}
    sets = ["reference_subsample", "astronaut_set"]
    for i, ds in enumerate(sets):
        p, c, cnt = recovery_panel(lay, s, ds, true, added, label_pt=5.5)
        colors |= c
        counts[ds] = cnt
        x0 = i * (pw + PANEL_GAP)
        body.append(panel_title(x0, 0, "ab"[i], DS_LABEL[ds], DS_N[ds]))
        body.append(p.svg(x0, TITLE_H))
    y0 = TITLE_H + ph + 10
    kb, kc, kh = recovery_key(0, y0, sets, counts, wide=True)
    body += kb
    return svg_doc(TEXT_W, y0 + kh + 6, body, colors | kc | {FAINT}, marker_size=5), counts


# ---------- Figure 6: where credit went ----------

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
    # white to deep blue through the Figure 2 blue, by square root so small shares stay visible
    f = min(1.0, math.sqrt(max(v, 0) / vmax)) if vmax else 0
    lo, mid, hi = (255, 255, 255), (76, 192, 232), (6, 54, 122)
    a, b, g = (lo, mid, f / 0.5) if f < 0.5 else (mid, hi, (f - 0.5) / 0.5)
    return "#%02x%02x%02x" % tuple(round(a[i] + (b[i] - a[i]) * g) for i in range(3))


def credit_r(v, vmax):
    return min(2.4 + 9.0 * math.sqrt(v / vmax), 11.5) if v > 0 else 2.0


def fig6(lay, dataset, label_min, edges):
    sh = credit_shares(dataset)
    vmax = 0.30
    pw = (TEXT_W - 2 * PANEL_GAP) / 3
    s = pw / (lay["width"] + 2 * PANEL_PAD)
    ph = lay["height"] * s
    row_h = TITLE_H + ph + PANEL_GAP
    body = []
    for i, (m, title) in enumerate(METHODS):
        p = Panel(lay, s)
        vals = {nid: (float(sh.loc[nid, m]) if nid in sh.index else 0.0) for nid in p.nodes}
        for nid, v in vals.items():
            p.rad[nid] = credit_r(v, vmax)
        for a, b in edges:
            p.edge(a, b, FAINT, 0.9, arrow=False)
        for nid, v in vals.items():
            stroke = INK if nid == OUTCOME else "#6b7686"
            p.node(nid, colour(v, vmax) if nid != OUTCOME else "#ffffff", stroke, p.rad[nid], ring=(nid == OUTCOME), width=0.8)
        p.place_labels(sorted(k for k, v in vals.items() if v >= label_min or k == OUTCOME))
        x0, y0 = (i % 3) * (pw + PANEL_GAP), (i // 3) * row_h
        body.append(panel_title(x0, y0, "abcde"[i], title))
        body.append(p.svg(x0, y0 + TITLE_H))
    # key in the sixth slot
    x0, y0 = 2 * (pw + PANEL_GAP), row_h
    body.append(text(x0, y0 + 9, "Share of credit", 7))
    for j, v in enumerate([0.01, 0.05, 0.10, 0.20, 0.30]):
        r = credit_r(v, vmax) * s
        x = x0 + 10 + j * 30
        body.append(f'<circle cx="{x:.1f}" cy="{y0 + 40:.1f}" r="{r:.2f}" fill="{colour(v, vmax)}" stroke="#6b7686" stroke-width="{0.8 * s:.2f}"/>')
        body.append(text(x, y0 + 58, f"{v:.0%}", 6, weight=500, anchor="middle"))
    notes = [f"{DS_LABEL[dataset]} ({DS_N[dataset]}), all-features runs,", "mean over 50 replicates.",
             f"Labels: nodes with at least {label_min:.0%} of credit.", "Ringed node: the outcome."]
    for j, t in enumerate(notes):
        body.append(text(x0, y0 + 80 + j * 9, t, 5.5, weight=500, color=GREY))
    return svg_doc(TEXT_W, 2 * row_h - PANEL_GAP, body, {FAINT})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--portrait", default=str(HERE / "fig2_portrait.json"))
    ap.add_argument("--added-threshold", type=float, default=0.10)
    ap.add_argument("--dataset", default="all", help="data set for Figure 6, or all (the other two are supplementary)")
    ap.add_argument("--label-min", type=float, default=0.03)
    ap.add_argument("--figs", default="2,3,6", help="which figures to draw")
    a = ap.parse_args()
    figs = set(a.figs.split(","))
    OUT.mkdir(parents=True, exist_ok=True)
    lay, edges = load_portrait(a.portrait), read_edges()
    if "2" in figs:
        (OUT / "fig2_true_graph.svg").write_text(fig2(lay, edges))
    if "3" in figs:
        svg3, counts = fig3(lay, a.added_threshold)
        (OUT / "fig3_recovery_grid.svg").write_text(svg3)
        svg3s, counts_s = fig3_supp(lay, a.added_threshold)
        (OUT / "supp_fig_recovery_grid.svg").write_text(svg3s)
        print(json.dumps({**counts, **counts_s}, indent=1))
    if "6" in figs:
        for ds in (DATASETS if a.dataset == "all" else [a.dataset]):
            (OUT / f"fig5_credit_dag_{ds}.svg").write_text(fig6(lay, ds, a.label_min, edges))


if __name__ == "__main__":
    main()
