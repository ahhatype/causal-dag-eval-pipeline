"""Node-label placement for the DAG figures: Inter Bold capitals, placed by a greedy search that avoids
other labels, nodes and edges and keeps each label nearest its own node."""

# Inter Bold advance widths, in em
ADV = {' ': 0.2368, 'A': 0.7466, 'B': 0.6616, 'C': 0.7397, 'D': 0.7222, 'E': 0.6074, 'F': 0.5869, 'G': 0.7505,
       'H': 0.7471, 'I': 0.2808, 'J': 0.5845, 'K': 0.7192, 'L': 0.5654, 'M': 0.9316, 'N': 0.7622, 'O': 0.7705,
       'P': 0.6479, 'Q': 0.7769, 'R': 0.6567, 'S': 0.6548, 'T': 0.6675, 'U': 0.7319, 'V': 0.7466, 'W': 1.0376,
       'X': 0.7383, 'Y': 0.731, 'Z': 0.6641, '0': 0.6743, '1': 0.4312, '2': 0.6299, '3': 0.6455, '4': 0.6763,
       '5': 0.6392, '6': 0.6494, '7': 0.5815, '8': 0.6509, '9': 0.6494, '-': 0.4678, '/': 0.3882, '+': 0.6787,
       '%': 1.0156, '.': 0.334, ',': 0.334, '(': 0.377, ')': 0.377, "'": 0.3389, '&': 0.6719}
CAP = 0.727
TRACK, LEAD, WEIGHT = 0.01, 1.18, 700

BASE = ["R", "L", "A", "B", "AR", "AL", "BR", "BL"]
PREF = {"R": 0, "L": 0, "A": 1.5, "B": 1.5, "AR": 2.5, "AL": 2.5, "BR": 2.5, "BL": 2.5}


def text_width(s, fs):
    return sum(ADV.get(c, 0.65) for c in s) * fs + TRACK * fs * (len(s) - 1)


def seg_len_in_box(p, q, box):
    x0, y0, x1, y1 = box
    t0, t1 = 0.0, 1.0
    dx, dy = q[0] - p[0], q[1] - p[1]
    for pp, qq in ((-dx, p[0] - x0), (dx, x1 - p[0]), (-dy, p[1] - y0), (dy, y1 - p[1])):
        if pp == 0:
            if qq < 0:
                return 0.0
            continue
        t = qq / pp
        if pp < 0:
            t0 = max(t0, t)
        else:
            t1 = min(t1, t)
        if t0 > t1:
            return 0.0
    return (t1 - t0) * (dx * dx + dy * dy) ** 0.5


def box_dist(b, p):
    dx = max(b[0] - p[0], 0, p[0] - b[2])
    dy = max(b[1] - p[1], 0, p[1] - b[3])
    return (dx * dx + dy * dy) ** 0.5


def overlap(a, b):
    w = min(a[2], b[2]) - max(a[0], b[0])
    h = min(a[3], b[3]) - max(a[1], b[1])
    return w * h if w > 0 and h > 0 else 0.0


def variants(label, first=None, max_w=None, fs=6.5):
    """Line-break options for a label written with | breaks; `first` is tried first."""
    words = label.replace("|", " ").upper().split()
    vs = [first] if first else []
    vs += [[t.strip().upper() for t in label.split("|")], [" ".join(words)], words]
    vs += [[" ".join(words[:i]), " ".join(words[i:])] for i in range(1, len(words))]
    out = []
    for v in vs:
        if v and v not in out and (max_w is None or max(text_width(t, fs) for t in v) < max_w):
            out.append(v)
    return out or [words]


def place(P, rad, labelled, options, segs, width, height, fs, fixed=None, gap=None, iters=12, preferred=None, clear=0.12):
    """Choose a position and line breaks for each labelled node.

    P: node -> (x, y) for every node drawn; rad: node -> radius; options: node -> list of line lists;
    segs: drawn edges as ((x, y), (x, y)); fixed: node -> ((position, dy), lines).
    Returns node -> dict(lines, anchor, lx, ly), with (lx, ly) the first baseline relative to the node.
    """
    fixed, preferred = fixed or {}, preferred or {}
    lh, g = fs * LEAD, gap if gap is not None else fs * 0.34
    k_off = fs / 6.5
    cands = [(c, dy) for c in BASE for dy in ((0, -3 * k_off, 3 * k_off) if c in ("R", "L") else (0,))]

    def sz(ls):
        return max(text_width(t, fs) for t in ls), CAP * fs + (len(ls) - 1) * lh

    def box(k, cand, ls):
        (cx, cy), r = P[k], rad[k]
        w, h = sz(ls)
        c, dy = cand
        cy += dy
        s = 2.5 * k_off
        return {"R": (cx + r + g, cy - h / 2, cx + r + g + w, cy + h / 2),
                "L": (cx - r - g - w, cy - h / 2, cx - r - g, cy + h / 2),
                "A": (cx - w / 2, cy - r - g - h, cx + w / 2, cy - r - g),
                "B": (cx - w / 2, cy + r + g, cx + w / 2, cy + r + g + h),
                "AR": (cx + s, cy - r - g - h, cx + s + w, cy - r - g),
                "AL": (cx - s - w, cy - r - g - h, cx - s, cy - r - g),
                "BR": (cx + s, cy + r + g, cx + s + w, cy + r + g + h),
                "BL": (cx - s - w, cy + r + g, cx - s, cy + r + g + h)}[c]

    def pref(cand, ls):
        return PREF[cand[0]] + (1.5 if cand[1] else 0) + 1.2 * (len(ls) - 1) + (3 if len(ls) > 2 else 0)

    def cost(k, cand, ls, chosen):
        b = box(k, cand, ls)
        e = 0.8 * k_off
        pad = (b[0] - e, b[1] - e, b[2] + e, b[3] + e)
        c = pref(cand, ls)
        if preferred.get(k) == cand[0] and cand[1] == 0 and ls == options[k][0]:
            c -= 6
        c += 400 * (max(0, -b[0] + 2) + max(0, b[2] - width + 2) + max(0, -b[1] + 2) + max(0, b[3] - height + 2))
        q = clear * fs
        lpad = (b[0] - q, b[1] - q, b[2] + q, b[3] + q)
        for j, (cj, lj) in chosen.items():
            if j != k:
                c += 30 * overlap(lpad, box(j, cj, lj)) / k_off ** 2
        for j, (x, y) in P.items():
            rr = rad[j] + e
            c += 40 * overlap(pad, (x - rr, y - rr, x + rr, y + rr)) / k_off ** 2
        for p, q in segs:
            c += 2.5 * seg_len_in_box(p, q, (b[0] - 2 * e, b[1] - 1.5 * e, b[2] + 2 * e, b[3] + 1.5 * e)) / k_off
        d0 = box_dist(b, P[k])
        for j, pj in P.items():
            if j != k:
                c += 25 * max(0, d0 + 2.5 * k_off - box_dist(b, pj)) / k_off
        return c

    chosen = {k: fixed[k] if k in fixed else (("R", 0), options[k][0]) for k in labelled}
    free = sorted((k for k in labelled if k not in fixed), key=lambda k: -sz(options[k][0])[0])
    for _ in range(iters):
        changed = False
        for k in free:
            best = min(((cd, ls) for cd in cands for ls in options[k]), key=lambda t: cost(k, t[0], t[1], chosen))
            if best != chosen[k]:
                chosen[k], changed = best, True
        if not changed:
            break
    out = {}
    for k, (cand, ls) in chosen.items():
        b = box(k, cand, ls)
        x, y = P[k]
        anchor, lx = "start", b[0] - x
        if cand[0] in ("A", "B"):
            anchor, lx = "middle", 0.0
        elif cand[0] in ("L", "AL", "BL"):
            anchor, lx = "end", b[2] - x
        out[k] = {"lines": ls, "anchor": anchor, "lx": round(lx, 2), "ly": round(b[1] + CAP * fs - y, 2), "place": cand[0],
                  "residual": round(cost(k, cand, ls, chosen) - pref(cand, ls), 1)}
    return out
