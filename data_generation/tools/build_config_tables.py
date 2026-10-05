#!/usr/bin/env python3
"""Build config/nodes.csv and config/edges.csv from Supplementary Tables 1 and 2.

Usage: python3 data_generation/tools/build_config_tables.py <appendix.md> <out_dir>
"""
import csv
import re
import sys

# Edges whose sign is unspecified in the source table. PROVISIONAL.
PROVISIONAL_SIGNS = {
    ("evacuation", "loss of crew life"): 0.8,         # evacuation risk (vs rescue)
    ("percutaneous nephrostomy", "medical illness"): -0.4,  # relief (vs complications)
}


def slug(s):
    s = s.lower().replace("k+", "k")
    return re.sub(r"[^a-z0-9]+", "_", s).strip("_")


def table_after(md, heading):
    lines = md[md.index(heading):].split("\n")[1:]
    rows, started = [], False
    for line in lines:
        if line.startswith("|"):
            started = True
            rows.append([c.strip() for c in line.strip().strip("|").split("|")])
        elif started:
            break
    return rows[2:]


def node_rows(md):
    for name, level, scale, observed, role, _parents, s4 in table_after(md, "## Supplementary Table 1"):
        if "Via Mission Era" in s4:
            feature_set = "context"
        elif "Ancestor set: Yes" in s4:
            feature_set = "ancestor"
        elif "all-features set: Yes" in s4:
            feature_set = "non_ancestor"
        else:
            feature_set = "none"
        if level.startswith("Context"):
            lv = "era_fixed"
        elif level.startswith("Era"):
            lv = "era_root"
        else:
            lv = level.lower()
        yield [slug(name), name, lv, "binary" if scale == "Binary" else "continuous",
               "TRUE" if observed.startswith("Yes") else "FALSE", role, feature_set,
               "FALSE" if level.startswith("Sampling") else "TRUE"]


def edge_rows(md):
    for parent, child, etype, cscale, _tier, coef, _anc, _notes in table_after(md, "## Supplementary Table 2"):
        value, status = "", "given"
        if (parent, child) in PROVISIONAL_SIGNS:
            value, status = PROVISIONAL_SIGNS[(parent, child)], "provisional_sign"
        elif re.fullmatch(r"-?\d*\.?\d+", coef):
            value = coef
        elif coef == "era-fixed":
            status = "era_fixed"
        elif coef.startswith("see Table 3"):
            status = "era_table"
        elif coef == "5% flip":
            status = "flip"
        elif coef == "[blank]":
            status = "blank"
        else:
            sys.exit(f"unrecognised coefficient for {parent} -> {child}: {coef!r}")
        yield [slug(parent), slug(child), etype,
               "binary" if cscale.startswith("Binary") else "continuous", value, status]


def write(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main(appendix, out_dir):
    md = open(appendix, encoding="utf-8").read()
    write(f"{out_dir}/nodes.csv",
          ["id", "name", "level", "scale", "observed", "role", "feature_set", "simulate"], node_rows(md))
    write(f"{out_dir}/edges.csv",
          ["parent", "child", "edge_type", "child_scale", "coef", "coef_status"], edge_rows(md))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
