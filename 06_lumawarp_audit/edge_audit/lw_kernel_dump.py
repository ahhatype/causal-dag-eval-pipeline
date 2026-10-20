#!/usr/bin/env python3
"""Read saved LumaWarp models (.lsmdl, container format 1) and write each kernel axis's
learned slice k = A^T A with its spectrum. Pure Python + numpy: no engine build needed.

Layout of the KERNEL section (worked out from the files, checked against the theta values
in the diagnostics logs):
  f64 amplitude | u64 n_axes | per axis: u8 space | string name (u64 len + bytes) |
  u32 dim | u32 rank | f32 vector A (u64 count = rank*dim) | f32 vector mean_phi (dim) |
  f32 theta | i32 source | then component weights / component_of_axis vectors.

Usage: lw_kernel_dump.py out_prefix model1.lsmdl [model2.lsmdl ...]
Writes <prefix>_axes.csv (one row per model, member, axis) and <prefix>_slices.npz
(full slices, key 'model|member|axis').
"""
import csv, struct, sys
import numpy as np
from lsmdl_container import sections

def read_kernel(k):
    off = 0
    amp, n = struct.unpack_from("<dQ", k, off); off += 16
    axes = []
    for _ in range(n):
        space = k[off]; off += 1
        (ln,) = struct.unpack_from("<Q", k, off); off += 8
        name = k[off:off+ln].decode(); off += ln
        dim, rank, na = struct.unpack_from("<IIQ", k, off); off += 16
        assert na == dim * rank, (name, dim, rank, na)
        A = np.frombuffer(k, "<f4", na, off).reshape(rank, dim).astype(float); off += 4 * na
        (nm,) = struct.unpack_from("<Q", k, off); off += 8
        mean_phi = np.frombuffer(k, "<f4", nm, off).astype(float); off += 4 * nm
        theta, source = struct.unpack_from("<fi", k, off); off += 8
        axes.append(dict(name=name, space=space, dim=dim, rank=rank, A=A, mean_phi=mean_phi,
                         theta=theta, source=source))
    return amp, axes

def members(path):
    buf = open(path, "rb").read()
    _, secs = sections(buf)
    tags = [t for t, _ in secs]
    if "EMETA" in tags:
        for t, p in secs:
            if t.startswith("MEMBER"):
                yield int(t[6:]), dict(sections(p)[1])
    else:
        yield 0, dict(secs)

def main():
    prefix, paths = sys.argv[1], sys.argv[2:]
    fa = open(prefix + "_axes.csv", "w", newline=""); wa = csv.writer(fa)
    slices = {}  # full slices go to one compressed .npz (the CSV form ran to ~700 MB)
    wa.writerow("model member axis_index axis space source dim rank theta amplitude trace effective_modes "
                "share_mode0 degree_centroid share_deg0 share_deg1 share_deg2 share_deg3plus "
                "eigenvalues_desc diagonal".split())
    for path in paths:
        try:
            for m, secs in members(path):
                if "KERNEL" not in secs:
                    continue
                amp, axes = read_kernel(secs["KERNEL"])
                for i, ax in enumerate(axes):
                    K = ax["A"].T @ ax["A"]
                    w = np.sort(np.clip(np.linalg.eigvalsh(K), 0, None))[::-1]
                    tr = w.sum()
                    d = np.diag(K)
                    ds = d / d.sum() if d.sum() > 0 else d
                    eff = tr * tr / (w * w).sum() if tr > 0 else 0.0
                    wa.writerow([path, m, i, ax["name"], ax["space"], ax["source"], ax["dim"], ax["rank"],
                                 "%.9g" % ax["theta"], "%.9g" % amp, "%.9g" % tr, "%.6g" % eff,
                                 "%.6g" % (w[0] / tr if tr > 0 else 0),
                                 "%.6g" % float((np.arange(len(ds)) * ds).sum()),
                                 "%.6g" % ds[0], "%.6g" % (ds[1] if len(ds) > 1 else 0),
                                 "%.6g" % (ds[2] if len(ds) > 2 else 0), "%.6g" % ds[3:].sum(),
                                 ";".join("%.6g" % x for x in w), ";".join("%.6g" % x for x in d)])
                    slices["%s|%d|%s" % (path, m, ax["name"])] = K.astype("f4")
        except Exception as e:
            print("SKIP", path, repr(e), file=sys.stderr)
    fa.close()
    np.savez_compressed(prefix + "_slices.npz", **slices)

if __name__ == "__main__":
    main()
