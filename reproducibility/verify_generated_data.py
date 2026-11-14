"""Usage: python3 reproducibility/verify_generated_data.py [--data-dir 02_data] [--json report.json]

Compares the SHA-256 of every generated replicate CSV with expected-data-sha256.csv (50 replicates, three data sets,
observed and latent files: 300 files). The data are only read. The exit status is zero only if every file matches.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=Path("02_data"))
    ap.add_argument("--manifest", type=Path, default=Path(__file__).with_name("expected-data-sha256.csv"))
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()
    data_dir = args.data_dir.resolve()
    with args.manifest.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    failures = []
    for row in rows:
        p = (data_dir / row["file"]).resolve()
        if not p.is_relative_to(data_dir):
            raise ValueError("manifest path outside the data directory")
        if not p.is_file():
            failures.append({"file": row["file"], "reason": "missing"})
        elif hashlib.sha256(p.read_bytes()).hexdigest() != row["sha256"]:
            failures.append({"file": row["file"], "reason": "SHA-256 mismatch"})
    report = {"expected_files": len(rows), "matching_files": len(rows) - len(failures), "passed": not failures,
              "failures": failures}
    output = json.dumps(report, indent=2) + "\n"
    if args.json:
        args.json.write_text(output, encoding="utf-8")
    print(output, end="")
    raise SystemExit(0 if report["passed"] else 1)


if __name__ == "__main__":
    main()
