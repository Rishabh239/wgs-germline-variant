#!/usr/bin/env python3
"""
check_samplesheet.py - validate and normalize the input samplesheet.

Expected columns (header required): sample,fastq_1,fastq_2
  - sample    : unique sample identifier (used as read-group SM and file prefix)
  - fastq_1   : path/URL to R1 (.fastq.gz / .fq.gz)
  - fastq_2   : path/URL to R2 (.fastq.gz / .fq.gz)

Writes a cleaned samplesheet on success; exits non-zero with a clear message
on any validation error so the pipeline fails fast at launch.
"""
import argparse
import csv
import sys

REQUIRED = ["sample", "fastq_1", "fastq_2"]
FASTQ_EXT = (".fastq.gz", ".fq.gz")


def die(msg):
    sys.stderr.write(f"[check_samplesheet] ERROR: {msg}\n")
    sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description="Validate the input samplesheet.")
    ap.add_argument("samplesheet", help="Input CSV")
    ap.add_argument("out", nargs="?", default="samplesheet.valid.csv",
                    help="Output normalized CSV")
    args = ap.parse_args()

    with open(args.samplesheet, newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            die("samplesheet is empty")
        header = [h.strip() for h in reader.fieldnames]
        missing = [c for c in REQUIRED if c not in header]
        if missing:
            die(f"missing required column(s): {', '.join(missing)}; found {header}")

        rows = []
        seen = set()
        for i, row in enumerate(reader, start=2):  # line 2 = first data row
            sample = (row.get("sample") or "").strip()
            r1 = (row.get("fastq_1") or "").strip()
            r2 = (row.get("fastq_2") or "").strip()
            if not sample:
                die(f"line {i}: empty 'sample'")
            if " " in sample:
                die(f"line {i}: 'sample' may not contain spaces: '{sample}'")
            if sample in seen:
                die(f"line {i}: duplicate sample id '{sample}'")
            seen.add(sample)
            for label, val in (("fastq_1", r1), ("fastq_2", r2)):
                if not val:
                    die(f"line {i}: empty '{label}'")
                if not val.endswith(FASTQ_EXT):
                    die(f"line {i}: '{label}' must end with one of {FASTQ_EXT}: '{val}'")
            if r1 == r2:
                die(f"line {i}: fastq_1 and fastq_2 are identical")
            rows.append({"sample": sample, "fastq_1": r1, "fastq_2": r2})

    if not rows:
        die("no data rows found")

    with open(args.out, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=REQUIRED)
        writer.writeheader()
        writer.writerows(rows)

    sys.stderr.write(f"[check_samplesheet] OK: {len(rows)} sample(s) validated -> {args.out}\n")


if __name__ == "__main__":
    main()
