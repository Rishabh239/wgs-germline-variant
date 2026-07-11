#!/usr/bin/env python3
"""
vcf_summary.py - lightweight, dependency-free VCF summary for QC reporting.

Reads a (optionally gzipped) VCF and reports:
  - total / PASS variant counts
  - SNP / INDEL / MNP / other counts
  - Ti/Tv ratio (biallelic SNPs)
  - het / hom-alt genotype counts and het:hom ratio (across all samples)
  - per-chromosome variant counts
  - mean QUAL

Outputs:
  <prefix>.summary.tsv              human-readable metric<TAB>value table
  <prefix>_vcf_summary_mqc.json     MultiQC custom-content (table) file

Intentionally uses only the Python standard library so it runs anywhere
(no pysam/cyvcf2), which keeps the container small and the step reproducible.
"""
import argparse
import gzip
import json
import sys
from collections import Counter

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
BASES = {"A", "C", "G", "T"}


def _open(path):
    """Open plain or gzipped text transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def classify(ref, alt):
    """Return 'SNP', 'INDEL', 'MNP', or 'OTHER' for a single REF/ALT pair."""
    if alt in (".", "*", "") or "<" in alt or "[" in alt or "]" in alt:
        return "OTHER"  # symbolic / breakend / spanning-deletion
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(ref) == len(alt):
        return "MNP"
    return "INDEL"


def parse_gt(sample_field):
    """Return the genotype string (first FORMAT sub-field) or None."""
    if sample_field in (".", "./.", ".|.", ""):
        return None
    gt = sample_field.split(":", 1)[0]
    return gt if gt not in (".", "./.", ".|.") else None


def main():
    ap = argparse.ArgumentParser(description="Dependency-free VCF summary.")
    ap.add_argument("--vcf", required=True, help="Input VCF (.vcf or .vcf.gz)")
    ap.add_argument("--prefix", default="cohort", help="Output file prefix")
    ap.add_argument("--outdir", default=".", help="Output directory")
    args = ap.parse_args()

    total = 0
    passed = 0
    typ = Counter()
    per_chrom = Counter()
    ti = tv = 0
    het = hom_alt = 0
    qual_sum = 0.0
    qual_n = 0
    n_samples = 0

    with _open(args.vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                n_samples = max(0, len(cols) - 9)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, _pos, _id, ref, alt, qual, flt = f[0], f[1], f[2], f[3], f[4], f[5], f[6]
            total += 1
            per_chrom[chrom] += 1
            if flt in ("PASS", ".", ""):
                passed += 1
            if qual not in (".", ""):
                try:
                    qual_sum += float(qual)
                    qual_n += 1
                except ValueError:
                    pass

            alts = alt.split(",")
            # Per-variant type: call it a SNP only if every ALT is a SNP.
            kinds = {classify(ref, a) for a in alts}
            if kinds == {"SNP"}:
                typ["SNP"] += 1
                # Ti/Tv only for clean biallelic SNPs.
                if len(alts) == 1 and ref.upper() in BASES and alts[0].upper() in BASES:
                    if (ref.upper(), alts[0].upper()) in TRANSITIONS:
                        ti += 1
                    else:
                        tv += 1
            elif "INDEL" in kinds:
                typ["INDEL"] += 1
            elif "MNP" in kinds:
                typ["MNP"] += 1
            else:
                typ["OTHER"] += 1

            # Genotype het/hom across all samples.
            if n_samples and len(f) >= 10:
                for s in f[9:]:
                    gt = parse_gt(s)
                    if gt is None:
                        continue
                    alleles = gt.replace("|", "/").split("/")
                    if len(alleles) != 2:
                        continue
                    a1, a2 = alleles
                    if a1 == "0" and a2 == "0":
                        continue  # hom-ref
                    if a1 == a2:
                        hom_alt += 1
                    else:
                        het += 1

    titv = round(ti / tv, 3) if tv else 0.0
    hethom = round(het / hom_alt, 3) if hom_alt else 0.0
    mean_qual = round(qual_sum / qual_n, 2) if qual_n else 0.0
    pass_rate = round(100.0 * passed / total, 2) if total else 0.0

    # ---- human-readable TSV ----
    tsv_path = f"{args.outdir}/{args.prefix}.summary.tsv"
    with open(tsv_path, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"total_variants\t{total}\n")
        out.write(f"pass_variants\t{passed}\n")
        out.write(f"pass_rate_pct\t{pass_rate}\n")
        out.write(f"snps\t{typ['SNP']}\n")
        out.write(f"indels\t{typ['INDEL']}\n")
        out.write(f"mnps\t{typ['MNP']}\n")
        out.write(f"other\t{typ['OTHER']}\n")
        out.write(f"ti_tv\t{titv}\n")
        out.write(f"het\t{het}\n")
        out.write(f"hom_alt\t{hom_alt}\n")
        out.write(f"het_hom_ratio\t{hethom}\n")
        out.write(f"mean_qual\t{mean_qual}\n")
        out.write(f"n_samples\t{n_samples}\n")
        for c in sorted(per_chrom):
            out.write(f"chrom:{c}\t{per_chrom[c]}\n")

    # ---- MultiQC custom-content (table) ----
    mqc = {
        "id": "vcf_summary",
        "section_name": "VCF Summary",
        "description": "Cohort variant summary computed by vcf_summary.py (no external deps).",
        "plot_type": "table",
        "pconfig": {"id": "vcf_summary_table", "title": "VCF Summary", "scale": False},
        "data": {
            args.prefix: {
                "Total": total,
                "PASS": passed,
                "PASS %": pass_rate,
                "SNPs": typ["SNP"],
                "INDELs": typ["INDEL"],
                "Ti/Tv": titv,
                "Het/Hom": hethom,
                "Mean QUAL": mean_qual,
            }
        },
    }
    mqc_path = f"{args.outdir}/{args.prefix}_vcf_summary_mqc.json"
    with open(mqc_path, "w") as out:
        json.dump(mqc, out, indent=2)

    sys.stderr.write(
        f"[vcf_summary] {total} variants "
        f"({typ['SNP']} SNP / {typ['INDEL']} INDEL), "
        f"Ti/Tv={titv}, Het/Hom={hethom}, PASS={pass_rate}%\n"
    )


if __name__ == "__main__":
    main()
