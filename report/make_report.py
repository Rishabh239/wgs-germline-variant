#!/usr/bin/env python3
"""
make_report.py - build the tables and figures behind report/README.md.

Two stages:
  collect   pipeline outputs under --results  ->  small tidy tables in report/tables/
  plot      report/tables/                    ->  PNG figures in report/figures/

    python report/make_report.py --results results   # after a pipeline run
    python report/make_report.py                     # re-plot from the committed tables

The collect stage only reads files the pipeline already writes (hap.py summary,
extended metrics, ROC tables and annotated VCF; vcf_summary; MarkDuplicates;
MultiQC data; the Nextflow trace), so every number in the report traces back to
a pipeline output. Requires Python 3.9+, pandas >= 1.5 and matplotlib.
"""
import argparse
import gzip
import re
import shutil
from collections import Counter
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"

SNP_FILTERS = ["QD2", "QUAL30", "SOR3", "FS60", "MQ40", "MQRankSum-12.5", "ReadPosRankSum-8"]
INDEL_FILTERS = ["QD2", "QUAL30", "FS200", "ReadPosRankSum-20", "SOR10"]
FILTER_LABELS = {
    "QD2": "QD < 2", "QUAL30": "QUAL < 30", "SOR3": "SOR > 3", "FS60": "FS > 60",
    "MQ40": "MQ < 40", "MQRankSum-12.5": "MQRankSum < -12.5",
    "ReadPosRankSum-8": "ReadPosRankSum < -8",
}
INDEL_SIZES = ["D16_PLUS", "D6_15", "D1_5", "I1_5", "I6_15", "I16_PLUS"]
SIZE_LABELS = {
    "D16_PLUS": "Del ≥16 bp", "D6_15": "Del 6–15 bp", "D1_5": "Del 1–5 bp",
    "I1_5": "Ins 1–5 bp", "I6_15": "Ins 6–15 bp", "I16_PLUS": "Ins ≥16 bp",
}
# Pipeline steps in execution order, as named in the Nextflow trace.
STEP_LABELS = {
    "FASTQC": "FastQC", "BWA_MEM": "BWA-MEM", "SAMTOOLS_SORT": "samtools sort",
    "MARKDUPLICATES": "MarkDuplicates", "BASERECALIBRATOR": "BaseRecalibrator",
    "APPLYBQSR": "ApplyBQSR", "HAPLOTYPECALLER": "HaplotypeCaller",
    "COMBINE_GVCFS": "CombineGVCFs", "GENOTYPE_GVCFS": "GenotypeGVCFs",
    "HARD_FILTER": "Hard filters", "HAPPY": "hap.py benchmark",
}
STEPS = list(STEP_LABELS)


def f1(precision, recall):
    return 2 * precision * recall / (precision + recall) if precision + recall else 0.0


# ============================================================================
# collect: pipeline outputs -> report/tables/
# ============================================================================
def collect(results, truth_bed, fasta):
    TABLES.mkdir(parents=True, exist_ok=True)
    bench = results / "benchmark"

    # hap.py's own summary + stratified metrics, copied verbatim.
    shutil.copy(bench / "happy.summary.csv", TABLES / "happy.summary.csv")
    shutil.copy(bench / "happy.extended.csv", TABLES / "happy.extended.csv")
    summary = pd.read_csv(bench / "happy.summary.csv")
    ext = pd.read_csv(bench / "happy.extended.csv")

    collect_pr_curves(bench)
    collect_filter_audit(bench, summary)
    collect_qual_comparison(summary)
    collect_indel_sizes(ext)
    collect_genotypes(ext)
    collect_qc(results, ext, truth_bed, fasta)
    collect_runtime(results / "pipeline_info" / "trace.txt")

    # The pipeline's own MultiQC report, kept alongside the tables.
    report = next((results / "multiqc").glob("*multiqc_report.html"))
    shutil.copy(report, HERE / "multiqc_report.html")


def collect_pr_curves(bench):
    """Precision/recall as the QUAL threshold rises (hap.py --roc QUAL), all calls."""
    parts = []
    for vtype in ("SNP", "INDEL"):
        roc = pd.read_csv(bench / f"happy.roc.Locations.{vtype}.csv.gz")
        roc = roc[(roc["Subtype"] == "*") & (roc["Subset"] == "*")]
        part = roc[["QQ", "METRIC.Recall", "METRIC.Precision"]].rename(columns={
            "QQ": "qual_threshold", "METRIC.Recall": "recall", "METRIC.Precision": "precision"})
        part.insert(0, "type", vtype)
        parts.append(part.sort_values("qual_threshold"))
    out = pd.concat(parts).dropna(subset=["precision"])
    out.to_csv(TABLES / "pr_curve_qual.csv", index=False, float_format="%.6g",
               lineterminator="\n")


def collect_filter_audit(bench, summary):
    """
    Which hard filters removed true variants vs false positives.

    hap.py annotates every query call with a decision (TP/FP/UNK). Its PASS
    metrics are the ALL metrics with filtered calls removed, so the effect of
    dropping one filter can be re-scored exactly from these annotations:
    calls whose only failing filter is X return to the callset.
    """
    fails = Counter()   # (type, filter, decision) -> calls failing that filter
    combos = Counter()  # (type, exact filter set, decision) -> calls
    with gzip.open(bench / "happy.vcf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[6] in ("PASS", "."):
                continue
            keys = f[8].split(":")
            query = dict(zip(keys, f[10].split(":")))
            decision, vtype = query.get("BD"), query.get("BVT")
            if decision not in ("TP", "FP"):
                continue  # outside the high-confidence regions (UNK) or no call
            flt = tuple(sorted(f[6].split(";")))
            combos[(vtype, flt, decision)] += 1
            for x in flt:
                fails[(vtype, x, decision)] += 1

    audit = []
    for vtype, filters in (("SNP", SNP_FILTERS), ("INDEL", INDEL_FILTERS)):
        for x in filters:
            tp, fp = fails[(vtype, x, "TP")], fails[(vtype, x, "FP")]
            audit.append({
                "type": vtype, "filter": x,
                "tp_removed": tp, "fp_removed": fp,
                "tp_removed_only_by_this": combos[(vtype, (x,), "TP")],
                "fp_removed_only_by_this": combos[(vtype, (x,), "FP")],
                "fp_share_of_removed": round(fp / (tp + fp), 4) if tp + fp else None,
            })
    audit = pd.DataFrame(audit)
    audit.to_csv(TABLES / "filter_audit.csv", index=False, lineterminator="\n")

    # Leave-one-filter-out re-scoring (SNPs; INDEL filters remove only a handful of calls).
    s = summary.set_index(["Type", "Filter"])
    p = s.loc[("SNP", "PASS")]
    truth_total, truth_tp, query_fp = p["TRUTH.TOTAL"], p["TRUTH.TP"], p["QUERY.FP"]
    query_tp = p["QUERY.TOTAL"] - p["QUERY.FP"] - p["QUERY.UNK"]

    def score(label, add_tp=0, add_fp=0):
        rec = (truth_tp + add_tp) / truth_total
        prec = (query_tp + add_tp) / (query_tp + add_tp + query_fp + add_fp)
        return {"callset": label, "tp_returned": int(add_tp), "fp_returned": int(add_fp),
                "recall": round(rec, 6), "precision": round(prec, 6),
                "f1": round(f1(prec, rec), 6)}

    snp = audit[audit["type"] == "SNP"]
    loo = [score("All hard filters (PASS)")]
    for r in snp.itertuples():
        if r.tp_removed or r.fp_removed:
            loo.append(score(f"Without {FILTER_LABELS[r.filter]}",
                             r.tp_removed_only_by_this, r.fp_removed_only_by_this))
    all_tp = sum(n for (v, _f, d), n in combos.items() if v == "SNP" and d == "TP")
    all_fp = sum(n for (v, _f, d), n in combos.items() if v == "SNP" and d == "FP")
    loo.append(score("No filters (ALL)", all_tp, all_fp))

    # Returning every filtered call must reproduce hap.py's own ALL row; if it
    # does not, the leave-one-out numbers above cannot be trusted.
    a = s.loc[("SNP", "ALL")]
    assert abs(loo[-1]["recall"] - a["METRIC.Recall"]) < 1e-5, "ALL recall mismatch"
    assert abs(loo[-1]["precision"] - a["METRIC.Precision"]) < 1e-5, "ALL precision mismatch"
    pd.DataFrame(loo).to_csv(TABLES / "filter_leave_one_out_snp.csv", index=False, lineterminator="\n")


def collect_qual_comparison(summary):
    """
    Could a single QUAL cut-off replace the hard filters? For each filtered
    operating point, find the QUAL threshold on the all-calls curve that reaches
    the same precision, and report the recall it keeps.
    """
    curve = pd.read_csv(TABLES / "pr_curve_qual.csv")
    loo = pd.read_csv(TABLES / "filter_leave_one_out_snp.csv").set_index("callset")
    s = summary.set_index(["Type", "Filter"])
    points = [
        ("SNP", "All hard filters (PASS)",
         s.loc[("SNP", "PASS"), "METRIC.Precision"], s.loc[("SNP", "PASS"), "METRIC.Recall"]),
        ("SNP", "Hard filters without SOR > 3",
         loo.loc["Without SOR > 3", "precision"], loo.loc["Without SOR > 3", "recall"]),
        ("INDEL", "All hard filters (PASS)",
         s.loc[("INDEL", "PASS"), "METRIC.Precision"], s.loc[("INDEL", "PASS"), "METRIC.Recall"]),
    ]
    rows = []
    for vtype, label, prec, rec in points:
        c = curve[(curve["type"] == vtype) & (curve["precision"] >= prec)]
        best = c.loc[c["recall"].idxmax()]
        rows.append({"type": vtype, "filtered_callset": label,
                     "precision": round(prec, 6), "recall": round(rec, 6),
                     "qual_cutoff_same_precision": round(best["qual_threshold"], 1),
                     "qual_cutoff_recall": round(best["recall"], 6),
                     "recall_difference_pts": round(100 * (best["recall"] - rec), 2)})
    pd.DataFrame(rows).to_csv(TABLES / "qual_cutoff_comparison.csv", index=False, lineterminator="\n")


def collect_indel_sizes(ext):
    sel = ext[(ext["Type"] == "INDEL") & (ext["Subset"] == "*") & (ext["Filter"] == "PASS")
              & (ext["Subtype"].isin(INDEL_SIZES))].set_index("Subtype").loc[INDEL_SIZES]
    out = pd.DataFrame({
        "size_class": [SIZE_LABELS[k] for k in sel.index],
        "truth_total": sel["TRUTH.TOTAL"].values,
        "tp": sel["TRUTH.TP"].values,
        "fn": sel["TRUTH.FN"].values,
        "fp": sel["QUERY.FP"].values,
        "recall": sel["METRIC.Recall"].round(6).values,
        "precision": sel["METRIC.Precision"].round(6).values,
        "f1": sel["METRIC.F1_Score"].round(6).values,
    })
    out.to_csv(TABLES / "indel_by_size.csv", index=False, lineterminator="\n")


def collect_genotypes(ext):
    rows = []
    for vtype in ("SNP", "INDEL"):
        r = ext[(ext["Type"] == vtype) & (ext["Subtype"] == "*") & (ext["Subset"] == "*")
                & (ext["Filter"] == "PASS")].iloc[0]
        for gt in ("het", "homalt"):
            rows.append({"type": vtype, "genotype": gt,
                         "truth_total": int(r[f"TRUTH.TOTAL.{gt}"]),
                         "tp": int(r[f"TRUTH.TP.{gt}"]), "fn": int(r[f"TRUTH.FN.{gt}"]),
                         "recall": round(r[f"TRUTH.TP.{gt}"] / r[f"TRUTH.TOTAL.{gt}"], 6),
                         "fp": int(r[f"QUERY.FP.{gt}"])})
        rows.append({"type": vtype, "genotype": "FP: wrong genotype", "fp": int(r["FP.gt"])})
        rows.append({"type": vtype, "genotype": "FP: wrong allele", "fp": int(r["FP.al"])})
    out = pd.DataFrame(rows)
    counts = ["truth_total", "tp", "fn", "fp"]
    out[counts] = out[counts].astype("Int64")   # keep counts integral alongside empty cells
    out.to_csv(TABLES / "benchmark_by_genotype.csv", index=False, lineterminator="\n")


def collect_qc(results, ext, truth_bed, fasta):
    mqc = next((results / "multiqc").glob("*_data"))
    fq = pd.read_csv(mqc / "multiqc_fastqc.txt", sep="\t")
    flag = pd.read_csv(mqc / "multiqc_samtools_flagstat.txt", sep="\t").iloc[0]
    md_file = next((results / "qc" / "markduplicates").glob("*.metrics.txt"))
    md = pd.read_csv(md_file, sep="\t", comment="#", nrows=1).iloc[0]
    vs = pd.read_csv(results / "qc" / "vcf_summary" / "cohort.summary.tsv", sep="\t",
                     index_col=0)["value"]
    snp = ext[(ext["Type"] == "SNP") & (ext["Subtype"] == "*") & (ext["Subset"] == "*")
              & (ext["Filter"] == "PASS")].iloc[0]

    reads = int(fq["Total Sequences"].sum())
    read_len = int(float(str(fq["Sequence length"].iloc[0]).split("-")[-1]))
    bases = int(reads * read_len)
    rows = [
        ("Read pairs", reads // 2, "FastQC"),
        ("Read length (bp)", read_len, "FastQC"),
        ("Sequenced bases", bases, "FastQC"),
    ]
    if fasta.exists():
        length = n = 0
        with open(fasta) as fh:
            for line in fh:
                if not line.startswith(">"):
                    s = line.strip()
                    length += len(s)
                    n += s.upper().count("N")
        rows += [("Reference length, non-N (bp)", length - n, fasta.name),
                 ("Raw depth (bases / non-N length)", round(bases / (length - n), 1), "derived")]
    if truth_bed.exists():
        conf = sum(int(c[2]) - int(c[1]) for c in
                   (l.split("\t") for l in truth_bed.read_text().splitlines() if l.strip()))
        rows.append(("High-confidence region (bp)", conf, truth_bed.name))
    rows += [
        ("Mapped reads (%)", flag["mapped_passed_pct"], "samtools flagstat"),
        ("Properly paired (%)", flag["properly paired_passed_pct"], "samtools flagstat"),
        ("Duplicate rate (%)", round(100 * md["PERCENT_DUPLICATION"], 2), "MarkDuplicates"),
        ("Variant records", int(vs["total_variants"]), "vcf_summary.py"),
        ("PASS records", int(vs["pass_variants"]), "vcf_summary.py"),
        ("PASS rate (%)", vs["pass_rate_pct"], "vcf_summary.py"),
        ("SNP records", int(vs["snps"]), "vcf_summary.py"),
        ("INDEL records", int(vs["indels"]), "vcf_summary.py"),
        ("Ti/Tv, all records", vs["ti_tv"], "vcf_summary.py"),
        ("Het / hom-alt ratio, all records", vs["het_hom_ratio"], "vcf_summary.py"),
        ("Ti/Tv, PASS callset", round(snp["QUERY.TOTAL.TiTv_ratio"], 2), "hap.py"),
        ("Ti/Tv, PASS true positives", round(snp["QUERY.TP.TiTv_ratio"], 2), "hap.py"),
        ("Ti/Tv, PASS false positives", round(snp["QUERY.FP.TiTv_ratio"], 2), "hap.py"),
        ("Ti/Tv, PASS calls outside high-confidence regions",
         round(snp["QUERY.UNK.TiTv_ratio"], 2), "hap.py"),
        ("Ti/Tv, GIAB truth", round(snp["TRUTH.TOTAL.TiTv_ratio"], 2), "hap.py"),
    ]
    # object dtype keeps integer counts and decimal metrics in their own formats.
    pd.DataFrame({"metric": [r[0] for r in rows],
                  "value": pd.Series([r[1] for r in rows], dtype=object),
                  "source": [r[2] for r in rows]}).to_csv(TABLES / "qc_metrics.csv", index=False, lineterminator="\n")


def to_minutes(text):
    """Nextflow trace durations look like '1h 2m 3s', '32m 21s', '17.3s' or '637ms'."""
    scale = {"d": 1440, "h": 60, "m": 1, "s": 1 / 60, "ms": 1 / 60000}
    return sum(float(num) * scale[unit]
               for num, unit in re.findall(r"([\d.]+)\s*(ms|d|h|m|s)", str(text)))


def collect_runtime(trace):
    t = pd.read_csv(trace, sep="\t")
    t = t[t["status"].isin(["COMPLETED", "CACHED"])].copy()
    t["step"] = t["name"].str.replace(r"\s*\(.*\)$", "", regex=True)
    t = t[t["step"].isin(STEPS)]
    out = pd.DataFrame({
        "step": t["step"],
        "wall_min": t["realtime"].map(to_minutes).round(2),
        "cpu_pct": t["%cpu"].str.rstrip("%").astype(float),
        "peak_rss": t["peak_rss"],
    })
    out = out.assign(order=out["step"].map(STEPS.index)).sort_values("order")
    out.drop(columns="order").to_csv(TABLES / "runtime.csv", index=False, lineterminator="\n")


# ============================================================================
# plot: report/tables/ -> report/figures/
# ============================================================================
# Validated categorical slots 1-3 (blue, orange, aqua) on a light chart surface.
SURFACE = "#fcfcfb"
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
PX = 1 / 100   # one display pixel, in inches: figures are laid out at 100 px/in, saved at 2x
WIDTH = 9.2    # inches -> ~920 px, about the width of a GitHub README column


def style(plt):
    import logging
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)  # quiet font fallback
    plt.rcParams.update({
        "font.family": ["Segoe UI", "Arial", "DejaVu Sans"],
        "font.size": 9.5, "axes.titlesize": 10.5, "axes.titleweight": "semibold",
        "axes.titlelocation": "left", "axes.titlecolor": INK, "axes.titlepad": 8,
        "axes.labelcolor": INK_2, "axes.labelsize": 9.5,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.75, "axes.facecolor": SURFACE,
        "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.75, "grid.linestyle": "-",
        "axes.axisbelow": True,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK_2,
        "ytick.labelcolor": INK_2, "xtick.major.size": 0, "ytick.major.size": 0,
        "xtick.major.pad": 5, "ytick.major.pad": 5,
        "legend.frameon": False, "legend.fontsize": 9, "legend.labelcolor": INK_2,
    })


def header(fig, title, subtitle):
    """Left-aligned title + subtitle above the plot area."""
    fig.text(0.012, 0.985, title, ha="left", va="top", fontsize=12, fontweight="semibold",
             color=INK)
    fig.text(0.012, 0.915, subtitle, ha="left", va="top", fontsize=9.5, color=INK_2)


def despine(ax, keep=("bottom",)):
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)


def hbar(ax, y, width, height, color, radius_px=4):
    """Horizontal bar from x=0 with a rounded data end and a square baseline end.

    Call after the layout is final: the corner radius is converted from display
    pixels to data units using the axes' current size.
    """
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as MPath
    bbox = ax.get_window_extent().transformed(ax.figure.dpi_scale_trans.inverted())
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    rx = min(radius_px * PX * abs(x1 - x0) / bbox.width, width / 2)
    ry = min(radius_px * PX * abs(y1 - y0) / bbox.height, height / 2)
    b, t = y - height / 2, y + height / 2
    verts = [(0, b), (width - rx, b), (width, b), (width, b + ry), (width, t - ry),
             (width, t), (width - rx, t), (0, t), (0, b)]
    codes = [MPath.MOVETO, MPath.LINETO, MPath.CURVE3, MPath.CURVE3, MPath.LINETO,
             MPath.CURVE3, MPath.CURVE3, MPath.LINETO, MPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MPath(verts, codes), facecolor=color, edgecolor="none"))


def dot(ax, x, y, color, **kw):
    ax.plot(x, y, "o", ms=8.5, color=color, mec=SURFACE, mew=1.5, zorder=4, **kw)


def pct_axis(axis, step, decimals):
    from matplotlib.ticker import MultipleLocator, PercentFormatter
    axis.set_major_locator(MultipleLocator(step))
    axis.set_major_formatter(PercentFormatter(1.0, decimals=decimals))


def plot_pr_curves(plt):
    curve = pd.read_csv(TABLES / "pr_curve_qual.csv")
    summ = pd.read_csv(TABLES / "happy.summary.csv").set_index(["Type", "Filter"])
    loo = pd.read_csv(TABLES / "filter_leave_one_out_snp.csv").set_index("callset")

    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 4.4))
    windows = {"SNP": ((0.95, 1.0, 0.01), (0.985, 1.0, 0.005)),
               "INDEL": ((0.90, 0.99, 0.02), (0.975, 1.0, 0.005))}
    for ax, vtype in zip(axes, ("SNP", "INDEL")):
        c = curve[curve["type"] == vtype]
        (xlo, xhi, xstep), (ylo, yhi, ystep) = windows[vtype]
        ax.plot(c["recall"], c["precision"], color=BLUE, lw=1.5, solid_capstyle="round",
                solid_joinstyle="round", zorder=2)
        a, p = summ.loc[(vtype, "ALL")], summ.loc[(vtype, "PASS")]
        ax_, ay = a["METRIC.Recall"], a["METRIC.Precision"]
        px_, py = p["METRIC.Recall"], p["METRIC.Precision"]
        dot(ax, ax_, ay, BLUE)
        dot(ax, px_, py, ORANGE)
        if vtype == "SNP":
            r = loo.loc["Without SOR > 3"]
            dot(ax, r["recall"], r["precision"], AQUA)
            ax.annotate("No filters", (ax_, ay), xytext=(0, -12), textcoords="offset points",
                        ha="center", va="top", color=INK_2, fontsize=9)
            ax.annotate("GATK hard filters", (px_, py), xytext=(0, 10),
                        textcoords="offset points", ha="center", va="bottom", color=INK,
                        fontsize=9, fontweight="semibold")
            ax.annotate("Without\nSOR > 3", (r["recall"], r["precision"]), xytext=(8, 5),
                        textcoords="offset points", ha="left", va="bottom", color=INK_2,
                        fontsize=9, linespacing=1.2)
        else:
            # Only a handful of benchmarked INDELs fail a hard filter, so the points coincide.
            ax.annotate("Hard filters ≈ no filters: only 6\nbenchmarked INDELs fail them",
                        (px_, py), xytext=(-10, -12), textcoords="offset points", ha="right",
                        va="top", color=INK, fontsize=9, linespacing=1.3)
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(ylo, yhi)
        pct_axis(ax.xaxis, xstep, 0)
        pct_axis(ax.yaxis, ystep, 1)
        ax.set_title(f"{vtype}s")
        ax.set_xlabel("Recall")
        despine(ax)
    axes[0].set_ylabel("Precision")
    handles = [plt.Line2D([], [], color=BLUE, lw=1.5),
               plt.Line2D([], [], marker="o", ls="", ms=8, color=BLUE, mec=SURFACE),
               plt.Line2D([], [], marker="o", ls="", ms=8, color=ORANGE, mec=SURFACE),
               plt.Line2D([], [], marker="o", ls="", ms=8, color=AQUA, mec=SURFACE)]
    fig.legend(handles, ["Raising a QUAL cut-off (all calls)", "No filters",
                         "GATK hard filters (reported callset)",
                         "Hard filters without SOR > 3 (re-scored)"],
               loc="lower left", ncol=4, bbox_to_anchor=(0.005, 0.0), handlelength=1.6,
               columnspacing=1.3, handletextpad=0.5)
    header(fig, "Precision and recall against GIAB HG002 v4.2.1 (chr20)",
           "Axes are zoomed to the high-accuracy corner; up and to the right is better.")
    fig.tight_layout(rect=(0, 0.07, 1, 0.86), w_pad=2.5)
    save(fig, "fig1_precision_recall.png")


def plot_filter_audit(plt):
    audit = pd.read_csv(TABLES / "filter_audit.csv")
    snp = audit[audit["type"] == "SNP"]
    idle = [FILTER_LABELS[f] for f in snp[(snp["tp_removed"] + snp["fp_removed"]) == 0]["filter"]]
    snp = snp[(snp["tp_removed"] + snp["fp_removed"]) > 0].sort_values("tp_removed")

    n = len(snp)
    fig, ax = plt.subplots(figsize=(WIDTH, 1.25 + 0.72 * n))
    xmax = snp[["tp_removed", "fp_removed"]].to_numpy().max() * 1.42
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_yticks(range(n))
    ax.set_yticklabels([FILTER_LABELS[f] for f in snp["filter"]], color=INK, fontsize=10)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Benchmarked SNP calls that failed the filter")
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    despine(ax)
    handles = [plt.Rectangle((0, 0), 1, 1, color=BLUE), plt.Rectangle((0, 0), 1, 1, color=ORANGE)]
    fig.legend(handles, ["True variants removed (recall lost)",
                         "False positives removed (precision gained)"],
               loc="lower left", bbox_to_anchor=(0.005, 0.0), ncol=2, handlelength=1.0,
               columnspacing=1.6)
    header(fig, "Which SNP hard filters earn their keep",
           f"No high-confidence calls failed {', '.join(idle[:-1])} or {idle[-1]}. "
           "A call can fail more than one filter.")
    fig.tight_layout(rect=(0, 0.11, 1, 0.80))

    bar_h, gap = 0.30, 0.02          # 2 px-ish surface gap between the pair
    for i, r in enumerate(snp.itertuples()):
        for value, y, color in ((r.tp_removed, i + (bar_h + gap) / 2, BLUE),
                                (r.fp_removed, i - (bar_h + gap) / 2, ORANGE)):
            hbar(ax, y, value, bar_h, color)
            ax.text(value + xmax * 0.01, y, f"{value:,}", va="center", color=INK_2, fontsize=9)
        share = r.fp_removed / (r.tp_removed + r.fp_removed)
        ax.text(xmax, i, f"{share:.0%} of the calls it\nremoved were errors", va="center",
                ha="right", color=INK, fontsize=9, linespacing=1.3)
    save(fig, "fig2_filter_audit.png")


def plot_indel_sizes(plt):
    d = pd.read_csv(TABLES / "indel_by_size.csv")
    fig, ax = plt.subplots(figsize=(WIDTH, 3.9))
    x = list(range(len(d)))
    ax.plot(x, d["recall"], color=BLUE, lw=0, marker="o", ms=8.5, mec=SURFACE, mew=1.5, zorder=3)
    ax.plot(x, d["precision"], color=ORANGE, lw=0, marker="o", ms=8.5, mec=SURFACE, mew=1.5,
            zorder=3)
    for i in (d["recall"].idxmin(), d["recall"].idxmax()):
        ax.annotate(f"{d.loc[i, 'recall']:.1%}", (i, d.loc[i, "recall"]), xytext=(9, 0),
                    textcoords="offset points", va="center", color=INK, fontsize=9)
    ax.axvline(2.5, color=AXIS, lw=0.75)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}\nn = {n:,}" for s, n in zip(d["size_class"], d["truth_total"])],
                       fontsize=9, linespacing=1.4)
    ax.set_xlim(-0.5, len(d) - 0.5)
    ax.set_ylim(0.92, 1.005)
    pct_axis(ax.yaxis, 0.02, 0)
    ax.grid(axis="x", visible=False)
    despine(ax)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=8, color=c, mec=SURFACE)
               for c in (BLUE, ORANGE)]
    fig.legend(handles, ["Recall", "Precision"], loc="lower left", ncol=2,
               bbox_to_anchor=(0.005, 0.0), handletextpad=0.3)
    header(fig, "INDEL accuracy by event size (hard-filtered calls)",
           "Deletions left of the divider, insertions right. "
           "n = GIAB truth INDELs in each size class.")
    fig.tight_layout(rect=(0, 0.08, 1, 0.84))
    save(fig, "fig3_indel_size.png")


def plot_runtime(plt):
    d = pd.read_csv(TABLES / "runtime.csv").iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(WIDTH, 0.9 + 0.34 * len(d)))
    xmax = d["wall_min"].max() * 1.15
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.6, len(d) - 0.4)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([STEP_LABELS[s] for s in d["step"]], color=INK, fontsize=9.5)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Wall-clock minutes per task")
    despine(ax)
    sequential = d[d["step"] != "FASTQC"]["wall_min"].sum()
    header(fig, "Where the time goes: one chr20 run on a laptop",
           "WSL2 + Docker, test profile capped at 4 CPUs / 6 GB. FastQC runs alongside "
           f"alignment; the other steps run in order and sum to {sequential:.0f} min.")
    fig.tight_layout(rect=(0, 0, 1, 0.85))

    for i, r in d.iterrows():
        hbar(ax, i, r["wall_min"], 0.56, BLUE)
        label = f"{r['wall_min']:.1f} min" if r["wall_min"] >= 1 else f"{r['wall_min'] * 60:.0f} s"
        ax.text(r["wall_min"] + xmax * 0.008, i, label, va="center", color=INK_2, fontsize=9)
    save(fig, "fig4_runtime.png")


def save(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / name, dpi=200)
    print(f"[make_report] wrote {FIGURES / name}")


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    style(plt)
    plot_pr_curves(plt)
    plot_filter_audit(plt)
    plot_indel_sizes(plt)
    plot_runtime(plt)
    plt.close("all")


def main():
    ap = argparse.ArgumentParser(description="Build report tables and figures.")
    ap.add_argument("--results", type=Path, default=None,
                    help="pipeline --outdir to collect from (default: re-plot committed tables)")
    ap.add_argument("--truth_bed", type=Path, default=Path("reference/HG002.chr20.bed"))
    ap.add_argument("--fasta", type=Path, default=Path("reference/chr20.fa"))
    args = ap.parse_args()
    if args.results:
        collect(args.results, args.truth_bed, args.fasta)
        print(f"[make_report] tables written to {TABLES}")
    plot()


if __name__ == "__main__":
    main()
