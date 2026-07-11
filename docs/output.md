# Outputs (under `results/`)

```
results/
├── alignment/            *.recal.bam (+ .bai)          final analysis-ready BAM
├── variants/
│   ├── cohort.vcf.gz                                    joint-genotyped calls
│   └── cohort.filtered.vcf.gz                           hard-filtered (FILTER column set)
├── annotation/
│   ├── cohort.ann.vcf                                   SnpEff-annotated
│   └── snpEff_summary.html / snpEff.csv / snpEff_genes.txt
├── qc/
│   ├── fastqc/            per-sample FastQC
│   ├── markduplicates/    duplicate metrics
│   ├── bcftools/          bcftools stats
│   └── vcf_summary/       cohort.summary.tsv (Ti/Tv, het/hom, PASS rate, ...)
├── benchmark/            happy.summary.csv (SNP/INDEL precision, recall, F1)   [--benchmark]
├── multiqc/              multiqc_report.html
└── pipeline_info/        timeline / report / trace / dag
```

**Reading the callset:** `FILTER == PASS` passed all hard filters; other tags
(e.g. `QD2`, `FS60`) name the filter(s) a record failed. The raw records are
retained (soft-filtered) so nothing is silently dropped.
