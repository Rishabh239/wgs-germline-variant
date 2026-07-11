# wgs-germline-nf

A reproducible **germline whole-genome short-variant pipeline** that takes raw reads
to an annotated, filtered VCF following **GATK4 best practices**, orchestrated in
**Nextflow DSL2** with per-process **Docker** containers and an **AWS Batch** profile.

Variant filtering is **rule-based** (GATK hard filters) — deliberately **no machine
learning** (no VQSR, which fits a Gaussian mixture model). Calls are benchmarked
against the **Genome in a Bottle HG002** truth set with **hap.py**.

```
FASTQ ─► FastQC ─► BWA-MEM ─► sort ─► MarkDuplicates ─► BQSR
      ─► HaplotypeCaller (GVCF) ─► CombineGVCFs ─► GenotypeGVCFs
      ─► hard-filter (SNP/INDEL) ─► SnpEff ─► bcftools stats + VCF summary ─► MultiQC
                                          └► hap.py vs GIAB HG002  (--benchmark)
```

See [`docs/architecture.md`](docs/architecture.md) for the full DAG.

## Quick start

```bash
# 1. reference + known-sites + GIAB truth (one-time, large)
./download_data.sh --refs --truth --reads     # or: --all

# 2. fast local end-to-end check on chr20 (WSL2 + Docker on Windows)
nextflow run . -profile test,docker --input assets/samplesheet.csv

# 3. full genome + benchmark (recommended on AWS Batch)
nextflow run . -profile awsbatch -w s3://BUCKET/work \
    --input assets/samplesheet.csv --outdir s3://BUCKET/results --benchmark
```

Full instructions, parameters, and caveats: [`docs/usage.md`](docs/usage.md).
Output layout: [`docs/output.md`](docs/output.md).

## What each stage does

| Stage | Tool | Purpose |
|-------|------|---------|
| Read QC | FastQC (+ optional fastp) | raw-read quality; optional trimming |
| Align | BWA-MEM (alt-aware) → samtools | map reads to GRCh38, sort, index |
| Dedup | GATK MarkDuplicates | flag PCR/optical duplicates |
| Recalibrate | BaseRecalibrator + ApplyBQSR | correct systematic base-quality error |
| Call | HaplotypeCaller (`-ERC GVCF`) | per-sample variant calling |
| Joint-genotype | CombineGVCFs + GenotypeGVCFs | cohort VCF |
| Filter | SelectVariants + VariantFiltration | **rule-based** SNP/INDEL hard filters |
| Annotate | SnpEff + SnpSift | gene / consequence / impact |
| QC | bcftools stats, `vcf_summary.py`, MultiQC | Ti/Tv, het/hom, PASS rate, one report |
| Benchmark | hap.py vs GIAB HG002 v4.2.1 | precision / recall / F1 (`--benchmark`) |

## Hard-filter thresholds (GATK germline)

**SNPs:** `QD<2 · QUAL<30 · SOR>3 · FS>60 · MQ<40 · MQRankSum<-12.5 · ReadPosRankSum<-8`
**INDELs:** `QD<2 · QUAL<30 · FS>200 · ReadPosRankSum<-20 · SOR>10`

## Data sources

- Reference + known-sites: GATK GRCh38 resource bundle (`gcp-public-data--broad-references/hg38/v0`)
- Sample: GIAB **HG002** (Ashkenazi son), Illumina 2×250
- Truth: GIAB **HG002 GRCh38 v4.2.1** benchmark VCF + high-confidence BED

## Repo layout

```
main.nf · nextflow.config · download_data.sh
modules/            fastqc · fastp · bwa_mem · samtools_sort · markduplicates · bqsr ·
                    haplotypecaller · joint_genotyping · hard_filter · snpeff ·
                    bcftools_stats · vcf_summary · happy · multiqc
bin/                vcf_summary.py · check_samplesheet.py   (dependency-free, unit-tested)
assets/             samplesheet.csv · multiqc_config.yaml
docs/               usage · output · architecture
```

Execution profiles (`docker`, `singularity`, `test`, `awsbatch`) live in `nextflow.config`.

## Requirements

Nextflow ≥ 23.04, Java 11+, and Docker (or Singularity). On Windows, run inside
**WSL2 + Docker Desktop**. Containers are pulled per process, so no manual tool
installs are needed.
