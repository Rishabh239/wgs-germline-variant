# Usage

## 0. Prerequisites

Nextflow needs a POSIX shell + Java, so **on Windows run everything inside WSL2
(Ubuntu)** with **Docker Desktop** (WSL2 backend) enabled.

```bash
# in WSL2 / Linux / macOS
curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/
nextflow -version          # needs >= 23.04
docker --version
```

## 1. Get the data

```bash
./download_data.sh --refs      # GRCh38 reference + BWA index + known-sites (~30 GB)
./download_data.sh --truth     # GIAB HG002 v4.2.1 truth set (for --benchmark)
./download_data.sh --reads     # HG002 2x250 reads, one lane (edit assets/samplesheet.csv)
```

The reference download is large but one-time. On AWS, stage `reference/` to S3 once
and point `--fasta` etc. at the S3 paths.

## 2. Fast local check on a chr20 sample dataset (WSL2 + Docker)

There is no tiny official HG002 sample — the reads are all full WGS. Build a
**chr20 slice** instead: it's ~1–3 GB, runs the whole pipeline (including the
hap.py benchmark) in ~15–30 min on a laptop, and still gives a real
precision/recall number. chr20 is the field-standard held-out benchmark
chromosome.

```bash
./download_data.sh --truth      # small; needed for the chr20 truth subset
./make_test_data.sh             # builds chr20 reference, known-sites, truth, and ~30x reads
nextflow run . -profile test,docker --input assets/samplesheet.chr20.csv
```

`make_test_data.sh` streams only chr20 (downsampled to ~30x) from the verified
GRCh38 300x HG002 BAM and subsets the reference/known-sites/truth to chr20, so no
20 GB reference download is required. Speed knobs: `REGION=chr20:1-15000000` or
`DOWNSAMPLE=0.05`.

For an even smaller smoke test, downsample the reads further before running.

## 3. Full genome + benchmark (recommended on AWS Batch)

```bash
nextflow run . -profile awsbatch \
    -w s3://YOUR-BUCKET/work \
    --input assets/samplesheet.csv \
    --outdir s3://YOUR-BUCKET/results \
    --benchmark
```

`--benchmark` runs **hap.py** against the GIAB HG002 truth set within the
high-confidence regions and writes SNP/INDEL precision, recall, and F1 to
`benchmark/happy.summary.csv`.

On the RTX laptop you *can* run the full genome via WSL2 + Docker, but it is slow
(GATK is CPU-bound and the GPU does not help); AWS Batch with spot instances is the
practical route and matches the EC2/S3/Batch workflow this project targets.

## Key parameters

| param | default | meaning |
|-------|---------|---------|
| `--input` | `assets/samplesheet.csv` | `sample,fastq_1,fastq_2` |
| `--outdir` | `results` | output location (local path or `s3://...`) |
| `--intervals` | `null` | restrict to a region, e.g. `chr20` |
| `--trim` | `false` | run fastp adapter/quality trimming |
| `--snpeff_db` | `GRCh38.105` | SnpEff database name |
| `--benchmark` | `false` | run hap.py vs GIAB truth |
| `--max_cpus` / `--max_memory` | `8` / `28.GB` | per-task ceilings |

## Notes / honest caveats

- Reference choice: this uses the GATK bundle `Homo_sapiens_assembly38.fasta`
  (alt+decoy+HLA analysis set). Its primary contigs share coordinates with the
  GIAB truth, so benchmarking restricted to the chr1–22 high-confidence BED is
  valid. For a purist match you can call against
  `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta` instead.
- SnpEff downloads its GRCh38 database on first run (needs internet); pre-cache it
  to make runs fully offline/reproducible.
- Filtering is intentionally rule-based (GATK hard filters), not VQSR.
