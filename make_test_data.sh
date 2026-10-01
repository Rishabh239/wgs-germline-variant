#!/usr/bin/env bash
#
# make_test_data.sh - build a REAL but small chr20 HG002 sample dataset.
#
# Produces a chr20-only reference, chr20 known-sites, chr20 GIAB truth, and
# ~30x HG002 chr20 reads (streamed from the GRCh38 300x BAM). Total ~1-3 GB;
# the whole pipeline incl. the hap.py benchmark then runs on a 4-CPU laptop in
# ~80 min. chr20 is the field-standard held-out benchmark chromosome.
#
#   ./make_test_data.sh
#   nextflow run . -profile test,docker --input assets/samplesheet.chr20.csv
#
# Requires: samtools, tabix, bgzip, bwa, gatk, bcftools, awk, wget (all present
# in the pipeline's containers, or install via conda/WSL2).
#
# Tuning:
#   REGION=chr20:1-15000000 ./make_test_data.sh   # even faster (partial chr20)
#   DOWNSAMPLE=0.05          ./make_test_data.sh   # ~15x instead of ~30x
set -euo pipefail

REGION="${REGION:-chr20}"
DOWNSAMPLE="${DOWNSAMPLE:-0.1}"     # fraction of the 300x BAM to keep (0.1 ~ 30x)
mkdir -p reference reads assets

GATK="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
UCSC="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz"
BAM="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam"

echo ">> [1/4] chr20 reference (UCSC hg38 chr20; 'chr20' matches GIAB/GATK naming)"
[[ -s reference/chr20.fa ]] || wget -qO- "$UCSC" | gunzip > reference/chr20.fa
samtools faidx reference/chr20.fa
[[ -s reference/chr20.dict ]] || gatk CreateSequenceDictionary -R reference/chr20.fa -O reference/chr20.dict
[[ -s reference/chr20.fa.bwt ]] || bwa index reference/chr20.fa

echo ">> [2/4] chr20 known-sites (tabix-stream only chr20 from the GATK bundle)"
tabix -h "$GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" chr20 | bgzip > reference/mills.chr20.vcf.gz
tabix -p vcf reference/mills.chr20.vcf.gz
# dbSNP: stream chr20 if the bundle's bgzipped dbSNP is remotely tabix-indexed;
# otherwise fall back to Mills (BQSR/HaplotypeCaller still run — the test just
# won't annotate rsIDs). If you already ran ./download_data.sh --refs, you can
# instead subset locally: bcftools view -r chr20 reference/Homo_sapiens_assembly38.dbsnp138.vcf ...
if tabix -h "$GATK/Homo_sapiens_assembly38.dbsnp138.vcf.gz" chr20 2>/dev/null | bgzip > reference/dbsnp138.chr20.vcf.gz && [[ -s reference/dbsnp138.chr20.vcf.gz ]]; then
    tabix -p vcf reference/dbsnp138.chr20.vcf.gz
    echo "   dbSNP chr20 streamed OK"
else
    echo "   dbSNP remote stream unavailable -> reusing Mills as the known-site/-D source for the test"
    cp reference/mills.chr20.vcf.gz     reference/dbsnp138.chr20.vcf.gz
    cp reference/mills.chr20.vcf.gz.tbi reference/dbsnp138.chr20.vcf.gz.tbi
fi

echo ">> [3/4] chr20 truth (subset the GIAB truth VCF + BED to chr20)"
if [[ -s reference/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz ]]; then
    bcftools view -r chr20 reference/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -Oz -o reference/HG002.chr20.truth.vcf.gz
    tabix -p vcf reference/HG002.chr20.truth.vcf.gz
    awk '$1=="chr20"' reference/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > reference/HG002.chr20.bed
    echo "   truth ready"
else
    echo "   NOTE: run ./download_data.sh --truth first to enable --benchmark"
fi

echo ">> [4/4] HG002 chr20 reads (stream ${REGION} from the 300x GRCh38 BAM, keep ${DOWNSAMPLE})"
# -F 0x900 drops secondary/supplementary so pairs convert cleanly to FASTQ.
samtools view -h -F 0x900 -s "$DOWNSAMPLE" "$BAM" "$REGION" \
    | samtools collate -O -u - \
    | samtools fastq -1 reads/HG002_chr20_R1.fastq.gz -2 reads/HG002_chr20_R2.fastq.gz \
                     -0 /dev/null -s /dev/null -n

cat > assets/samplesheet.chr20.csv <<CSV
sample,fastq_1,fastq_2
HG002,reads/HG002_chr20_R1.fastq.gz,reads/HG002_chr20_R2.fastq.gz
CSV

echo ""
echo "Done. chr20 sample dataset is in reference/ and reads/. Run it with:"
echo "    nextflow run . -profile test,docker --input assets/samplesheet.chr20.csv"
