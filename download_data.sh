#!/usr/bin/env bash
#
# download_data.sh - fetch reference, known-sites, GIAB truth, and HG002 reads.
#
#   ./download_data.sh --refs     # GRCh38 reference + BWA index + known-sites (GATK bundle)
#   ./download_data.sh --truth    # GIAB HG002 v4.2.1 truth VCF + high-confidence BED
#   ./download_data.sh --reads    # HG002 Illumina 2x250 reads (one lane; large)
#   ./download_data.sh --all      # everything
#
# All files land in ./reference (refs/truth) and ./reads (fastq).
set -euo pipefail

REF_DIR="reference"
READS_DIR="reads"
mkdir -p "$REF_DIR" "$READS_DIR"

GATK="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
GIAB_TRUTH="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
GIAB_READS="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads"

fetch() {  # fetch <url> <dest>
    local url="$1" dest="$2"
    if [[ -s "$dest" ]]; then
        echo "  [skip] $(basename "$dest") already present"
        return
    fi
    echo "  [get ] $(basename "$dest")"
    wget -q --show-progress -c -O "$dest" "$url"
}

refs() {
    echo ">> Reference + BWA index + known-sites (GATK GRCh38 bundle)"
    fetch "$GATK/Homo_sapiens_assembly38.fasta"                        "$REF_DIR/Homo_sapiens_assembly38.fasta"
    fetch "$GATK/Homo_sapiens_assembly38.fasta.fai"                    "$REF_DIR/Homo_sapiens_assembly38.fasta.fai"
    fetch "$GATK/Homo_sapiens_assembly38.dict"                         "$REF_DIR/Homo_sapiens_assembly38.dict"
    for ext in alt amb ann bwt pac sa; do
        fetch "$GATK/Homo_sapiens_assembly38.fasta.64.$ext"           "$REF_DIR/Homo_sapiens_assembly38.fasta.64.$ext"
    done
    # BWA looks for <fasta>.{amb,ann,bwt,pac,sa}; the bundle ships them as
    # .64.* so we symlink to the names BWA expects.
    ( cd "$REF_DIR"
      for ext in alt amb ann bwt pac sa; do
          ln -sf "Homo_sapiens_assembly38.fasta.64.$ext" "Homo_sapiens_assembly38.fasta.$ext"
      done )
    fetch "$GATK/Homo_sapiens_assembly38.dbsnp138.vcf"                "$REF_DIR/Homo_sapiens_assembly38.dbsnp138.vcf"
    fetch "$GATK/Homo_sapiens_assembly38.dbsnp138.vcf.idx"            "$REF_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
    fetch "$GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"    "$REF_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    fetch "$GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" "$REF_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
    fetch "$GATK/Homo_sapiens_assembly38.known_indels.vcf.gz"         "$REF_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz"
    fetch "$GATK/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"     "$REF_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
}

truth() {
    echo ">> GIAB HG002 v4.2.1 truth set (for hap.py benchmarking)"
    fetch "$GIAB_TRUTH/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"                 "$REF_DIR/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    fetch "$GIAB_TRUTH/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"             "$REF_DIR/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
    fetch "$GIAB_TRUTH/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"     "$REF_DIR/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
}

reads() {
    echo ">> HG002 Illumina 2x250 reads (one lane ~ a few GB)"
    echo "   NOTE: this single lane is low-coverage. For a real ~30x run use the"
    echo "   HiSeq300x FASTQs, or add more lanes; see docs/usage.md."
    fetch "$GIAB_READS/D1_S1_L001_R1_001.fastq.gz" "$READS_DIR/HG002_R1.fastq.gz"
    fetch "$GIAB_READS/D1_S1_L001_R2_001.fastq.gz" "$READS_DIR/HG002_R2.fastq.gz"
    echo ""
    echo "   Samplesheet row:"
    echo "     HG002,$(pwd)/$READS_DIR/HG002_R1.fastq.gz,$(pwd)/$READS_DIR/HG002_R2.fastq.gz"
}

case "${1:-}" in
    --refs)  refs ;;
    --truth) truth ;;
    --reads) reads ;;
    --all)   refs; truth; reads ;;
    *) echo "usage: $0 [--refs|--truth|--reads|--all]"; exit 1 ;;
esac

echo ""
echo "Done. snpEff will download its GRCh38 database on first run; to pre-cache:"
echo "    docker run --rm -v \$PWD:/data staphb/snpeff:5.1 snpEff download GRCh38.105"
