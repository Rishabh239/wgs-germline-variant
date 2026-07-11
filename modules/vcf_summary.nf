process VCF_SUMMARY {
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/qc/vcf_summary", mode: 'copy'

    input:
    path vcf

    output:
    path "*.summary.tsv"           , emit: tsv
    path "*_vcf_summary_mqc.json"  , emit: mqc

    script:
    // vcf_summary.py is auto-staged from bin/ onto PATH by Nextflow.
    """
    vcf_summary.py --vcf ${vcf} --prefix cohort
    """
}
