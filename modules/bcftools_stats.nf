process BCFTOOLS_STATS {
    label 'process_low'
    container 'staphb/bcftools:1.17'
    publishDir "${params.outdir}/qc/bcftools", mode: 'copy'

    input:
    path vcf

    output:
    path "${vcf.baseName}.stats.txt", emit: stats

    script:
    """
    bcftools stats ${vcf} > ${vcf.baseName}.stats.txt
    """
}
