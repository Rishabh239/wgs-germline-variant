process HARD_FILTER {
    label 'process_low'
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple path("cohort.filtered.vcf.gz"), path("cohort.filtered.vcf.gz.tbi"), emit: vcf

    script:
    // GATK-recommended germline hard filters (rule-based; NOT VQSR/ML).
    """
    # --- SNPs ---
    gatk SelectVariants -R ${fasta} -V ${vcf} --select-type-to-include SNP -O snps.vcf.gz
    gatk VariantFiltration -R ${fasta} -V snps.vcf.gz -O snps.filt.vcf.gz \\
        --filter-expression "QD < 2.0"              --filter-name "QD2" \\
        --filter-expression "QUAL < 30.0"           --filter-name "QUAL30" \\
        --filter-expression "SOR > 3.0"             --filter-name "SOR3" \\
        --filter-expression "FS > 60.0"             --filter-name "FS60" \\
        --filter-expression "MQ < 40.0"             --filter-name "MQ40" \\
        --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

    # --- INDELs ---
    gatk SelectVariants -R ${fasta} -V ${vcf} --select-type-to-include INDEL -O indels.vcf.gz
    gatk VariantFiltration -R ${fasta} -V indels.vcf.gz -O indels.filt.vcf.gz \\
        --filter-expression "QD < 2.0"               --filter-name "QD2" \\
        --filter-expression "QUAL < 30.0"            --filter-name "QUAL30" \\
        --filter-expression "FS > 200.0"             --filter-name "FS200" \\
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \\
        --filter-expression "SOR > 10.0"             --filter-name "SOR10"

    # --- merge back into one callset ---
    gatk MergeVcfs -I snps.filt.vcf.gz -I indels.filt.vcf.gz -O cohort.filtered.vcf.gz
    """
}
