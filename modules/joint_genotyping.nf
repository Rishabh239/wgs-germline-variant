process COMBINE_GVCFS {
    label 'process_medium'
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    path gvcfs
    path tbis
    tuple path(fasta), path(fai), path(dict)
    val intervals

    output:
    tuple path("cohort.g.vcf.gz"), path("cohort.g.vcf.gz.tbi"), emit: gvcf

    script:
    def variants = gvcfs.collect { "-V ${it}" }.join(' ')
    def ival = intervals ? "-L ${intervals}" : ''
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" CombineGVCFs \\
        -R ${fasta} ${variants} ${ival} \\
        -O cohort.g.vcf.gz
    """
}

process GENOTYPE_GVCFS {
    label 'process_medium'
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple path(gvcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)
    tuple path(dbsnp), path(dbsnp_idx)
    val intervals

    output:
    tuple path("cohort.vcf.gz"), path("cohort.vcf.gz.tbi"), emit: vcf

    script:
    def ival = intervals ? "-L ${intervals}" : ''
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenotypeGVCFs \\
        -R ${fasta} -V ${gvcf} -D ${dbsnp} ${ival} \\
        -O cohort.vcf.gz
    """
}
