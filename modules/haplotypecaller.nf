process HAPLOTYPECALLER {
    tag   "$meta.id"
    label 'process_high'
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/variants/gvcf", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)
    tuple path(dbsnp), path(dbsnp_idx)
    val intervals

    output:
    tuple val(meta), path("${meta.id}.g.vcf.gz"), path("${meta.id}.g.vcf.gz.tbi"), emit: gvcf

    script:
    def ival = intervals ? "-L ${intervals}" : ''
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller \\
        -R ${fasta} -I ${bam} \\
        -ERC GVCF \\
        -D ${dbsnp} \\
        ${ival} \\
        -O ${meta.id}.g.vcf.gz
    """
}
