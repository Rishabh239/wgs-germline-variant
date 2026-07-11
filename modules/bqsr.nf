process BASERECALIBRATOR {
    tag   "$meta.id"
    label 'process_medium'
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)
    path known_vcfs
    path known_idxs
    val intervals

    output:
    tuple val(meta), path("${meta.id}.recal.table"), emit: table

    script:
    def ks = known_vcfs.collect { "--known-sites ${it}" }.join(' ')
    def ival = intervals ? "-L ${intervals}" : ''
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" BaseRecalibrator \\
        -I ${bam} -R ${fasta} \\
        ${ks} \\
        ${ival} \\
        -O ${meta.id}.recal.table
    """
}

process APPLYBQSR {
    tag   "$meta.id"
    label 'process_medium'
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(table)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.id}.recal.bam"), path("${meta.id}.recal.bam.bai"), emit: bam

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" ApplyBQSR \\
        -I ${bam} -R ${fasta} \\
        --bqsr-recal-file ${table} \\
        -O ${meta.id}.recal.bam
    mv ${meta.id}.recal.bai ${meta.id}.recal.bam.bai
    """
}
