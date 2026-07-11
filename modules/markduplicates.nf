process MARKDUPLICATES {
    tag   "$meta.id"
    label 'process_medium'
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "*.md.bam*"
    publishDir "${params.outdir}/qc/markduplicates", mode: 'copy', pattern: "*.metrics.txt"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.md.bam"), path("${meta.id}.md.bam.bai"), emit: bam
    path "${meta.id}.md.metrics.txt", emit: metrics

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" MarkDuplicates \\
        -I ${bam} \\
        -O ${meta.id}.md.bam \\
        -M ${meta.id}.md.metrics.txt \\
        --CREATE_INDEX true
    mv ${meta.id}.md.bai ${meta.id}.md.bam.bai
    """
}
