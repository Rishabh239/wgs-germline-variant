process SAMTOOLS_SORT {
    tag   "$meta.id"
    label 'process_medium'
    container 'staphb/samtools:1.17'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam
    path "${meta.id}.flagstat.txt", emit: flagstat

    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam ${sam}
    samtools index ${meta.id}.sorted.bam
    samtools flagstat ${meta.id}.sorted.bam > ${meta.id}.flagstat.txt
    """
}
