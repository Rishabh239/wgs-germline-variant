process FASTP {
    tag   "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'
    publishDir "${params.outdir}/qc/fastp", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: reads
    path "*.fastp.json", emit: json
    path "*.fastp.html", emit: html

    script:
    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${meta.id}_R1.trim.fastq.gz -O ${meta.id}_R2.trim.fastq.gz \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    """
}
