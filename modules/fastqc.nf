process FASTQC {
    tag   "$meta.id"
    label 'process_low'
    container 'staphb/fastqc:0.12.1'
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "*.zip" , emit: zip
    path "*.html", emit: html

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}
