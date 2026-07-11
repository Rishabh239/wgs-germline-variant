process MULTIQC {
    label 'process_low'
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path '*'
    path config

    output:
    path "*.html", emit: report
    path "*_data", emit: data, optional: true

    script:
    """
    multiqc -f -c ${config} .
    """
}
