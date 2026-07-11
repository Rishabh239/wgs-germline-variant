process HAPPY {
    label 'process_medium'
    container 'jmcdani20/hap.py:v0.3.12'
    publishDir "${params.outdir}/benchmark", mode: 'copy'

    input:
    path query_vcf
    tuple path(truth_vcf), path(truth_tbi)
    path truth_bed
    tuple path(fasta), path(fai), path(dict)

    output:
    path "happy.*"          , emit: results
    path "happy.summary.csv", emit: summary

    script:
    // Precision/recall/F1 for SNPs and INDELs, restricted to the GIAB
    // high-confidence regions (-f). Default xcmp engine only needs the FASTA.
    """
    /opt/hap.py/bin/hap.py \\
        ${truth_vcf} ${query_vcf} \\
        -f ${truth_bed} \\
        -r ${fasta} \\
        -o happy \\
        --threads ${task.cpus}
    """
}
