process SNPEFF {
    label 'process_medium'
    container 'quay.io/biocontainers/snpeff:5.1--hdfd78af_2'
    // When a pre-downloaded cache is provided, bind-mount it into the container
    // at the same path so snpEff can read it without any in-container networking.
    containerOptions { params.snpeff_cache ? "-v ${params.snpeff_cache}:${params.snpeff_cache}" : '' }
    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    val db

    output:
    path "cohort.ann.vcf", emit: vcf
    path "snpEff.csv"    , emit: csv
    path "*.html"        , emit: report, optional: true
    path "*genes.txt"    , emit: genes,  optional: true

    script:
    // With --snpeff_cache set, use the mounted DB and -nodownload (no network).
    // Without it, snpEff downloads the DB into a writable workdir path (needs
    // in-container internet).
    def cache = params.snpeff_cache ? "-dataDir ${params.snpeff_cache} -nodownload" : '-dataDir \$PWD/snpeff_data'
    """
    mkdir -p snpeff_data
    snpEff -Xmx${task.memory.toGiga()}g ${cache} \\
        -csvStats snpEff.csv \\
        ${db} ${vcf} > cohort.ann.vcf
    """
}
