process BWA_MEM {
    tag   "$meta.id"
    label 'process_high'
    container 'staphb/bwa:0.7.17'

    input:
    tuple val(meta), path(reads)
    path fasta            // reference FASTA (index prefix)
    path index            // bwa index files: .amb .ann .bwt .pac .sa (.alt)

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam

    script:
    // Read group is required by GATK. Alt-aware mapping is automatic when
    // ${fasta}.alt is present alongside the index.
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}"
    """
    bwa mem -t ${task.cpus} -M -R "${rg}" ${fasta} ${reads[0]} ${reads[1]} > ${meta.id}.sam
    """
}
