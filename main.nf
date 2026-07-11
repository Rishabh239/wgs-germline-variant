#!/usr/bin/env nextflow
/*
========================================================================================
    wgs-germline-nf : GATK4 best-practices germline short-variant pipeline
    FASTQ -> aligned BAM -> joint-genotyped, hard-filtered, annotated VCF -> QC/benchmark
    (rule-based hard-filtering; no machine learning)
========================================================================================
*/
nextflow.enable.dsl = 2

include { FASTQC          } from './modules/fastqc.nf'
include { FASTP           } from './modules/fastp.nf'
include { BWA_MEM         } from './modules/bwa_mem.nf'
include { SAMTOOLS_SORT   } from './modules/samtools_sort.nf'
include { MARKDUPLICATES  } from './modules/markduplicates.nf'
include { BASERECALIBRATOR; APPLYBQSR } from './modules/bqsr.nf'
include { HAPLOTYPECALLER } from './modules/haplotypecaller.nf'
include { COMBINE_GVCFS; GENOTYPE_GVCFS } from './modules/joint_genotyping.nf'
include { HARD_FILTER     } from './modules/hard_filter.nf'
include { SNPEFF          } from './modules/snpeff.nf'
include { BCFTOOLS_STATS  } from './modules/bcftools_stats.nf'
include { VCF_SUMMARY     } from './modules/vcf_summary.nf'
include { HAPPY           } from './modules/happy.nf'
include { MULTIQC         } from './modules/multiqc.nf'

// ----------------------------------------------------------------------------
// Helper: fail early with a readable message if a reference file is missing.
// ----------------------------------------------------------------------------
def req(path, what) {
    def f = file(path)
    if (!f.exists()) {
        exit 1, "Missing ${what}: ${path}\n  -> run ./download_data.sh or set the matching --param."
    }
    return f
}

// Pick the right index extension: .tbi for bgzipped VCFs, .idx for plain VCFs.
def vcfIndex(vcf) {
    return vcf.toString().endsWith('.gz') ? "${vcf}.tbi" : "${vcf}.idx"
}

workflow {

    // ---- reference channels (staged together where GATK needs companions) ----
    ch_ref = Channel.value([
        req(params.fasta, 'reference FASTA'),
        req(params.fasta + '.fai', 'FASTA .fai index'),
        req(params.dict, 'sequence dictionary (.dict)')
    ])
    ch_fasta = Channel.value(req(params.fasta, 'reference FASTA'))
    ch_bwa   = Channel.value(files("${params.fasta}.{amb,ann,bwt,pac,sa,alt}"))
    ch_dbsnp = Channel.value([
        req(params.dbsnp, 'dbSNP VCF'),
        req(vcfIndex(params.dbsnp), 'dbSNP index')
    ])
    // Known-sites for BQSR. Deduplicate by path so the same file is never
    // staged twice (Nextflow rejects that as a name collision); the test
    // profile may legitimately reuse one file for several roles.
    known_paths   = [ params.dbsnp, params.mills, params.known_indels ].unique()
    ch_known_vcfs = Channel.value( known_paths.collect { req(it, 'known-sites VCF') } )
    ch_known_idxs = Channel.value( known_paths.collect { req(vcfIndex(it), 'known-sites index') } )
    // Value channel for the optional -L interval. Must never be null (a null
    // process input is rejected), so fall back to '' when no interval is set.
    ch_intervals = Channel.value(params.intervals ?: '')

    // ---- input samplesheet -> [ meta, [r1, r2] ] ----
    ch_reads = Channel.fromPath(req(params.input, 'samplesheet'))
        .splitCsv(header: true)
        .map { row ->
            tuple(
                [ id: row.sample ],
                [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]
            )
        }

    // ---- read QC (+ optional trimming) ----
    FASTQC(ch_reads)
    ch_align_in = params.trim ? FASTP(ch_reads).reads : ch_reads

    // ---- align -> sort -> mark duplicates -> BQSR ----
    BWA_MEM(ch_align_in, ch_fasta, ch_bwa)
    SAMTOOLS_SORT(BWA_MEM.out.sam)
    MARKDUPLICATES(SAMTOOLS_SORT.out.bam)

    BASERECALIBRATOR(MARKDUPLICATES.out.bam, ch_ref, ch_known_vcfs, ch_known_idxs, ch_intervals)
    ch_bqsr_in = MARKDUPLICATES.out.bam.join(BASERECALIBRATOR.out.table)   // join on meta
    APPLYBQSR(ch_bqsr_in, ch_ref)

    // ---- per-sample calling -> joint genotyping ----
    HAPLOTYPECALLER(APPLYBQSR.out.bam, ch_ref, ch_dbsnp, ch_intervals)
    ch_gvcfs = HAPLOTYPECALLER.out.gvcf.map { it[1] }.collect()
    ch_tbis  = HAPLOTYPECALLER.out.gvcf.map { it[2] }.collect()
    COMBINE_GVCFS(ch_gvcfs, ch_tbis, ch_ref, ch_intervals)
    GENOTYPE_GVCFS(COMBINE_GVCFS.out.gvcf, ch_ref, ch_dbsnp, ch_intervals)

    // ---- rule-based hard-filtering -> (optional) annotation ----
    HARD_FILTER(GENOTYPE_GVCFS.out.vcf, ch_ref)

    if (params.skip_snpeff) {
        ch_final_vcf  = HARD_FILTER.out.vcf.map { it[0] }   // just the .vcf.gz
        ch_snpeff_mqc = Channel.empty()
    } else {
        SNPEFF(HARD_FILTER.out.vcf, params.snpeff_db)
        ch_final_vcf  = SNPEFF.out.vcf
        ch_snpeff_mqc = SNPEFF.out.csv
    }

    // ---- QC on the final callset ----
    BCFTOOLS_STATS(ch_final_vcf)
    VCF_SUMMARY(ch_final_vcf)

    // ---- optional GIAB truth-set benchmarking (hap.py) ----
    if (params.benchmark) {
        ch_truth = Channel.value([
            req(params.truth_vcf, 'GIAB truth VCF'),
            req(params.truth_vcf + '.tbi', 'GIAB truth .tbi')
        ])
        ch_truth_bed = Channel.value(req(params.truth_bed, 'GIAB high-confidence BED'))
        HAPPY(HARD_FILTER.out.vcf.map { it[0] }, ch_truth, ch_truth_bed, ch_ref)
    }

    // ---- aggregate everything into one MultiQC report ----
    ch_multiqc = FASTQC.out.zip
        .mix(SAMTOOLS_SORT.out.flagstat)
        .mix(MARKDUPLICATES.out.metrics)
        .mix(BCFTOOLS_STATS.out.stats)
        .mix(ch_snpeff_mqc)
        .mix(VCF_SUMMARY.out.mqc)
        .collect()
    MULTIQC(ch_multiqc, file(params.multiqc_config))
}
