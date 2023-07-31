/*
 * Candida auris WGS pipeline
 * ==========================
 */

/*
 * Sequence QC
 * ===========
 */
include { ILLUMINA_FASTP } from './modules/quality_control'
include { ILLUMINA_FASTQC } from './modules/quality_control'
include { ILLUMINA_SEQKIT_STATS as ILLUMINA_SEQKIT_STATS_BEFORE } from './modules/quality_control'
include { ILLUMINA_CHECK_DEPTH } from './modules/assembly'

/*
 * Genome assembly
 * ===============
 */
include { PILON } from './modules/assembly'
include { CLEANUP_CONSENSUS as CLEANUP_CONSENSUS_HYBRID_PILON } from './modules/assembly'
include { QUAST_HYBRID } from './modules/assembly'

/*
 * Logging
 * =======
 */

log.info """\
    C. auris WGS PIPELINE
    ==========================
    nanopore_contigs                : ${params.nanopore_contigs}
    nanopore_read_for_downstream_ch : ${params.nanopore_reads}
    ilumina_reads_1                 : ${params.illumina_reads_1}
    ilumina_reads_2                 : ${params.illumina_reads_2}
    outdir                          : ${params.outdir}
    """
    .stripIndent()

workflow {

    // Input files
    // ===========
    nanopore_read_for_downstream_ch = Channel.fromPath(params.nanopore_reads, checkIfExists: true)
    nanopore_contig_ch = Channel.fromPath(params.nanopore_contigs, checkIfExists: true)
    read_ch1 = Channel.fromPath(params.illumina_reads_1, checkIfExists: true)
    read_ch2 = Channel.fromPath(params.illumina_reads_2, checkIfExists: true)

    // Adapter trimming and sequencing QC
    // ==================================
    (fastp_json, fastp_html, read_trim_qc_ch1, read_trim_qc_ch2) = ILLUMINA_FASTP(read_ch1, read_ch2)
    ILLUMINA_FASTQC(read_trim_qc_ch1, read_trim_qc_ch2)
    
    // Decontamination
    // ===============
    illumina_stats_ch = ILLUMINA_SEQKIT_STATS_BEFORE(read_trim_qc_ch1, read_trim_qc_ch2, "no_decont")

    // Subsampling reads to a specific depth
    // =====================================
    (illumina_check_depth_log, read_for_downstream_ch1, read_for_downstream_ch2) = ILLUMINA_CHECK_DEPTH(read_trim_qc_ch1, read_trim_qc_ch2, illumina_stats_ch)

    // Genome polishing
    // ===============================
    (pilon_ch, pilon_change_ch) = PILON(nanopore_contig_ch, read_for_downstream_ch1, read_for_downstream_ch2)

    // Assembly QC with QUAST
    nanopore_quast_ch = QUAST_HYBRID(pilon_ch, nanopore_read_for_downstream_ch, read_for_downstream_ch1, read_for_downstream_ch2, "hybrid.pilon")

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Check the result at $params.outdir\n" : "Oops .. something went wrong" )
}
