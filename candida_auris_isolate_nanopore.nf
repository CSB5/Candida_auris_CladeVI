/*
 * Candida auris WGS pipeline
 * ==============================
 */

/*
 * Common parameters
 */
params.outdir = "/home/ubuntu/data/sgh_candida_auris/20230604_analysis_results"
params.sample_id = "999999999"
params.min_contig_length = "1000"
params.genome_size = "12500000"
params.genus = "Candida"

params.nanopore_reads = "/home/ubuntu/data/sgh_candida_auris/N00466_BC01.fastq.gz"
params.medaka_model = "r941_min_sup_g507"
params.nanopore_min_target_depth = 30
params.nanopore_max_target_depth = 150
params.nanopore_min_target_coverage = 95

/*
 * Sequence QC
 * ===========
 */
include { NANOPORE_CHOPPER } from './modules/quality_control'
include { NANOPORE_NANOSTAT } from './modules/quality_control'

/*
 * Decontamination
 * ===============
 */
include { NANOPORE_DECONTAMINATION_MINIMAP2 } from './modules/quality_control'
include { NANOPORE_SEQKIT_STATS as NANOPORE_SEQKIT_STATS_BEFORE} from './modules/quality_control'
include { NANOPORE_SEQKIT_STATS as NANOPORE_SEQKIT_STATS_AFTER} from './modules/quality_control'

/*
 * Genome assembly
 * ===============
 */
include { NANOPORE_CHECK_DEPTH } from './modules/assembly'
include { FLYE_ONT_HQ } from './modules/assembly'
include { MEDAKA_X2 } from './modules/assembly'
include { CLEANUP_CONSENSUS as CLEANUP_CONSENSUS_NANOPORE } from './modules/assembly'
include { QUAST_NANOPORE } from './modules/assembly'

/*
 * Logging
 * =======
 */

log.info """\
    C. auris WGS PIPELINE
    ==========================
    nanopore_reads         : ${params.nanopore_reads}
    ilumina_reads_1        : ${params.illumina_reads_1}
    ilumina_reads_2        : ${params.illumina_reads_2}
    outdir                 : ${params.outdir}
    """
    .stripIndent()

workflow {

    // Input files
    // ===========
    read_ch = Channel.fromPath(params.nanopore_reads, checkIfExists: true)

    read_trim_qc_ch = NANOPORE_CHOPPER(read_ch)
    nanostat_report = NANOPORE_NANOSTAT(read_trim_qc_ch)
    
    // Decontamination
    // ===============
    nanopore_stats_ch = NANOPORE_SEQKIT_STATS_BEFORE(read_trim_qc_ch, "no_decont")

    // Subsampling reads to a specific depth
    // =====================================
    (nanopore_check_depth_log, read_for_downstream_ch) = NANOPORE_CHECK_DEPTH(read_trim_qc_ch, nanopore_stats_ch)

    // Genome assembly and assembly QC
    // ===============================
    flye_assembly_ch = FLYE_ONT_HQ(read_for_downstream_ch)
    (nanopore_polish_ch, nanopore_polish_round1, nanopore_polish_round2, nanopore_flye_medakax2_ch) = MEDAKA_X2(read_for_downstream_ch, flye_assembly_ch)
    nanopore_assembly_ch = CLEANUP_CONSENSUS_NANOPORE(nanopore_flye_medakax2_ch, "nanopore.flye.medaka_x2")

    // Assembly QC with QUAST
    nanopore_quast_ch = QUAST_NANOPORE(nanopore_assembly_ch, read_for_downstream_ch, "nanopore.flye.medaka_x2")
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Check the result at $params.outdir\n" : "Oops .. something went wrong" )
}
