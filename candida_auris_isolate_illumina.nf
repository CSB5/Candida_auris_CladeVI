/*
 * Candida auris WGS pipeline
 * ==========================
 */

params.snippy_ref = "/mnt/data/discovery_of_cauris/data/C_auris_B8441_current_chromosomes.fasta"
params.snippy_cauris_mito_ref = "/mnt/data/discovery_of_cauris/data/MT849287.1.fasta"
params.snpeff_db = "candida_auris_GCA_002759435.2"
params.snpeff_path = "/home/download/snpEff/snpEff.jar"
params.sra_accession = false

/*
 * Sequence QC
 * ===========
 */
include { ILLUMINA_FASTP } from './modules/quality_control'
include { ILLUMINA_FASTQC } from './modules/quality_control'

/*
 * Decontamination
 * ===============
 */
include { ILLUMINA_SEQKIT_STATS } from './modules/quality_control'

/*
 * Genome assembly
 * ===============
 */
include { ILLUMINA_CHECK_DEPTH } from './modules/assembly'
include { SPADES } from './modules/assembly'
include { CLEANUP_CONSENSUS as CLEANUP_CONSENSUS_ILLUMINA } from './modules/assembly'
include { QUAST_ILLUMINA } from './modules/assembly'

/*
 * Genome characterization
 * =======================
 */
include { FUNGAL_SPECIES_IDEN_MASH_READ } from './modules/species_identification'
include { FUNGAL_SPECIES_IDEN_MASH_CONTIG } from './modules/species_identification'

include { VARIENT_CALLING_SNIPPY as VARIENT_CALLING_SNIPPY_GENOME } from './modules/genome_characterization/snps_analysis'
include { VARIENT_CALLING_SNIPPY as VARIENT_CALLING_SNIPPY_CAURIS_MITO } from './modules/genome_characterization/snps_analysis'
include { VARIENT_ANNOTATION_SNPEFF } from './modules/genome_characterization/snps_analysis'

include { GET_FASTQ_FROM_SRA } from './modules/utils'

/*
 * Logging
 * =======
 */

log.info """\
    C. auris WGS PIPELINE
    ==========================
    sra_accession          : ${params.sra_accession}
    outdir                 : ${params.outdir}
    """
    .stripIndent()

workflow {

    // Input files
    // ===========
    (read_ch1, read_ch2) = GET_FASTQ_FROM_SRA(params.sra_accession)

    // Adapter trimming and sequencing QC
    (fastp_json, fastp_html, read_trim_qc_ch1, read_trim_qc_ch2) = ILLUMINA_FASTP(read_ch1, read_ch2)
    ILLUMINA_FASTQC(read_trim_qc_ch1, read_trim_qc_ch2)
    
    // Decontamination
    // ===============
    illumina_stats_ch = ILLUMINA_SEQKIT_STATS(read_trim_qc_ch1, read_trim_qc_ch2, "no_decont")

    // Species identificiation
    // =======================
    species_illumina_read_ch = FUNGAL_SPECIES_IDEN_MASH_READ(read_trim_qc_ch1, read_trim_qc_ch2)

    // Subsampling reads to a specific depth
    // =====================================
    (illumina_check_depth_log, read_for_downstream_ch1, read_for_downstream_ch2) = ILLUMINA_CHECK_DEPTH(read_trim_qc_ch1, read_trim_qc_ch2, illumina_stats_ch)

    // Genome assembly and assembly QC
    // ===============================
    (spade_ch, illumina_spades_ch) = SPADES(read_for_downstream_ch1, read_for_downstream_ch2)
    illumina_assembly_ch = CLEANUP_CONSENSUS_ILLUMINA(illumina_spades_ch, "illumina.spades")

    // Genome annotation and QC
    // ========================
    illumina_quast_ch = QUAST_ILLUMINA(illumina_assembly_ch, read_for_downstream_ch1, read_for_downstream_ch2, 'illumina.spades')
    species_illumina_contig_ch = FUNGAL_SPECIES_IDEN_MASH_CONTIG(illumina_assembly_ch)

    snippy_vcf_ch = VARIENT_CALLING_SNIPPY_GENOME(params.snippy_ref, read_for_downstream_ch1, read_for_downstream_ch2, 'snippy')
    VARIENT_ANNOTATION_SNPEFF(snippy_vcf_ch, params.snpeff_db, 'snippy')

    snippy_cauris_mito_vcf_ch = VARIENT_CALLING_SNIPPY_CAURIS_MITO(params.snippy_cauris_mito_ref, read_for_downstream_ch1, read_for_downstream_ch2, 'cauris_mito.snippy')

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Check the result at $params.outdir\n" : "Oops .. something went wrong" )
}
