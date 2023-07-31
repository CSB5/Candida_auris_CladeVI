/*
 * A simple stat for nanpore fastq file
 */
process NANOPORE_SEQKIT_STATS {
    tag "Calculate stats for $read_ch"
    publishDir params.outdir, mode:'copy'
    container 'quality_control:0.0.1'

    input:
    path read_ch
    val suffix

    output:
    path "${params.sample_id}.nanopore.${suffix}.seqkit-stats.txt"

    script:
    """
    seqkit stats $read_ch -a -T > ${params.sample_id}.nanopore.${suffix}.seqkit-stats.txt
    """
    
}

/*
 * A simple stat for illumina paired-end fastq file
 */
process ILLUMINA_SEQKIT_STATS {
    tag "Calculate stats for $read_ch1 and $read_ch2"
    publishDir params.outdir, mode:'copy'

    input:
    path read_ch1
    path read_ch2
    val suffix

    output:
    path "${params.sample_id}.illumina-pe.${suffix}.seqkit-stats.txt"

    script:
    """
    seqkit stats $read_ch1 $read_ch2 -a -T > ${params.sample_id}.illumina-pe.${suffix}.seqkit-stats.txt
    """

}

/*
 * Adapter trimming and quality filtering for illumina paired-end fastq file
 */
process ILLUMINA_FASTP {
    tag "Run adapter trimming and quality filtering for Illumina"
    publishDir params.outdir, mode:'copy'

    input:
    path read_ch1
    path read_ch2

    output:
    path "${params.sample_id}.illumina-pe.fastp.json"
    path "${params.sample_id}.illumina-pe.fastp.html"
    path "${params.sample_id}.trim.qc.1.fastq"
    path "${params.sample_id}.trim.qc.2.fastq"

    script:
    """
    fastp --in1 $read_ch1 --in2 $read_ch2 \
    --out1 ${params.sample_id}.trim.qc.1.fastq --out2 ${params.sample_id}.trim.qc.2.fastq \
    --qualified_quality_phred 18 --unqualified_percent_limit 40 --average_qual 18 --n_base_limit 5 \
    --length_required 100 \
    --low_complexity_filter --complexity_threshold 30 \
    --json ${params.sample_id}.illumina-pe.fastp.json \
    --html ${params.sample_id}.illumina-pe.fastp.html
    """

}

/*
 * Sequencing QC for illumina paired-end fastq file after adapter trimming and quality filtering.
 */
process ILLUMINA_FASTQC {
    tag "Run FastQC for Illumina"
    publishDir params.outdir, mode:'copy'

    input:
    path read_ch1
    path read_ch2

    output:
    path "${params.sample_id}.illumina-pe.fastqc"

    script:
    """
    mkdir ${params.sample_id}.illumina-pe.fastqc
    fastqc $read_ch1 $read_ch2 --outdir ${params.sample_id}.illumina-pe.fastqc
    """

}

/*
 * Quality filtering for illumina paired-end fastq file.
 */
process NANOPORE_CHOPPER {
    tag "Run Chopper for Nanopore"
    publishDir params.outdir, mode:'copy'
    container 'quality_control_chopper:0.0.1'

    input:
    path read_ch

    output:
    path "${params.sample_id}.trim.qc.fastq"

    script:
    """
    gunzip -c $read_ch | chopper --quality 10 --minlength 500 --threads 16 > ${params.sample_id}.trim.qc.fastq
    """

}

/*
 * Sequencing QC for nanopore fastq file after adapter trimming and quality filtering
 */
process NANOPORE_NANOSTAT {
    tag "Run NanoStat for Nanopore"
    publishDir params.outdir, mode:'copy'
    container 'quality_control:0.0.1'

    input:
    path read_ch

    output:
    path "${params.sample_id}.nanopore.nanostat.txt"

    script:
    """
    NanoStat --fastq $read_ch --name ${params.sample_id}.nanopore.nanostat.txt
    """

}
