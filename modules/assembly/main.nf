/*
 * Check sequencing depth and perform subsampling if needed (for nanopore fastq)
 */
process NANOPORE_CHECK_DEPTH {
    tag "Check sequencing depth and perform subsampling if needed"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )
    container 'quality_control:0.0.1'

    input:
    path read_ch
    path nanopore_after_decont_stats_ch

    output:
    path "${params.sample_id}.nanopore.check_depth.log"
    path "${params.sample_id}_for_downstream.fastq"

    script:
    """
    TOTAL_BASES=\$(tail -n +2 $nanopore_after_decont_stats_ch | cut -f 5 | awk -v OFMT='%d' '{s+=\$1} END {print s}')

    CHECK_DEPTH_OUTPUT=\$(python /home/scripts/check_depth.py --total_bases \$TOTAL_BASES --genome_size ${params.genome_size} \
    --min_target_depth ${params.nanopore_min_target_depth} --max_target_depth ${params.nanopore_max_target_depth})

    echo -e \$CHECK_DEPTH_OUTPUT > "${params.sample_id}.nanopore.check_depth.log"
    PROPORTION=\$(cat ${params.sample_id}.nanopore.check_depth.log | tail -n 1)
    cat $read_ch | seqkit sample -s 100 -p \$PROPORTION -o ${params.sample_id}_for_downstream.fastq
    """
}

/*
 * Check sequencing depth and perform subsampling if needed (for illumina paired-read fastq)
 */
process ILLUMINA_CHECK_DEPTH {
    tag "Check sequencing depth and perform subsampling if needed"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )

    input:
    path read_ch1
    path read_ch2
    path illumina_after_decont_stats_ch

    output:
    path "${params.sample_id}.illumina-pe.check_depth.log"
    path "${params.sample_id}.for_downstream.1.fastq"
    path "${params.sample_id}.for_downstream.2.fastq"

    script:
    """
    TOTAL_BASES=\$(tail -n +2 $illumina_after_decont_stats_ch | cut -f 5 | awk -v OFMT='%d' '{s+=\$1} END {print s}')

    CHECK_DEPTH_OUTPUT=\$(python /home/scripts/check_depth.py --total_bases \$TOTAL_BASES --genome_size ${params.genome_size} \
    --min_target_depth ${params.illumina_min_target_depth} --max_target_depth ${params.illumina_max_target_depth})

    echo -e \$CHECK_DEPTH_OUTPUT > "${params.sample_id}.illumina-pe.check_depth.log"
    PROPORTION=\$(cat ${params.sample_id}.illumina-pe.check_depth.log | tail -n 1)
    cat $read_ch1 | seqkit sample -s 100 -p \$PROPORTION -o ${params.sample_id}.for_downstream.1.fastq
    cat $read_ch2 | seqkit sample -s 100 -p \$PROPORTION -o ${params.sample_id}.for_downstream.2.fastq
    """
    
}

/*
 * Genome assembly using Fyle for nanopore reads
 */
process FLYE_ONT_HQ {
    tag "Using flye to assemble $read_ch"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )
    container 'assembly:0.0.1'

    input:
    path read_ch

    output:
    path "${params.sample_id}.nanopore.flye"

    script:
    """
    flye --nano-hq $read_ch --out-dir ${params.sample_id}.nanopore.flye --threads $task.cpus --genome-size $params.genome_size
    """
}

/*
 * Genome polishing for nanopore assembly (from Flye) using medaka (run twice)
 */
process MEDAKA_X2 {
    tag "Using medaka to polish $read_ch"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )
    container 'assembly:0.0.1'

    input:
    path read_ch
    path flye_assembly_ch

    output:
    path "${params.sample_id}.nanopore.medaka_x2"
    path "${params.sample_id}.nanopore.medaka_x2/medaka_1"
    path "${params.sample_id}.nanopore.medaka_x2/medaka_2"
    path "${params.sample_id}.nanopore.medaka_x2/medaka_2/consensus.fasta"

    script:
    """
    mkdir -p ${params.sample_id}.nanopore.medaka_x2/medaka_1
    mkdir -p ${params.sample_id}.nanopore.medaka_x2/medaka_2
    medaka_consensus -i $read_ch -d $flye_assembly_ch/assembly.fasta -o ${params.sample_id}.nanopore.medaka_x2/medaka_1 -t $task.cpus -m $params.medaka_model
    medaka_consensus -i $read_ch -d ${params.sample_id}.nanopore.medaka_x2/medaka_1/consensus.fasta -o ${params.sample_id}.nanopore.medaka_x2/medaka_2 -t $task.cpus -m $params.medaka_model
    """
}

/*
 * Genome assembly using SPAdes for illumina paired reads
 */
process SPADES {
    tag "Genome assembly using SPAdes for illumina paired reads"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )

    input:
    path read_ch1
    path read_ch2

    output:
    path "${params.sample_id}.illumina-pe.spades"
    path "${params.sample_id}.illumina-pe.spades/scaffolds.fasta"

    script:
    """
    /home/download/SPAdes-3.15.5-Linux/bin/spades.py \
    --pe1-1 $read_ch1 --pe1-2 $read_ch2 \
    -o ${params.sample_id}.illumina-pe.spades --isolate -t $task.cpus
    """
    
}

/*
 * Genome polishing
 */
process PILON {
    tag "Using pilon to polish long read assembly with illumina reads"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )
    container 'assembly:0.0.1'

    input:
    path assembly_ch
    path read_ch1
    path read_ch2

    output:
    path "${params.sample_id}.hybrid.pilon.fasta"
    path "${params.sample_id}.hybrid.pilon.changes"

    script:
    """
    bwa index $assembly_ch
    bwa mem -t $task.cpus $assembly_ch $read_ch1 $read_ch2 | samtools view -b - | samtools sort > ${params.sample_id}.hybrid.pilon.sorted.bam
    samtools index ${params.sample_id}.hybrid.pilon.sorted.bam
    java -Xmx16G -jar /home/download/pilon-1.24.jar --genome $assembly_ch --changes \
    --frags ${params.sample_id}.hybrid.pilon.sorted.bam --threads $task.cpus --output ${params.sample_id}.hybrid.pilon
    """
}

/*
 * Clean up the consensus by removing contigs that are shorter than "min_contig_length"
 */
process CLEANUP_CONSENSUS {
    tag "Removing contigs that have length < $params.min_contig_length"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )

    input:
    path assembly_ch
    val suffix

    output:
    path "${params.sample_id}.${suffix}.fasta"

    script: 
    """
    seqkit seq --min-len $params.min_contig_length $assembly_ch > ${params.sample_id}.${suffix}.fasta
    """

}

/*
 * QC for genome assembly. We will also evaluate genomic features (Prokka's output)
 */
process QUAST_NANOPORE {
    tag "Checking assembly QC with QUAST"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )
    container 'staphb/quast:5.0.2'

    input:
    path assembly_ch
    path read_ch
    val suffix

    output:
    path "${params.sample_id}.${suffix}.quast"

    script:
    """
    quast.py $assembly_ch --nanopore $read_ch \
    --threads $task.cpus --labels ${params.sample_id} -o ${params.sample_id}.${suffix}.quast
    """
}

/*
 * QC for genome assembly. We will also evaluate genomic features (Prokka's output)
 */
process QUAST_ILLUMINA {
    tag "Checking assembly QC with QUAST"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )

    input:
    path assembly_ch
    path read_ch1
    path read_ch2
    val suffix

    output:
    path "${params.sample_id}.${suffix}.quast"

    script:
    """
    quast.py $assembly_ch --pe1 $read_ch1 --pe2 $read_ch2 \
    --threads $task.cpus --labels ${params.sample_id} -o ${params.sample_id}.${suffix}.quast
    """
}

/*
 * QC for genome assembly. We will also evaluate genomic features (Prokka's output)
 */
process QUAST_HYBRID {
    tag "Checking assembly QC with QUAST"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )
    container 'staphb/quast:5.0.2'

    input:
    path assembly_ch
    path read_ch
    path read_ch1
    path read_ch2
    val suffix

    output:
    path "${params.sample_id}.${suffix}.quast"

    script:
    """
    quast.py $assembly_ch --pe1 $read_ch1 --pe2 $read_ch2 --nanopore $read_ch --k-mer-stats \
    --threads $task.cpus --labels ${params.sample_id} -o ${params.sample_id}.${suffix}.quast
    """
}
