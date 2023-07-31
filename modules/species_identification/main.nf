process FUNGAL_SPECIES_IDEN_MASH_READ {
    tag "Using mash to identify species from illumina paired reads"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )
    containerOptions "--volume ${params.fungal_mash_db_path}:${params.fungal_mash_db_path}"

    input:
    path read_ch1
    path read_ch2

    output:
    path "${params.sample_id}.read.k21.msh"
    path "${params.sample_id}.read.k21.mash.tsv"

    script:
    """
    mash sketch -o ${params.sample_id}.read.k21.msh -m 2 -s 10000 -g ${params.genome_size} -r ${read_ch1} ${read_ch2}
    ls  ${params.fungal_mash_db_path}/*.k21.msh > k21_mash_list.txt
    mash dist ${params.sample_id}.read.k21.msh -l k21_mash_list.txt > ${params.sample_id}.read.k21.mash.tsv
    """

}

process FUNGAL_SPECIES_IDEN_MASH_CONTIG {
    tag "Using mash to identify species from illumina contigs"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )
    containerOptions "--volume ${params.fungal_mash_db_path}:${params.fungal_mash_db_path}"

    input:
    path assembly_ch

    output:
    path "${params.sample_id}.contig.k21.msh"
    path "${params.sample_id}.contig.k21.mash.tsv"

    script:
    """
    mash sketch -o ${params.sample_id}.contig.k21.msh -s 10000 -g ${params.genome_size} ${assembly_ch}
    ls  ${params.fungal_mash_db_path}/*.k21.msh > k21_mash_list.txt
    mash dist ${params.sample_id}.contig.k21.msh -l k21_mash_list.txt > ${params.sample_id}.contig.k21.mash.tsv
    """

}
