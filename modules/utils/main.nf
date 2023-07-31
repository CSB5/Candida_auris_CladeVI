
process GET_FASTQ_FROM_SRA {

    input:
    val sra_accession

    output:
    path "${sra_accession}_1.fastq"
    path "${sra_accession}_2.fastq"

    script: 
    """
    prefetch ${sra_accession}
    fasterq-dump --split-files ${sra_accession}
    """
}