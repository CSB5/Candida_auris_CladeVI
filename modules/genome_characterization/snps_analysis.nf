process VARIENT_ANNOTATION_SNPEFF {
    tag "Annotating varient using snpEff"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )

    input:
    path vcf_ch
    val db_name
    val suffix

    output: 
    path "${params.sample_id}.${suffix}.annot.vcf"

    script:
    """
    sed -i 's/_C_auris_B8441/.1/g' ${vcf_ch}

    java -Xmx8g -jar ${params.snpeff_path} -v ${db_name} \
    ${vcf_ch} > ${params.sample_id}.${suffix}.annot.vcf
    """

}

process VARIENT_CALLING_SNIPPY {
    tag "Annotating varient using snpEff"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )

    input:
    path ref_ch
    path illumina_reads_1_ch
    path illumina_reads_2_ch
    val suffix

    output: 
    path "${params.sample_id}.${suffix}/snps.vcf"

    script:
    """
    snippy --cpus ${task.cpus} --outdir ${params.sample_id}.${suffix} --ref ${ref_ch} \
    --R1 ${illumina_reads_1_ch} --R2 ${illumina_reads_2_ch}
    """

}
