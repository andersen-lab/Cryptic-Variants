/*
 *  Process definitions
 */

process SORT { 
    input:
    path input_bam

    output:
    path "${input_bam.baseName}.sorted.bam"

    script:
    """
    samtools sort -o ${input_bam.baseName}.sorted.bam ${input_bam}
    """
}

process TRIM {
    input:
    path sorted_bam
    path primer_bed

    output:
    path "${sorted_bam.baseName}.trimmed.bam"

    script:
    """
    ivar trim -x 5 -e -i ${sorted_bam} -b ${primer_bed} -p ${sorted_bam.baseName}.trimmed.bam
    """
}

process COVARIANTS {
    publishDir "data/output", mode: 'copy'
    input:
    path trimmed_bam
    path ref
    path gff_file

    output:
    path "${trimmed_bam.baseName}.covariants.tsv"

    script:
    """
    samtools index ${trimmed_bam}
    freyja covariants ${trimmed_bam} ${params.min_site} ${params.max_site} --gff-file ${gff_file} --ref-genome  ${ref} --output ${trimmed_bam.baseName}.covariants.tsv 
    """
}

process DETECT_CRYPTIC {
    publishDir "data/output", mode: 'copy'

    input:
    path covariants
    path detect_cryptic_script

    output:
    path "${covariants.baseName}.cryptic.tsv"

    script:
    """
    python ${detect_cryptic_script} --max_clinical_count ${params.max_gisaid_count} --location_id ${params.location_id} --output ${covariants.baseName}.cryptic.tsv ${covariants}
    """
}