/*
 * Cryptic variant detection steps
 */

 process FREYJA_COVARIANTS {
    publishDir "${params.output_dir}/covariants", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path trimmed_bam
    path ref
    path gff_file

    output:
    path "${trimmed_bam.baseName}.covariants.tsv"

    script:
    """
    samtools sort -o ${trimmed_bam.baseName}.sorted.bam ${trimmed_bam}
    samtools index ${trimmed_bam.baseName}.sorted.bam
    freyja covariants ${trimmed_bam.baseName}.sorted.bam ${params.min_site} ${params.max_site} --min_count ${params.min_WW_count} --gff-file ${gff_file} --ref-genome  ${ref} --output ${trimmed_bam.baseName}.covariants.tsv 
    """
}

process DETECT_CRYPTIC {

    publishDir "${params.output_dir}/cryptic", mode: 'copy'

    input:
    path covariants
    path detect_cryptic_script

    output:
    path "${covariants.baseName}.cryptic.tsv"

    script:
    """
    python ${detect_cryptic_script} --max_clinical_count ${params.max_clinical_count} --location_id ${params.location_id} --output ${covariants.baseName}.cryptic.tsv ${covariants}
    """
}