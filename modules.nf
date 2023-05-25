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

process INDEX {
    input:
    path sorted_bam

    output:
    path "${sorted_bam}.bai"

    script:
    """
    samtools index ${sorted_bam}
    """
}

process TRIM {
    input:
    path sorted_bam

    output:
    path "${sorted_bam.baseName}.trimmed.bam"

    script:
    """
    ivar trim -x 5 -e -i ${sorted_bam} -b ${params.primer_bed} -p ${sorted_bam.baseName}.trimmed.bam
    """
}

process COVARIANTS {
    input:
    path trimmed_bam

    output:
    path "${trimmed_bam.baseName}.covariants.tsv"

    script:
    """
    freyja covariants ${trimmed_bam} ${params.min_site} ${params.max_site} --ref-genome ${params.ref} --output ${trimmed_bam.baseName}.covariants.tsv 
    """
}