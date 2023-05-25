#!/usr/bin/env nextflow

/*
 *  SARS-CoV-2 Cryptic Variant Detection
 */

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

// Define default parameters
params.input = "$baseDir/data/input_bam/*.bam"
params.is_trimmed = true
params.ref = "$baseDir/data/NC_045512_Hu-1.fasta"
params.primer_bed = "$baseDir/data/nCov-2019_v3.primer.bed"

// Freyja covariants parameters
params.min_site = 21563
params.max_site = 25384

// Cryptic variant detection parameters
params.min_WW_count = 10
params.max_gisaid_count = 10
params.location_id = "USA"

// Import modules
include {
    SORT;
    INDEX;
    TRIM;
    COVARIANTS;
} from "./modules.nf"

Channel 
    .fromPath(params.input)
    .set { input_bam_ch }

workflow {
    SORT(input_bam_ch)
    INDEX(SORT.out)

    if (!params.is_trimmed) {
        TRIM(SORT.out)
        COVARIANTS(TRIM.out)
    } else {
        COVARIANTS(SORT.out)
    }
}