#!/usr/bin/env nextflow

/*
 *  Wastewater Cryptic Variant Detection
 */

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

// Define input/output directories
params.input = "$baseDir/data/input/*.bam"
params.output = "$baseDir/data/output"

params.is_trimmed = true

// Sars-Cov-2 specific parameters
params.ref = "$PWD/data/NC_045512_Hu-1.fasta"
params.gff_file = "$PWD/data/NC_045512_Hu-1.gff"
params.primer_bed = "$PWD/data/nCov-2019_v3.primer.bed"

// Freyja covariants parameters (S-gene by default)
params.min_site = 21563
params.max_site = 25384

// Cryptic variant detection parameters
params.min_WW_count = 50
params.max_gisaid_count = 10
params.location_id = "USA"

ref = file(params.ref)
gff_file = file(params.gff_file)
primer_bed = file(params.primer_bed)

detect_cryptic_script = file("$PWD/scripts/detect_cryptic.py")

// Import modules
include {
    SORT;
    TRIM;
    COVARIANTS;
    DETECT_CRYPTIC;
} from "./modules.nf"

Channel 
    .fromPath(params.input)
    .set { input_bam_ch }

workflow {
    SORT(input_bam_ch)

    if (!params.is_trimmed) {
        TRIM(SORT.out)
        COVARIANTS(TRIM.out, ref, gff_file)
    } else {
        COVARIANTS(SORT.out, ref, gff_file)
    }
    DETECT_CRYPTIC(COVARIANTS.out, detect_cryptic_script)
}