#!/usr/bin/env nextflow

/*
 *  Wastewater Cryptic Variant Detection
 */


nextflow.enable.dsl = 2

// Define default input/output directories
params.input_dir = "$PWD/data/input"
params.output_dir = "$PWD/data/output"

// If input bam file is trimmed, set to true
params.skip_trimming = false

// Sars-Cov-2 specific parameters
params.ref = "$PWD/data/NC_045512_Hu-1.fasta"
params.gff_file = "$PWD/data/NC_045512_Hu-1.gff"
params.primer_bed = "$PWD/data/nCov-2019_v3.primer.bed"

// Freyja covariants parameters (SARS-CoV-2 RBD by default)
params.min_site = 22556
params.max_site = 23156

// Cryptic variant detection parameters
params.min_WW_count = 30
params.max_clinical_count = 5
params.location_id = 'global'

log.info """\
    ====================================================================
    W A S T E W A T E R  C R Y P T I C  V A R I A N T  D E T E C T I O N
    ====================================================================
    input dir          : ${params.input_dir}
    output dir         : ${params.output_dir}
    reference          : ${params.ref}
    gff file           : ${params.gff_file}
    primer bed file    : ${params.primer_bed}
    skip trimming      : ${params.skip_trimming}
    min genomic site   : ${params.min_site}
    max genomic site   : ${params.max_site}
    min WW count       : ${params.min_WW_count}
    max clinical count : ${params.max_clinical_count}
    location id        : ${params.location_id}
    """
    .stripIndent()

ref = file(params.ref)
gff_file = file(params.gff_file)
primer_bed = file(params.primer_bed)
detect_cryptic_script = file("$PWD/scripts/detect_cryptic.py")

// Import modules
include {
    FASTERQ_DUMP;
    GET_AMPLICON_SCHEME;
    GET_ACCESSIONS;
} from "./modules/sra.nf"


include {
    MINIMAP2;
    IVAR_TRIM;
} from "./modules/preprocessing.nf"

include {
    FREYJA_COVARIANTS;
    DETECT_CRYPTIC
} from "./modules/variant_detection.nf"

workflow from_fastq {
    Channel
        .fromFilePairs(params.input_dir + "/*{1,2}.fastq*", checkIfExists: true, size:2)
        .set { input_fastq_ch }
    MINIMAP2(input_fastq_ch, ref)
    input_bam_ch = MINIMAP2.out

    detect_cryptic(input_bam_ch)
}

workflow from_bam {
    Channel
        .fromPath(params.input_dir + "/*.bam", checkIfExists: true)
        .set { input_bam_ch }

    detect_cryptic(input_bam_ch)
}

workflow from_sra {
    Channel
        .fromPath(params.input_dir)
        .set { input_ch }

    GET_ACCESSIONS(input_ch)
        .splitCsv()
        .map { line -> line.join('') }
        .set { acc_ch }

    acc_ch.view()
    FASTERQ_DUMP(acc_ch)
    MINIMAP2(FASTERQ_DUMP.out, ref)
    input_bam_ch = MINIMAP2.out

    detect_cryptic(input_bam_ch)
}

workflow detect_cryptic {
    take:
    input_bam_ch

    main:
    if (!params.skip_trimming) {
        IVAR_TRIM(input_bam_ch, primer_bed)
        FREYJA_COVARIANTS(IVAR_TRIM.out, ref, gff_file)
    } else {
        FREYJA_COVARIANTS(input_bam_ch, ref, gff_file)
    }
    DETECT_CRYPTIC(FREYJA_COVARIANTS.out, detect_cryptic_script)
}