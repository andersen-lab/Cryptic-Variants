# Cryptic Variant Detection

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](http://nextflow.io)

SARS-CoV-2 wastewater cryptic variant detection pipeline.

## Installation
---

### Install via Git
```
git clone https://github.com/dylanpilz/cryptic-variants.git
cd cryptic-variants
```
## Environment setup
---

#### Via conda:
```
conda env create -f environment.yml
conda activate cryptic-variants
```
#### Via mamba:
```
mamba create -n cryptic-variants
mamba env update -n cryptic-variants --file environment.yml
mamba activate cryptic-variants
```

## Usage
---
### GISAID Authentication
To access data from GISAID, you must have a GISAID account. In order to authenticate your credentials, run the following command:
```
python -c "from outbreak_data import authenticate_user; authenticate_user.authenticate_new_user()"
```
This will provide a link prompting you to enter your GISAID username and password. Once completed, your token will be saved and you will not need to authenticate again.

### Running the pipeline
```
nextflow run main.nf -entry [from_fastq|from_bam] --input_dir <path/to/input/dir>  --output_dir <path/to/output/dir> 
```
Set `-entry` to `from_fastq` if you are providing paired fastq files, or `from_bam` if you are providing aligned bam files. In the output directory, you will find a `cryptic_variants` directory containing potential cryptic variants, as well as a `covariants` directory containing the output from `freyja covariants` for the provided samples (see [freyja](https://github.com/andersen-lab/Freyja)).

### Optional parameters
```
--ref <path/to/reference.fasta>
            Reference genome to use for alignment and covariant detection
            (default: data/NC_045512.2_Hu-1.fasta)
--gff_file <path/to/gff_file.gff>
            GFF file containing gene annotations
            (default: data/NC_045512.2_Hu-1.gff)
--primer_bed <path/to/primer.bed>
            BED file containing primer locations for primer trimming
            (default: data/nCoV-2019_v3.primer.bed)
--skip_trimming <true|false>
            Whether or not to trim primer sequences from reads. If true, primer
            trimming will be skipped.
            (default: false)
--min_site <int>
            Minimum genomic site to consider for cryptic variant detection
            (default: 22556) (RBD start)
--max_site <int>
            Maximum genomic site to consider for cryptic variant detection
            (default: 23156) (RBD end)
--min_WW_count <int>
            Minimum number of wastewater hits to consider a cluster of
            variants in a given sample
            (default: 30)
--max_clinical_count <int>
            Maximum number of clinical hits for a variant to be considered
            cryptic
            (default: 5)
--location_id <str>
            Location ID to query from GISAID
            (default: 'global')
```
### Output

The pipeline produces two output directories: `cryptic_variants` and `covariants`. The `cryptic_variants` directory will contain a `{sample}.covariants.cryptic.tsv` for each sample in the input directory. This file contains the following columns:

`Covariants`
    Mutation cluster detected
`WW_Count`
    Number of wastewater hits for the mutation cluster in this sample
`Clinical_Count`
    Number of clinical hits for the mutation cluster in this sample
`Lineages`
    Lineages associated with the mutation cluster (if any)

```
Covariants      WW_Count        Clinical_Count  Lineages
['S:G416E', 'S:K417N', 'S:N440K', 'S:L452Q']    12      0       NA
['S:K417N', 'S:S438P', 'S:N440K', 'S:L452Q']    14      0       NA
['S:K417N', 'S:Y421H', 'S:N440K', 'S:L452Q']    11      0       NA
['S:G416G', 'S:K417N', 'S:N440K', 'S:L452Q']    11      2       ['ba.2.12.1']
['S:K417N', 'S:N440K', 'S:G446G', 'S:L452Q']    14      2       ['bg.5']
['S:K417N', 'S:N440K', 'S:L452Q', 'S:F456F']    15      2       ['ba.2.12.1']
```

The `covariants` directory contains the raw output from [freyja covariants](https://github.com/andersen-lab/Freyja).