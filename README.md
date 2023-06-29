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
nextflow run cryptic-variants.nf \
-entry [from_fastq|from_bam] \
--input_dir <path/to/input/dir> \
--output_dir <path/to/output/dir> \
--ref <path/to/reference.fasta> \
--gff_file <path/to/gff_file.gff> \
--primer_bed <path/to/primer.bed> \
--min_site
--max_site
--min_WW_count 
--max_clinical_count
--location_id
```
The main workflow will accept either paired fastq files, or aligned bam files from `input_dir` (specified by setting `-entry` to either `from_fastq` or `from_bam`). Note that the parameters `--ref`, `--gff_file`, and `--primer_bed` are optional. If not provided, the pipeline will use the default SARS-CoV-2 reference, gff file, and primer bed file located in the `data` directory. `--input` and `--output` will default to the respective files in the `data` directory if not provided, and `--min_site` and `--max_site` default to the SARS-CoV-2 RBD.

### Example output
```
