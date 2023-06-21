# Cryptic Variant Detection

SARS-CoV-2 wastewater cryptic lineage detection pipeline.

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
To access from GISAID, you must have a GISAID account. In order to authenticate your credentials, run the following command:
```
python -c "from outbreak_data import authenticate_user; authenticate_user.authenticate_new_user()"
```
This will provide a link prompting you to enter your GISAID username and password. Once completed, the Outbreak.info API will create and store your API key.

### Running the pipeline
```
nextflow run cryptic-variants.nf \
--input <path/to/input/dir/*.bam> \
--output <path/to/output/dir> \
[--ref <path/to/reference.fasta> --gff_file <path/to/gff_file.gff> --primer_bed <path/to/primer.bed>]
```
Note that the parameters `--ref`, `--gff_file`, and `--primer_bed` are optional. If not provided, the pipeline will use the default SARS-CoV-2 reference, gff file, and primer bed file located in the `data` directory. `--input` and `--output` will default to the respective files in the `data` directory if not provided.