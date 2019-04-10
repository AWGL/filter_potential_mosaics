## Expanded filter for potential mosiac cases

DISCLAIMER - Neither the GermlineEnrichment pipeline or this script have been validated to detect mosaic variants.
The purpose of this script is to pull out variants previously filtered as part of the pipeline, so that the 
clinical scientists have an expanded list of variants that they can manually investigate in the case that there is a 
suspected mosaic variant. This is not a script to detect mosiac variants, it only parses VCF files that have been filtered. 
Therefore, there is a high chance that the variants are artefacts and each variant must be manually inspected in IGV and 
confirmed using another technique. 

---
### Program output

The program outputs a TSV file containing variants that were filtered out at each filtering step, in reverse order of when 
they were filtered. The variants in filter step 1 were filtered out just before uploading to the variant database, 
the variants in filter step 3 were filtered out much earlier and are therefore more likely to be artefacts.

Each variant should be checked in IGV and confirmed by another method.

The data included in the columns for each variant are:
- `chr`: The chromosome the variant is located on
- `pos`: The genomic co-ordinate of the variant
- `ref`: The reference allele
- `alt`: The alternate allele, there may be more than one, seperated by a comma
- `AD`: The number of reads supporting each allele seperated by commas, the first number is the ref allele and each 
subsequent number corresponds to the ref allele(s) in the same order as they are listed in the ref column
- `DP`: The total number of reads for this position

---
### Setup

1. Clone the repository: `git clone https://github.com/erikwaskiewicz/mosiac.git`

2. Setup the Conda environment:

  - If not already installed, install miniconda from [here](https://conda.io/en/latest/miniconda.html)  
  - Make the conda environment: `conda env create -f mosiac_env.yaml`  
  - Activate the environment: `conda activate mosaic_env.yaml`  

3. Edit the config file:

`bcftools_path`: The path to the executable file for the bcftools program  
`output_path`: The path to the folder where the output files will be saved  
`genes`: A dictionary of genes that can be applied, any additional genes added in the future should be added here.  
Each entry must be in the format:  
`'<gene_name: string>': ['<chrom: string>', <start_pos: int>,  <end_pos: int>]`

### Running the program

The script is designed to work with TruSightCancer data from the Germline Enrichment pipeline versions 2.5.0 onwards only. 
The script will also work for TruSightOne data run through the same pipeline, but additional functionality for applying 
BED files will likely have to be developed, due to the large number of genes in TruSightOne.

To get help: `python mosaic.py -h`

**Argument descriptions**

`run_folder`: path the the run folder, where GermlineEnrichment script 2 runs. e.g. 
/data/results/<run_id>/IlluminaTruSightCancer/

`sample_id`: The patient ID, e.g. 19M123456

**gene flags**: Each gene within the config file will have it's own flag, to include the gene in the output, use the flag for 
that gene. To include multiple genes, use multiple flags.  
e.g. `python mosaic.py <run_folder> <sample_id> --NF1` for NF1 only,  
or `python mosaic.py <run_folder> <sample_id> --TSC1 --TSC2` for both TSC1 and TSC2.
