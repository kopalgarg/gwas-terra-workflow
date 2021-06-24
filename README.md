# gwas-terra-workflow
Genome-Wide Association Study (GWAS) 

Required files:
  * PLINK binary files (.bim, .bed, .fam)
  * list of samples to use (samples.tsv)
  * covariates file (columns: IID, trait, covariates)

1. `recode_plink_split_vcf.wdl`

input: 
- basename of PLINK files \
- PLINK binary files \

output: per chromosome VCF and index files \

Download VCF (.vcf.gz) and Index (.vcf.gz.csi) files from samples table

2. Upload VCFs to TOPMed and submit for QC, phasing and imputation.

`download_unzip_topmed.wdl`

input: link and password for download all files from TOPMed \

output: per chromosome VCF and index files

3. `SAIGE_association_testing.wdl`

inputs: 
  - VCF , Index files and chromosome number (from samples table)
  - trait type: "quantitative" or "binary"
  - covariates file
  - genotype: "GT" or "DS"
  - ID column name
  - phenotype column name
  - covariate column name(s)
  - file with list of samples
  - PLINK binary files \
  
output: summary statistics from the association tests saved in the samples table
