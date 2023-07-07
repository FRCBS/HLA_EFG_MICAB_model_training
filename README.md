# Training_and_validation_of_HLA_EFG_MICAB_imp_models

Training and validation of imputation reference panels for HLA-E, -F and -G and MICA/B using Finnish and 1000G reference data

## code (./src)

`extract_MIC.sh` extraction of reads in gene area from 1000 Genomes whole exome sequencing data

`HLA_G_haplotype_processing.R` inference of HLA-G UTR reference haplotypes from SNP data

`Model_training.R` training models in data compositions I-VII using training and whole reference data

`Validation.R` validation of models 

`imputation.R` example of running imputation using models VI for MICA, MICB, HLA-E, HLA-F and models I for HLA-G, HLA-G 3'UTR and HLA-G 5'UTR

`plot_results.R` plotting results
