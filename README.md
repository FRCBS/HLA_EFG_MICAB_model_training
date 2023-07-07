# HLA-E, HLA-F, HLA-G, MICA and MICB imputation reference panels

Training and validation of imputation reference panels for HLA-E, -F and -G and MICA/B using Finnish and 1000G reference data

## code (./src)

`extract_MIC.sh` extraction of reads in gene area from 1000 Genomes whole exome sequencing data

`HLA_G_haplotype_processing.R` inference of HLA-G UTR reference haplotypes from SNP data

`Model_training.R` training models in data compositions I-VII using training and whole reference data

`Validation.R` validation of models 

`plot_results.R` plotting results
