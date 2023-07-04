### ................................................................................................................................................ ###

### fit HIBAG models for MICA, MICB, HLA-E, HLA-F and HLA-G in different reference data combinations (i - vii)

### ................................................................................................................................................ ###

# libraries and functions

library(tidyverse)
library(HIBAG)
library(parallel)
library(readxl)
library(openxlsx)

source('./src/FINALS/functions.R')

### ................................................................................................................................................ ###

# Fit models and test imputation accuracies using different reference data combinations:

# (I) Finnish reference - no intersect (all SNPs)
# (II) Finnish reference - intersect SNPS
# (III) 1000G EUR reference - intersect SNPs
# (IV) 1000G EUR + Finnish reference - intersect SNPs
# (V) 1000G all populations reference - intersect SNPs
# (VI) 1000G all populations + Finnish reference - intersect SNPs
# (VII) 1000G all populations - no intersect (all SNPs) 

### ................................................................................................................................................ ###

### LOAD GENOTYPE AND PHENOTYPE DATA 

### ................................................................................................................................................ ###

## Finnish reference 

# Load Plink genotype data

geno_FIN_I <- hlaBED2Geno(bed.fn = "data/Genotype_data/MHC_selected_samples_MICAB.bed", 
                       fam.fn = "data/Genotype_data/MHC_selected_samples_MICAB.fam", 
                       bim.fn = "data/Genotype_data/MHC_selected_samples_MICAB.bim", assembly="hg38") # used for model I (HLA-E, -F, MICA and MICB)
geno_FIN_II <- hlaBED2Geno(bed.fn = "data/Genotype_data/MHCmergedunion_donors_hg38.bed", 
                         fam.fn = "data/Genotype_data/MHCmergedunion_donors_hg38.fam", 
                         bim.fn = "data/Genotype_data/MHCmergedunion_donors_hg38.bim", assembly="hg38") # used for model I (HLA-F)

# remove duplicated sample ID from FIN II genotype data set
geno_FIN_II$sample.id <- str_split_fixed(geno_FIN_II$sample.id, '-', 2)[, 1]

## 1000G reference

# Load Plink genotype data and convert to hg38

Genotypedata_1000G <- hlaBED2Geno(bed.fn = "data/Tienari/plink_SNP_data/chr6_phase3_MHC_rem_dup.bed", 
                                  fam.fn = "data/Tienari/plink_SNP_data/chr6_phase3_MHC_rem_dup.fam", 
                                  bim.fn = "data/Tienari/plink_SNP_data/chr6_phase3_MHC_rem_dup.bim", assembly="hg19") # used for model VII

# hg37 -> hg38 conversion done using https://www.ncbi.nlm.nih.gov/genome/tools/remap
hg38_converted_1000G_rem_dup_snps <- read.table("./data/Tienari/report_1000Gsnps_rem_dup_hg37b.txt", sep = "\t", header = TRUE)
Genotypedata_1000G$snp.position <- hg38_converted_1000G_rem_dup_snps$mapped_start
Genotypedata_1000G$assembly <- 'hg38'

# Combine Finnish I or II with 1000G genotype data by SNP position
geno_comb_1000G_I <-hlaGenoCombine(Genotypedata_1000G, geno_FIN_I, match.type="Position") # used for models IV and VI
geno_comb_1000G_II <-hlaGenoCombine(Genotypedata_1000G, geno_FIN_II, match.type="Position") # used for models IV and VI

# extract shared SNPs from Finnish reference data and 1000G data
intersect_FIN_I_comb <- intersect(geno_FIN_I$snp.position, geno_comb_1000G_FIN_I$snp.position)
geno_FIN_I_shared <- hlaGenoSubset(geno_FIN_I, snp.sel= geno_FIN_I$snp.position %in% intersect_FIN_I_comb) # used for model II (HLA-E, -F, MICA, MICB)

intersect_FIN_II_comb <- intersect(geno_FIN_II$snp.position, geno_comb_1000G_FIN_II$snp.position)
geno_FIN_II_shared <- hlaGenoSubset(geno_FIN_II, snp.sel= geno_FIN_II$snp.position %in% intersect_FIN_II_comb) # used for model II (HLA-F)

intersect_1000G_comb_FIN_I <- intersect(Genotypedata_1000G$snp.position, geno_comb_1000G_FIN_I$snp.position)
geno_1000G_shared <- hlaGenoSubset(Genotypedata_1000G, snp.sel= Genotypedata_1000G$snp.position %in% intersect_1000G_comb_FIN_I) # used for models III and V (HLA-E, -F, MICA, MICB)

intersect_1000G_comb_FIN_II <- intersect(Genotypedata_1000G$snp.position, geno_comb_1000G_HSCT$snp.position)
geno_1000G_shared_FIN_II <- hlaGenoSubset(Genotypedata_1000G, snp.sel= Genotypedata_1000G$snp.position %in% intersect_1000G_comb_FIN_II) # used for models III and V (HLA-F)


## PHENOTYPE data 1000G

# read in 1000G phenotype data (all populations)
pheno_1000G_MICA <- read_excel("./data/Tienari/1000G_results_all_processed_pop_info.xlsx", sheet = "MICA") %>% as.data.frame()
pheno_1000G_MICB <- read_excel("./data/Tienari/1000G_results_all_processed_pop_info.xlsx", sheet = "MICB") %>% as.data.frame()
pheno_1000G_E <- read_excel("./data/Tienari/1000G_results_all_processed_pop_info.xlsx", sheet = "HLA-E") %>% as.data.frame()
pheno_1000G_F <- read_excel("./data/Tienari/1000G_results_all_processed_pop_info.xlsx", sheet = "HLA-F") %>% as.data.frame()

# get samples that are in genotype data
pheno_1000G_MICA_intersect <- filter(pheno_1000G_MICA, Sample_ID %in% Genotypedata_1000G$sample.id)
pheno_1000G_MICB_intersect <- filter(pheno_1000G_MICB, Sample_ID %in% Genotypedata_1000G$sample.id)
pheno_1000G_E_intersect <- filter(pheno_1000G_E, Sample_ID %in% Genotypedata_1000G$sample.id)
pheno_1000G_F_intersect <- filter(pheno_1000G_F, Sample_ID %in% Genotypedata_1000G$sample.id)

# Convert phenotype data into hlaAlleleClass object
mica_1000G <- hlaAllele(pheno_1000G_MICA_intersect$Sample_ID, H1=pheno_1000G_MICA_intersect[, "allele1"], H2=pheno_1000G_MICA_intersect[, "allele2"], 
                        locus = "MICA", assembly = "hg38")
micb_1000G <- hlaAllele(pheno_1000G_MICB_intersect$Sample_ID, H1=pheno_1000G_MICB_intersect[, "allele1"], H2=pheno_1000G_MICB_intersect[, "allele2"], 
                        locus = "MICB", assembly = "hg38")
hlae_1000G <- hlaAllele(pheno_1000G_E_intersect$Sample_ID, H1=pheno_1000G_E_intersect[, "allele1"], H2=pheno_1000G_E_intersect[, "allele2"], 
                        locus = "E", assembly = "hg38")
hlaf_1000G <- hlaAllele(pheno_1000G_F_intersect$Sample_ID, H1=pheno_1000G_F_intersect[, "allele1"], H2=pheno_1000G_F_intersect[, "allele2"], 
                        locus = "F", assembly = "hg38")

# Split samples randomly into 2/3 training and 1/3 test
set.seed(100)
mica_1000G_divide <- hlaSplitAllele(mica_1000G, train.prop=0.666)
set.seed(100)
micb_1000G_divide <- hlaSplitAllele(micb_1000G, train.prop=0.666)
set.seed(100)
hlae_1000G_divide <- hlaSplitAllele(hlae_1000G, train.prop=0.666)
set.seed(100)
hlaf_1000G_divide <- hlaSplitAllele(hlaf_1000G, train.prop=0.666)

# check that alleles are evenly divided between sets
hlaUniqueAllele(mica_1000G_divide$training)
hlaUniqueAllele(mica_1000G_divide$validation)
hlaUniqueAllele(micb_1000G_divide$training)
hlaUniqueAllele(micb_1000G_divide$validation)
hlaUniqueAllele(hlae_1000G_divide$training)
hlaUniqueAllele(hlae_1000G_divide$validation)
hlaUniqueAllele(hlaf_1000G_divide$training)
hlaUniqueAllele(hlaf_1000G_divide$validation)

# subset superpopulation validation samples 
EUR_id_mica <- filter(pheno_1000G_MICA_intersect, Superpopulation_code == "EUR")
EUR_mica <- hlaAlleleSubset(mica_1000G_divide$validation, samp.sel=mica_1000G_divide$validation$value$sample.id %in% EUR_id_mica$Sample_ID)
#EUR_test_geno_mica <- hlaGenoSubset(Genotypedata_1000G, samp.sel = match(EUR_mica$value$sample.id, Genotypedata_1000G$sample.id))

EUR_id_micb <- filter(pheno_1000G_MICB_intersect, Superpopulation_code == "EUR")
EUR_micb <- hlaAlleleSubset(micb_1000G_divide$validation, samp.sel=micb_1000G_divide$validation$value$sample.id %in% EUR_id_micb$Sample_ID)
#EUR_test_geno_micb <- hlaGenoSubset(Genotypedata_1000G, samp.sel = match(EUR_micb$value$sample.id, Genotypedata_1000G$sample.id))

EUR_id_e <- filter(pheno_1000G_E_intersect, Superpopulation_code == "EUR")
EUR_e <- hlaAlleleSubset(hlae_1000G_divide$validation, samp.sel=hlae_1000G_divide$validation$value$sample.id %in% EUR_id_e$Sample_ID)
#EUR_test_geno_e <- hlaGenoSubset(Genotypedata_1000G, samp.sel = match(EUR_e$value$sample.id, Genotypedata_1000G$sample.id))

EUR_id_f <- filter(pheno_1000G_F_intersect, Superpopulation_code == "EUR")
EUR_f <- hlaAlleleSubset(hlaf_1000G_divide$validation, samp.sel=hlaf_1000G_divide$validation$value$sample.id %in% EUR_id_f$Sample_ID)
#EUR_test_geno_f <- hlaGenoSubset(Genotypedata_1000G, samp.sel = match(EUR_f$value$sample.id, Genotypedata_1000G$sample.id))

# repeat same for all superpopulations (EUR, AFR, EAS, SAS, AMR)


## PHENOTYPE data Finnish reference (FIN_I and FIN_II)

MICAB_pheno_FIN_I <- read.table("./data/Phenotype_data/MICAB_alleles_merged_Histog") 
HLA_E_pheno_FIN_I <- read.table("./data/Phenotype_data/HLA_E_alleles_Histog")
HLA_F_pheno_FIN_II <- read.table("./data/Phenotype_data/HLA_F_alleles_HSCT") 
HLA_G_pheno_FIN_I <- read.table("./data/Phenotype_data/HLA_G_alleles_Histog")
HLA_G_pheno_UTR_FIN_I <- read.table("./data/Phenotype_data/HLA_G_ext_haplotypes_UTR")

mica_FIN_I <- hlaAllele(MICAB_pheno_FIN_I$SampleID, H1=MICAB_pheno_FIN_I[, "MICA_Allele1"], H2=MICAB_pheno_FIN_I[, "MICA_Allele2"], 
                         locus = "MICA", assembly = "hg38")
micb_FIN_I <- hlaAllele(MICAB_pheno_FIN_I$SampleID, H1=MICAB_pheno_FIN_I[, "MICB_Allele1"], H2=MICAB_pheno_FIN_I[, "MICB_Allele2"],
                         locus = "MICB", assembly = "hg38")
HLA_E_FIN_I <- hlaAllele(HLA_E_pheno_FIN_I$SampleID, H1=HLA_E_pheno_FIN_I[, "HLA_E_Allele1"], H2=HLA_E_pheno_FIN_I[, "HLA_E_Allele2"], 
                          locus = "E", assembly = "hg38")
HLA_F_FIN_II <- hlaAllele(HLA_F_pheno_FIN_II$SampleID, H1=HLA_F_pheno_FIN_II[, "HLA_F_Allele1"], H2=HLA_F_pheno_FIN_II[, "HLA_F_Allele2"], 
                        locus = "F", assembly = "hg38")
HLA_G_FIN_I <- hlaAllele(HLA_G_pheno_FIN_I$SampleID, H1=HLA_G_pheno_FIN_I[, "HLA_G_Allele1"], H2=HLA_G_pheno_FIN_I[, "HLA_G_Allele2"], 
                          locus = "G", assembly = "hg38")
HLA_G_5UTR_FIN_I <- hlaAllele(HLA_G_pheno_UTR_FIN_I$SampleID, H1=HLA_G_pheno_UTR_FIN_I[, "UTR5_Allele1"], H2=HLA_G_pheno_UTR_FIN_I[, "UTR5_Allele2"],
                        locus = "G", assembly = "hg38")
HLA_G_3UTR_FIN_I <- hlaAllele(HLA_G_pheno_UTR_FIN_I$SampleID, H1=HLA_G_pheno_UTR_FIN_I[, "UTR3_Allele1"], H2=HLA_G_pheno_UTR_FIN_I[, "UTR3_Allele2"], 
                        locus = "G", assembly = "hg38")

# Extract one of the samples with alleles present in <3 samples  (MICA*011:01, MICA*009:02, MICA*001, MICB*005:08, MICB*005:06, MICB*005:01, HLA-G*01:01:01:09, HLA-G*01:21N)
mica_FIN_I_rare <- hlaAlleleSubset(mica_FIN_I, c(241,277,412))
micb_FIN_I_rare <- hlaAlleleSubset(micb_FIN_I, c(82,153,216))
hla_G_FIN_I_rare <- hlaAlleleSubset(HLA_G_FIN_I, c(117,257))

# subset rest of the samples
mica_FIN_I_rare_excluded <- hlaAlleleSubset(mica_FIN_I, samp.sel = !(mica_FIN_I$value$sample.id %in% mica_FIN_I_rare$value$sample.id))
micb_FIN_I_rare_excluded <- hlaAlleleSubset(micb_FIN_I, samp.sel = !(micb_FIN_I$value$sample.id %in% micb_FIN_I_rare$value$sample.id))
hla_G_FIN_I_rare_excluded <- hlaAlleleSubset(HLA_G_FIN_I, samp.sel = !(HLA_G_FIN_I$value$sample.id %in% hla_G_FIN_I_rare$value$sample.id))

# Split samples randomly into 2/3 training and 1/3 test
set.seed(100)
micatab <- hlaSplitAllele(mica_FIN_I_rare_excluded, train.prop=0.666)
set.seed(100)
micbtab <- hlaSplitAllele(micb_FIN_I_rare_excluded, train.prop=0.666)
set.seed(100)
hlaetab <- hlaSplitAllele(HLA_E_FIN_I, train.prop=0.666)
set.seed(100)
hlaftab <- hlaSplitAllele(HLA_F_FIN_II, train.prop=0.666)
set.seed(100)
hlagtab <- hlaSplitAllele(hla_G_FIN_I_rare_excluded, train.prop=0.666)
set.seed(100)
hlag5UTRtab <- hlaSplitAllele(HLA_G_5UTR_FIN_I, train.prop=0.666)
set.seed(100)
hlag3UTRtab <- hlaSplitAllele(HLA_G_3UTR_FIN_I, train.prop=0.666)

# add samples of rare alleles to validation sets
micatab$validation<- hlaCombineAllele(micatab$validation, mica_FIN_I_rare)
micbtab$validation<- hlaCombineAllele(micbtab$validation, micb_FIN_I_rare)
hlagtab$validation<- hlaCombineAllele(hlagtab$validation, hla_G_FIN_I_rare)

# check that rare alleles are evenly divided between sets
hlaUniqueAllele(micatab$training)
hlaUniqueAllele(micatab$validation)
hlaUniqueAllele(micbtab$training)
hlaUniqueAllele(micbtab$validation)
hlaUniqueAllele(hlaetab$training)
hlaUniqueAllele(hlaetab$validation)
hlaUniqueAllele(hlaftab$training)
hlaUniqueAllele(hlaftab$validation)
hlaUniqueAllele(hlagtab$training)
hlaUniqueAllele(hlagtab$validation)
hlaUniqueAllele(hlag5UTRtab$training)
hlaUniqueAllele(hlag5UTRtab$validation)
hlaUniqueAllele(hlag3UTRtab$training)
hlaUniqueAllele(hlag3UTRtab$validation)

### ................................................................................................................................................ ###

# phenotype sets for training

# (I) & (II): sets created previously in this script: 
# micatab$training, micb = micbtab$training, hlae = hlaetab$training, hlaf= hlaftab$training, hlag= hlagtab$training, hlag5UTRtab$training, hlag3UTRtab$training)

# (III): 1000G EUR:
EUR_iii_mica <- hlaAlleleSubset(mica_1000G_divide$training, samp.sel=mica_1000G_divide$training$value$sample.id %in% EUR_id_mica$Sample_ID)
EUR_iii_micb <- hlaAlleleSubset(micb_1000G_divide$training, samp.sel=micb_1000G_divide$training$value$sample.id %in% EUR_id_micb$Sample_ID)
EUR_iii_e <- hlaAlleleSubset(hlae_1000G_divide$training, samp.sel=hlae_1000G_divide$training$value$sample.id %in% EUR_id_e$Sample_ID)
EUR_iii_f <- hlaAlleleSubset(hlaf_1000G_divide$training, samp.sel=hlaf_1000G_divide$training$value$sample.id %in% EUR_id_f$Sample_ID)

# (IV) 1000G EUR + Finnish reference:
mica_iv <- hlaCombineAllele(EUR_iii_mica, micatab$training)
micb_iv <- hlaCombineAllele(EUR_iii_micb, micbtab$training)
hlae_iv <- hlaCombineAllele(EUR_iii_e, hlaetab$training)
hlaf_iv <- hlaCombineAllele(EUR_iii_f, hlaftab$training)

# (V) and (VII): sets created previously in this script: mica_1000G_divide$training, micb_1000G_divide$training, hlae_1000G_divide$training, hlaf_1000G_divide$training

# (VI) 1000G all populations + Finnish reference:
mica_vi <- hlaCombineAllele(mica_1000G_divide$training, micatab$training)
micb_vi <- hlaCombineAllele(micb_1000G_divide$training, micbtab$training)
hlae_vi <- hlaCombineAllele(hlae_1000G_divide$training, hlaetab$training)
hlaf_vi <- hlaCombineAllele(hlaf_1000G_divide$training, hlaftab$training)

### ................................................................................................................................................ ###

###           FLANKING REGION TESTINGS

### ................................................................................................................................................ ###

# mica
# loop to fit models with different flanking regions in FIN and 1000G data
loop_model_FIN_I_mica <- map(c(seq(1000, 15000, by=2000), 50000), function(x) fitHLAmodel('MICA', micatab$training, geno_BB, x, 50))
loop_model_1000G_mica <- map(c(seq(1000, 15000, by=2000), 50000), function(x) fitHLAmodel('MICA', mica_1000G_divide$training, Genotypedata_1000G, x, 50))

loop_model_FIN_I_mica_out <- map_dfr(1:length(loop_model_FIN_I_mica), function(x) summary(loop_model_FIN_I_mica[[x]]) %>% .$info %>% .[3, ]) # out of bag accuracies
loop_model_1000G_mica_out <- map_dfr(1:length(loop_model_1000G_mica), function(x) summary(loop_model_1000G_mica[[x]]) %>% .$info %>% .[3, ]) # out of bag accuracies

# pred results  
loop_results_model_FIN_I_mica <- map(1:length(loop_model_FIN_I_mica), function (x) {
  PredHlaModel(loop_model_FIN_I_mica[[x]], geno_FIN_I, micatab$validation)
})

loop_results_model_1000G_mica <- map(1:length(loop_model_1000G_mica), function (x) {
  PredHlaModel(loop_model_1000G_mica[[x]], Genotypedata_1000G, mica_1000G_divide$validation)
})

# get overall accuracies from the loop´
results_FIN_I_loop_mica <- accuracy_loop(loop_results_model_FIN_I_mica, 1)
results_1000G_loop_mica <- accuracy_loop(loop_results_model_1000G_mica, 1)

# micb
# loop to fit models with different flanking regions in FIN and 1000G data
loop_model_FIN_I_micb <- map(c(seq(1000, 15000, by=2000), 50000), function(x) fitHLAmodel('MICB', micbtab$training, geno_FIN_I, x, 50))
loop_model_1000G_micb <- map(c(seq(1000, 15000, by=2000), 50000), function(x) fitHLAmodel('MICB', micb_1000G_divide$training, Genotypedata_1000G, x, 50))

loop_model_FIN_I_micb_out <- map_dfr(1:length(loop_model_FIN_I_micb), function(x) summary(loop_model_FIN_I_micb[[x]]) %>% .$info %>% .[3, ]) # out of bag accuracies
loop_model_1000G_micb_out <- map_dfr(1:length(loop_model_1000G_micb), function(x) summary(loop_model_1000G_micb[[x]]) %>% .$info %>% .[3, ]) # out of bag accuracies

# pred results 
loop_results_model_FIN_I_micb <- map(1:length(loop_model_FIN_I_micb), function (x) {
  PredHlaModel(loop_model_FIN_I_micb[[x]], geno_FIN_I, micbtab$validation)
})

loop_results_model_1000G_micb <- map(1:length(loop_model_1000G_micb), function (x) {
  PredHlaModel(loop_model_1000G_micb[[x]], Genotypedata_1000G, micb_1000G_divide$validation)
})

# get overall accuracies from the loop
results_FIN_I_loop_micb <- accuracy_loop(loop_results_model_FIN_I_micb, 1)
results_1000G_loop_micb <- accuracy_loop(loop_results_model_1000G_micb, 1)

# HLA-E
# loop to fit models with different flanking regions in FIN and 1000G data
loop_model_FIN_I_hlae <- map(c(seq(1000, 15000, by=2000), 50000), function(x) fitHLAmodel('E', hlaetab$training, geno_FIN_I, x, 50))
loop_model_1000G_hlae <- map(c(seq(1000, 15000, by=2000), 50000), function(x) fitHLAmodel('E', hlae_1000G_divide$training, Genotypedata_1000G, x, 50))

loop_model_FIN_I_hlae_out <- map_dfr(1:length(loop_model_FIN_I_hlae), function(x) summary(loop_model_FIN_I_hlae[[x]]) %>% .$info %>% .[3, ]) # out of bag accuracies
loop_model_1000G_hlae_out <- map_dfr(1:length(loop_model_1000G_hlae), function(x) summary(loop_model_1000G_hlae[[x]]) %>% .$info %>% .[3, ]) # out of bag accuracies

# pred results
loop_results_model_FIN_I_hlae <- map(1:length(loop_model_FIN_I_hlae), function (x) {
  PredHlaModel(loop_model_FIN_I_hlae[[x]], geno_FIN_I, hlaetab$validation)
})

loop_results_model_1000G_hlae <- map(1:length(loop_model_1000G_hlae), function (x) {
  PredHlaModel(loop_model_1000G_hlae[[x]], Genotypedata_1000G, hlae_1000G_divide$validation)
})

# get overall accuracies from the loop
results_FIN_I_loop_hlae <- accuracy_loop(loop_results_model_FIN_I_hlae, 1)
results_1000G_loop_hlae <- accuracy_loop(loop_results_model_1000G_hlae, 1)

# HLA-F
# loop to fit models with different flanking regions in FIN and 1000G data
loop_model_FIN_II_hlaf <- map(c(seq(1000, 15000, by=2000), 50000), function(x) fitHLAmodel('F', hlaftab$training, geno_FIN_II, x, 50))
loop_model_1000G_hlaf <- map(c(seq(1000, 15000, by=2000), 50000), function(x) fitHLAmodel('F', hlaf_1000G_divide$training, Genotypedata_1000G, x, 50))

loop_model_FIN_II_hlaf_out <- map_dfr(1:length(loop_model_FIN_II_hlaf), function(x) summary(loop_model_FIN_II_hlaf[[x]]) %>% .$info %>% .[3, ]) # out of bag accuracies
loop_model_1000G_hlaf_out <- map_dfr(1:length(loop_model_1000G_hlaf), function(x) summary(loop_model_1000G_hlaf[[x]]) %>% .$info %>% .[3, ]) # out of bag accuracies

# pred results
loop_results_model_FIN_II_hlaf <- map(1:length(loop_model_FIN_II_hlaf), function (x) {
  PredHlaModel(loop_model_FIN_II_hlaf[[x]], geno_FIN_II, hlaftab$validation)
})

loop_results_model_1000G_hlaf <- map(1:length(loop_model_1000G_hlaf), function (x) {
  PredHlaModel(loop_model_1000G_hlaf[[x]], Genotypedata_1000G, hlaf_1000G_divide$validation)
})

# get overall accuracies from the loop
results_FIN_II_loop_hlaf <- accuracy_loop(loop_results_model_FIN_II_hlaf, 1)
results_1000G_loop_hlaf <- accuracy_loop(loop_results_model_1000G_hlaf, 1)

# HLA-G
# loop to fit models with different flanking regions in FIN data
loop_model_FIN_I_hlag <- map(c(seq(1000, 15000, by=2000), 50000), function(x) fitHLAmodel('G', hlagtab$training, geno_FIN_I, x, 50))
loop_model_FIN_I_hlag_out <- map_dfr(1:length(loop_model_FIN_I_hlag), function(x) summary(loop_model_FIN_I_hlag[[x]]) %>% .$info %>% .[3, ]) # out of bag accuracies

# pred results
loop_results_model_FIN_I_hlag <- map(1:length(loop_model_FIN_I_hlag), function (x) {
  PredHlaModel(loop_model_FIN_I_hlag[[x]], geno_FIN_I, hlagtab$validation)
})

# get overall accuracies from the loop
results_FIN_I_loop_hlag <- accuracy_loop(loop_results_model_FIN_I_hlag, 1)

### ................................................................................................................................................ ###

###             TRAIN MODELS USING TRAINING PROPORTION (2/3) OF DATA

### ................................................................................................................................................ ###

# train models with different data combinations, using training sets, 100 classifiers and 10 kb flanking region
# (I) FIN (no intersect)
fitHLAmodel("MICA", micatab$training, geno_FIN_I, 10000, 100) # MICA,MICB,HLA-E, HLA-G: geno_FIN_I, HLA-F: geno_FIN_II
# (II) FIN
fitHLAmodel("MICA", micatab$training, geno_FIN_I_shared, 10000, 100) # MICA,MICB,HLA-E: geno_FIN_I, HLA-F: geno_FIN_II_shared
# (III) EUR
fitHLAmodel("MICA", EUR_iii_mica, geno_1000G_shared, 10000, 100) # MICA, MICB, HLA-E: geno_1000G_shared, HLA-F: geno_1000G_shared_FIN_II
# (IV) EUR + FIN
fitHLAmodel("MICA", mica_iv, geno_comb_1000G_FIN_I, 10000, 100) # MICA, MICB, HLA-E: geno_comb_1000G_FIN_I, HLA-F: geno_comb_1000G_FIN_II
# (V) 1000G all
fitHLAmodel("MICA", mica_1000G_divide$training, geno_1000G_shared, 10000, 100) # MICA, MICB, HLA-E: geno_1000G_shared, HLA-F: geno_1000G_shared_FIN_II
# (VI) 1000G all + FIN
fitHLAmodel("MICA", mica_vi, geno_comb_1000G_FIN_I, 10000, 100) # MICA, MICB, HLA-E: geno_comb_1000G_FIN_I, HLA-F: geno_comb_1000G_FIN_II 
# (VII) 1000G all (no intersect)
fitHLAmodel("MICA", mica_1000G_divide$training, Genotypedata_1000G, 10000, 100) 

# HLA-G (training only with FIN I reference)
fitHLAmodel("G", hlagtab$training, geno_FIN_I, 10000,100)
fitHLAmodel("G", hlag3UTRtab$training, geno_FIN_I, 10000,100)
fitHLAmodel("G", hlag5UTRtab$training, geno_FIN_I, 10000,100)

### ................................................................................................................................................ ###

### FIT MODELS WITH ALL DATA

### ................................................................................................................................................ ###

# phenotype sets for all data

# (I) & (II): sets created previously in this script: mica_Histog, micb_Histog, HLA_E_Histog, HLA_F_HSCT

# (III): 1000G EUR:
EUR_mica_all <- hlaAlleleSubset(mica_1000G, samp.sel=mica_1000G$value$sample.id %in% EUR_id_mica$Sample_ID) 
EUR_micb_all <- hlaAlleleSubset(micb_1000G, samp.sel=micb_1000G$value$sample.id %in% EUR_id_micb$Sample_ID) 
EUR_e_all <- hlaAlleleSubset(hlae_1000G, samp.sel=hlae_1000G$value$sample.id %in% EUR_id_e$Sample_ID)
EUR_f_all <- hlaAlleleSubset(hlaf_1000G, samp.sel=hlaf_1000G$value$sample.id %in% EUR_id_f$Sample_ID)

# (IV) 1000G EUR + Finnish reference:
mica_iv_all <- hlaCombineAllele(EUR_mica_all, mica_FIN_I) 
micb_iv_all <- hlaCombineAllele(EUR_micb_all, micb_FIN_I) 
e_iv_all <- hlaCombineAllele(EUR_e_all, HLA_E_FIN_I)
f_iv_all <- hlaCombineAllele(EUR_f_all, HLA_F_FIN_II)

# (V) & (VII: sets created previously in this script: mica_1000G, micb_1000G, hlae_1000G, hlaf_1000G

# (VI) 1000G all populations + Finnish reference:
mica_vi_all <- hlaCombineAllele(mica_1000G, mica_FIN_I) 
micb_vi_all <- hlaCombineAllele(micb_1000G, micb_FIN_I) 
e_vi_all <- hlaCombineAllele(hlae_1000G, HLA_E_FIN_I)
f_vi_all <- hlaCombineAllele(hlaf_1000G, HLA_F_FIN_II)

# train models with different data combinations, using whole data sets, 100 classifiers and 10 kb flanking region
# (I) FIN (no intersect)
fitHLAmodel("MICA", mica_FIN_I, geno_FIN_I, 10000, 100) # MICA,MICB,HLA-E: geno_FIN_I, HLA-F: geno_FIN_II
# (II) FIN
fitHLAmodel("MICA", mica_FIN_I, geno_FIN_I_shared, 10000, 100) # MICA,MICB,HLA-E: geno_FIN_I, HLA-F: geno_FIN_II_shared
# (III) EUR
fitHLAmodel("MICA", EUR_mica_all, geno_1000G_shared, 10000, 100) # MICA, MICB, HLA-E: geno_1000G_shared, HLA-F: geno_1000G_shared_FIN_II
# (IV) EUR + FIN
fitHLAmodel("MICA", mica_iv_all, geno_comb_1000G_FIN_I, 10000, 100) # MICA, MICB, HLA-E: geno_comb_1000G_FIN_I, HLA-F: geno_comb_1000G_FIN_II
# (V) 1000G all
fitHLAmodel("MICA", mica_1000G, geno_1000G_shared, 10000, 100) # MICA, MICB, HLA-E: geno_1000G_shared, HLA-F: geno_1000G_shared_FIN_II
# (VI) 1000G all + FIN
fitHLAmodel("MICA", mica_vi_all, geno_comb_1000G_FIN_I, 10000, 100) # MICA, MICB, HLA-E: geno_comb_1000G_FIN_I, HLA-F: geno_comb_1000G_FIN_II 
# (VII) 1000G all (no intersect)
fitHLAmodel("MICA", mica_1000G, Genotypedata_1000G, 10000, 100) 

# HLA-G (training only with FIN I reference) 
fitHLAmodel("G", HLA_G_FIN_I, geno_FIN_I, 10000,100)
fitHLAmodel("G", HLA_G_3UTR_FIN_I, geno_FIN_I, 10000,100)
fitHLAmodel("G", HLA_G_5UTR_FIN_I, geno_FIN_I, 10000,100)
