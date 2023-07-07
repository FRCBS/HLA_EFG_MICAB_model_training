### ................................................................................................................................................ ###

#                              CROSS-VALIDATION

### ................................................................................................................................................ ###

# libraries and functions

library(HIBAG)
library(tidyverse)

source('./src/functions.R')

### ................................................................................................................................................ ###

#                   IMPORT MODEL FILES

# list model files (created in Model_training.R)
model_files_mica <- list.files(path="./data/Imputation_models/Combination_testing_models/MICA", full.names=F, recursive=F)
model_files_micb <- list.files(path="./data/Imputation_models/Combination_testing_models/MICB", full.names=F, recursive=F)
model_files_hlae <- list.files(path="./data/Imputation_models/Combination_testing_models/HLA_E", full.names=F, recursive=F)
model_files_hlaf <- list.files(path="./data/Imputation_models/Combination_testing_models/HLA_F", full.names=F, recursive=F)

# load models
models_mica <- load_models(model_files_mica, "MICA")
models_micb <- load_models(model_files_micb, "MICB")
models_hlae <- load_models(model_files_hlae, "HLA_E")
models_hlaf <- load_models(model_files_hlaf, "HLA_F")

# change the order of MICA, MICB, HLA-E and HLA-F models to correspond i - vii
models_mica <- models_mica[c("model_MICA_geno_FIN_I_micatab$training.RData","model_MICA_geno_comb_1000G_I_micatab$training.RData", 
                             "model_MICA_geno_comb_1000G_I_EUR_iii_mica.RData", "model_MICA_geno_comb_1000G_I_mica_iv.RData",
                             "model_MICA_geno_comb_1000G_I_mica_1000G_divide$training.RData", "model_MICA_geno_comb_1000G_I_mica_vi.RData",
                             "model_MICA_Genotypedata_1000G_mica_1000G_divide$training.RData")]
models_micb <- models_micb[c("model_MICB_geno_FIN_I_micbtab$training.RData","model_MICB_geno_comb_1000G_I_micbtab$training.RData", 
                             "model_MICB_geno_comb_1000G_I_EUR_iii_micb.RData", "model_MICB_geno_comb_1000G_I_micb_iv.RData",
                             "model_MICB_geno_comb_1000G_I_micb_1000G_divide$training.RData", "model_MICB_geno_comb_1000G_I_micb_vi.RData",
                             "model_MICB_Genotypedata_1000G_micb_1000G_divide$training.RData")]
models_hlae <- models_hlae[c("model_E_geno_FIN_I_hlaetab$training.RData","model_E_geno_comb_1000G_I_hlaetab$training.RData", 
                             "model_E_geno_comb_1000G_I_EUR_iii_e.RData", "model_E_geno_comb_1000G_I_hlae_iv.RData",
                             "model_E_geno_comb_1000G_I_hlae_1000G_divide$training.RData", "model_E_geno_comb_1000G_I_hlae_vi.RData",
                             "model_E_Genotypedata_1000G_hlae_1000G_divide$training.RData")]
models_hlaf <- models_hlaf[c("model_F_geno_FIN_II_hlaftab$training.RData","model_F_geno_comb_1000G_II_hlaftab$training.RData", 
                             "model_F_geno_comb_1000G_II_HSCT_EUR_iii_f.RData", "model_F_geno_comb_1000G_II_hlaf_iv.RData",
                             "model_F_geno_comb_1000G_II_hlaf_1000G_divide$training.RData", "model_F_geno_comb_1000G_II_hlaf_vi.RData",
                             "model_F_Genotypedata_1000G_hlaf_1000G_divide$training.RData")]

# load HLA-G models (created in Model_training.R)
model_hlag <- hlaModelFromObj(get(load("./data/Imputation_models/Combination_testing_models/HLA_G/model_G_geno_FIN_I_hlagtab$training.RData")))
model_hlag_3utr <- hlaModelFromObj(get(load("./data/Imputation_models/Combination_testing_models/HLA_G/model_G_geno_FIN_I_hlag3UTRtab$training.RData")))
model_hlag_5utr <- hlaModelFromObj(get(load("./data/Imputation_models/Combination_testing_models/HLA_G/model_G_geno_FIN_I_hlag5UTRtab$training.RData")))

### ............................................................................................................................................... ###

#         LOAD FIN AND 1000G PLINK GENOTYPE DATA


geno_FIN_I <- hlaBED2Geno(bed.fn = "data/Genotype_data/MHC_selected_samples_MICAB.bed", 
                       fam.fn = "data/Genotype_data/MHC_selected_samples_MICAB.fam", 
                       bim.fn = "data/Genotype_data/MHC_selected_samples_MICAB.bim", assembly="hg38")
geno_FIN_II <- hlaBED2Geno(bed.fn = "data/Genotype_data/MHCmergedunion_donors_hg38.bed", 
                         fam.fn = "data/Genotype_data/MHCmergedunion_donors_hg38.fam", 
                         bim.fn = "data/Genotype_data/MHCmergedunion_donors_hg38.bim", assembly="hg38")

# remove duplicated sample ID from FIN II genotype data set
geno_FIN_II$sample.id <- str_split_fixed(geno_FIN_II$sample.id, '-', 2)[, 1]

# Load 1000G plink genotype data and convert positions hg37 -> hg38

Genotypedata_1000G <- hlaBED2Geno(bed.fn = "data/Genotype_data/plink_SNP_data/chr6_phase3_MHC_rem_dup.bed", 
                                  fam.fn = "data/Genotype_data/plink_SNP_data/chr6_phase3_MHC_rem_dup.fam", 
                                  bim.fn = "data/Genotype_data/plink_SNP_data/chr6_phase3_MHC_rem_dup.bim", assembly="hg19")
hg38_converted_1000G_rem_dup_snps <- read.table("./data/Genotype_data/report_1000Gsnps_rem_dup_hg37b.txt", sep = "\t", header = TRUE)

# hg37 -> hg38 conversion done using https://www.ncbi.nlm.nih.gov/genome/tools/remap
hg38_converted_1000G_rem_dup_snps <- read.table("./data/Genotype_data/plink_SNP_data/report_1000Gsnps_rem_dup_hg37b.txt", sep = "\t", header = TRUE)
Genotypedata_1000G$snp.position <- hg38_converted_1000G_rem_dup_snps$mapped_start
Genotypedata_1000G$assembly <- 'hg38'

### ................................................................................................................................................ ###

#               IMPUTE FIN AND 1000G TEST SETS USING MODELS I - VII


## impute FIN reference test samples using all models I - VII
pred_FIN_mica <- pred_FIN(models_mica, geno_FIN_I, micatab$validation)
pred_FIN_micb <- pred_FIN(models_micb, geno_FIN_I, micbtab$validation)
pred_FIN_hlae <- pred_FIN(models_hlae, geno_FIN_I, hlaetab$validation)
pred_FIN_hlaf <- pred_FIN(models_hlaf, geno_FIN_II, hlaftab$validation)

## impute 1000G superpopulation (EUR, AFR, EAS, SAS, AMR) test sets

# list superpopulation test sets
val_superpop_mica <- list(EUR_mica, AFR_mica, EAS_mica, SAS_mica, AMR_mica) # phenotype sets from Model_training.R
val_superpop_micb <- list(EUR_micb, AFR_micb, EAS_micb, SAS_micb, AMR_micb)
val_superpop_hlae <- list(EUR_e, AFR_e, EAS_e, SAS_e, AMR_e)
val_superpop_hlaf <- list(EUR_f, AFR_f, EAS_f, SAS_f, AMR_f)

# impute 1000G test samples using all models i - vii
pred_1000G_mica <- pred_1000G_superpop(models_mica, val_superpop_mica)
pred_1000G_micb <- pred_1000G_superpop(models_micb, val_superpop_micb)
pred_1000G_hlae <- pred_1000G_superpop(models_hlae, val_superpop_hlae)
pred_1000G_hlaf <- pred_1000G_superpop(models_hlaf, val_superpop_hlaf)

# combine FIN and 1000G results into same list
pred_1000G_FIN_mica <- lapply(1:length(pred_1000G_mica),function(i) append(pred_1000G_mica[[i]],pred_FIN_mica[i]))
pred_1000G_FIN_micb <- lapply(1:length(pred_1000G_micb),function(i) append(pred_1000G_micb[[i]],pred_FIN_micb[i]))
pred_1000G_FIN_hlae <- lapply(1:length(pred_1000G_hlae),function(i) append(pred_1000G_hlae[[i]],pred_FIN_hlae[i]))
pred_1000G_FIN_hlaf <- lapply(1:length(pred_1000G_hlaf),function(i) append(pred_1000G_hlaf[[i]],pred_FIN_hlaf[i]))

# add model names
names(pred_1000G_FIN_mica) <- c("i", "ii", "iii", "iv", "v", "vi", "vii")
names(pred_1000G_FIN_micb) <- c("i", "ii", "iii", "iv", "v", "vi", "vii")
names(pred_1000G_FIN_hlae) <- c("i", "ii", "iii", "iv", "v", "vi", "vii")
names(pred_1000G_FIN_hlaf) <- c("i", "ii", "iii", "iv", "v", "vi", "vii")

# add population names
for (i in 1:length(pred_1000G_FIN_mica)) { names(pred_1000G_FIN_mica[[i]]) <-  c("EUR", "AFR", "EAS", "SAS", "AMR", "FIN")}

for (i in 1:length(pred_1000G_FIN_micb)) { names(pred_1000G_FIN_micb[[i]]) <-  c("EUR", "AFR", "EAS", "SAS", "AMR", "FIN")}

for (i in 1:length(pred_1000G_FIN_hlae)) { names(pred_1000G_FIN_hlae[[i]]) <-  c("EUR", "AFR", "EAS", "SAS", "AMR", "FIN")}

for (i in 1:length(pred_1000G_FIN_hlaf)) { names(pred_1000G_FIN_hlaf[[i]]) <-  c("EUR", "AFR", "EAS", "SAS", "AMR", "FIN")}

# save results
saveRDS(pred_1000G_FIN_mica, "./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/MICA_all_models_1000G_BB_pred_compare")
saveRDS(pred_1000G_FIN_micb, "./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/MICB_all_models_1000G_BB_pred_compare")
saveRDS(pred_1000G_FIN_hlae, "./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/HLA_E_all_models_1000G_BB_pred_compare")
saveRDS(pred_1000G_FIN_hlaf, "./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/HLA_F_all_models_1000G_BB_pred_compare")

### impute HLA-G and HLA-G 3' and 5'UTRs (only model I) and save results
pred_hla_G <- PredHlaModel(model_hlag, geno_FIN_I, hlagtab$validation)
pred_hla_G_3UTR <- PredHlaModel(model_hlag_3utr, geno_FIN_I, hlag3UTRtab$validation)
pred_hla_G_5UTR <- PredHlaModel(model_hlag_5utr, geno_FIN_I, hlag5UTRtab$validation)

saveRDS(pred_hla_G, "./data/Output/Data_comb_test_results/HLA_G/HLA_G_BB_pred_compare")
saveRDS(pred_hla_G_3UTR, "./data/Output/Data_comb_test_results/HLA_G/HLA_G_3UTR_BB_pred_compare")
saveRDS(pred_hla_G_5UTR, "./data/Output/Data_comb_test_results/HLA_G/HLA_G_5UTR_BB_pred_compare")

###

# Impute all populations test sets (FIN + 1000G) together

# combine FIN and 1000G test sets
mica_1000G_BB_test <- hlaCombineAllele(mica_1000G_divide$validation, micatab$validation)
micb_1000G_BB_test <- hlaCombineAllele(micb_1000G_divide$validation, micbtab$validation)
hlae_1000G_BB_test <- hlaCombineAllele(hlae_1000G_divide$validation, hlaetab$validation)
hlaf_1000G_BB_test <- hlaCombineAllele(hlaf_1000G_divide$validation, hlaftab$validation)

# impute
pred_all_1000G_mica <- pred_FIN(models_mica, geno_comb_1000G_FIN_I, mica_1000G_BB_test)
pred_all_1000G_micb <- pred_FIN(models_micb, geno_comb_1000G_FIN_I, micb_1000G_BB_test)
pred_all_1000G_hlae <- pred_FIN(models_hlae, geno_comb_1000G_FIN_I, hlae_1000G_BB_test)
pred_all_1000G_hlaf <- pred_FIN(models_hlaf, geno_comb_1000G_FIN_II, hlaf_1000G_BB_test)

# add model names
names(pred_all_1000G_mica) <- c("i", "ii", "iii", "iv", "v", "vi", "vii")
names(pred_all_1000G_micb) <- c("i", "ii", "iii", "iv", "v", "vi", "vii")
names(pred_all_1000G_hlae) <- c("i", "ii", "iii", "iv", "v", "vi", "vii")
names(pred_all_1000G_hlaf) <- c("i", "ii", "iii", "iv", "v", "vi", "vii")

# save results
saveRDS(pred_all_1000G_mica, "./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb/MICA_1000G_superpop_comb_pred_compare")
saveRDS(pred_all_1000G_micb, "./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb/MICB_1000G_superpop_comb_pred_compare")
saveRDS(pred_all_1000G_hlae, "./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb/HLA_E_1000G_superpop_comb_pred_compare")
saveRDS(pred_all_1000G_hlaf, "./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb/HLA_F_1000G_superpop_comb_pred_compare")

### ................................................................................................................................................ ###

#           CROSS-VALIDATION OF 1000G REFERENCE WITH FINNISH REFERENCE


# impute all FIN refrence using 1000G model (VII)
# phenotype sets for all data created in Model_training.R

pred_all_FIN_I_mica <- PredHlaModel(models_mica$`model_MICA_Genotypedata_1000G_mica_1000G_divide$training.RData` , geno_FIN_I, mica_FIN_I)
pred_all_FIN_I_micb <- PredHlaModel(models_micb$`model_MICB_Genotypedata_1000G_micb_1000G_divide$training.RData` , geno_FIN_I, micb_FIN_I)
pred_all_FIN_I_hlae <- PredHlaModel(models_hlae$`model_E_Genotypedata_1000G_hlae_1000G_divide$training.RData` , geno_FIN_I, HLA_E_FIN_I)
pred_all_FIN_II_hlaf <- PredHlaModel(models_hlaf$`model_F_Genotypedata_1000G_hlaf_1000G_divide$training.RData` , geno_FIN_II, HLA_F_FIN_II)

# saveRDS(pred_all_FIN_I_mica, "./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/FIN_all/MICA_BB_all_pred_compare")

