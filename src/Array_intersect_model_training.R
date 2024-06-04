### ................................................................................................................................................ ###

### fit HIBAG models for MICA, MICB, HLA-E, HLA-F and HLA-G for GSA and PMRA array SNPs

### ................................................................................................................................................ ###

# libraries and functions

library(tidyverse)
library(HIBAG)
library(parallel)
library(data.table)
library(readxl)

# load functions
source('./src/functions.R')

# load genotype and phenotype data
source('./src/Geno_pheno_loading.R')

### ................................................................................................... # 

### FIT MODELS 

# .................................................................................................... #

# ILLUMINA GSA

# .................................................................................................... #

# GSA SNP list
GSA_snps <- read.table("./data/Kekalainen/KEK_MHC_pos.txt")

# remove 6_ from GSA snp list and choose the SNPs that are common in GSA and FinnGen
GSA_snps <- GSA_snps %>% mutate(V1=gsub("6_","", V1))
nrow(GSA_snps) # 8635

# check how many intersect snps there are
length(which(GSA_snps$V1 %in% geno_comb_1000G_I$snp.position)) # 6955
length(which(GSA_snps$V1 %in% geno_comb_1000G_II$snp.position)) # 7040
length(which(GSA_snps$V1 %in% geno_FIN_I$snp.position)) # 7002

# extract 1000G/FIN reference - array combination intersect snps from genotype data
geno_GSAintersect_I <- hlaGenoSubset(geno_comb_1000G_I, snp.sel= geno_comb_1000G_I$snp.position %in% GSA_snps$V1) # used for model II (HLA-E, -F, MICA, MICB)
geno_GSAintersect_II <- hlaGenoSubset(geno_comb_1000G_II, snp.sel= geno_comb_1000G_II$snp.position %in% GSA_snps$V1) # used for model II (HLA-E, -F, MICA, MICB)
geno_GSAintersect_FIN <- hlaGenoSubset(geno_FIN_I, snp.sel= geno_FIN_I$snp.position %in% GSA_snps$V1)


# train MICA, MICB, HLA-E and HLA-F GSA models with training data
map2(c('MICA', 'MICB', 'E', 'F'), names(divide_comb_train), function(x,y) {
   if(x == 'F') {
     fitHLAmodel(x, divide_comb_train[[y]], geno_GSAintersect_II, 500000, 100)
   } else {
     fitHLAmodel(x, divide_comb_train[[y]], geno_GSAintersect_I, 500000, 100)
   }
})

# train HLA-G, HLA-G 3'UTR and HLA-G 5'UTR GSA models with training data
map(c('HLA_G', 'HLA_G_3UTR', 'HLA_G_5UTR'), function(x) {
    fitHLAmodel('G', divide_FIN[[x]], geno_GSAintersect_FIN, 500000, 100)
})

# ................................................................................................ #

###  validation of GSA models

# import models for MICA, MICB, HLA-E and HLA-F
model_list <- list.files("./data/Imputation_models/Combination_testing_models/Array_intersect/GSA/Training", pattern = 'divide_comb', full.names = F)
model_list_g <- list.files("./data/Imputation_models/Combination_testing_models/Array_intersect/GSA/Training", pattern = 'divide_FIN', full.names = F)

models_GSAint <- sapply(model_list, function(x) {
  list(hlaModelFromObj(get(load(paste0("./data/Imputation_models/Combination_testing_models/Array_intersect/GSA/Training/", x)))))
})

models_GSAint_g <- sapply(model_list_g, function(x) {
  list(hlaModelFromObj(get(load(paste0("./data/Imputation_models/Combination_testing_models/Array_intersect/GSA/Training/", x)))))
})

names(superpop_val) <- names(superpop_val) %>% gsub('HLA-', '', .)
names <- names(models_GSAint) %>% gsub('model_|_geno_GSAintersect_I_divide_comb_train_E.RData', '', .) %>%
  gsub('_geno_GSAintersect_II_divide_comb_train_F.RData|_geno_GSAintersect_I_divide_comb_train_MICA.RData', '', .) %>%
  gsub('_geno_GSAintersect_I_divide_comb_train_MICB.RData', '', .)
names_g <- names(models_GSAint_g) %>% gsub('model_|_geno_GSAintersect_FIN_divide_FIN_G.RData', '', .) %>%
  gsub('G_geno_GSAintersect_FIN_divide_FIN_G_', '', .) %>% gsub('.RData', '', .)

names(divide_FIN) <- gsub('HLA_', '', names(divide_FIN))
names(divide_comb_val) <- gsub('HLA_', '', names(divide_comb_val))

# combine FIN and ALL to superpop_val validation sets
for (i in names(superpop_val)) {
  superpop_val[[i]]$FIN <- divide_FIN[[i]]$validation
  superpop_val[[i]]$ALL <- divide_comb_val[[i]]
}

# pred all populations using all models
pred_superpop <- map2(1:length(models_GSAint), names, function(x,y) {
  map(1:length(superpop_val[[1]]), function(z) {
    if(y == 'F') {
      PredHlaModel(models_GSAint[[x]], geno_GSAintersect_II, superpop_val[[x]][[z]])
    } else {
      PredHlaModel(models_GSAint[[x]], geno_GSAintersect_I, superpop_val[[x]][[z]])
    }
  })
})

# set names
all_pop <- c(superpop, 'FIN', 'ALL')
names(pred_superpop) <- names
for (i in 1:length(pred_superpop)) {
  names(pred_superpop[[i]]) <- all_pop
}

# save results
# saveRDS(pred_superpop, "./data/Imputation_models/Combination_testing_models/Array_intersect/GSA/Training/pred_superpop.RDS")

# hla-G
names_g <- c(paste0('G_', names_g[1:2]), names_g[3])
pred_fin_g <- map2(1:length(models_GSAint_g), names_g, function(x, y) {
  PredHlaModel(models_GSAint_g[[x]], geno_GSAintersect_FIN, divide_FIN[[y]]$validation)
})

names(pred_fin_g) <- names_g
# saveRDS(pred_fin_g, "./data/Imputation_models/Combination_testing_models/Array_intersect/pred_superpop_hlag.RDS")

# .................................................................................................. #

# fit GSA intersect models with all data

# .................................................................................................. #

# train MICA, MICB, HLA-E and HLA-F GSA models with all data
map2(c('MICA', 'MICB', 'E', 'F'), names(divide_comb_train), function(x,y) {
   if(x == 'F') {
     fitHLAmodel(x, divide_comb[[y]], geno_GSAintersect_II, 500000, 100)
   } else {
     fitHLAmodel(x, divide_comb[[y]], geno_GSAintersect_I, 500000, 100)
   }
})

# train HLA-G, HLA-G 3'UTR and HLA-G 5'UTR GSA models with all data
map(c('HLA_G', 'HLA_G_3UTR', 'HLA_G_5UTR'), function(x) {
  fitHLAmodel('G', divide_FIN_all[[x]], geno_GSAintersect_FIN, 500000, 100)
})

# ................................................................................................ #

# Create validation result and model parameter tables

# GSA models trained with training data
model_list <- list.files("./data/Imputation_models/Combination_testing_models/Array_intersect/GSA/Training", pattern = '.RData', full.names = F)

models_GSAint <- sapply(model_list, function(x) {
  list(hlaModelFromObj(get(load(paste0("./data/Imputation_models/Combination_testing_models/Array_intersect/GSA/Training/", x)))))
})

names <- names(models_GSAint) %>% gsub('model_|_geno_GSAintersect_I_divide_comb_train_E.RData', '', .) %>%
  gsub('_geno_GSAintersect_II_divide_comb_train_F.RData|_geno_GSAintersect_I_divide_comb_train_MICA.RData', '', .) %>%
  gsub('_geno_GSAintersect_I_divide_comb_train_MICB.RData', '', .) %>% gsub('model_|_geno_GSAintersect_FIN_divide_FIN_G.RData', '', .) %>%
  gsub('G_geno_GSAintersect_FIN_divide_FIN_G_', '', .) %>% gsub('.RData', '', .)


# GSA models trained with all data
model_list <- list.files("./data/Imputation_models/Combination_testing_models/Array_intersect/GSA/All_data/Final", pattern = '.RData', full.names = F)

models_GSAint_all <- sapply(model_list, function(x) {
  list(hlaModelFromObj(get(load(paste0("./data/Imputation_models/Combination_testing_models/Array_intersect/GSA/All_data/Final/", x)))))
})

# extract model parameter info
model_param_gsa <- list()
for (model in names(models_GSAint)) {
  overall.temp.p <-data.frame(Training_n=models_GSAint[[model]]$n.samp, 
                              SNPs_total=models_GSAint[[model]]$n.snp, 
                              SNPs_used=summary(models_GSAint[[model]])$num.snp,
                              Alleles=length(models_GSAint[[model]]$hla.allele),
                              OOB_train=summary(models_GSAint[[model]])$info[3,1]) %>% mutate(Model=model)
  model_param_gsa[[length(model_param_gsa)+1]] <- overall.temp.p
}

# add all data OOB
for (i in 1:length(model_param_gsa)) {
  model_param_gsa[[i]] <- mutate(model_param_gsa[[i]], OOB_All_data=summary(models_GSAint_all[[i]])$info[3,1])
  model_param_gsa[[i]] <- mutate(model_param_gsa[[i]], All_data_n=models_GSAint_all[[i]]$n.samp)
  }

# combine and reorder columns
param.list.gsa <-do.call("rbind", model_param_gsa) 
param.list.gsa <- param.list.gsa[, c('Model', 'Training_n', 'All_data_n', 'SNPs_total', 'SNPs_used', 'Alleles', 'OOB_train', 'OOB_All_data')]
#fwrite(param.list.gsa, row.names = F, sep = "\t", "~/imputation/Non_classical_HLA_MICAB_imputation/data/Imputation_models/Combination_testing_models/Array_intersect/GSA/Training/Model_prop_GSA_cv_all.txt")

# .....

# overall accuracy table GSA

# validation results
pred_superpop_gsa <- readRDS("./data/Imputation_models/Combination_testing_models/Array_intersect/GSA/Training/pred_superpop_GSA.RDS")
pred_fin_g_gsa <- readRDS("./data/Imputation_models/Combination_testing_models/Array_intersect/GSA/Training/pred_superpop_hlag_GSA.RDS")

# extract population results
all_pop <- c('EUR', 'AFR', 'EAS', 'SAS', 'AMR', 'FIN')
overalls_gsa <- list()
for (gene in names(pred_superpop_gsa)) {
  for (pop in all_pop) {
    overall.temp <-pred_superpop_gsa[[gene]][[pop]]$overall %>% mutate(Gene=gene) %>% mutate(Pop=pop)
    overalls_gsa[[length(overalls_gsa)+1]] <- overall.temp
  }
}

# combine and reorder columns
overall.list.gsa <-do.call("rbind", overalls_gsa) 
overall.list.gsa <- overall.list.gsa[, c(ncol(overall.list.gsa)-1, ncol(overall.list.gsa), 5, 4, 1:3, 6:(ncol(overall.list.gsa)-2))]
# fwrite(overall.list.gsa, row.names = F, sep = "\t", "~/imputation/Non_classical_HLA_MICAB_imputation/data/Imputation_models/Combination_testing_models/Array_intersect/Overall_acc_GSA_cv.txt")

# hla-g
overalls_g_gsa <- list()
for (gene in names(pred_fin_g_gsa)) {
  overall.temp.g <-pred_fin_g_gsa[[gene]]$overall %>% mutate(Gene=gene) %>% mutate(Pop='FIN')
  overalls_g_gsa[[length(overalls_g_gsa)+1]] <- overall.temp.g
}

overall.list.g.gsa <- do.call("rbind", overalls_g_gsa)
overall.list.g.gsa <- overall.list.g.gsa[, c(ncol(overall.list.g.gsa)-1, ncol(overall.list.g.gsa), 5, 4, 1:3, 6:(ncol(overall.list.g.gsa)-2))]
# fwrite(overall.list.g.gsa, row.names = F, sep = "\t", "~/imputation/Non_classical_HLA_MICAB_imputation/data/Imputation_models/Combination_testing_models/Array_intersect/Overall_acc_GSA_cv_hlag.txt")

# combine all loci - overall acc results
overall.list.gsa.all <- rbind(overall.list.gsa, overall.list.g.gsa)
# fwrite(overall.list.gsa.all, row.names = F, sep = "\t", "~/imputation/Non_classical_HLA_MICAB_imputation/data/Imputation_models/Combination_testing_models/Array_intersect/GSA/Training/Overall_acc_GSA_cv_all.txt")


### .................................................................................................... #

### PMRA 

### .................................................................................................... #

PMRA_snps <- fread("./data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training/MHC_SNPs_PMRA.tsv")
length(PMRA_snps$hg38) # 9563

# check how many intersect snps there are
length(which(PMRA_snps$hg38 %in% geno_comb_1000G_I$snp.position)) # 6871
length(which(PMRA_snps$hg38 %in% geno_comb_1000G_II$snp.position)) # 7296
length(which(PMRA_snps$hg38 %in% geno_FIN_I$snp.position)) # 6930

# extract 1000G/FIN reference - array combination intersect snps from genotype data
geno_PMRAintersect_I <- hlaGenoSubset(geno_comb_1000G_I, snp.sel= geno_comb_1000G_I$snp.position %in% PMRA_snps$hg38) # used for HLA-E, MICA, MICB
geno_PMRAintersect_II <- hlaGenoSubset(geno_comb_1000G_II, snp.sel= geno_comb_1000G_II$snp.position %in% PMRA_snps$hg38) # used for HLA-F
geno_PMRAintersect_FIN <- hlaGenoSubset(geno_FIN_I, snp.sel= geno_FIN_I$snp.position %in% PMRA_snps$hg38)


# fit MICA, MICB, HLA-E and HLA-F PMRA models with training data
map2(c('MICA', 'MICB', 'E', 'F'), names(divide_comb_train), function(x,y) {
  if(x == 'F') {
    fitHLAmodel(x, divide_comb_train[[y]], geno_PMRAintersect_II, 500000, 100)
  } else {
    fitHLAmodel(x, divide_comb_train[[y]], geno_PMRAintersect_I, 500000, 100)
  }
})

# train HLA-G, HLA-G 3'UTR and HLA-G 5'UTR PMRA models with training data
map(c('HLA_G', 'HLA_G_3UTR', 'HLA_G_5UTR'), function(x) {
  fitHLAmodel('G', divide_FIN[[x]], geno_PMRAintersect_FIN, 500000, 100)
})

##### ....................................................................................

###  validation of PMRA models

# import models for MICA, MICB, HLA-E and HLA-F

model_list_pmra <- list.files("./data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training", pattern = 'divide_comb', full.names = F)
model_list_pmra_g <- list.files("./data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training", pattern = 'divide_FIN', full.names = F)

models_PMRAint <- sapply(model_list_pmra, function(x) {
  list(hlaModelFromObj(get(load(paste0("./data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training/", x)))))
})

models_PMRAint_g <- sapply(model_list_pmra_g, function(x) {
  list(hlaModelFromObj(get(load(paste0("./data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training/", x)))))
})

# predict and cross-validate all population test sets
names(superpop_val) <- names(superpop_val) %>% gsub('HLA-', '', .)
names <- names(models_PMRAint) %>% gsub('model_|_geno_PMRAintersect_I_divide_comb_train_E.RData', '', .) %>%
  gsub('_geno_PMRAintersect_II_divide_comb_train_F.RData|_geno_PMRAintersect_I_divide_comb_train_MICA.RData', '', .) %>%
  gsub('_geno_PMRAintersect_I_divide_comb_train_MICB.RData', '', .)
names_g <- names(models_PMRAint_g) %>% gsub('model_G_geno_PMRAintersect_FIN_divide_FIN_', '', .) %>%
  gsub('_train.RData', '', .) %>% gsub('.RData', '', .)

names(divide_FIN) <- gsub('HLA_', '', names(divide_FIN))
names(divide_comb_val) <- gsub('HLA_', '', names(divide_comb_val))

# combine FIN and ALL to superpop_val validation sets
for (i in names(superpop_val)) {
  superpop_val[[i]]$FIN <- divide_FIN[[i]]$validation
  superpop_val[[i]]$ALL <- divide_comb_val[[i]]
}

# pred all populations using all models
pred_superpop_pmra <- map2(1:length(models_PMRAint), names, function(x,y) {
  
  map(1:length(superpop_val[[1]]), function(z) {
    if(y == 'F') {
      PredHlaModel(models_PMRAint[[x]], geno_PMRAintersect_II, superpop_val[[x]][[z]])
    } else {
      PredHlaModel(models_PMRAint[[x]], geno_PMRAintersect_I, superpop_val[[x]][[z]])
    }
  })
})

# set names
all_pop <- c(superpop, 'FIN', 'ALL')
names(pred_superpop_pmra) <- names
for (i in 1:length(pred_superpop_pmra)) {
  names(pred_superpop_pmra[[i]]) <- all_pop
}

# save results
saveRDS(pred_superpop_pmra, "./data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training/pred_superpop_PMRA.RDS")

# hla-G
#names_g <- c(paste0('G_', names_g[1:2]), names_g[3])
pred_fin_g_pmra <- map2(1:length(models_PMRAint_g), names_g, function(x, y) {
  PredHlaModel(models_PMRAint_g[[x]], geno_PMRAintersect_FIN, divide_FIN[[y]]$validation)
})

names(pred_fin_g_pmra) <- names_g
saveRDS(pred_fin_g_pmra, "./data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training/pred_superpop_hlag_PMRA.RDS")

# ................................................................................................ #

# Create validation result and model parameter tables

# Model parameter table

# PMRA models trained with training data
model_list <- list.files("./data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training", pattern = '.RData', full.names = F)

models_PMRAint <- sapply(model_list, function(x) {
  list(hlaModelFromObj(get(load(paste0("./data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training/", x)))))
})

names <- names(models_PMRAint) %>% gsub('model_|_geno_PMRAintersect_I_divide_comb_train_E.RData', '', .) %>%
  gsub('_geno_PMRAintersect_II_divide_comb_train_F.RData|_geno_PMRAintersect_I_divide_comb_train_MICA.RData', '', .) %>%
  gsub('_geno_PMRAintersect_I_divide_comb_train_MICB.RData', '', .) %>% gsub('model_|_geno_PMRAintersect_FIN_divide_FIN_G.RData', '', .) %>%
  gsub('G_geno_PMRAintersect_FIN_divide_FIN_G_', '', .) %>% gsub('.RData', '', .)


# PMRA models trained with all data
model_list <- list.files("./data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/All_data/Final", pattern = '.RData', full.names = F)

models_PMRAint_all <- sapply(model_list, function(x) {
  list(hlaModelFromObj(get(load(paste0("./data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/All_data/Final/", x)))))
})

# extract model parameter info
model_param_pmra <- list()
for (model in names(models_PMRAint)) {
  overall.temp.p <-data.frame(Training_n=models_PMRAint[[model]]$n.samp, 
                              SNPs_total=models_PMRAint[[model]]$n.snp, 
                              SNPs_used=summary(models_PMRAint[[model]])$num.snp,
                              Alleles=length(models_PMRAint[[model]]$hla.allele),
                              OOB_train=summary(models_PMRAint[[model]])$info[3,1]) %>% mutate(Model=model)
  model_param_pmra[[length(model_param_pmra)+1]] <- overall.temp.p
}

# add all data OOB
for (i in 1:length(model_param_pmra)) {
  model_param_pmra[[i]] <- mutate(model_param_pmra[[i]], OOB_All_data=summary(models_PMRAint_all[[i]])$info[3,1])
  model_param_pmra[[i]] <- mutate(model_param_pmra[[i]], All_data_n=models_PMRAint_all[[i]]$n.samp)
}

# combine and reorder columns
param.list.pmra <-do.call("rbind", model_param_pmra) 
param.list.pmra <- param.list.pmra[, c('Model', 'Training_n', 'All_data_n', 'SNPs_total', 'SNPs_used', 'Alleles', 'OOB_train', 'OOB_All_data')]
# fwrite(param.list.pmra, row.names = F, sep = "\t", "~/imputation/Non_classical_HLA_MICAB_imputation/data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training/Model_prop_PMRA_cv_all.txt")

# ...

# overall accuracy table
pred_superpop_pmra <- readRDS("./data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training/pred_superpop_PMRA.RDS")
pred_fin_g_pmra <- readRDS("./data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training/pred_superpop_hlag_PMRA.RDS")

# extract population results
overalls_pmra <- list()
for (gene in names(pred_superpop_pmra)) {
  for (pop in all_pop) {
    overall.temp <-pred_superpop_pmra[[gene]][[pop]]$overall %>% mutate(Gene=gene) %>% mutate(Pop=pop)
    overalls_pmra[[length(overalls_pmra)+1]] <- overall.temp
  }
}

# combine and reorder columns
overall.list.pmra <-do.call("rbind", overalls_pmra) 
overall.list.pmra <- overall.list.pmra[, c(ncol(overall.list_pmra)-1, ncol(overall.list.pmra), 5, 4, 1:3, 6:(ncol(overall.list_pmra)-2))]
# fwrite(overall.list.pmra, row.names = F, sep = "\t", "~/imputation/Non_classical_HLA_MICAB_imputation/data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training/Overall_acc_PMRA_cv.txt")

# hla-g
overalls_g_pmra <- list()
for (gene in names(pred_fin_g_pmra)) {
  overall.temp.g <-pred_fin_g_pmra[[gene]]$overall %>% mutate(Gene=gene) %>% mutate(Pop='FIN')
  overalls_g_pmra[[length(overalls_g_pmra)+1]] <- overall.temp.g
}


overall.list.g.pmra <- do.call("rbind", overalls_g_pmra)
overall.list.g.pmra <- overall.list.g[, c(ncol(overall.list.g.pmra)-1, ncol(overall.list.g.pmra), 5, 4, 1:3, 6:(ncol(overall.list.g.pmra)-2))]
# fwrite(overall.list.g.pmra, row.names = F, sep = "\t", "~/imputation/Non_classical_HLA_MICAB_imputation/data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training/Overall_acc_PMRA_cv_hlag.txt")


# combine all loci - overall acc results
overall.list.all.pmra <- rbind(overall.list.pmra, overall.list.g.pmra)
# fwrite(overall.list.all.pmra, row.names = F, sep = "\t", "~/imputation/Non_classical_HLA_MICAB_imputation/data/Imputation_models/Combination_testing_models/Array_intersect/PMRA/Training/Overall_acc_PMRA_cv_all.txt")

# .................................................................................................. #

# fit PMRA intersect models with all data

# .................................................................................................. #


# train MICA, MICB, HLA-E and HLA-F PMRA models with all data
map2(c('MICA', 'MICB', 'E', 'F'), names(divide_comb_train), function(x,y) {
   if(x == 'F') {
     fitHLAmodel(x, divide_comb[[y]], geno_PMRAintersect_II, 500000, 100)
   } else {
     fitHLAmodel(x, divide_comb[[y]], geno_PMRAintersect_I, 500000, 100)
   }
})

# train HLA-G, HLA-G 3'UTR and HLA-G 5'UTR PMRA models with all data
map(c('HLA_G', 'HLA_G_3UTR', 'HLA_G_5UTR'), function(x) {
  fitHLAmodel('G', divide_FIN_all[[x]], geno_PMRAintersect_FIN, 500000, 100)
})
