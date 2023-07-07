#####-------------------------------------------------------------------------------------------------------------------#####

# FINALIZE MODELS TRAINED WITH ALL DATA

#####-------------------------------------------------------------------------------------------------------------------#####

library(HIBAG)

### CHANGE NAMES OF NOVEL MICA AND MICB ALLELES TO MODELS TRAINED WITH ALL DATA

# Novel alleles

# MICA*R29H -> MICA*285
# MICA*D226N = D249N in the model -> MICA*193
# MICB*G198R -> no submission, remaines named as it is in the model
# MICB*A337V = A314V in the model -> MICB*053
# MICB*R29C = R6C in the model -> MICB*052

# novels alleles only present in Finnish reference -> change allele names to models (i), (ii) and (iv) and (vi), here MICA and MICB model i as an example
# mica
model_mica_all <- hlaModelFromObj(get(load("./data/Imputation_models/Combination_testing_models/All_data/MICA/model_MICA_geno_FIN_I_mica_Histog.RData")))
model_mica_all$hla.allele
model_mica_all$hla.freq

model_mica_all$hla.allele[model_mica_all$hla.allele == 'R6H'] <- '285'
names(model_mica_all$hla.freq)[names(model_mica_all$hla.freq) == 'R6H'] <- '285'
model_mica_all$hla.allele[model_mica_all$hla.allele == 'D226N'] <- '193'
names(model_mica_all$hla.freq)[names(model_mica_all$hla.freq) == 'D226N'] <- '193'
model_mica_all$hla.allele
model_mica_all$hla.freq

model_mica_all <- hlaModelToObj(model_mica_all)
save(model_mica_all, file="./data/Imputation_models/Combination_testing_models/All_data/MICA/model_MICA_geno_FIN_I_mica_Histog_final.RData")
model_test <- hlaModelFromObj(get(load("./data/Imputation_models/Combination_testing_models/All_data/MICA/model_MICA_geno_FIN_I_mica_Histog_final.RData")))
model_test$hla.allele
model_test$hla.freq

# micb

model_micb_all <- hlaModelFromObj(get(load("./data/Imputation_models/Combination_testing_models/All_data/MICB/model_MICB_geno_FIN_I_micb_Histog.RData")))
summary(model_micb_all)
model_micb_all$hla.allele
model_micb_all$hla.freq

model_micb_all$hla.allele[model_micb_all$hla.allele == 'A314V'] <- '053'
names(model_micb_all$hla.freq)[names(model_micb_all$hla.freq) == 'A314V'] <- '053'
model_micb_all$hla.allele[model_micb_all$hla.allele == 'R6C'] <- '052'
names(model_micb_all$hla.freq)[names(model_micb_all$hla.freq) == 'R6C'] <- '052'
model_micb_all$hla.allele[model_micb_all$hla.allele == '009N'] <- '009:01N'
names(model_micb_all$hla.freq)[names(model_micb_all$hla.freq) == '009N'] <- '009:01N'
model_micb_all$hla.allele
model_micb_all$hla.freq

model_micb_all <- hlaModelToObj(model_micb_all)
save(model_micb_all, file="./data/Imputation_models/Combination_testing_models/All_data/MICB/model_MICB_geno_FIN_I_micb_Histog_final.RData")
model_test <- hlaModelFromObj(get(load("./data/Imputation_models/Combination_testing_models/All_data/MICB/model_MICB_geno_FIN_I_micb_Histog_final.RData")))
model_test$hla.allele

##### ------------------------------------------------------------------------------------------------------------------------------------------- #####

### FINALIZE WITH HIBAG hlaPublish (delete sample IDs and unused SNPs and add change model name)


model_to_finalize <- hlaModelFromObj(get(load("./data/Imputation_models/Combination_testing_models/All_data/MICA/model_mica_geno_FIN_I_mica_Histog_final.RData")))
model_to_finalize$hla.allele
model_to_finalize$hla.freq
summary(model_to_finalize)
model_to_finalize$n.snp

mobj <- hlaPublish(model_to_finalize, rm.unused.snp=T, anonymize=T)
mobj$n.snp
mobj$sample.id
summary(mobj)

save(mobj, file="./data/Imputation_models/Combination_testing_models/All_data/MICA/finals/MICA_model_I.RData")
model_test <- hlaModelFromObj(get(load("./data/Imputation_models/Combination_testing_models/All_data/MICA/finals/MICA_model_I.RData")))
summary(model_test)
model_test$sample.id
model_test$hla.allele
model_test$hla.freq



