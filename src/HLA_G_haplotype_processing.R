##############################################################################################################################
## HLA-G 3'UTR and 5'UTR haplotype imputation
##############################################################################################################################

library(data.table)
library(readxl)
library(tidyverse)
library(HIBAG)
library(parallel)

##############################################################################################################################
### Conversion of binary haplotype data to literal
##############################################################################################################################

# read phased haplotype data and sampleID information from .haps and .sample files
haplo_3UTR <- read.table("./data/Genotype_data/3UTR_SNPs_haplo.haps")
haplo_5UTR <- read.table("./data/Genotype_data/5UTR_SNPs_haplo.haps")
samples_UTRs <- read.table("./data/Genotype_data/3UTR_SNPs_haplo.sample")

# get ref and alt alleles of haplotype SNPs
letters_3 <- haplo_3UTR[, 1:5] #data.frame(haplo_3UTR$V1, haplo_3UTR$V2, haplo_3UTR$V3, haplo_3UTR$V4, haplo_3UTR$V5)
colnames(letters_3) <- c("chrom","snp_ID", "pos", "A1", "A2")
letters_5 <- haplo_5UTR[, 1:5] #data.frame(haplo_5UTR$V1, haplo_5UTR$V2, haplo_5UTR$V3, haplo_5UTR$V4, haplo_5UTR$V5)
colnames(letters_5) <- c("chrom","snp_ID", "pos", "A1", "A2")

# extract columns that contain actual haplotypes
haplo_3UTR <- haplo_3UTR %>% select(-(1:5))
haplo_5UTR <- haplo_5UTR %>% select(-(1:5))

# convert binary haplotype data into literal
# 3' UTR
haplos_literal_3UTR <- sapply(1:nrow(haplo_3UTR), function(i) {
  ifelse(unlist(haplo_3UTR[i, ])==0, letters_3[i, 'A1'], letters_3[i, 'A2'])
})

# 5' UTR
haplos_literal_5UTR <- sapply(1:nrow(haplo_5UTR), function(i) {
  ifelse(unlist(haplo_5UTR[i, ])==0, letters_5[i, 'A1'], letters_5[i, 'A2'])
})

# double sampleID rows to mach haplotype data before combining with haplotype data
samples_UTRs_doubled <- samples_UTRs[rep(seq_len(nrow(samples_UTRs)), each=2), ]

# extract rows with actual sampleIDs
samples_UTRs_doubled <- samples_UTRs_doubled[c(5:(nrow(samples_UTRs_doubled))), ]

# Combine sampleIDs and converted haplotypes
combine_3 <- cbind(samples_UTRs_doubled, haplos_literal_3UTR)
combine_5 <- cbind(samples_UTRs_doubled, haplos_literal_5UTR)

# extract samples that are in phenotype data
HLA_G_pheno_Histog <- read.table("./data/Phenotype_data/HLA_G_alleles_Histog")
chosen_samples_3 <- subset(combine_3, combine_3$V2 %in% HLA_G_pheno_Histog$SampleID)
chosen_samples_5 <- subset(combine_5, combine_5$V2 %in% HLA_G_pheno_Histog$SampleID)

# select columns with sample ID and haplotype information
chosen_samples_3 <- chosen_samples_3 %>% select(2,5:ncol(chosen_samples_3))
chosen_samples_5 <- chosen_samples_5 %>% select(2,5:ncol(chosen_samples_5))

# change column names (positions are ordered)
colnames(chosen_samples_3) <- c("SampleID", letters_3$snp_ID %>% sort %>% as.character)
colnames(chosen_samples_5) <- c("SampleID", letters_5$snp_ID %>% sort %>% as.character)

# create 3'UTR haplotype file and continue to process 5'UTR multiallelic SNP
write.table(chosen_samples_3, sep = "\t", "./data/Genotype_data/HLA_G_3_UTR_haplotypes")

#########################################################################################################################################################################
### processing of multiallelic 5'UTR SNP (pos 29827120)
#########################################################################################################################################################################

# make vectors of 29827120_G_C and 29827120_G_T
a <- vector()
for (i in chosen_samples_5$chr6_29827120_G_C) {
  ifelse (i=="G", a <- c(a, paste0("G")), a <- c(a, paste0("C")))
}

b <- vector()
for (i in chosen_samples_5$chr6_29827120_G_T) {
  ifelse (i=="T", b <- c(b, paste0("T")), b <- c(b, paste0("G")))
}

# combine vectors and replace with correct genotype
c <- paste0(a,b)
d <- vector()
for (i in c) {
  if (i=="GT") {
    d <- c(d, paste0("T"))
  } else if (i=="CG") {
    d <- c(d, paste0("C"))
  } else {
    d <- c(d, paste0("G"))
  }
}

## Add processed multiallelic SNP to data
# remove two columns of 29827120
chosen_samples_5_without_multiall_SNP <- chosen_samples_5[, -c(10:11)]

# insert "d" containing correct data (G,C or T) for pos 29827120
chosen_samples_5_corrected <- cbind(chosen_samples_5_without_multiall_SNP, d)

#rearrange columns, add column name and create 5'UTR haplotype file
chosen_samples_5_corrected <- chosen_samples_5_corrected[, c(1:9, 28, 10:27)]
names(chosen_samples_5_corrected)[names(chosen_samples_5_corrected)=="d"] <- "29827120_G_C_T"
write.table(chosen_samples_5_corrected, sep = "\t", "./data/Genotype_data/HLA_G_5_UTR_haplotypes")

##############################################################################################################################################
### Inference of 3'UTR and 5'UTR haplotypes
##############################################################################################################################################

# Read haplotype models (according to Castelli et al. 2014), remove unnecessary columns, name columns and change columns into rows
model_3UTR <- read_excel("./data/Genotype_data/HLA_G_UTR_haplotyypit_imputaatio_malli.xlsx", sheet = "3'UTR")
model_3UTR <- model_3UTR[, 2:ncol(model_3UTR)]
model_3UTR <- rbind(names(model_3UTR), model_3UTR)
model_3UTR <- as.data.frame(t(model_3UTR))

model_5UTR <- read_excel("./data/Genotype_data/HLA_G_UTR_haplotyypit_imputaatio_malli.xlsx", sheet = "5'UTR")
model_5UTR <- model_5UTR[, 2:ncol(model_5UTR)]
model_5UTR <- rbind(names(model_5UTR), model_5UTR)
model_5UTR <- as.data.frame(t(model_5UTR))

# convert 3'UTR and 5'UTR haplotypes into strings in haplotype models
haplotype_strings_3 <- unite(model_3UTR, Haplotype, -c("V1"), sep = "")
haplotype_strings_5 <- unite(model_5UTR, Haplotype, -c("V1"), sep = "")

# convert 3'UTR and 5'UTR haplotypes into strings in sample data
chosen_samples_3_haplotype_strings <- unite(chosen_samples_3, Haplotype, -c("SampleID"), sep = "")
chosen_samples_5_haplotype_strings <- unite(chosen_samples_5_corrected, Haplotype, -c("SampleID"), sep = "")

# infer 3'UTR and 5'UTR haplotypes from haplotype strings
utr_name3 <- rep('other', length(chosen_samples_3_haplotype_strings$Haplotype))
for (i in 1:length(utr_name3)) {
  hap <- chosen_samples_3_haplotype_strings$Haplotype[i]
  if (hap %in% haplotype_strings_3$Haplotype) {
    utr_name3[i] <- haplotype_strings_3$V1[haplotype_strings_3$Haplotype==hap]
  } 
}

utr_name5 <- rep('other', length(chosen_samples_5_haplotype_strings$Haplotype))
for (i in 1:length(utr_name5)) {
  hap <- chosen_samples_5_haplotype_strings$Haplotype[i]
  if (hap %in% haplotype_strings_5$Haplotype) {
    utr_name5[i] <- haplotype_strings_5$V1[haplotype_strings_5$Haplotype==hap]
  } 
}

# add haplotype names to sampleID and haplotype information and create result files
utr_haplos_3 <- cbind(chosen_samples_3_haplotype_strings, utr_name3)
utr_haplos_5 <- cbind(chosen_samples_5_haplotype_strings, utr_name5)
write.table(utr_haplos_3, sep = "\t", "./data/Genotype_data/HLA_G_3_UTR_haplotypes_final")
write.table(utr_haplos_5, sep = "\t", "./data/Genotype_data/HLA_G_5_UTR_haplotypes_final")

##################################################################################################################################
### HLA-G 3'UTR and 5'UTR haplotype frequencies
##################################################################################################################################

# count 3'UTR frequencies and create output file
counts_3UTR <- sort(table(utr_haplos_3$utr_name3), decreasing = TRUE)
freq_3UTR <- prop.table(counts_3UTR)
count_freq_3UTR <- cbind(as.data.frame(counts_3UTR), as.data.frame(freq_3UTR))
count_freq_3UTR <- count_freq_3UTR[,-3] 
colnames(count_freq_3UTR) <- c("3'UTR_Haplotype", "Count", "Freq") 
# write.table(count_freq_3UTR, sep = "\t", "./data/Output/HLA_G/HLA_G_UTR/HLA_G_3UTR_freq_Histog")

# count 5'UTR frequencies and create output file
counts_5UTR <- sort(table(utr_haplos_5$utr_name5), decreasing = TRUE)
freq_5UTR <- prop.table(counts_5UTR)
count_freq_5UTR <- cbind(as.data.frame(counts_5UTR), as.data.frame(freq_5UTR))
count_freq_5UTR <- count_freq_5UTR[,-3]
colnames(count_freq_5UTR) <- c("5'UTR_Haplotype", "Count", "Freq")
# write.table(count_freq_5UTR, sep = "\t", "./data/Output/HLA_G/HLA_G_UTR/HLA_G_5UTR_freq_Histog")

####################################################################################################################################
### Combine HLA-G alleles and 5'/3'UTR haplotypes to get extended haplotypes
####################################################################################################################################

## order alleles vertically into two columns
# extract every other line 
firstallele_5UTR <- utr_haplos_5[seq(1, nrow(utr_haplos_5), 2), ]
secondallele_5UTR <- utr_haplos_5[seq(2, nrow(utr_haplos_5), 2), ]

firstallele_3UTR <- utr_haplos_3[seq(1, nrow(utr_haplos_3), 2), ]
secondallele_3UTR <- utr_haplos_3[seq(2, nrow(utr_haplos_3), 2), ]

# bind alleles horizontally
HLA_G_haplotypes_5UTR <- cbind(firstallele_5UTR[, c(1,3)], secondallele_5UTR[,3])
HLA_G_haplotypes_3UTR <- cbind(firstallele_3UTR[, c(1,3)], secondallele_3UTR[,3])
colnames(HLA_G_haplotypes_5UTR) <- c("SampleID", "UTR5_Allele1", "UTR5_Allele2")
colnames(HLA_G_haplotypes_3UTR) <- c("SampleID", "UTR3_Allele1", "UTR3_Allele2")

# combine 5'UTR haplotype, HLA-G allele and 3'UTR haplotype data and write a table
HLA_G_ext_haplotypes_UTR <- merge(HLA_G_haplotypes_5UTR, HLA_G_pheno_Histog, by.x = "SampleID")
HLA_G_ext_haplotypes_UTR <- merge(HLA_G_ext_haplotypes_UTR, HLA_G_haplotypes_3UTR, by.x = "SampleID")
write.table(HLA_G_ext_haplotypes_UTR, sep = "\t", "./data/Phenotype_data/HLA_G_ext_haplotypes_UTR")