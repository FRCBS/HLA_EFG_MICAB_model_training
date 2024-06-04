# ................................................................................................. #

# LOAD AND PROCESS FIN AND 1000G REFERENCE GENOTYPE AND PHENOTYPE DATA

# ................................................................................................. #

# GENOTYPING DATA

# ................................................................................................. #

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

Genotypedata_1000G <- hlaBED2Geno(bed.fn = "data/Genotype_data/plink_SNP_data/chr6_phase3_MHC_rem_dup.bed", 
                                  fam.fn = "data/Genotype_data/plink_SNP_data/chr6_phase3_MHC_rem_dup.fam", 
                                  bim.fn = "data/Genotype_data/plink_SNP_data/chr6_phase3_MHC_rem_dup.bim", assembly="hg19") # used for model VII

# hg37 -> hg38 conversion done using https://www.ncbi.nlm.nih.gov/genome/tools/remap
hg38_converted_1000G_rem_dup_snps <- read.table("./data/Genotype_data/plink_SNP_data/report_1000Gsnps_rem_dup_hg37b.txt", sep = "\t", header = TRUE)
Genotypedata_1000G$snp.position <- hg38_converted_1000G_rem_dup_snps$mapped_start
Genotypedata_1000G$assembly <- 'hg38'

# Combine Finnish I or II with 1000G genotype data by SNP position
geno_comb_1000G_I <-hlaGenoCombine(Genotypedata_1000G, geno_FIN_I, match.type="Position") # used for MICA, MICB, HLA-E
geno_comb_1000G_II <-hlaGenoCombine(Genotypedata_1000G, geno_FIN_II, match.type="Position") # used for HLA-F

# the combined genotyping data are:
# MICA, MICB, HLA-E: geno_comb_1000G_I
# HLA-F: geno_comb_1000G_II
# HLA-G: geno_FIN_I

# .................................................................................................. #

# PHENOTYPE DATA

# .................................................................................................. #

## 1000G

# read in 1000G phenotype data (all populations)
sheetlist <- excel_sheets(path="./data/Phenotype_data/1000G_results_all_processed_pop_info.xlsx")

sheets <- list()
for (i in 1:length(sheetlist)) {
  tempdf <- read_excel(path="./data/Phenotype_data/1000G_results_all_processed_pop_info.xlsx", sheet = sheetlist[i])
  sheets[[i]] <- tempdf
  # get samples that are in genotype data
  sheets[[i]] <- sheets[[i]] %>% as.data.frame() %>% filter(Sample_ID %in% Genotypedata_1000G$sample.id)
}
names(sheets) <- sheetlist
loci <- gsub('HLA-', '', sheetlist)

# split 1000G all samples into 2/3 train and 1/3 test
divide_1000G <- list()
for (i in 1:length(sheets)) { # Convert phenotype data into hlaAlleleClass object
  divide_1000G[[i]] <- hlaAllele(sheets[[i]]$Sample_ID, H1=sheets[[i]][, "allele1"], H2=sheets[[i]][, "allele2"], 
                                 locus = loci[i], assembly = "hg38")
  set.seed(100)
  divide_1000G[[i]] <- hlaSplitAllele(divide_1000G[[i]], train.prop=0.666) 
}
names(divide_1000G) <- sheetlist

# subset superpopulation validation samples 
superpop <- c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')
superpop_val <- list()
for(x in names(divide_1000G))  {
  for(y in superpop) {
    superpop_val[[x]][[y]] <- hlaAlleleSubset(divide_1000G[[x]]$validation, samp.sel=divide_1000G[[x]]$validation$value$sample.id %in% filter(sheets[[x]], Superpopulation_code == y)$Sample_ID)
  }
}

# 1000G all samples, no splitting into train and test
divide_1000G_all <- list() 
loci <- gsub('HLA-', '', sheetlist)
for (i in 1:length(sheets)) { 
  divide_1000G_all[[i]] <- hlaAllele(sheets[[i]]$Sample_ID, H1=sheets[[i]][, "allele1"], H2=sheets[[i]][, "allele2"], 
                                     locus = loci[i], assembly = "hg38")
}
names(divide_1000G_all) <- sheetlist
names(divide_1000G_all)[[3]] <- 'HLA_E'
names(divide_1000G_all)[[4]] <- 'HLA_F'

# .....

## Finnish reference (FIN_I and FIN_II)

pheno_FIN <- list.files('./data/Phenotype_data/', pattern='Histog$|F_alleles|haplotypes_UTR', full.names = T) %>% map(read.table)
names(pheno_FIN) <- c(list.files('./data/Phenotype_data/', 'Histog$|F_alleles|haplotypes_UTR') %>% 
                        gsub('_alleles_Histog|_alleles_HSCT|_alleles_merged_Histog|_ext_haplotypes', '', .))

# separate MICA, MICB, HLA-G 3'UTR and 5'UTR into own dataframes
pheno_FIN[['MICA']] <- pheno_FIN[['MICAB']] %>% select(c(SampleID, MICA_Allele1, MICA_Allele2))
pheno_FIN[['MICB']] <- pheno_FIN[['MICAB']] %>% select(c(SampleID, MICB_Allele1, MICB_Allele2))
pheno_FIN[['HLA_G_3UTR']] <- pheno_FIN[['HLA_G_UTR']] %>% select(c(SampleID, UTR3_Allele1, UTR3_Allele2)) %>% 
  rename(HLA_G_3UTR_Allele1=UTR3_Allele1) %>% rename(HLA_G_3UTR_Allele2=UTR3_Allele2)
pheno_FIN[['HLA_G_5UTR']] <- pheno_FIN[['HLA_G_UTR']] %>% select(c(SampleID, UTR5_Allele1, UTR5_Allele2)) %>%
  rename(HLA_G_5UTR_Allele1=UTR5_Allele1) %>% rename(HLA_G_5UTR_Allele2=UTR5_Allele2)
pheno_FIN[['HLA_G_UTR']] <- NULL
pheno_FIN[['MICAB']] <- NULL

divide_FIN <- list()
loci <- gsub('HLA_|_3UTR|_5UTR', '', names(pheno_FIN))

# FIN all samples as hlaAllele object
for (i in 1:length(pheno_FIN)) { # Convert phenotype data into hlaAlleleClass object
  divide_FIN[[i]] <- hlaAllele(pheno_FIN[[i]]$SampleID, H1=pheno_FIN[[i]][, paste0(names(pheno_FIN)[i],"_Allele1")], H2=pheno_FIN[[i]][, paste0(names(pheno_FIN)[i],"_Allele2")], 
                               locus = loci[i], assembly = "hg38")
}
names(divide_FIN) <- names(pheno_FIN)

# Extract one of the samples with alleles present in <3 samples  (MICA*011:01, MICA*009:02, MICA*001, MICB*005:08, MICB*005:06, MICB*005:01, HLA-G*01:01:01:09, HLA-G*01:21N)
mica_FIN_rare <- hlaAlleleSubset(divide_FIN[['MICA']], c(241,277,412))
micb_FIN_rare <- hlaAlleleSubset(divide_FIN[['MICB']], c(82,153,216))
hla_G_FIN_rare <- hlaAlleleSubset(divide_FIN[['HLA_G']], c(117,257))

# subset rest of the samples
divide_FIN[['MICA']] <- hlaAlleleSubset(divide_FIN[['MICA']], samp.sel = !(divide_FIN[['MICA']]$value$sample.id %in% mica_FIN_rare$value$sample.id))
divide_FIN[['MICB']] <- hlaAlleleSubset(divide_FIN[['MICB']], samp.sel = !(divide_FIN[['MICB']]$value$sample.id %in% micb_FIN_rare$value$sample.id))
divide_FIN[['HLA_G']] <- hlaAlleleSubset(divide_FIN[['HLA_G']], samp.sel = !(divide_FIN[['HLA_G']]$value$sample.id %in% hla_G_FIN_rare$value$sample.id))

# Split samples randomly into 2/3 training and 1/3 test
for (i in names(divide_FIN)) { # Convert phenotype data into hlaAlleleClass object
  set.seed(100)
  divide_FIN[[i]] <- hlaSplitAllele(divide_FIN[[i]], train.prop=0.666) # Split samples randomly into 2/3 training and 1/3 test
}

# add samples of rare alleles to validation sets
divide_FIN[['MICA']]$validation <- hlaCombineAllele(divide_FIN[['MICA']]$validation, mica_FIN_rare)
divide_FIN[['MICB']]$validation <- hlaCombineAllele(divide_FIN[['MICB']]$validation, micb_FIN_rare)
divide_FIN[['HLA_G']]$validation <- hlaCombineAllele(divide_FIN[['HLA_G']]$validation, hla_G_FIN_rare)
names(divide_1000G)[[3]] <- 'HLA_E'
names(divide_1000G)[[4]] <- 'HLA_F'

# all FIN samples, no splitting into train and test
divide_FIN_all <- list() # all FIN samples, no split between train and test
loci <- gsub('HLA_|_3UTR|_5UTR', '', names(pheno_FIN))
for (i in 1:length(divide_FIN)) { 
  divide_FIN_all[[i]] <- hlaAllele(pheno_FIN[[i]]$SampleID, H1=pheno_FIN[[i]][, paste0(names(pheno_FIN)[i],"_Allele1")], H2=pheno_FIN[[i]][, paste0(names(pheno_FIN)[i],"_Allele2")], 
                                   locus = loci[i], assembly = "hg38")
}
names(divide_FIN_all) <- names(pheno_FIN)

# .....

# combine 1000G all + FIN all samples
divide_comb <- list()
for (i in names(divide_1000G_all)) {
  divide_comb[[i]] <- hlaCombineAllele(divide_1000G_all[[i]], divide_FIN_all[[i]])
}
names(divide_comb) <- names(divide_1000G_all)

# combine 1000G + FIN training
divide_comb_train <- list()
for (i in names(divide_1000G)) {
  divide_comb_train[[i]] <- hlaCombineAllele(divide_1000G[[i]]$training, divide_FIN[[i]]$training)
}

# combine 1000G + FIN validation
divide_comb_val <- list()
for (i in names(divide_1000G)) {
  divide_comb_val[[i]] <- hlaCombineAllele(divide_1000G[[i]]$validation, divide_FIN[[i]]$validation)
}