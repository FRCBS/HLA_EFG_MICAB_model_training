# ************************************************************************************************** #

# Dosage correlation vs allele frequency

# ************************************************************************************************** #
library(tidyverse)
library(data.table)
library(patchwork)

source("./src/FINALS/dosage_function.R")

# import imputation results
# ***
# vi
results <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb', 'compare$', full.names = T) %>% map(readRDS)
names(results) <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb', 'compare$') %>% 
  gsub('_1000G_superpop_comb_pred_compare', '', .)

# hla-g
results_g <- list.files('data/Output/Data_comb_test_results/HLA_G/', 'compare$', full.names = T) %>% map(readRDS)
names(results_g) <- list.files('data/Output/Data_comb_test_results/HLA_G/', 'compare$') %>% 
  gsub('_BB_pred_compare', '', .)

# convert to dosage
loci <- as.list(c('E', 'F', 'MICA', 'MICB')) # huom! check this!
loci_g <- as.list(c('HLA_G_3UTR', 'HLA_G_5UTR', 'HLA_G')) # huom! check this!

# convert imputation posterior probabilities to dosages
sapply(1:length(results), function(x) {
  imputation2dosage(results[[x]][[x]], loci[[x]], paste0('./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb/dosage/', loci[[x]], "_dosage"))
}) 

# import dosages 
dosages <- list.files('./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb/dosage/test2', 'raw$', full.names = T) %>% map(function(x) {
  fread(x)[, -c(1, 3:6)] 
}) 

dosages_g <- list.files('./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb/dosage/test2/HLA_G', 'raw$', full.names = T) %>% map(function(x) {
  fread(x)[, -c(1, 3:6)] 
}) 

# *****
# 1000G-FIN cv
results_fin <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/FIN_all', 'compare$', full.names = T) %>% map(readRDS)
names(results_fin) <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/FIN_all', 'compare$') %>% 
  gsub('_1000G_superpop_comb_pred_compare', '', .)

loci <- as.list(c('E', 'F', 'MICA', 'MICB')) # huom! check this!
names(results_fin) <- loci

# convert imputation posterior probabilities to dosages
sapply(1:length(results), function(x) {
  imputation2dosage(results[[x]], loci[[x]], paste0('./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/FIN_all/dosage/', loci[[x]], "_dosage"))
}) 

# import dosages 
dosages_fin <- list.files('./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/FIN_all/dosage', 'raw$', full.names = T) %>% map(function(x) {
  fread(x)[, -c(1, 3:6)] %>% as.data.frame()
}) 
# *****

# saple.id:s beginning with a number have 'X' in front of them in dosages, remove that so that id:s  match
# and remove '<present>'
for (i in 1:length(dosages)) {
  dos <- as.data.frame(dosages[[i]])
  for (row in 1:nrow(dos)) {
    a <- dos[row, 'IID']
    dos[row, 'IID'] <- ifelse(nchar(a) == 13, gsub('^.', '', a), a)
    dos <- dos %>% `colnames<-`(gsub('_<present>', '', fixed =T, colnames(.)))
    dosages[[i]] <- dos
  }
  names(dosages) <- loci
}

for (i in 1:length(dosages_g)) {
  dos <- as.data.frame(dosages_g[[i]])
  for (row in 1:nrow(dos)) {
    a <- dos[row, 'IID']
    dos[row, 'IID'] <- ifelse(nchar(a) == 13, gsub('^.', '', a), a)
    dos <- dos %>% `colnames<-`(gsub('_<present>', '', fixed =T, colnames(.)))
    dosages_g[[i]] <- dos
  }
  names(dosages_g) <- loci_g
}

for (i in 1:length(dosages_fin)) {
  dos <- as.data.frame(dosages_fin[[i]])
  for (row in 1:nrow(dos)) {
    a <- dos[row, 'IID']
    dos[row, 'IID'] <- ifelse(nchar(a) == 13, gsub('^.', '', a), a)
    dos <- dos %>% `colnames<-`(gsub('_<present>', '', fixed =T, colnames(.)))
    dosages_fin[[i]] <- dos
  }
  names(dosages_fin) <- loci
}

# extract dosages by locus
dos_e <- dosages[['E']]
dos_f <- dosages[['F']]
dos_mica <- dosages[['MICA']]
dos_micb <- dosages[['MICB']]

dos_g_3UTR <- dosages_g[['HLA_G_3UTR']]
dos_g_5UTR <- dosages_g[['HLA_G_5UTR']]
dos_g <- dosages_g[['HLA_G']]

dos_e_fin <- dosages_fin[['E']]
dos_f_fin <- dosages_fin[['F']]
dos_mica_fin <- dosages_fin[['MICA']]
dos_micb_fin <- dosages_fin[['MICB']]

# function to extract and process true alleles
true_function <- function(df, locus) {
  res <- df[[locus]]$vi$individual %>% select(c(sample.id, true.hla))
  res$true.hla <- str_replace_all(res$true.hla, '008:01/04', '008:01or04')
  res <- res %>% separate_wider_delim(true.hla, delim = "/", names = c("allele1", "allele2")) %>% as.data.frame()
  return(res)
}

true_function_fin <- function(df, locus) {
  res <- df[[locus]]$individual %>% select(c(sample.id, true.hla))
  res$true.hla <- str_replace_all(res$true.hla, '008:01/04', '008:01or04')
  res <- res %>% separate_wider_delim(true.hla, delim = "/", names = c("allele1", "allele2")) %>% as.data.frame()
  return(res)
}

true_function_g <- function(df, locus) {
  res <- df[[locus]]$individual %>% select(c(sample.id, true.hla))
  res$true.hla <- str_replace_all(res$true.hla, '01:01:01:01/02', '01:01:01:01or02')
  res <- res %>% separate_wider_delim(true.hla, delim = "/", names = c("allele1", "allele2")) %>% as.data.frame()
  return(res)
}

# apply function
res_mica <- true_function(results, 'MICA')
res_micb <- true_function(results, 'MICB')
res_e <- true_function(results, 'HLA_E')
res_f <- true_function(results, 'HLA_F')

res_g <- true_function_g(results_g, 'HLA_G')
res_g_3UTR <- true_function_g(results_g, 'HLA_G_3UTR')
res_g_5UTR <- true_function_g(results_g, 'HLA_G_5UTR')

res_mica_fin <- true_function_fin(results_fin, 'MICA')
res_micb_fin <- true_function_fin(results_fin, 'MICB')
res_e_fin <- true_function_fin(results_fin, 'E')
res_f_fin <- true_function_fin(results_fin, 'F')


# add true allele info to dosage tables
dos_res_mica <- inner_join(dos_mica, res_mica, by= c('IID' = 'sample.id')) %>% 
  mutate(allele1 = paste0('MICA*', allele1)) %>% mutate(allele2 = paste0('MICA*', allele2))
dos_res_micb <- inner_join(dos_micb, res_micb, by= c('IID' = 'sample.id')) %>% 
  mutate(allele1 = paste0('MICB*', allele1)) %>% mutate(allele2 = paste0('MICB*', allele2))
dos_res_e <- inner_join(dos_e, res_e, by= c('IID' = 'sample.id')) %>% 
  mutate(allele1 = paste0('E*', allele1)) %>% mutate(allele2 = paste0('E*', allele2))
dos_res_f <- inner_join(dos_f, res_f, by= c('IID' = 'sample.id')) %>% 
  mutate(allele1 = paste0('F*', allele1)) %>% mutate(allele2 = paste0('F*', allele2))

dos_res_g <- inner_join(dos_g, res_g, by= c('IID' = 'sample.id')) %>% 
  mutate(allele1 = paste0('G*', allele1)) %>% mutate(allele2 = paste0('G*', allele2))
dos_res_g3UTR <- inner_join(dos_g_3UTR, res_g_3UTR, by= c('IID' = 'sample.id')) %>% 
  mutate(allele1 = paste0('G*', allele1)) %>% mutate(allele2 = paste0('G*', allele2))
dos_res_g5UTR <- inner_join(dos_g_5UTR, res_g_5UTR, by= c('IID' = 'sample.id')) %>% 
  mutate(allele1 = paste0('G*', allele1)) %>% mutate(allele2 = paste0('G*', allele2))

dos_res_mica_fin <- inner_join(dos_mica_fin, res_mica_fin, by= c('IID' = 'sample.id')) %>% 
  mutate(allele1 = paste0('MICA*', allele1)) %>% mutate(allele2 = paste0('MICA*', allele2))
dos_res_micb_fin <- inner_join(dos_micb_fin, res_micb_fin, by= c('IID' = 'sample.id')) %>% 
  mutate(allele1 = paste0('MICB*', allele1)) %>% mutate(allele2 = paste0('MICB*', allele2))
dos_res_e_fin <- inner_join(dos_e_fin, res_e_fin, by= c('IID' = 'sample.id')) %>% 
  mutate(allele1 = paste0('E*', allele1)) %>% mutate(allele2 = paste0('E*', allele2))
dos_res_f_fin <- inner_join(dos_f_fin, res_f_fin, by= c('IID' = 'sample.id')) %>% 
  mutate(allele1 = paste0('F*', allele1)) %>% mutate(allele2 = paste0('F*', allele2))

# function for true dosages
dosage_calc <- function(df) {
  df <- df[, c(1, ncol(df)-1, ncol(df), 2:(ncol(df)-2))]
  for (row in 1:nrow(df)) {
    true <- unlist(c(df[row, 'allele1'], df[row, 'allele2']), use.names = F) #%>% str_replace('G_3UTR', 'G*UTR') %>% str_replace('G_5UTR_', 'G*')
    for (col in colnames(df[, 4:ncol(df)])) {
      ntrue <- sum(col %in% true[1], col %in% true[2])
      df[row, col] <- ifelse((col %in% true & ntrue == 2), 2, 
             ifelse((col %in% true & ntrue == 1), 1, 0))
    }
  }
  return(df[, c(1, 4:ncol(df))])
}

# apply function
true_dos_mica <- dosage_calc(dos_res_mica)
true_dos_micb <- dosage_calc(dos_res_micb)
true_dos_e <- dosage_calc(dos_res_e)
true_dos_f <- dosage_calc(dos_res_f)

true_dos_g <- dosage_calc(dos_res_g)
true_dos_g3UTR <- dosage_calc(dos_res_g3UTR)
true_dos_g5UTR <- dosage_calc(dos_res_g5UTR)

true_dos_mica_fin <- dosage_calc(dos_res_mica_fin)
true_dos_micb_fin <- dosage_calc(dos_res_micb_fin)
true_dos_e_fin <- dosage_calc(dos_res_e_fin)
true_dos_f_fin <- dosage_calc(dos_res_f_fin)

# function to calculate imp dosage vs true dosage correlation 
cor_function <- function(df_imp, df_true) {
  # df_imp = imputed allele dosages
  # df_true = true allele dosages
  df_cor <- data.frame()
  for (col in colnames(df_imp[-1])) {
     df_cor[1, col] <- cor(df_imp[, col], df_true[, col])
  }
  return(df_cor)
}

# apply function
cor_mica <- cor_function(dos_mica, true_dos_mica) %>% t() %>% as.data.frame() %>% 
  rownames_to_column() %>% `colnames<-`(c("allele", "cor"))
cor_micb <- cor_function(dos_micb, true_dos_micb) %>% t() %>% as.data.frame() %>% 
  rownames_to_column() %>% `colnames<-`(c("allele", "cor"))
cor_e <- cor_function(dos_e, true_dos_e) %>% t() %>% as.data.frame() %>% 
  rownames_to_column() %>% `colnames<-`(c("allele", "cor"))
cor_f <- cor_function(dos_f, true_dos_f) %>% t() %>% as.data.frame() %>% 
  rownames_to_column() %>% `colnames<-`(c("allele", "cor"))

cor_g <- cor_function(dos_g, true_dos_g) %>% t() %>% as.data.frame() %>% 
  rownames_to_column() %>% `colnames<-`(c("allele", "cor"))
cor_g3UTR <- cor_function(dos_g_3UTR, true_dos_g3UTR) %>% t() %>% as.data.frame() %>% 
  rownames_to_column() %>% `colnames<-`(c("allele", "cor")) 
cor_g5UTR <- cor_function(dos_g_5UTR, true_dos_g5UTR) %>% t() %>% as.data.frame() %>%
  rownames_to_column() %>% `colnames<-`(c("allele", "cor")) 

cor_mica_fin <- cor_function(dos_mica_fin[dos_mica_fin$IID %in% true_dos_mica_fin$IID, ], true_dos_mica_fin) %>% t() %>% as.data.frame() %>% 
  rownames_to_column() %>% `colnames<-`(c("allele", "cor"))
cor_micb_fin <- cor_function(dos_micb_fin[dos_micb_fin$IID %in% true_dos_micb_fin$IID, ], true_dos_micb_fin) %>% t() %>% as.data.frame() %>% 
  rownames_to_column() %>% `colnames<-`(c("allele", "cor"))
cor_e_fin <- cor_function(dos_e_fin, true_dos_e_fin) %>% t() %>% as.data.frame() %>% 
  rownames_to_column() %>% `colnames<-`(c("allele", "cor"))
cor_f_fin <- cor_function(dos_f_fin, true_dos_f_fin) %>% t() %>% as.data.frame() %>% 
  rownames_to_column() %>% `colnames<-`(c("allele", "cor"))


# function for allele frequencies
freq_function <- function(df_res, df_cor, locus, name) {
  # locus = name in results
  # gene_name = gene name to be removed from df_cor$allele
  freq_table <- df_res[[locus]]$vi$detail
  for (row in 1:nrow(df_cor)) {
    allele <- df_cor[row, 'allele'] 
    allele <- gsub(paste0(name, "*"), "", fixed=T, allele)
    allele_freq <- freq_table[freq_table$allele == allele, ]
    df_cor[row, 'freq'] <- allele_freq$train.freq
  }
  return(df_cor)
}

freq_function_g <- function(df_res, df_cor, locus, name) {
  # locus = name in results
  # gene_name = gene name to be removed from df_cor$allele
  freq_table <- df_res[[locus]]$detail
  for (row in 1:nrow(df_cor)) {
    allele <- df_cor[row, 'allele'] 
    allele <- gsub(paste0(name, "*"), "", fixed=T, allele)
    allele_freq <- freq_table[freq_table$allele == allele, ]
    df_cor[row, 'freq'] <- allele_freq$train.freq
  }
  return(df_cor)
}

# apply function
cor_mica <- freq_function(results, cor_mica, 'MICA', 'MICA') %>% mutate(locus = 'MICA')
cor_micb <- freq_function(results, cor_micb, 'MICB', 'MICB') %>% mutate(locus = 'MICB')
cor_e <- freq_function(results, cor_e, 'HLA_E', 'E') %>% mutate(locus = 'E')
cor_f <- freq_function(results, cor_f, 'HLA_F', 'F') %>% mutate(locus = 'F')

cor_g <- freq_function_g(results_g, cor_g, 'HLA_G', 'G') %>% mutate(locus = 'G')
cor_g3UTR <- freq_function_g(results_g, cor_g3UTR, 'HLA_G_3UTR', 'G') %>% mutate(locus = 'G_3UTR') 
cor_g3UTR$allele <- str_replace_all(cor_g3UTR$allele, fixed('G*UTR'), 'G_3UTR_') %>% str_replace_all(fixed('G*'), 'G_3UTR_')
cor_g5UTR <- freq_function_g(results_g, cor_g5UTR, 'HLA_G_5UTR', 'G') %>% mutate(locus = 'G_5UTR')
cor_g5UTR$allele <- str_replace_all(cor_g5UTR$allele, fixed('G*'), 'G_5UTR_')

cor_mica_fin <- freq_function_g(results_fin, cor_mica_fin, 'MICA', 'MICA') %>% mutate(locus = 'MICA')
cor_micb_fin <- freq_function_g(results_fin, cor_micb_fin, 'MICB', 'MICB') %>% mutate(locus = 'MICB')
cor_e_fin <- freq_function_g(results_fin, cor_e_fin, 'E', 'E') %>% mutate(locus = 'E')
cor_f_fin <- freq_function_g(results_fin, cor_f_fin, 'F', 'F') %>% mutate(locus = 'F')

# combine all
cor_all <- do.call("rbind", list(cor_mica, cor_micb, cor_e, cor_f))
cor_all_g <- do.call("rbind", list(cor_g, cor_g3UTR, cor_g5UTR))
cor_all_fin <- do.call("rbind", list(cor_mica_fin, cor_micb_fin, cor_e_fin, cor_f_fin))

# plot correlation vs allele frequency
cor_vs_freq <- ggplot(cor_all, aes(x=freq, y=cor, color=locus, label=allele)) +
  geom_point(size=4) +
  ylab('Dosage r') +
  xlab('Allele frequency') +
  theme_classic() +
  guides(color=guide_legend(title="Locus")) +
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18)) +
  theme(axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=18)) +
  theme(legend.position = c(0.9, 0.2)) +
  geom_text(aes(label=ifelse(cor<0.75, as.character(allele), '')), hjust=-0.2, vjust=0.4, size = 7, fontface='bold') 
  

cor_vs_freq_g <- ggplot(cor_all_g, aes(x=freq, y=cor, color=locus, label = allele)) +
  geom_point(size=4) +
  ylab('Dosage r') +
  xlab('Allele frequency') +
  theme_classic() +
  guides(color=guide_legend(title="Locus")) +
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18)) +
  theme(axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=18)) +
  theme(legend.position = c(0.9, 0.2)) +
  geom_text(aes(label=ifelse(cor<0.9, as.character(allele), '')), hjust=-0.2, vjust=0.4, size = 7, fontface= 'bold') 
 

cor_vs_freq_fin <- ggplot(cor_all_fin, aes(x=freq, y=cor, color=locus, label=allele)) +
  geom_point(size=4) +
  ylab('Dosage r') +
  xlab('Allele frequency') +
  theme_classic() +
  guides(color=guide_legend(title="Locus")) +
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18)) +
  theme(axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=18)) +
  theme(legend.position = c(0.9, 0.2)) +
  geom_text(aes(label=ifelse(cor<0.75, as.character(allele), '')), hjust=-0.2, vjust=0.4, size = 7, fontface= 'bold') 

dos_cor_plot <- ((cor_vs_freq + cor_vs_freq_g + cor_vs_freq_fin + 
  plot_layout(axis_titles = "collect") +
  plot_annotation(tag_levels = list(c('A', 'B', 'C'))) & 
  theme(plot.tag = element_text(size = 22))))
  
# save 
ggsave(dos_cor_plot, file="results/Plots/Manuscript/Supplementary/SF_5_dosage_correlation.jpeg", width=18, height=6, dpi = 600)

