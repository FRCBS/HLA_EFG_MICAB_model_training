# ************************************************************************************************** #

# Dosage correlation vs allele frequency (Supplementary Figure 5)

# ************************************************************************************************** #
library(tidyverse)
library(data.table)
library(patchwork)

source("./src/dosage_function.R")
source("./src/functions_for_dosage_cor.R")

# import imputation results
# ***
# HLA-E, -F, MICA and MICA (model vi)
results <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb', 'compare$', full.names = T) %>% map(readRDS)
names(results) <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb', 'compare$') %>% 
  gsub('_1000G_superpop_comb_pred_compare', '', .)

# HLA-G (model i)
results_g <- list.files('data/Output/Data_comb_test_results/HLA_G/', 'compare$', full.names = T) %>% map(readRDS)
names(results_g) <- list.files('data/Output/Data_comb_test_results/HLA_G/', 'compare$') %>% 
  gsub('_BB_pred_compare', '', .)

# convert imputation posterior probabilities to dosages
loci <- names(results) %>% gsub('HLA_', '', .)
loci_names_g <- names(results_g) %>% gsub('HLA_', '', .)
loci_g <- loci_names_g %>% gsub('_3UTR|_5UTR', '', .)

#sapply(1:length(results), function(x) {
#  imputation2dosage(results[[x]][[x]], loci[[x]], paste0('./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb/dosage/', loci[[x]], "_dosage"))
#}) 

#sapply(1:length(results_g), 1:length(loci_names_g), function(x) {
#  imputation2dosage(results_g[[x]], loci_g[[x]], paste0('./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb/dosage/', loci_names_g[[x]], "_dosage"))
#}) 

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
  gsub('_BB_all_pred_compare', '', .) %>% gsub('_HSCT_all_pred_compare', '', .)

loci <- names(results_fin) %>% gsub('HLA_', '', .) 

# convert imputation posterior probabilities to dosages
#sapply(1:length(results_fin), function(x) {
#  imputation2dosage(results_fin[[x]], loci[[x]], paste0('./data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/FIN_all/dosage/', loci[[x]], "_dosage"))
#}) 

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
  names(dosages_g) <- loci_names_g
}

for (i in 1:length(dosages_fin)) {
  dos <- as.data.frame(dosages_fin[[i]]) %>% 
    filter(IID %in% results_fin[[i]]$individual$sample.id) # only alleles common in FIN and 1000G reference
  for (row in 1:nrow(dos)) {
    a <- dos[row, 'IID']
    dos[row, 'IID'] <- ifelse(nchar(a) == 13, gsub('^.', '', a), a)
    dos <- dos %>% `colnames<-`(gsub('_<present>', '', fixed =T, colnames(.)))
    dosages_fin[[i]] <- dos
  }
  names(dosages_fin) <- loci
}

# apply true function
res_1000G <- list()
for (i in 1:length(results)) {
  res_1000G[[i]] <- true_function(results[[i]], names(results)[[i]], '1000G')
  names(res_1000G)[i] <- names(results)[[i]]
}

res_g <- list()
for (i in 1:length(results_g)) {
  res_g[[i]] <- true_function(results_g[[i]], loci_g[[i]], 'FIN')
  names(res_g)[i] <- names(results_g)[[i]]
}

res_fin <- list()
for (i in 1:length(results_fin)) {
  res_fin[[i]] <- true_function(results_fin[[i]], names(results_fin)[[i]], 'FIN')
  names(res_fin)[i] <- names(results_fin)[[i]]
}

# combine imp dosages and true alleles
dos_res_1000G <- sapply(1:length(dosages), function(x) {
  inner_join(dosages[[x]], res_1000G[[x]], by= c('IID' = 'sample.id'))
})
names(dos_res_1000G) <- names(dosages)

dos_res_g <- sapply(1:length(dosages_g), function(x) {
  inner_join(dosages_g[[x]], res_g[[x]], by= c('IID' = 'sample.id'))
})
names(dos_res_g) <- names(dosages_g)

dos_res_fin <- sapply(1:length(dosages_fin), function(x) {
  inner_join(dosages_fin[[x]], res_fin[[x]], by= c('IID' = 'sample.id'))
})
names(dos_res_fin) <- names(dosages_fin)


# apply dosage_calc function to get true allele dosages
true_dos_1000G <- list()
for (i in 1:length(dos_res_1000G)) {
  true_dos_1000G[[i]] <- dosage_calc(dos_res_1000G[[i]])
  names(true_dos_1000G)[i] <- names(dos_res_1000G)[[i]]
}

true_dos_g <- list()
for (i in 1:length(dos_res_g)) {
  true_dos_g[[i]] <- dosage_calc(dos_res_g[[i]])
  names(true_dos_g)[i] <- names(dos_res_g)[[i]]
}

true_dos_fin <- list()
for (i in 1:length(dos_res_fin)) {
  true_dos_fin[[i]] <- dosage_calc(dos_res_fin[[i]])
  names(true_dos_fin)[i] <- names(dos_res_fin)[[i]]
}


# apply cor_function to calculate imputed vs. true allele dosages
cor_1000G <- list()
for (x in 1:length(dosages)) {
  cor_1000G[[x]] <- cor_function(dosages[[x]], true_dos_1000G[[x]]) %>% t() %>% as.data.frame() %>% 
    rownames_to_column() %>% `colnames<-`(c("allele", "cor"))
}
names(cor_1000G) <- names(dosages)

cor_g <- list()
for (x in 1:length(dosages_g)) {
  cor_g[[x]] <- cor_function(dosages_g[[x]], true_dos_g[[x]]) %>% t() %>% as.data.frame() %>% 
    rownames_to_column() %>% `colnames<-`(c("allele", "cor"))
}
names(cor_g) <- names(dosages_g)

cor_fin <- list()
for (x in 1:length(dosages_fin)) {
  cor_fin[[x]] <- cor_function(dosages_fin[[x]], true_dos_fin[[x]]) %>% t() %>% as.data.frame() %>% 
    rownames_to_column() %>% `colnames<-`(c("allele", "cor"))
}
names(cor_fin) <- names(dosages_fin)



# add allele frequencies 
for (i in 1:length(cor_1000G)) {
  cor_1000G[[i]] <- freq_function(results, cor_1000G[[i]], names(results)[[i]], loci[[i]], '1000G') %>% mutate(locus= loci[[i]])
}

for (i in 1:length(cor_g)) {
  cor_g[[i]] <- freq_function(results_g, cor_g[[i]], names(results_g)[[i]], loci_g[[i]], 'FIN') %>% mutate(locus= sub('HLA_', '', names(results_g))[[i]])
}

cor_g[['G_3UTR']]$allele <- str_replace_all(cor_g[['G_3UTR']]$allele, fixed('G*UTR'), 'G_3UTR') %>% str_replace_all(fixed('G*'), 'G_3UTR_')
cor_g[['G_5UTR']]$allele <- str_replace_all(cor_g[['G_5UTR']]$allele, fixed('G*'), 'G_5UTR_')

for (i in 1:length(cor_fin)) {
  cor_fin[[i]] <- freq_function(results_fin, cor_fin[[i]], names(results_fin)[[i]], loci[[i]], 'FIN') %>% mutate(locus= loci[[i]])
}

# combine all
cor_all_1000G <- do.call(rbind, cor_1000G)
cor_all_g <- do.call(rbind, cor_g)
cor_all_fin <- do.call(rbind, cor_fin)

# plot correlation vs allele frequency
cor_vs_freq <- ggplot(cor_all_1000G, aes(x=freq, y=cor, color=locus, label=allele)) +
  geom_point(size=4) +
  ylab('Dosage r') +
  xlab('Allele frequency') +
  theme_classic() +
  guides(color=guide_legend(title="Locus")) +
  theme(legend.title=element_text(size=16), legend.text=element_text(size=14)) +
  theme(axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=18)) +
  theme(legend.position = c(0.9, 0.2)) +
  geom_text(aes(label=ifelse(cor<0.75, as.character(allele), '')), hjust=-0.2, vjust=0.4, size = 7, fontface='bold') 
  

cor_vs_freq_g <- ggplot(cor_all_g, aes(x=freq, y=cor, color=locus, label = allele)) +
  geom_point(size=4) +
  ylab('Dosage r') +
  xlab('Allele frequency') +
  theme_classic() +
  guides(color=guide_legend(title="Locus")) +
  theme(legend.title=element_text(size=16), legend.text=element_text(size=14)) +
  theme(axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=18)) +
  theme(legend.position = c(0.9, 0.2)) +
  geom_text(aes(label=ifelse(cor<0.9, as.character(allele), '')), hjust=-0.2, vjust=0.4, size = 7, fontface= 'bold') 
 

cor_vs_freq_fin <- ggplot(cor_all_fin, aes(x=freq, y=cor, color=locus, label=allele)) +
  geom_point(size=4) +
  ylab('Dosage r') +
  xlab('Allele frequency') +
  theme_classic() +
  guides(color=guide_legend(title="Locus")) +
  theme(legend.title=element_text(size=16), legend.text=element_text(size=14)) +
  theme(axis.title=element_text(size=18,face="bold"), axis.text=element_text(size=18)) +
  theme(legend.position = c(0.9, 0.2)) +
  geom_text(aes(label=ifelse(cor<0.75, as.character(allele), '')), hjust=-0.2, vjust=0.4, size = 7, fontface= 'bold') 

dos_cor_plot <- ((cor_vs_freq + cor_vs_freq_g + cor_vs_freq_fin + 
  plot_layout(axis_titles = "collect") +
  plot_annotation(tag_levels = list(c('A', 'B', 'C'))) & 
  theme(plot.tag = element_text(size = 22))))
  
# save 
ggsave(dos_cor_plot, file="results/Plots/Manuscript/Supplementary/SF_5_dosage_correlation.jpeg", width=18, height=6, dpi = 600)

