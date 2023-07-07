### ------------------------------------------------------------------------------------------------------------------------------------------------ ###

### PLOTTING RESULTS

### ------------------------------------------------------------------------------------------------------------------------------------------------ ###

# libraries and functions

library(tidyverse)
library(readODS)
library(data.table)
library(readxl)
library(gridExtra) 
library(cowplot) 
library(ggpubr)
library(flextable) 
library(patchwork)

source('./src/functions_for_plotting.R')

### ------------------------------------------------------------------------------------------------------------------------------------------------ ###

###        FIGURE 2 (overall accuracy heatmap, limit = NULL)

### ------------------------------------------------------------------------------------------------------------------------------------------------ ###

# Read in data

dat.null <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/NULL', 'null$', full.names = T) %>% map(readRDS)
names(dat.null) <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/NULL', 'null$') %>% 
  gsub('_all_models_1000G_BB_pred_compare_null', '', .)

gene.names       <- names(dat.null)
model.names      <- dat.null[['HLA_E']] %>% names
pop.names        <- dat.null[['HLA_E']][['i']] %>% names

# get overall accuracies $overall for all genes and models
dat.overalls <- list()
for (gene in gene.names) {
  for (model in model.names) {
    for (pop in pop.names) {
      overall.temp <-dat.null[[gene]][[model]][[pop]]$overall %>% mutate(Gene=gene) %>% mutate(Model=model) %>% mutate(Pop=pop)
      dat.overalls[[length(dat.overalls)+1]] <- overall.temp
    }
  }
}

# rbind into table and melt for plotting
overall.list <-do.call("rbind", dat.overalls)
overall.list.melt <- overall.list[,c(5, 9:11)] %>% melt(id=c("Model", "Gene", "Pop")) %>% 
  mutate(Model = recode(Model, 'i'='I','ii'='II','iii'='III','iv'='IV','v'='V','vi'='VI','vii'='VII'))

# plot accuracy heatmap
accuracy_plot_all_null <- accuracy_plot_facet(overall.list.melt)

# save
ggsave(accuracy_plot_all_null, file='./results/Plots/Figure_2.jpeg', height=7.5, width = 20, dpi=600)


### ------------------------------------------------------------------------------------------------------------------------------------------------ ###

###         FIGURE 3 (Confusion matrices for model VI (MICA, MICB, HLA-E, HLA-F) and model I (HLA-G, 3'UTR, 5'UTR))

### ------------------------------------------------------------------------------------------------------------------------------------------------ ###

# read in data (results for model VI in 1000G + FIN test set)
dat.1000G.FIN.all <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb/NULL', 'null$', 
                                full.names=T) %>% map(readRDS)
dat.hla.G <- list.files('data/Output/Data_comb_test_results/HLA_G/NULL', 'null$', full.names=T) %>% 
  map(readRDS)

# name list elements
names(dat.1000G.FIN.all) <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/1000G_all_superpop_comb/NULL', 
                                       'null$') %>% gsub('_1000G_superpop_comb_pred_compare_null', '', .)

names(dat.hla.G) <- list.files('data/Output/Data_comb_test_results/HLA_G/NULL', 'null$') %>% 
  gsub('_BB_pred_compare_null', '', .)

gene.names  <- c(names(dat.1000G.FIN.all))
model.names <- "vi"
names.hla.G <- names(dat.hla.G)


### Plotting

# get confusion matrices ($confusion) for model VI for MICA, MICB, HLA-E and HLA-F 
# y = title string in the plot
dat.confusion <- map2(gene.names, c('HLA-E', 'HLA-F', 'MICA', 'MICB'), function(gene, y) {
  dat.1000G.FIN.all[[gene]][['vi']]$confusion %>% as.data.frame() %>% 
    confusion_processing %>% conf_table_plot_3(title=y)
})
names(dat.confusion) <- gene.names
dat.confusion <- dat.confusion[c('MICA', 'MICB', 'HLA_E', 'HLA_F')] # change order 

# plot hla-e and hla-f in same plot

# HLA-E and HLA-F get $confusion and add gene names for combining
dat.hla.E <- dat.1000G.FIN.all[['HLA_E']][['vi']]$confusion %>% as.data.frame() %>% confusion_processing() %>%
  mutate(name = gsub("^*", "E*", name)) %>% mutate(rowname = gsub("^*", "E*", rowname))

dat.hla.F <- dat.1000G.FIN.all[['HLA_F']][['vi']]$confusion %>% as.data.frame() %>% confusion_processing() %>% 
  mutate(name = gsub("^*", "F*", name)) %>% mutate(rowname = gsub("^*", "F*", rowname))

# combine & plot
dat.hla.EF <- rbind(dat.hla.E, dat.hla.F)

dat.confusion.hlaef <- dat.hla.EF %>% conf_table_plot_3(title="HLA-E and HLA-F")

# HLA-G, 3'UTR and 5'UTR
dat.confusion.G <- map2(1:length(dat.hla.G), 
                        c("HLA-G 3'UTR", "HLA-G 5'UTR", 'HLA-G'), function(x, y) {
  dat.hla.G[[x]]$confusion %>% as.data.frame() %>%
    confusion_processing  %>%
    conf_table_plot_3(title=y) %>% return
})

# list element names
names(dat.confusion.G) <- names.hla.G

# reorder
dat.confusion.G <- dat.confusion.G[c('HLA_G', 'HLA_G_3UTR', 'HLA_G_5UTR')] 


# plot composition

panel.1 <- dat.confusion[['MICA']] & theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

panel.2 <- (dat.confusion[['MICB']] / dat.confusion.hlaef) + plot_layout(heights=c(2.5,2)) & theme(plot.margin = unit(c(0, 1, 0, 1), "cm")) 

panel.3 <- dat.confusion.G[[1]] + dat.confusion.G[[2]] + dat.confusion.G[[3]] + 
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

confusion_plot <- ((panel.1 + panel.2 + plot_layout(widths = c(2.3, 1))) / panel.3) +
  plot_layout(heights = c(2.3, 1)) + 
  plot_annotation(tag_levels = list(c('A', '', '', 'B'))) & 
  theme(plot.tag = element_text(size = 40)) 

# save 
ggsave(confusion_plot, file="results/Plots/Figure_3.jpeg", width=18, height=20, dpi = 600)


### ------------------------------------------------------------------------------------------------------------------------------------------------ ###

###         FIGURE 4 (allelic accuracy, sensitivity and specificity vs. allele frequency)

### ------------------------------------------------------------------------------------------------------------------------------------------------ ###

# read in data
dat.nonnull <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G', 'compare$', full.names = T) %>% 
  map(readRDS)
names(dat.nonnull) <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G', 'compare$') %>% 
  gsub('_all_models_1000G_BB_pred_compare', '', .)

gene.names       <- names(dat.nonnull)
model.names.pick <- c('ii', 'v', 'vi') # model 'ii' = FIN reference, model 'v' = 1000G reference, model 'vi' = combined reference
# otan tässä vain nämä kolme mallia edustamaan referenssejä
pop.names        <- dat.nonnull[['HLA_E']][['i']] %>% names

# get allelic results ($detail) for all genes and models
dat.temp <- list()
for (gene in gene.names) {
  for (model in model.names.pick) {
    for (pop in pop.names) {
      dat.detail <-dat.nonnull[[gene]][[model]][[pop]]$detail %>% 
      mutate(Pop=pop) %>% mutate(Gene=gene) %>% mutate(Ref= ifelse(model=='ii', 'FIN', 
                                                                   ifelse(model=='v', '1000G', 'combined')))
      dat.temp[[length(dat.temp)+1]] <- dat.detail
    }
  }
}

# rbind into table and melt for plotting
detail.list <-do.call("rbind", dat.temp)

# plot accuracy, sensitivity, specificity vs. allele frequency
allele_accuracy_all <- plot_allele_prop(detail.list, accuracy, 'Accuracy' ,'')
allele_sens_all     <- plot_allele_prop(detail.list, sensitivity, 'Sensitivity', 'Allele frequency')
allele_spec_all     <- plot_allele_prop(detail.list, specificity, 'Specificity', '')


allele_freq_all <- (allele_accuracy_all + plot_spacer() +
                      allele_sens_all +  plot_spacer() +
                      allele_spec_all) + 
  plot_layout(guides = "collect", widths = c(1, .15, 1, .15, 1)) & 
  theme(legend.position = "top", 
        legend.box = 'horizontal', legend.box.just = 'left', 
        legend.margin = margin(1, 35, 20, 35)) 

ggsave(allele_freq_all, file='./results/Plots/Manuscript/Figure_4.jpeg', 
       height = 14, width = 18)


### ------------------------------------------------------------------------------------------------------------------------------------------------ ###

###         FIGURE 5 (Confusion matrices for 1000G - FIN crossvalidation)

### ------------------------------------------------------------------------------------------------------------------------------------------------ ###

# plot confusion matrices for 1000G reference - FIN reference cross-validation (all FIN imputed with 1000G model (VII))

# data
dat.FIN.all <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/FIN_all', 'compare$', 
                          full.names=T) %>% map(readRDS)

# name list elements
names(dat.FIN.all) <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/FIN_all', 'compare$') %>% 
  gsub('_BB_all_pred_compare', '', .) %>% gsub('_HSCT_all_pred_compare', '', .)

gene.names  <- names(dat.FIN.all)
accuracies <- c("0.998", "1.000", "0.996", "0.998")

# extract confusion matrices, process and plot
dat.confusion.FIN <- map2(gene.names, c('HLA-E', 'HLA-F', 'MICA', 'MICB'), function(gene, y) {
  dat.FIN.all[[gene]]$confusion %>% as.data.frame() %>% 
    confusion_processing %>% conf_table_plot_3(title=y) 
})

dat.confusion.FIN <- map2(gene.names, c('HLA-E', 'HLA-F', 'MICA', 'MICB'), function(gene, y) {
  dat.FIN.all[[gene]]$confusion %>% as.data.frame() %>% 
    confusion_processing %>% conf_table_plot_3(title=y)
})

# HLA-E and HLA-F combined 

# HLA-E and HLA-F get $confusion and add gene names for combining
dat.hla.E.FIN <- dat.FIN.all[['HLA_E']]$confusion %>% as.data.frame() %>% confusion_processing() %>%
  mutate(name = gsub("^*", "E*", name)) %>% mutate(rowname = gsub("^*", "E*", rowname))

dat.hla.F.FIN <- dat.FIN.all[['HLA_F']]$confusion %>% as.data.frame() %>% confusion_processing() %>% 
  mutate(name = gsub("^*", "F*", name)) %>% mutate(rowname = gsub("^*", "F*", rowname))

# combine
dat.hla.EF.FIN <- rbind(dat.hla.E.FIN, dat.hla.F.FIN)

# plot hla-e and hla-f in same plot
dat.confusion.hlaef.FIN <- dat.hla.EF.FIN %>% conf_table_plot_3(title="HLA-E and HLA-F")

panel.1.FIN <- dat.confusion.FIN[[3]] & theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

panel.2.FIN <- dat.confusion.FIN[[4]] / dat.confusion.hlaef.FIN #& theme(plot.margin = unit(c(0, 1, 0, 1), "cm"))

conf_plot_FIN <- (panel.1.FIN + panel.2.FIN + plot_layout(widths=c(2,1)))

ggsave(conf_plot_FIN, file = './results/Plots/Plot_5.jpeg', height = 14, width =20, dpi=600)



# **************************************************************************************************************************************************** #

#         SUPPLEMENTARY FIGURE 1 (flanking region testings)

# **************************************************************************************************************************************************** #

# read in results and add locus and dataset info
# mica
fl_tests_FIN_mica <- read_ods('data/Output/1000G/fl_region_test_BB_1000G_summary.ods', sheet='FIN_reference_mica') %>% 
    setnames(old= c('Mean', 'Acc_test'), new= c('OOB', 'Test')) %>% mutate(locus= rep("MICA", nrow(fl_tests_BB_mica))) %>% 
    mutate(dataset=rep('FIN', nrow(fl_tests_BB_mica)))
fl_tests_1000G_mica <- read_ods('data/Output/1000G/fl_region_test_BB_1000G_summary.ods', sheet='1000G_reference_mica') %>% 
    setnames(old= c('Mean', 'Acc_test'), new= c('OOB', 'Test')) %>% mutate(locus= rep("MICA", nrow(fl_tests_1000G_mica))) %>% 
    mutate(dataset=rep('1000G', nrow(fl_tests_1000G_mica)))
fl_tests_1000G_FIN_mica <- read_ods('data/Output/1000G/fl_region_test_BB_1000G_summary.ods', sheet='FIN_1000G_intersect_ref_mica') %>% 
    setnames(old= c('Mean', 'Acc_test'), new= c('OOB', 'Test')) %>% mutate(locus= rep("MICA", nrow(fl_tests_1000G_BB_mica))) %>% 
    mutate(dataset= rep('Combined', nrow(fl_tests_1000G_BB_mica)))

# micb
fl_tests_FIN_micb <- read_ods('data/Output/1000G/fl_region_test_BB_1000G_summary.ods', sheet='FIN_reference_micb') %>% 
    setnames(old= c('Mean', 'Acc_test'), new= c('OOB', 'Test')) %>% mutate(locus= rep("MICB", nrow(fl_tests_BB_micb))) %>% 
    mutate(dataset=rep('FIN', nrow(fl_tests_BB_micb)))
fl_tests_1000G_micb <- read_ods('data/Output/1000G/fl_region_test_BB_1000G_summary.ods', sheet='1000G_reference_micb') %>% 
    setnames(old= c('Mean', 'Acc_test'), new= c('OOB', 'Test')) %>% mutate(locus= rep("MICB", nrow(fl_tests_1000G_micb))) %>% 
    mutate(dataset=rep('1000G', nrow(fl_tests_1000G_micb)))
fl_tests_1000G_FIN_micb <- read_ods('data/Output/1000G/fl_region_test_BB_1000G_summary.ods', sheet='FIN_1000G_intersect_ref_micb') %>% 
    setnames(old= c('Mean', 'Acc_test'), new= c('OOB', 'Test')) %>% mutate(locus= rep("MICB", nrow(fl_tests_1000G_BB_micb))) %>% 
    mutate(dataset= rep('Combined', nrow(fl_tests_1000G_BB_micb)))

# hla-e
fl_tests_FIN_e <- read_ods('data/Output/1000G/fl_region_test_BB_1000G_summary.ods', sheet='FIN_reference_hla_e') %>% 
    setnames(old= c('Mean', 'Acc_test'), new= c('OOB', 'Test')) %>% mutate(locus= rep("HLA-E", nrow(fl_tests_BB_e))) %>% 
    mutate(dataset=rep('FIN', nrow(fl_tests_BB_e)))
fl_tests_1000G_e <- read_ods('data/Output/1000G/fl_region_test_BB_1000G_summary.ods', sheet='1000G_reference_hla_e') %>% 
    setnames(old= c('Mean', 'Acc_test'), new= c('OOB', 'Test')) %>% mutate(locus= rep("HLA-E", nrow(fl_tests_1000G_e))) %>% 
    mutate(dataset=rep('1000G', nrow(fl_tests_1000G_e)))
fl_tests_1000G_FIN_e <- read_ods('data/Output/1000G/fl_region_test_BB_1000G_summary.ods', sheet='FIN_1000G_intersect_ref_hla_e') %>% 
    setnames(old= c('Mean', 'Acc_test'), new= c('OOB', 'Test')) %>% mutate(locus= rep("HLA-E", nrow(fl_tests_1000G_BB_e))) %>% 
    mutate(dataset= rep('Combined', nrow(fl_tests_1000G_BB_e)))

# hla-f
fl_tests_FIN_f <- read_ods('data/Output/1000G/fl_region_test_BB_1000G_summary.ods', sheet='FIN_reference_hla_f') %>% 
    setnames(old= c('Mean', 'Acc_test'), new= c('OOB', 'Test')) %>% mutate(locus= rep("HLA-F", nrow(fl_tests_hsct_f))) %>% 
    mutate(dataset=rep('FIN', nrow(fl_tests_hsct_f)))
fl_tests_1000G_f <- read_ods('data/Output/1000G/fl_region_test_BB_1000G_summary.ods', sheet='1000G_reference_hla_f') %>% 
    setnames(old= c('Mean', 'Acc_test'), new= c('OOB', 'Test')) %>% mutate(locus= rep("HLA-F", nrow(fl_tests_1000G_f))) %>% 
    mutate(dataset=rep('1000G', nrow(fl_tests_1000G_f)))
fl_tests_1000G_FIN_f <- read_ods('data/Output/1000G/fl_region_test_BB_1000G_summary.ods', sheet='FIN_1000G_intersect_ref_hla_f') %>% 
    setnames(old= c('Mean', 'Acc_test'), new= c('OOB', 'Test')) %>% mutate(locus= rep("HLA-F", nrow(fl_tests_1000G_hsct_f))) %>% 
    mutate(dataset= rep('Combined', nrow(fl_tests_1000G_hsct_f)))

# hla-G
fl_tests_FIN_g <- read_ods('data/Output/1000G/fl_region_test_BB_1000G_summary.ods', sheet='FIN_reference_hla_G') %>% 
    setnames(old= c('Mean', 'Acc_test'), new= c('OOB', 'Test')) %>% mutate(locus = rep('HLA-G', nrow(fl_tests_BB_g))) %>% 
    mutate(dataset=rep('FIN', nrow(fl_tests_BB_g)))

# combine data
acc_FIN <- rbind(fl_tests_FIN_mica, fl_tests_FIN_micb, fl_tests_FIN_e, fl_tests_FIN_f, fl_tests_FIN_g)
acc_1000G <- rbind(fl_tests_1000G_mica, fl_tests_1000G_micb, fl_tests_1000G_e, fl_tests_1000G_f) %>% select(c(1:8,14:15))
acc_1000G_FIN <- rbind(fl_tests_1000G_FIN_mica, fl_tests_1000G_FIN_micb, fl_tests_1000G_FIN_e, fl_tests_1000G_FIN_f) %>% select(c(1:8, 16:17))

acc_all <- rbind(acc_FIN, acc_1000G, acc_1000G_FIN)

# facet wrap all loci in the same plot
fl_plot_oob_test <- ggplot(acc_all, aes(x= kb_flanking, group=locus)) +
    geom_point(aes(y=OOB), shape=2) + 
    geom_line(aes(y=OOB, colour='OOB'), linetype=1) + 
    geom_point(aes(y=Test), shape=8) + 
    geom_line(aes(y=Test, colour='Test'), linetype=1) +
    scale_color_manual(values=c("black","dark grey")) +
    facet_wrap(~factor(locus, levels= c('MICA', 'MICB', 'HLA-E', 'HLA-F', 'HLA-G')) + factor(dataset, levels=c('FIN', '1000G', 'Combined')), ncol=3, 
               labeller = label_wrap_gen(multi_line=FALSE)) +
    ylim(95, 100) +
    xlim(0, 100) +
    xlab('Flanking region kb') +
    ylab('Accuracy (%)') +
    theme_bw() +
    labs(colour="") +
    scale_x_continuous(breaks = c(0,5, 10, 15, 50)) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.background=element_rect(fill='white', colour='black'),strip.background=element_rect(fill='white', colour='white')) +
    theme(strip.text = element_text(size = 20, margin = margin())) +
    theme(axis.text.x = element_text(size = 18)) +
    theme(axis.text.y = element_text(size = 18)) +
    theme(legend.text = element_text(size=22)) +
    theme(legend.position = c(0.7, 0.05)) +
    theme(axis.title=element_text(size=24))

ggsave(fl_plot_oob_test, file='./results/Plots/Supplementary/SF_1.jpeg', height=18, width = 20, dpi=600)


# ************************************************************************************************************************************************** ###

#         SUPPLEMENTARY FIGURE 2 (model properties)

# ************************************************************************************************************************************************** ###

### plot train and all data OOB accuracies 

model_properties <- read_excel("./data/Output/1000G/tests_model_properties.xlsx") %>% mutate(model = recode(model, 'i'='I','ii'='II','iii'='III','iv'='IV','v'='V','vi'='VI','vii'='VII'))

# plot train OOB accuracies

model_properties_train <- model_properties %>% select(c(1:2,5:6)) %>% melt(id.vars = c("locus", "model")) %>% filter(variable=="oob_train")

plot_train_oob <- ggplot(model_properties_train, aes(x=factor(locus, level=c('MICA', 'MICB', 'HLA-E', 'HLA-F')), y = value, fill = model)) +
  geom_col(position = position_dodge(0.9)) +
  geom_text(aes(label = value), vjust = 0.5, hjust = -0.1, position = position_dodge(width = 0.9), size=4.5, angle = 90) +
  ylab('Train OOB \naccuracy (%)') +
  #scale_fill_manual(values=c("lightskyblue2", "lightskyblue3", alpha("orange1", 0.6), "orange2", "palegreen2", "palegreen3","palegreen4")) +
  scale_fill_brewer(palette ='Paired') +
  coord_cartesian(ylim=c(95, 100.5)) +
  #ggtitle('B') +
  guides(fill = guide_legend(nrow = 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16, color='black'),
        axis.text.y = element_text(size = 14, color='black')) +
  theme(axis.ticks.x =  element_blank()) +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0, size=22, face = 'bold')) +
  theme(legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=14),
        legend.title = element_text(size=18)) +
  theme(axis.title.y = element_text(size=16),
        axis.title.x = element_blank())

plot_train_oob

# plot all data OOB accuracies

model_properties_alldata <- model_properties %>% select(c(1:2,6:7)) %>% melt(id.vars = c("locus", "model")) %>% filter(variable=="oob_all_data")

plot_alldata_oob <- ggplot(model_properties_alldata, aes(x=factor(locus, level=c('MICA', 'MICB', 'HLA-E', 'HLA-F')), y = value, fill = model)) +
  geom_col(position = position_dodge(0.9)) +
  geom_text(aes(label = value), vjust = 0.5, hjust = -0.1, position = position_dodge(width = 0.9), size=4.5, angle = 90) +
  xlab('Locus') +
  ylab('All data OOB \naccuracy (%)') +
  scale_fill_brewer(palette = 'Paired') +
  coord_cartesian(ylim=c(95, 100.5)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16, color='black'),
        axis.text.y = element_text(size = 14, color='black')) +
  theme(axis.ticks.x =  element_blank()) +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0, size=22, face = 'bold')) +
  theme(legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=14),
        legend.title = element_text(size=18)) +
  theme(axis.title.y = element_text(size=16),
        axis.title.x = element_blank())

plot_alldata_oob

### Plot SNP and allele numbers (panels D + E)

# plot SNP numbers in models

model_properties_snps <- model_properties %>% select(c(1:2,4)) %>% melt(id.vars = c("locus", "model"))

SNP_plot <- ggplot(model_properties_snps, aes(x=factor(locus, level=c('MICA', 'MICB', 'HLA-E', 'HLA-F')), y = value, fill = model)) +
  geom_col(position = position_dodge(0.9)) +
  geom_text(aes(label = value), vjust = 0.5, hjust = -0.1, position = position_dodge(width = 0.9), size=4.5, angle = 90) +
  xlab('Locus') +
  ylab('Number of SNPs') +
  scale_fill_brewer(palette = 'Paired') +
  scale_y_continuous(limits=c(0, 1550), breaks = seq(0, 1400, by = 200), expand = c(0, 0)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16, color='black'),
        axis.text.y = element_text(size = 14, color='black')) +
  theme(axis.ticks.x =  element_blank()) +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0, size=22, face = 'bold')) +
  theme(legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=14),
        legend.title = element_text(size=18)) +
  theme(axis.title.y = element_text(size=16),
        axis.title.x = element_blank())

SNP_plot

# plot number of alleles in models

model_properties_alleles <- model_properties %>% select(c(1:3)) %>% melt(id.vars = c("locus", "model"))

Allele_plot <- ggplot(model_properties_alleles, aes(x=factor(locus, level=c('MICA', 'MICB', 'HLA-E', 'HLA-F')), y = value, fill = model)) +
  geom_col(position = position_dodge(0.9)) +
  geom_text(aes(label = value), vjust = -0.7, position = position_dodge(width = 0.9), size=4.5) +
  xlab('Locus') +
  ylab('Number of alleles') +
  scale_fill_brewer(palette = 'Paired') +
  scale_y_continuous(limits=c(0, 50), breaks = seq(0, 50, by = 5), expand = c(0, 0)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16, color='black'),
        axis.text.y = element_text(size = 14, color='black')) +
  theme(axis.ticks.x =  element_blank()) +
  theme(legend.position="bottom") +
  theme(plot.title = element_text(hjust = 0, size=22, face = 'bold')) +
  theme(legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=14),
        legend.title = element_text(size=18)) +
  theme(axis.title.y = element_text(size=16),
        axis.title.x = element_blank())

Allele_plot

model_properties_plot <- ( plot_train_oob + plot_alldata_oob + SNP_plot + Allele_plot + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')) + 
  plot_annotation(tag_levels = 'A')  & theme(plot.tag = element_text(size = 22))

ggsave(model_properties_plot, file='./results/Plots/Supplementary/SF_2.jpeg', height=10, width = 16)

# **************************************************************************************************************************************************** #

#         SUPPLEMENTARY FIGURE 3 (Overall heatmap, limit = model)

# **************************************************************************************************************************************************** #

dat.nonnull <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G', 'compare$', full.names = T) %>% map(readRDS)
names(dat.nonnull) <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G', 'compare$') %>% 
  gsub('_all_models_1000G_BB_pred_compare', '', .)

gene.names       <- names(dat.nonnull)
model.names      <- dat.nonnull[['HLA_E']] %>% names
pop.names        <- dat.nonnull[['HLA_E']][['i']] %>% names

# get overall accuracies $overall for all genes and models
dat.overalls.nonnull <- list()
for (gene in gene.names) {
  for (model in model.names) {
    for (pop in pop.names) {
      overall.temp.nonnull <-dat.nonnull[[gene]][[model]][[pop]]$overall %>% mutate(Gene=gene) %>% mutate(Model=model) %>% mutate(Pop=pop)
      dat.overalls.nonnull[[length(dat.overalls.nonnull)+1]] <- overall.temp.nonnull
    }
  }
}

# rbind into table and melt for plotting
overall.list.nonnull <-do.call("rbind", dat.overalls.nonnull)
overall.list.melt.nonnull <- overall.list.nonnull[,c(5, 9:11)] %>% melt(id=c("Model", "Gene", "Pop")) %>% 
  mutate(Model = recode(Model, 'i'='I','ii'='II','iii'='III','iv'='IV','v'='V','vi'='VI','vii'='VII'))


# plot using accuracy_plot_facet function 
accuracy_plot_all <- accuracy_plot_facet(overall.list.melt.nonnull)
accuracy_plot_all

ggsave(accuracy_plot_all, file = './results/Plots/Supplementary/SF_3.jpeg', height = 7.5, width =20, dpi=600)

# **************************************************************************************************************************************************** #

#         SUPPLEMENTARY FIGURE 4 (Overall delta)

# **************************************************************************************************************************************************** #

# overall.list.melt.nonnull from Supplementary Figure 3 & overall.list.melt from Figure 2
comb.overall.list <- cbind(overall.list.melt.nonnull, overall.list.melt$value) %>% rename(value_nonnull = "value") %>% 
  rename(value_null = "overall.list.melt$value")

# delta per pop and model

plot_overall_delta <-
  ggplot(comb.overall.list, aes(x= value_null, y=value_nonnull)) +
  geom_point(aes(shape=Model, color=Pop), size=3, show.legend = T) +
  scale_shape_manual(values=c(1,2,5,6,7,8,9)) +
  geom_abline(linetype=3) +
  xlab("Accuracy (%)/no mask") +
  ylab("Accuracy (%)/untrained alleles masked") +
  ylim(0.8,1) +
  xlim(0.8,1) +
  facet_wrap(vars(factor(Gene, levels=c("MICA", "MICB", "HLA_E", "HLA_F")))) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size=12), axis.title=element_text(size=14,)) +
  theme(strip.text.x = element_text(size = 14)) +
  theme(legend.title = element_text(size=12), legend.text = element_text(size=10)) +
  theme(panel.background=element_rect(fill='white', colour='black'),strip.background=element_rect(fill='white', colour='white'))


# ggsave(plot_overall_delta, file='./results/Plots/overall_delta_point.jpeg', height=10, width = 12)


# delta heatmap

# overall.list.melt.nonnull from Supplementary Figure 3 & overall.list.melt from Figure 2
comb.overall.list$delta <- comb.overall.list$value_nonnull - comb.overall.list$value_null

gene.names <- c("MICA", "MICB", "HLA-E", "HLA-F")
names(gene.names) <- c("MICA", "MICB", "HLA_E", "HLA_F")

accuracy_delta <- ggplot(comb.overall.list, aes(factor(Pop, levels=c('FIN','EUR','AFR','EAS','SAS','AMR')), Model, fill=delta, group=Gene)) + 
  geom_tile(color = "black", show.legend = T, linewidth = .2) + 
  guides(colour = "none") +
  geom_text(aes(label = round(delta, digits=3)), size=6) +
  scale_fill_viridis_c(alpha=0.4) + 
  scale_y_discrete(limits=rev) +
  facet_wrap(vars(factor(Gene, levels=c("MICA", "MICB", "HLA_E", "HLA_F"))), ncol=4,
             labeller = as_labeller(gene.names)) +
  labs(fill = "delta") +
  ylab(" ") +
  xlab(" ") +
  theme(panel.grid.major.x=element_line(), 
        panel.grid.minor.x=element_line(), 
        panel.grid.major.y=element_line(), 
        panel.grid.minor.y=element_line(),
        panel.background=element_rect(fill="white"), 
        axis.text.x = element_text(hjust = 0.5, vjust = 1,size = 18, face="bold"),
        plot.title = element_text(size = 26, hjust = 0, face = 'bold'),
        axis.text.y = element_text(size = 18, face="bold"),
        axis.ticks = element_blank()) + 
  theme(strip.text.x = element_text(size = 20),
        strip.background = element_blank()) +
  theme(legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20, margin = margin(0, 0, 9, 0))) +
  theme(plot.title = element_text(hjust = 0, size=26)) +
  theme(legend.position = "bottom") +
  theme(legend.key.width= unit(2, 'cm'))


ggsave(accuracy_delta, file='./results/Plots/Supplementary/SF_4.jpeg', height=7.5, width = 20, dpi=600)