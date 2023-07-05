
## ===============================================================
##
## Confusion matrices (all models and genes, supplementary file 2)
##
## ===============================================================

# libraries and functions

library(tidyverse)
library(HIBAG)

source('./src/FINALS/functions_for_plotting.R')

## ---------------------------------------------------------------
## data read-in
## ---------------------------------------------------------------

# list and read rds files
dat.nonnull <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G', 'compare$', full.names = T) %>% map(readRDS)
dat.null    <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/NULL', 'null$', full.names = T) %>% map(readRDS)
dat.hla.G.nonnull <- list.files('data/Output/Data_comb_test_results/HLA_G', 'compare$', full.names=T) %>% map(readRDS)
dat.hla.G.null <- list.files('data/Output/Data_comb_test_results/HLA_G/NULL', 'null$', full.names=T) %>% map(readRDS)

# name list elements
names(dat.nonnull) <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G', 'compare$') %>% 
  gsub('_all_models_1000G_BB_pred_compare', '', .)
names(dat.null) <- list.files('data/Output/Data_comb_test_results/BB_1000G/Combined_BB_1000G/NULL', 'null$') %>% 
  gsub('_all_models_1000G_BB_pred_compare_null', '', .)
names(dat.hla.G.nonnull) <- list.files('data/Output/Data_comb_test_results/HLA_G', 'compare$') %>% 
  gsub('_BB_pred_compare', '', .)
names(dat.hla.G.null) <- list.files('data/Output/Data_comb_test_results/HLA_G/NULL', 'null$') %>% 
  gsub('_BB_pred_compare_null', '', .)

# create name variables
gene.names       <- c(names(dat.nonnull))#, names(dat.null))
model.null.names <- c(rep('model', 4), rep('NULL', 4))
model.names      <- dat.nonnull[['HLA_E']] %>% names
pop.names        <- dat.nonnull[['HLA_E']][['i']] %>% names
gene.names.G     <- c(names(dat.hla.G.nonnull))

## ---------------------------------------------------------------
## produce plots into a single pdf file
## ---------------------------------------------------------------

pdf('results/Plots/confusions.pdf', onefile = T, paper='A4', width = 7, height = 8) 
map(gene.names, function(gene) {
  map(model.names, function(model) {
    map(pop.names, function(pop) {
      list(
      dat.nonnull[[gene]][[model]][[pop]]$confusion %>% 
        confusion_processing %>% conf_table_plot(title=paste0('gene = ', gene, 
                                                              '\nmodel = ', model, 
                                                              '\nmodel limit = model', 
                                                              '\npop = ', pop)),
      dat.null[[gene]][[model]][[pop]]$confusion %>% 
        confusion_processing %>% conf_table_plot(title=paste0('gene = ', gene, 
                                                              '\nmodel = ', model, 
                                                              '\nmodel limit = NULL', 
                                                              '\npop = ', pop))
      )
    })
  })
}) 
dev.off()

# HLA-G

pdf('results/Plots/confusions_hlaG.pdf', onefile = T, paper='A4', width = 7, height = 8) 
map(gene.names.G, function(gene) {
      list(
        dat.hla.G.nonnull[[gene]]$confusion %>% 
          confusion_processing %>% conf_table_plot(title=paste0('gene = ', gene, 
                                                                '\nmodel = ', "i", 
                                                                '\nmodel limit = model',
                                                                '\npop = ', 'FIN')),
        dat.hla.G.null[[gene]]$confusion %>% 
          confusion_processing %>% conf_table_plot(title=paste0('gene = ', gene, 
                                                                '\nmodel = ', "i", 
                                                                '\nmodel limit = NULL',
                                                                '\npop = ', 'FIN'))
      )
}) 
dev.off()





