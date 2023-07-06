# function for plotting accuracy heatmap

accuracy_plot_facet <- function(i) {
  # i = data.frame
  i$value <- i$value*100
  
  gene.names <- c("MICA", "MICB", "HLA-E", "HLA-F")
  names(gene.names) <- c("MICA", "MICB", "HLA_E", "HLA_F")
  
  a <- ggplot(i, aes(factor(Pop, levels=c('FIN','EUR','AFR','EAS','SAS','AMR')), Model, fill=value, group=Gene)) + 
    geom_tile(color = "black", show.legend = T, linewidth = .2) + 
    guides(colour = "none") +
    geom_text(aes(label = round(value, digits=1)), size=6.5) +
    scale_fill_viridis_c(alpha=0.8, begin=0.2, label = function(x) sprintf("%.0f", x), limits=c(80,100), breaks=c(80,85,90,95,100)) + 
    scale_y_discrete(limits=rev) +
    facet_wrap(vars(factor(Gene, levels=c("MICA", "MICB", "HLA_E", "HLA_F"))), ncol=4,
               labeller = as_labeller(gene.names)) +
    labs(fill = "Accuracy (%)") +
    ylab(" ") +
    xlab(" ") +
    theme(panel.grid.major.x=element_line(), 
          panel.grid.minor.x=element_line(), 
          panel.grid.major.y=element_line(), 
          panel.grid.minor.y=element_line(),
          panel.background=element_rect(fill="white"), 
          axis.text.x = element_text(hjust = 0.5, vjust = 1,size = 18, face="bold"),
          plot.title = element_text(size = 24, hjust = 0, face = 'bold'),
          axis.text.y = element_text(size = 18, face="bold"),
          axis.ticks = element_blank()) + 
    theme(strip.text.x = element_text(size = 20),
          strip.background = element_blank()) +
    theme(legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size=18),
          legend.title = element_text(size=20, margin = margin(0, 0, 9, 0))) +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(hjust = 0, size=26)) +
    theme(legend.key.width= unit(1.5, 'cm'))
} 

# function to process confusion matrices for plotting

confusion_processing <- function(x) {
  # matrix to long format
  x <- tibble::rownames_to_column(x %>% data.frame)
  x <- x %>% slice(1:(n()-1))
  x <- pivot_longer(x, cols=(!rowname))
  x$name <- gsub("X","", x$name)
  x$name <- str_replace_all(x$name, "\\.", ":")
  x$name <- str_replace_all(x$name, "UTR:", "UTR-")
  x$rowname <- str_replace_all(x$rowname, "UTR:", "UTR-")
  x <- x %>% mutate(name = sub("^((?:.*?:){3}.*?):", "\\1/", name))
  x$name <- str_replace_all(x$name, "008:01:04", "008:01/04") # correct allele naming in MICA
  x$name <- str_replace_all(x$name, c("009N" = "009:01N", "A314V" = "053", "R6C" = "052", "R6H" = "285", "D226N" = "193")) # change novel allele names and correct 009N in MICB
  x$rowname <- str_replace_all(x$rowname, c("009N" = "009:01N", "A314V" = "053", "R6C" = "052", "R6H" = "285", "D226N" = "193"))
  # add goodbad col for color fills
  x <- x %>% mutate(goodbad = ifelse(x$rowname == x$name, "good", "bad")) %>%
    group_by(rowname)
  x[x$value==0, 'goodbad'] <- 'absent'
  
  return(x)
}

# function for confusion matrix plotting (supplementary confusion matrices)
# input is a processed confusion data

conf_table_plot <- function(x, title=NA) {
  ggplot(x, aes(x = name, y = rowname, fill = goodbad)) +
    geom_tile(show.legend = FALSE, color="black") + 
    theme(axis.text.x = element_text(angle = 90)) +
    geom_text(aes(label=value), 
              color= ifelse(x$goodbad=="absent", "white", "black"), 
              position = position_dodge(width=0), 
              size=3.5, 
              hjust=0.7, 
              vjust=0.5) +
    scale_fill_manual(values = c(good = "green", bad = "red", absent="white")) +
    xlab("True") +
    ylab("Predicted") +
    ggtitle(title) %>% return
}

# function v3 for confusion matrix plotting (Figure 3 & 5 confusion matrices)
# input is a processed confusion data

conf_table_plot_3 <- function(x, title=NA) {
  ggplot(x, aes(x = name, y = rowname, fill = goodbad)) +
    geom_tile(show.legend = FALSE, color="black") + 
    geom_text(aes(label=value),
              color= ifelse(x$goodbad=="absent", "white", "black"),
              position = position_dodge(width=0),
              size=5,
              #fontface="bold",
              hjust=0.5,
              vjust=0.5) +
    scale_fill_manual(values = c(good = alpha("green", 0.5), bad = alpha("red", 0.4), absent="white")) +
    xlab("True") +
    ylab("Predicted") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 24),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 0.95)) %>% return
}


# function for plotting allelic results (accuracy, sensitivity, specificity etc)

plot_allele_prop <- function(results, parameter, ytitle, xtitle) { 
  gene.names <- c("MICA", "MICB", "HLA-E", "HLA-F")
  names(gene.names) <- c("MICA", "MICB", "HLA_E", "HLA_F")
  ggplot(results, aes(x= train.freq, y={{parameter}})) +
    geom_point(aes(shape=factor(Pop, levels=c("FIN", "EUR", "AFR", "EAS", "SAS", "AMR")), 
                   colour=factor(Ref, levels=c("FIN", "1000G", "combined")), stroke=1), 
               size = 3.7, alpha = 0.8, show.legend = T) +
    xlab(xtitle) + ylab(ytitle) +
    facet_wrap(vars(factor(Gene, levels=c("MICA", "MICB", "HLA_E", "HLA_F"))), 
               ncol = 1,
               labeller = as_labeller(gene.names)) +
    labs(shape = 'Population') +
    scale_shape_manual(values=c(0,1,2,4,5,6)) +
    scale_color_manual(values=c('#046C9A','#FF0000','#FDD262'), name = 'Reference') +
    theme_classic() +
    theme(axis.text.x = element_text(size = 18, angle = 0, vjust = 0.5), 
          axis.text.y = element_text(size=18), axis.title=element_text(size=18)) +
    theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size=18), axis.title=element_text(size=18)) +
    theme(strip.text.x = element_text(size = 24)) +
    theme(legend.title = element_text(size=24), legend.text = element_text(size=18)) +
    theme(panel.background=element_rect(fill='white', colour='black'),
          strip.background=element_rect(fill='white', colour='white')) +
    theme(plot.title = element_text(size = 22, hjust = 0, face = 'bold')) +
    theme(legend.position = "top") +
    guides(shape = guide_legend(override.aes = list(stroke = 1)))
}
