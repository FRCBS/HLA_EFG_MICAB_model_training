# function to extract and process true alleles

true_function <- function(df, locus, pop) {
  if (pop == 'FIN') {
    res <- df$individual %>% select(c(sample.id, true.hla))
  }
  else {
    res <- df$vi$individual %>% select(c(sample.id, true.hla))
  }
  if (locus == 'G'){
    res$true.hla <- str_replace_all(res$true.hla, '01:01:01:01/02', '01:01:01:01or02')
  }
  else {
    res$true.hla <- str_replace_all(res$true.hla, '008:01/04', '008:01or04')
  }
  res <- res %>% separate_wider_delim(true.hla, delim = "/", names = c("allele1", "allele2")) %>% as.data.frame()
  res <- res %>% mutate(allele1 = paste0(gsub('HLA_', '', locus), '*', allele1)) %>% mutate(allele2 = paste0(gsub('HLA_', '', locus),'*', allele2))
  return(res)
}

# function to get true allele dosages

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

# function for allele frequencies

freq_function <- function(df_res, df_cor, locus, name, pop) {
  
  if (name == 'G'| pop == 'FIN') {
    freq_table <- df_res[[locus]]$detail
  }
  else {
    freq_table <- df_res[[locus]]$vi$detail
  }
  for (row in 1:nrow(df_cor)) {
    allele <- df_cor[row, 'allele'] 
    allele <- gsub(paste0(name, "*"), "", fixed=T, allele)
    allele_freq <- freq_table[freq_table$allele == allele, ]
    df_cor[row, 'freq'] <- allele_freq$train.freq
  }
  return(df_cor)
}