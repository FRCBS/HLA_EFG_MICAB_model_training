### ------------------------------------------------------------------------------------------------------------------------------ ### 

#             functions for model fitting and cross-validation

### ------------------------------------------------------------------------------------------------------------------------------ ###

### function to load models listed in model_files

load_models <- function(filelist, locus) {
  # filelist = list of model files
  sapply(filelist, function (x) {
    models <- list(hlaModelFromObj(get(load(paste0("./data/Imputation_models/Combination_testing_models/",locus,"/",x)))))
    return(models)
  } )
} 

### function for model fitting

fitHLAmodel <- function(gene, pheno, genotype, region, n) {
  # gene = gene name
  # pheno = training phenotypes, an object of hlaAlleleClass
  # genotype = training SNP genotypes, an object of hlaSNPGenoClass
  # region = flanking bp
  # n = total number of individual classifiers
  
  oma.cl <- makeCluster(7)
  tmp.snps <- hlaFlankingSNP(genotype$snp.id, genotype$snp.position, gene, 
                             flank.bp = region, assembly = "hg38")
  genotypes <- hlaGenoSubset(genotype, snp.sel = match(tmp.snps, genotype$snp.id),
                             samp.sel = genotype$sample.id %in% pheno$value$sample.id)
  file_path <- paste0("./data/Imputation_models/Combination_testing_models/model_", gene, "_", deparse(substitute(genotype)),
                      "_", deparse(substitute(pheno)),".RData")
  set.seed(100)
  out <- hlaParallelAttrBagging(cl=oma.cl, pheno, genotypes, nclassifier = n, 
                                auto.save =file_path)
  
}

### function for prediction and calculation of accuracy

PredHlaModel <- function(model, genotypes, validation) {
  # model = imputation model used in the prediction
  # genotypes = test SNP genotypes, an object of hlaSNPGenoClass
  # validation = test phenotypes, an object of hlaAlleleClass
  
  genotype.set <- hlaGenoSubset(genotypes, samp.sel=match(validation$value$sample.id, genotypes$sample.id))
  pred.results <- hlaPredict(model, genotype.set, type="response+prob", match.type="Position")
  accuracy <- hlaCompareAllele(validation, pred.results, allele.limit = NULL, output.individual = TRUE)
  return(c(accuracy, pred.results))
  
}

### function to get overall accuracies from flanking region test prediction loops
accuracy_loop <- function(lst, n){
  sapply(lst, `[`, n)
}


### function to impute FIN test set (FIN I/FIN II) using all models (models I-VII)

pred_FIN <- function(model.list, geno.data, validation.set) {
  # model.list = list of objects of hlaAttrBagClass
  # geno.data =  test SNP genotypes, an object of hlaSNPGenoClass
  # validation.set = test phenotypes, an object of hlaAlleleClass
  
  map(1:length(model.list), function (x) {
    a <-PredHlaModel(model.list[[x]], geno.data, validation.set)
    return(a)
  })
}

### function to predict all superpopulations (EUR, AFR, EAS, SAS, AMR) using all models (I-VII)
pred_1000G_superpop <- function(model.list, val.superpop) {
  # model.list = list of objects of hlaAttrBagClass
  # val.superpop = list of objects of hlaAlleleClass
  
  map(1:length(model.list), function (x) {
    map(1:length(val.superpop), function (y) {
      out <- PredHlaModel(model.list[[x]], Genotypedata_1000G, val.superpop[[y]])
      return(out)
    })
  })
}
