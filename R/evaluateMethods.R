#' @title Evaluated multiple normalization methods and multiple DE methods
#' 
#' @param counts Gene by sample expression count matrix (G by N). 
#'               Use raw count data before filtering.
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param null binary vector of TRUE/FALSE the gene a null gene.
#' @param methodsNormalize Chararacter vector of evaluted methods. To run all methods, use
#'                            c("normalize.cpm", "normalize.tmm", "normalize.rle",
#'                             "normalize.census", "normalize.scnorm", "normalize.scran")
#' @param methodsMeanExpression Chararacter vector of evaluted methods. To run all methods, use
#'                             c("DESeq2", "limmaVoom", "edgeR","BPSC", "MAST", "ROTS")
#' 
#' @return 
#'     \code{pvals_longformat} data.frame of pvals.
#'
#' @export
query.pipeline <- function(counts, condition, null,
                           methodsNormalize, methodsMeanExpression) {

  #----- filtering
  counts_filtered <- filter.excludeAllZeros(counts)
  featuresToInclude <- filterFeatures.fractionExpressed(counts_filtered,
                                                       thresholdDetection = 1,
                                                       fractionExpressed = .01)$index_filter
  
  samplesToInclude <-  filterSamples.fractionExpressed(counts_filtered,
                                                       thresholdDetection = 1,
                                                       fractionExpressed = .01)$index_filter
  
  counts_filtered <- counts_filtered[featuresToInclude, samplesToInclude]
  
  #----- normalization
  methodsNormalize <- c("TMM", "RLE", "census","scran")
  libsize_factors_list <- query.methodsNormalization(counts = counts_filtered, 
                                            condition = condition,
                                            methodsNormalize = methodsNormalize)
  
  #---- run DE methods
  pvals_list <- vector("list", length = length(methodsNormalize))
  names(pvals_list) <- methodsNormalize
  for (index in 1:length(methodsNormalize)) {
    counts_normed <- normalize.cpm(counts_filtered, libsize_factors)$cpm
    libsize_factors <- libsize_factors_list[[index]]
    pvals_list[[index]] <- query.methodsMeanExpression(
                                  counts = counts_filtered,
                                  counts_normed = counts_normed,
                                  condition = condition,
                                  libsize_factors = libsize_factors,
                                  methodsMeanExpression = c("limmaVoom", "edgeR", "MAST"))
  }
  
  
  # make list into a long-format data.frame
  pvals_longformat <- do.call(rbind, lapply(1:length(pvals_list), function(index) {
    pvals_list[[index]] <- data.frame(pvals_list[[index]])
    pvals_list[[index]]$methodsNormalize <- names(pvals_list)[index]
    return(pvals_list[[index]])
  }) )
  
  return(pvals_longformat)
  
}
  

#' @title Evaluated multiple normalization methods
#' 
#' @param counts Gene by sample expression count matrix (G by N). Use filtered count data.
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param methodsNormalize Chararacter vector of evaluted methods. To run all methods, use
#'                            c("normalize.cpm", "normalize.tmm", "normalize.rle",
#'                             "normalize.census", "normalize.scnorm", "normalize.scran")
#'
#' @return 
#'    \code{libsize_factors} List of multiple size factors.

#' @example
#' ipsc_eset <- get(load(system.file("testdata", "HumanTungiPSC.rda", package = "ashbun")))
#' counts <- exprs(ipsc_eset)[sample(nrow(exprs(ipsc_eset)), 500), ]
#' condition <- pData(ipsc_eset)$replicate
#' 
#' ----- Step 1: filtering
#' counts_filtered <- filter.excludeAllZeros(counts)
#' featuresToInclude <- filterFeatures.fractionExpressed(counts_filtered, 
#'                                                      thresholdDetection = 1,
#'                                                      fractionExpressed = .01)$index_filter
#'                                                    
#' samplesToInclude <-  filterSamples.fractionExpressed(counts_filtered, 
#'                                                      thresholdDetection = 1,
#'                                                      fractionExpressed = .01)$index_filter
#'                                                      
#' counts_filtered <- counts_filtered[featuresToInclude, samplesToInclude]
#'
#' ---- Step 2: compute library size factors
#' sizefactors <- query.methodsNormalization(counts_filtered, condition = condition,
#'                                           methodsNormalize = c("TMM", "RLE",
#'                                                                 "census","scran"))
#' 
#' @author Chiaowen Joyce Hsiao
#'
#' @export
query.methodsNormalization <- function(counts, condition,
                                       methodsNormalize = c("TMM", "RLE",
                                                            "census", "SCnorm", "scran")) {
  
  # SCnorm is the only method that requires the input of condition/group variable
  # and is the only method that def. has different scale factors for different
  # gene groups
  methods_same_for_genes <- c("TMM", "RLE", "census", "scran")
  methods_same_for_genes_function <- c( "normalize.tmm",
                                        "normalize.rle",
                                        "normalize.census",
                                        "normalize.scran")
  which_methods_same_for_genes <- which(methods_same_for_genes %in% methodsNormalize)
  
  methods_diff_for_genes <- c("SCnorm")
  methods_diff_for_genes_function <- "normalize.scnorm"  
  which_methods_diff_for_genes <- which(methods_diff_for_genes %in% methodsNormalize)
  
  # run methods that apply the same scale factors for every gene
  scalefactors_same_for_genes <- vector("list", 
                                        length = length(which_methods_same_for_genes))
  names(scalefactors_same_for_genes) <- methods_same_for_genes[which_methods_same_for_genes]
  for (index in seq_along(which_methods_same_for_genes)) {
    index_method <- which_methods_same_for_genes[index]
    cat(methods_same_for_genes[index_method], "\n")
    output <- do.call(methods_same_for_genes_function[index_method], 
                      list(counts = counts))
    scalefactors_same_for_genes[[index]] <- output$libsize_factors
  }
  
  # run methods that apply different scale factors for every gene
  scalefactors_diff_for_genes <- vector("list", 
                                        length = length(which_methods_diff_for_genes))
  names(scalefactors_diff_for_genes) <- methods_diff_for_genes[which_methods_diff_for_genes]
  for (index in seq_along(which_methods_diff_for_genes)) {
    index_method <- which_methods_diff_for_genes[index]
    cat(methods_diff_for_genes[index_method], "\n")
    output <- do.call(methods_diff_for_genes_function[index_method], 
                      list(counts = counts, condition = condition))
    scalefactors_diff_for_genes[[index]] <- output
  }
  
  
  libsize_factors <- list(scalefactors_same_for_genes = scalefactors_same_for_genes,
                          scalefactors_diff_for_genes = scalefactors_diff_for_genes)
  libsize_factors <- c(scalefactors_same_for_genes,
                       scalefactors_diff_for_genes)
  
  return(libsize_factors)
}





#' @title Evaluated multiple DE methods
#' 
#' @param counts Gene by sample expression count matrix (G by N). Use filtered count data.
#' @param counts_normed Normalized expression count matrix 
#'                      (typicall CPM with normlized library size).
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param libsize_factors Numeric vector of scale factors for library size factors.
#' @param methodsMeanExpression Chararacter vector of evaluted methods. To run all methods, use
#'                             c("DESeq2", "limmaVoom", "edgeR","BPSC", "MAST", "ROTS")
#'
#' @param control
#'   \code{pseudocount} Default .5. For limmaVoom and MAST. If NULL, then do not use pseudocount.
#'
#' @return 
#'    \code{pvalues} data.frame of significance values. Columns corresond to input methods.

#' @example
#' ipsc_eset <- get(load(system.file("testdata", "HumanTungiPSC.rda", package = "ashbun")))
#' counts <- exprs(ipsc_eset)[sample(nrow(exprs(ipsc_eset)), 500), ]
#' condition <- pData(ipsc_eset)$replicate
#' 
#' ----- Step 1: filtering
#' counts_filtered <- filter.excludeAllZeros(counts)
#' featuresToInclude <- filterFeatures.fractionExpressed(counts_filtered, 
#'                                                      thresholdDetection = 1,
#'                                                      fractionExpressed = .01)$index_filter
#'                                                    
#' samplesToInclude <-  filterSamples.fractionExpressed(counts_filtered, 
#'                                                      thresholdDetection = 1,
#'                                                      fractionExpressed = .01)$index_filter
#'                                                      
#' counts_filtered <- counts_filtered[featuresToInclude, samplesToInclude]
#'
#' ---- Step 2: compute library size factors
#' libsize_factors <- normalize.scran(counts = counts_filtered)$libsize_factors
#' counts_normed <- normalize.cpm(counts_filtered, libsize_factors)$cpm
#' 
#' ---- Step 3: run DE methods
#' pvals_list <- query.methodsMeanExpression(counts = counts_filtered,
#'                                           counts_normed = counts_normed,
#'                                           condition = condition,
#'                                           libsize_factors = libsize_factors,
#'                                           methodsMeanExpression = c("limmaVoom", "edgeR",
#'                                                                     "MAST"))
#' @author Chiaowen Joyce Hsiao
#'
#' @export
query.methodsMeanExpression <- function(counts, counts_normed, condition,
                                  libsize_factors, 
                                  methodsMeanExpression = c("DESeq2", "limmaVoom", "edgeR",
                                                            "BPSC", "MAST", "ROTS")) {
  # match methods that input raw counts
  methods_raw <- c("DESeq2", "edgeR")
  methods_raw_function <- c("methodWrapper.DESeq2",
                             "methodWrapper.edgeR")
  which_methods_raw <- which(methods_raw %in% methodsMeanExpression)
  
  # match methods that input normalied expression count
  methods_normed <- c("limmaVoom", "BPSC", "MAST", "ROTS")
  methods_normed_function <- c("methodWrapper.limmaVoom",
                               "methodWrapper.bpsc",
                               "methodWrapper.mast",
                               "methodWrapper.rots")
  which_methods_normed <- which(methods_normed %in% methodsMeanExpression)
  
  # run methods for raw counts
  pvals_raw <- vector("list", length = length(which_methods_raw))
  names(pvals_raw) <- methods_raw[which_methods_raw]
  for (index in seq_along(which_methods_raw)) {
    index_method <- which_methods_raw[index]
    cat(methods_raw[index_method], "\n")
    output <- do.call(methods_raw_function[index_method], 
                      list(counts = counts, 
                           condition = condition,
                           libsize_factors = libsize_factors))
    pvals_raw[[index]] <- output$pvalue
  }

  # methods for normalized counts
  pvals_normed <- vector("list", length = length(which_methods_normed))
  names(pvals_normed) <- methods_normed[which_methods_normed]
  for (index in seq_along(which_methods_normed)) {
    index_method <- which_methods_normed[index]
    cat(methods_normed[index_method], "\n")
    # pseudcount for mast and limma is set at .5.
    output <- do.call(methods_normed_function[index_method], 
                      list(counts = counts_normed, 
                           condition = condition))
    pvals_normed[[index]] <- output$pvalue
  }
  
  pvalues <- cbind(do.call(cbind, pvals_raw), 
                   do.call(cbind, pvals_normed))
  
  return(pvalues)
}

  



  
  