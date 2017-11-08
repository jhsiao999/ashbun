# evaluateMethods.R
#
# Description:
#   This file contains high-level functions that combines and performs multiple methods,
#   including methods for DE, normalization, and evaluation.
#
# Content (from higher to lower-level functions)
#
# query.evaluation: inputs a single count table and outputs evaluation results of
#                   multiple normalization methods and DE methods
#
# query.pipeline: inputs a single count table, and outputs processed data and p-values of
#                  multiple normalization methods and DE methods
#
#
# query.methodsNormalization: input a single count table and outputs library size factors of
#                   multiple normalization methods.
#
# query.methodsMeanExpression: inputs a single count table and outputs p-values of
#                   multiple DE methods
#

#' @title Evaluate multiple normalization methods and multiple DE methods
#'
#' @param thresholdDetection minimum count per gene/sample. Default value = 1.
#' @param fractionExpressed fraction of samples expressed (above thresholdDetection).
#'        Default value = .01.
#'
#' @examples
#' ipsc_eset <- get(load(system.file("testdata", "HumanTungiPSC.rda", package = "ashbun")))
#' counts <- exprs(ipsc_eset)[sample(nrow(exprs(ipsc_eset)), ), ]
#'
#' #---- generat simulated datasets
#' library(ashbun)
#' simdata_list <- simulationWrapper(counts, Nsim = 2,
#'                                  Ngenes = 100,
#'                                  Nsam = 20,
#'                                  sample_method = "all_genes",
#'                                  pi0 = .5,
#'                                  beta_args = args.big_normal(betapi = 1,
#'                                                              betamu = 0, betasd = .8))
#'
#' #---- extract a single dataset as an example
#' #---- take pi0 = .9, the first simulated data
#' simdata <- simdata_list[[1]]
#'
#' # ---- gather evaluation results
#' eval_ouptut <- query.evaluation(counts = simdata$counts,
#'                                 condition = simdata$condition,
#'                                 is_nullgene = simdata$is_nullgene,
#'                                 methodsNormalize = c("TMM", "RLE"),
#'                                 methodsMeanExpression = c("DESeq2", "limmaVoom"),
#'                                 report = "fdr_cutoff_summary")
#'
#' @author Chiaowen Joyce Hsiao
#'
#' @export
query.evaluation <- function(counts, condition, is_nullgene,
                             methodsNormalize = c("LIB", "TMM", "RLE","census","SCnorm","scran"),
                             methodsMeanExpression = c("DESeq2", "limmaVoom"),
#                             thresholdDetection = 1, fractionExpressed = .01,
                             report.control = list(fdr_cutoff = .05), nsim = NULL) {

  results <- query.pipeline(counts, condition, is_nullgene,
                  methodsNormalize = methodsNormalize,
                  methodsMeanExpression = methodsMeanExpression)
  message("Analysis done!", "\n")

  num_evals <- dim(results$pvals_longformat)[1]/length(results$data$is_nullgene)
  df_summarize <- results$pvals_longformat
  df_summarize$is_nullgene <- rep(results$data$is_nullgene, num_evals)

  message("Evaluating results", "\n")

  output <- vector("list", 2)
  names(output) <- c("fdr_cutoff", "roc")

  # find true positive rate given false discovery rate .05,
    suppressPackageStartupMessages(library(dplyr))
    output$fdr_cutoff <- df_summarize %>%
      group_by(methodsNormalize, methodsMeanExpression) %>%
      summarise(tpr = getTPR.pROC(response = is_nullgene,
                                  predictor = pvalues,
                                  fdr_cutoff = .05))

  # compute receiver operating curve
    list_methodsNormalize <- unique(results$pvals_longformat$methodsNormalize)
    list_methodsMeanExpression <- unique(results$pvals_longformat$methodsMeanExpression)
    is_nullgene <- results$data$is_nullgene

    output$roc <- do.call(rbind, lapply(seq_along(list_methodsNormalize),
                                    function(index_normalize) {

      one_methodsNormalize <- do.call(rbind, lapply(seq_along(list_methodsMeanExpression),
                                                    function(index_meanExpression) {
            df_sub <- subset(results$pvals_longformat,
                           methodsNormalize == list_methodsNormalize[index_normalize] &
                           methodsMeanExpression == list_methodsMeanExpression[index_meanExpression])
            pvals_sub <- df_sub$pvalues
            roc_output <- getROC(response = is_nullgene, predictor = pvals_sub)
            foo <- data.frame(TPR = roc_output$sensitivities,
                              FPR = 1- roc_output$specificities,
                              methodsMeanExpression = list_methodsMeanExpression[index_meanExpression],
                              methodsNormalize = list_methodsNormalize[index_normalize])
            foo <- foo[foo$FPR < .2, ]
            return(foo)
      }) )

    }) )
    message("DONE!", "\n")

    output$fdr_cutoff$nsim <- nsim
    output$roc$nsim <- nsim

    return(output)
  }



#' @title Run multiple normalization methods and multiple DE methods
#'
#' @param counts Gene by sample expression count matrix (G by N).
#'               Use raw count data before filtering.
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param null binary indicator of true/false. True = Non-null gene and FALSE = Null gene.
#' @param methodsNormalize Chararacter vector of evaluted methods. To run all methods, use
#'                            c("normalize.cpm", "normalize.tmm", "normalize.rle",
#'                             "normalize.census", "normalize.scnorm", "normalize.scran")
#' @param methodsMeanExpression Chararacter vector of evaluted methods. To run all methods, use
#'                             c("DESeq2", "limmaVoom", "edgeR","BPSC", "MAST", "ROTS")
#'
#' @examples
#' ipsc_eset <- get(load(system.file("testdata", "HumanTungiPSC.rda", package = "ashbun")))
#' counts <- exprs(ipsc_eset)[sample(nrow(exprs(ipsc_eset)), 500), ]
#' condition <- pData(ipsc_eset)$replicate
#'
#' results <- query.pipeline(counts = counts,
#'                           condition = condition,
#'                           is_nullgene = NULL,
#'                           methodsNormalize = c("TMM", "RLE", "census","scran"),
#'                           methodsMeanExpression = c("DESeq2", "limmaVoom"))
#'
#' @return
#'     \code{data} List of filtered data, including count matrix, sample condition vector,
#'                  and logical vector for null gene status (TRUE if null).
#'     \code{pvals_longformat} data.frame of pvals.
#'
#' @export
query.pipeline <- function(counts, condition, is_nullgene = NULL,
                           methodsNormalize = c("LIB", "TMM", "RLE", "census","SCnorm","scran"),
                           methodsMeanExpression = c("DESeq2", "limmaVoom", "edgeR",
                                                     "BPSC", "MAST", "ROTS", "scde")) {
#                           thresholdDetection = 1, fractionExpressed = .01) {

  data_filtered <- list(counts = counts,
                        condition = condition,
                        is_nullgene = is_nullgene)

  #----- normalization
  libsize_factors_list <- with(data_filtered,
                               query.methodsNormalization(counts = counts,
                                                          condition = condition,
                                                          methodsNormalize = methodsNormalize))

  #---- run DE methods
  counts_filtered <- data_filtered$counts
  condition_filtered <- data_filtered$condition
  pvals_list <- vector("list", length = length(methodsNormalize))
  names(pvals_list) <- methodsNormalize

  # consider methods that use that same scaling factors for each gene
  for (index in 1:length(methodsNormalize)) {

    counts_normed <- normalize.cpm(counts_filtered, libsize_factors)$cpm
    libsize_factors <- libsize_factors_list$scalefactors_same_for_genes[[index]]
    pvals_list[[index]] <- query.methodsMeanExpression(
                                  counts = counts_filtered,
                                  counts_normed = counts_normed,
                                  condition = condition_filtered,
                                  libsize_factors = libsize_factors,
                                  methodsMeanExpression = methodsMeanExpression)
  }


  # make list into a long-format data.frame
  pvals_longformat <- do.call(rbind, lapply(1:length(pvals_list), function(index) {
    obj <- data.frame(pvals_list[[index]])
    obj_name <- names(pvals_list)[[index]]

    obj_transformed <- as.data.frame(do.call(rbind,
        lapply(1:dim(pvals_list[[index]])[2], function(index_sub) {
              foo <- cbind(pvalues = as.numeric(obj[[index_sub]]),
                    methodsMeanExpression = colnames(pvals_list[[index]])[[index_sub]])
              return(foo)
      })),
      stringsAsFactors = FALSE)
    obj_transformed$pvalues <- as.numeric(obj_transformed$pvalues)

    obj_transformed$methodsNormalize <- obj_name
    return(obj_transformed)
  }) )

  return(list(data = data_filtered,
              pvals_longformat = pvals_longformat))
}


#' @title Run multiple normalization methods
#'
#' @param counts Gene by sample expression count matrix (G by N). Use filtered count data.
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param methodsNormalize Chararacter vector of evaluted methods. To run all methods, use
#'                            c("normalize.cpm", "normalize.tmm", "normalize.rle",
#'                             "normalize.census", "normalize.scnorm", "normalize.scran")
#'
#' @return
#'    \code{libsize_factors} List of multiple size factors.

#' @examples
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
#' condition_filtered <- condition[samplesToInclude]
#'
#' ---- Step 2: compute library size factors
#' sizefactors <- query.methodsNormalization(counts_filtered,
#'                                           condition = condition_filtered,
#'                                           methodsNormalize = c("TMM", "RLE",
#'                                                                 "census","scran"))
#'
#' @author Chiaowen Joyce Hsiao
#'
#' @export
query.methodsNormalization <- function(counts, condition,
                                       methodsNormalize = c("LIB", "TMM", "RLE",
                                                            "census", "SCnorm", "scran")) {

  # SCnorm is the only method that requires the input of condition/group variable
  # and is the only method that def. has different scale factors for different
  # gene groups
  methods_same_for_genes <- c("LIB", "TMM", "RLE", "census", "scran")
  methods_same_for_genes_function <- c( "normalize.lib",
                                        "normalize.tmm",
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

  if (length(which_methods_diff_for_genes) == 0) {
    libsize_factors <- list(scalefactors_same_for_genes = scalefactors_same_for_genes)
   }

  if (length(which_methods_diff_for_genes) != 0) {
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
  }

  return(libsize_factors)
}





#' @title Runs multiple DE methods
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

#' @examples
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
#'                                           condition = condition_filtered,
#'                                           libsize_factors = libsize_factors,
#'                                           methodsMeanExpression = c("limmaVoom",
#'                                                                     "DESeq2",
#'                                                                     "edgeR",
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
