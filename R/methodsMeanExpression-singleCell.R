# This document contains code for implementing DE methods developed for single cell RNA-seq data

#' @title BPSC
#'
#' @description Implement BPSC 0.99.1 (not available via Bioconductor). This method
#'              requires that the input data is normalized (counts-per-million or
#'              fragments per-kilo-base per million).
#'
#' @param counts Gene by sample expression count matrix (G by N). BPSC requires this input
#'               to be normalized expression matrix, such as FPKM or CPM.
#' @param condition Binary vector of length N indicating sample biological condition.
#' #' @param libsize_factors Numeric vector of scale factors for library sizes.
#'                       Default NULL multiplies all library sizes by 1.

#' @param control A list with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete DESeq2 fit results.
#'   \code{estIntPar} BPSC internal argument (TRUE/FALSE) for using only the expressed
#'                    genes to compute parameter estimates.
#'   \code{useParalle} BPSC internal argument (TRUE/FALSE) for using code written for
#'                    parallel computing.
#'
#' @return List with the following objects
#'    \code{pvalue} significance value for the group effect.
#'    \code{fit} BPSC complete output of the model fit.

#' @author Chiaowen Joyce Hsiao
#'
#' @export
methodWrapper.bpsc <- function(counts, condition,
                               default = FALSE,
                               control = list(save_modelFit = FALSE,
                                              estIntPar = FALSE,
                                              useParallel = TRUE)) {

  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(dim(counts)[2] == length(condition))

  counts <- as.matrix(counts)
  if (!is.factor(condition)) {condition <- factor(condition)}

  #<--------------------------------------
  # First, choose IDs of all cells of the control group for estimating
  # parameters of BP models
  controlIds <- which(as.integer(condition) == 1)

  # Create a design matrix including hte group labels, batch effectd
  # can be added here if they are available
  design <- model.matrix(~condition)
  
  if (default == FALSE) {
    # Model fitting; estIntPar requires to use only the expressed genes
    # to compute estiamtes, according to the vignette, this option requires longer
    # computational time
    fit <- BPSC::BPglm(data = counts,
                       controlIds = controlIds,
                       design = design,
                       coef = 2,
                       estIntPar = control$estIntPar,
                       useParallel = control$useParallel)
  }

  if (default == TRUE) {
    # Model fitting; estIntPar requires to use only the expressed genes
    # to compute estiamtes, according to the vignette, this option requires longer
    # computational time
    counts.cpm <- normalize.cpm(counts)
    fit <- BPSC::BPglm(data = counts.cpm$cpm,
                       controlIds = controlIds,
                       design = design,
                       coef = 2,
                       estIntPar = TRUE,
                       useParallel = control$useParallel)
  }
  
  # extract results
  pvalue <- fit$PVAL

  # if save_modelFit, then output will include the original model fit
  if (control$save_modelFit) {
    fit <- fit
  } else {
    fit <- NULL
  }

  return(list(betahat=rep(NA, length(pvalue)),
              sebetahat=rep(NA, length(pvalue)),
              df=rep(NA, length(pvalue)),
              pvalue = pvalue,
              fit=fit))
}


#' @title MAST
#'
#' @describeIn Implement MAST.
#'
#' @param counts Gene by sample expression count matrix (G by N). MAST requires the default
#'                input to be log2CPM. In this function, the required input is
#'                normalized expression counts, such as CPM. log2CPM is computed internally.
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param pseudocount Default .5.
#' @param control A list with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete DESeq2 fit results.
#'   \code{include_cdr} MAST internal argument. TRUE if include cellular detection
#'                      rate (fraction of genes detected per sample) as a covariate
#'                      (scaled to mean 0 and standard deviation 1).
#'
#' @return List with the following objects
#'    \code{betahat} The estimate effect size of condition for all genes.
#'    \code{sebetahat} The standard errors of the effect sizes.
#'    \code{df} The degrees of freedom associated with the effect sizes.
#'    \code{pvalue} P-values of the effect sizes from the likelihood ratio test of
#'                  the hurdle model, which combines the continous and the discrete
#'                  component (the above betahat, sebetahat, and df are extracted from
#'                  the continuous component.)
#'    \code{fit} MAST complete output of the model fit.

#' @author Chiaowen Joyce Hsiao
#'
#' @export
methodWrapper.mast <- function(counts, condition, default = FALSE,
                               pseudocount = .5,
                               control = list(save_modelFit = FALSE,
                                              include_cdr = TRUE)) {
  
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(dim(counts)[2] == length(condition))
  
  counts <- as.matrix(counts)
  if (!is.factor(condition)) {condition <- factor(condition)}
  
  suppressPackageStartupMessages(library(MAST))
  
  #<--------------------------------------
  if (default == TRUE) {
    # compute log2CPM
    counts.cpm <- normalize.cpm(counts)$cpm
    
    # make data.frame into singleCellAssay object
    colData <- data.frame(condition = condition)
    if (!is.null(rownames(counts))) {
      rowData <- data.frame(gene = rownames(counts))
    }
    if (is.null(rownames(counts))) {
      rowData <- data.frame(gene = paste0("gene.", c(1:nrow(counts))))
    }
    rowData$gene <- rowData$gene
    sca <- MAST::FromMatrix(counts.cpm, colData, rowData)
    
    # adaptive threshold in MAST
    thres <- thresholdSCRNACountMatrix(assay(sca), data_log=FALSE)
    
    # MAST vignette default
    freq_expressed <- .2
    
    #    assays(sca) <- list(thresh=thres$counts_threshold, cpm=assay(sca))
    expressed_genes <- freq(sca) > freq_expressed
    sca <- sca[expressed_genes,]
    
    # calculate cellualr detection rate; normalized to mean 0 and sd 1
    colData(sca)$cdr.normed <- scale(colSums(assay(sca) > 0))
    
    # the default method for fitting is bayesGLM
    fit <- MAST::zlm(~ condition + cdr.normed, sca)
  }
  
  if (default == FALSE) {
    # compute log2CPM
    log2CPM <- log2(counts + pseudocount)
    
    # make data.frame into singleCellAssay object
    colData <- data.frame(condition = condition)
    rowData <- data.frame(gene = rownames(counts))
    sca <- MAST::FromMatrix(log2CPM, colData, rowData)
    
    if (control$include_cdr) {
      # calculate cellualr detection rate; normalized to mean 0 and sd 1
      colData(sca)$cdr.normed <- scale(colSums(assay(sca) > 0))
      # the default method for fitting is bayesGLM
      fit <- MAST::zlm(~ condition + cdr.normed, sca)
    } else {
      fit <-  MAST::zlm(~ condition, sca)
    }
  }
  
  # LRT test for the significance of the condition effect
  lrt <-  MAST::lrTest(fit, "condition")
  
  # extract p.value from their "hurdle" model
  pvalue <- lrt[,3,3]
  
  # extract effect size, standard error, and df from the non-zero component
  betahat <- fit@coefC[,2]
  setbetahat <- sqrt(sapply(1:dim(fit@vcovC)[3], function(i) {
    diag(fit@vcovC[,,i])[2]} ) )
  df <- fit@df.resid[,1]
  
  # if save_modelFit, then output will include the original model fit
  if (control$save_modelFit) {
    fit <- fit
  } else {
    fit <- NULL
  }
  
  return(list(betahat=betahat,
              sebetahat=setbetahat,
              df=df,
              pvalue=pvalue,
              fit=fit))
}








#' @title ROTS
#'
#' @describeIn Implement ROTS 1.2.0.
#'
#' @param counts Gene by sample expression count matrix (G by N). ROTS requires this input
#'               to be normalized expression matrix, such as FPKM or CPM.
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param control List with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete DESeq2 fit results.
#'   \code{B} ROTS internal argument: Number of bootstraps. Default 100.
#'   \code{seed} ROTS internal argument: An integer seed for the random number generator.
#'               Default NULL.
#'
#' @return A list with the following objects
#'    \code{sig_order} Z-socres of expression difference, ordered from the most
#'                     significant to the least significant.
#'    \code{fit} ROTS complete output of the model fit.

#' @author Chiaowen Joyce Hsiao
#'
#' @export
methodWrapper.rots <- function(counts, condition,
                               control = list(save_modelFit = FALSE,
                                              seed = NULL,
                                              B = 100)) {

  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(dim(counts)[2] == length(condition))
  assertthat::assert_that(is.numeric(condition))

  counts <- as.matrix(counts)

  #<--------------------------------------
  suppressPackageStartupMessages(library(ROTS))
  
  # fitting model
  fit <- ROTS::ROTS(data = counts,
                     groups = condition,
                      B = control$B,
                      K = dim(counts)[1],
                      seed = NULL,
                      log = FALSE) 

  # extract results
  pvalue <- fit$pvalue

  # if save_modelFit, then output will include the original model fit
  if (control$save_modelFit) {
    fit <- fit
  } else {
    fit <- NULL
  }

  return(list(pvalue = pvalue,
              fit=fit))
}




#' @title SCDE
#'
#' @description Implement SCDE 1.99.2. This older version of SCDE is need for computers
#'               that are compatible with an older version of flexmix (flexmix_2.3-13;
#'               for related discussions on this, see https://github.com/hms-dbmi/scde/issues/40).
#'               Currently, this function uses the SCDE method for normalizatino and
#'               does not allow using other normalization factors.
#'
#' @param counts Gene by sample expression count matrix (G by N).
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param control List with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete DESeq2 fit results.
#'   \code{n.randomization} SCDE internal argument: Nubmer of bootstraps.
#'   \code{n.cores} SCDE internal argument: Nubmer of cores.
#'   \code{min.size.entries} SCDE internal argument: Mininum number of genes to use when
#'                            determining expected expression magnitude during model fitting.
#'
#' @return List with the following objects
#'    \code{sig_order} Z-socres of expression difference, ordered from the most
#'                     significant to the least significant.
#'    \code{fit} SCDE complete output of the model fit.

#' @author Chiaowen Joyce Hsiao
#'
#' @export
methodWrapper.scde <- function(counts, condition,
                               control = list(save_modelFit = FALSE,
                                              n.randomizations = 100,
                                              n.cores = 4,
                                              min.size.entries = ncol(counts))) {

  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(counts))
  assertthat::assert_that(dim(counts)[2] == length(condition))

  # convert condition to factor
  if (!is.factor(condition)) {condition <- factor(condition)}

  # convert count table to integer
  if (!is.integer(counts)) {
    counts <- apply(counts, 2, function(x) { storage.mode(x) <- 'integer'; x})
  }
  #<--------------------------------------
  # fitting error models, this step takes a considerable amount of time...
  o.ifm <- scde::scde.error.models(counts = counts,
                                   groups = condition,
                                   n.cores = control$n.cores,
                                   min.count.threshold = 1,
                                   threshold.segmentation = TRUE,
                                   save.crossfit.plots = FALSE,
                                   save.model.plots = FALSE,
                                   min.size.entries = control$min.size.entries,
                                   max.pairs = (min(table(condition)))^2,
                                   verbose = 0)

  # filter out cells that dont't show positive correlation with
  # the expected expression magnitudes
  o.ifm <- o.ifm[o.ifm$corr.a > 0, ]

  # estimate gene expression prior
  o.prior <- scde::scde.expression.prior(models = o.ifm,
                                         counts = counts,
                                         length.out = 400, show.plot = FALSE)

  # DE analysis
  fit <- scde::scde.expression.difference(o.ifm, counts, o.prior,
                                          groups  =  as.factor(condition),
                                          n.randomizations  =  control$n.randomizations,
                                          # number of bootstrap randomizations to be performed
                                          n.cores  =  control$n.cores, verbose  =  0)

  # order results by magnitude of statistical signifcance
  sig_order <- fit[order(abs(fit$Z), decreasing  =  TRUE), ]

  # if save_modelFit, then output will include the original model fit
  if (control$save_modelFit) {
    fit <- fit
  } else {
    fit <- NULL
  }

  return(list(sig_order = sig_order,
              fit=fit))
}
