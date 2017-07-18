# This document contains code for implementing DE methods developed for single cell RNA-seq data



#' @title MAST
#' 
#' @describeIn Implement MAST.
#' 
#' @param count_matrix Gene by sample expression count matrix (G by N).
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param control A list with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete DESeq2 fit results.
#'   \code{include_cdr} MAST internal argument. TRUE if include cellular detection 
#'                      rate (fraction of genes detected per sample) as a covariate 
#'                      (scaled to mean 0 and standard deviation 1).
#'
#' @return A list with the following objects
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
methodWrapper.mast <- function(count_matrix, condition,
                               control = list(save_modelFit = FALSE,
                                              include_cdr = TRUE)) {
  
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(count_matrix))
  assertthat::assert_that(is.integer(count_matrix), 
                          msg = "count_matrix is not integer-values")
  assertthat::assert_that(dim(count_matrix)[2] == length(condition))
  
  # convert condition to factor
  if (!is.factor(condition)) {condition <- factor(condition)}

  #<--------------------------------------
  
  # compute log2CPM
  v <- limma::voom(count_matrix, design = model.matrix(~condition), plot=FALSE)
  log2CPM <- v$E
  
  # make data.frame into singleCellAssay object 
  colData <- data.frame(condition = condition)
  rowData <- data.frame(gene = rownames(count_matrix))
  sca <- MAST::FromMatrix(count_matrix, colData, rowData)
  
  if (control$include_cdr) {
    # calculate cellualr detection rate; normalized to mean 0 and sd 1
    colData(sca)$cdr.normed <- scale(colSums(assay(sca) > 0))
    # the default method for fitting is bayesGLM
    fit <- zlm.SingleCellAssay(~ condition + cdr.normed, sca)
  } else {
    fit <- zlm.SingleCellAssay(~ condition, sca)
  }
  
  # LRT test for the significance of the condition effect
  lrt <- lrTest(fit, "condition")
  
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






#' @title SCDE
#' 
#' @description Implement SCDE 1.99.1. This older version of SCDE is required for systems using
#'               an older version of flexmix (flexmix_2.3-13) (such as mine). (For related 
#'               discussions: https://github.com/hms-dbmi/scde/issues/40) 
#' 
#' @param count_matrix Gene by sample expression count matrix (G by N).
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param control A list with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete DESeq2 fit results.
#'   \code{filter_by_scde} SCDE internal argument. TRUE if apply SCDE built-in filtering criteria.
#'   \code{bootstrap_times} SCDE internal argument. Default 100 bootstraps.
#'
#' @return A list with the following objects
#'    \code{sig_order} Z-socres of expression difference, ordered from the most 
#'                     significant to the least significant.
#'    \code{fit} SCDE complete output of the model fit.

#' @author Chiaowen Joyce Hsiao
#'
#' @export
methodWrapper.scde <- function(count_matrix, condition,
                               control = list(save_modelFit = FALSE,
                                              n_cores = 1,
                                              filter_by_scde = TRUE,
                                              bootstrap_times = 100)) {
  
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(count_matrix))
  assertthat::assert_that(is.integer(count_matrix), 
                          msg = "count_matrix is not integer-values")
  assertthat::assert_that(dim(count_matrix)[2] == length(condition))
  
  # convert condition to factor
  if (!is.factor(condition)) {condition <- factor(condition)}
  
  #<--------------------------------------
  # apply SCDE default filtering criteria
  if (control$filter_by_scde == TRUE) {
    cd <- clean.counts(count_matrix, min.lib.size=1000, min.reads = 1, min.detected = 1)
  } else {
    cd <- count_matrix
  }
  
  # fitting error models, this step takes a considerable amount of time...
  o.ifm <- scde::scde.error.models(counts = cd, 
                                   groups = condition, 
                                   n.cores = control$n_cores, 
                                   threshold.segmentation = TRUE, 
                                   save.crossfit.plots = FALSE, 
                                   save.model.plots = FALSE, verbose = 1)
  
  # filter out cells that dont't show positive correlation with
  # the expected expression magnitudes 
  o.ifm <- o.ifm[o.ifm$corr.a > 0, ]
  
  # estimate gene expression prior
  o.prior <- scde::scde.expression.prior(models = o.ifm, 
                                         counts = cd, 
                                         length.out = 400, show.plot = FALSE)
  
  # DE analysis
  fit <- scde::scde.expression.difference(o.ifm, cd, o.prior, 
                                          groups  =  condition, 
                                          n.randomizations  =  control$bootstrap_times, 
                                          # number of bootstrap randomizations to be performed
                                          n.cores  =  control$n_cores, verbose  =  1)
  
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





#' @title ROTS
#' 
#' @describeIn Implement ROTS 1.2.0.
#' 
#' @param count_matrix Gene by sample expression count matrix (G by N).
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param control A list with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete DESeq2 fit results.
#'   \code{bootstrap_times} SCDE internal argument. Default 100 bootstraps.
#'
#' @return A list with the following objects
#'    \code{sig_order} Z-socres of expression difference, ordered from the most 
#'                     significant to the least significant.
#'    \code{fit} SCDE complete output of the model fit.

#' @author Chiaowen Joyce Hsiao
#'
#' @export
methodWrapper.rots <- function(count_matrix, condition,
                               control = list(save_modelFit = FALSE,
                                              bootstrap_times = 100)) {
  
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(count_matrix))
  assertthat::assert_that(is.integer(count_matrix), 
                          msg = "count_matrix is not integer-values")
  assertthat::assert_that(dim(count_matrix)[2] == length(condition))
  
  # convert condition to factor
  if (!is.factor(condition)) {condition <- factor(condition)}
  
  #<--------------------------------------
  # fitting model
  fit <- ROTS::ROTS(data = count_matrix,
                    groups = condition,
                    B = control$bootstrap_times,
                    K = dim(count_matrix)[1],
                    seed = 1234,
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




#' @title BPSC
#' 
#' @description Implement BPSC 0.99.1 (not available via Bioconductor). This method
#'              requires that the input data is normalized (counts-per-million or
#'              fragments per-kilo-base per million).
#' 
#' @param count_matrix Gene by sample expression count matrix (G by N).
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param control A list with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete DESeq2 fit results.
#'   \code{estIntPar} BPSC internal argument (TRUE/FALSE) for using only the expressed
#'                    genes to compute parameter estimates.
#'
#' @return A list with the following objects
#'    \code{pvalue} significance value for the group effect.
#'    \code{fit} BPSC complete output of the model fit.

#' @author Chiaowen Joyce Hsiao
#'
#' @export
methodWrapper.bpsc <- function(count_matrix, condition,
                               control = list(save_modelFit = FALSE,
                                              bootstrap_times = 100,
                                              option_estIntPar = FALSE)) {
  
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(count_matrix))
  assertthat::assert_that(is.integer(count_matrix), 
                          msg = "count_matrix is not integer-values")
  assertthat::assert_that(dim(count_matrix)[2] == length(condition))
  
  # convert condition to factor
  if (!is.factor(condition)) {condition <- factor(condition)}
  
  #<--------------------------------------
  # First, choose IDs of all cells of the control group for estimating
  # parameters of BP models
  controlIds <- which(as.integer(condition) == 1)
  
  # Create a design matrix including hte group labels, batch effectd
  # can be added here if they are available
  design <- model.matrix(~condition)
  
  # Model fitting; estIntPar requires to use only the expressed genes
  # to compute estiamtes, according to the vignette, this option requires longer
  # computational time
  fit <- BPSC::BPglm(data = count_matrix,
                     controlIds = controlIds,
                     design = design,
                     coef = 2, 
                     estIntPar = control$option_estIntPar, useParallel = FALSE)

  # extract results
  pvalue <- fit$PVAL
  
  # if save_modelFit, then output will include the original model fit
  if (control$save_modelFit) {
    fit <- fit
  } else {
    fit <- NULL
  }
  
  return(list(pvalue = pvalue,
              fit=fit))
}




