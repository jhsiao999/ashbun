# This document contains code for implementing DE methods originally developed for bulk RNA-seq data

#' @title limma + voom
#' 
#' @description Implement limma with voom and apply the empirical Bayes method 
#'    in limma. 
#'
#' @param count_matrix Gene by sample expression count matrix (G by N).
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param control A list with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete limmaVoom fit results.
#'
#' @return A list with the following objects
#'    \code{betahat} The estimate effect size of condition for all genes.
#'    \code{sebetahat} The standard errors of the effect sizes.
#'    \code{df} The degrees of freedom associated with the effect sizes.
#'    \code{pvalue} P-values of the effect sizes.
#'    \code{fit} limmaVoom complete output of the model fit.
#'
#' @author Chiaowen Joyce Hsiao
#'
#' @export
methodWrapper.limmaVoom <- function(count_matrix, condition,
                                    control = list(save_modelFit = FALSE)){

  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(count_matrix))
  assertthat::assert_that(is.integer(count_matrix), 
                          msg = "count_matrix is not integer-values")
  assertthat::assert_that(dim(count_matrix)[2] == length(condition))
  
  # This gives erros... don't know why...
  # assertthat::assert_that(length(unique(condition)) != 2,
  #                         msg = "condition vector only allows 2 groups")
  
  # convert condition to factor
  if (!is.factor(condition)) {condition <- factor(condition)}
  
  #<--------------------------------------
  # create design matrix
  design <- model.matrix(~condition)
  
  # apply voom: outputs log2CPM values, 
  # does not apply normalization methods (eg., TMM, quantile)
  v <- limma::voom(count_matrix, design, plot=FALSE)
  fit <- limma::lmFit(v)
  fit.ebayes <- limma::eBayes(fit)
  
  # given that the condition is a binary vector
  # extract the coefficient corresponds to the difference between the two conditions
  betahat <- fit.ebayes$coefficients[,2]
  sebetahat <- with(fit.ebayes, stdev.unscaled[,2]*sigma)
  pvalue <- fit.ebayes$p.value[,2]
  df <- fit.ebayes$df.total

  # if save_modelFit, then output will include the original model fit
  if (control$save_modelFit) {
    fit <- fit.ebayes
  } else {
    fit <- NULL
  }
  
  return(list(betahat=betahat, sebetahat=sebetahat, 
              df=df, pvalue = pvalue, fit = fit))
}



#' @title DESeq2
#' 
#' @description Implement DESeq2.
#'
#' @param count_matrix Gene by sample expression count matrix (G by N).
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param control A list with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete DESeq2 fit results.
#'   \code{independentFiltering} DESeq2 internal argument. Default = FALSE to not apply gene filtering.
#'   \code{cooksCutoff} DESeq2 internal argument. Default = FALSE to not exclude samples.
#'
#' @return A list with the following objects
#'   \code{betahat} The estimate effect size of condition for all genes.
#'   \code{sebetahat} The standard errors of the effect sizes.
#'   \code{pvalue} P-values of the effect sizes.
#'   \code{fit} DESeq2 complete output of the model fit. NULL if \code{save_modelFit = FALSE}.
#'
#' @author Chiaowen Joyce Hsiao
#'
#' @export
methodWrapper.DESeq2 <- function(count_matrix, condition,
                                 control = list(save_modelFit = FALSE,
                                                independentFiltering = FALSE,
                                                cooksCutoff = FALSE)){
  
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(count_matrix))
  assertthat::assert_that(is.integer(count_matrix), 
                          msg = "count_matrix is not integer-values")
  assertthat::assert_that(dim(count_matrix)[2] == length(condition))
  
  # this gives errors... dont' know why
  # assertthat::assert_that(length(unique(condition)) != 2,
  #                         msg = "condition vector only allows 2 groups")
  
  # convert condition to factor
  if (!is.factor(condition)) {condition <- factor(condition)}
  
  #<--------------------------------------
  # Make "DESeqDataSet" object
  dds <- DESeq2::DESeqDataSetFromMatrix(count_matrix, DataFrame(condition), ~ condition)
  
  # Run DE analysis
  dds <- DESeq(dds)
  
  # Call results table without any arguments
  # this will extract the estimated log2 fold changes and p values for the
  # last variable in the design formula
  res <- results(dds, 
                 cooksCutoff = control$cooksCutoff,
                 independentFiltering = control$independentFiltering)
  
  betahat <- res$log2FoldChange
  sebetahat <- res$lfcSE
  pvalue <- res$pvalue
  
  # if save_modelFit, then output will include the original model fit
  if (control$save_modelFit) {
    fit <- res
  } else {
    fit <- NULL
  }
  
  return(list(betahat=betahat, sebetahat=sebetahat, 
              pvalue = pvalue, fit = fit))
}





#' @title edgeR
#' 
#' @description Implement edgeR.
#'
#' @param count_matrix Gene by sample expression count matrix (G by N).
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param control A list with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete DESeq2 fit results.
#'   \code{libsize_factors} Numeric vector of library size normalization factor.
#'
#' @return A list with the following objects
#'   \code{betahat} The estimate effect size of condition for all genes.
#'   \code{pvalue} P-values of the effect sizes.
#'   \code{fit} edgeR complete output of the model fit. NULL if \code{save_modelFit = FALSE}.
#'
#' @author Chiaowen Joyce Hsiao
#'
#' @export
methodWrapper.edgeR <- function(count_matrix, condition,
                                 control = list(save_modelFit = FALSE,
                                                libsize_factors = rep(1:ncol(count_matrix)))){
  
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(count_matrix))
  assertthat::assert_that(is.integer(count_matrix), 
                          msg = "count_matrix is not integer-values")
  assertthat::assert_that(dim(count_matrix)[2] == length(condition))
  
  # This gives erros... don't know why...
  # assertthat::assert_that(length(unique(condition)) != 2,
  #                         msg = "condition vector only allows 2 groups")
  
  # convert condition to factor
  if (!is.factor(condition)) {condition <- factor(condition)}
  
  #<--------------------------------------
  # Make "DGEList" object
  dge <- edgeR::DGEList(counts = count_matrix, 
                        group = condition,
                        genes = rownames(count_matrix),
                        norm.factors = libsize_factors)
  
  # estimate dispersion
  dge <- estimateDisp(dge, design = model.matrix(~condition), robust = TRUE)
  
  # Run DE analysis; dispersion = NULL will extract tagwise (genewise) dispersion estimates
  # for DE analysis
  fit <- edgeR::glmFit(dge, dispersion = NULL)
  
  # Run LRT test
  lrt <- edgeR::glmLRT(fit, coef = 2)

  betahat <- lrt$coefficients[,2]
  pvalue <- lrt$table$PValue
  
  # if save_modelFit, then output will include the original model fit
  if (control$save_modelFit) {
    fit <- lrt
  } else {
    fit <- NULL
  }
  
  return(list(betahat=betahat, 
              pvalue = pvalue, fit = fit))
}




