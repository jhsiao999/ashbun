# This document contains code for implementing DE methods originally developed for bulk RNA-seq data

#' @title DESeq2
#'
#' @description Implement DESeq2.
#'
#' @param counts Gene by sample expression count matrix (G by N).
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param libsize_factors Numeric vector of scale factors for library sizes.
#'                       Default NULL computes libsize_factors using DESeq2.
#' @param control A list with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete DESeq2 fit results.
#'   \code{independentFiltering} DESeq2 internal argument. Default = FALSE to not apply gene filtering.
#'   \code{cooksCutoff} DESeq2 internal argument. Default = FALSE to not exclude samples.
#'
#' @return List with the following objects
#'   \code{betahat} Estimate effect size of condition for all genes.
#'   \code{sebetahat} Standard errors of the effect sizes.
#'   \code{pvalue} P-values of the effect sizes.
#'   \code{fit} DESeq2 complete output of the model fit. NULL if \code{save_modelFit = FALSE}.
#'
#' @author Chiaowen Joyce Hsiao
#'
#' @export
methodWrapper.DESeq2 <- function(counts, condition, libsize_factors = NULL,
                                 default = FALSE,
                                 control = list(save_modelFit = FALSE,
                                                independentFiltering = FALSE,
                                                cooksCutoff = FALSE)){

  suppressPackageStartupMessages(library(DESeq2))
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(counts))
  assertthat::assert_that(dim(counts)[2] == length(condition))

  # this gives errors... dont' know why
  # assertthat::assert_that(length(unique(condition)) != 2,
  #                         msg = "condition vector only allows 2 groups")

  # convert condition to factor
  if (!is.factor(condition)) {condition <- factor(condition)}

  #<--------------------------------------
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, S4Vectors::DataFrame(condition), ~ condition)

  # input pre-defined size factor    
  if (default == FALSE) {
    dds <- DESeq2::estimateSizeFactors(dds)
    names(libsize_factors) <- colnames(counts)
    dds@colData@listData$sizeFactor <- libsize_factors
  }
  
  # use DESeq size factor
  if (default == TRUE) {
    dds <- DESeq2::estimateSizeFactors(dds)
  }
  
  # Run DE analysis
  dds <- DESeq(dds, quiet = TRUE)

  # Call results table without any arguments
  # this will extract the estimated log2 fold changes and p values for the
  # last variable in the design formula
  res <- results(dds,
                 cooksCutoff = control$cooksCutoff,
                 independentFiltering = control$independentFiltering)

  betahat <- res$log2FoldChange
  sebetahat <- res$lfcSE
  df <- nrow(colData(dds)) - ncol(model.matrix(design(dds)))
  pvalue <- res$pvalue
  names(pvalue) <- rownames(dds)

  # if save_modelFit, then output will include the original model fit
  if (control$save_modelFit) {
    fit <- res
  } else {
    fit <- NULL
  }

  return(list(betahat=betahat, sebetahat=sebetahat,
              df=df,
              pvalue = pvalue, fit = fit))
}





#' @title edgeR
#'
#' @description Implement edgeR.
#'
#' @param counts Gene by sample expression count matrix (G by N).
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param libsize_factors Numeric vector of scale factors for library sizes.
#'                       Default NULL computes libsize_factors using edgeR.
#' @param control List with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete DESeq2 fit results.
#'
#' @return List with the following objects
#'   \code{betahat} Estimate effect size of condition for all genes.
#'   \code{pvalue} P-values of the effect sizes.
#'   \code{fit} edgeR complete output of the model fit. NULL if \code{save_modelFit = FALSE}.
#'
#' @author Chiaowen Joyce Hsiao
#'
#' @export
methodWrapper.edgeR <- function(counts, condition, libsize_factors = NULL,
                                default = FALSE,
                                control = list(save_modelFit = FALSE)) {

  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(dim(counts)[2] == length(condition))

  counts <- as.matrix(counts)
  if (!is.factor(condition)) {condition <- factor(condition)}

  #<--------------------------------------
  # Make "DGEList" object
  if (default == FALSE) {
    dge <- edgeR::DGEList(counts = counts,
                          group = condition,
                          genes = rownames(counts),
                          norm.factors = libsize_factors)
  }
  if (default == TRUE) {
    dge <- edgeR::DGEList(counts = counts,
                          group = condition,
                          genes = rownames(counts))
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
  }

  # estimate dispersion
  dge <- edgeR::estimateDisp(dge, design = model.matrix(~condition), robust = TRUE)

  # Run DE analysis; dispersion = NULL will extract tagwise (genewise) dispersion estimates
  # for DE analysis
  fit <- edgeR::glmFit(dge, dispersion = NULL)

  # Run LRT test
  lrt <- edgeR::glmLRT(fit, coef = 2)

  betahat <- lrt$coefficients[,2]
  df <- lrt$df.residual
  pvalue <- lrt$table$PValue
  names(pvalue) <- rownames(fit)
  
  # if save_modelFit, then output will include the original model fit
  if (control$save_modelFit) {
    fit <- lrt
  } else {
    fit <- NULL
  }

  return(list(betahat=betahat,
              sebetahat=rep(NA, length(betahat)),
              df=df,
              pvalue = pvalue, fit = fit))
}




#' @title limma + voom
#'
#' @description Implement limma with voom and apply the empirical Bayes method
#'    in limma.
#'
#' @param counts normalized gene by sample expression count matrix. (CPM)
#' @param condition binary vector of length N indicating sample biological condition.
#' @param libsize_factors Numeric vector of scale factors for library sizes.
#'                       Default NULL multiplies all library sizes by 1.
#' @param pseudocount default .5.
#' @param control List with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete limmaVoom fit results.
#'
#' @return List with the following objects
#'    \code{betahat} Estimate effect size of condition for all genes.
#'    \code{sebetahat} Standard errors of the effect sizes.
#'    \code{df} Degrees of freedom associated with the effect sizes.
#'    \code{pvalue} P-values of the effect sizes.
#'    \code{fit} limmaVoom complete output of the model fit.
#' @author Chiaowen Joyce Hsiao
#'
#' @export
methodWrapper.limmaVoom <- function(counts, condition, pseudocount = .5,
                                    default = FALSE,
                                    control = list(save_modelFit = FALSE)){

  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(dim(counts)[2] == length(condition))

  counts <- as.matrix(counts)
  if (!is.factor(condition)) {condition <- factor(condition)}

  # convert condition to factor
  if (!is.factor(condition)) {condition <- factor(condition)}

  #<--------------------------------------
  # create design matrix
  design <- model.matrix(~condition)

  suppressPackageStartupMessages(library(limma))

  if (default==FALSE) {
    # compute log2CPM
    if (is.null(pseudocount)) {
      log2CPM <- log2(counts)
    } else {
      log2CPM <- log2(counts + pseudocount)
    }
    
    # don't apply normalization methods (eg., TMM, quantile)
    # the default setting of voom.controlPseudocount is the same as
    # in voom package; but here we made pseudocount and pseudo library size
    # adjustable by users
    weights <- voom.controlPseudocount(counts, design)
    fit <- limma::lmFit(log2CPM, design, weights = weights)
    fit.ebayes <- limma::eBayes(fit)
  }
  
  if (default==TRUE) {
    v <- voom(counts, design, plot = FALSE)
    fit <- lmFit(v, design)
    fit.ebayes <- limma::eBayes(fit)
  }
  
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





#' @title limma + voom
#'
#' @description Implement limma with voom and apply the empirical Bayes method
#'    in limma. This is a modified version of voom which allows to specify pseudocount and
#'    pseudo library size.
#'    Hence, voom: y <- t(log2(t(counts + pseudocount)/(lib.size + pseudo_libsizes) * 1e+06)).
#'
#'
#' @param counts gene by sample expression count matrix.
#' @param design design matrix, generated by R function model.matrix()
#' @param pseudocount default .5
#' @param pseudo_libsizes default 1.
#'
#' @return
#'    \code{w} Weights of dimension G by N.
#' @author Chiaowen Joyce Hsiao
#'
#' @export
voom.controlPseudocount <- function(counts, design,
                                    pseudocount = .5,
                                    pseudo_libsizes = 1,
                                    span = .5) {
  # this function wad adpated from the limma voom function
  # prior count and library size adjustment are set to be argument
  lib.size <- colSums(counts)

  y <- t(log2(t(counts + pseudocount)/(lib.size + pseudo_libsizes) * 1e+06))
  fit <- lmFit(y, design)
  if (is.null(fit$Amean))  fit$Amean <- rowMeans(y, na.rm = TRUE)

  # compute variance-mean dependency
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  f <- approxfun(l, rule = 2)
  if (fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[,
                                                                  j, drop = FALSE])
  }
  else {
    fitted.values <- fit$coef %*% t(fit$design)
  }
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + pseudo_libsizes))
  fitted.logcount <- log2(fitted.count)
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)

  return(w)
}
