## This document contains the code for implementing the various existing methods for single cell differential abundance analysis



#' @title voom + limma
#'
#' @param object gene by sample expression matrix. Can be
#' @param log2counts G genes by N samples log2count matrix.
#' @param condition Binary vector of length N indicating sample biological condition.
#'
#' @example
#' ## objectNormed <- normalize.TMM(zeisel_cortex_subset)
#' @export
methodWrapper.limmaVoom <- function(object, condition, W=NULL){

  if (is.null(W)){
    design <- model.matrix(~condition)
  }else{
    design <- model.matrix(~condition+W)
  }

  v <- limma::voom(object, design, plot=FALSE)
  fit <- limma::lmFit(v)
  fit.ebayes <- limma::eBayes(fit)

  betahat <- fit.ebayes$coefficients[,2]
  sebetahat <- with(fit.ebayes, stdev.unscaled[,2]*sigma)
  p.value <- fit.ebayes$p.value
  df.total <- fit.ebayes$df.total

  return(list(betahat=betahat, sebetahat=sebetahat, df=df.total, p.value = p.value,
  fit.model = fit.ebayes))
}




#' @title MAST
#'
#' @param object gene by sample expression matrix containing log2 adjusted counts.
#' @param condition vector specifying sample treatment condition.
#' @param cdr cellular detection rate, computed as \code{scale(colSums(counts > 0))} where counts is raw expression count matrix.
#'
#' @return p.value significance value from likelihood ratio test.
#' @export
methodWrapper.mast <- function(object, condition, cdr) {

  # Check column names of log2counts;
  # if repeated, then change to unique
  if (sum(duplicated(colnames(object))) > 0) {
    colnames(object) <- paste0("sample.", c(1:ncol(object)))
  }

  # make data.frame into singleCellAssay object for MAST computing
  library(MAST)
  colData <- data.frame(condition = condition)
  rowData <- data.frame(gene = rownames(object))
  sca <- FromMatrix(object, colData, rowData)

  if (!is.null(cdr)) {
    # calculate cellualr detection rate
    # normalized to mean 0 and sd 1
    colData(sca)$cdr.normed <- cdr
    # the default method for fitting is bayesGLM
    fit <- zlm.SingleCellAssay(~ condition + cdr.normed, sca)
  }
  if (is.null(cdr)) {
    fit <- zlm.SingleCellAssay(~ condition, sca)
  }

  # LRT test for the significance of the condition effect
  lrt <- lrTest(fit, "condition")

  # extract p.value from their hurdle model
  p.value <- lrt[,3,3]

    # extract effect size, standard error, and df from the non-zero component
  betahat <- fit@coefC[,2]
  setbetahat <- sqrt(sapply(1:dim(fit@vcovC)[3], function(i) {
    diag(fit@vcovC[,,i])[2]} ) )
  df <- fit@df.resid[,1]

  return(list(betahat=betahat,
              sebetahat=setbetahat,
              df=df,
              p.value = p.value,
              fit.model = fit))
}




#methodWrapper.deseq2


#methodWrapper.deseq


#methodWrapper.scde
