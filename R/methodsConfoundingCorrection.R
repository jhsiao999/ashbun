## This document contains the code for confounding correction


#' @title Estimate number of surrogate variables
#'
#' @param object gene by sample expression matrix. Can be
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param num_sv number of surrogate variables
#'
#' @export
methodWrapper.num_sv <- function(object, condition) {
  Y0 <- log2(object + 1)
  X0 <- model.matrix(~ condition)
  sva::num.sv(dat = Y0, mod = X0, method = "be")
}


#' @title Surrogate Variable Analysis
#'
#' @param count_matrix gene by sample expression matrix. Can be
#' @param design_matrix Binary vector of length N indicating sample biological condition.
#' @param num_sv number of surrogate variables
#'
#' @export
methodWrapper.sva <- function(count_matrix, design_matrix, num_sv){

  Y <- count_matrix
  X <- design_matrix

  trash     <- capture.output(sva_out <- sva::sva(dat = Y, mod = X, n.sv = num_sv))
  X.sv      <- cbind(X, sva_out$sv)
  limma_out <- limma::lmFit(object = Y, design = X.sv)
  limma.ebayes <- limma::eBayes(limma_out)
  betahat <- limma.ebayes$coefficients[,2]
  sebetahat <- with(limma.ebayes, stdev.unscaled[,2]*sigma)
  p.value <- limma.ebayes$p.value
  df.total <- limma.ebayes$df.total
  # betahat   <- limma_out$coefficients[, 2]
  # sebetahat <- limma_out$stdev.unscaled[, 2] * limma_out$sigma
  # df        <- limma_out$df.residual[1]
  # tstats    <- betahat / sebetahat
  # pvalues   <- 2 * pt(-abs(tstats), df = df)

  return(list(betahat = betahat, sebetahat = sebetahat, df = df.total,
              pvalues = p.value))
}
