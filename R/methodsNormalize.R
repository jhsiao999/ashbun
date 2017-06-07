## This document contains for the various methods for normalizing
## single cell expression data


#' @title TMM normalization
#'
#' @param object Gene by sample expression count matrix.
#'
#' @return an edgeR \code{DGE} object containing the raw count expression matrix, the raw library size (\code{lib.size}) and the library size normalizing factor (\code{norm.factors}).
#'
#' @export
normalize.TMM <- function(object) {
    object.integer <- as.matrix(as.integer(object), dim = dim(object))
    dgecounts <- edgeR::calcNormFactors(DGEList(object.integer), method = "TMM")
    return(dgecounts)
}
