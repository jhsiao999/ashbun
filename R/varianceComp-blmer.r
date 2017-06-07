## The code here was used to perform
## variance component analysis in the 2017 single-cell batch paper
## from the Gilad lab DOI: 10.1038/srep39921


#' Per gene variance component model
#'
#' @description Wrapper of the blmer function. Fits a bayesian nested model for one gene at a time.
#'
#' @param xx Matrix of expression measurements on log scale.
#' @param annotation Meta-data matrix of each column of xx.
#'
#' @examples
#' molecules_final <- get(load("molecules-final.txt"))
#' anno_filter <- get(load("annotation-filter.txt"))
#' blme_final <- gene_variation(counts = molecules_final,
#'                              annotation = anno_filter)
#' @export
gene_variation <- function(counts, annotation) {

  individual <- as.factor(annotation$individual)
  replicate <- as.factor(annotation$replicate)

  ## fit bayesian GLM one gene at a time

  blme_fit <- lapply( 1:NROW(counts), function(i) {

      value <- unlist(counts[i,])

      fit_try <- tryCatch(
        fit <- blme::blmer(value ~ 1|individual/replicate,
                           cov.prior = gamma(shape = 2),
                           resid.prior = gamma(shape = 2)),
                             condition = function(c) c)
      if(inherits(fit_try, "condition")){
      var_foo <- rep(NA, 3)
      return(var_foo)
      }
      if(!inherits(fit_try, "condition")){
        var_foo <- as.data.frame(VarCorr(fit_try))[,4]
        var_foo <- var_foo[c(2,1,3)]
        var_foo
      }
    })
  blme_fit <- do.call(rbind, blme_fit)
  rownames(blme_fit) <- rownames(counts)
  colnames(blme_fit) <- c("individual","replicate","residual")
  blme_fit
}


###<--- Compute the proportion of variance explained
# the above analysis produces variance component estimates (e.g., $\sigma^2_b$ for batch effect) that are based on a penalized maximum likelihood approach. We compute naive approximation of sum of squared variation for the individual effect and for the batch effect, and their proportions of variation. Specifically, to simplify the computations of degrees of freedom for each factor, we approximate a balanced nested design and compute estiamted number of levels of each factor as the average of the observed number of levels of each factor: the approximate number of batches is 2.67 (i.e., (2+3+3)/3) and and the approximate number of cell is 70.5 (i.e., average number of cell samples per batch).


##load("../data/blme-variance.rda")
#res <- blme_final
#ms_ind <- (res[,1]*2.67*70.5) + (res[,2]*70.5) + res[,3]
#ms_batch <- (res[,2]*70.5) + res[,3]
#ms_resid <- res[,3]
#ss_ind <- ms_ind*(3-1)
#ss_batch <- ms_batch*3*(2.67-1)
#ss_resid <- ms_resid*3*2.67*(70.5-1)
## individual component
#prop_ind <- ss_ind/(ss_ind + ss_batch + ss_resid)
## batch component
#prop_batch <- ss_batch/(ss_ind + ss_batch + ss_resid)
