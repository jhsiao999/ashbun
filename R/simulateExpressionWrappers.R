###########################################################
#  Code for simulating single cell expression count data
#
#  Wrapper for functions in simulateExpression.R
###########################################################


#' @title Wrapper for simulating datasets
#'
#' @param counts gene by sample count matrix
#' @param Nsim number of simulated datasets
#' @param Nsample number of samples per biological condition
#' @param Ngene number of genes. Defaults to include all genes in the input data.
#' 
#' @examples 
#' ipsc_eset <- get(load(system.file("testdata", "HumanTungiPSC.rda", package = "ashbun")))
#' counts <- exprs(ipsc_eset)[sample(nrow(exprs(ipsc_eset)), ), ]
#' 
#' simdata_list <- simulationWrapper(counts, Nsim = 5, Nsample = 80, Ngene = 500)
#' 
#' @return List of data generated under three different fractions of null genes.
#'   \code{null} pi0 = 0. List of Nsim simulated datasets.
#'   \code{normal_5} pi0 = 0.5. List of Nsim simulated datasets. For each simulated dataset, I store the count table and a logical vector indicating TRUE = null gene, and FALSE = true DE gene.
#'   \code{norma_9} pi0 = 0.9. List of Nsim simulated datasets. List of Nsim simulated datasets. For each simulated dataset, I store the count table and a logical vector indicating TRUE = null gene, and FALSE = true DE gene.
#' 
#' @export
simulationWrapper <- function(counts,
                              Nsim, Nsample, Ngene = NULL) {
  # null dataset
  counts_null <- lapply(1:Nsim, function(i) {
    foo <- ashbun::makeSimCount2groups(counts = counts,
                              Nsam = Nsample, Ngene = Ngene,
                              sample_method = "all_genes",
                              control = list(pi0 = NULL))
    condition <- rep(1:2, each = Nsample)
    condition <- condition[sample(1:(Nsample*2), size = Nsample*2)]

    geneToInclude <- which(rowSums(foo$counts) != 0)
    sampleToInclude <- which(colSums(foo$counts) != 0)

    return(list(counts = foo$counts[geneToInclude, sampleToInclude],
                condition = condition[sampleToInclude]))
  })

  counts_bignormal_5 <- lapply(1:Nsim, function(i) {
  #  set.seed(999*i)
    foo <- ashbun::makeSimCount2groups(counts = counts,
                               Nsamp = Nsample, Ngene = Ngene,
                               sample_method = "all_genes")
    foo2 <- ashbun::non_null_sim(counts = foo$counts,
                         args = args.big_normal(nsam = Nsample, pi0 = .5))

    geneToInclude <- which(rowSums(foo2$counts) != 0)
    sampleToInclude <- which(colSums(foo2$counts) != 0)
    foo2$counts <- foo2$counts[geneToInclude, sampleToInclude]
    foo2$condition <- foo2$condition[sampleToInclude]
    foo2$null <- foo2$null[geneToInclude]
    return(foo2)
  })

  counts_bignormal_9 <- lapply(1:Nsim, function(i) {
  #  set.seed(999*i)
    foo <- ashbun::makeSimCount2groups(counts = counts,
                               Nsamp = Nsample, Ngene = Ngene,
                               sample_method = "all_genes")
    foo2 <- ashbun::non_null_sim(counts = foo$counts,
                         args = args.big_normal(nsam = Nsample, pi0 = .9))

    geneToInclude <- which(rowSums(foo2$counts) != 0)
    sampleToInclude <- which(colSums(foo2$counts) != 0)
    foo2$counts <- foo2$counts[geneToInclude, sampleToInclude]
    foo2$condition <- foo2$condition[sampleToInclude]
    foo2$null <- foo2$null[geneToInclude]
    return(foo2)
  })

  return(list(null = counts_null,
              normal_5 = counts_bignormal_5,
              normal_9 = counts_bignormal_9))
}
