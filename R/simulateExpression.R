###############################################################
#  Simulate single cell expression count data the Stephens way..
#
#  Date:2017-07-26
#  Adapted from Mengyin Lu's work: Variable and function names
#  may be different from the original code.
#
#  The main top-level functions are
#  makeSimCount2groups(...): generate null gene expression set
#     from experimental data,
#  non_null_sim(...): sample a random set of null expression data
#     to be signal gene and determine signal based on a normal
#     distribution
###############################################################

#' @title Wrapper for simulating M datasets
#'
#' @description This different version of simulationWrapper first selects
#'              single cell samples, performs filtering, selects a subset of genes,
#'              and then performs  permutation at the gene-level or sample-level.
#'              This is in contrast with the previous approach, where filtering is performed
#'              after select a subset of genes and permuting labels.
#' @param counts gene by sample count matrix
#' @param Nsim number of simulated datasets
#' @param Nsample number of samples per biological condition
#' @param Ngenes number of genes. Defaults to include all genes in the input data.
#'
#' @examples
#' ipsc_eset <- get(load(system.file("testdata", "HumanTungiPSC.rda", package = "ashbun")))
#' counts <- exprs(ipsc_eset)[sample(nrow(exprs(ipsc_eset)), ), ]
#'
#' simdata_list <- simulationWrapper(counts, Nsim = 5, Nsample = 80, Ngenes = 500)
#'
#' @return List of data generated under three different fractions of null genes.
#'   \code{null} pi0 = 0. List of Nsim simulated datasets.
#'   \code{normal_5} pi0 = 0.5. List of Nsim simulated datasets. For each simulated dataset, I store the count table and a logical vector indicating TRUE = null gene, and FALSE = true DE gene.
#'   \code{norma_9} pi0 = 0.9. List of Nsim simulated datasets. List of Nsim simulated datasets. For each simulated dataset, I store the count table and a logical vector indicating TRUE = null gene, and FALSE = true DE gene.
#'
#' @examples
#' library(singleCellRNASeqHumanTungiPSC)
#' eset <- HumanTungiPSC
#' counts <- exprs(eset)[,pData(eset)$individual == "NA19101"]
#'
#' sim_data_null <- simulationWrapper(counts,
#'                               Ngenes = 100,
#'                               Nsamples = 20,
#'                               sample_method = "all_genes",
#'                               pi0 = 1)
#'
#' sim_data_nonnull <- simulationWrapper(counts, Nsim = 2,
#'                               Ngenes = 100,
#'                               Nsam = 20,
#'                               sample_method = "all_genes",
#'                               pi0 = .5,
#'                               beta_args = args.big_normal(betapi = 1,
#'                                                           betamu = 0, betasd = .8))
#'
#' @export
simulationWrapper.filter <- function(counts,
                              Nsim = 1, Nsamples = 50, Ngenes = NULL,
                              pi0 = NULL,
                              sample_method = c("all_genes", "per_gene"),
                              samplesFractionExpressed=.25,
                              featuresFractionExpressed=.25,
                              thresholdDetection=1,
                              beta_args = args.big_normal(betapi = c(1),
                                                          betamu = c(0),
                                                          betasd = c(1))) {
  output <- list(counts = counts)

  Nsamples_total <- 2*Nsamples

  simdata <- lapply(1:Nsim, function(i) {
    # filter samples
    samplesToInclude.filter <-  filterSamples.fractionExpressed(counts,
                                    thresholdDetection=thresholdDetection,
                                    fractionExpressed=samplesFractionExpressed)$index_filter
    # subset samples
    counts_subset <- counts[, samplesToInclude.filter]

    # choose a random set of samples
    samplesToInclude.permute <- sample((1:NCOL(counts_subset)), Nsamples_total, replace = FALSE)
    counts_subset <- counts_subset[, samplesToInclude.permute]

    # permute sample labels
    if (sample_method == "per_gene") {
     #For each gene, randomly select 2*Nsamp samples from counts
      counts_perm <- t(apply(counts_subset, 1, sampleingene, Nsamples_total=Nsamples_total))
      }
    if (sample_method == "all_genes") {
      counts_perm <- counts_subset
    }

    # assign samples to 2 arbitrary conditions
    condition <- rep(1:2, each = Nsamples)
    condition <- condition[sample(1:(Nsamples_total), size = Nsamples_total)]
    output$condition <- condition

    # reorder sample columns
    output$condition <- output$condition[order(output$condition)]
    output$counts <- counts_perm[, order(output$condition)]

    # make colnames be unique
    colnames(output$counts) <- paste0("sample_",c(1:dim(output$counts)[2]))

    # filter genes
    featuresToInclude.filter <- filterFeatures.fractionExpressed(output$counts,
                              condition=output$condition,
                              thresholdDetection=thresholdDetection,
                              fractionExpressed=featuresFractionExpressed)$index_filter
    output$counts <- output$counts[featuresToInclude.filter, ]

    # select random subset of genes
    featuresToInclude.permute <- sample(1:NROW(output$counts), Ngenes, replace = FALSE)
    output$counts <- output$counts[featuresToInclude.permute, ]

    # add rownames if null
    # make pseudo rownames
    if (is.null(rownames(output$counts))) {
    rownames(output$counts) <-paste0("feature_", c(1:dim(output$counts)[1]))
    }

    #  set.seed(999*i)
    if (pi0 < 1) {
      data_signal <- non_null_sim(counts = output$counts,
                           condition = output$condition,
                           pi0,
                           beta_args = beta_args)
      return(data_signal)
    }
  })
}

#' @title Generate count matrix of all null genes
#'
#' @description This is an updated version of makeSimCount2groups which takes in filtered data
#'               and performs permutation at gene-level and at the sample level.
#' @param counts Gene expression count matrix from a dataset.
#' @param Ngenes Number of genes in the simulated dataset.
#' @param Nsample Number of samples in each condition.
#' @param pi0 Proportion of null genes. Default to be 1.
#'
#' @examples
#' library(singleCellRNASeqHumanTungiPSC)
#' eset <- HumanTungiPSC
#' counts <- exprs(eset)[,pData(eset)$individual == "NA19101"]
#'
#' sim_counts <- makeSimCount2groups(counts,
#'                                   Ngenes = 100,
#'                                   Nsample = 20,
#'                                   sample_method = "all_genes")
#'
#' @export
makeSimCount2groups.filter <- function(counts, Ngenes = NULL,
                                sample_method = c("per_gene", "all_genes")){
  output <- list(counts = counts)
  Nsamples_total <- ncol(counts)
  Nsamples <- ncol(counts)/2

  # select a subset of genes
  if (!is.null(Ngenes)) {
    genes_to_include <- sample(1:NROW(counts), Ngenes, replace = FALSE)
    counts_subset <- counts[genes_to_include, ]
  } else {
    counts_subset <- counts
  }

  if (sample_method == "per_gene") {
    # For each gene, randomly select 2*Nsamp samples from counts
    counts_subset <- t(apply(counts_subset, 1, sampleingene, Nsamples_total=Nsamples_total))
    output$counts <- counts_subset
#    # Remove genes without any reads
#    if ( sum(apply(counts_subset,1,sum) == 0) > 0) {
#      counts_subset <- counts_subset[which(apply(counts_subset,1,sum)>0), ]
#    } else {
#      counts_subset <- counts_subset
#    }
#    # #Ngenes <- NROW(counts_subset)
#    # output$counts <- counts_subset
  }

  if (sample_method == "all_genes") {
    output$counts <- counts_subset
#    samples_to_include <- sample((1:NCOL(counts)), Nsamples_total, replace = TRUE)
#    counts_subset <- counts[, samples_to_include]
#    # Remove genes without any reads
#    counts_subset <- counts_subset[which(apply(counts_subset,1,sum)>0), ]
#    # #Ngenes <- NROW(counts_subset)
#    # output$counts <- counts_subset
  }


  # assign samples to 2 arbitrary conditions
  condition <- rep(1:2, each = Nsamples)
  condition <- condition[sample(1:(Nsamples_total), size = Nsamples_total)]
  output$condition <- condition

  # reorder sample columns
  output$condition <- output$condition[order(output$condition)]
  output$counts <- output$counts[, order(output$condition)]

  # make colnames be unique
  colnames(output$counts) <- paste0("sample_",c(1:dim(output$counts)[2]))

  return(output)
}


#---------------------------------------------------------------

#' @title Wrapper for simulating M datasets
#'
#' @param counts gene by sample count matrix
#' @param Nsim number of simulated datasets
#' @param Nsample number of samples per biological condition
#' @param Ngenes number of genes. Defaults to include all genes in the input data.
#'
#' @examples
#' ipsc_eset <- get(load(system.file("testdata", "HumanTungiPSC.rda", package = "ashbun")))
#' counts <- exprs(ipsc_eset)[sample(nrow(exprs(ipsc_eset)), ), ]
#'
#' simdata_list <- simulationWrapper(counts, Nsim = 5, Nsample = 80, Ngenes = 500)
#'
#' @return List of data generated under three different fractions of null genes.
#'   \code{null} pi0 = 0. List of Nsim simulated datasets.
#'   \code{normal_5} pi0 = 0.5. List of Nsim simulated datasets. For each simulated dataset, I store the count table and a logical vector indicating TRUE = null gene, and FALSE = true DE gene.
#'   \code{norma_9} pi0 = 0.9. List of Nsim simulated datasets. List of Nsim simulated datasets. For each simulated dataset, I store the count table and a logical vector indicating TRUE = null gene, and FALSE = true DE gene.
#'
#' @examples
#' library(singleCellRNASeqHumanTungiPSC)
#' eset <- HumanTungiPSC
#' counts <- exprs(eset)[,pData(eset)$individual == "NA19101"]
#'
#' sim_data_null <- simulationWrapper(counts,
#'                               Ngenes = 100,
#'                               Nsamples = 20,
#'                               sample_method = "all_genes",
#'                               pi0 = 1)
#'
#' sim_data_nonnull <- simulationWrapper(counts, Nsim = 2,
#'                               Ngenes = 100,
#'                               Nsam = 20,
#'                               sample_method = "all_genes",
#'                               pi0 = .5,
#'                               beta_args = args.big_normal(betapi = 1,
#'                                                           betamu = 0, betasd = .8))
#'
#' @export
simulationWrapper <- function(counts,
                              Nsim = 1, Nsamples = 50, Ngenes = NULL,
                              pi0 = NULL,
                              sample_method = c("all_genes", "per_gene"),
                              beta_args = args.big_normal(betapi = c(1),
                                                          betamu = c(0),
                                                          betasd = c(1))) {
  simdata <- lapply(1:Nsim, function(i) {
    #  set.seed(999*i)
    if (pi0 == 1) {
      foo <- makeSimCount2groups(counts = counts,
                                 Nsamples = Nsamples, Ngenes = Ngenes,
                                 sample_method = sample_method)
      return(foo)
    }

    if (pi0 < 1) {
      foo <- makeSimCount2groups(counts = counts,
                                 Nsamples = Nsamples, Ngenes = Ngenes,
                                 sample_method = sample_method)
      foo2 <- non_null_sim(counts = foo$counts,
                           condition = foo$condition,
                           pi0,
                           beta_args = beta_args)
      return(foo2)
    }
  })
}

#' @title simulation wrapper with flexible gene sampling method
#' 
#' @export
simulationWrapper.sampleGenes <- function(counts,
                              Nsim = 1, Nsamples = 50, Ngenes = NULL,
                              pi0 = NULL,
                              sample_method = c("all_genes", "per_gene"),
                              sampleGenes_method = c("high", "medium", "low"),
                              beta_args = args.big_normal(betapi = c(1),
                                                          betamu = c(0),
                                                          betasd = c(.8))) {
  simdata <- lapply(1:Nsim, function(i) {
    #  set.seed(999*i)
    if (pi0 == 1) {
      foo <- makeSimCount2groups.sampleGenes(counts = counts,
                                 Nsamples = Nsamples, Ngenes = Ngenes,
                                 sample_method = sample_method)
      rownames(foo$counts) <- paste0("gene.",c(1:nrow(foo$counts)))
      names(foo$is_nullgene) <- rownames(foo$counts)
      return(foo)
    }
    
    if (pi0 < 1) {
      foo <- makeSimCount2groups.sampleGenes(counts = counts,
                                 Nsamples = Nsamples, Ngenes = Ngenes,
                                 sample_method = sample_method,
                                 sampleGenes_method = sampleGenes_method)
      foo2 <- non_null_sim(counts = foo$counts,
                           condition = foo$condition,
                           pi0,
                           beta_args = beta_args)
      rownames(foo2$counts) <- paste0("gene.",c(1:nrow(foo2$counts)))
      names(foo2$is_nullgene) <- rownames(foo2$counts)
      return(foo2)
    }
  })

}


#' Sample genes
#' 
#' @description Methods for sampling genes
#' 
#' @export
sampleGenes <- function(counts, Ngenes, method = c("high", "medium", "low")) {
  quants <- quantile(rowMeans(counts))
  if (method == "high") {
    ind <- which(rowMeans(counts) > quants[4])
  }
  
  if (method == "medium") {
    ind <- which(rowMeans(counts) > quants[2] &
                   (rowMeans(counts) < quants[4]))
  }
  if (method == "low") {
    ind <- which(rowMeans(counts) < quants[2])
  }
  
  if (length(ind) > Ngenes) {
    ind <- sample(ind, Ngenes, replace = FALSE)
  }  
  return(ind)
}


#' @title Generate count matrix for genes of different expression levels
#'
#' @param counts Gene expression count matrix from a dataset.
#' @param Ngenes Number of genes in the simulated dataset.
#' @param Nsample Number of samples in each condition.
#' @param pi0 Proportion of null genes. Default to be 1.
#'
#' @export
makeSimCount2groups.sampleGenes <- function(counts, Ngenes = NULL, Nsamples = 50,
                                sampleGenes_method = c("high", "medium", "low"),
                                sample_method = c("per_gene", "all_genes")){
  output <- list(counts = counts)
  Nsamples_total <- 2*Nsamples
  
  if (sample_method == "per_gene") {
    # For each gene, randomly select 2*Nsamp samples from counts
    counts_subset <- t(apply(counts, 1, sampleingene, Nsamples_total=Nsamples_total))
    # Remove genes without any reads
    if ( sum(apply(counts_subset,1,sum) == 0) > 0) {
      counts_subset <- counts_subset[which(apply(counts_subset,1,sum)>0), ]
    } else {
      counts_subset <- counts_subset
    }
    # #Ngenes <- NROW(counts_subset)
    # output$counts <- counts_subset
  }
  
  if (sample_method == "all_genes") {
    samples_to_include <- sample((1:NCOL(counts)), Nsamples_total, replace = TRUE)
    counts_subset <- counts[, samples_to_include]
    # Remove genes without any reads
    counts_subset <- counts_subset[which(apply(counts_subset,1,sum)>0), ]
    # #Ngenes <- NROW(counts_subset)
    # output$counts <- counts_subset
  }
  
  if (!is.null(Ngenes)) {
    # genes_to_include <- sample(1:NROW(counts_subset), Ngenes, replace = FALSE)
    genesToInclude <- sampleGenes(counts_subset, Ngenes, method = sampleGenes_method)
    counts_subset <- counts_subset[genesToInclude, ]
    output$counts <- counts_subset
  } else {
    output$counts <- counts_subset
  }
  
  # assign samples to 2 arbitrary conditions
  condition <- rep(1:2, each = Nsamples)
  condition <- condition[sample(1:(Nsamples_total), size = Nsamples_total)]
  output$condition <- condition
  
  # reorder sample columns
  output$condition <- output$condition[order(output$condition)]
  output$counts <- output$counts[, order(output$condition)]
  
  # make colnames be unique
  colnames(output$counts) <- paste0("sample_",c(1:dim(output$counts)[2]))
  
  return(output)
}



#' @title Generate count matrix of all null genes
#'
#' @param counts Gene expression count matrix from a dataset.
#' @param Ngenes Number of genes in the simulated dataset.
#' @param Nsample Number of samples in each condition.
#' @param pi0 Proportion of null genes. Default to be 1.
#'
#' @examples
#' library(singleCellRNASeqHumanTungiPSC)
#' eset <- HumanTungiPSC
#' counts <- exprs(eset)[,pData(eset)$individual == "NA19101"]
#'
#' sim_counts <- makeSimCount2groups(counts,
#'                                   Ngenes = 100,
#'                                   Nsample = 20,
#'                                   sample_method = "all_genes")
#'
#' @export
makeSimCount2groups <- function(counts, Ngenes = NULL, Nsamples = 50,
                                sample_method = c("per_gene", "all_genes")){
  output <- list(counts = counts)
  Nsamples_total <- 2*Nsamples

  if (sample_method == "per_gene") {
    # For each gene, randomly select 2*Nsamp samples from counts
    counts_subset <- t(apply(counts, 1, sampleingene, Nsamples_total=Nsamples_total))
    # Remove genes without any reads
    if ( sum(apply(counts_subset,1,sum) == 0) > 0) {
      counts_subset <- counts_subset[which(apply(counts_subset,1,sum)>0), ]
    } else {
      counts_subset <- counts_subset
    }
    # #Ngenes <- NROW(counts_subset)
    # output$counts <- counts_subset
  }

  if (sample_method == "all_genes") {
    samples_to_include <- sample((1:NCOL(counts)), Nsamples_total, replace = TRUE)
    counts_subset <- counts[, samples_to_include]
    # Remove genes without any reads
    counts_subset <- counts_subset[which(apply(counts_subset,1,sum)>0), ]
    # #Ngenes <- NROW(counts_subset)
    # output$counts <- counts_subset
  }

  if (!is.null(Ngenes)) {
    genes_to_include <- sample(1:NROW(counts_subset), Ngenes, replace = FALSE)
    counts_subset <- counts_subset[genes_to_include, ]
    output$counts <- counts_subset
  } else {
    output$counts <- counts_subset
  }

  # assign samples to 2 arbitrary conditions
  condition <- rep(1:2, each = Nsamples)
  condition <- condition[sample(1:(Nsamples_total), size = Nsamples_total)]
  output$condition <- condition

  # reorder sample columns
  output$condition <- output$condition[order(output$condition)]
  output$counts <- output$counts[, order(output$condition)]

  # make colnames be unique
  colnames(output$counts) <- paste0("sample_",c(1:dim(output$counts)[2]))

  return(output)
}



#' @title Randomisation of sample at each gene
#'
#' @param count_onegene Count vector of one gene.
#' @param Nsamp_total Total number of samples in the simulated data.
#'
#' @export
sampleingene <- function(count_onegene, Nsamples_total){
  sample_subset <- sample(length(count_onegene), Nsamples_total)
  return(c(count_onegene[sample_subset]))
}

#' @title Simulate count matrix
#'
#' @param counts G*N null count matrix, rows are genes and columns are samples
#' @param args List of arguments, including
#'   \code{beta_args}: parameters for normal mixture distributions.
#'        \code{betapi} Probability vector for the k components in the normal mixture.
#'        \code{betamu} Mean vector for k mixture components.
#'        \code{betasd} Standard deviation vector for k mixture components.
#'   \code{pi0}: null proportion. If pi0=="random" then pi0 will be randomly selected from U(0,1)
#'
#' @examples
#' library(singleCellRNASeqMouseZeiselBrain)
#' eset <- get(data(MouseZeiselBrain))
#' counts <- exprs(eset)
#'
#' counts_null <- makeSimCount2groups(counts,
#'                                    Ngenes = 1000,
#'                                    Nsamples = 20,
#'                                    sample_method = "all_genes")
#' counts_sim <- non_null_sim(counts_null$counts,
#'                            counts_null$condition,
#'                            pi0 = .5,
#'                            beta_args = args.big_normal())
#'
#' @export
non_null_sim <- function(counts, condition, pi0,
                         beta_args = list(betapi, betamu, betasd) ){
  # Thinned effect sizes generated from normal mixture prior
  Ngenes <- dim(counts)[1]
  be <- make_normalmix(Ngenes,
                       beta_args = beta_args,
                       pi0)
  is_nullgene <- be$is_nullgene # null gene indicators
  beta <- be$beta
  # Use Poisson thinning to add effects to null data
  sim_list <- pois_thinning(counts, condition, beta)
  return(list(counts = sim_list$counts,
              condition = sim_list$condition,
              beta = beta,
              is_nullgene = is_nullgene))
}

#' @title Generate beta (effects) from normal mixture prior
#'
#' @param Ngenes Total number of genes (null + non-null in the data).
#' @param pi K-component probability vector. K = number of components in the mixture prior.
#' @param mu Length-K mean vector of the mixture components.
#' @param sd Length-K standard error vector of the mixture components.
#' @param pi0 Proportion of null genes.Default = "random", which selects a random number between 0 and 1 from a uniform distribution.
#'
#' @return null Binary vector of length Ngenes. 1 = null, 0 = no null.
#' @return scalar value pi0 proportion of null
#'
#' @examples
#' betas <- make_normalmix(100, pi = .5,
#'                         beta_args = args.big_normal())
#' @export
make_normalmix <- function(Ngenes, pi0,
                          beta_args = list(betapi, betamu, betasd)){
  if (pi0=="random"){
    pi0 = runif(1,0,1) #generate the proportion of true nulls randomly
  }
  k <- length(beta_args$betapi) # number of components
  comp <- sample(1:k, Ngenes, beta_args$betapi,replace=TRUE) #randomly draw a component
  is_nullgene <- (runif(Ngenes,0,1) < pi0)
  beta <- ifelse(is_nullgene, 0, rnorm(Ngenes, beta_args$betamu[comp], beta_args$betasd[comp]))
  return(list(beta = beta, pi0 = pi0, is_nullgene = is_nullgene))
}

#' @title Poisson thinning
#'
#' @param counts G genes by N samples null count matrix.
#' @param log2foldchanges Prior distribution of the betas (true effects).
#'
#' @return counts A simulated count matrix of non-null and null genes.
#'
#' @export
pois_thinning <- function(counts, condition, log2foldchanges){
  nsamples_per_group <- length(condition)/2
  is_nullgene <- log2foldchanges == 0
  log2foldchanges <- log2foldchanges[!is_nullgene]
  foldchanges <- 2^log2foldchanges

  # thin group A
  num_positive_foldchange <- sum(log2foldchanges > 0)
  which_positive_foldchange <- which(log2foldchanges > 0)

  counts[which(!is_nullgene)[which_positive_foldchange], which(condition == 1)] <-
    matrix(rbinom(num_positive_foldchange*nsamples_per_group,
            size = c(as.matrix(counts[which(!is_nullgene)[which_positive_foldchange],
                                      which(condition == 1)]) ),
            prob = rep(1/foldchanges[which_positive_foldchange], nsamples_per_group)),
            ncol = nsamples_per_group)

  # thin group B
  num_negative_foldchange <- sum(log2foldchanges < 0)
  which_negative_foldchange <- which(log2foldchanges < 0)

  counts[which(!is_nullgene)[which_negative_foldchange], which(condition == 2)] <-
    matrix(rbinom(num_negative_foldchange*nsamples_per_group,
            size = c(as.matrix(counts[which(!is_nullgene)[which_negative_foldchange],
                                      which(condition == 2)])),
            prob = rep(foldchanges[which_negative_foldchange], nsamples_per_group)),
            ncol = nsamples_per_group)

  return(list(counts = counts, condition = condition))
}


#' @title spiky prior
#'
#' @export
args.spiky <- function(...) {
  list(betapi=c(.4,.2,.2,.2),betamu=c(0,0,0,0),betasd=c(.25,.5,1,2)/sqrt(2*Nsamp-2))
}

#' @title near normal prior
#' @export
args.near_normal <- function(...) {
  list(betapi=c(2/3,1/3),betamu=c(0,0),betasd=c(1,2)/sqrt(2*Nsamp-2))
}

#' @title flat top prior
#'
#' @export
args.flat_top <- function(...) {
  list(betapi=rep(1/7,7),
       betamu=c(-1.5,-1,-0.5,0,0.5,1,1.5),
       betasd=rep(0.5,7)/sqrt(2*Nsamp-2))
}

#' @title Normal mixture prior. Default 1 component.
#'
#' @param Nsamp Number of samples per group
#' @param betapi Probabliy vector of k components in the normal mixture.
#' @param betamu Mean vector of k components in the normal mixture.
#' @param betasd Standard deviation of the effect sizes.
#' @param pi0 Fraction of null genes. Default = Random.
#'
#' @export
args.big_normal <- function(betapi = c(1), betamu = c(0), betasd = c(1)) {
  list(betapi = betapi,
       betamu = betamu,
       betasd = betasd )
}

#' @title bimodal prior
#'
#' @export
args.bimodal <- function(...) {
  list(betapi = c(0.5,0.5),
       betamu = c(-2,2),
       betasd=c(1,1)/sqrt(2*Nsamp-2))
}
