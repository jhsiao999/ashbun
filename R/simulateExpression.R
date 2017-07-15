###########################################################
#  Code for simulating single cell expression count data
#
#  Date:2017-04-03
#  Code adapted from Mengyin Lu's work
###########################################################


#' @title Generate count matrix of all null genes
#'
#' @param counts Gene expression count matrix from a dataset (usually consisting of one biological condition).
#' @param Ngene Number of genes in the simulated dataset.
#' @param Nsamp Number of samples per condition in the simulated data.
#' @param nullpi Proportion of null genes.
#'
#' @export
makeSimCount2groups <- function(counts, Ngene = NULL, Nsamp,
                              sample_method = c("per_gene", "all_genes"),
                              control = list(pi0 = NULL)){
   output <- list(counts = counts)

  if (!is.null(Ngene)) {
    genes_to_include <- sample(1:NROW(counts), Ngene, replace = FALSE)
    counts_subset <- counts[genes_to_include, ]
  } else {
    counts_subset <- counts
  }

  if (sample_method == "per_gene") {
    # For each gene, randomly select 2*Nsamp samples from counts
    counts_subset <- t(apply(counts_subset, 1, sampleingene, Nsamp=2*Nsamp))
    # Remove genes without any reads
    if ( sum(apply(counts_subset,1,sum) == 0) > 0) {
      counts_subset <- counts_subset[which(apply(counts_subset,1,sum)>0), ]
    } else {
      counts_subset <- counts_subset
    }
    #Ngene <- NROW(counts_subset)
    output$counts <- counts_subset
  }

  if (sample_method == "all_genes") {
    samples_to_include <- sample((1:NCOL(counts)), 2*Nsamp, replace = TRUE)
    counts_subset <- counts_subset[, samples_to_include]
    # Remove genes without any reads
    counts_subset <- counts_subset[which(apply(counts_subset,1,sum)>0), ]
    #Ngene <- NROW(counts_subset)
    output$counts <- counts_subset
  }

  if (!is.null(control$pi0)) {
    if (control$pi0 == 1) {
    # assign samples to 2 arbitrary conditions
    condition <- rep(1:2, each = Nsamp)
    condition <- condition[sample(1:(Nsamp*2), size = Nsamp*2)]
    output$condition <- condition
    output$pi0 <- 1
    }
  }

  return(output)
}

#' @title Randomisation of sample at each gene
#'
#' @param count_onegene Count vector of one gene.
#' @param Nsamp Total number of samples in the simulated data (both conditions included).
#'
#' @export
sampleingene <- function(count_onegene, Nsamp){
  sample_subset <- sample(length(count_onegene), Nsamp)
  return(c(count_onegene[sample_subset]))
}

#' @title Simulate count matrix
#'
#' @param counts G*N null count matrix, rows are genes and columns are samples
#' @param args List of arguments, including
#'  - betaargs: parameters for the normal mixture distribution to generate signals
#'  - pi0: null proportion. If pi0=="random" then pi0 will be randomly selected from U(0,1)
#'
#' @export
non_null_sim <- function(counts, args){
  # Thinned effect sizes generated from normal mixture prior
  ngene <- dim(counts)[1]
  be <- make_normalmix(ngene,
                     args$betaargs$betapi, args$betaargs$betamu, args$betaargs$betasd,
                     args$pi0)
  null <- be$null # null gene indicators
  beta <- be$beta
  # Use Poisson thinning to add effects to null data
  sim_list = pois_thinning(counts, beta)
  return(list(counts=sim_list$counts,
              condition = sim_list$condition,
              beta = beta,
              null=null))
}

#' @title Generate beta (effects) from normal mixture prior
#'
#' @param ngene Total number of genes (null + non-null in the data).
#' @param pi K-component probability vector. K = number of components in the mixture prior.
#' @param mu Length-K mean vector of the mixture components.
#' @param sd Length-K standard error vector of the mixture components.
#' @param pi0 Proportion of null genes.Default = "random", which selects a random number between 0 and 1 from a uniform distribution.
#'
#' @return null Binary vector of length ngene. 1 = null, 0 = no null.
#' @return scalar value pi0 proportion of null
#'
#' @export
make_normalmix = function(ngene, pi, mu, sd, pi0){
  if (pi0=="random"){
    pi0 = runif(1,0,1) #generate the proportion of true nulls randomly
  }
  k = length(pi) # number of components
  comp = sample(1:k,ngene,pi,replace=TRUE) #randomly draw a component
  isnull = (runif(ngene,0,1) < pi0)
  beta = ifelse(isnull, 0, rnorm(ngene,mu[comp],sd[comp]))
  return(list(beta=beta, pi0=pi0, null=isnull))
}

#' @title Poisson thinning
#'
#' @param counts G genes by N samples null count matrix.
#' @param log2foldchanges Prior distribution of the betas (true effects).
#'
#' @return counts A simulated count matrix of non-null and null genes.
#'
#' @export
pois_thinning = function(counts, log2foldchanges){
  nsamp = dim(counts)[2]/2
  null = (log2foldchanges==0)
  log2foldchanges = log2foldchanges[!null]
  foldchanges = 2^log2foldchanges

  # thin group A
  counts[which(!null)[log2foldchanges>0],1:nsamp] =
    matrix(rbinom(sum(log2foldchanges>0)*nsamp,
                  size=c(as.matrix(counts[which(!null)[log2foldchanges>0],1:nsamp])),
                  prob=rep(1/foldchanges[log2foldchanges>0],nsamp)),ncol=nsamp)
  # thin group B
  counts[which(!null)[log2foldchanges<0],(nsamp+1):(2*nsamp)] =
    matrix(rbinom(sum(log2foldchanges<0)*nsamp,
                  size=c(as.matrix(counts[which(!null)[log2foldchanges<0],(nsamp+1):(2*nsamp)])),
                  prob=rep(foldchanges[log2foldchanges<0],nsamp)),ncol=nsamp)

  condition <- c(rep(1, nsamp), rep(2, nsamp))

  return(list(counts=counts, condition = condition))
}


#' @title spiky prior
#'
#' @export
args.spiky <- function(nsamp, pi0 = "random") {
  list(pi0=pi0,
       betaargs=list(betapi=c(.4,.2,.2,.2),betamu=c(0,0,0,0),betasd=c(.25,.5,1,2)/sqrt(2*nsamp-2)))
}

#' @title near normal prior
#' @export
args.near_normal <- function(nsamp, pi0 = "random") {
  list(pi0=pi0,
       betaargs=list(betapi=c(2/3,1/3),betamu=c(0,0),betasd=c(1,2)/sqrt(2*nsamp-2)))
}

#' @title flat top prior
#'
#' @export
args.flat_top <- function(nsamp, pi0 = "random") {
  list(pi0=pi0,
       betaargs=list(betapi=rep(1/7,7),betamu=c(-1.5,-1,-0.5,0,0.5,1,1.5),betasd=rep(0.5,7)/sqrt(2*nsamp-2)))
}

#' @title big normal prior
#'
#' @param nsamp Number of samples per group
#' @param betasd Standard deviation of the effect sizes. Input is used to compute standard error: betasd/sqrt(2*nsamp-2).
#' @param pi0 Fraction of null genes. Default = Random. 
#'
#' @export
args.big_normal <- function(nsamp, betasd = NULL, pi0 = "random") {
  if (betasd == NULL) { betasd.0 <- 4/sqrt(2*nsamp-2) }
  if (betasd != NULL) { betasd.0 <- betasd}

  list(pi0=pi0,
       betaargs=list(betapi=c(1),betamu=c(0),betasd=betasd.0 ))
}

#' @title bimodal prior
#'
#' @export
args.bimodal <- function(nsamp, pi0 = "random") {
  list(pi0=pi0,
       betaargs=list(betapi=c(0.5,0.5),betamu=c(-2,2),betasd=c(1,1)/sqrt(2*nsamp-2)))
}
