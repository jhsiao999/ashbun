source("simulateExpression.R")
simulationWrapper <- function(counts,
                              Nsample, Ngenes = NULL,
                              pi0 = NULL,
                              sample_method = c("all_genes", "per_gene"),
                              beta_args = args.big_normal(betapi = c(1),
                                                          betamu = c(0),
                                                          betasd = c(1))) {
    if (pi0 == 1) {
      foo <- makeSimCount2groups(counts = counts,
                                 Nsamp = Nsample, Ngenes = Ngenes,
                                 sample_method = sample_method)
      return(foo)
    }

    if (pi0 < 1) {
      foo <- makeSimCount2groups(counts = counts,
                                 Nsamp = Nsample, Ngenes = Ngenes,
                                 sample_method = sample_method)
      foo2 <- non_null_sim(counts = foo$counts,
                           condition = foo$condition,
                           pi0,
                           beta_args = beta_args)
      return(foo2)
    }
}
beta_args = do.call(beta_function, list(betapi = betapi, betamu = betamu, betasd = betasd))
output = simulationWrapper(counts, Nsample, Ngenes, pi0, sample_method, beta_args)
