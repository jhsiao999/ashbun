#################################################################
#       ash-Poisson related code
#
# Objective:
# We consider several scenario of modality:
# - uni_natural: unimodal distribution with a non-zero mode (unimodal nonzero)
# - uni_zero: unimodal distribution with a zero mode (unimodal zero)
# - uni_natural_pm_zero: bimodal distribution with two components, including a zero component of point mass at zero, and a non-zero component of a mixture of distributions with non-zero mode
# - uni_natural_uni_zero: bimodal distribution with two components, including a mixture of distribution with mode at zero and a mixture of distribution with mode at non-zero count.
#################################################################

# Arguments for the scenarios
#' @title unimodal distribution with a non-zero mode
#' @export
args.uni_natural <- function(...) { list(mixcompdist = "halfuniform",
                                         mode = "estimate",
                                         prior = "uniform") }

#' @title unimodal distribution with mode at zero
#' @export
args.uni_zero <- function(...) { list(mixcompdist = "+uniform",
                                      mode = 0,
                                      prior = "uniform") }

#' @title bi-modal distribution, with a non-zero mode and point mass at zero
#' @export
args.uni_natural_pm_zero <- function(ash_out_uni_natural) {
  g_uni_natural <- get_fitted_g(ash_out_uni_natural)
  k <- length(g_uni_natural$pi)
  g_bimodal <- with(g_uni_natural,
                    unimix(pi=rep(1,k+1)/(k+1), a=c(0, a), b=c(0, b)) )
  list(g = g_bimodal,
       prior = "uniform") }

#' @title bi-modal distribution, with a non-zero mode and a mode at zero
#' @export
args.uni_natural_uni_zero <- function(ash_out_uni_zero) {
  g_uni_zero <- ash_out_uni_zero$fitted_g
  k <- length(g_uni_zero$pi)
  g_uni_zero$pi <- rep(1/k, k)

  list(gpart = g_uni_zero,
       mixcompdist = "halfuniform",
       prior = "uniform")
  }


#' @title  ash-Poisson wrapper for evaluating different bimodality scenarios
#'
#' @param counts gene by sample count matrix
#' @param args pre-defined arguments for the four modality scenarios: args_uninonzero for unimodal distribution modeling after a mixture of half uniform distributions with mode at non-zero count, args_zero fo runimodal distribution modeling after a mixture of positive uniform distributionw ith mode at zero, and finally args_bimodal1 for a two-component situation where there's a zero-component modeling point mass at 0 and a non-zero component modeling a mixture of half uniform distribution with mode at 0.
#'
#' @export
ashPoissonWrapper <- function(counts, args = c(args.uni_natural,
                                               args.uni_zero,
                                               args.uni_natural_pm_zero,
                                               args.uni_natural_uni_zero)) {
  library(ashr)
  ngenes <- NROW(counts)
  foo <- lapply(1:ngenes, function(i) {
    y <- as.numeric(counts[i,])

    fit.uni_natural <- ash(rep(0, length(y)), 1, lik = lik_pois(y),
                          mixcompdist = args.uni_natural()$mixcompdist,
                          mode = args.uni_natural()$mode,
                          prior = args.uni_natural()$prior,
                          optmethod = "mixIP",
                          control = list(rtol = 1e-15))

    fit.uni_zero <- ash(rep(0, length(y)), 1, lik = lik_pois(y),
                       mixcompdist = args.uni_zero()$mixcompdist,
                       mode = args.uni_zero()$mode,
                       prior = args.uni_zero()$prior,
                       optmethod = "mixIP",
                       control = list(rtol = 1e-15))

    input.uni_natural_pm_zero <- args.uni_natural_pm_zero(fit.uni_natural)

    input.uni_natural_uni_zero <- args.uni_natural_uni_zero(fit.uni_zero)

    fit.uni_natural_pm_zero <-ash(rep(0, length(y)), 1, lik = lik_pois(y),
                       g = input.uni_natural_pm_zero$g,
                       prior = input.uni_natural_pm_zero$prior,
                       optmethod = "mixIP",
                       control = list(rtol = 1e-15))

    fit.uni_natural_uni_zero <-ash(rep(0, length(y)), 1, lik = lik_pois(y),
                                  gpart = input.uni_natural_uni_zero$gpart,
                                  mixcompdist = input.uni_natural_uni_zero$mixcompdist,
                                  prior = input.uni_natural_uni_zero$prior,
                                  optmethod = "mixIP",
                                  control = list(rtol = 1e-15))

    loglik <- data.frame(uni_natural = get_loglik(fit.uni_natural),
                         uni_zero = get_loglik(fit.uni_zero),
                         uni_natural_pm_zero = get_loglik(fit.uni_natural_pm_zero),
                         uni_natural_uni_zero = get_loglik(fit.uni_natural_uni_zero))

    return(list(g_fitted.uni_natural = get_fitted_g(fit.uni_natural),
                g_fitted.uni_zero = get_fitted_g(fit.uni_zero),
                g_fitted.uni_natural_pm_zero = get_fitted_g(fit.uni_natural_pm_zero),
                g_fitted.uni_natural_uni_zero = get_fitted_g(fit.uni_natural_uni_zero),
                loglik = loglik) )
  })
  return(list(g_fitted.uni_natural = lapply(foo, "[[", 1) ,
              g_fitted.uni_zero = lapply(foo, "[[", 2) ,
              g_fitted.uni_natural_pm_zero = lapply(foo, "[[", 3) ,
              g_fitted.uni_natural_uni_zero = lapply(foo, "[[", 4) ,
              loglik = do.call(rbind, lapply(foo, "[[", 5) ) ) )
}
