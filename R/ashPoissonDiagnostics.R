#################################################
#     ash Poisson diagnostics
#################################################


# Utility functions from Mengyin

#' @title compute density for a single component in the mixture prior
#'
#' @export
dens_unimix_sing <- function(x,pi,a,b){
  sum((x>=a & x<b)/(b-a)*pi, na.rm=TRUE)
}

#' @title compute density for all components in the mixture prior
#'
#' @export
dens_unimix <- function(g, x){
  sapply(x, dens_unimix_sing, pi=g$pi, a=g$a, b=g$b)
}
