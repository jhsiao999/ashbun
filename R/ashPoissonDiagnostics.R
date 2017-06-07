#################################################
#     ash Poisson diagnostics
#################################################



# Some functions for ash Poisson diagnostics


# Utility functions from Mengyin

dens_unimix_sing <- function(x,pi,a,b){
  sum((x>=a & x<b)/(b-a)*pi, na.rm=TRUE)
}

dens_unimix <- function(g, x){
  sapply(x, dens_unimix_sing, pi=g$pi, a=g$a, b=g$b)
}
