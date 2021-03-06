% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methodsNormalize.R
\name{normalize.scnorm}
\alias{normalize.scnorm}
\title{SCnorm 1.1.0}
\usage{
normalize.scnorm(counts, condition, control = list(save_modelFit = FALSE,
  FilterCellNum = 10, K = NULL, NCores = 4, useSpikes = FALSE, reportSF = TRUE))
}
\arguments{
\item{counts}{Gene by sample expression count matrix (G by N).}

\item{condition}{Binary vector of length N indicating sample biological condition.}

\item{control}{A list with control arguments, including
\code{save_modelFit} TRUE to output the complete SCnorm output.
\code{FilterCellNum} SCnorm argument. Minimum number of cells required to
                     estimate the count-depth relationship.
\code{K} SCnorm argument. Number of quantile groups.
         Default NULL: letting SCnorm find the optimal K.
\code{Ncores} SCnorm argument. Number of cores to use. Default 4.
\code{useSpikes} SCnorm argument. Whether to use spike-ins to perform the between-condition
                 scaling using spike-ins. Assume that the spike-in labels start with ERCC
\code{reportSF} SCnorm argument. Whehter to provide a matrix of scaling factors in the output.}
}
\value{
List of the following objects
   \code{count_normed} numeric vector of the scale factors for library size.
   \code{model_output} SCnorm complete output.
}
\description{
The method forces to output diagnoistic plots. The PLOT = T/F option is obsolete.
}
