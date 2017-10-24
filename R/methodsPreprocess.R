# typical steps of preprocessing, includign prior count and filtering


#----- filtering criteria

#' @title Wrapper for all filtering steps
#'
#' @param counts Gene by sample expression count matrix (G by N).
#' @param thresholdDetection minimum count per gene/sample. Default value = 1.
#' @param fractionExpressed fraction of samples expressed (above thresholdDetection).
#'        Default value = .01.
#'
#' @return
#'    \code{counts} filtered count matrix.
#'
#' @export
filter.Wrapper <- function(counts, condition,
                           thresholdDetection = 1, fractionExpressed = .01,
                           is_nullgene = NULL) {

  counts_filtered <- filter.excludeAllZeros(counts)
  featuresToInclude <- filterFeatures.fractionExpressed(counts_filtered,
                            condition=condition,
                            thresholdDetection=thresholdDetection,
                            fractionExpressed=fractionExpressed)$index_filter

  samplesToInclude <-  filterSamples.fractionExpressed(counts_filtered,
                            thresholdDetection=thresholdDetection,
                            fractionExpressed=fractionExpressed)$index_filter

  counts_filtered <- counts_filtered[featuresToInclude, samplesToInclude]
  condition_filtered <- condition[samplesToInclude]

  if (!is.null(is_nullgene)) {
    is_nullgene_filtered <- is_nullgene[featuresToInclude]
    } else {
    is_nullgene_filtered <- NULL
    }

  return(list(counts = counts_filtered,
              condition = condition_filtered,
              is_nullgene = is_nullgene_filtered))
}


#' @title Filter all-zero samples and features.
#'
#' @param counts Gene by sample expression count matrix (G by N).
#'
#' @return
#'    \code{counts} filtered count matrix.
#'
#' @export
filter.excludeAllZeros <- function(counts) {

  samplesToInclude <- which(colSums(counts > 0) > 0)
  genesToInclude <- which(rowSums(counts > 0) > 0)

  counts <- counts[genesToInclude, samplesToInclude]

  return(counts)
}



#' @title Filter genes by fraction of samples detected as expressed
#'
#'
#' @param counts Gene by sample expression count matrix (G by N).
#' @param thresholdDetection minimum count (or tpm or cpm, depends on the data) required to
#'                           qualified as "detected as expressed".
#' @param fractionExpressed fraction of samples qualified as detected as expressed.
#' @param control List with the following arguments
#'   \code{plot} TRUE/FALSE to output diagnositic plot.
#'
#' @return List of the following objects
#'    \code{fraction_expressed_cond} fraction of samples qualified as "expressed" in each condition.
#'    \code{index_filter} TRUE/FALSE passing the filter.
#'
#' @export
filterFeatures.fractionExpressed <- function(counts, condition,
                                             thresholdDetection = 1,
                                             fractionExpressed = .01,
                                             control = list(plot=FALSE)) {
 condition_vals <- unique(condition)
 fraction_expressed_cond <- lapply(unique(condition), function(index) {
          tmp <- rowMeans(counts[, condition==condition_vals[index]] > thresholdDetection)
          return(tmp) })
 index_filter <- which(fraction_expressed_cond[[1]] > fractionExpressed & fraction_expressed_cond[[2]] > fractionExpressed)

 if (control$plot) {

   ggplot2::qplot(fraction_expressed, geom = "histogram",
                  ylab = "Frequency",
                  xlab = paste("Fraction of samples detected as expressed (>",
                                thresholdDetection, ")")) +
     geom_vline(xintercept = fractionExpressed, col = "red")
 }

 return(list(fraction_expressed_cond = fraction_expressed_cond,
             index_filter = index_filter))
}



#' @title Filter samples by fraction of genes detected as expressed
#'
#'
#' @param counts Gene by sample expression count matrix (G by N).
#' @param thresholdDetection minimum count (or tpm or cpm, depends on the data) required to
#'                           qualified as "detected as expressed".
#' @param fractionExpressed fraction of genes qualified as detected as expressed.
#' @param control List with the following arguments
#'   \code{plot} TRUE/FALSE to output diagnositic plot.
#'
#' @return List of the following objects
#'    \code{fraction_expressed} fraction of genes qualified as "expressed".
#'    \code{index_filter} TRUE/FALSE passing the filter.
#'
#' @export
filterSamples.fractionExpressed <- function(counts,
                                             thresholdDetection = 1,
                                             fractionExpressed = .25,
                                             control = list(plot=FALSE)) {

  fraction_expressed <- colMeans(counts > thresholdDetection)
  index_filter <- which(fraction_expressed > fractionExpressed)

  if (control$plot) {

    ggplot2::qplot(fraction_expressed, geom = "histogram",
                   ylab = "Frequency",
                   xlab = paste("Fraction of genes detected as expressed (>",
                                thresholdDetection, ")")) +
      geom_vline(xintercept = fractionExpressed, col = "red")
  }

  return(list(fraction_expressed = fraction_expressed,
              index_filter = index_filter))
}




#----- pseudocount
