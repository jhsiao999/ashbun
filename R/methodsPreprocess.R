# typical steps of preprocessing, includign prior count and filtering


#----- filtering criteria

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
#'    \code{fraction_expressed} fraction of samples qualified as "expressed".
#'    \code{index_filter} TRUE/FALSE passing the filter.
#'
#' @export
filterFeatures.fractionExpressed <- function(counts, 
                                             thresholdDetection = 1,
                                             fractionExpressed = .25,
                                             control = list(plot=FALSE)) {

 fraction_expressed <- rowMeans(counts > thresholdDetection)
 index_filter <- which(fraction_expressed > fractionExpressed)
 
 if (control$plot) {
  
   ggplot2::qplot(fraction_expressed, geom = "histogram",
                  ylab = "Frequency", 
                  xlab = paste("Fraction of samples detected as expressed (>",
                                thresholdDetection, ")")) +
     geom_vline(xintercept = fractionExpressed, col = "red")
 }

 return(list(fraction_expressed = fraction_expressed,
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
