# This document contains codes for method performance evaluation

#' @title Compute sensitivity or true positive rate given fixed FDR
#'
#' @param true_positive_rate numeric vector of sensitivity.
#' @param false_positive_rate numeric vector of 1 - specificity.
#' @param fdr_cutoff default .05
#'
#' @return
#'    \code{tpr} true postive rate at the given FDR threshold.
#' @export
getTPR <- function(true_positive_rate, false_positive_rate,
                   fdr_cutoff = .05) {

  df <- data.frame(true_positive_rate = true_positive_rate,
                   false_positive_rate = false_positive_rate)
  df <- df[order(df$false_positive_rate, df$true_positive_rate), ]
  which_max_fdr <- max(which(df$false_positive_rate < .05))
  tpr <- with(df, true_positive_rate[which_max_fdr])
  return(tpr)
}



#' @title Compute sensitivity or true positive rate given fixed FDR using pROC output
#'
#' @param pROC_output output from pROC package
#' @param fdr_cutoff default .05
#'
#' @return
#'    \code{tpr} true postive rate at the given FDR threshold.
#'
#' @export
getTPR.pROC <- function(response, predictor,
                   fdr_cutoff = .05) {

  if (length(unique(response)) == 1) {
    tpr <- NA
  } else {
  pROC_output <- getROC(response, predictor)
  tpr <- with(pROC_output, getTPR(sensitivities,
                                  1- specificities, fdr_cutoff = .05))
  }
  return(tpr)
}




#' Compute AUC
#'
#' @param response Binary indicator of true/false
#' @param predictor Numeric and continous variable.
#'
#' @export
getAUC <- function(response, predictor) {
  library(pROC)
  if (length(unique(response)) == 1) {
    auc <- NA
  } else {
    auc <- pROC::roc(response, predictor)$auc
  }
  return(auc)
}

#' @title compute ROC related informatino using pROC package
#'
#' @export
getROC <- function(response, predictor) {
  library(pROC)

  if (length(unique(response)) == 1) {
    auc <- NA
  } else {
    roc <- pROC::roc(response, predictor)
  }
  return(roc)
}


#' @title Plot ROC curve average
#'
#' @param roc_output data.frame of roc output, including columns of true positive rate, 
#'                  false positive rate, and simulation dataset index.
#' @param length.out default 20. FPR is sorted and grouped into \code{length.out} internals.
#'
#' @return 
#'    \code{output} data.frame of average true positive rates given false positive rates.
#'    
#' @export
getROC.average <- function(roc_output, length.out = 20) {
    range_fpr <- range(roc_output$FPR)
    seq_fpr <- seq.int(from=range_fpr[1], to = range_fpr[2], length.out = length.out+1)
    avg_tpr <- rep(NA, length.out+1)
    for (index in 1:(length.out)) {
      avg_tpr[index] <- with(roc_output, 
                             mean( TPR[which(FPR >= seq_fpr[index] & FPR < seq_fpr[index+1])] ) ) 
    }
    output <- data.frame(TPR = avg_tpr, FPR = seq_fpr)
    output <- output[!is.na(output$TPR),]
    
    return(output)
}
  
# @param response Binary indicator of true/false. True = Non-null gene and FALSE = Null gene.
# @param predictor Significance value (e.g., lfsr or p-value).
# @param fdr_cutoff False Discovery Rate threshold.
# getTPR <- function(response, predictor, fdr_cutoff = .05) {
#   if (length(unique(response)) == 1) {
#     tpr <- NA
#   } else {
#     ordered_list <- data.frame(response, predictor)
#     ordered_list <- ordered_list[order(ordered_list$predictor), ]
#
#     # cumulative number of null genes
#     cc <- cumsum(ordered_list$response == 0)
#
#     # fraction of false discoveries at each cutoff
#     ff <- cc/sum(ordered_list$response == 0)
#
#     # find the cutoff at which fdr is no more than threshold
#     cutoff_index <- sum(ff < fdr_cutoff)
#
#     # find the significance value corresponding to the fdr cutoff
#     cutoff_sigvalue <- ordered_list$predictor[cutoff_index]
#
#     # compute true positive rate
#     tpr <- sum(response == 1 & predictor < cutoff_sigvalue)/sum(response == 1)
#   }
#   return(tpr)
# }
