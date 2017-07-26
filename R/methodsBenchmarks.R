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

  if (length(unique(response)) == 1) {
    tpr <- NA
  } else {
    df <- data.frame(true_positive_rate = true_positive_rate,
                     false_positive_rate = false_positive_rate)
    df <- df[order(df$false_positive_rate, df$true_positive_rate), ]
    which_max_fdr <- max(which(df$false_positive_rate < .05))
    tpr <- with(df, true_positive_rate[which_max_fdr])
  }
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
    pROC_output <- ashbun::getROC(response, predictor)
    tpr <- with(pROC_output, ashbun::getTPR(sensitivities,
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
  if (length(unique(response)) == 1) {
    auc <- NA
  } else {
    roc <- pROC::roc(response, predictor)
  }
  return(roc)
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
