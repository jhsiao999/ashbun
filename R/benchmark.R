# This document contains codes for method performance evaluation

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

#' @title compute ROC
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
