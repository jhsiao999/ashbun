# This document contains codes for method performance evaluation

#' Compute AUC
#'
#' @param response Binary indicator of true/false
#' @param predictor Numeric and continous variable.
#'
getAUC <- function(response, predictor) {
  if (length(unique(response)) == 1) {
    auc <- NA
  } else {
    auc <- pROC::roc(response, predictor)$auc
  }
  return(auc)
}


getROC <- function(response, predictor) {
  if (length(unique(response)) == 1) {
    auc <- NA
  } else {
    roc <- pROC::roc(response, predictor)
  }
  return(roc)
}
