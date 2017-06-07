#' @title Hurdle negative binomial
#' @export
negbin_hurdle <- function(df, norm_factor) {
  fit <- try(pscl::hurdle(counts ~ group - 1 + offset(log(norm_factor)) | offset(log(norm_factor)),
                    data = df,
                    dist = "negbin"),
             silent = TRUE)
  if ("hurdle" %in% attr(fit, "class")) {
    # resid <- residuals(fit)
    # mse <- sum(resid^2)/NCOL(counts)
    exp_zero <- sum(predict(fit, type = "prob")[,1])
    exp_coef <- coef(fit, "count")[1]
    exp_theta <- fit$theta
#    exp_coef_se <- vcov(fit)[2,2]
    logLL <- logLik(fit)
    data.frame(exp_zero = exp_zero,
               exp_coef = exp_coef,
               exp_disp = exp_theta,
#               exp_coef_se = exp_coef_se,
               logLL = logLL)
  } else {
    return(rep(NA, 4))
  }
}


#' @title Zero-inflated negative binomial
#'
#' @export
negbin_zif <- function(counts, group, norm_factor, control = list(EM = FALSE)) {

  # check control list
  control_default <- list(EM = FALSE)
  control_names <- names(control)
  if (!all(control_names %in% names(control_default)))
    stop("unknown names in control: ",
         control_names[!(control_names %in% names(control_default))])
  control <- modifyList(control_default, control)
  EM <- control$EM

  # fit model
  fit <- try(pscl::zeroinfl(counts ~ group + offset(log(norm_factor)),
                       dist = "negbin", EM = EM),
             silent = TRUE)
  if ("zeroinfl" %in% attr(fit, "class")) {
    # resid <- residuals(fit)
    # mse <- sum(resid^2)/NCOL(counts)
    est_zero <- sum(predict(fit, type = "prob")[,1])
    # est_coef_zero <- coef(fit, "zero")[2]
    # est_coef_exprs <- coef(fit, "count")[2]
    # est_theta <- fit$theta

#    exp_coef_se <- vcov(fit, "count")[2,2]
#    logLL <- logLik(fit)
    # data.frame(est_count_zero = est_count_zero,
    #            est_coef_zero = est_coef_zero,
    #            est_coef_exprs = est_coef_exprs,
    #            est_disp = est_theta)
#               exp_coef_se = exp_coef_se,
#               logLL = logLL)
    data.frame(est_zero = est_zero)

  } else {

    return(rep(NA, 1))

  }

}





#' @title Fit Negative binomial to gene expression counts (univariate)
#'
#' @export
#'
negbin <- function(counts, group, norm_factor = NULL) {

  if (!is.null(norm_factor)) {
    offset_constant <- log(norm_factor)
    } else {
    offset_constant <- 0
    }

  fit <- try(MASS::glm.nb(counts ~ group + offset(offset_constant)),
             silent = TRUE)

  if ("glm" %in% attr(fit, "class")) {
    est_zero <- sum(dnbinom(0, mu = fitted(fit), size = fit$theta))
    # est_coef_exprs <- coef(fit)[2]
    # est_theta <- fit$theta

    #    exp_coef_se <- vcov(fit)[2,2]
#    logLL <- logLik(fit)
    data.frame(est_zero = est_zero)
    #               exp_coef_se = exp_coef_se,
#               logLL = logLL)
  } else {

    return(rep(NA, 1))

  }

}
