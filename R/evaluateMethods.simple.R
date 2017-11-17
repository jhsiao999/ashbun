# evaluateMethods.simple.R
#
# Description:
#   This file is for running defautl setting of each method.
#
# Content (from higher to lower-level functions)
#
# query.evaluation: inputs a single count table and outputs evaluation results of
#                   multiple normalization methods and DE methods
#
# query.methodsMeanExpression: inputs a single count table and outputs p-values of
#                   multiple DE methods
#

#' @title Evaluate multiple normalization methods and multiple DE methods
#'
#' @author Chiaowen Joyce Hsiao
#'
#' @export
query.evaluation.simple <- function(counts, condition, is_nullgene,
                                    methodsMeanExpression,
                                    report.control = list(fdr_control_threshold = .05), 
                                    nsim = NULL) {
  
  library(ashr)
  
  if (is.null(rownames(counts))) {
    names(is_nullgene) <- paste0("gene.", c(1:nrow(counts)))
  }
  if (!is.null(rownames(counts))) {
    names(is_nullgene) <- rownames(counts)
  }
  
  message("Simulation data ", nsim, "\n")

  results <- query.methodsMeanExpression.simple(counts, condition, is_nullgene,
                  default = TRUE,
                  methodsMeanExpression = methodsMeanExpression)
  message("Analysis done!", "\n")
  

  message("Evaluating results", "\n")
  
  output <- vector("list", 3)
  names(output) <- c("fdr_control", "roc", "pi0")
  
  list_methodsMeanExpression <- as.character(unique(results$res_longformat$methodsMeanExpression))
  
  # how well the method calibrates (FDR control) at pvalue < .05
  pval_summarize <- results$res_longformat
  pval_summarize$is_nullgene <- sapply(pval_summarize$genes, function(x) 
    is_nullgene[which(names(is_nullgene) %in% x)] )
  
  output$fdr_control <- do.call(rbind, lapply(seq_along(list_methodsMeanExpression),
    function(index_meanExpression) {
        df_sub <- subset(pval_summarize,
                         methodsMeanExpression == list_methodsMeanExpression[index_meanExpression])
        df_sub <- df_sub[!is.na(df_sub$pvalues),]
        df_sig <- df_sub[df_sub$pvalues < report.control$fdr_control_threshold, ]
        
        foo <- data.frame(fdr_control = sum(!df_sig$is_nullgene)/nrow(df_sig),
                          methodsMeanExpression = list_methodsMeanExpression[index_meanExpression])
        return(foo)
      }) )

  # compute receiver operating curve
  output$roc <- do.call(rbind, lapply(seq_along(list_methodsMeanExpression),
    function(index_meanExpression) {
        df_sub <- subset(pval_summarize,
                  methodsMeanExpression == list_methodsMeanExpression[index_meanExpression])
        pvals_sub <- df_sub$pvalues
        is_nullgene_sub <- df_sub$is_nullgene
        
        roc_output <- getROC(response = is_nullgene_sub, predictor = pvals_sub)
        foo <- data.frame(TPR = roc_output$sensitivities,
                FPR = 1- roc_output$specificities,
                methodsMeanExpression = list_methodsMeanExpression[index_meanExpression])
        return(foo)
      }) )
       
  
  # null proportion estimates
  output$pi0 <- do.call(rbind, lapply(seq_along(list_methodsMeanExpression),
    function(index_meanExpression) {
      df_sub <- subset(pval_summarize,
                       methodsMeanExpression == list_methodsMeanExpression[index_meanExpression])
      df_sub <- df_sub[!is.na(df_sub$pvalues),]
      
      if (mean(!is.na(df_sub$betahat) & !is.na(df_sub$sebetahat) & !is.na(df_sub$df)) > .9) {
        ash_res <- ash(betahat=df_sub$betahat, 
                       sebetahat=df_sub$sebetahat)
        foo <- data.frame(pi0 = get_pi0(ash_res),
                          methodsMeanExpression = list_methodsMeanExpression[index_meanExpression])
      } else {
        foo <- data.frame(pi0 = NA,
                          methodsMeanExpression = list_methodsMeanExpression[index_meanExpression])
      }
      return(foo)
    }) )
  
  
  message("DONE!", "\n")
  
  output$fdr_control$nsim <- nsim
  output$roc$nsim <- nsim
  output$pi0$nsim <- nsim
  
  return(output)
}


#' @title Runs multiple DE methods under their default setting
#'
#' @param counts Gene by sample expression count matrix (G by N). Use filtered count data.
#' @param counts_normed Normalized expression count matrix
#'                      (typicall CPM with normlized library size).
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param libsize_factors Numeric vector of scale factors for library size factors.
#' @param methodsMeanExpression Chararacter vector of evaluted methods. To run all methods, use
#'                             c("DESeq2", "limmaVoom", "edgeR","BPSC", "MAST", "ROTS")
#'
#' @author Chiaowen Joyce Hsiao
#'
#' @export
query.methodsMeanExpression.simple <- function(counts, condition, 
                                  is_nullgene, default = TRUE,
                                  methodsMeanExpression) {
  
  names(is_nullgene) <- rownames(counts)
  
  methods <-  c("DESeq2", "edgeR", "limmaVoom",
                "BPSC", "MAST")
  methods_function <- c("methodWrapper.DESeq2",
                        "methodWrapper.edgeR",
                        "methodWrapper.limmaVoom",
                        "methodWrapper.bpsc",
                        "methodWrapper.mast")
  which_methods <- which(methods %in% methodsMeanExpression)

  # run methods for raw counts
  output <- vector("list", length = length(which_methods))
  names(output) <- methods[which_methods]
  
  print(names(output))

  if ("DESeq2" %in% names(output)) {
    cat("DESeq2", "\n")
    output[["DESeq2"]] <- methodWrapper.DESeq2(counts, condition, libsize_factors = NULL,
                             default = TRUE,
                             control = list(save_modelFit = FALSE,
                                            independentFiltering = TRUE,
                                            cooksCutoff = TRUE))
  }

  if ("edgeR" %in% names(output)) {
    cat("edgeR", "\n")
    output[["edgeR"]] <- methodWrapper.edgeR(counts, condition, libsize_factors = NULL,
                             default = TRUE,
                             control = list(save_modelFit = FALSE))
  }

  if ("limmaVoom" %in% names(output)) {
    cat("limmaVoom", "\n")
    output[["limmaVoom"]] <- methodWrapper.limmaVoom(counts, condition,
                                default = TRUE,
                                control = list(save_modelFit = FALSE))
  }

  if ("BPSC" %in% names(output)) {
    cat("BPSC", "\n")
    output[["BPSC"]] <- methodWrapper.bpsc(counts, condition,
                                           default=TRUE,
                                           control = list(save_modelFit = FALSE,
                                                          estIntPar = TRUE,
                                                          useParallel = TRUE))
  }

  if ("MAST" %in% names(output)) {
    cat("MAST", "\n")
    output[["MAST"]] <- methodWrapper.mast(counts, condition,
                                           default=TRUE,
                                           control = list(save_modelFit = FALSE,
                                                          include_cdr = TRUE))
  }

  res_longformat <- vector("list", length(methodsMeanExpression))
  names(res_longformat) <- methodsMeanExpression

  for (index in 1:length(methodsMeanExpression)) {
    res_longformat[[index]] <- data.frame(genes=paste0("gene.",c(1:length(output[[index]]$pvalue))),
                                pvalues=output[[index]]$pvalue,
                                betahat=output[[index]]$betahat,
                                sebetahat=output[[index]]$sebetahat,
                                df=output[[index]]$df,
                                methodsMeanExpression=rep(methodsMeanExpression[index],
                                                          length(output[[index]]$pvalue)))
  }
  res_longformat <- do.call(rbind, res_longformat)

  return(list(data = list(counts=counts,
                          condition=condition,
                          is_nullgene=is_nullgene),
               res_longformat = res_longformat))
  
}
