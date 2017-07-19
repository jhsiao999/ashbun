# inclue this in the simulation code
# Check column names of log2counts; if repeated, then rename
if (sum(duplicated(colnames(count_matrix))) > 0) {
  colnames(count_matrix) <- paste0("sample.", c(1:ncol(object)))
}




fit <- scde::scde.expression.difference(o.ifm, cd, o.prior, 
                                        groups  =  condition, 
                                        n.randomizations  =  control$bootstrap_times, 
                                        # number of bootstrap randomizations to be performed
                                        n.cores  =  control$n_cores, verbose  =  0)


scde::scde.expression.difference
function (models = o.ifm, counts = count_matrix, prior = o.prior, groups = condition, batch = NULL, 
          n.randomizations = 150, n.cores = 4, batch.models = models, 
          return.posteriors = FALSE, expectation = 0, verbose = 0) 
{
  if (!all(rownames(models) %in% colnames(counts))) {
    stop("ERROR: provided count data does not cover all of the cells specified in the model matrix")
  }
  ci <- match(rownames(models), colnames(counts))
  counts <- as.matrix(counts[, ci])
  if (is.null(groups)) {
    groups <- as.factor(attr(models, "groups"))
    if (is.null(groups)) 
      stop("ERROR: groups factor is not provided, and models structure is lacking groups attribute")
    names(groups) <- rownames(models)
  }
  if (length(levels(groups)) != 2) {
    stop(paste("ERROR: wrong number of levels in the grouping factor (", 
               paste(levels(groups), collapse = " "), "), but must be two.", 
               sep = ""))
  }
  correct.batch <- FALSE
  if (!is.null(batch)) {
    if (length(levels(batch)) > 1) {
      correct.batch <- TRUE
    }
    else {
      if (verbose) {
        cat("WARNING: only one batch level detected. Nothing to correct for.")
      }
    }
  }
  if (correct.batch) {
    batch <- as.factor(batch)
    bgti <- table(groups, batch)
    bgti.ft <- fisher.test(bgti)
    if (verbose) {
      cat("controlling for batch effects. interaction:\n")
      print(bgti)
    }
    if (bgti.ft$p.value < 0.001) {
      cat("WARNING: strong interaction between groups and batches! Correction may be ineffective:\n")
      print(bgti.ft)
    }
    if (verbose) {
      cat("calculating batch posteriors\n")
    }
    batch.jpl <- tapply(seq_len(nrow(models)), groups, function(ii) {
      scde.posteriors(models = batch.models, counts = counts, 
                      prior = prior, batch = batch, composition = table(batch[ii]), 
                      n.cores = n.cores, n.randomizations = n.randomizations, 
                      return.individual.posteriors = FALSE)
    })
    if (verbose) {
      cat("calculating batch differences\n")
    }
    batch.bdiffp <- calculate.ratio.posterior(batch.jpl[[1]], 
                                              batch.jpl[[2]], prior, n.cores = n.cores)
    batch.bdiffp.rep <- quick.distribution.summary(batch.bdiffp)
  }
  else {
    if (verbose) {
      cat("comparing groups:\n")
      print(table(as.character(groups)))
    }
  }
  
  jpl <- tapply(seq_len(nrow(models)), groups, function(ii) {
    scde.posteriors(models = models[ii, , drop = FALSE], 
                    counts = counts[, ii, drop = FALSE], prior = prior, 
                    n.cores = n.cores, n.randomizations = n.randomizations)
  })
  
  
  if (verbose) {
    cat("calculating difference posterior\n")
  }
  bdiffp <- calculate.ratio.posterior(jpl[[1]], jpl[[2]], prior, 
                                      n.cores = n.cores)
  if (verbose) {
    cat("summarizing differences\n")
  }
  bdiffp.rep <- quick.distribution.summary(bdiffp, expectation = expectation)
  if (correct.batch) {
    if (verbose) {
      cat("adjusting for batch effects\n")
    }
    a.bdiffp <- calculate.ratio.posterior(bdiffp, batch.bdiffp, 
                                          prior = data.frame(x = as.numeric(colnames(bdiffp)), 
                                                             y = rep(1/ncol(bdiffp), ncol(bdiffp))), skip.prior.adjustment = TRUE, 
                                          n.cores = n.cores)
    a.bdiffp.rep <- quick.distribution.summary(a.bdiffp, 
                                               expectation = expectation)
    if (return.posteriors) {
      return(list(batch.adjusted = a.bdiffp.rep, results = bdiffp.rep, 
                  batch.effect = batch.bdiffp.rep, difference.posterior = bdiffp, 
                  batch.adjusted.difference.posterior = a.bdiffp, 
                  joint.posteriors = jpl))
    }
    else {
      return(list(batch.adjusted = a.bdiffp.rep, results = bdiffp.rep, 
                  batch.effect = batch.bdiffp.rep))
    }
  }
  else {
    if (return.posteriors) {
      return(list(results = bdiffp.rep, difference.posterior = bdiffp, 
                  joint.posteriors = jpl))
    }
    else {
      return(bdiffp.rep)
    }
  }
}

