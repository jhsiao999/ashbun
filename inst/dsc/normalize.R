source("methodsNormalize.R")
if (is.null(condition)) {
  output = do.call(paste0('normalize.', method), list(counts = counts))
} else {
  output = do.call(paste0('normalize.', method), list(counts = counts, condition = condition))
}
