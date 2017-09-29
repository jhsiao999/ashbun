source("methodsMeanExpression-bulk.R")
source("methodsMeanExpression-singleCell.R")
if (is.null(libsize)) {
  output = do.call(paste0('methodWrapper.', method), list(counts = counts, condition = condition))
} else {
  output = do.call(paste0('methodWrapper.', method), list(counts = counts, condition = condition, libsize_factor = libsize))
}
