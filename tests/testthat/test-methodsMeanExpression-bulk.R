
test_that("Check input format", {
  ipsc_eset <- get(load(system.file("testdata", "HumanTungiPSC.rda", package = "ashbun")))
  count_matrix <- exprs(ipsc_eset)
  condition <- pData(ipsc_eset)$replicate
  
  ee <- methodWrapper.edgeR(count_matrix, condition = condition)
  de <- methodWrapper.DESeq2(count_matrix, condition = condition)
  ll <- methodWrapper.limmaVoom(count_matrix, condition = condition)
  mm <- methodWrapper.mast(count_matrix, condition = condition)
  ss <- methodWrapper.scde(count_matrix, condition = condition)
  rr <- methodWrapper.rots(count_matrix, condition = condition)
  
})




# inclue this in the simulation code
# Check column names of log2counts; if repeated, then rename
if (sum(duplicated(colnames(count_matrix))) > 0) {
  colnames(count_matrix) <- paste0("sample.", c(1:ncol(object)))
}



