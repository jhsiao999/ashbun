
test_that("", {
  ispc_eset <- system.file("HumanTungiPSC", package = "ashbun")
  count_matrix <- exprs(ipsc_eset)
})




condition <- pData(eset_subset)$replicate
count_matrix <- exprs(eset_subset)




ed <- methodWrapper.edgeR(count_matrix, condition = condition)
de <- methodWrapper.DESeq2(count_matrix, condition = condition)
ll <- methodWrapper.limmaVoom(count_matrix, condition = condition)
mm <- methodWrapper.mast(count_matrix, condition = condition)
ss <- methodWrapper.scde(count_matrix, condition = condition)
rr <- methodWrapper.rots(count_matrix, condition = condition)

str(ll)
str(de)
str(ed)


# inclue this in the simulation code
# Check column names of log2counts; if repeated, then rename
if (sum(duplicated(colnames(count_matrix))) > 0) {
  colnames(count_matrix) <- paste0("sample.", c(1:ncol(object)))
}



