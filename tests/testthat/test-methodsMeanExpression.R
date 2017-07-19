
test_that("Check input format", {
  ipsc_eset <- get(load(system.file("testdata", "HumanTungiPSC.rda", package = "ashbun")))
  count_matrix <- exprs(ipsc_eset)[sample(nrow(exprs(ipsc_eset)), 1000), ]
  condition <- pData(ipsc_eset)$replicate
  
  methods <- c("methodWrapper.DESeq2",
               "methodWrapper.edgeR",
               "methodWrapper.limmaVoom",
               "methodWrapper.bpsc",
               "methodWrapper.mast",
               "methodWrapper.rots")
"methodWrapper.scde", 

  for (index in 1:length(methods)) {
    output <- do.call(methods[index], list(count_matrix, condition))
    source("common-methodsMeanExpression.R", local = TRUE)
  }

})






