ipsc_eset <- get(load(system.file("testdata", "HumanTungiPSC.rda", package = "ashbun")))
counts <- exprs(ipsc_eset)[sample(nrow(exprs(ipsc_eset)), 500), ]
condition <- pData(ipsc_eset)$replicate

query.pipeline <- function(counts, condition, null,
                           methodsNormalize, methodsMeanExpression) {
  