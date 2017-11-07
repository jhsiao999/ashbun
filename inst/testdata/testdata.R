# This script documents the code used to generate test datasets

# HumanTungiPSC (UMI)
# Take the raw counts stored in singleCellRNASeqHumanTungiPSC
# and subset samples belong to two conditions 
# Also, exclude ERCC genes
library(singleCellRNASeqHumanTungiPSC)
eset <- get(data("HumanTungiPSC"))
samplesToInclude <- which( pData(eset)$individual == "NA19239" & pData(eset)$replicate != "r3" )
genesToInclude <- grep("ERCC", rownames(exprs(eset)), invert = TRUE)

HumanTungiPSC <- HumanTungiPSC[genesToInclude, samplesToInclude] 

save(HumanTungiPSC, file = "inst/testdata/HumanTungiPSC.rda")

