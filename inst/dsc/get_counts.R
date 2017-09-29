# HumanTungiPSC (UMI)
# Take the raw counts stored in singleCellRNASeqHumanTungiPSC
# and subset samples belong to two conditions
# Also, exclude ERCC genes
library(singleCellRNASeqHumanTungiPSC)
eset <- get(data("HumanTungiPSC"))
samplesToInclude <- which(Biobase::pData(eset)$individual == "NA19239" & Biobase::pData(eset)$replicate != "r3" )
genesToInclude <- grep("ERCC", rownames(Biobase::exprs(eset)), invert = TRUE)
ipsc_eset <- HumanTungiPSC[genesToInclude, samplesToInclude]
counts <- Biobase::exprs(ipsc_eset)[sample(nrow(Biobase::exprs(ipsc_eset)), ), ]
