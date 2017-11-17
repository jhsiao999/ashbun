library(ashbun)

df <- readRDS("tests/mouseengeltcell.rds")
dim(df)

counts <- round(df)

simdata <- simulationWrapper.sampleGenes(counts, Nsim = 2,
                Ngenes = 1000,
                Nsamples = 50,
                sample_method = "all_genes",
                sampleGenes_method = "high",
                pi0 = .9,
                beta_args = args.big_normal(betapi = 1,
                                            betamu = 0, betasd = .8))


output <- vector("list", length(simdata))
for (index in 1:length(output)) {
  output[[index]]  <- query.evaluation.simple(counts = simdata[[index]]$counts,
                                       condition = simdata[[index]]$condition,
                                       is_nullgene = simdata[[index]]$is_nullgene,
                                       methodsMeanExpression = c("DESeq2", "limmaVoom",
                                                                 "MAST", "BPSC"),
                                       report.control = list(fdr_control_threshold = .05),
                                       nsim = index)
}

