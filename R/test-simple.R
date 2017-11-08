# 
# library(ashbun)
# 
# df <- readRDS("tests/humantungipsc.rds")
# simdata_list <- simulationWrapper(df, Nsim = 2,
#                                   Ngenes = 100,
#                                   Nsamples = 40,
#                                   sample_method = "all_genes",
#                                   pi0 = .9,
#                                   beta_args = args.big_normal(betapi = 1,
#                                                               betamu = 0, betasd = .8))
# 
# output <- vector("list", length(simdata_list))
# for (index in 1:length(output)) {
#   output[[index]]  <- query.evaluation.simple(counts = simdata_list[[index]]$counts,
#                                        condition = simdata_list[[index]]$condition,
#                                        is_nullgene = simdata_list[[index]]$is_nullgene,
#                                        methodsMeanExpression = c("DESeq2", "limmaVoom",
#                                                                  "MAST"),
#                                        report.control = list(fdr_control_threshold = .05),
#                                        nsim = index)
# }
# results.multiple.data <- saveRDS(output, file = "tests.results.multiple.data.rds")


