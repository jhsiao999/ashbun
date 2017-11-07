test_that("check to see if method runs w/o erros", {
  library(Biobase)
  ipsc_eset <- get(load(system.file("testdata", "HumanTungiPSC.rda", package = "ashbun")))
  counts <- exprs(ipsc_eset)[sample(nrow(exprs(ipsc_eset)), 500), ]
  condition <- pData(ipsc_eset)$replicate
  
  # methods that don't use condition information
  methods_1 <- c("normalize.tmm",
                 "normalize.rle",
                 "normalize.census",
                 "normalize.scran")
  
  # methods that perform normalization within group
  methods_2 <- c("normalize.scnorm")
  
  for (index in 1:length(methods_1)) {
    # check if libsize_factors is null
    output <- do.call(methods_1[index], list(counts))
    libsize_factors <- libsize_factors
    cpm <- normalize.cpm(counts, libsize_factors)$cpm
    
    source("common-methodsNormalize.R", local = TRUE)
  }
  
  for (index in 1:length(methods_2)) {
    # check if libsize_factors is null
    output <- do.call(methods_2[index], list(counts, condition))
    cpm <- output$counts_normed

    source("common-methodsNormalize.R", local = TRUE)
  }
  
} )




