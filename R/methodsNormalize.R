# Code for all evalulated between-sample normalization methods,
# including ones originally designed for bulk RNA-seq and for scRNA-seq



#--------------- bulk RNA-seq

#' @title Counts per million
#'
#' @description Compute CPM using the scale factors for library sizes.
#'
#' @param counts gene by sample expression count matrix (G by N).
#' @param libsize_factors numreic vector of the scale factors for library sizes.
#'                        Default to 1 - no adjustment.
#'
#' @param control List with control arguments, including
#'
#' @return
#'    \code{cpm} numeric matrix of counts per million.
#'
#' @export
normalize.cpm <- function(counts, libsize_factors = rep(1, dim(counts)[2]),
                          control = list(NULL)) {
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(counts))

  #<--------------------------------------
  libsize <- colSums(counts)
  libsize_factors <- rep(1, dim(counts)[2])
  libsize_normed <- libsize*libsize_factors
  cpm <- t((t(counts)/libsize_normed)*(10^6))

  return(list(cpm = cpm))
}




#' @title Library size normalization without adjustment
#'
#' @description Scale factors equal 1 for all libraries
#'
#' @param counts gene by sample expression count matrix (G by N).
#' @param libsize_factors numreic vector of the scale factors for library sizes.
#'                        Default to 1 - no adjustment.
#'
#' @param control List with control arguments, including
#'
#' @return
#'    \code{libsize_factors} numeric vector of the scale factors for library size.
#'
#' @export
normalize.lib <- function(counts, libsize_factors = rep(1, dim(counts)[2]),
                          control = list(NULL)) {
  return(list(libsize_factors = rep(1,ncol(counts))))
}



#' @title TMM
#'
#' @param counts Gene by sample expression count matrix (G by N).
#'
#' @return
#'    \code{libsize_factors} numeric vector of the scale factors for library size.
#'
#' @export
normalize.tmm <- function(counts, control = list(NULL)) {
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(counts))

  #<--------------------------------------
  libsize_factors <- edgeR::calcNormFactors(counts, method = "TMM")

  return(list(libsize_factors = libsize_factors))
  }




#' @title RLE (relative log expression)
#'
#' @description DESeq was used to implement RLE. Note that we chose to use the implementation by
#'              the developer of the method. edgeR also implments a version of RLE, which gives
#'              different results: (calcNormFactors(method = "RLE"))
#' @param counts Gene by sample expression count matrix (G by N).
#'
#' @return
#'    \code{libsize_factors} numeric vector of the scale factors for library size.
#'
#' @export
normalize.rle <- function(counts, control = list(NULL)) {
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(counts))

  #<--------------------------------------
  libsize_factors <- DESeq2::estimateSizeFactorsForMatrix(counts)

  return(list(libsize_factors = libsize_factors))
}


#--------------- single cell RNA-seq

#' @title census
#'
#' @description census is implemented in Monocle.
#'
#' @param counts Gene by sample expression count matrix (G by N). Has to be TPM or FPKM.
#' @param control List with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete census output.
#'
#' @return List of the following objects
#'    \code{count_normed} numeric matrix of normalized count matrix, or also know as true
#'                        abundance by monocle.
#'    \code{libsize_factors} numeric matrix of scale factors for library sizes.
#'    \code{model_output} monocle complete output.
#'
#' @export
normalize.census <- function(counts,
                             control = list(save_modelFit = FALSE)) {
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(counts))

  #  assertthat::assert_that(dim(counts)[2] == length(condition))

  # convert condition to factor
  #  if (!is.factor(condition)) {condition <- factor(condition)}

  #<--------------------------------------
  # construct Monocle object
  suppressPackageStartupMessages(library(monocle))
  # phenoData <- data.frame(sampleID = paste0("sample_", c(1:length(condition))),
  #                         condition = condition)
  phenoData <- data.frame(sampleID = paste0("sample_", c(1:dim(counts)[2])))
  rownames(phenoData) <- phenoData$sampleID
  colnames(counts) <- phenoData$sampleID


  featureData <- data.frame(geneID = rownames(counts),
                            row.names = rownames(counts))

  pd <- new("AnnotatedDataFrame", data = phenoData)
  fd <- new("AnnotatedDataFrame", data = featureData)

  # First create a CellDataSet from the relative expression levels
  cds <- newCellDataSet(as.matrix(counts),
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit=0.1,
                        expressionFamily=tobit(Lower=0.1))

  # Next, use it to estimate RNA counts
  counts_normed <- relative2abs(cds)

  # Now, make a new CellDataSet using the RNA counts
  cds_new <- newCellDataSet(as(as.matrix(counts_normed), "sparseMatrix"),
                            phenoData = pd,
                            featureData = fd,
                            lowerDetectionLimit=0.5,
                            expressionFamily=negbinomial.size())

  # estimate size factors and extract size factors
  libsize_factors <- sizeFactors(estimateSizeFactors(cds_new))

  # if save_modelFit, then output will include the original model fit
  if (control$save_modelFit) {
    model_output <- cds_new
  } else {
    model_output <- NULL
  }

  return(list(counts_normed = counts_normed,
              libsize_factors = libsize_factors,
              model_output = model_output))
}




#' @title SCnorm 1.1.0
#'
#' @description The method forces to output diagnoistic plots. The PLOT = T/F option is obsolete.
#'
#' @param counts Gene by sample expression count matrix (G by N).
#' @param condition Binary vector of length N indicating sample biological condition.
#' @param control A list with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete SCnorm output.
#'   \code{FilterCellNum} SCnorm argument. Minimum number of cells required to
#'                        estimate the count-depth relationship.
#'   \code{K} SCnorm argument. Number of quantile groups.
#'            Default NULL: letting SCnorm find the optimal K.
#'   \code{Ncores} SCnorm argument. Number of cores to use. Default 4.
#'   \code{useSpikes} SCnorm argument. Whether to use spike-ins to perform the between-condition
#'                    scaling using spike-ins. Assume that the spike-in labels start with ERCC
#'   \code{reportSF} SCnorm argument. Whehter to provide a matrix of scaling factors in the output.
#'
#' @return List of the following objects
#'    \code{count_normed} numeric vector of the scale factors for library size.
#'    \code{model_output} SCnorm complete output.
#'
#' @export
normalize.scnorm <- function(counts, condition,
                             control = list(save_modelFit = FALSE,
                                            FilterCellNum = 10,
                                            K = NULL,
                                            NCores = 4,
                                            useSpikes = FALSE,
                                            reportSF = TRUE)) {
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(counts))

  assertthat::assert_that(dim(counts)[2] == length(condition))

  # convert condition to factor
  if (!is.factor(condition)) {condition <- factor(condition)}

  #<--------------------------------------
  suppressPackageStartupMessages(library(SummarizedExperiment))
  suppressPackageStartupMessages(library(SCnorm))
  scnorm_output <- suppressMessages(
                      SCnorm(Data = counts, Conditions = condition,
                                      PLOT = F,
                                      FilterExpression = 1,
                                      FilterCellNum = control$FilterCellNum,
                                      withinSample = NULL,
                                      NCores = control$NCores,
                                      K = control$K,
                                      useSpikes = control$useSpikes,
                                      reportSF = control$reportSF) )

  # if save_modelFit, then output will include the original model fit
  if (control$save_modelFit) {
    model_output <- scnorm_output
  } else {
    model_output <- NULL
  }

  return(list(counts_normed = counts_normed,
              libsize_factors = libsize_factors,
              model_output = model_output))
}





#' @title scran
#'
#' @description The most recent version of scran depends on R > 3.4.
#'              I use scran 1.2.2, which depeneds on R.3.3.0. One difference that I have found
#'              so far is about SCESet - in the older version, the spike-in variable is a binary
#'              vector of yes/no, and in the newer version, the spike-in variable is a numeric
#'              vector specifying which genes is spike-in controls.
#'
#' @param counts Gene by sample expression count matrix (G by N).
#' @param control List with control arguments, including
#'   \code{save_modelFit} TRUE to output the complete SCnorm output.
#'   \code{get_cell_clusters} TRUE/FALSE to use scran to cluster cells and then estimate
#'                            size factors by cell clusters.
#'   \code{min.size} scran internal argument. Integer scalar specifying the minimum size
#'                   of each cluster. Default to be 50.
#'
#' @return
#'    \code{libsize_factors} numeric vector of the scale factors for library size.
#'
#' @examples
#' ipsc_eset <- get(load(system.file("testdata", "HumanTungiPSC.rda", package = "ashbun")))
#' counts <- exprs(ipsc_eset)[sample(nrow(exprs(ipsc_eset)), ), ]
#'
#' #---- generat simulated datasets
#' simdata_list <- simulationWrapper(counts, Nsim = 5, Nsample = 100, Ngene = 500)
#'
#' #---- extract a single dataset as an example
#' #---- take pi0 = .9, the first simulated data
#' simdata <- simdata_list[[3]][[1]]
#'
#' #---- normalize
#' output <- normalize.scran(simdata$counts)
#'
#' @export
normalize.scran <- function(counts,
                             control = list(save_modelFit = FALSE,
                                            get_cell_clusters = FALSE,
                                            min.size = 50
                                            )) {
  library(scran)
  #--------------------------
  # Make sure input format is correct
  assertthat::assert_that(is.matrix(counts))

  #<--------------------------------------
  # construct SCE object
  suppressPackageStartupMessages(library(scater))
  phenoData <- data.frame(sampleID = paste0("sample_", c(1: dim(counts)[2])))
  rownames(phenoData) <- phenoData$sampleID
  colnames(counts) <- phenoData$sampleID

  if (!is.integer(counts)) {
    counts <- apply(counts, 2, function(x) { storage.mode(x) <- 'integer'; x})
  }

  sce <- SingleCellExperiment(list(counts=counts))
  sce <- computeSumFactors(sce)

  # cluster cells
  if (control$get_cell_clusters) {
    cell_clusters <- scran::quickCluster(sce, min.size = control$min.size)

    # scran requires the number of cells in each cluster should at leaset be
    # twice that of the larges pool size
#    control$sizes <- seq(10, round(min(table(cell_clusters))/2 ), by = 20)
    # re-compute size factors
    scran_output <- scran::computeSumFactors(sce,
                                             clusters = cell_clusters)
  } else {
    control$cell_clusters <- NULL
    scran_output <- sce
  }

  # extract size factors
  libsize_factors <- sizeFactors(scran_output)

  # normalization
  # dont' use their method for computing CPM, I don't know how they compute CPM
  # it seems that some sort of pseudocount procedure is used, but the code
  # is not easy to decipher...
  # scran_output <- scater::normalize(scran_output,
  #                                   recompute_cpm = TRUE)

  # if save_modelFit, then output will include the original model fit
  if (control$save_modelFit) {
    model_output <- scran_output
  } else {
    model_output <- NULL
  }
  #
  return(list(libsize_factors = libsize_factors,
              model_output = model_output))
}









#' #<-------------- Utility functions
#'
#' #' @title Get genomic features using biomaRT
#' #'
#' #' @description Use /code{biomaRt} package to get some gene-leve features, including
#' #'        transcript count and length (max & median)
#' #'
#' #' @param genome name of genome, can be hg19 (human), mm9 or mm10 (mice genome)
#' #' @param ensembl_gene_id Character vector of ensembl gene IDs.
#' #' @param entrez_gene_id Character vector of entrez gene IDs. Default NULL.
#' #' @param hgnc_symbol Character vector of hgnc symbol. Default NULL.
#' #'
#' #' @example
#' #' library(singleCellRNASeqHumanTungiPSC)
#' #' gene_names <- fData(HumanTungiPSC)$gene_name
#' #' gene_features <- getTranscriptInfo("hg19", ensembl_gene_id = gene_names)
#' #'
#' #' @return data.frame of
#' #'    \code{ensembl_gene_id}
#' #'    \code{transcript_count}
#' #'    \code{max_transript_length}
#' #'    \code{median_transcript_length}
#' #'
#' #' @export
#' #'
#' getTranscriptInfo <- function(genome,
#'                               ensembl_gene_id = NULL,
#'                               entrez_gene_id = NULL,
#'                               hgnc_symbol = NULL) {
#'   # # use GenomicFeatures to get more informatin on the transcriptome
#'   # # such as length of coding region, length of 5'UTR and 3'UTR
#'   # # gene_id in TxDb are Entrez IDs
#'   # if (genome == "hg19") {
#'   #   library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#'   #   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'   # }
#'   # if (genome == "mm9") {
#'   #   library("TxDb.Mmusculus.UCSC.mm9.knownGene")
#'   #   txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#'   # }
#'   # if (genome == "mm10") {
#'   #   library("TxDb.Mmusculus.UCSC.mm10.knownGene")
#'   #   txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#'   # }
#'   #
#'   # # find the length of non-overlapping exons
#'   # exon_list_per_gene <- GenomicFeatures::transcriptLengths(txdb,
#'   #                                                          with.cds_len = TRUE,
#'   #                                                          with.utr5_len = TRUE,
#'   #                                                          with.utr3_len = TRUE)
#'   # exon_list <- exonsBy(txdb, "gene")
#'   # exon_overlap_list <- exonsByOverlaps(txdb, exon_list_per_gene)
#'
#'   library(biomaRt)
#'   if (genome == "hg19") {
#'     ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
#'                        host = "grch37.ensembl.org",
#'                        dataset="hsapiens_gene_ensembl")
#'   }
#'   if (genome == "mm9") {
#'     ensembl <- useMart(biomart='ENSEMBL_MART_ENSEMBL',
#'                        host = 'may2012.archive.ensembl.org',
#'                        dataset = "mmusculus_gene_ensembl")
#'   }
#'   if (genome == "mm10") {
#'     ensembl <- useMart(biomart='ENSEMBL_MART_ENSEMBL',
#'                        host = 'jul2012.archive.ensembl.org',
#'                        dataset = "mmusculus_gene_ensembl")
#'   }
#'
#'   gene_methods <- list(ensembl_gene_id, entrez_gene_id, hgnc_symbol)
#'   method_lookup <- which(!is.null( gene_id_methods) )
#'
#'   # get attributes
#'   gene_info <- getBM(attributes = c("ensembl_gene_id",
#'                                     "hgnc_symbol",
#'                                     "entrezgene",
#'                                     "transcript_count",
#'                                     "transcript_length"),
#'                      filters = "ensembl_gene_id",
#'                      values = gene_methods[method_lookup],
#'                      mart = ensembl)
#'
#'   # compute summary information
#'   library(dplyr)
#'   gene_summary <- gene_info %>%
#'     group_by(ensembl_gene_id) %>%
#'     summarise(entrez_gene_id = first(entrezgene),
#'               hgnc_symbol = first(hgnc_symbol),
#'               transcript_count = first(transcript_count),
#'               max_transcript_length = max(transcript_length),
#'               median_transcript_length = median(transcript_length))
#'
#'   gene_summary <- as.data.frame(gene_summary)
#'
#'   return(gene_summary)
#' }
#'
#'
#'
#'
#'
#' #<-------------- Within-sample normalization
#'
#' #' @title RPKM (Reads Per Kilobase of transcript per Million mapped reads)
#' #'
#' #' @param counts Gene by sample expression count matrix (G by N).
#' #' @param gene_length Numeric vector of gene length.
#' #' @param control List with control arguments, including
#' #'
#' #' @return
#' #'    \code{rpkm} numeric matrix of
#' #'
#' #' @export
#' normalize.rpkm <- function(counts,
#'                            gene_length, control= list(NULL)) {
#'
#'   #--------------------------
#'   # Make sure input format is correct
#'   assertthat::assert_that(is.matrix(counts))
#'   assertthat::assert_that(is.integer(counts),
#'                           msg = "counts is not integer-values")
#'   assertthat::assert_that(length(gene_length) == nrow(counts))
#'
#'   #<--------------------------------------
#'   rpk <- counts/gene_length
#'   rpkm <- t(t(rpk)/colSums(rpk))*(10^6)
#'
#'   return(rpkm)
#' }



#######<------  TBD
#' #' @title BASiCS 0.7.28
#' #'
#' #' @description downloaed from https://github.com/catavallejos/BASiCS
#' #'
#' #' @param counts Gene by sample expression count matrix (G by N). Has to be RPKM or FPKM.
#' #' @param condition Binary vector of length N indicating sample biological condition.
#' #' @param control List with control arguments, including
#' #'   \code{save_modelFit} TRUE to output the complete SCnorm output.
#' #'
#' #' @return List of the following objects
#' #'    \code{count_normed} numeric matrix of normalized count matrix, or also know as true
#' #'                        abundance by monocle.
#' #'    \code{model_output} monocle complete output.
#' #'
#' #' @export
#' normalize.basics <- function(counts, condition,
#'                              control = list(save_modelFit = FALSE)) {
#'   #--------------------------
#'   # Make sure input format is correct
#'   assertthat::assert_that(is.matrix(counts))
#'   assertthat::assert_that(is.integer(counts),
#'                           msg = "counts is not integer-values")
#'
#'   assertthat::assert_that(dim(counts)[2] == length(condition))
#'
#'   # convert condition to factor
#'   if (!is.factor(condition)) {condition <- factor(condition)}
#'
#'   #<--------------------------------------
#'   # construct Monocle object
#'   library(BASiCS)
#'   bsd <- newBASiCS_Data(count_matrix, tech = spikes, SpikeInfo)
#'
#'
#'   return(list(counts_normed = counts_normed,
#'               libsize_factors = libsize_factors,
#'               model_output = model_output))
#' }




#######<------  TBD
# some scran code for using spike-in expression to compute normalization factor
# if (control$get.spikes) {
#     assertthat::assert_that(length(unique(spikes)) == 2)
#
#     scran::isSpike(sce) <- "spikes"
#
#     # scran option to normalize data by the scale factors computed from spike-ins
#     control$general.use <- TRUE
#
#     # compute normailzation factor based on spike-in
#     # scran doesn't do clustering based on spike-in?
#     scran_output <- scran::computeSpikeFactors(sce,
#                                                general.use = control$general.use)
# } else {
#     # first group cells into clusters of similar expression
#     if (control$get_cell_clusters) {
#       cell_clusters <- scran::quickCluster(sce,
#                                            min.size = control$min.size,
#                                            get.spikes = control$get.spikes)
#
#       # scran requires the number of cells in each cluster should at leaset be
#       # twice that of the larges pool size
#       control$sizes <- seq(10, round(max(table(cell_clusters))/2 ), by = 20)
#     } else {
#       control$cell_clusters <- NULL
#     }
#
#     # compute scale factor
#     scran_output <- scran::computeSumFactors(sce,
#                                              sizes = control$sizes,
#                                              clusters = cell_clusters)
# }
