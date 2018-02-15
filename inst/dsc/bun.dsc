# RNA-seq analysis
# Joice and Gao (c) 2017 - 2018

#
# Data
#

get_counts: get_counts.R
  $counts: counts

simulate: simulate.R
    counts: $counts
    Nsample: 100
    Ngenes: 1000
    sample_method: all_genes, per_genes
    pi0: 0.5, 0.9
    betapi: 1
    betamu: 0
    betasd: 0.8
    beta_function: args.big_normal
    counts: output$counts
    condition: output$condition
    null_gene: output$is_nullgene

#
# Normalization
# this is a "lazy" style that I adopted from ashbun codes
# Can be made more concise / flexible
#

rle: normalize.R
  counts: $counts
  condition: NULL
  method: rle
  $model: output$model_output
  $libsize: output$libsize_factors

tmm(rle):
  method: tmm

cpm(rle):
  method: cpm
  $counts: output$cpm

census(rle):
  method: census
  $counts: output$counts_normed

scnorm(rle):
  method: scnorm
  condition: $condition
  $counts: output$counts_normed

scran(rle):
  method: scran

#
# Mean expression
#
rots: mean.R
  counts: $counts
  condition: $condition
  libsize: NULL
  method: rots
  $pvalue: output$pvalue
  $fit: output$fit

bpsc(rots):
  method: bpsc

mast(rots):
  method: mast

scde(rots):
  method: scde

limmaVoom(rots):
  method: limmaVoom

DESeq2(rots):
  method: DESeq2
  libsize: $libsize

edgeR(DESeq2):
  method: edgeR

DSC:
  define:
    data: get_counts * simulate
    norm: cpm, tmm, rle, census, scnorm
    de: DESeq2, edgeR, limmaVoom, bpsc, mast, rots, scde
  run: data * norm * de
  R_libs: singleCellRNASeqHumanTungiPSC@jhsiao999/singleCellRNASeqHumanTungiPSC, Biobase, assertthat,
          SCnorm@rhondabacher/SCnorm, scater, scran, edgeR, DESeq2, monocle, MAST, BPSC@nghiavtr/BPSC
  lib_path: ../../R