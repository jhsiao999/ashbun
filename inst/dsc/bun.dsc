# RNA-seq analysis
# Joice and Gao (c) 2017

#
# Data
#

get_counts:
  exec: get_counts.R
  return:
    counts

simulate:
  exec: simulate.R
  params:
    counts: $counts
    Nsample: 100
    Ngenes: 1000
    sample_method: all_genes, per_genes
    pi0: 0.5, 0.9
    betapi: 1
    betamu: 0
    betasd: 0.8
    beta_function: args.big_normal
  return:
    counts = R(output$counts),
    condition = R(output$condition),
    null_gene = R(output$is_nullgene)

#
# Normalization
# this is a "lazy" style that I adopted from ashbun codes
# Can be made more concise / flexible
#

rle:
  exec: normalize.R
  .alias: rle
  params:
    counts: $counts
    condition: NULL
    method: rle
  return:
    model = R(output$model_output),
    libsize = R(output$libsize_factors)

tmm(rle):
  .alias: tmm
  params:
    method: tmm

cpm(rle):
  .alias: cpm
  params:
    method: cpm
  return:
    counts = R(output$cpm),
    model = R(output$model_output),
    libsize = R(output$libsize_factors)

census(rle):
  .alias: census
  params:
    method: census
  return:
    counts = R(output$counts_normed),
    model = R(output$model_output),
    libsize = R(output$libsize_factors)

scnorm(rle):
  .alias: scnorm
  params:
    method: scnorm
    condition: $condition
  return:
    counts = R(output$counts_normed),
    model = R(output$model_output),
    libsize = R(output$libsize_factors)

scran(rle):
  .alias: scran
  params:
    method: scran

#
# Mean expression
#
rots:
  exec: mean.R
  .alias: rots
  params:
    counts: $counts
    condition: $condition
    libsize: NULL
    method: rots
  return:
    pvalue = R(output$pvalue),
    fit = R(output$fit)

bpsc(rots):
  .alias: bpsc
  params:
    method: bpsc

mast(rots):
  .alias: mast
  params:
    method: mast

scde(rots):
  .alias: scde
  params:
    method: scde

limmaVoom(rots):
  .alias: limmaVoom
  params:
    method: limmaVoom

DESeq2(rots):
  .alias: DESeq2
  params:
    method: DESeq2
    libsize: $libsize

edgeR(rots):
  .alias: edgeR
  params:
    method: edgeR
    libsize: $libsize

DSC:
  run:
    all: get_counts * simulate * (cpm, tmm, rle, census, scnorm) * (DESeq2, edgeR, limmaVoom, bpsc, mast, rots, scde)
  R_libs: jhsiao999/singleCellRNASeqHumanTungiPSC, Biobase, assertthat,
          rhondabacher/SCnorm, scater, scran, edgeR, DESeq2, monocle, MAST, nghiavtr/BPSC
  lib_path: ../../R