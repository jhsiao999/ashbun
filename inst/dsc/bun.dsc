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

cpm:
  exec: normalize.R
  .alias: cpm
  params:
    counts: $counts
    condition: NULL
    method: cpm
  return:
    counts = R(output$counts_normed),
    model = R(output$model_output),
    libsize = R(output$libsize_factors)

tmm(cpm):
  .alias: tmm
  params:
    method: tmm

rle(cpm):
  .alias: rle
  params:
    method: rle

census(cpm):
  .alias: census
  params:
    method: census

scnorm(cpm):
  .alias: scnorm
  params:
    method: scnorm
    condition: $condition

#
# Mean expression
#
rots:
  exec: mean.R
  .alias: rots
  params:
    counts: $counts
    condition:  $condition
    libsize: NULL
    method: rots
  return:
    pvalue = R(output$pvalue)
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
          rhondabacher/SCnorm, scater, scran, edgeR, DESeq, monocle, MAST
  lib_path: ../../R