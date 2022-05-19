#!/usr/bin/env Rscript

# Read command line arguments
args = commandArgs(trailingOnly=TRUE)

# load libraries
library(fastbaps)
library(ape)

# import fasta alignment
sparse.data <- import_fasta_sparse_nt(args[1])
multi <- multi_res_baps(sparse.data)

# write to output file
write.csv(multi, args[2], quote = FALSE)

