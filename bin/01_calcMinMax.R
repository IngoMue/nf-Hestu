#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

bamlist <- read.table(args[1], header = T, sep = "\t")

bamlist$minDoC <- round(bamlist$DoC/3, 0)
bamlist$maxDoC <- round(bamlist$DoC*2, 0)

write.table(bamlist, "bams_DoCs.list", row.names = F, quote = F, sep = "\t")
