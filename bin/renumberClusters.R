#!/usr/bin/env Rscript

library(data.table)

cl <- commandArgs(trailingOnly = TRUE)

infile <- cl[1]
outfile <- cl[2]

clusters <- fread(infile, check.names=FALSE)

if (min(clusters[,c(-1,-2)]) == 0){
    clusters <- cbind(clusters[,c(1,2)], clusters[,c(-1,-2)]+1)
}

fwrite(clusters, file=outfile, sep="\t", quote=FALSE, row.names = FALSE)
