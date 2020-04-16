#!/usr/bin/env Rscript

library(DropletUtils)

cl <- commandArgs(trailingOnly = TRUE)

indir <- cl[1]
outdir <- cl[2]

mat <- read10xCounts(indir)

colnames(mat) <- colData(mat)$Barcode
rownames(mat) <- rowData(mat)$ID

if (all(order(colnames(mat)) == 1:ncol(mat))){
    print("Columns didn't need reordering")
}
if (all(order(rownames(mat)) == 1:nrow(mat))){
    print("Rows didn't need reordering")
}

mat <- mat[order(rowData(mat)$ID), order(colData(mat)$Barcode)]
write10xCounts(outdir, assays(mat)[[1]]) 
