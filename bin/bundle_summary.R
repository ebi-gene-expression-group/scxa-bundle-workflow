#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-c", "--counts-dir"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to old-school 10X dir with barcodes.tsv, genes.tsv and matrix.mtx, for count data"
  ),
  make_option(
    c("-e", "--experiment-id"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Atlas experiment ID these markers were derived for"
  ),
  make_option(
    c("-t", "--tpm-dir"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to old-school 10X dir with barcodes.tsv, genes.tsv and matrix.mtx, for tpm data"
  ),
  make_option(
    c("-j", "--clusters-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File containing cluster definitions as downloaded from a Scanpy workflow with scanpy-scripts"
  ),
  make_option(
    c("-d", "--cluster-markers-dir"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Directory path in which to find marker files as downloaded from a Scanpy workflow with scanpy-scripts"
  ),
  make_option(
    c("-m", "--cluster-markers-file-pattern"),
    action = "store",
    default = "markers_[0-9]*.tsv",
    type = 'character',
    help = "A pattern used to match marker files as downloaded from a Scanpy workflow with scanpy-scripts"
  ),
  make_option(
    c("-g", "--cellgroups-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A tabular file in which additional cell groups column can be found"
  ),
  make_option(
    c("-f", "--celltype-fields"),
    action = "store",
    default = 'inferred_cell_type',
    type = 'character',
    help = "Field name(s) (comma separated if multiple) in the cell type file from which cell type groupings can be derived"
  ),
  make_option(
    c("-y", "--celltype-markers-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A file containing cell type marker statistics, to be included with the processing of the cluster markers"
  ),
  make_option(
    c("-o", "--output-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to CSV file for output"
  )
)

opt <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)
saveRDS(opt, file = "opt.rds")

# Argument checking

compulsory = c('counts_dir', 'experiment_id', 'clusters_file', 'cluster_markers_dir', 'output_file', 'cluster_markers_file_pattern')
for (c in compulsory){
  if(is.na(opt[[c]])){
    write(paste('Required argument', c, 'missing'), stderr())
    q(status = 1)
  }
}

suppressPackageStartupMessages(require(DropletUtils))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(R.utils))
suppressPackageStartupMessages(require(sparseMatrixStats))
suppressPackageStartupMessages(require(tidyr))

# Parse the cluster definitions

clusters <- fread(opt$clusters_file, check.names=FALSE)
ks <- clusters$K
clusters <- t(clusters[,c(-1,-2)])
colnames(clusters) <- as.character(ks)

# Add other cell groupings, where provided

if ((! is.na(opt$celltype_markers_file)) && (! is.na(opt$cellgroups_file)) && (! is.na(opt$celltype_fields))){
  cellgroups <- fread(opt$cellgroups_file, select = c('id', unlist(strsplit(opt$celltype_fields, ','))))
  for (cg in names(cellgroups)[-1]){
    cellgroups[[cg]] <- sub('^$', 'None', cellgroups[[cg]])
  }
  
  clusters <- merge(clusters, cellgroups, by.x = 'row.names', by.y = 'id', all.x = TRUE, sort = FALSE)
  rownames(clusters) <- clusters$Row.names
  clusters <- clusters[,-1]
}

# Read the count-based expression values

# (first check if we need to gunzip)
print("Checking expression data...")
matfiles <- file.path(opt$counts_dir, c('barcodes.tsv', 'genes.tsv', 'matrix.mtx'))
gzipped_matfiles <- paste0(matfiles, '.gz')

for (i in 1:length(matfiles)){
  if (! file.exists(matfiles[i])){
    if (file.exists(gzipped_matfiles[i])){
      gunzip(gzipped_matfiles[i])
    }else{
      write(paste(matfiles[i], 'missing'), stderr())
      q(status = 1)
    }
  }
}

# ... then do the actual read

print("Loading expression data...")
counts <- read10xCounts(samples = opt$counts_dir)
colnames(counts) <- colData(counts)$Barcode

# Now read in the markers

print("Loading markers...")
cluster_marker_files <- list.files(path = opt$cluster_markers_dir, pattern = basename(opt$cluster_markers_file_pattern), full.names = TRUE)

# Order cluster markers by k 

k_vals <- sub('markers_([0-9]+).tsv', '\\1', basename(cluster_marker_files))
marker_files <- structure(cluster_marker_files[order(as.numeric(k_vals))], names = sort(as.numeric(k_vals)))

# Add in cell markers where provided

if (! is.na(opt$celltype_markers_file)){
  print("Cell type markers present...")
  marker_files['inferred_cell_type'] <- opt$celltype_markers_file
}

cluster_markers <- do.call(rbind, lapply(names(marker_files), function(x) cbind(exp_id = opt$experiment_id, variable = x, fread(marker_files[x], select = c('cluster', 'genes', 'pvals_adj'), colClasses = c('character', 'character', 'integer', 'character', 'numeric', 'numeric', 'numeric', 'numeric')))))
cluster_markers$cluster <- sub('^nan$', 'None', cluster_markers$cluster)

# We're ultimately populating a table like:
#
# Materialized view "atlasprd3.scxa_marker_gene_stats"
# Column          |          Type          | Collation | Nullable | Default 
# -------------------------+------------------------+-----------+----------+---------
#   experiment_accession    | character varying(255) |           |          | 
#   gene_id                 | character varying(255) |           |          | 
#   grouping_where_marker          | integer                |           |          | 
#   group_where_marker | integer                |           |          | 
#   cluster_id              | integer                |           |          | 
#   marker_p_value          | double precision       |           |          | 
#   mean_expression         | double precision       |           |          | 
#   median_expression       | double precision       |           |          | 
#
# ... so do some column renaming

cluster_markers <- cluster_markers[,c('exp_id',  'genes', 'variable', 'cluster','pvals_adj')]
colnames(cluster_markers) <- c('experiment_accession', 'gene_id', 'grouping_where_marker', 'group_where_marker', 'marker_p_value')

# Remove clusterings with no markers

clusters <- clusters[,colnames(clusters) %in% cluster_markers$grouping_where_marker, drop = FALSE]

# Match the matrix order to the clustering file, and remove any genes not present in the markers for speed

counts <- counts[rownames(counts) %in% cluster_markers$gene_id, match(rownames(clusters), colnames(counts))]

# Now we need mean and median expression values per k and cluster, using the sparseMatrixStats package

print("Calculating cluster-wise summary stats...")
cluster_stats <- do.call(rbind, lapply(colnames(clusters), function(x){

  # We're only interested in genes that were actually markers at this K
  k_genes <- unique(subset(cluster_markers, grouping_where_marker == x)$gene_id)
    
  cells_by_cluster <- split(rownames(clusters), factor(clusters[,x]))
    
  do.call(rbind, lapply(
    names(cells_by_cluster), 
    function(y){
      data.table(
        experiment_accession = opt$experiment_id,
        gene_id = k_genes,
        grouping_where_marker = x,
        cluster_id = y,
        mean_expression = structure(rowMeans2(assays(counts[k_genes,cells_by_cluster[[y]]])[[1]]), names = k_genes),
        median_expression = structure(rowMedians(assays(counts[k_genes,cells_by_cluster[[y]]])[[1]]), names = k_genes)
      )
    }
  ))
}))
  
# Now we have mean and median for each gene in each cluster in each k where it's
# a marker. But if it's a marker in multiple of those, then we have to duplicate
# the row accordingly

cluster_stats$key <- paste(cluster_stats$grouping_where_marker, cluster_stats$gene_id)
cluster_markers$key <- paste(cluster_markers$grouping_where_marker, cluster_markers$gene_id)

if (! all(cluster_markers$key %in% cluster_stats$key)){
  write("Not all clusters and genes have summary stats", stderr())
  q(status = 1)
}

print("Merging summary stats and marker info...")
final <- merge(cluster_stats, cluster_markers[,c('key', 'group_where_marker', 'marker_p_value')], by = 'key', all.x = TRUE, allow.cartesian=TRUE)

# Write output

print("Writing output...")
fwrite(final[,c('experiment_accession', 'gene_id', 'grouping_where_marker', 'group_where_marker', 'cluster_id', 'marker_p_value', 'mean_expression', 'median_expression')], file = opt$output_file, sep=",", quote = TRUE)
