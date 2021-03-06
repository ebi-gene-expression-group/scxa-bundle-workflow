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
    c("-y", "--meta-markers-dir"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Directory path in which to find metadata marker files as downloaded from a Scanpy workflow with scanpy-scripts"
  ),
  make_option(
    c("-n", "--meta-markers-file-pattern"),
    action = "store",
    default = "markers_.*.tsv",
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
    c("-s", "--select-top"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "Select the top n markers in each clustering"
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

# Argument checking

compulsory = c('counts_dir', 'clusters_file', 'cluster_markers_dir', 'output_file', 'cluster_markers_file_pattern')
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
all_marker_files <- unlist(lapply(c('cluster', 'meta'), function(markerType){
    
    if (! is.na(opt[[paste0(markerType, '_markers_dir')]])){
        marker_files <- list.files(path = paste0(markerType, '_markers'), pattern = opt[[paste0(markerType, '_markers_file_pattern')]], full.names = TRUE)
        
        if (length(marker_files) > 0){

            marker_cell_groupings <- gsub('_', ' ', sub('markers_(.+).tsv', '\\1', basename(marker_files)))

            # Load metadata groupings and subset to those relevant in markers 
            
            if (markerType == 'meta'){
                if ((! is.na(opt$cellgroups_file))){
                  cellgroups <- fread(opt$cellgroups_file)
                  colnames(cellgroups) <- gsub('_', ' ', colnames(cellgroups))
                  for (cg in names(cellgroups)[-1]){
                    cellgroups[[cg]] <- sub('^$', 'None', cellgroups[[cg]])
                  }
                  if (any(! marker_cell_groupings %in% colnames(cellgroups))){
                    write(paste("Can't find cell groups annotations for all of fields", paste(marker_cell_groupings, collapse=', '), 'available fields are:', paste(colnames(cellgroups), collapse=',')), stderr())
                    q(status = 1)
                  }

                  cols <- c('id', marker_cell_groupings)
                  clusters <- merge(clusters, cellgroups[, ..cols], by.x = 'row.names', by.y = 'id', all.x = TRUE, sort = FALSE)
                  rownames(clusters) <- clusters$Row.names
                  clusters <<- clusters[,-1]
                }
            }

            missing <- marker_cell_groupings[! marker_cell_groupings %in% colnames(clusters)]
            if (length(missing) > 0){
                write(paste('Marker file for', paste(missing, collapse=', '), 'missing its cell groups/ clusters'), stderr())
                q(status = 1)
            }

            marker_files <- structure(marker_files, names = marker_cell_groupings)
            if (markerType == 'cluster'){
                marker_files <- marker_files[order(as.numeric(names(marker_files)))]
            }
            marker_files
        }else{
            NULL
        }
    }else{
        NULL
    }
}))

if (length(all_marker_files) == 0){
    write('Found no valid marker files', stderr())
    q(status = 1) 
}

cluster_markers <- do.call(rbind, lapply(names(all_marker_files), function(x) cbind(variable = x, fread(all_marker_files[x], select = c('cluster', 'genes', 'pvals_adj'), colClasses = c('character', 'character', 'integer', 'character', 'numeric', 'numeric', 'numeric', 'numeric')))))
cluster_markers$cluster <- sub('^nan$', 'None', cluster_markers$cluster)

# Rank by padj within each clustering (and filter)

if (! is.na(opt$select_top)){
  cluster_markers <- do.call(rbind, lapply(split(cluster_markers, paste(cluster_markers$variable, cluster_markers$cluster)), function(x){
    x[order(x$pvals_adj) <= opt$select_top, ]
  }))
}

# Do some column renaming

cluster_markers <- cluster_markers[,c('genes', 'variable', 'cluster','pvals_adj')]
colnames(cluster_markers) <- c('gene_id', 'grouping_where_marker', 'group_where_marker', 'marker_p_value')

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
  missing <- cluster_markers$key[! cluster_markers$key %in%  cluster_stats$key]
  write(paste("Not all clusters and genes have summary stats, missing:",paste(missing, collapse=',')), stderr())
  q(status = 1)
}

print("Merging summary stats and marker info...")
final <- merge(cluster_stats, cluster_markers[,c('key', 'group_where_marker', 'marker_p_value')], by = 'key', all.x = TRUE, allow.cartesian=TRUE)

# Write output

print("Writing output...")
fwrite(final[,c('gene_id', 'grouping_where_marker', 'group_where_marker', 'cluster_id', 'marker_p_value', 'mean_expression', 'median_expression')], file = opt$output_file, sep=",", quote = TRUE)
