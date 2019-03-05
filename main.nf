#!/usr/bin/env nextflow

resultsRoot = params.resultsRoot

RAW_FILTERED_MATRIX = Channel.fromPath( "$resultsRoot/${params.rawFilteredMatrix}", checkIfExists: true)
NORMALISED_MATRIX = Channel.fromPath( "$resultsRoot/${params.normalisedMatrix}", checkIfExists: true)
RAW_TPM_MATRIX = Channel.fromPath( "$resultsRoot/${params.tpmMatrix}", checkIfExists: true)
SCANPY_CLUSTERS = Channel.fromPath( "$resultsRoot/${params.clusters}", checkIfExists: true)
SCANPY_TSNE = Channel.fromPath( "$resultsRoot/${params.tsneDir}/embeddings*.csv", checkIfExists: true )
SCANPY_MARKERS = Channel.fromPath( "$resultsRoot/${params.markersDir}/markers_*.csv", checkIfExists: true )
softwareTemplate = params.softwareTemplate

// Send channels to different processes

RAW_FILTERED_MATRIX.into{
    RAW_FILTERED_MATRIX_FOR_TPM_FILTERING
    RAW_FILTERED_MATRIX_FOR_MTX
    RAW_FILTERED_MATRIX_FOR_TSV
}

NORMALISED_MATRIX.into{
    NORMALISED_MATRIX_FOR_MTX
    NORMALISED_MATRIX_FOR_TSV
}

// We need a TPM matrix for the same cell and gene sets as the Scanpy-filtered
// count matrices. Filter the abundance matrix to match the count matrices

process filter_tpms {

    conda 'bioconductor-dropletutils zip'
    
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 10
    
    input:
        file(filteredCountsMatrix) from RAW_FILTERED_MATRIX_FOR_TPM_FILTERING
        file(tpmMatrix) from RAW_TPM_MATRIX

    output:
        set val("${tpmMatrix.getBaseName()}"), file("${tpmMatrix.getBaseName()}_filter_cells_genes/matrix.mtx") into TPM_FILTER_CELLS_MTX
        set val("${tpmMatrix.getBaseName()}"), file("${tpmMatrix.getBaseName()}_filter_cells_genes/genes.tsv") into TPM_FILTER_CELLS_MTX_ROWS
        set val("${tpmMatrix.getBaseName()}"), file("${tpmMatrix.getBaseName()}_filter_cells_genes/barcodes.tsv") into TPM_FILTER_CELLS_MTX_COLS

    """
        #!/usr/bin/env Rscript
        
        suppressPackageStartupMessages(require(DropletUtils))
   
        unzip('$filteredCountsMatrix')
        filt_counts_sce <- read10xCounts(sub('.zip', '', '$filteredCountsMatrix'))
        
        unzip('$tpmMatrix')
        tpm_sce <- read10xCounts(sub('.zip', '', '$tpmMatrix'))

        tpm_sce <- tpm_sce[rownames(filt_counts_sce), which(colData(tpm_sce)\$Barcode %in% colData(filt_counts_sce)\$Barcode ) ]
        write10xCounts(assays(tpm_sce)[[1]], path = '${tpmMatrix.getBaseName()}_filter_cells_genes', barcodes = colData(tpm_sce)\$Barcode, gene.id = rownames(tpm_sce))
    """
}

// Compress the mtx files into a zip. Couldn't get R to do this directly in the
// rule above

process compress_filtered_tpms {
    
    publishDir "$resultsRoot/bundle", mode: 'copy', overwrite: true
    
    input:
        set val(matName), file(mtx) from TPM_FILTER_CELLS_MTX
        set val(matName), file(genes) from TPM_FILTER_CELLS_MTX_ROWS
        set val(matName), file(barcodes) from TPM_FILTER_CELLS_MTX_COLS

    output:
        file("${matName}_filter_cells_genes.zip") into RAW_FILTERED_TPM_MATRIX

    """
        mkdir -p ${matName}_filter_cells_genes
        mv matrix.mtx genes.tsv barcodes.tsv ${matName}_filter_cells_genes 
        zip -r ${matName}_filter_cells_genes.zip ${matName}_filter_cells_genes
    """
}


// Report software versions

process make_software_report {

    publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true
    
    output:
        file "software.tsv" into SOFTWARE

    """
        generateSoftwareReport.sh ${softwareTemplate} software.tsv
    """
}

// Find out what perplexities are represented by the t-SNE files

process mark_perplexities {

    input:
        file tSNE from SCANPY_TSNE

    output:
        set stdout, file (tSNE) into EMBEDDINGS_BY_PERPLEXITY 

    """
       echo $tSNE | grep -o -E '[0-9]+' | tr '\\n' '\\0' 
    """
}

// Convert the t-SNE files to tsv

process tsne_to_tsv {
    
    publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true
    
    input:
        set val(perplexity), file(embeddings) from EMBEDDINGS_BY_PERPLEXITY

    output:
        set val(perplexity), file ("tsne_*.tsv") into TSV_EMBEDDINGS 

    """
    cat $embeddings | sed 's/,/\t/g' > tsne_${perplexity}.tsv
    """
}

// Combine the listing of t-SNE files for the manifest

process tsne_lines {

    input:
        set val(perplexity), file(embeddings) from TSV_EMBEDDINGS

    output:
        stdout TSNE_MANIFEST_LINES 

    """
    echo -e "tnse_embeddings\t${embeddings}\t$perplexity"
    """
}

// Transform the clusters into the format expected. Currently we need to change
// matrix orientation, and insert a logical initial column pointing at the
// default parameterisation (resolution = 1) 

process transform_clusters{

    conda 'r-base'
    
    publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true
        
    input:
        file clustersFile from SCANPY_CLUSTERS

    output:
        file 'clusters_for_bundle.txt' into BUNDLE_CLUSTERS

    """
        #!/usr/bin/env Rscript
        
        clusters <- t(read.delim("$clustersFile", row.names=1))
        resolutions <- sub('louvain_r', '', rownames(clusters))
        clusters <- cbind(sel.resolution = as.character(resolutions == '1.0'), resolution = resolutions, clusters)
        write.table(clusters, file="clusters_for_bundle.txt", sep='\t', quote=FALSE, row.names= FALSE)
    """        
}

// Repackage the matrices 


Channel.from( 'raw_filtered', 'normalised', 'tpm_filtered' ).into{
    EXPRESSION_TYPES_FOR_MTX
    EXPRESSION_TYPES_FOR_TSV
}

RAW_FILTERED_TPM_MATRIX.into{
    RAW_FILTERED_TPM_MATRIX_FOR_MTX
    RAW_FILTERED_TPM_MATRIX_FOR_TSV
}

process repackage_matrices {

    publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true
    
    input:
        file expressionMatrix from RAW_FILTERED_MATRIX_FOR_MTX.concat(NORMALISED_MATRIX_FOR_MTX).concat(RAW_FILTERED_TPM_MATRIX_FOR_MTX)
        val expressionType from EXPRESSION_TYPES_FOR_MTX

    output:
        set val(expressionType), file("${expressionType}/genes.tsv.gz") into MTX_MATRIX_ROWNAMES
        set val(expressionType), file("${expressionType}/barcodes.tsv.gz") into MTX_MATRIX_COLNAMES
        set val(expressionType), file("${expressionType}/matrix.mtx.gz") into MTX_MATRIX_CONTENT


    """
        zipdir=\$(unzip -qql ${expressionMatrix.getBaseName()}.zip | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/||')
        unzip ${expressionMatrix.getBaseName()}        
    
        mkdir -p ${expressionType}
        mv \${zipdir}/matrix.mtx ${expressionType} && gzip ${expressionType}/matrix.mtx
        mv \${zipdir}/genes.tsv ${expressionType} && gzip ${expressionType}/genes.tsv
        mv \${zipdir}/barcodes.tsv ${expressionType} && gzip ${expressionType}/barcodes.tsv
    """        

}

// Make tsv-format matrices

process mtx_to_tsv {
    
    conda 'bioconductor-dropletutils r-data.table'

    publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true
    
    memory { 5.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        file expressionMatrix from RAW_FILTERED_MATRIX_FOR_TSV.concat(NORMALISED_MATRIX_FOR_TSV).concat(RAW_FILTERED_TPM_MATRIX_FOR_TSV)
        val expressionType from EXPRESSION_TYPES_FOR_TSV
        
    output:
        set val(expressionType), file("${expressionType}/${expressionType}.tsv") into TSV_MATRICES

    """
        #!/usr/bin/env Rscript
        
        suppressPackageStartupMessages(require(DropletUtils))
        suppressPackageStartupMessages(require(data.table))
        source(file.path(Sys.getenv(c("SCRIPTS_DIR")), "utils.R"))    
   
        unzip('$expressionMatrix')
        sce <- read10xCounts(sub('.zip', '', '$expressionMatrix'))
        colnames(sce) <- colData(sce)\$Barcode 

        dir.create('${expressionType}')
        write.tsv(as.data.frame(cbind(Feature = rownames(sce), as.matrix(assays(sce)[[1]])), col.names = c('Feature', colData(sce)\$Barcode)), "${expressionType}/${expressionType}.tsv")      
    """
}

// Make manifest lines for matrices

process matrix_lines {

    input:
        set val(expressionType), file(matrixRows) from MTX_MATRIX_ROWNAMES
        set val(expressionType), file(matrixCols) from MTX_MATRIX_COLNAMES
        set val(expressionType), file(matrixContent) from MTX_MATRIX_CONTENT
        set val(expressionType), file(tsvMatrix) from TSV_MATRICES

    output:
        stdout MATRIX_MANIFEST_LINES 

    """
    echo -e "mtx_matrix_rows\t$matrixRows\t$expressionType"
    echo -e "mtx_matrix_cols\t$matrixCols\t$expressionType"
    echo -e "mtx_matrix_content\t$matrixContent\t$expressionType"
    echo -e "tsv_matrix\t$tsvMatrix\t$expressionType"
    """
}

// Find out what resolutions are represented by the marker files

process mark_marker_resolutions {

    input:
        file markersFile from SCANPY_MARKERS

    output:
        set stdout, file (markersFile) into MARKERS_BY_RESOLUTION 

    """
       echo $markersFile | grep -o -E '[0-9]+' | tr '\\n' '\\0' 
    """
}

// Convert the marker files to tsv

process markers_to_tsv {
    
    publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true
    
    input:
        set val(resolution), file(markersFile) from MARKERS_BY_RESOLUTION

    output:
        set val(resolution), file("markers_*.tsv") into TSV_MARKERS

    """
    cat $markersFile | sed 's/,/\t/g' > markers_${resolution}.tsv
    """
}

// Combine the listing of markers files for the manifest

process markers_lines {

    input:
        set val(resolution), file(markersFile) from TSV_MARKERS

    output:
        stdout MARKER_MANIFEST_LINES 

    """
    echo -e "cluster_markers\t${markersFile}\t$resolution"
    """
}

// Build a manifest containing bundle components

TSNE_MANIFEST_LINES
    .collectFile(name: 'tsnes.txt', newLine: false, sort: true)
    .set{ TSNE_MANIFEST_CONTENT }

MATRIX_MANIFEST_LINES
    .collectFile(name: 'matrixes.txt',  newLine: false, sort: false )
    .set{ MATRIX_MANIFEST_CONTENT }

MARKER_MANIFEST_LINES
    .collectFile(name: 'markers.txt',  newLine: false, sort: true )
    .set{ MARKER_MANIFEST_CONTENT }

// Build the manifest from collected components

process manifest {

    publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true
    
    input:
        file software from SOFTWARE
        file tsne from TSNE_MANIFEST_CONTENT
        file matrices from MATRIX_MANIFEST_CONTENT
        file clusters from BUNDLE_CLUSTERS
        file markers from MARKER_MANIFEST_CONTENT

    output:
        file "MANIFEST"

    """
        echo -e "Description\tFile\tParameterisation" > MANIFEST
        echo -e "software_versions_file\t${software}\t" >> MANIFEST
        cat ${tsne} >> MANIFEST
        cat ${matrices} >> MANIFEST
        cat ${markers} >> MANIFEST
        echo -e "cluster_memberships\t${clusters}" >> MANIFEST
    """

}
