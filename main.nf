#!/usr/bin/env nextflow

dropletProtocols = [ '10xv1', '10xv1a', '10xv1i', '10xv2', 'drop-seq' ]
smartProtocols = [ 'smart-seq', 'smart-seq2', 'smarter', 'smart-like' ]

resultsRoot = params.resultsRoot
masterWorkflow = params.masterWorkflow

if ( params.containsKey('tertiaryWorkflow' )){
    tertiaryWorkflow = params.tertiaryWorkflow
}else{
    tertiaryWorkflow = 'none'
}

PROTOCOL_LIST = Channel.from(params.protocolList.split(','))
RAW_MATRIX = Channel.fromPath( "$resultsRoot/${params.rawMatrix}", checkIfExists: true)
REFERENCE_FASTA = Channel.fromPath( "${params.referenceFasta}", checkIfExists: true )
REFERENCE_GTF = Channel.fromPath( "${params.referenceGtf}", checkIfExists: true )

if ( tertiaryWorkflow == 'scanpy-workflow'){
    RAW_FILTERED_MATRIX = Channel.fromPath( "$resultsRoot/${params.rawFilteredMatrix}", checkIfExists: true)
    NORMALISED_MATRIX = Channel.fromPath( "$resultsRoot/${params.normalisedMatrix}", checkIfExists: true)
    SCANPY_CLUSTERS = Channel.fromPath( "$resultsRoot/${params.clusters}", checkIfExists: true)
    SCANPY_TSNE = Channel.fromPath( "$resultsRoot/${params.tsneDir}/embeddings*.csv", checkIfExists: true )
    SCANPY_MARKERS = Channel.fromPath( "$resultsRoot/${params.markersDir}/markers_*.csv", checkIfExists: true )
}else{
    RAW_FILTERED_MATRIX = Channel.empty()
    NORMALISED_MATRIX = Channel.empty()
    RAW_TPM_MATRIX = Channel.empty()
    SCANPY_CLUSTERS = Channel.empty()
    SCANPY_TSNE = Channel.empty()
    SCANPY_MARKERS = Channel.empty()
}

if ( params.containsKey('tpmMatrix') ){
    RAW_TPM_MATRIX = Channel.fromPath( "$resultsRoot/${params.tpmMatrix}", checkIfExists: true)
}else{
    RAW_TPM_MATRIX = Channel.empty()
}

// Send channels to different processes

RAW_MATRIX.into{
    RAW_MATRIX_FOR_MTX
    RAW_MATRIX_FOR_TSV
}

RAW_FILTERED_MATRIX.into{
    RAW_FILTERED_MATRIX_FOR_TPM_FILTERING
    RAW_FILTERED_MATRIX_FOR_MTX
    RAW_FILTERED_MATRIX_FOR_TSV
}

NORMALISED_MATRIX.into{
    NORMALISED_MATRIX_FOR_MTX
    NORMALISED_MATRIX_FOR_TSV
}

RAW_TPM_MATRIX.into{
    RAW_TPM_MATRIX_FOR_FILTERING
    RAW_TPM_MATRIX_FOR_MTX
    RAW_TPM_MATRIX_FOR_TSV
}

// Make manifest lines for references

process reference_lines {

    input:
        file(referenceFasta) from REFERENCE_FASTA
        file(referenceGtf) from REFERENCE_GTF

    output:
        stdout REFERENCE_MANIFEST_LINES 

    """
    echo -e "reference_transcriptome\t$referenceFasta\t"
    echo -e "reference_annotation\t$referenceGtf\t"
    """
}

// We need a TPM matrix for the same cell and gene sets as the Scanpy-filtered
// count matrices. Filter the abundance matrix to match the count matrices

process filter_tpms {

    conda "${workflow.projectDir}/envs/bioconductor-dropletutils.yml"
    
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 10
    
    input:
        file(filteredCountsMatrix) from RAW_FILTERED_MATRIX_FOR_TPM_FILTERING
        file(tpmMatrix) from RAW_TPM_MATRIX_FOR_FILTERING

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
        set val(matName), file(mtx), file(genes), file(barcodes) from TPM_FILTER_CELLS_MTX.join(TPM_FILTER_CELLS_MTX_ROWS).join(TPM_FILTER_CELLS_MTX_COLS)

    output:
        file("${matName}_filter_cells_genes.zip") into RAW_FILTERED_TPM_MATRIX

    """
        mkdir -p ${matName}_filter_cells_genes
        mv matrix.mtx genes.tsv barcodes.tsv ${matName}_filter_cells_genes 
        zip -r ${matName}_filter_cells_genes.zip ${matName}_filter_cells_genes
    """
}

// Report sotware from master workflow

process master_workflow_software {

    output:
        file('master.software.tsv') into MASTER_SOFTWARE

    """
        generateSoftwareReport.sh ${masterWorkflow} master.software.tsv
    """        
}

// Report software versions

process make_base_software_report {

    input:
        val(protocol) from PROTOCOL_LIST 

    output:
        file "${protocol}.software.tsv" into BASE_SOFTWARE

    script:

        if ( dropletProtocols.contains(protocol) ) {
            baseWorkflow = 'scxa-droplet-quantification-workflow'
        } else if ( smartProtocols.contains(protocol) ) {
            baseWorkflow = 'scxa-smartseq-quantification-workflow'
        } else {
            baseWorkflow = ''
        }
    

    """
        if [ -z "$baseWorkflow" ]; then
            echo "I don't know what workflow would have been used for $protocol" 1>&2
            exit 1
        fi

        generateSoftwareReport.sh ${baseWorkflow} ${protocol}.software.tsv
    """
}

MASTER_SOFTWARE
    .concat(BASE_SOFTWARE)
    .collectFile(name: 'software.tsv', newLine: true, keepHeader: true )
    .set { ALL_BASE_SOFTWARE }

if ( tertiaryWorkflow == 'scanpy-workflow'){

    process make_tertiary_software_report {

        publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true
        
        output:
            file "${tertiaryWorfklow}.software.tsv" into TERTIARY_SOFTWARE

        """
            generateSoftwareReport.sh ${tertiaryWorkflow} software.tsv
        """
    }

    ALL_BASE_SOFTWARE
        .concat(TERTIARY_SOFTWARE)
        .set{ SOFTWARE }

}else{
    ALL_BASE_SOFTWARE
        .set{ SOFTWARE }
}

SOFTWARE
    .collectFile(name: 'software.tsv', storeDir: "$resultsRoot/bundle", keepHeader: true)
    .set{ MERGED_SOFTWARE }


// Make sure we have no duplicate lines from merging software reports

process finalise_software {

    publishDir "$resultsRoot/bundle", mode: 'copy', overwrite: true

    input:
        file(software) from MERGED_SOFTWARE

    output:
        file('output/software.tsv') into SOFTWARE_FOR_MANIFEST

    """
    mkdir -p output
    head -n 1 $software > output/software.tsv
    tail -n +2 $software | sort | uniq >> output/software.tsv
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


Channel.from( 'raw', 'tpm', 'raw_filtered', 'normalised', 'tpm_filtered' ).into{
    EXPRESSION_TYPES_FOR_MTX
    EXPRESSION_TYPES_FOR_TSV
}

RAW_FILTERED_TPM_MATRIX.into{
    RAW_FILTERED_TPM_MATRIX_FOR_MTX
    RAW_FILTERED_TPM_MATRIX_FOR_TSV
}

// Collect a list of matrices to repackage 

RAW_MATRIX_FOR_MTX
    .concat(RAW_TPM_MATRIX_FOR_MTX)
    .concat(RAW_FILTERED_MATRIX_FOR_MTX)
    .concat(NORMALISED_MATRIX_FOR_MTX)
    .concat(RAW_FILTERED_TPM_MATRIX_FOR_MTX)
    .merge(EXPRESSION_TYPES_FOR_MTX)
    .set{
        MATRICES_TO_REPACKAGE
    }

process repackage_matrices {

    publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true
    
    input:
        set file(expressionMatrix), val(expressionType) from MATRICES_TO_REPACKAGE

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

// Collect a list of matrices to convert to tsv

RAW_MATRIX_FOR_TSV
    .concat(RAW_TPM_MATRIX_FOR_TSV)
    .concat(RAW_FILTERED_MATRIX_FOR_TSV)
    .concat(NORMALISED_MATRIX_FOR_TSV)
    .concat(RAW_FILTERED_TPM_MATRIX_FOR_TSV)
    .set{
        MATRICES_FOR_TSV
    }

// Count the number of cells

process cell_count {

    input:
        file(expressionMatrix) from MATRICES_FOR_TSV

    output:
        set stdout, file("out/${expresssionMatrix}") into MATRICES_FOR_TSV_WITH_COUNT

    """
        zipdir=\$(unzip -qql ${expressionMatrix.getBaseName()}.zip | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/||')
        unzip ${expressionMatrix.getBaseName()}        
    
        zcat \${zipdir}/barcodes.tsv | wc

        mkdir -p out
        cp -p $expressionMatrix out/${expressionMatrix}
    """
}

// Only make TSVs for small enough matrices

SMALL_MATRICES = Channel.create()
BIG_MATRICES = Channel.create()

MATRICES_FOR_TSV_WITH_COUNT
    .merge( EXPRESSION_TYPES_FOR_TSV)
    .choice( SMALL_MATRICES, BIG_MATRICES ) {a -> 
    a[0] < params.largeMatrixThreshold ? 0 : 1
}

// Make tsv-format matrices

process mtx_to_tsv {
    
    conda "${workflow.projectDir}/envs/bioconductor-dropletutils.yml"

    publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true
    
    memory { 5.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(cellCount), file(expressionMatrix), val(expressionType) from SMALL_MATRICES
        
    output:
        set val(expressionType), file("${expressionType}/${expressionType}.tsv") into TSV_MATRICES

    """
        #!/usr/bin/env Rscript
        
        suppressPackageStartupMessages(require(DropletUtils))
        suppressPackageStartupMessages(require(data.table))
        source(file.path(Sys.getenv(c("SCXA_BIN")), "utils.R"))
   
        unzip('$expressionMatrix')
        sce <- read10xCounts(sub('.zip', '', '$expressionMatrix'))
        colnames(sce) <- colData(sce)\$Barcode 

        dir.create('${expressionType}')
        write.tsv(as.data.frame(cbind(Feature = rownames(sce), as.matrix(assays(sce)[[1]])), col.names = c('Feature', colData(sce)\$Barcode)), "${expressionType}/${expressionType}.tsv")      
    """
}

// Mark the big matrices as not having TSVs, but include them in the data
// structure for the matrix_lines process

BIG_MATRICES
    .map{ row-> tuple( row[1], file('NOTSV')) }        
    .concat( TSV_MATRICES)
    .set { TSV_AND_NOTSV_MATRICES }

// Make manifest lines for matrices

process matrix_lines {

    input:
        set val(expressionType), file(matrixRows), file(matrixCols), file(matrixContent), file(tsvMatrix) from MTX_MATRIX_ROWNAMES.join(MTX_MATRIX_COLNAMES).join(MTX_MATRIX_CONTENT).join(TSV_AND_NOTSV_MATRICES)

    output:
        stdout MATRIX_MANIFEST_LINES 

    """
    echo -e "mtx_matrix_rows\t$expressionType/$matrixRows\t$expressionType"
    echo -e "mtx_matrix_cols\t$expressionType/$matrixCols\t$expressionType"
    echo -e "mtx_matrix_content\t$expressionType/$matrixContent\t$expressionType"
    if [ "${tsvMatrix.getName}" != 'NOTSV' ]; then 
        echo -e "tsv_matrix\t$expressionType/$tsvMatrix\t$expressionType"
    fi
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

process base_manifest {

    input:
        file matrices from MATRIX_MANIFEST_CONTENT
        file software from SOFTWARE_FOR_MANIFEST 
        file reference from REFERENCE_MANIFEST_LINES

    output:
        file "BASE_MANIFEST" into BASE_MANIFEST

    """
        echo -e "Description\tFile\tParameterisation" > BASE_MANIFEST
        echo -e "software_versions_file\t\$(basename ${software})\t" >> BASE_MANIFEST
        cat ${matrices} >> BASE_MANIFEST
        cat ${reference} >> BASE_MANIFEST
        echo -e protocol\t\t${params.protocolList} >> BASE_MANIFEST
    """

}

// Add in any tertiary data to the bundle. If there's no teriary data, just
// copy the base manifest

if ( tertiaryWorkflow == 'scanpy-workflow'){

    process tertiary_manifest {

        publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true
        
        input:
            file baseManifest from BASE_MANIFEST
            file tsne from TSNE_MANIFEST_CONTENT
            file clusters from BUNDLE_CLUSTERS
            file markers from MARKER_MANIFEST_CONTENT

        output:
            file "MANIFEST"

        """
            cp $baseManifest MANIFEST
            cat ${tsne} >> MANIFEST
            cat ${markers} >> MANIFEST
            echo -e "cluster_memberships\t${clusters}" >> MANIFEST
        """

    }

}else{

    process publish_manifest {
        
        publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true
        
        input:
            file baseManifest from BASE_MANIFEST

        output:
            file "MANIFEST"

        """
            cp $baseManifest MANIFEST
        """
    } 

}

