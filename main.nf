#!/usr/bin/env nextflow

dropletProtocols = [ '10xv1', '10xv1a', '10xv1i', '10xv2', '10xv3', 'drop-seq' ]
smartProtocols = [ 'smart-seq', 'smart-seq2', 'smarter', 'smart-like' ]
expressionTypes = [ 'raw' ]

resultsRoot = params.resultsRoot
masterWorkflow = params.masterWorkflow
tertiarySoftwareReport = 'None'

if ( params.containsKey('tertiaryWorkflow' )){
    tertiaryWorkflow = params.tertiaryWorkflow
    if ( tertiaryWorkflow == 'scanpy-galaxy' ){
        tertiarySoftwareReport = "$resultsRoot/${params.tertiarySoftwareReport}"
    }      
}else{
    tertiaryWorkflow = 'none'
}

// Test protocols

protocolList = params.protocolList.split(',')
isDroplet = protocolList.any { dropletProtocols.contains( it ) }
isSmart = protocolList.any { smartProtocols.contains( it ) }

if (isDroplet && isSmart){
    println("Error: Protocol list ${params.protocolList} for this experiment contains both droplet and smart protocols.")
    System.exit(1)
}else if(isDroplet){
    println("Droplet experiment from ${params.protocolList}")
}else if (isSmart){
    println("Smart experiment from ${params.protocolList}")
}else{
    println("Cannot determine experiment type for ${params.protocolList}")
    System.exit(1)
}

// See what other inputs are provided

RAW_MATRIX = Channel.fromPath( "$resultsRoot/${params.rawMatrix}", checkIfExists: true)
REFERENCE_FASTA = Channel.fromPath( "$resultsRoot/${params.referenceFasta}", checkIfExists: true ).first()
REFERENCE_GTF = Channel.fromPath( "$resultsRoot/${params.referenceGtf}", checkIfExists: true ).first()

if ( tertiaryWorkflow == 'scanpy-workflow' || tertiaryWorkflow == 'scanpy-galaxy' ){
    expressionTypes = expressionTypes + [ 'raw_filtered', 'filtered_normalised' ]

    RAW_FILTERED_MATRIX = Channel.fromPath( "$resultsRoot/${params.rawFilteredMatrix}", checkIfExists: true)
    NORMALISED_MATRIX = Channel.fromPath( "$resultsRoot/${params.normalisedMatrix}", checkIfExists: true)
    SCANPY_CLUSTERS = Channel.fromPath( "$resultsRoot/${params.clusters}", checkIfExists: true)
    SCANPY_TSNE = Channel.fromPath( "$resultsRoot/${params.tsneDir}/tsne_perplexity*.csv", checkIfExists: true )
    SCANPY_CLUSTER_MARKERS = Channel.fromPath( "$resultsRoot/${params.markersDir}/markers_*.csv" )
    SCANPY_META_MARKERS = Channel.fromPath( "$resultsRoot/${params.markersDir}/*_markers.csv" )
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
    expressionTypes = expressionTypes + [ 'tpm', 'tpm_filtered' ]
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

process reference_manifest_lines {

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

// Publish reference files to bundle

process publish_reference {
    
    publishDir "$resultsRoot/bundle", mode: 'copy', overwrite: true
    
    input:
        file(referenceFasta) from REFERENCE_FASTA
        file(referenceGtf) from REFERENCE_GTF

    output:
        file("reference/$referenceFasta")
        file("reference/$referenceGtf")
        
    """
    mkdir -p reference
    cp -P $referenceFasta reference
    cp -P $referenceGtf reference
    """
}

// Reference lines for supplementary/ software table. Putting these data there
// is a bit of a hack until we figure out better ways of flagging the reference
// used

process reference_supplementary_lines {

    input:
        file(referenceFasta) from REFERENCE_FASTA
        file(referenceGtf) from REFERENCE_GTF

    output:
        file("software_reference.tsv") into REFERENCE_SOFTWARE 

    """
    ensembl_version=\$(echo $referenceGtf | cut -d '.' -f 3)
    echo "Analysis\tSoftware\tVersion\tCitation" > software_reference.tsv
    echo -e "Reference\tEnsembl\t\$ensembl_version\t$referenceFasta, $referenceGtf" >> software_reference.tsv
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

    output:
        file "base.software.tsv" into BASE_SOFTWARE

    script:

        if ( isDroplet ) {
            baseWorkflow = 'scxa-workflows/w_droplet_quantification'
        } else if ( isSmart ) {
            baseWorkflow = 'scxa-workflows/w_smart-seq_quantification'
        } else {
            baseWorkflow = ''
        }
    

    """
        generateSoftwareReport.sh ${baseWorkflow} base.software.tsv
    """
}

MASTER_SOFTWARE
    .concat(BASE_SOFTWARE)
    .concat(REFERENCE_SOFTWARE)
    .collectFile(name: 'software.tsv', newLine: true, keepHeader: true )
    .set { ALL_BASE_SOFTWARE }

if ( tertiaryWorkflow == 'scanpy-workflow' || tertiaryWorkflow == 'scanpy-galaxy'){

    process make_tertiary_software_report {

        output:
            file "${tertiaryWorkflow}.software.tsv" into TERTIARY_SOFTWARE

        script:

        if ( tertiaryWorkflow == 'scanpy-workflow' )

            """
                generateSoftwareReport.sh ${tertiaryWorkflow} ${tertiaryWorkflow}.software.tsv
            """

        else
            """
               cp ${tertiarySoftwareReport} ${tertiaryWorkflow}.software.tsv
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
    .collectFile(name: 'collected_software.tsv', keepHeader: true, newLine: false)
    .set{ MERGED_SOFTWARE }


// Make sure we have no duplicate lines from merging software reports

process finalise_software {

    publishDir "$resultsRoot/bundle", mode: 'copy', overwrite: true

    input:
        file(software) from MERGED_SOFTWARE

    output:
        file('software.tsv') into SOFTWARE_FOR_MANIFEST

    """
    head -n 1 $software > software.tsv
    tail -n +2 $software | sort | uniq >> software.tsv

    # Fetch the current SHA of the config repo, for reproducibility

    pushd $SCXA_WORKFLOW_ROOT/workflow/scxa-workflows > /dev/null
    current_sha=\$(git rev-parse --short HEAD)
    popd > /dev/null
    echo -e "Configuration\tscxa-workflows\t\$current_sha\thttps://github.com/ebi-gene-expression-group/scxa-workflows" >> software.tsv

    # Remove any empty lines
    sed -i '/^\$/d' software.tsv
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

// Repackage the matrices 

Channel.from( expressionTypes ).into{
    EXPRESSION_TYPES_FOR_MTX
    EXPRESSION_TYPES_FOR_TSV
}

RAW_FILTERED_TPM_MATRIX.into{
    RAW_FILTERED_TPM_MATRIX_FOR_MTX
    RAW_FILTERED_TPM_MATRIX_FOR_TSV
}

// Collect a list of matrices to repackage 

RAW_MATRIX_FOR_MTX
    .concat(RAW_FILTERED_MATRIX_FOR_MTX)
    .concat(NORMALISED_MATRIX_FOR_MTX)
    .concat(RAW_TPM_MATRIX_FOR_MTX)
    .concat(RAW_FILTERED_TPM_MATRIX_FOR_MTX)
    .merge(EXPRESSION_TYPES_FOR_MTX)
    .set{
        MATRICES_TO_REPACKAGE
    }

process repackage_matrices {

    cache 'deep'

    publishDir "$resultsRoot/bundle", mode: 'copy', overwrite: true
    
    input:
        set file(expressionMatrix), val(expressionType) from MATRICES_TO_REPACKAGE

    output:
        set val(expressionType), file("${expressionType}/genes.tsv.gz") into MTX_MATRIX_ROWNAMES
        set val(expressionType), file("${expressionType}/barcodes.tsv.gz") into MTX_MATRIX_COLNAMES
        set val(expressionType), file("${expressionType}/matrix.mtx.gz") into MTX_MATRIX_CONTENT


    """
        zipdir=\$(unzip -qql ${expressionMatrix.getBaseName()}.zip | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/||')
        unzip ${expressionMatrix.getBaseName()}        
   
        if [ "\$zipdir" != ${expressionType} ]; then
            ln -s \$zipdir ${expressionType}            
        fi
 
        gzip ${expressionType}/matrix.mtx
        gzip ${expressionType}/genes.tsv
        gzip ${expressionType}/barcodes.tsv
    """        

}

// Make cell - library mappings for droplet experiments 

MTX_MATRIX_COLNAMES
    .into{
        MTX_MATRIX_COLNAMES_FOR_MANIFEST_LINES
        MTX_MATRIX_COLNAMES_FOR_CELLMAPPING
    }

process cell_library_mappings {

    publishDir "$resultsRoot/bundle/$expressionType", mode: 'move'

    input:
        set val(expressionType), file(barcodesFile) from MTX_MATRIX_COLNAMES_FOR_CELLMAPPING

    output:
        set val(expressionType), file('cell_to_library.txt') optional true

    script:

        def sampleField = params.fields.run
        if ( params.fields.containsKey('techrep') ){
            sampleField = params.fields.techrep
        }

    """
        if [ "$isDroplet" = 'true' ]; then
            echo "# $sampleField" > cell_to_library.txt

            zcat $barcodesFile | while read -r b; do 
                barcode=\${b##*-} 
                run=\${b/-\$barcode/''}
                echo -e "\$b\t\$run" 
            done >> cell_to_library.txt
        fi
    """
}

// Collect a list of matrices to convert to tsv

RAW_MATRIX_FOR_TSV
    .concat(RAW_FILTERED_MATRIX_FOR_TSV)
    .concat(NORMALISED_MATRIX_FOR_TSV)
    .concat(RAW_TPM_MATRIX_FOR_TSV)
    .concat(RAW_FILTERED_TPM_MATRIX_FOR_TSV)
    .merge( EXPRESSION_TYPES_FOR_TSV)
    .set{
        MATRICES_FOR_TSV
    }

// Count the number of cells

process cell_count {

    input:
        set file(expressionMatrix), val(expressionType) from MATRICES_FOR_TSV

    output:
        set stdout, file("out/${expressionMatrix}"), val(expressionType) into MATRICES_FOR_TSV_WITH_COUNT

    """
        zipdir=\$(unzip -qql ${expressionMatrix} | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/||')
        unzip -p ${expressionMatrix} \${zipdir}/barcodes.tsv | wc -l | tr -d \'\\n\'  
        
        mkdir -p out
        cp -p $expressionMatrix out/${expressionMatrix}
    """
}

// Only make TSVs for small enough matrices

SMALL_MATRICES = Channel.create()
BIG_MATRICES = Channel.create()
 
MATRICES_FOR_TSV_WITH_COUNT
    .choice( SMALL_MATRICES, BIG_MATRICES ) {a -> a[0].toInteger() < params.largeMatrixThreshold ? 0 : 1 }

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
    .map{ row-> tuple( row[2], file('NOTSV')) }        
    .concat( TSV_MATRICES)
    .set { 
        TSV_AND_NOTSV_MATRICES 
    }

// Make manifest lines for matrices

process matrix_lines {

    input:
        set val(expressionType), file(matrixRows), file(matrixCols), file(matrixContent), file(tsvMatrix) from MTX_MATRIX_ROWNAMES.join(MTX_MATRIX_COLNAMES_FOR_MANIFEST_LINES).join(MTX_MATRIX_CONTENT).join(TSV_AND_NOTSV_MATRICES)

    output:
        stdout MATRIX_MANIFEST_LINES 

    """
    echo -e "mtx_matrix_rows\t$expressionType/$matrixRows\t$expressionType"
    echo -e "mtx_matrix_cols\t$expressionType/$matrixCols\t$expressionType"
    echo -e "mtx_matrix_content\t$expressionType/$matrixContent\t$expressionType"
    if [ "${tsvMatrix.name}" != 'NOTSV' ]; then 
        echo -e "tsv_matrix\t$expressionType/$tsvMatrix\t$expressionType"
    fi
    """
}

// Renumber clusters where numbering starts at 0

process renumber_clusters {
    
    memory { 5.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20

    input:
        file 'possibly_misnumbered_clusters.txt' from SCANPY_CLUSTERS
    
    output:
        file 'clusters_for_bundle.txt' into FINAL_CLUSTERS

    """
        #!/usr/bin/env Rscript

        clusters <- read.delim('possibly_misnumbered_clusters.txt', check.names=FALSE)

        if (min(clusters[,c(-1,-2)]) == 0){
            clusters[,c(-1,-2)] <- clusters[,c(-1,-2)]+1
        }

        write.table(clusters, file='clusters_for_bundle.txt', sep="\t", quote=FALSE, row.names = FALSE)
    """
}

// Find out what resolutions are represented by the marker files

process mark_marker_resolutions {

    input:
        file markersFile from SCANPY_CLUSTER_MARKERS

    output:
        stdout, file (markersFile) into CLUSTER_MARKERS_BY_RESOLUTION 

    """
        echo $markersFile | grep -o -E '[0-9]+' | tr '\\n' '\\0' 
    """
}

// Find out what resolutions are represented by the marker files

process mark_marker_meta {

    input:
        file markersFile from SCANPY_META_MARKERS

    output:
        set val('meta_markers'), stdout, file (markersFile) into META_MARKERS_BY_VAR 

    """
        echo "$markersFile" | rev | cut -d"_" -f2-  | rev
    """
}

// Convert the marker files to tsv

process renumber_cluster_markers_to_tsv {
    
    publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true

    memory { 5.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(resolution), file(markersFile) from CLUSTER_MARKERS_BY_RESOLUTION

    output:
        set set val('cluster_markers'), val(resolution), file("out/${markersFile.baseName}") into RENUMBERED_CLUSTER_MARKERS_BY_RESOLUTION

    """
    #!/usr/bin/env Rscript

    markers <- read.csv('${markersFile}', check.names = FALSE)

    if ('groups' %in% names(markers) && min(markers\$groups) == 0){
        markers\$groups <- markers\$groups + 1
    }else if ('cluster' %in% names(markers) && min(markers\$cluster) == 0){
        markers\$cluster <- markers\$cluster + 1
    }
    dir.create('out')
    write.csv(markers, file='out/${markersFile.baseName}.csv', sep="\t", quote=FALSE, row.names=FALSE)
    """
}

process markers_to_tsv {
    
    input:
        set val(markerType), val(markerVal), file(markersFile) from RENUMBERED_CLUSTER_MARKERS_BY_RESOLUTION.concat(META_MARKERS_BY_VAR)

    output:
        set val(markerType), val(markerVal), file("${markersFile.baseName}.tsv") into TSV_MARKERS
        

    """
    #!/usr/bin/env Rscript

    markers <- read.csv('${markersFile}', check.names = FALSE)
    write.table(markers, file='${markersFile.baseName}.tsv', sep="\t", quote=FALSE, row.names=FALSE)
    """
}

// Combine the listing of markers files for the manifest

process markers_lines {

    input:
        set val(markerType), val(markerVal), file(markersFile) from TSV_MARKERS

    output:
        stdout MARKER_MANIFEST_LINES 

    """
    echo -e "$markerType\t${markersFile}\t$markerVal"
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

if ( tertiaryWorkflow == 'scanpy-workflow' || tertiaryWorkflow == 'scanpy-galaxy'){


    BASE_MANIFEST
        .concat(TSNE_MANIFEST_CONTENT)
        .concat(MARKER_MANIFEST_CONTENT)
        .collectFile(name: 'manifest_lines.tsv', newLine: false, sort: 'index' )
        .set { STARTING_MANIFEST }

    process tertiary_manifest {

        publishDir "$resultsRoot/bundle", mode: 'move', overwrite: true
        
        input:
            file startingManifest from STARTING_MANIFEST
            file clusters from FINAL_CLUSTERS

        output:
            file "MANIFEST"
            file(clusters)

        """
            cp $startingManifest MANIFEST
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

