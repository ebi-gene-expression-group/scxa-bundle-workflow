# Single-cell Expression Atlas bundling workflow

For digestion by our production processes, pipeline outputs must be organised into a specific structure we call a 'bundle'. Bundles have files organised in specific ways, and structured in manners not necessarily the same as the primary workflow outputs. This workflow carries out the appropriate restructuring and reformatting.

## Setup

### Conda/ Bioconda

Workflow dependencies are managed via Conda and Bioconda, so you'll need to set that up, see instructions [here](https://bioconda.github.io/#install-conda). 

### Nextflow

Obviously you'll need Nexflow itself. If you don't have it already you can install via Conda:

```
conda install nextflow
```

You may well want to do this within a Conda environment you create for the purpose.

## Run the workflow

### Inputs

Expected inputs are:

 * A raw count matrix, filtered by e.g. Scanpy to remove bad cells and genes.
 * A normalised count matrix, e.g. as produced by Scanpy, with the same dimensions as the filtered matrix.
 * A matrix of abundances, in transcripts per million (TPM). The rows and columns of this matrix will be filtered to match the count matrices.
 * t-SNE embeddings as output by scanpy-scripts
 * Cell clusters as output by scanpy-scripts
 * Cell cluster markers as output by scanpy-scripts
 
### Parameters
 
You can copy the default configuration, edit the Scanpy and other parameters, and provide it to Nextflow to override any of the settings. See the [Nexflow documentation](https://www.nextflow.io/docs/latest/executor.html) for executor settings.
 
### Execution

The workflow can be run directly from the repository:

```
nextflow run -config <your nextflow.config> ebi-gene-expression-group/scanpy-bundle-workflow --resultsRoot <results dir> --rawFilteredMatrix <raw matrix file name> --normalisedMatrix <normalised matrix file name> --tpmMatrix <matrix of TPM values> --clusters <clusers file> --tsneDir <subdirectory with t-SNE embedings> --markersDir <subdirectory with cluster markers> 
```

This will download the workflow, create any necessary environments, and run the workflow with the specified inputs. Input files and subdirectories are expected to exist in the specified results directory. Output will be to a directory called 'bundle' in the same location.


Future executions will use a cached copy of the pipeline, should you wish to update the code in future, you can do so like:

```
nextflow pull ebi-gene-expression-group/scanpy-bundle-workflow
```

