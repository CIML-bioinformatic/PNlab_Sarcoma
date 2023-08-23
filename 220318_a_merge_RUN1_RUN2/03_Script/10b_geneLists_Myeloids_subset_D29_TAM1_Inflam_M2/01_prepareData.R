# ##########################################
# This script reads a precomputed Seurat object, and eventual external data from
# TSV files for clustering results and 2D representation of cells such as result
# of dimensionality reduction algorithms. 
#
# ##########################################




## @knitr loadData


#### Seurat object

# Load full session object from previous step.
# Only on first run as this report is designed to be started twice, see Rmd...
previousResultsEnv = new.env() # Load in separate environment to prevent overwriting critical variables specific to current analysis step (e.g. output path defined in analysisParams)
load(PATH_RDATA_PREVIOUS_RESULTS, envir = previousResultsEnv)

# Extract variable of interest to current environment for easy access
sc10x = previousResultsEnv[['sc10x']]
#geneExpByHTO = previousResultsEnv[['geneExpByHTO']]
#geneExpByClusterAndHTO = previousResultsEnv[['geneExpByClusterAndHTO']]
clustersColor = previousResultsEnv[['clustersColor']]
HTOsColor = previousResultsEnv[['HTOsColor']]

# Remove loaded environment to prevent inflating current analysis environment 
# with variables of previous analysis step when saving in final step
rm(previousResultsEnv)

cat( paste0( "\n<br>Successfuly loaded Seurat object: ", ncol( sc10x), " Cells x ", nrow( sc10x)," Genes."));




# CHECK GENES LISTS
###################

## @knitr loadGenesLists

### Remove eventual NULL (empty) list elements from list of genes to be monitored
monitoredGroupEmpty = sapply( GENES_LIST, is.null);
if(any( monitoredGroupEmpty)) warning( paste("Following group(s) of genes to be monitored will be ignored because empty:", paste( names(monitoredGroupEmpty)[monitoredGroupEmpty], collapse=" - ")));
GENES_LIST = GENES_LIST[! monitoredGroupEmpty];

# Check whether genes in GENES_LIST are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
#matchMonitoredGenes = match( toupper( unlist( GENES_LIST)), toupper( rownames( GetAssayData( sc10x))));
matchMonitoredGenes = match( ( unlist( GENES_LIST)), ( rownames( GetAssayData( sc10x))));
monitoredGenesNotFound = unique( unlist( GENES_LIST)[is.na( matchMonitoredGenes)]);
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( paste0("'", monitoredGenesNotFound, "'"), collapse=" - ")));
# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
GENES_LIST = relist( rownames( GetAssayData( sc10x))[ matchMonitoredGenes ], skeleton = GENES_LIST); # Does not work with NULL list elements (removed earlier)
# Finally remove names that did not match
GENES_LIST = lapply( GENES_LIST, na.omit);




# FILTERINGS
############

## @knitr filtering

# Addition in this script version only to subset populations for heatmaps (this version only)
if( exists( "HTO_SELECT") && ( ! is.null( HTO_SELECT)))
{
  warning( paste( "Filtering cells population of interest based on HTOs: "), paste(HTO_SELECT, collapse = " - "))
  selection = NULL;
  selection = Reduce(function(x,y){grepl(y, sc10x[["factorHTO", drop = TRUE]], perl = TRUE) | x}, HTO_SELECT, init = logical(ncol(sc10x)))
  if(! is.null( selection))
  {
    warning( paste( "Number of cells selected from current object:", sum(selection)))
  }
  if(sum(selection) == 0) stop("No cells remaining after selection")
  sc10x = sc10x[, selection]

  # Remove eventual filtered levels from factor
  sc10x[["factorHTO"]] = factor(sc10x[["factorHTO", drop = TRUE]])
}

if( exists( "CLUSTER_SELECT") && ( ! is.null( CLUSTER_SELECT)))
{
  warning( paste( "Filtering cells population of interest based on Clusters identity: "), paste(CLUSTER_SELECT, collapse = " - "))
  selection = NULL;
  selection = Reduce(function(x,y){grepl(y, Idents(sc10x), perl = TRUE) | x}, CLUSTER_SELECT, init = logical(ncol(sc10x)))
  if(! is.null( selection))
  {
    warning( paste( "Number of cells selected from current object:", sum(selection)))
  }
  if(sum(selection) == 0) stop("No cells remaining after selection")
  sc10x = sc10x[, selection]

}
