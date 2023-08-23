# ##########################################
# This script reads a precomputed Seurat object, and eventual external data from
# TSV files for clustering results and 2D representation of cells such as result
# of dimensionality reduction algorithms. 
#
# ##########################################




## @knitr loadData




#### Seurat object

# Load Seurat from previously saved binary RDS file (must contain numID, )
sc10x = readRDS(PATH_RDS_SEURAT_OBJECT);

cat( paste0( "\n<br>Successfuly loaded Seurat object: ", ncol( sc10x), " Cells x ", nrow( sc10x)," Genes."));




#### External data files

# Read eventual TSV file to define groups of interest (and override eventual 
# clustering contained in seurat object).
if(file.exists( EXTERNAL_CLUSTERING_PATH))
{
  externalClusters = read.csv( EXTERNAL_CLUSTERING_PATH, sep = "\t", row.names = 1);

  if(!any( rownames( externalClusters) %in% colnames( sc10x)))
  {
    stop("Cells barcode from initial clustering TSV file do not match any barcode in Seurat object...");
  }

  if(!all( rownames( externalClusters) %in% colnames( sc10x))) 
  {
    warning( "Some cells barcode from initial clustering TSV file are not found in Seurat object...");
  }

  if(!all( colnames( sc10x) %in% rownames( externalClusters) ))
  {
    warning( "Some barcodes from Seurat object are not found in TSV file for initial clustering, they will be grouped as 'unknown' class...");
    # Restitute all barcodes of seurat object in 'externalClusters' and atribute them level 'unknown'
    externalClusters = externalClusters[colnames( sc10x),]; # Not found are NA
    rownames(externalClusters) = colnames( sc10x);
    externalClusters[["dentity"]] = fct_explicit_na( externalClusters[["identity"]], 
                                                      na_level = "unknown"); # Replace NA factor values by an actual level
  }

  # Sort and select cells identiity according to seurat object content
  externalClusters = externalClusters[colnames( sc10x),];

  # Set cells identity directly into Seurat object 
  Idents( sc10x) = externalClusters[["identity"]];
}




# Read eventual TSV file defining a color for each cluster (attribute ggplot defaults otherwise)
clustersColor = NULL;
if(file.exists( EXTERNAL_CLUSTERSCOLOR_PATH))
{
  clustersColor =  as.matrix( read.table( EXTERNAL_CLUSTERSCOLOR_PATH, 
                                          sep = "\t", 
                                          header = FALSE, 
                                          comment.char= "", 
                                          row.names = 1, 
                                          stringsAsFactors=FALSE))[,1]; # Convert to matrix to get a named vector when extracting column
  

  # Check consistency between colors file and clusters
  clusterInColorFile = levels( Idents( sc10x)) %in% names( clustersColor);
  if(!all( clusterInColorFile)) stop( paste0( "Following cluster name(s) could not be found in file defining clusters color: ", paste( levels( Idents( sc10x))[!clusterInColorFile], collapse = " - "), "."));

  clusterNameExists = names( clustersColor) %in% levels( Idents( sc10x));
  if(!all( clusterNameExists)) warning( paste0( "Following cluster name from file defining colors could not be found in loaded data: ", paste( names( clustersColor)[!clusterNameExists], collapse = " - "), "."));

} else
{
  # Define a set of colors for clusters (based on ggplot default)
  clustersColor = hue_pal()( nlevels( Idents( sc10x)));
  names( clustersColor) = levels( Idents( sc10x));
}




# Read eventual TSV file to define groups of interest (and override clustering)
# in first level.
cellsCoordinates = NULL;
if(file.exists( CELLS_COORDINATES))
{
  cellsCoordinates = read.csv( CELLS_COORDINATES, sep = "\t", row.names = 1);

  if(!all( colnames( sc10x) %in% rownames( cellsCoordinates) ))
  {
    stop( "Some barcodes from Seurat object are not found in TSV file for reference coordinates...");
  }

  if(!all( rownames( cellsCoordinates) %in% colnames( sc10x))) 
  {
    warning( "Some cells barcode from reference coordinates TSV file are not found in Seurat object...");
  }

  # Sort and select cells identiity according to seurat object content
  cellsCoordinates = cellsCoordinates[colnames( sc10x),];
  
} else if( CELLS_COORDINATES %in% Reductions( sc10x) ) 
{
  cellsCoordinates = Embeddings(Reductions(sc10x, slot = CELLS_COORDINATES))

  # PCA can have arbitrary number of dimensions, use requested ones (2d plot)
  if(CELLS_COORDINATES == "pca")
  {
    cellsCoordinates = cellsCoordinates[, PCA_DIMS[1:2]];
  }

} else stop("Parameter 'CELLS_COORDINATES' must refer to a DimReduc stored in loaded Seurat object, or contain a path to a valid file containing 2D cells coordinates...");



