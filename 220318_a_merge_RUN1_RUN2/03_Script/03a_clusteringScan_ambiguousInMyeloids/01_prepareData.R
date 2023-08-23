# ##########################################
# This script reads and normalize sc10x data
#
# It DOES NOT filter cells from loaded object based on #UMIs, #genes, %ribo, and
# %mito content, that should be done earlier. However, it assumes that these 
# values are contained in Seurat bject for plotting statistics.
# It also DOES NOT attribute a numeric ID to cells, this should have been done
# in a previous step (QC) before saving Seurat object as RDS file.
# ##########################################




# READ DATA
###########

## @knitr loadData

# Load Seurat from previously saved binary RDS file (must contain numID, )
sc10x = readRDS(PATH_RDS_SEURAT_OBJECT);

# TO REMOVE !!! For tests only
#sc10x = sc10x[,1:300]

cat( paste0( "\n<br>Successfuly loaded Seurat object: ", ncol( sc10x), " Cells x ", nrow( sc10x)," Genes."));

# Create a copy that will not be touched until subsetting for subclustering
sc10xBackupOriginal = sc10x;


# Read eventual TSV file to define groups of interest (and override clustering)
# in first level.
initialClustering = NULL;
if(file.exists( INITIAL_CLUSTERING_PATH) && (SUBCLUSTERING_RECURSION_CURRENTLEVEL==1))
{
  initialClustering = read.csv( INITIAL_CLUSTERING_PATH, sep = "\t", row.names = 1, stringsAsFactors = TRUE);

  if(!any( rownames( initialClustering) %in% colnames( sc10x)))
  {
    stop("Cells barcode from initial clustering TSV file do not match any barcode in Seurat object...");
  }

  if(!all( rownames( initialClustering) %in% colnames( sc10x))) 
  {
    warning( "Some cells barcode from initial clustering TSV file are not found in Seurat object...");
  }

  if(!all( colnames( sc10x) %in% rownames( initialClustering) ))
  {
    warning( "Some barcodes from Seurat object are not found in TSV file for initial clustering, they will be grouped as 'unknown' class...");
    # Restitute all barcodes of seurat object in 'initialClustering' and atribute them level 'unknown'
    initialClustering = initialClustering[colnames( sc10x), , drop = FALSE]; # Not found are NA
    rownames(initialClustering) = colnames( sc10x); # Restore correct row names for 'not found' items
    initialClustering[["identity"]] = fct_explicit_na( initialClustering[["identity"]], 
                                                       na_level = "unknown"); # Replace NA factor values by an actual level
  }

  # Sort and select cells identity according to seurat object content
  initialClustering = initialClustering[colnames( sc10x), , drop = FALSE]; # Not found are NA
}


# Read eventual TSV file to define cells coordinates
originaCellsCoordinates = NULL;
if(file.exists( ORIGINAL_CELLS_COORDINATES_PATH))
{
  originaCellsCoordinates = read.csv( ORIGINAL_CELLS_COORDINATES_PATH, sep = "\t", row.names = 1, stringsAsFactors = TRUE);

  if(!all( colnames( sc10x) %in% rownames( originaCellsCoordinates) ))
  {
    stop( "Some barcodes from Seurat object are not found in TSV file for reference coordinates, they will be removed from analysis...");
  }

  if(!all( rownames( originaCellsCoordinates) %in% colnames( sc10x))) 
  {
    warning( "Some barcodes from cells coordinates TSV file are not found in Seurat object...");
  }

  # Sort and select cells identity according to seurat object content (and ignore extra columns)
  originaCellsCoordinates = originaCellsCoordinates[colnames( sc10x), 1:2];
}

## MONITORED GENES


### Remove eventual NULL (empty) list elements from list of genes to be monitored
monitoredGroupEmpty = sapply( MONITORED_GENES, is.null);
if(any( monitoredGroupEmpty)) warning( paste("Following group(s) of genes to be monitored will be ignored because empty:", paste( names(monitoredGroupEmpty)[monitoredGroupEmpty], collapse=" - ")));
MONITORED_GENES = MONITORED_GENES[! monitoredGroupEmpty];

# Check whether genes in MONITORED_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
matchMonitoredGenes = match( toupper( unlist( MONITORED_GENES)), toupper( rownames( GetAssayData( sc10x))));
#matchMonitoredGenes = match( ( unlist( MONITORED_GENES)), ( rownames( GetAssayData( sc10x))));
monitoredGenesNotFound = unique( unlist( MONITORED_GENES)[is.na( matchMonitoredGenes)]);
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( paste0("'", monitoredGenesNotFound, "'"), collapse=" - ")));
# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
MONITORED_GENES = relist( rownames( GetAssayData( sc10x))[ matchMonitoredGenes ], skeleton = MONITORED_GENES); # Does not work with NULL list elements (removed earlier)
# Finally remove names that did not match
MONITORED_GENES = lapply( MONITORED_GENES, na.omit);


### Remove eventual NULL (empty) list elements from list of genes in modules
modulesGroupEmpty = sapply( MODULES_GENES, is.null);
if(any( modulesGroupEmpty)) warning( paste("Following module(s) of genes will be ignored because empty:", paste( names(modulesGroupEmpty)[modulesGroupEmpty], collapse=" - ")));
MODULES_GENES = MODULES_GENES[! modulesGroupEmpty];

# Check whether genes in MODULES_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
#matchModulesGenes = match( toupper( unlist( MODULES_GENES)), toupper( rownames( GetAssayData( sc10x))));
matchModulesGenes = match( ( unlist( MODULES_GENES)), ( rownames( GetAssayData( sc10x))));
modulesGenesNotFound = unique( unlist( MODULES_GENES)[is.na( matchModulesGenes)]);
if(any( is.na( matchModulesGenes))) warning( paste( "Following gene(s) from modules list could not be found in experimental data and will be ignored:", paste( paste0("'", modulesGenesNotFound, "'"), collapse=" - ")));
# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
MODULES_GENES = relist( rownames( GetAssayData( sc10x))[ matchModulesGenes ], skeleton = MODULES_GENES); # Does not work with NULL list elements (removed earlier)
# Finally remove names that did not match
MODULES_GENES = lapply( MODULES_GENES, na.omit);


### Transfer genes in very small modules (MODULES_GENES) to be analyzed individually (MONITORED_GENES)
modulesToTransfer = sapply(MODULES_GENES, length) < MONITORED_GENES_SMALL_MODULE
if(any(modulesToTransfer))
{
    warning( paste0( "Following Module(s) contained very few genes (<", MONITORED_GENES_SMALL_MODULE, "). These genes were transfered to 'Monitored genes' to be analyzed individually: ", paste( names(modulesToTransfer)[modulesToTransfer], collapse = " - ")));
    # Alter name so they can be recognized in list of Monitored genes
    names( MODULES_GENES)[modulesToTransfer] = paste0( "MOD_", names( MODULES_GENES)[modulesToTransfer]);
    MONITORED_GENES = c( MONITORED_GENES, MODULES_GENES[modulesToTransfer]);
    MODULES_GENES = MODULES_GENES[!modulesToTransfer];
}




# NORMALIZE DATA
################

## @knitr normalizeData

sc10x = NormalizeData( object = sc10x,
                       normalization.method = DATA_NORM_METHOD,
                       scale.factor = DATA_NORM_SCALEFACTOR,
                       verbose = .VERBOSE);

sc10x = ScaleData( object    = sc10x,
                   do.center = DATA_CENTER,
                   do.scale  = DATA_SCALE,
                   vars.to.regress = DATA_VARS_REGRESS,
                   verbose = .VERBOSE);
