###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#
ANALYSIS_STEP_NAME = "13c_genesListUMAP_APC"
PATH_ANALYSIS_OUTPUT = file.path( PATH_PROJECT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "APC - Genes List UMAP"


# !!! THIS ANALYSIS STEP NEEDS TO BE STARTED TWICE !!!
# On first time, it does all computation and generates figures as external files
# (using png/pdf and dev.off) which causes the rmarkdown layout to fail (figures
# in wrong tabset).
# On second execution, the RDATA result of previous run is loaded and previously
# generated external figures are not rendered again, but just integrated to
# rmarkdown using computed file path. It also skips time consuming executions
# (compute DEG analyses and enrichments).
# See chunk 'rmd_loadData' in RMD


#Path to Seurat object (from previous analysis steps)
PATH_RDS_SEURAT_OBJECT = dir( file.path( PATH_PROJECT_OUTPUT, 
                                         "03a_clusteringScan_ambiguousInMyeloids",
                                         "initialClustering",
                                         "subGroup_APC", 
                                         "R0.4"),
                              pattern = ".*_seuratObject_allClusters.RDS",
                              full.names = TRUE); # Take the subset object from global experiment

# TSV file listing cell identities to be used for clustering representation.
# Format: barcodes as rownames and one column named "identity" (other columns
# allowed but ignored).
# Barcodes not found in Seurat object are ignored (with a warning), missing ones
# are grouped as class 'unknown' (with warning).
# Invalid path (or empty string) to use eventual 'Idents' from Seurat object.
#EXTERNAL_CLUSTERING_PATH = "" 
EXTERNAL_CLUSTERING_PATH = file.path( PATH_PROJECT_OUTPUT, 
                                      "04c_groupClusters_APC", 
                                      "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_groupClusters_APC_cellsClusterIdentity.tsv");
                                      # With renamed and eventually grouped clusters

# TSV file giving an hexadecimal-coded color for each cluster.
# Format: no header, cluster names as first column, colors in second column. All
# existing clusters must be defined in file (match by name, error otherwise). A
# warning is raised if more colors are declared than than existing clusters. 
# Empty string (or non-existing file path) for automatic coloring. 
EXTERNAL_CLUSTERSCOLOR_PATH = ""
#EXTERNAL_CLUSTERSCOLOR_PATH = file.path( PATH_PROJECT_EXTERNALDATA, 
#                                         "03_clusters_for_formattedDimReduc", 
#                                          "CpG", 
#                                          "total",
#                                          "clustersColor.tsv");

# TSV file giving an hexadecimal-coded color for each HTO.
# Format: no header, HTOs names as first column, colors in second column. All
# existing HTOs must be defined in file (match by name, error otherwise). A
# warning is raised if more colors are declared than than existing HTOs. 
# Empty string (or non-existing file path) for automatic coloring. 
EXTERNAL_HTOSCOLOR_PATH = ""


# Dimreduc coordinates from Seurat object to be used for plot: "umap" "tsne" or 
# "pca". Must already be computed and stored in loaded Seurat object, an error 
# is raised otherwise.
# Alternatively, a valid path to a TSV file containing cells coordinates from a
# previously computed 2-dimensions representation (dimensionality reduction). An
# error occurs for cells of Seurat object not found in file. A warning is raised
# for cells of file not found in Seurat object (cells ignored). As a consequence
# eventual subsetting MUST be made on Seurat object.
# Format: barcodes as rownames and two (named) columns for x and y coordinates.
CELLS_COORDINATES = dir( file.path( PATH_PROJECT_OUTPUT, 
                                    "03a_clusteringScan_ambiguousInMyeloids",
                                    "initialClustering",
                                    "subGroup_APC"),
                                    pattern = ".*_cellsCoordinates_umap.tsv",
                                    full.names = TRUE); # Use the umap computed from subclustering

# If 'CELLS_COORDINATES' is "pca", define which dimensions are used for 2D plots
PCA_DIMS = c( 1, 2);



#### HTO parameters

HTO_METADATA_COLUMN = "HTO_classification" # Name of the metadata column to use for HTO
HTO_FACTOR_LEVELS = c("D21-PBS", "D21-Panth", "D29-PBS", "D29-Panth"); # Order of levels to use for HTOs (for plots, must match with values in selected metadata column). NULL to ignore.



#### Filtering cells (to be excluded before current and further analyses)

# Cluster(s) name to remove, NULL to ignore
FILTER_CLUSTERS = NULL;

# File(s) (csv) containing cells barcode to remove (in a column named 'Cell', as 
# in 'cellsExplorer' exports), NULL to ignore
FILE_FILTER_CELLS = NULL;




#### General

# Seed for pseudo-random numbers
SEED = 42;

# Number of cores to use when possible (using 'future' for Seurat3, mclapply for
# other loops)
NBCORES = 4;

# Set the max global amount of "shared" memory for paralellization (future)
options(future.globals.maxSize= 891289600)

# Number of cells above which use ggplot instead of interactive plotly
PLOT_RASTER_NBCELLS_THRESHOLD = 10000;




#### Normalization and analysis parameters 




#### Plot options

# Scatterplot for quantitative (expression/score) values on dimreduc (ggplot)
PLOT_DIMREDUC_EXPRESSION_MAXONTOP  = TRUE;  # TRUE to plot most expressed cells on top, FALSE to plot least expressed cells on top, NULL to not reorder cells before plotting
PLOT_DIMREDUC_EXPRESSION_ALPHA     = 1;     # Alpha value (0-1) for points, 1 to mimic FeaturePlot and get maximum contrasts, less to prevent masking overlaid cells
PLOT_DIMREDUC_EXPRESSION_POINTSIZE = 0      # Size of points (0 for using Seurat:::AutoPointSize internal function)

# Scatterplot for group of cells (clusters/highlights) on dimreduc (ggplot)
PLOT_DIMREDUC_GROUPS_ALPHA     = 0.4;  # Alpha value (0-1) for points, 1 to mimic FeaturePlot and get maximum contrasts, less to prevent masking overlaid cells
PLOT_DIMREDUC_GROUPS_POINTSIZE = 1.5     # Size of points (0 for using Seurat:::AutoPointSize internal function)






#### Lists of genes of interest

# Number of genes under which genes of a module (MODULES_GENES) are transfered to be analyzed individually (MONITORED_GENES)
MONITORED_GENES_SMALL_MODULE = 5;
MODULES_CONTROL_SIZE = 100;


## Contamination-related genes (txt file, one gene name by line)
#CONTAMINATION_GENES = readLines( file.path( PATH_PROJECT_EXTERNALDATA, "ContaminationGenes.txt"))
CONTAMINATION_GENES=NULL

## Genes monitored individually (tsv file, one column for each group of genes)
#MONITORED_GENES = list()
MONITORED_GENES = as.list( read.table( file.path( PATH_PROJECT_REFERENCEDATA,
                                                  "GenesList_For_13_UMAPs",
                                                  "220916_GenesList_ForUmaps.tsv"),
                                       sep = "\t",
                                       header = TRUE,
                                       stringsAsFactors = FALSE,
                                       row.names = NULL, fill = TRUE));
MONITORED_GENES = Map('[', MONITORED_GENES, lapply(MONITORED_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

# MONITORED_GENES = c(list(Contamination = CONTAMINATION_GENES), MONITORED_GENES) # Add contamination genes as a group of genes to be monitored


## Genes monitored as modules (tsv file, one column for each group of genes)
# MODULES_GENES = as.list( read.table( file.path( PATH_PROJECT_REFERENCEDATA, "Modules.csv"),
#                                        sep = "\t",
#                                        header = TRUE,
#                                        stringsAsFactors = FALSE,
#                                        row.names = NULL, fill = TRUE));
# MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings
MODULES_GENES = MONITORED_GENES




