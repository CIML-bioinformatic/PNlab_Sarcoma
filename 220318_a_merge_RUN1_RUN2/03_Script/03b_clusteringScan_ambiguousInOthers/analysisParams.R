###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "03b_clusteringScan_ambiguousInOthers"
PATH_ANALYSIS_OUTPUT = file.path( PATH_PROJECT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Clustering scan"

# Path to Seurat object from which to start analysis
PATH_RDS_SEURAT_OBJECT = dir( file.path( PATH_PROJECT_OUTPUT, "01a_QC"),
                              pattern = ".*seuratObject_final.RDS",
                              full.names = TRUE);


# TSV file listing cell identities to be used as clustering result for level 1.
# Useful when groups of interest (to process separately) were identified from a
# previously made analysis. Any cells left unidentified will be ignored (with a
# warning). Format: barcodes as rownames and one column named "identity" (other 
# columns ignored).
# Empty string to ignore and start computing 'new' clusters from level 1.
INITIAL_CLUSTERING_PATH = file.path( PATH_PROJECT_OUTPUT, 
                                     "02b_groupClusters_TwoOptionsForAmbiguousCells", 
                                     "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_ambiguousInOthers_cellsClusterIdentity.tsv");

# TSV file containing cells coordinates from a previously computed 2-dimension
# representation (e.g. dimensionality reduction). Will be used as reference for
# comparison with dimensionality reduction results at each level.
# Empty string to ignore. A warning is raised for cells not found in file.
ORIGINAL_CELLS_COORDINATES_PATH = dir( file.path( PATH_PROJECT_OUTPUT, "01a_QC"),
                                       pattern = ".*_cellsCoordinates_umap.tsv",
                                       full.names = TRUE);


# Scanning parameters
SUBCLUSTERING_RECURSION_NBLEVELS = 2; # How many recursive subclustering levels should be explored (Exponential: KEEP IT LOW !!!). 1 = no recursion/subclustering.
SUBCLUSTERING_RECURSION_MINCELLS = 50; # Number of cells in a cluster under which to stop subclustering 
FINDCLUSTERS_RESOLUTION_VALUES   = c(0.1, 0.2, 0.4, 0.6); # Clustering resolutions to be tested


# Parameters for re-running analyses (for levels > 1, use parameters above for level 1)
# NOTE: must set 'preventOverwrite=FALSE' in 'launch_reports_compilation.R' (or delete main html file) for being able to re-rnu an anlysis
USE_PRECOMPUTED_DIMREDUC_CLUSTERING = TRUE; # For reproduction of results, should subclustering use precomputed dimreduc and clusters when it finds corresponding 'csv' files in output folder ? To avoid inconsistencies, user must make sure that all parameters are identical between runs, and that csv files are provided for all levels.




#### Show distribution of categorical variable in clusters (mostly batch effect)

# Name of categorical (will be converted to factor) metadata column(s)
CATEGORICAL_IN_CLUSTERS = c( Batch = "orig.ident", 
                             HTO   = "HTO_classification")



#### General

# Seed for pseudo-random numbers
SEED = 42;

# Number of cores to use when possible (using 'future' for Seurat3, mclapply for
# recursion)
NBCORES = 2; # With R < 3.6.3 parallelization of umap crashes, not possible to paralellize recursive renderings
NBCORES_PLOT_EXTERNAL_FIGURES = 3; # Number of cores used in each reporting for plotting external 'png' figures (marker genes). Recommended to set to 1 when NBCORES>1.

# Set the max global amount of "shared" memory for paralellization (future)
options(future.globals.maxSize= 891289600)

# Number of cells above which use ggplot instead of interactive plotly
PLOT_RASTER_NBCELLS_THRESHOLD = 10000;

RENDEREXTERNALFIGURES_MARKERS   = TRUE; # Export as PNG files a larger number of dimreduc and violinplot for marker genes ?
RENDEREXTERNALFIGURES_MONITORED = TRUE; # Export as PNG files the monitored genes (not included in report anyway) ?
RENDEREXTERNALFIGURES_MODULES   = TRUE; # Same for modules

#### Filtering / Normalization
# Loading prefiltered object 

# Normalization parameters (see Seurat::NormalizeData())
DATA_NORM_METHOD      = "LogNormalize";
DATA_NORM_SCALEFACTOR = 10000;

# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE;
DATA_SCALE        = FALSE;
DATA_VARS_REGRESS = NULL;  # c("nCount_RNA") for UMIs (NULL to ignore)




#### Analysis parameters

# Maximum number of variable features to keep
VARIABLE_FEATURES_MAXNB   = 2000;  # For analysis (PCA)
VARIABLE_FEATURES_SHOWTOP = 200;   # For table in report

# Cluster identification parameters
FINDCLUSTERS_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results
FINDCLUSTERS_ALGORITHM      = 2;   # 1 = Louvain; 2 = Louvain with multilevel refinement; 3 = SLM; 4 = Leiden

# PCA parameters
PCA_NPC              = 50;  # Default number of dimensions to use for PCA (see Seurat::RunPCA())
PCA_PLOTS_NBDIMS     = 3;   # Number of dimensions to show in PCA-related plots
PCA_PLOTS_NBFEATURES = 15;  # Number of'top' features to show when plotting PCA loadings

# Dimensionality reduction parameters (UMAP)
DIMREDUC_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results


## Marker genes

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD           = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS          = TRUE;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT           = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR        = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR         = 0.001;    # PValue threshold for identification of significative markers


FINDMARKERS_SHOWTOP          = 100;      # Number of top marker genes to show in table (NULL for all)
FINDMARKERS_SHOWTOP_FIGS     = 10;       # Number of top marker genes to show in figures (dimreduc+violin). Smaller or equal to 'FINDMARKERS_SHOWTOP'.
FINDMARKERS_SHOWTOP_FIGS_EXT = 100;      # Number of top marker genes to plos as external figure files (dimreduc+violin). Smaller or equal to 'FINDMARKERS_SHOWTOP'.




#### Lists of genes of interest

# Number of genes under which genes of a module (MODULES_GENES) are transfered to be analyzed individually (MONITORED_GENES)
MONITORED_GENES_SMALL_MODULE = 1; # 1=disable
MODULES_CONTROL_SIZE = 100;


## Contamination-related genes (txt file, one gene name by line)
#CONTAMINATION_GENES = readLines( file.path( PATH_PROJECT_EXTERNALDATA, "ContaminationGenes.txt"))
CONTAMINATION_GENES=NULL

## Genes monitored individually (tsv file, one column for each group of genes)
#MONITORED_GENES = list()
MONITORED_GENES = as.list( read.table( file.path( PATH_PROJECT_REFERENCEDATA, "220331_GeneListCellClusterClean_plusCathepsin.tsv"),
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





#### Plot options

# Scatterplot for quantitative (expression/score) values on dimreduc (ggplot)
PLOT_DIMREDUC_EXPRESSION_MAXONTOP  = TRUE;  # TRUE to plot most expressed cells on top, FALSE to plot least expressed cells on top, NULL to not reorder cells before plotting
PLOT_DIMREDUC_EXPRESSION_ALPHA     = 1;     # Alpha value (0-1) for points, 1 to mimic FeaturePlot and get maximum contrasts, less to prevent masking overlaid cells
PLOT_DIMREDUC_EXPRESSION_POINTSIZE = 0      # Size of points (0 for using Seurat:::AutoPointSize internal function)

# Scatterplot for group of cells (clusters/highlights) on dimreduc (ggplot)
PLOT_DIMREDUC_GROUPS_ALPHA     = 0.4;  # Alpha value (0-1) for points, 1 to mimic FeaturePlot and get maximum contrasts, less to prevent masking overlaid cells
PLOT_DIMREDUC_GROUPS_POINTSIZE = 0     # Size of points (0 for using Seurat:::AutoPointSize internal function)


