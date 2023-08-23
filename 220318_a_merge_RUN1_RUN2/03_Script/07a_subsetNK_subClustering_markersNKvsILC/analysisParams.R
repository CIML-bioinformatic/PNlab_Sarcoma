###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "07a_subsetNK_subClustering_markersNKvsILC"
PATH_ANALYSIS_OUTPUT = file.path( PATH_PROJECT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "subClustering NK + markers NK vs ILC"

#Path to Seurat object (from previous analysis steps). Take the NK subset object from Lymphoid set.
PATH_RDS_SEURAT_OBJECT = file.path( PATH_PROJECT_OUTPUT, 
                                    "05a_subsetSeurat_Lymphoids_CD4_CD8_NK",
                                    "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_seuratObject_subset_NK.RDS")

# TSV file listing cell identities to be used for clustering representation.
# Format: barcodes as rownames and one column named "identity" (other columns
# allowed but ignored).
# Barcodes not found in Seurat object are ignored (with a warning), missing ones
# are grouped as class 'unknown' (with warning).
# Invalid path (or empty string) to use eventual 'Idents' from Seurat object.
#EXTERNAL_CLUSTERING_PATH = "" 
EXTERNAL_CLUSTERING_PATH = file.path( PATH_PROJECT_OUTPUT, 
                                      "04a_groupClusters_Lymphoids", 
                                      "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_groupClusters_Lymphoids_cellsClusterIdentity.tsv");
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
                                    "subGroup_Lymphoids"),
                                    pattern = ".*_cellsCoordinates_umap.tsv",
                                    full.names = TRUE); # Use the umap computed from subclustering

# If 'CELLS_COORDINATES' is "pca", define which dimensions are used for 2D plots
PCA_DIMS = c( 1, 2);



#### HTO parameters

HTO_METADATA_COLUMN = "HTO_classification" # Name of the metadata column to use for HTO
HTO_FACTOR_LEVELS = c("D21-PBS", "D21-Panth", "D29-PBS", "D29-Panth"); # Order of levels to use for HTOs (for plots, must match with values in selected metadata column). NULL to ignore.

## HTO differential expression analysis
# For each cluster, a differential analysis is made on conditions defined by HTOs.

# Define the list of DE comparisons to be made (each being a 2-elements list of individual conditions names to compare)
HTO_DIFFEXP_COMPARISONLIST = list("D21_PBS_vs_Panth" = list("D21_PBS" = c("D21-PBS"), "D21_Panth" = c("D21-Panth")),
                               "D29_PBS_vs_Panth" = list("D29_PBS" = c("D29-PBS"), "D29_Panth" = c("D29-Panth")),
                               "D21D29_PBS_vs_Panth" = list("D21D29_PBS" = c("D21-PBS", "D29-PBS"), "D21D29_Panth" = c("D21-Panth", "D29-Panth")),
                               "PBS_D21_vs_D29"   = list("D21_PBS" = c("D21-PBS"), "D29_PBS" = c("D29-PBS")),
                               "Panth_D21_vs_D29" = list("D21_Panth" = c("D21-Panth"), "D29_Panth" = c("D29-Panth")))
# Need a proper model (edger ?) to get the how Panth vs PBS is different between D21 and D29 ?? 

# Parameters for differential expression analysis (see Seurat::FindMarkers())
HTO_DIFFEXP_METHOD    = "wilcox"  # Method
HTO_DIFFEXP_ONLYPOS   = FALSE;    # Only consider overexpressed annotations ? (if FALSE downregulated genes can also be markers)
HTO_DIFFEXP_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
HTO_DIFFEXP_LOGFC_THR = 0.1;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
HTO_DIFFEXP_PVAL_THR  = 0.05;    # PValue threshold for identification of significative differentially expressed genes
# There is no 'double-dipping' problem here, we observe much more selective p-values than for cluster markers

HTO_DIFFEXP_SELECT_TOP = NULL; # Select top nb of differentially expressed genes for each condition and selected cluster(s)/identity (sorted by adjusted PValue as FindMarkers output). NULL for no selection.




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
# Just for recomputing marker genes of clusters (+HTOs) defined with input files. 
# There is no reclustering and no recomputation of dimensionality reduction here.
# In theory, if the Seurat object comes from a previously computed report or a
# 'clustering scan' subset, normalization, variable genes and PCA should already
# be computed, but we do it again here in case subset was made manually or one
# wants to change following parameters (again, just for searching markers...).

# Normalization parameters (see Seurat::NormalizeData())
DATA_NORM_METHOD      = "LogNormalize";
DATA_NORM_SCALEFACTOR = 10000;

# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE;
DATA_SCALE        = FALSE;
DATA_VARS_REGRESS = NULL;  # c("nCount_RNA") for UMIs (NULL to ignore)


# Maximum number of variable features to keep
VARIABLE_FEATURES_MAXNB   = 2000;  # For analysis (PCA)
VARIABLE_FEATURES_SHOWTOP = 200;   # For table in report

# Cluster identification parameters
FINDCLUSTERS_RESOLUTION     = 0.5;
FINDCLUSTERS_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results
FINDCLUSTERS_ALGORITHM      = 1;   # 1 = Louvain; 2 = Louvain with multilevel refinement; 3 = SLM; 4 = Leiden

# PCA parameters
PCA_NPC              = 50;  # Default number of dimensions to use for PCA (see Seurat::RunPCA())
PCA_PLOTS_NBDIMS     = 3;   # Number of dimensions to show in PCA-related plots
PCA_PLOTS_NBFEATURES = 15;  # Number of'top' features to show when plotting PCA loadings

# No recomputation of dimensionality reduction here, we use what is given as input (from object or tsv file)
# # Dimensionality reduction parameters (TSNE/UMAP)
# DIMREDUC_USE_PCA_NBDIMS = 20;  # Number of dimensions to use from PCA results


#### Plot options

# Scatterplot for quantitative (expression/score) values on dimreduc (ggplot)
PLOT_DIMREDUC_EXPRESSION_MAXONTOP  = TRUE;  # TRUE to plot most expressed cells on top, FALSE to plot least expressed cells on top, NULL to not reorder cells before plotting
PLOT_DIMREDUC_EXPRESSION_ALPHA     = 1;     # Alpha value (0-1) for points, 1 to mimic FeaturePlot and get maximum contrasts, less to prevent masking overlaid cells
PLOT_DIMREDUC_EXPRESSION_POINTSIZE = 0      # Size of points (0 for using Seurat:::AutoPointSize internal function)

# Scatterplot for group of cells (clusters/highlights) on dimreduc (ggplot)
PLOT_DIMREDUC_GROUPS_ALPHA     = 0.4;  # Alpha value (0-1) for points, 1 to mimic FeaturePlot and get maximum contrasts, less to prevent masking overlaid cells
PLOT_DIMREDUC_GROUPS_POINTSIZE = 1.5     # Size of points (0 for using Seurat:::AutoPointSize internal function)




#### Lists of genes from NK vs ILC1 analysis (EVlab) to be visualized

TOPGENES_NK_ILC = 20; # Number of top genes for modules (in addition to total lists)

## List of genes separating NK and ILCs (from diff expression)
# DEG NK vs ILCs from a previous analysis (EVlab)
# Cutoff 0.25 on avg_logFC. Contains only (very) significative genes. 
genes_NKvsILC = read.csv( file.path( PATH_PROJECT_REFERENCEDATA, 
                                       "subsetNK_subClustering_markersNKvsILC",
                                       "01_NK_vs_ILC1_tableDEG.csv"))

# Sort by PValue (in case it is not already)
genes_NKvsILC = genes_NKvsILC[order(genes_NKvsILC[["p_val_adj"]]),]

# Separate NK and ILC rows
genenamesNK = genes_NKvsILC[genes_NKvsILC[["avg_logFC"]]>0, "X"]
genenamesILC = genes_NKvsILC[genes_NKvsILC[["avg_logFC"]]<0, "X"]

# Create list of gene names
degList_NK_ILC = list("degEV_NK"   = genenamesNK, 
                      "degEV_ILC1" = genenamesILC)

# Also create list of top genes
degList_NK_ILC_top = lapply(degList_NK_ILC, head, TOPGENES_NK_ILC)
names(degList_NK_ILC_top) = paste0(names(degList_NK_ILC), "_top", TOPGENES_NK_ILC)


## List of genes classifying ILCs (markers from clustering)
# ILC (nad NK) From a previous analysis (EVlab).
# Contains only (very) significative genes as of adjusted P-value. 
markersTable_NK_ILC = read.csv( file.path( PATH_PROJECT_REFERENCEDATA, 
                                      "subsetNK_subClustering_markersNKvsILC",
                                      "02_markers_NK_0_2_3_ILC_1_4.csv"))

# Rename cluster numbers
markersTable_NK_ILC[["cluster"]] = forcats::fct_recode( factor( markersTable_NK_ILC[["cluster"]]), 
                                           "markersEV_NK_Dim"     = "0", 
                                           "markersEV_ILC1_1"     = "1", 
                                           "markersEV_NK_Bright"  = "2", 
                                           "markersEV_NK_1"       = "3", 
                                           "markersEV_ILC1_2"     = "4");

# Sort each cluster by PValue in case not already sorted and axtract as list of gene names
markersList_NK_ILC = by( markersTable_NK_ILC, 
                         markersTable_NK_ILC[["cluster"]], 
                         function(x){ return(x[order(x[["p_val_adj"]]),"gene"]) } )

# Also create list of top genes
markersList_NK_ILC_top = lapply(markersList_NK_ILC, head, TOPGENES_NK_ILC)
names(markersList_NK_ILC_top) = paste0(names(markersList_NK_ILC), "_top", TOPGENES_NK_ILC)




### Lists of genes of interest

# Number of genes under which genes of a module (MODULES_GENES) are transfered to be analyzed individually (MONITORED_GENES)
MONITORED_GENES_SMALL_MODULE = 1; # 1 = disable
MODULES_CONTROL_SIZE = 100;

# Replace typical list of monitored genes by ILC/NK DEG list and markers from EV
MONITORED_GENES = c(degList_NK_ILC, markersList_NK_ILC)
MONITORED_GENES = MONITORED_GENES[order(names(MONITORED_GENES))] # Order by name

# For modules also add the top genes selection (see 'TOPGENES_NK_ILC') 
MODULES_GENES = c(degList_NK_ILC, degList_NK_ILC_top, markersList_NK_ILC, markersList_NK_ILC_top)
MODULES_GENES = MODULES_GENES[order(names(MODULES_GENES))] # Order by name

