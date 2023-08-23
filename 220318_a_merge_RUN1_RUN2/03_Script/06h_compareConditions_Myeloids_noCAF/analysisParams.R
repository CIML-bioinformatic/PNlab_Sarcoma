###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "06h_compareConditions_Myeloids_noCAF"
PATH_ANALYSIS_OUTPUT = file.path( PATH_PROJECT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Myeloids no CAF - Compare conditions"


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
                                         "03b_clusteringScan_ambiguousInOthers",
                                         "initialClustering",
                                         "subGroup_Myeloids", 
                                         "R0.6"),
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
                                      "04h_groupClusters_Myeloids_noAmbiguous", 
                                      "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_groupClusters_Myeloids_cellsClusterIdentity.tsv");
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
                                    "03b_clusteringScan_ambiguousInOthers",
                                    "initialClustering",
                                    "subGroup_Myeloids"),
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
# HTO_DIFFEXP_COMPARISONLIST = list("D21_PBS_vs_Panth" = list("D21_PBS" = c("D21-PBS"), "D21_Panth" = c("D21-Panth")),
#                                   "D29_PBS_vs_Panth" = list("D29_PBS" = c("D29-PBS"), "D29_Panth" = c("D29-Panth")),
#                                   "D21D29_PBS_vs_Panth" = list("D21D29_PBS" = c("D21-PBS", "D29-PBS"), "D21D29_Panth" = c("D21-Panth", "D29-Panth")),
#                                   "PBS_D21_vs_D29"   = list("D21_PBS" = c("D21-PBS"), "D29_PBS" = c("D29-PBS")),
#                                   "Panth_D21_vs_D29" = list("D21_Panth" = c("D21-Panth"), "D29_Panth" = c("D29-Panth")))
# Reversed model for a better logic in figures representation (+reordered)
HTO_DIFFEXP_COMPARISONLIST = list("D21D29_Panth_vs_PBS" = list("D21D29_Panth" = c("D21-Panth", "D29-Panth"), "D21D29_PBS" = c("D21-PBS", "D29-PBS")),
                                  "D21_Panth_vs_PBS" = list("D21_Panth" = c("D21-Panth"), "D21_PBS" = c("D21-PBS")),
                                  "D29_Panth_vs_PBS" = list("D29_Panth" = c("D29-Panth"), "D29_PBS" = c("D29-PBS")),
                                  "PBS_D29_vs_D21"   = list("D29_PBS" = c("D29-PBS"), "D21_PBS" = c("D21-PBS")),
                                  "Panth_D29_vs_D21" = list("D29_Panth" = c("D29-Panth"), "D21_Panth" = c("D21-Panth")))
# Need a proper model (edger ?) to get the how Panth vs PBS is different between D21 and D29 ?? 

# Parameters for differential expression analysis (see Seurat::FindMarkers())
HTO_DIFFEXP_METHOD    = "wilcox"  # Method
HTO_DIFFEXP_ONLYPOS   = FALSE;    # Only consider overexpressed annotations ? (if FALSE downregulated genes can also be markers)
HTO_DIFFEXP_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
HTO_DIFFEXP_LOGFC_THR = 0.01;    # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
HTO_DIFFEXP_PVAL_THR  = 0.05;    # PValue threshold for identification of significative differentially expressed genes
# There is no 'double-dipping' problem here, we observe much more selective p-values than for cluster markers
HTO_DIFFEXP_FORCE_ADJUST_PVAL_BH = FALSE # Default is Bonferroni (conservative method), TRUE to recompute using BH

# Selection filters for DEGs
HTO_DIFFEXP_SELECT_TOP = NULL; # Select top nb of differentially expressed genes for each condition and selected cluster(s)/identity (sorted by adjusted PValue as FindMarkers output). NULL for no selection. Subset selection used for summary datatable (prior to dedicated filter 'HTO_DIFFEXP_TOP_DATATABLE'), enrichment analyses, heatmaps (prior to dedicated filter 'PLOT_HEATMAPS_MAXGENES'), NOT for volcano/MA plots, NOT for external csv file output.
HTO_DIFFEXP_TOP_DATATABLE = 100; # Select top nb of differentially expressed genes reported for each cluster/identity in report datatable. NULL for no selection. Filter applied on eventually pre-filtered top DEGs (see 'HTO_DIFFEXP_SELECT_TOP')




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

# No reclustering here, we use what is given as input (from object or tsv file)
# # Cluster identification parameters 
# FINDCLUSTERS_USE_PCA_NBDIMS = 20;  # Number of dimensions to use from PCA results
# FINDCLUSTERS_ALGORITHM      = 2;   # 1 = Louvain; 2 = Louvain with multilevel refinement; 3 = SLM; 4 = Leiden

# PCA parameters
PCA_NPC              = 50;  # Default number of dimensions to use for PCA (see Seurat::RunPCA())
PCA_PLOTS_NBDIMS     = 3;   # Number of dimensions to show in PCA-related plots
PCA_PLOTS_NBFEATURES = 15;  # Number of'top' features to show when plotting PCA loadings

# No recomputation of dimensionality reduction here, we use what is given as input (from object or tsv file)
# # Dimensionality reduction parameters (TSNE/UMAP)
# DIMREDUC_USE_PCA_NBDIMS = 20;  # Number of dimensions to use from PCA results

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP   = 15;       # Number of marker genes to show in report and tables (NULL for all)




#### Plot options

# Scatterplot for quantitative (expression/score) values on dimreduc (ggplot)
PLOT_DIMREDUC_EXPRESSION_MAXONTOP  = TRUE;  # TRUE to plot most expressed cells on top, FALSE to plot least expressed cells on top, NULL to not reorder cells before plotting
PLOT_DIMREDUC_EXPRESSION_ALPHA     = 1;     # Alpha value (0-1) for points, 1 to mimic FeaturePlot and get maximum contrasts, less to prevent masking overlaid cells
PLOT_DIMREDUC_EXPRESSION_POINTSIZE = 0      # Size of points (0 for using Seurat:::AutoPointSize internal function)

# Scatterplot for group of cells (clusters/highlights) on dimreduc (ggplot)
PLOT_DIMREDUC_GROUPS_ALPHA     = 0.4;  # Alpha value (0-1) for points, 1 to mimic FeaturePlot and get maximum contrasts, less to prevent masking overlaid cells
PLOT_DIMREDUC_GROUPS_POINTSIZE = 1.5     # Size of points (0 for using Seurat:::AutoPointSize internal function)


PLOT_HEATMAPS_MAXGENES = 50 # Limit the number of genes to show on heatmaps (top sorted by adjusted PValue). NULL for no selection.




#### Lists of genes of interest

# Number of genes under which genes of a module (MODULES_GENES) are transfered to be analyzed individually (MONITORED_GENES)
MONITORED_GENES_SMALL_MODULE = 5;
MODULES_CONTROL_SIZE = 100;


## Contamination-related genes (txt file, one gene name by line)
#CONTAMINATION_GENES = readLines( file.path( PATH_PROJECT_EXTERNALDATA, "ContaminationGenes.txt"))
CONTAMINATION_GENES=NULL

## Genes monitored individually (tsv file, one column for each group of genes)
#MONITORED_GENES = list()
MONITORED_GENES =  as.list( read.table( file.path( PATH_PROJECT_REFERENCEDATA, 
                                                   "geneLists",
                                                   "220725_TAM_signature.tsv"),
                                        sep = "\t",
                                        header = TRUE,
                                        stringsAsFactors = FALSE,
                                        row.names = NULL, 
                                        fill = TRUE));

MONITORED_GENES = Map('[', MONITORED_GENES, lapply(MONITORED_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

# MONITORED_GENES = c(list(Contamination = CONTAMINATION_GENES), MONITORED_GENES) # Add contamination genes as a group of genes to be monitored


## Genes monitored as modules (tsv file, one column for each group of genes)
# MODULES_GENES = as.list( read.table( file.path( PATH_PROJECT_REFERENCEDATA, "Modules.csv"),
#                                        sep = "\t",
#                                        header = TRUE,
#                                        stringsAsFactors = FALSE,
#                                        row.names = NULL, fill = TRUE));
# MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings
MODULES_GENES = list()




#### Functional enrichment analyses

FUNCTIONAL_ORGANISM_LIBRARY = "org.Mm.eg.db" # "org.Hs.eg.db" for Homo Sapiens. Bioconductor annotation for current species (convert gene symbol to entrez id)
FUNCTIONAL_ORGANISM_CODE = "mmu" # "hsa" for Homo Sapiens (see http://www.genome.jp/kegg/catalog/org_list.html) 

TOPTERMS_SUMMARY = 3    # Number of top terms (based on 'p.adjust') to show for each 'sample' (summary figures with genes names)
SEPARATE_ALLCELLS_ENRICHMENTS = TRUE; # Logical. When plotting functionnal enrichent summary, do a separate plot for 'AllCells' and cluster groups
SEPARATE_ALLCELLS_TOPTERMS = 10; # If 'AllCells' is plotted separately (see SEPARATE_ALLCELLS_ENRICHMENTS), eventually specify a different number for top terms to be shown
TOPTERMS_FIGURE    = 20   # Number of top terms (based on 'p.adjust') to show for each 'sample' (individual sample/cluster figures)
TOPTERMS_DATATABLE = 20   # Number of top terms for datatable in report (full table exported as separate file)

## KEGG functionnal analysis (clusterprofiler)
KEGG_PADJUST_METHOD   = "BH"
KEGG_PVALUE_CUTOFF    = 0.05

## GO functionnal analysis (clusterprofiler)
GO_CATEGORIES         = list( "GO-Bio.Process"      = "BP", # GO categories to check enrichments for (remove unused elements).
                              "GO-Cell.Compartment" = "CC", # Elements names shown on the report (can be changed).
                              "GO-Molec.Function"   = "MF") # Elements values used for enrichment function (do not change).
GO_PADJUST_METHOD     = "BH"
GO_PVALUE_CUTOFF      = 0.05
GO_SIMPLIFY_CUTOFF    = NULL # Cutoff for eventual simplification of redundant GO terms in results (NULL to ignore simplification)

#UNIVERSE_IS_UNIONGENES = TRUE # Controlled by Rmd file (both values). Define if the background for enrichment analysis must be defined as the union of genes from all clusters (all genes otherwise)

## GSEA

GSEA_PVALUE_CUTOFF = 0.05

## GSEA: Mitocarta pathways from downloaded excel file to perform GSEA enrichment analysis against mitochondrial functions
MITOCARTA_DATA = read.csv( file.path( PATH_PROJECT_REFERENCEDATA, "Mitocarta", "Mouse.MitoCarta3.csv" ),
                           sep = "\t",
                           quote = '"')
MITOCARTA_PATHWAYS = strsplit( MITOCARTA_DATA[["Genes"]], ", ")
names(MITOCARTA_PATHWAYS) = MITOCARTA_DATA[["MitoPathway"]]

## GSEA: gmt files to be tested
MSIGDB_GMT_PATHS = dir( file.path( PATH_PROJECT_REFERENCEDATA, 
                                  "msigdb_gmt",
                                  "mus_musculus"),
                       full.names = TRUE)

## KEGG SELECTED PATHWAYS OF INTEREST

# IDs of manually selected pathways of interest (in addition to enriched ones)
KEGG_MANUAL_PATHWAY_IDS = c("04650", "04612", "00010", "00020", "01212", "00190", "00071", "00770", "04064", "04668", "04066", "04350", "04620", "04624", "04621", "04622", "04623", "04625", "04650", "04612")


