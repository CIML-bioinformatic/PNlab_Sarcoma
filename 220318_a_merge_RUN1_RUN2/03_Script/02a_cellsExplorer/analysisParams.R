###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "02a_cellsExplorer"
PATH_ANALYSIS_OUTPUT = file.path( PATH_PROJECT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Cells explorer"

#Path to Seurat object from which to start analysis
PATH_RDS_SEURAT_OBJECT = dir( file.path( PATH_PROJECT_OUTPUT, "01a_QC"),
                              pattern = ".*seuratObject_final.RDS",
                              full.names = TRUE);

# TSV file listing cell identities to be used for clustering representation.
# Format: barcodes as rownames and one column named "identity" (other columns
# allowed but ignored).
# Empty string to ignore and use eventual 'Idents' from Seurat object.
EXTERNAL_CLUSTERING_PATH = "";

# TSV file giving an hexadecimal-coded color for each cluster.
# Format: no header, cluster names as first column, colors in second column. All
# existing clusters must be defined in file (match by name, error otherwise). A
# warning is raised if more colors are declared than than existing clusters. 
# Empty string (or non-existing file path) for automatic coloring. 
EXTERNAL_CLUSTERSCOLOR_PATH = "";

# Dimreduc coordinates from Seurat object to be used for plot: "umap" "tsne" or 
# "pca". Must already be computed and stored in loaded Seurat object, an error 
# is raised otherwise.
# Alternatively, a valid path to a TSV file containing cells coordinates from a
# previously computed 2-dimensions representation (dimensionality reduction). A
# warning is raised for cells of Seurat object not found in file.
# Format: barcodes as rownames and two columns for x and y coordinates.
CELLS_COORDINATES = dir( file.path( PATH_PROJECT_OUTPUT, "01a_QC"),
                                    pattern = ".*_cellsCoordinates_umap.tsv",
                                    full.names = TRUE);

# If 'CELLS_COORDINATES' is "pca", define which dimensions are used for 2D plots
PCA_DIMS = c( 1, 2);




#### General

# Seed for pseudo-random numbers
SEED = 42;

# Number of cores to use when possible (using 'future' for Seurat3, mclapply for
# other loops)
NBCORES = 4;




#### Analysis parameters



