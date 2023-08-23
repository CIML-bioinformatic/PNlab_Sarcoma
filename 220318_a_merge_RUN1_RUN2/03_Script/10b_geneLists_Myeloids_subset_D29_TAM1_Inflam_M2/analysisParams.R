###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "10b_geneLists_Myeloids_subset_D29_TAM1_Inflam_M2"
PATH_ANALYSIS_OUTPUT = file.path( PATH_PROJECT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Myeloids - Genes lists heatmaps"


# !!! THIS ANALYSIS STEP NEEDS TO BE STARTED TWICE !!!
# On first time, it does all computation and generates figures as external files
# (using png/pdf and dev.off) which causes the rmarkdown layout to fail (figures
# in wrong tabset).
# On second execution, the RDATA result of previous run is loaded and previously
# generated external figures are not rendered again, but just integrated to
# rmarkdown using computed file path. It also skips time consuming executions
# (compute DEG analyses and enrichments).
# See chunk 'rmd_loadData' in RMD


# Path to R session object (from "compare conditions" previous analysis step).
# Will be loaded in a separate env to extract objects required for current step.
PATH_RDATA_PREVIOUS_RESULTS = file.path( PATH_PROJECT_OUTPUT, 
                                         "08b_compareConditionsBH_Myeloids", 
                                         "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_sessionImage_final.RDATA");

# List of genes to be analyzed / plotted
GENES_LIST = as.list( read.table( file.path( PATH_PROJECT_REFERENCEDATA, 
                                            "geneLists",
                                            "220725_TAM_signature.tsv"),
                                  sep = "\t",
                                  header = TRUE,
                                  stringsAsFactors = FALSE,
                                  row.names = NULL, 
                                  fill = TRUE));

# Remove empty elements 
GENES_LIST = Map('[', GENES_LIST, lapply(GENES_LIST, function(x){ which( nchar( x)>0)})); # Remove empty strings

HTO_FACTOR_LEVELS = c("D21-PBS", "D21-Panth", "D29-PBS", "D29-Panth"); # Order of levels to use for HTOs (for plots, must match with actual values).


# Addition in this script version only to subset populations for heatmaps (this version only)
HTO_SELECT = c("D29")
CLUSTER_SELECT = c("Mono.*TAM 1", ".*Inflammatory.*", "M2.*")



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

#### Analyses parameters

MODULES_CONTROL_SIZE = 100 # Number of genes to sample as control for each group of 'module score'

# Factor levels used for ordering clusters in heatmaps
HEATMAPS_CLUSTERS_ORDERING = c("Neutrophils", "Monocytes", "TAM Inflammatory", "Mono-TAM 1", "Mono-TAM 2", "Mono-TAM 3", "TAM TNF NK-kb", "M2-like TAM", "Proliferating", "MoDC", "CAF")


