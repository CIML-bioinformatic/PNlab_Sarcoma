# Merge clusters from TSV file (taken from 01_QC) to be used in downstream analyses

library( funr)
library(forcats);


# Get path of current script
WORKING_DIR = dirname( sys.script())

# Get globalParams variables (project paths)
source( file.path( WORKING_DIR, "..", "globalParams.R"));

# Get path to original clustering output and read as data.frame
inputTSV = file.path( PATH_PROJECT_OUTPUT,
                      "03a_clusteringScan_ambiguousInMyeloids",
                      "initialClustering",
                      "subGroup_Myeloids",
                      "R0.6",
                      "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_initialClustering_Myeloids_R0.6_cellsClusterIdentities.tsv")
idents = read.csv( inputTSV, sep = "\t", row.names = 1);

# Recode factor levels
idents[["identity"]] = fct_recode( factor( idents[["identity"]]), 
                                   "M2-like TAM" = "0", 
                                   "Mono-TAM 3" = "1", 
                                   "Proliferating" = "2", 
                                   "Mono-TAM 2" = "3", 
                                   "Mono-TAM 1" = "4", 
                                   "Monocytes" = "5", 
                                   "TAM Inflammatory" = "6", 
                                   "TAM TNF NK-kb" = "7",
                                   "CAF" = "8",
                                   "MoDC" = "9",
                                   "Neutrophils" = "10");

# Write results in folder "external" with expected filename
outputTSV = file.path( PATH_PROJECT_OUTPUT,
                       "04b_groupClusters_Myeloids",
                       "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_groupClusters_Myeloids_cellsClusterIdentity.tsv")

# Make sure output directory exists
dir.create( dirname( outputTSV))

# Save cells cluster identity as determined with 'FindClusters'
write.table( idents, 
             file = outputTSV, 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");




