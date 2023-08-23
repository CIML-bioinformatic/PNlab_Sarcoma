# Merge clusters from TSV file (taken from 01_QC) to be used in downstream analyses

library( funr)
library(forcats);


# Get path of current script
WORKING_DIR = dirname( sys.script())

# Get globalParams variables (project paths)
source( file.path( WORKING_DIR, "..", "globalParams.R"));

# Get path to original clustering output and read as data.frame
inputTSV = file.path( PATH_PROJECT_OUTPUT,
                      "03b_clusteringScan_ambiguousInOthers",
                      "initialClustering",
                      "subGroup_Myeloids",
                      "R0.6",
                      "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_initialClustering_Myeloids_R0.6_cellsClusterIdentities.tsv")
idents = read.csv( inputTSV, sep = "\t", row.names = 1);

# Recode factor levels
idents[["identity"]] = fct_recode( factor( idents[["identity"]]), 
                                   "Mono-TAM 0" = "0", 
                                   "Mono-TAM TNF NFkb" = "1", 
                                   "M2-like TAM" = "2", 
                                   "Proliferating TAM" = "3", 
                                   "Mono-TAM 4" = "4", 
                                   "Inflammatory Mac" = "5", 
                                   "HSP - M2-Like pheno" = "6", 
                                   "Monocytes" = "7",
                                   "Prolif. Inflammatory Mac" = "8",
                                   "MoDC" = "9",
                                   "Neutrophils" = "10");

# Write results in folder "external" with expected filename
outputTSV = file.path( PATH_PROJECT_OUTPUT,
                       "04h_groupClusters_Myeloids_noAmbiguous",
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




