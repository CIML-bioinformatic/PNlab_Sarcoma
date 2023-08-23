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
                      "subGroup_Others",
                      "R0.2",
                      "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_initialClustering_Others_R0.2_cellsClusterIdentities.tsv")
idents = read.csv( inputTSV, sep = "\t", row.names = 1);

# Recode factor levels
idents[["identity"]] = fct_recode( factor( idents[["identity"]]), 
                                   "Tumor" = "0", 
                                   "Fibroblasts" = "1", 
                                   "Mastocytes" = "2", 
                                   "CAF" = "3");

# Write results in folder "external" with expected filename
outputTSV = file.path( PATH_PROJECT_OUTPUT,
                       "04d_groupClusters_Others",
                       "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_groupClusters_Others_cellsClusterIdentity.tsv")

# Make sure output directory exists
dir.create( dirname( outputTSV))

# Save cells cluster identity as determined with 'FindClusters'
write.table( idents, 
             file = outputTSV, 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");




