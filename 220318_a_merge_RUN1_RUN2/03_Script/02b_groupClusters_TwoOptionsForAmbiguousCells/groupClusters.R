# Merge clusters from TSV file (taken from 01_QC) to be used in downstream analyses (subclustering)

library( funr)
library(forcats);


# Get path of current script
WORKING_DIR = dirname( sys.script())

# Get globalParams variables (project paths)
source( file.path( WORKING_DIR, "..", "globalParams.R"));

# Get path to original clustering output and read as data.frame
inputTSV = file.path( PATH_PROJECT_OUTPUT,
                      "01a_QC",
                      "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_cellsClusterIdentity.tsv")
idents = read.csv( inputTSV, sep = "\t", row.names = 1);

#table(idents[["identity"]])

# Get path to the CSV file containing identity of ambiguous cells (from cellsExplorer, see readme.md)
ambigousCells = read.csv( file.path( PATH_PROJECT_REFERENCEDATA,
                                     "cellsExplorer_manualSelectAmbiguousCells.csv"))
# Define ambiguous cells in original identity
idents[ambigousCells[["Cell"]], "identity"] = "ambiguous"

#table(idents[["identity"]])


# Recode factor levels (first option: ambiguous->Myeloids)
identity_ambiguousInMyeloids = fct_recode( factor( idents[["identity"]]), 
                                     "Lymphoids" = "2", 
                                     "Myeloids"  = "0", 
                                     "Myeloids"  = "1", 
                                     "Myeloids"  = "6", 
                                     "Myeloids"  = "ambiguous", 
                                     "APC"       = "4", 
                                     "APC"       = "5", 
                                     "Others"    = "3");


# Recode factor levels (second option: ambiguous->Others)
identity_ambiguousInOthers = fct_recode( factor( idents[["identity"]]), 
                                     "Lymphoids" = "2", 
                                     "Myeloids"  = "0", 
                                     "Myeloids"  = "1", 
                                     "Myeloids"  = "6", 
                                     "APC"       = "4", 
                                     "APC"       = "5", 
                                     "Others"    = "3",
                                     "Others"    = "ambiguous");


# Write results in folder "external" with expected filename
outputTSV_ambiguousInMyeloids = file.path( PATH_PROJECT_OUTPUT,
                                         "02b_groupClusters_TwoOptionsForAmbiguousCells",
                                         "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_ambiguousInMyeloids_cellsClusterIdentity.tsv")


# Write results in folder "external" with expected filename
outputTSV_ambiguousInOthers = file.path( PATH_PROJECT_OUTPUT,
                                       "02b_groupClusters_TwoOptionsForAmbiguousCells",
                                       "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_ambiguousInOthers_cellsClusterIdentity.tsv")


# Make sure output directory exists
dir.create( dirname( outputTSV_ambiguousInMyeloids))


# Update identity in 'idents' and save
idents[["identity"]] = identity_ambiguousInMyeloids;
write.table( idents, 
             file = outputTSV_ambiguousInMyeloids, 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

# Update identity in 'idents' and save
idents[["identity"]] = identity_ambiguousInOthers;
write.table( idents, 
             file = outputTSV_ambiguousInOthers, 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");




