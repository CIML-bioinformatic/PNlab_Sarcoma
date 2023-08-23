# Load seurat object and subset it based on cell identity as specified in tsv file
# Used for downstream analysis when clustering is not taken as resulting from 
# subclustering script (which gives a subset Seurat object too).


library( funr)
library( forcats);
library( Seurat);


# Get path of current script
WORKING_DIR = dirname( sys.script())
# debug: WORKING_DIR = getwd()


# Get globalParams variables (project paths)
source( file.path( WORKING_DIR, "..", "globalParams.R"));

# Get path to original seurat object (take it from subclustering results)
inputSeurat = file.path( PATH_PROJECT_OUTPUT,
                         "03a_clusteringScan_ambiguousInMyeloids",
                         "initialClustering",
                         "subGroup_Lymphoids",
                         "R0.6",
                         "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_initialClustering_Lymphoids_R0.6_seuratObject_allClusters.RDS") # Could have taken the one with all cells from 01_QC but this one has new umap coordinates inside (no need to use external csv file for coordinates)

# Load Seurat object
sc10x = readRDS(inputSeurat)


# Get path to TSV file containing corresponding cells identities and load it
inputTSV = file.path( PATH_PROJECT_OUTPUT,
                      "04a_groupClusters_Lymphoids",
                      "scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_groupClusters_Lymphoids_cellsClusterIdentity.tsv")
cellsIdentity = read.csv( inputTSV, sep = "\t", row.names = 1);




### Define subgroups of interest
identsList = list();

identsList[["CD4"]] = (cellsIdentity[grepl("Treg|CD4", cellsIdentity[["identity"]]),])
identsList[["CD8"]] = (cellsIdentity[grepl("CD8", cellsIdentity[["identity"]]),])
identsList[["NK"]] = (cellsIdentity[grepl("NK", cellsIdentity[["identity"]]),])

# DEBUG: lapply(identsList, function(x){ table(x[["identity"]]) })




###
# Subset seurat object and save result

# Make sure output directory exists
outputFolder = file.path( PATH_PROJECT_OUTPUT,
                          "05a_subsetSeurat_Lymphoids_CD4_CD8_NK")
dir.create( outputFolder)

for(currentSetName in names(identsList))
{
  message(currentSetName)
  # Subset the original object using cell names
  resultSubset = sc10x[, rownames(identsList[[currentSetName]])]
  # Update identity (in case downstream don't want to load identity csv file separately)
  Idents(resultSubset) = identsList[[currentSetName]][["identity"]]
  
  # Save as binary file for downstream analyses
  saveRDS( object = resultSubset,
           file = file.path(outputFolder, paste0("scRNAseq_sarcoma_220318_a_merge_RUN1_RUN2_seuratObject_subset_", currentSetName, ".RDS")))
  
}



