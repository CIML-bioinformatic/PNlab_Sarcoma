###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples. 
#


#### General

GLOBAL_DESCRIPTION = "Sarcoma scRNAseq - RUN2 CDNA1"

SCIENTIFIC_GROUP = "PNLAB"
SCIENTIFIC_PROJECT_NAME = "scRNAseq_sarcoma"
EXPERIMENT_PROJECT_NAME = "220214_VH00228_85_AAAT2JWM5_SARC_CD45_cDNA1"



#### Additional custom variables (for a specific analysis step, tools, ...)

#SAMPLE_ID = ""



#### Input / Output

# Output folder name in data folder (for R session object, lists of cells/genes) 
PATH_PROJECT = file.path( "/mnt/DOSI", 
                                        SCIENTIFIC_GROUP,
                                        "BIOINFO", 
                                        "Project",
                                        SCIENTIFIC_PROJECT_NAME,
                                        EXPERIMENT_PROJECT_NAME)

PATH_PROJECT_RAWDATA       = file.path( PATH_PROJECT, "00_RawData")
PATH_PROJECT_REFERENCEDATA = file.path( PATH_PROJECT, "01_Reference")
PATH_PROJECT_OUTPUT        = file.path( PATH_PROJECT, "05_Output")



# Create a 'safe' unique prefix for output files
outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME, "_",     
                            EXPERIMENT_PROJECT_NAME, "_"
                            #startTimeFileName, "_",
                            #paramsHash, "_"
                            )



#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



