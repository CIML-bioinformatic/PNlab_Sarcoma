#!/bin/bash

# This script execute the CellRanger count analysis on JPGlab data
# with mode 'feature barcodes' to also assign HTO counts.
# It requires following informations:
# TAGLIST_FILE, LIBS_CSV_FILE, OUTPUT_FOLDER
# Path to reference and cellranger command (--id=) must fit correct genome (currently mm10)



# Path to the taglist file (CSV with HTO barcode sequences and corresponding names)
# MUST BE MADE WITH NEW FORMAT for "feature barcode" mode with columns (id, name, read, pattern, sequence, feature_type)
# See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis
TAGLIST_FILE="/mnt/DOSI/PNLAB/BIOINFO/Projet/scRNAseq_sarcoma/210921_VH00228_67_AAAMYKYM5_SARC_CD45_cDNA2/01_Reference/HTO_taglist.csv"

# Path to the file describing libraries (paths to fastq files and sample/HTO names)
LIBS_CSV_FILE="/mnt/DOSI/PNLAB/BIOINFO/Projet/scRNAseq_sarcoma/210921_VH00228_67_AAAMYKYM5_SARC_CD45_cDNA2/01_Reference/LibrariesDescription.csv"

# Number of expected cells (optionnal)
#EXPECTED_CELLS=3000

# Path to the CellRanger singularity image
SINGULARITY_IMAGE="/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/01_COMMON_DATA/02_CONTAINER/02_SINGULARITY/CellRanger/6_0_1/cellranger_601.simg"

# Folder where the CellRanger reference file of the genome are stored
#FOLDER_TO_REFERENCE="/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/01_COMMON_DATA/01_REFERENCE/01_GENOME/CellRanger/human/2020-A/refdata-gex-GRCh38-2020-A/"
FOLDER_TO_REFERENCE="/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/01_COMMON_DATA/01_REFERENCE/01_GENOME/CellRanger/mouse/2020-A/refdata-gex-mm10-2020-A"


# Output folder for the CellRanger output files
OUTPUT_FOLDER="//mnt/DOSI/PNLAB/BIOINFO/Project/scRNAseq_sarcoma/210921_VH00228_67_AAAMYKYM5_SARC_CD45_cDNA2/05_Output/01_CellRanger_FeatureBarcoding"
mkdir -p $OUTPUT_FOLDER



#### Execution of the analysis
echo FOLDER_TO_REFERENCE = $FOLDER_TO_REFERENCE
echo TAGLIST_FILE = $TAGLIST_FILE 
echo LIBS_CSV_FILE = $LIBS_CSV_FILE
echo OUTPUT_FOLDER = $OUTPUT_FOLDER


cd $OUTPUT_FOLDER


# --localcores 20
echo -e "\n\nCommand:\nsingularity exec -B /mnt/DOSI:/mnt/DOSI $SINGULARITY_IMAGE cellranger count --id=mm10 --libraries $LIBS_CSV_FILE --transcriptome=$FOLDER_TO_REFERENCE --feature-ref $TAGLIST_FILE --chemistry SC3Pv3\n\n"
singularity exec -B /mnt/DOSI:/mnt/DOSI $SINGULARITY_IMAGE cellranger count --id=mm10 --libraries $LIBS_CSV_FILE --transcriptome=$FOLDER_TO_REFERENCE --feature-ref $TAGLIST_FILE --chemistry SC3Pv3
