# CDNA1

This is an analysis of sample CDNA1 from second RUN.

## What are 'globalParams.R' and 'analysisParams.R' file ?

Both files are R scripts that aim to be sourced by other (analysis) scripts to provide pre-configured variables.

The file 'globalParams.R' defines general project-specific variables (mostly paths to input data and result folders).
The file 'analysisParams.R' is located in the folder of an analysis step and defines parameters specific to this analytic step (e.g. 01_QC).

## What is 'launch_reports_compilation.R' file ?

This file is a helper that starts the rendering of associated Rmd files. It is found next to Rmd files.
It is in charge of loading variables from 'globalParams.R' and 'analysisParams.R' before starting to render the associated Rmd file.

## Main commands

Most scripts reload data from results of previous steps using RDS files (and sometimes csv files).

Commands used to start each processing step (using container), using following environment variables:

NBPROC=4 # Number of parallel processing (when 'find' matches several scripts)


These command search for a matching script(s) to execute, initialize the container (mount volumes and set user using S6 bundled with Rstudio) and execute found script(s).

### 01_CellRanger_FeatureBarcoding

Started manually (bash) on a node with singularity and access to /mnt/DOSI shared folder.

### 02_QC

find -path '*02_QC*' -not -path '*archive*' -iname 'launch_reports_compilation.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c \
'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


## Helper commands

Delete temporary Rmd copies that can accumulate when process does not end and clean them on its own.
```
find -type f -iname 'tempRmdCopy_*'        -delete
```


