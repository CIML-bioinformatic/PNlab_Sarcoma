# 220318_a_merge_RUN1_RUN2

This is a merged reanalysis of all samples (CDNA1 and CDNA2 = same biological material, separate libraries) for both RUN1 and RUN2 (separate experiments).

## What are 'globalParams.R' and 'analysisParams.R' file ?

Both files are R scripts that aim to be sourced by other (analysis) scripts to provide pre-configured variables.

The file 'globalParams.R' defines general project-specific variables (mostly paths to input data and result folders).
The file 'analysisParams.R' is located in the folder of an analysis step and defines parameters specific to this analytic step (e.g. 01_QC).

## What is 'launch_reports_compilation.R' file ?

This file is a helper that starts the rendering of associated Rmd files. It is found next to Rmd files.
It is in charge of loading variables from 'globalParams.R' and 'analysisParams.R' before starting to render the associated Rmd file.

## Main commands

Most scripts reload data from results of previous steps using RDS files (and sometimes csv files.

Commands used to start each processing step (using container), use following environment variable:
NBPROC=4 # Number of parallel processing (when 'find' matches several scripts)

These command search for a matching script(s) to execute, initialize the container (mount volumes and set user using S6 bundled with Rstudio) and execute found script(s).


### 01a_QC

find -path '*01a_QC*' -not -path '*archive*' -iname 'launch_reports_compilation.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c \
'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


### 01b_QC_HighRes

find -path '*01b_QC_HighRes*' -not -path '*archive*' -iname 'launch_reports_compilation.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c \
'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


### 02a_cellsExplorer

find -path '*02a_cellsExplorer*' -not -path '*archive*' -iname 'launch_reports_compilation.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c \
'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


### 02b_groupClusters_TwoOptionsForAmbiguousCells

find -path '*02b_groupClusters_TwoOptionsForAmbiguousCells*' -not -path '*archive*' -iname 'groupClusters.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c \
'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


### 03*_clusteringScan

find -path '*03*_clusteringScan*' -not -path '*archive*' -iname 'launch_reports_compilation.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c \
'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


### 04*_groupClusters_*

find -path '*04*_groupClusters_*' -not -path '*archive*' -iname 'groupClusters.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c \
'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


### 06*_compareConditions_*

find -path '*06*_compareConditions_*' -not -path '*archive*' -iname 'launch_reports_compilation.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c \
'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


### 08*_compareConditionsBH_*

find -path '*08*_compareConditionsBH_*' -not -path '*archive*' -iname 'launch_reports_compilation.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c \
'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


### 09*_geneLists_combinatorialExpression_*
find -path '*09*' -not -path '*archive*' -iname 'launch_reports_compilation.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c 'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


### 10b_geneLists_Myeloids
find -path '*10b*' -not -path '*archive*' -iname 'launch_reports_compilation.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c 'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


### 11_velocyto (Bash script, executed in container with velocito !)
find -path '*11_velocyto*' -not -path '*archive*' -iname 'execute_velocyto.sh' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c \
'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4_scvelo /init s6-setuidgid $(whoami) \
            ./$(readlink -m $0)'


### 12*_RNAvelocity*
# Only generates 'h5ad' files from loom and selected cells from seurat object (one for each HTO).
# The scVelo analysis and report is done in Jupyterlab using corresponding docker image.
find -path '*12*_RNAvelocity*' -not -path '*archive*' -iname 'launch_reports_compilation.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c 'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4_scvelo /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


### Interactive jupyterlab
# Start a webserver with jupyterlab to load corresponding notebooks
docker run -d -u $(id -u ${USER}):$(id -g ${USER}) \
           -p 8888:8888 \
           -e TOKEN=myPass \
           -v /mnt:/mnt \
           rfenouil/r3.6.3_seurat_scvelo0.2.1_jupyter:latest


### 13*_genesListUMAP*
find -path '*13*_genesListUMAP*' -not -path '*archive*' -iname 'launch_reports_compilation.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c 'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


### 14j_UMAP*
find -path '*14j_UMAP*' -not -path '*archive*' -iname 'launch_reports_compilation.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c 'docker run --rm -v /mnt:/mnt \
            -e PASSWORD=Pass -e USER=$(whoami) \
            -e USERID=$(id -u) -e GROUPID=$(id -g) \
            rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
            Rscript "$(readlink -m $0)"'


### 99 tests heatmaps
find -path '*99_tests_DE_heatmaps*' -not -path '*archive*' -iname 'launch_reports_compilation.R' -print0 | xargs -0 -t -L1 -P$NBPROC sh -c \
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


