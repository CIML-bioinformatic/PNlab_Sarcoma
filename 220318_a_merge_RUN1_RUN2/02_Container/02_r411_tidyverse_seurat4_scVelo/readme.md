
# r411_tidyverse_seurat4_scvelo

Based on rocker/tidyverse distribution (includes Rstudio) from previous container.
Add:
 - scVelo (latest from git) + dependencies
 - velocyto
 - RNA velocity R wrappers and utils
 - 'matplotlib-venn' and 'tabulate' (python libs)



## Build

docker build . -t rfenouil/r411_tidyverse_seurat4_scvelo



## Save

docker save rfenouil/r411_tidyverse_seurat4_scvelo:latest | gzip > /mnt/DOSI/PNLAB/BIOINFO/Project/scRNAseq_sarcoma/220318_a_merge_RUN1_RUN2/02_Container/02_r4.1.1_tidyverse_seurat4_scvelo.tar.gz


## Run Rstudio

Give user details to internal script which sets user and permissions:

```
docker run --rm -d \
           --name containerName \
           -p 8787:8787 \
           -e PASSWORD=yourPass \
           -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) \
           -v /mnt:/mnt \
           rfenouil/r411_tidyverse_seurat4_scvelo
```

Then connect to the machine running docker (localhost) on mapped port (8787):
http://127.0.0.1:8787



## Run a command as specified user (not starting Rstudio)

docker run --rm \
           -e PASSWORD=yourPass \
           -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) \
           -v /mnt:/mnt \
           rfenouil/r411_tidyverse_seurat4_scvelo \
           /init s6-setuidgid $(whoami) command

