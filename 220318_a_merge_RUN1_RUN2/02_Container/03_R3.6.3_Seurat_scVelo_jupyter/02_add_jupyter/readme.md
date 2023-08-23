
# Jupyter lab

This dockerfile imports the docker image of interest (FROM) and adds jupyter lab
 to it for ineractive tests.

This version is added to a R-based image so it also adds:
 - IRkenel for running R sessions in jupyter
 - rpy2 for running R code from python session

## Build

Here, jupyter layer is added to existing 'rfenouil/r3.6.3_seurat_scvelo0.2.1' Dockerfile
docker build . -t rfenouil/r3.6.3_seurat_scvelo0.2.1_jupyter

## Save

docker save rfenouil/r3.6.3_seurat_scvelo0.2.1_jupyter:latest | gzip > /mnt/DOSI/PNLAB/BIOINFO/Project/scRNAseq_sarcoma/220318_a_merge_RUN1_RUN2/02_Container/03_R3.6.3_Seurat_scVelo0.2.1_jupyter.tar.gz

## Run

Must forward a port for browser connection, and a token used for initial login:

```
docker run -d -u $(id -u ${USER}):$(id -g ${USER}) \
           -p 8888:8888 \
           -e TOKEN=myPass \
           -v /mnt:/mnt \
           rfenouil/r3.6.3_seurat_scvelo0.2.1_jupyter:latest
```

Then a browser must be used to connect locally: `https://127.0.0.1:8888`
