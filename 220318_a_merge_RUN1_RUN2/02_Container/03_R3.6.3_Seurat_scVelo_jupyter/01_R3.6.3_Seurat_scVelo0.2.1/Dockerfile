# Based on rocker - https://github.com/rocker-org/rocker-versioned
FROM rocker/r-ver:3.6.3

MAINTAINER RFenouil (fenouil@ciml.univ-mrs.fr)



# ##############################
# R packages system dependencies
# ##############################

RUN apt-get update && apt-get install -y optipng \
                                         pandoc \
                                         libpng-dev \
                                         zlib1g-dev \
                                         libssl-dev \
                                         libxml2-dev \
                                         libcurl4-openssl-dev


# ##########
# R PACKAGES 
# ##########

#### BiocManager
RUN Rscript -e 'install.packages( "BiocManager")'


#### Shiny
RUN Rscript -e 'install.packages( "shiny")'          # Various dependencies (htmlwidgets...)


#### RMarkdown & Reporting
RUN Rscript -e 'install.packages( "knitr")'
RUN Rscript -e 'install.packages( "rmarkdown")'
RUN Rscript -e 'install.packages( "bookdown")'
RUN Rscript -e 'install.packages( "htmltools")'      # browsable


#### Figures & layout management
# ggplot2
RUN Rscript -e 'install.packages( "ggplot2")'
RUN Rscript -e 'install.packages( "cowplot")'        # plot_grid, themes, ...
RUN Rscript -e 'install.packages( "ggpubr")'         # add_summary, geom_signif, ...
RUN Rscript -e 'install.packages( "ggrepel")'        # geom_text_repel, geom_label_repel
RUN Rscript -e 'install.packages( "gridExtra")'      # grid.arrange, ...
RUN Rscript -e 'BiocManager::install( "patchwork")'  # +/ operators for ggplots

# plotly
RUN Rscript -e 'install.packages( "plotly")'

# general
RUN Rscript -e 'install.packages( "gplots")'         # heatmap.2
RUN Rscript -e 'install.packages( "heatmaply")'      # heatmaply (interactive)
RUN Rscript -e 'BiocManager::install( "iheatmapr")'  # iheatmap (interactive, uses plotly), dependencies OK with BiocManager
RUN Rscript -e 'install.packages( "pheatmap")'       # pheatmap


#### Reporting
RUN Rscript -e 'install.packages( "DT")'             # datatable
RUN Rscript -e 'install.packages( "htmlwidgets")'    # JS (for datatable)
RUN Rscript -e 'install.packages( "kableExtra")'     # kable_styling
RUN Rscript -e 'install.packages( "pander")'         # pander


#### General
RUN Rscript -e 'install.packages( "digest")'         # digest (hashing)
RUN Rscript -e 'install.packages( "dplyr")'
RUN Rscript -e 'install.packages( "forcats")'        # fct_inorder (factors)
RUN Rscript -e 'install.packages( "funr")'           # get_script_path
RUN Rscript -e 'install.packages( "future")'         # plan (multicore)
RUN Rscript -e 'install.packages( "fs")'             # path_sanitize
RUN Rscript -e 'install.packages( "magrittr")'
RUN Rscript -e 'install.packages( "reshape")'        # melt
RUN Rscript -e 'install.packages( "scales")'         # hue_pal


#### Seurat
RUN Rscript -e 'BiocManager::install( "Seurat")'     # Dependencies OK with BiocManager (https://github.com/satijalab/seurat/issues/2409)
RUN Rscript -e 'BiocManager::install( "limma")'      # As recommended by Seurat for efficient Wilcoxon implementation 
RUN Rscript -e 'install.packages( "umap")'


#### scVelo
# Python3 for scVelo
RUN apt-get update && apt-get install -y bison \
                                         flex \
                                         git \
                                         libtool \
                                         libxml2-dev \
                                         python3 \
                                         python3-pip
# Optional libraries for scVelo
RUN pip3 install python-igraph louvain
RUN pip3 install pybind11 hnswlib
# Dependency (specific version as latest published one is not functional yet)
RUN pip3 install llvmlite==0.31.0
RUN pip3 install numba==0.50.1
# scVelo
RUN pip3 install -U scvelo==0.2.1


#### Velocyto
# samtools
RUN apt-get update && apt-get install -y samtools
# Fix permissions for numba (velocyto dependency) when container started with a specified user (https://github.com/numba/numba/issues/4032)
RUN mkdir /tmp/numba_cache \
    && chmod 777 /tmp/numba_cache 
ENV NUMBA_CACHE_DIR /tmp/numba_cache  
# Install command line velocyto for loom object creation
RUN pip3 install --upgrade cython
RUN pip3 install velocyto


#### Additional Seurat tools for interactions with Python and scVelo
RUN Rscript -e 'install.packages( "remotes")'        # Installing R packages from github
# seurat-wrappers (https://github.com/satijalab/seurat-wrappers)
RUN Rscript -e 'remotes::install_github( "satijalab/seurat-wrappers", upgrade = "never")'
# seurat-disk (https://github.com/mojaveazure/seurat-disk)
RUN apt-get update && apt-get install -y libhdf5-dev
RUN Rscript -e 'remotes::install_github( "mojaveazure/seurat-disk", upgrade = "never")'
# velocyto.R (https://github.com/velocyto-team/velocyto.R)
RUN Rscript -e 'BiocManager::install( "pcaMethods")' # Bioconductor dependency for velocyto.R
RUN apt-get update && apt-get install -y libboost-filesystem-dev \
                                         libboost-system-dev
RUN Rscript -e 'remotes::install_github( "velocyto-team/velocyto.R", upgrade = "never")'


