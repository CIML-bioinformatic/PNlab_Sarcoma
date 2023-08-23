# Based on rocker - https://github.com/rocker-org/rocker-versioned
FROM rfenouil/r411_tidyverse_seurat4

MAINTAINER RFenouil (fenouil@ciml.univ-mrs.fr)


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
RUN pip3 install numba==0.51.2
# scVelo
RUN pip3 install -U git+https://github.com/theislab/scvelo.git


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
RUN Rscript -e 'BiocManager::install( "R.utils")' # dependency for seurat-wrappers
RUN Rscript -e 'remotes::install_github( "satijalab/seurat-wrappers", upgrade = "never")'
# seurat-disk (https://github.com/mojaveazure/seurat-disk)
RUN apt-get update && apt-get install -y libhdf5-dev
RUN Rscript -e 'remotes::install_github( "mojaveazure/seurat-disk", upgrade = "never")'
# velocyto.R (https://github.com/velocyto-team/velocyto.R)
RUN Rscript -e 'BiocManager::install( "pcaMethods")' # Bioconductor dependency for velocyto.R
RUN apt-get update && apt-get install -y libboost-filesystem-dev \
                                         libboost-system-dev
RUN Rscript -e 'remotes::install_github( "velocyto-team/velocyto.R", upgrade = "never")'


#### Add matplotlib-venn

RUN pip3 install matplotlib-venn
RUN pip3 install tabulate


