FROM ubuntu:20.04

MAINTAINER Lionel Spinelli (lionel.spinelli@univ-amu.fr)

# CellRanger source have been downloaded from https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
# on April 14, 2021 - version Cell Ranger - 6.0.1 (April 13, 2021)
# Copy the cell ranger tarball into the image
COPY cellranger-6.0.1.tar.gz /opt

# Extract the tarball
WORKDIR /opt
RUN tar -xzvf cellranger-6.0.1.tar.gz
RUN rm cellranger-6.0.1.tar.gz

# Modify the PATH to get access to CellRanger
RUN echo 'export PATH=/opt/cellranger-6.0.1:$PATH' | tee -a /etc/profile
RUN ln -s /opt/cellranger-6.0.1/cellranger /usr/local/bin/cellranger

