#### Add jupyter lab to any base container for interactive tests
# See:
# - github.com/jupyter/docker-stacks/blob/master/base-notebook/Dockerfile
# - github.com/jupyter/docker-stacks/blob/master/base-notebook/start.sh
#
# Start command example:
# docker run -d -u $(id -u ${USER}):$(id -g ${USER}) \
#            -p 8888:8888 \
#            -e TOKEN=myPass \
#            -v /mnt:/mnt \
#            rfenouil/jpglab_peyerspatch_phagocytes_r3.6.3_seurat_scvelo_jupyter:latest
#


#### Select base image to which we add jupyter lab
FROM rfenouil/r3.6.3_seurat_scvelo0.2.1

LABEL maintainer="RFenouil (fenouil@ciml.univ-mrs.fr)"



#### Make sure Python and pip are installed
RUN apt-get update && apt-get install -y python3 \
                                         python3-pip



#### Create a 'dummy' user for installing jupyterlab that will be overriden with
# actual user at runtime (docker run -u).
# We run container as user with correct uid/gid to respect NFS permissions.
# Jupyterlab cannot be run as a user if it has been installed by 'root'.
ENV NB_USER="jupyter_user"
ENV NB_UID="1000"
ENV NB_GID="100"
ENV HOME=/home/$NB_USER

RUN useradd -m -s /bin/bash -N -u $NB_UID $NB_USER



#### Switch to created user, install jupyter lab and update PATH
USER $NB_UID
RUN pip3 install jupyterlab -U
ENV PATH="${PATH}:$HOME/.local/bin/"



#### Add custom jupyterlab modules/packages here

# IRkernel for R integration in jupyter lab (temp switch to root for R install)
USER root
RUN Rscript -e 'install.packages( "IRkernel")'
RUN Rscript -e 'IRkernel::installspec(user = FALSE)'
USER $NB_UID

# Add rpy2 for running R code from python session (and give location to R shared libraries)
RUN pip3 install rpy2
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib/R/lib/



#### Make jupyterlab usable by all (nothing safety-critical in this folder)
# It prevents the need of a script to switch users with UIG/GID declared in 
# environment variables at runtime, which is very complex and has same 
# consequences on security anyway (i.e. no security !).
RUN chmod -R 777 $HOME/.



# Declare web port used by jupyter lab
EXPOSE 8888

# Set default command
CMD $HOME/.local/bin/jupyter-lab --port=8888 --no-browser --ip=0.0.0.0 --LabApp.token="${TOKEN}"

