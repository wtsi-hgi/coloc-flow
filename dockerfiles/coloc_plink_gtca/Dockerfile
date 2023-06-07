FROM ubuntu:jammy
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install --yes -y --no-install-recommends \
    wget unzip \
    dirmngr apt-transport-https ca-certificates software-properties-common gnupg2 \
    build-essential curl git \
    g++ gcc gfortran make autoconf automake libtool \
    zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
    libblas-dev libx11-dev libxml2 libxml2-dev \
    libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*

ENV NB_USER container_user
ENV HOME /home/${NB_USER}
RUN mkdir ${HOME}

# install R and packages
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/'

RUN apt-get install --yes --no-install-recommends r-base

RUN apt-get install -y libpng-dev libssl-dev liblapack-dev libopenblas-dev

RUN R -e 'install.packages(c("BiocManager", "remotes", "R.utils", "dplyr", "optparse", "tidyr"), repos="https://cloud.r-project.org", Ncpus=4)'
# RUN R -e 'remotes::install_github("mrcieu/ieugwasr")'
RUN R -e 'remotes::install_github("chr1swallace/coloc@main",build_vignettes=FALSE)'
RUN R -e 'BiocManager::install("biomaRt")'
RUN R -e 'install.packages("R.utils", repos="https://cloud.r-project.org")'
RUN which conda
# SHELL ["conda", "run", "-n", "container_user", "/bin/bash", "-c"]
SHELL ["conda", "run", "-n", "root", "/bin/bash", "-c"]
RUN conda install -c bioconda plink
RUN conda install -c bioconda gcta
RUN conda install -c conda-forge r-r.utils
RUN conda install -c conda-forge unzip 
RUN which gcta64
RUN which unzip
RUN which gzip
SHELL ["/bin/sh", "-c"]
COPY gcta_v1.94.1Beta_linux_kernel_5_x86_64_static ${HOME}
RUN ln -s ${HOME}/gcta_v1.94.1Beta_linux_kernel_5_x86_64_static /usr/bin/gcta
RUN which gcta
