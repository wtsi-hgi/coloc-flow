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

RUN R -e 'install.packages(c("BiocManager", "remotes", "R.utils", "dplyr", "optparse"), repos="https://cloud.r-project.org", Ncpus=4)'
# RUN R -e 'remotes::install_github("mrcieu/ieugwasr")'
RUN R -e 'remotes::install_github("chr1swallace/coloc@main",build_vignettes=FALSE)'
RUN R -e 'BiocManager::install(c("biomaRt", "liftOver), Ncpus=4)'

# install plink2
RUN wget -q https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20221024.zip -O /tmp/plink2.zip && \
    unzip /tmp/plink2.zip && mv plink2 /usr/bin/ && rm /tmp/plink2.zip

# install gcta
RUN wget -q https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip -O /tmp/gcta.zip && \
    unzip /tmp/gcta.zip gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 && \
    mv gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 /usr/bin/gcta && rm /tmp/gcta.zip

SHELL ["/bin/sh", "-c"]
