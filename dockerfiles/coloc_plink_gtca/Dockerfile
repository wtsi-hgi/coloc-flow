FROM ubuntu
RUN cat /etc/os-release
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install --yes -y --no-install-recommends dirmngr apt-transport-https ca-certificates software-properties-common gnupg2 \
    build-essential curl git \ 
    g++ gcc gfortran make autoconf automake libtool \
    zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
    libblas-dev libx11-dev libxml2 libxml2-dev \
    libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*
ENV NB_USER container_user
ENV HOME /home/${NB_USER}
RUN mkdir ${HOME}
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
RUN apt update && apt-get install wget
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh -O /tmp/miniconda.sh && chmod +x  /tmp/miniconda.sh && bash /tmp/miniconda.sh -b -p ${HOME}/conda
# # RUN bash Miniconda3-latest-Linux-x86_64.sh
ENV PATH ${HOME}/conda/bin:$PATH
COPY environment.yml ${HOME}/tmp2/environment.yml
RUN conda config --set always_yes yes --set changeps1 no && \
    conda env update -v -n root --file ${HOME}/tmp2/environment.yml && \
    conda list --name root && \
    conda clean -atipy && \
    rm ${HOME}/tmp2/environment.yml
WORKDIR ${HOME}/tmp2
RUN which R
RUN R -e 'install.packages("remotes", repos="https://cloud.r-project.org")' 
RUN R -e 'remotes::install_github("mrcieu/ieugwasr")' 
RUN R -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")'
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