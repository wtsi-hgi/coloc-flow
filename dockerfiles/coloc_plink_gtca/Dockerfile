FROM genomicdatascience/coloc-susie
SHELL ["conda", "run", "-n", "coloc-susie", "/bin/bash", "-c"]
RUN conda install -c bioconda bioconductor-biomart
RUN conda install -c anaconda libxml2 
SHELL ["/bin/sh", "-c"]
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install --yes -y --no-install-recommends dirmngr apt-transport-https ca-certificates software-properties-common gnupg2 \
    build-essential curl git \ 
    g++ gcc gfortran make autoconf automake libtool \
    zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
    libblas-dev libx11-dev libxml2 libxml2-dev \
    libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*

RUN /bin/anaconda/envs/coloc-susie/bin/R -e 'install.packages("remotes", repos="https://cloud.r-project.org")' 
RUN /bin/anaconda/envs/coloc-susie/bin/R -e 'remotes::install_github("mrcieu/ieugwasr")' 
RUN /bin/anaconda/envs/coloc-susie/bin/R -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")'


RUN apt-get install libcurl4-openssl-dev
RUN apt-get install libxml2
RUN which conda
RUN R -e 'BiocManager::install("XML")'
RUN R -e 'BiocManager::install("biomaRt")'

# RUN yes | yum -y install libxml2 libxml2-devel
# RUN wget https://github.com/jianyangqt/gcta/releases/download/v1.94.1beta/gcta_v1.94.1Beta_linux_kernel_5_x86_64_static | ls
# RUN unzip gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip
# RUN https://github.com/jianyangqt/gcta/releases/download/v1.94.1beta/gcta_v1.94.1Beta_linux_kernel_5_x86_64_static
# RUN chmod a+x gcta_v1.94.1Beta_linux_kernel_5_x86_64_static
COPY gcta_v1.94.1Beta_linux_kernel_5_x86_64_static ${HOME}
# LIBXML_INCDIR
# LIBXML_LIBDIR
RUN /bin/anaconda/envs/coloc-susie/bin/R -e 'requireNamespace("BiocManager"); BiocManager::install("biomaRt")'
SHELL ["conda", "run", "-n", "root", "/bin/bash", "-c"]
RUN conda install -c conda-forge r-r.utils
RUN conda install -c bioconda bioconductor-biomart
RUN conda install -c bioconda/label/gcc7 bioconductor-biomart
RUN conda install -c conda-forge r-r.utils
# RUN conda install -c bioconda r-biomartr
# RUN conda install -c conda-forge r-biocmanager
# RUN /bin/anaconda/envs/coloc-susie/bin/R -e 'requireNamespace("BiocManager", quietly = FALSE)'
# RUN /bin/anaconda/envs/coloc-susie/bin/R -e 'BiocManager::install("biomaRt"), repos="https://cloud.r-project.org")'
# RUN conda install -c bioconda r-biomartr
ENV LIBXML_INCDIR=/usr/bin/xmllint
RUN which xmllint
SHELL ["/bin/sh", "-c"]

RUN R -e 'install.packages("XML", repos="https://cloud.r-project.org")' 
SHELL ["conda", "run", "-n", "coloc-susie", "/bin/bash", "-c"]

# RUN which R
