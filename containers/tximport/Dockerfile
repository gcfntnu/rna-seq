FROM bioconductor/release_base2:R3.5.2_Bioc3.8 
# use C.UTF-8 locale 
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# R pre-requisites ++
RUN apt-get update --fix-missing && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    zip \
    unzip \
    wget \
    ca-certificates \
    fonts-dejavu \
    build-essential \
    gfortran \
    libxml2 \
    libudunits2-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libcairo2-dev \
    libssl-dev \
    libbz2-dev \
    liblzma-dev \
    git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN R --slave -e "Sys.setenv(TAR = '/bin/tar')"
RUN R --slave -e "options(unzip = '/usr/bin/unzip')"

# RUN R --slave -e "devtools::install_github('user/example')"

RUN R --slave -e "install.packages( c('data.table', 'argparse'), repos='https://cloud.r-project.org')"

RUN R --slave -e "BiocManager::install('tximport')"