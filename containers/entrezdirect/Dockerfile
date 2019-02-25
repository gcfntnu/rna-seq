FROM continuumio/miniconda:4.5.12


ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV MPLBACKEND TkAgg

RUN apt-get update --fix-missing && \
    apt-get install -y \
    wget \
    curl \
    bzip2 \
    ca-certificates \
    build-essential \
    git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN conda update conda
RUN conda config --add channels defaults
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda update --all

RUN conda install --quiet --yes entrez-direct  perl-net-ssleay 

RUN conda clean -tipy
RUN conda build purge-all