FROM ubuntu:18.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV STAR_VERSION 2.6.1d
ENV SAMTOOLS_VERSION 1.9

RUN set -ex

RUN apt-get update --fix-missing && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y \
  wget \
  ca-certificates \
  git \
  unzip \
  bzip2 \
  build-essential \
  autoconf \
  gcc \
  g++ \
  make \
  wget \
  flex \
  bison \
  autoconf \
  cmake \
  libboost-all-dev \
  libjemalloc-dev \
  zlib1g-dev \
  libbz2-dev \
  libcereal-dev \
  liblzma-dev \
  libncurses5 \
  libncurses5-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  curl && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/COMBINE-lab/salmon.git --branch develop  && \
  mkdir -p salmon/build && \
  mkdir -p /opt/salmon/ && \
  cd salmon/build && \
  cmake -DCMAKE_INSTALL_PREFIX=/opt/salmon .. && \
  make && make install && \
  cd && \
  rm -rf salmon/
ENV PATH /opt/salmon/bin:$PATH
ENV LD_LIBRARY_PATH /opt/salmon/lib:$LD_LIBRARY_PATH

RUN \
  mkdir -p /opt/gffread && \
  cd /opt/gffread && \
  git clone https://github.com/gpertea/gclib && \
  git clone https://github.com/gpertea/gffread && \
  cd gffread && \
  make release
ENV PATH /opt/gffread/gffread:$PATH

RUN git clone https://github.com/lh3/bioawk.git && \
  cd bioawk && \
  make && \
  cp bioawk /usr/local/bin/ && \
  cd .. && rm -rf bioawk

RUN cd /opt && \
  wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz && \
  tar xzf ${STAR_VERSION}.tar.gz && \
  cd STAR-${STAR_VERSION}/source && \
  make STARstatic && \
  mkdir -p /opt/STAR/bin && \
  cp STAR /opt/STAR/bin && \
  cd /opt && \
  'rm' -rf STAR-${STAR_VERSION}

ENV PATH /opt/STAR/bin:${PATH}

RUN curl -L https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2 > /tmp/samtools-$SAMTOOLS_VERSION.tar.bz2 && \
  cd /tmp && tar xfjv samtools-$SAMTOOLS_VERSION.tar.bz2 && rm samtools-$SAMTOOLS_VERSION.tar.bz2 && \
  cd samtools-$SAMTOOLS_VERSION && ./configure && make && make prefix=/opt install && \
  cd /tmp && rm -rf samtools-$SAMTOOLS_VERSION

ENV PATH /opt/bin:$PATH