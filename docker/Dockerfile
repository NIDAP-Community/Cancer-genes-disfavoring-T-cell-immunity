## Emacs, make this -*- mode: sh; -*-
FROM rocker/r-ubuntu:18.04
#FROM debian:testing

LABEL org.label-schema.license="GPL-2.0" \
      org.label-schema.vcs-url="https://github.com/rocker-org/r-base" \
      org.label-schema.vendor="Rocker Project" \
      maintainer="Dirk Eddelbuettel <edd@debian.org>"

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		fonts-texgyre libx11-dev\
        vim \
	&& rm -rf /var/lib/apt/lists/*

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		software-properties-common \
                dirmngr \
                ed \
		less \
		locales \
		vim-tiny \
		wget curl libcurl4-openssl-dev \
		ca-certificates \
        && add-apt-repository --enable-source --yes "ppa:marutter/rrutter4.0" \
        && add-apt-repository --enable-source --yes "ppa:c2d4u.team/c2d4u4.0+"

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libx11-6 \
    libxss1 \
    libxt6 \
    libxext6 \
    libsm6 \
    libice6 \
    r-cran-png \
    xdg-utils libxt-dev xorg-dev libcairo2 libcairo2-dev libpango1.0-dev firefox \
  && rm -rf /var/lib/apt/lists/*


RUN apt-get update && apt -y install zlib1g-dev
RUN apt -y install libcurl4-openssl-dev
RUN apt -y install libxml2-dev
RUN apt -y install libssl-dev
RUN apt -y install libpng-dev
RUN apt -y install libhdf5-dev

RUN pwd
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-py37_4.8.2-Linux-x86_64.sh  -b \
    && rm -f Miniconda3-py37_4.8.2-Linux-x86_64.sh

RUN conda --version

RUN conda config --append channels  anaconda
RUN conda config --append channels  bioconda
RUN conda config --append channels conda-forge
RUN conda config --append channels bu_cnio
RUN conda config --append channels  esgf

RUN conda install -c conda-forge  mamba


#Download l2p-0.0_1
RUN wget https://github.com/CCBR/l2p/raw/master/r-l2p-0.0_1-r35_0.tar.bz2 -O /tmp/r-l2p-0.0_1-r35_0.tar.bz2 
# Install in a conda enviroment
RUN mamba install /tmp/r-l2p-0.0_1-r35_0.tar.bz2

COPY conda-packages.txt .

RUN PKGS=`cat conda-direct-packages.txt |grep -P -v "r-l2p|^vector-|^r-foundr|spark|snappy|python-file-paths|foundry|networkx|decorator|umap-learn" | tr " " "="  | tr "\n" " "`; mamba install $PKGS


CMD ["bash"]
