FROM continuumio/miniconda3:4.8.2

LABEL description="Autoseq Snakemake workflow"

COPY env/* /env/

# 1. base env for autoseq-snakemake to run
RUN conda env update -f /env/base.yml && conda clean -a
# 2. conda env for python 2.7 and gatk3
RUN conda env create -f /env/gatk_3.yml && conda clean -a
# 3. conda env for somaticseq
RUN conda env create -f /env/somaticseqenv.yml && conda clean -a

RUN conda env create -f /env/ensemblvep.yml && conda clean -a

RUN export PATH=/opt/conda/bin:$PATH; \
    conda env create -f /env/purecn-env.yml && conda clean -a && \
    conda init bash

ENV BASH_ENV ~/.bashrc
SHELL ["/bin/bash", "-c"]

RUN conda activate purecn-env && \
    R --quiet -e 'install.packages("BiocManager", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'BiocManager::install("PureCN")' \
    -e 'install.packages("optparse", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'install.packages("futile.logger", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' 

RUN conda env create -f /env/liqbiocna-env.yml

RUN conda activate liqbiocna-env && \
    R --quiet -e 'install.packages("BiocManager", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'BiocManager::install("VariantAnnotation")' 

RUN conda env create -f /env/flanken.yml && 

RUN conda activate flanken && \
    R --quiet -e 'install.packages("BiocManager", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'BiocManager::install("VariantAnnotation")'  \
    -e 'install.packages("getopt", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'install.packages("gridExtra", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' 

ENV PATH /opt/conda/envs/autoseq-base/bin:$PATH