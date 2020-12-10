FROM continuumio/miniconda3:4.8.2

LABEL description="Autoseq Snakemake workflow"

COPY env/conda_base.yml /

RUN conda env create -f /conda_base.yml && conda clean -a

ENV PATH /opt/conda/envs/autoseq-base/bin:$PATH