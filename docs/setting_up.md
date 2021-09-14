# Setting up environment for autoseq pipeline

You can set up conda environment or singularity container to run autoseq pipeline. Recommended option would be sigularity container.

## conda envs

autoseq pipeline required multiple tools with specific versions, To avoid conflicts between tools, we have multiple conda envs. Conda env yaml files are stored in `autoseq-snakemake/env/*.yml`.

| Conda envs          |
|---------------------|
| base.yml            |
| ensemblvep.yml      |
| gatk_3.yml          |
| liqbiocna-env.yml   |
| purecn-env.yml      |
| somaticseqenv.yml   |
| svcallerenv.yml     |


Creating conda environments

```sh
# you must have installed conda, to update base env with required tools

$ conda env update -f env/base.yml && conda clean -a

# gatk_3 and python 2.7 conda env 

$ conda env create -f env/gatk_3.yml && conda clean -a

# somaticseq - env

$ conda env create -f env/somaticseqenv.yml && conda clean -a

# ensembl-vep - variant annotation

$ conda env create -f env/ensemblvep.yml && conda clean -a

# conda env for svcallerenv 

$ conda env create -f /env/svcallerenv.yml && \ 
    conda activate svcallerenv && \
    pip install git+https://github.com/tomwhi/svcaller.git

# create purecn-env and install required R packages

$ conda env create -f env/purecn-env.yml && conda clean -a && \
    conda activate purecn-env && \
    conda install boost && \
    conda install -c bioconda bioconductor-rhdf5lib && \
    R --quiet -e 'install.packages("BiocManager", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'BiocManager::install("PureCN")' \
    -e 'install.packages("optparse", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'install.packages("futile.logger", repos = "http://ftp.acc.umu.se/mirror/CRAN/")'

# create liqbiocna-env and install required packages

$ conda env create -f /env/liqbiocna-env.yml && \
    conda activate liqbiocna-env && \
    conda install boost && \
    conda install -c bioconda bioconductor-rhdf5lib && \
    R --quiet -e 'install.packages("BiocManager", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'BiocManager::install("VariantAnnotation")' \
    -e 'install.packages("getopt", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'install.packages("gridExtra", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'install.packages("RJSONIO", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'install.packages("ggplot2", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'install.packages("plotly", repos = "http://ftp.acc.umu.se/mirror/CRAN/")'

```

---

**Note: Conda installation might take few minutes depending on tools**

---

## Other bioinformatic tools 

Tools need to be installed in autoseq pipeline. (Not available in conda or modified specific to autoseq pipeline)

* strelka
* msings
* somaticseq

Installation of msings

```sh
# install  msings v3.6
# for customised script run_msings.sh in pipeline/scripts/ 

$ git clone git@gitlab.com:sheenamt/msings.git

# virutualenv required 

$ conda activate gatk_3 && \
    conda install virtualenv && \
    cd msings/ && \
    bash dev/bootstrap.sh

# git clone somaticseq 

$ git clone git@github.com:bioinform/somaticseq.git
```

Modified strelka available on `/nfs/PROBIO/tools/strelka`. 


## ENV_PATH and ENV_VAR

You need to set environment variables and paths to run autoseq-pipeline with conda environment.

```sh
# need to include in your ~/.bashrc file

$ export MSINGSENV=/path/to/msings/msings-env

$ export GRISS_JAR=/path/to/autoseq-snakemake/pipeline/scripts/gridss-2.10.2-gridss-jar-with-dependencies.jar

$ export PATH=$PATH:/path/to/autoseq-snakemake/pipeline/scripts/:/path/to/tools/strelka/bin:/path/to/tools/somaticseq/somaticseq:
```


## Singularity containers

Snakemake uses singularity as primary container system. We created dockerfile to docker container and convert it into singularity. 

```sh
# building docker container

$ cd autoseq-snakemake

$ docker build -t autoseq-smk .

# converting into singularity

$ sudo singularity build container/autoseq-smk.sif docker-daemon://autoseq-smk:latest

# gridss container 

$ sudo singularity build container/gridss.sif docker-daemon://gridss/gridss:latest

```





