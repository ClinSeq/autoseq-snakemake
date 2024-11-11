FROM continuumio/miniconda3:4.8.3

LABEL description="Autoseq Snakemake workflow"

COPY env/* /env/

ENV BASH_ENV ~/.bashrc
SHELL ["/bin/bash", "-c"]

RUN conda init bash
# 1. base env for autoseq-snakemake to run
RUN conda env update -f /env/base.yml && conda clean -a
# 2. conda env for python 2.7 and gatk3
RUN conda env create -f /env/gatk_3.yml && conda clean -a
# 3. conda env for somaticseq
RUN conda env create -f /env/somaticseqenv.yml && conda clean -a

RUN conda env create -f /env/ensemblvep.yml && conda clean -a

RUN conda env create -f /env/purecn-env.yml && conda clean -a && \
    conda activate purecn-env && \
    conda install boost && \
    conda install -c bioconda bioconductor-rhdf5lib && \
    R --quiet -e 'install.packages("BiocManager", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'BiocManager::install("PureCN")' \
    -e 'install.packages("optparse", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'install.packages("futile.logger", repos = "http://ftp.acc.umu.se/mirror/CRAN/")' \
    -e 'install.packages("rjson", repos = "http://ftp.acc.umu.se/mirror/CRAN/")'

RUN ln -s /opt/conda/lib/libreadline.so.7 /opt/conda/lib/libreadline.so.6 && \
    ln -s /opt/conda/lib/libncurses.so.6 /opt/conda/lib/libncurses.so.5

RUN apt-get --allow-releaseinfo-change update && \
    apt-get -y install build-essential --fix-missing && \
    apt-get -y install zlib1g-dev && \
    apt-get -y install libncurses-dev && \
    apt-get install -y gcc --fix-missing && \
    conda env create -f /env/svcallerenv.yml && \
    conda activate svcallerenv && \
    pip install git+https://github.com/tomwhi/svcaller.git 

COPY tools/msings /tools/msings

RUN apt-get update && \
    conda activate gatk_3 && \
    conda install virtualenv && \
    cd /tools/msings/ && \
    bash dev/bootstrap.sh 

COPY . /autoseq-snakemake/
COPY tools/strelka /tools/strelka
COPY tools/somaticseq /tools/somaticseq
COPY tools/vcf2maf /tools/vcf2maf
COPY tools/VarDictJava /tools/VarDictJava

RUN pip install /autoseq-snakemake/

RUN ln -s /opt/conda/envs/gatk_3/lib/libcrypto.so.1.1 /opt/conda/envs/gatk_3/lib/libcrypto.so.1.0.0 && \
    ln -s /opt/conda/lib/libcrypto.so.1.1 /opt/conda/lib/libcrypto.so.1.0.0

ENV MSINGSENV=/tools/msings/msings-env
ENV GRISS_JAR=/autoseq-snakemake/pipeline/scripts/gridss-2.10.2-gridss-jar-with-dependencies.jar
ENV PATH=$PATH:/autoseq-snakemake/pipeline/scripts/:/tools/strelka/bin:/tools/somaticseq/somaticseq:/tools/vcf2maf:/tools/VarDictJava/build/install/VarDict/bin
CMD [ "/bin/bash", "-c" ]