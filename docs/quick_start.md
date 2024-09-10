# 1 Quick start

## 1.1 Requirements

###### System Recommendation:

|||
|----------------------|----------------|
| **Operating System** |  Any Linux Distrubution (Tested on Ubuntu 22.04)  |
| **RAM**              |  minimum of 64 GB for small dataset. 128G or more RAM is preferred.|
| **Disk Space**       |  500G is required for complete installation and supported files. But disk space for data is additional.|
| **CPU cores**        |  minimum of 8 CPU cores is required. |

<br>
**Note:** These are the minimum system requirements and may vary with pipeline version and size of input files.

###### Software Requirements:

**Required Dependencies:**

```
* python =3.8.12
* Singularity > 3.0 
* conda
* snakemake==6.2.1  # This will be installed automatically
* click             # This will be installed automatically
* pyyaml            # This will be installed automatically
* pandas            # This will be installed automatically
* rich              # This will be installed automatically
* loguru            # This will be installed automatically
```

**Optional Dependencies:**

```
* slurm
```

**Note:** It is highly recommended to install slurm as a workload manager if you are planning to run your samples on any HPC.

## 1.2 Installation

#### 1.2.1 Installing Singularity

Singularity is a container platform. It allows us to create and run containers that package up pieces of software in a way that is portable and reproducible. To install singularity, first, we need to download all the prerequisites using the following command.

```
sudo apt-get update
sudo apt-get install -y build-essential libssl-dev uuid-dev libgpgme11-dev \
    squashfs-tools libseccomp-dev wget pkg-config git cryptsetup debootstrap \
    libglib2.0-dev runc
```

The above command will download all the prerequisites to install singularity. Since most of the singularity code is written with a programming language called GO, we need to install it first inorder to run singularity. But first, we need ensure that GO language is not already installed in our system. We can delete any pre-existing GO with the following command.

```
rm -rf /usr/local/go
```

Then, download and install GO with the following command.

```
wget https://go.dev/dl/go1.20.5.linux-amd64.tar.gz
sudo tar --directory=/usr/local -xzvf go1.20.5.linux-amd64.tar.gz
export PATH=/usr/local/go/bin:$PATH
```

Then you can download the latest version of singularity from [singularity download page](https://github.com/sylabs/singularity/releases). Install it with the following set of command.

```
wget https://github.com/sylabs/singularity/releases/download/v3.11.4/singularity-ce-3.11.4.tar.gz  # paste the latest singularity link here
tar -xzvf singularity-ce-3.11.4.tar.gz  # Change the singularity version according to the version you have downloaded.

cd singularity-ce-3.11.4  # Change directory name as per your downloaded version.
./mconfig && sudo make -C builddir && sudo make -C builddir install
```

To check if singularity is installed correctly, type the following command. If singularity is installed correctly, it should print help page of singularity

```
singularity --help
```

To know more about singularity installation, you can visit [Singularity's Installation Page](https://sylabs.io/guides/3.0/user-guide/installation.html)

#### 1.2.2 Installing Conda

Conda offers two distributions: Anaconda and Miniconda. Anaconda is a full-featured installer with numerous data science packages and Anaconda Navigator, while Miniconda is a minimal installer with just the Conda package manager and Python. We will use Miniconda because it consumes less memory and disk space compared to Anaconda.

To install miniconda, download the latest version of conda specific to your system from [this link](https://docs.anaconda.com/miniconda/)

Download and install miniconda with the following set of commands.

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Accept the licencing agreement and select the default directory for installing miniconda when asked for it.
```

To ensure conda is installed successfully, type the following command which should list conda help page.

```
conda --help
```

To know more about conda installation, visit [Conda Home Page](https://docs.anaconda.com/miniconda/#quick-command-line-install)


#### 1.2.3 Autoseq Pipeline Installation

Download and install autoseq-snakemake and the requirements using pip.

```sh
conda create --name autoseq_pipeline python=3.8.12
conda activate autoseq_pipeline

git clone https://github.com/Clinseq/autoseq-snakemake.git 
pip install -e autoseq-snakemake/
cd autoseq-snakemake

# To install all required tools in current environment, run the following command.
conda env update --file env/base.yml
```

```sh
$ autoseq --help
     _         _       ____             
    / \  _   _| |_ ___/ ___|  ___  __ _ 
   / _ \| | | | __/ _ \___ \ / _ \/ _` |
  / ___ \ |_| | || (_) |__) |  __/ (_| |
 /_/   \_\__,_|\__\___/____/ \___|\__, |
                                     |_| 🐍
                         version: 3.3.2


Usage: autoseq [OPTIONS] COMMAND [ARGS]...

  Autoseq - pipeline

  Autoseq consists of a custom-pipeline with additional support modules
  aimed  primarily for the analysis of data from high-throughput sequencing
  of liquid biopsies.

Options:
  --loglevel TEXT  level of logging
  -v, --verbose    Print verbose output to the console.
  --help           Show this message and exit.

Commands:
  config  Create sample json for given clinseq barcodes
  launch  launch the respective pipeline with samples json
  list    List autoseq available pipelines with version
```

##### Downloading Reference Genome:

Our developers have customized the reference genome specifically for liquid biopsy analysis which is available in our S3 cloud. Kindly contact any of our [developers](index.md) if you need access to our reference genome. 

##### Creating docker containers:

Our developers have created containerization script to create docker containers for all the tools used in the pipeline and it is available in our [github](https://github.com/ClinSeq/autoseq-docker) page. You can build the containers with the following set of commands.

```
git clone https://github.com/ClinSeq/autoseq-docker.git
cd autoseq-docker/

sudo systemctl start docker
sudo systemctl enable docker

docker build -t autoseq-base -f autoseq-base.Dockerfile .
docker build -t autoseq-ensemblvep -f autoseq-ensemblvep.Dockerfile .
docker build -t autoseq-franken -f autoseq-franken.Dockerfile .
docker build -t autoseq-gatk3 -f autoseq-gatk3.Dockerfile .
docker build -t autoseq-gridss -f autoseq-gridss.Dockerfile .
docker build -t autoseq-jumble -f autoseq-jumble.Dockerfile .
docker build -t autoseq-purecn -f autoseq-purecn.Dockerfile .
docker build -t autoseq-somaticseq -f autoseq-somaticseq.Dockerfile .
docker build -t autoseq-svcaller -f autoseq-svcaller.Dockerfile .
```
The above set of command will create all the docker containers required for running this tool.

Once you have successfully installed the pipeline, you can check if the installation was successful with the following command.

## 1.3 Launching autoseq pipeline

Autoseq pipeline requires two config files (autoseq-genome and input samples), path to input and output directories and singularity containers, number of cores and profile options. Please make sure that the input directory name is consistent with [barcode](barcodes.md) format.

**Input config file**

The input config json file has to be in the following format.
```
{
    "sdid": "P-*******",
    "T": "",
    "N": ["PB-P-*******-N-********-TD1-C31"],
    "CFDNA": ["PB-P-*******-CFDNA-********-TD1-C31"]
}
```

**Reference config file**

Similary, in reference genome config file, you need to specify the following parameters.
```
{
    "1KG": "variants/1000G_phase1.indels.b37.vcf.gz", 
    "Mills_and_1KG_gold_standard": "variants/Mills_and_1000G_gold_standard.indels.b37.vcf.gz", 
    "ar_regions": "intervals/ar_regions.bed", 
    "brca_exchange": "variants/BrcaExchangeClinvar_15Jan2019_v26_hg19.vcf.gz", 
    "bwaIndex": "bwa/human_g1k_v37_decoy.fasta", 
    "dbSNP": "variants/dbsnp142-germline-only.vcf.gz", 
    "fusion_regions": "intervals/fusion_regions.bed", 
    "genes_bed": "genes/human_grch37_87.bed", 
    "oncokb": "variants/OncoKB_6Mar19_v1.9.txt", 
    "reference_dict": "genome/human_g1k_v37_decoy.dict", 
    "reference_genome": "genome/human_g1k_v37_decoy.fasta", 
    "swegene_common": "variants/swegen_common.vcf.gz",
    "pop_vcf": "variants/swegen_common_pop.vcf.gz", 
    "pondir":"pondir",*
    "no_chr":"genome/human_g1k_v37_decoy_non_chr.bed",
    "gene_track": "gene_track_output.csv",
    "small_design": {
        "probio_biomarkersignature2": {
            "targets-bed": "intervals/targets/probio_biomarkersignature2.bed",
            "targets-bed-gz": "intervals/targets/probio_biomarkersignature2.bed.gz",
            "targets-interval_list": "intervals/targets/probio_biomarkersignature2.interval_list",
            "snvindel": {
                "targets-bed": "intervals/targets/probio_snvindel2.bed",
                "targets-bed-gz": "intervals/targets/probio_snvindel2.bed.gz",
                "targets-interval_list": "intervals/targets/probio_snvindel2.interval_list"
            },
            "baseline": {
                "targets-bed": "intervals/targets/probio_baseline2.bed",
                "targets-bed-gz": "intervals/targets/probio_baseline2.bed.gz",
                "targets-interval_list": "intervals/targets/probio_baseline2.interval_list"
            }
        }        
    },
    "sv_filter": "variants/svfilter_2022-07-11.json",
    "cgcann": "genes/cgcann_filter_v1.txt",*
    "exons_gtf": "genes/Homo_sapiens.GRCh37.87.exons-only.gtf",
    "targets": {  
        "probio_comprehensive3": {
            "blacklist-bed": null, 
            "cnvkit-ref": {
                "KAPA_HYPERPREP": {
                    "CFDNA": "intervals/targets/probio_comprehensive3.KAPA_HYPERPREP.CFDNA.cnn", 
                    "N": "intervals/targets/probio_comprehensive3.KAPA_HYPERPREP.N.cnn", 
                    "T": "intervals/targets/probio_comprehensive3.KAPA_HYPERPREP.T.cnn"
                }
            }, 
            "msings-baseline": "intervals/targets/probio_comprehensive3.msings.baseline", 
            "msings-bed": "intervals/targets/probio_comprehensive3.msings.bed", 
            "msings-msi_intervals": "intervals/targets/probio_comprehensive3.msings.msi_intervals", 
            "msisites": "intervals/targets/probio_comprehensive3.slopped20.msisites.tsv", 
            "purecn_targets": null, 
            "targets-bed-slopped20": "intervals/targets/probio_comprehensive3.slopped20.bed", 
            "targets-bed-slopped20-gz": "intervals/targets/probio_comprehensive3.slopped20.bed.gz", 
            "targets-interval_list": "intervals/targets/probio_comprehensive3.interval_list", 
            "targets-interval_list-slopped20": "intervals/targets/probio_comprehensive3.slopped20.interval_list",
            "jumble-ref": "intervals/targets/comprehensive3_baits_twist.bed.reference.RDS"
        }
    },
    "wgs": {
        "cnvkit-ref": "wgs/cnvkit/wgs_reference.N.cnn",
        "jumble-ref": "intervals/targets/wgs.reference.RDS",
        "targets": {
                "bed": "wgs/targets/",
                "interval_list": "wgs/intervals/"
        },
        "hartwig": {
                "ensembl-dir": "wgs/hartwig/ensembl/",
                "actionable-somatic-panel-bed": "wgs/hartwig/ActionableCodingPanel.somatic.37.bed.gz",
                "known-hotspots-somatic-vcf": "wgs/hartwig/KnownHotspots.somatic.37.vcf.gz",
                "NA12878-highconf-bed": "wgs/hartwig/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz"
        }
    }, 
    "ts_regions": "intervals/ts_regions.bed",
    "vep_dir": "vep"
}
```

Once you have prepared your config files as mentioned above, you can run autoseq pipeline with the following command.
```
autoseq launch -r autoseq-genome/autoseq-genome.json --samples /path/to/sample.json --outdir /path/to/autoseq-output/ --libdir /path/to/INBOX/ --use-singularity --singularity /path/to/container_dir --umi --cores 8 --profile slurm --smk-opt " --singularity-args '--bind /path/to/autoseq-snakemake/:/path/to/autoseq-snakemake/'"
```

## 1.4 Results

Once the pipeline gets completed, the resuts folder will appear in the following structure.

```
.
|-- analysis_finished
|-- config_PB-P-*******-CFDNA-********-KH20240305-C420240306_PB-P-*******-N-********-KH20240305-C420240306.yml
|-- msisensor-PB-P-*******-N-04244257-KH-C4-PB-P-*******-CFDNA-********-KH-C4.tsv
|-- PB-P-*******-CFDNA-********-KH20240305-C420240306_PB-P-*******-N-********-KH20240305-C420240306_jobdb.json
|-- PB-P-*******-CFDNA-********-KH-C4-PB-P-*******-N-********-KH-C4-igvnav-input.txt
|-- PB-P-*******-N-********-KH-C4-igvnav-input.txt
|-- bams
|-- cnv
|-- contamination
|-- IGVnav
|-- logs
|-- msings-PB-P-*******-CFDNA-********-KH-C4
|-- multiqc
|-- purecn
|-- qc
|-- svs
|-- variants
```

**Analysis finished file:**
A file named "analysis_finished" will be generated once the pipeline is finished completely. It contains date and time at which the analysis got completed. If any error occurs during any stage of the pipeline, this file will not be generated.

**Config file:**
Depending on the input parameters for each sample, a configuration file will be generated automatically by modifying `autoseq-snakemake/config.yml` as per the input parameters. The `autoseq-snakemake/config.yml` file contains default information such as additional parameters for each tool, gnomad location etc. On top of this default information additional information such as location of each singularity containers, input and output directory, reference config file, sample config file etc will be added based on the parameters provided by the users. All of these information will be stored in `config_PB-P-*****-T-*****-KH20241026-C20241026_PB-P-*****-N-*****-KH20241026-C20241026.yml` file.

**JobDB file:**
`PB-P-*****-T-*****-KH20241026-C20241026_PB-P-*****-N-*****-KH20241026-C20241026_jobdb.json` file will contain ID of each submitted job and their corresponding log file.

**IGVNavigator input files:**
`PB-P-*****-N-*****-KH-WG-igvnav-input.txt` and `PB-P-*****-T-*****-KH-WG-PB-P-*****-N-*****-KH-WG-igvnav-input.txt` files contains information that can be visualized when using the tool IGVNavigator. IGVNavigator is an add-on tool for IGV Visualizer which can assist in analyzing variants during manual review.

**Bams:**
This directory contains all the BAM files that were generated during analysis; however, intermediate BAM files will be removed.

**cnv:**
This directory contains all the files generated by copy number calling tools (JUMBLE).

**svs:**
This directory contains all the files generated by structural variant callers (svcaller and GRIDSS). 

**Variants**
This directory contains all the files generated during variant calling stage. 

**QC**
This directory contains all the files generated by quality check tools such as picard, fastqc, samtools, etc. It also contains QC overview plot and liqbio-cna plot.

**Multiqc**
This directory contains the results of multiqc tool.

**Logs**
This directory contains all the log file of each tool that were created while running the pipeline. Log file for alignment process will be stored directly in this directory; however, the log files of skewer, svs, variants, samtools, picard, fastqc, contamination, and cluster will be stored in seperate directory.

## 1.5 Workflow Structure
The autoseq-snakemake directory is organized as per the following structure.

```
.
├── docs
├── env
├── pipeline
│   ├── autoseq
│   │   ├── rules
│   │   │   ├── alignment.smk
│   │   │   ├── cnvcalling.smk
│   │   │   ├── germline.smk
│   │   │   ├── __init__.py
│   │   │   ├── pre_processing.smk
│   │   │   ├── qc.smk
│   │   │   ├── somatic.smk
│   │   │   ├── split_targets.smk
│   │   │   ├── structuralvariants.smk
│   │   │   ├── umi_processing.smk
│   │   │   └── vep.smk
│   │   └── Snakefile
│   ├── autoseq-rerun
│   │   ├── rules
│   │   │   ├── cnvcalling.smk
│   │   │   ├── germline.smk
│   │   │   ├── qc.smk
│   │   │   ├── somatic.smk
│   │   │   ├── structuralvariants.smk
│   │   │   └── vep.smk
│   │   └── Snakefile
│   ├── autoseq-sd
│   │   ├── rules
│   │   │   ├── qc.smk
│   │   │   ├── split_targets.smk
│   │   │   ├── structuralvariants.smk
│   │   │   ├── umi_processing.smk
│   │   │   └── variant_calling.smk
│   │   └── Snakefile
│   ├── autoseq-wgs
│   │   ├── rules
│   │   │   ├── alignment.smk
│   │   │   ├── cnvcalling.smk
│   │   │   ├── germline.smk
│   │   │   ├── pre_processing.smk
│   │   │   ├── qc.smk
│   │   │   ├── somatic.smk
│   │   │   ├── structuralvariants.smk
│   │   │   └── vep.smk
│   │   └── Snakefile
└── tests
```

**docs**
This directory contains all the documentations related to running autoseq.

**env**
This directory contains all the conda environments associated with running autoseq

**pipeline**
This directory contains a python code called `cli.py` which acts as a entry point to the pipeline. Depending on the input parameters, this code will run specific pipeline. The source code for all the pipelines will be stored in this directory. Each pipeline is store in a seperate subdirectory which are autoseq, autoseq-rerun, autoseq-sd, autoseq-wgs, tumor_only. Inside each of this directory, there will be a Snakefile and a directory named rule which contains all the associated rules for a specific pipeline. Scripts directory contains all the python, shell, and R codes used in the analysis. Utils directory contains all the basic utilities used in the pipeline.

**tests**
This directory contains sample dataset and dummy reference genome and associated files to test the pipeline. 