# Quick start

Installation
-------------

###### System Recommendation:

|||
|----------------------|----------------|
| **Operating System** |  Ubuntu 22.04  |
| **RAM**              |  minimum of 64 GB for small dataset. 128G or more RAM is preferred.|
| **Disk Space**       |  500G is required for complete installation and supported files. But disk space for data is additional.|
| **CPU cores**        |  minimum of 8 CPU cores is required. |

<br>

**Note:** These are the minimum system requirements and may vary with pipeline version and size of input files.
### Software Requirements:

**Required Dependencies:**

```
* python =3.8.12
* Singularity > 3.0 
* conda (package manager)
```

**Optional Dependencies:**

```
* slurm
```

#### Installing Singularity

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

#### Installing Conda

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


### Downloading Reference Genome:

```
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
ar_regions.bed
BrcaExchangeClinvar_15Jan2019_v26_hg19.vcf.gz
human_g1k_v37_decoy.fasta

```

### Autoseq Pipeline Installation

Download and install autoseq-snakemake and the requirements using pip.

```sh
conda create --name autoseq_doc python=3.8.12
conda activate autoseq_doc

git clone https://github.com/Clinseq/autoseq-snakemake.git 
pip install -e autoseq-snakemake/
cd autoseq-snakemake

# To install all required tools in current environment, run the following command.
conda env update --file env/base.yml

# Setting Environmental Variable
GRIDSS_JAR=/path/to/autoseq-snakemake/pipeline/scripts/gridss-2.10.2-gridss-jar-with-dependencies.jar
GRIDSS_SCRIPT = /path/to/tools/gridss-2.13.2/
SAGE_JAR = /path/to/tools/sage_v3.2.3.jar
```

Pull singularity containers 

```sh
singularity pull --arch amd64 library://imsarath/default/autoseq-smk:latest 
singularity pull library://imsarath/default/gridss:latest
# There are additional singularity images present in Ravenclaw server. They are autoseq-ensemblvep.sif, autoseq-franken.sif, autoseq-gatk3.sif, autoseq-jumble.sif, autoseq-purecn.sif, autoseq-somaticseq.sif, autoseq-svcaller.sif.
# Make sure you have copied all these images before launching the pipeline.
```

Once you have successfully installed the pipeline, you can check if the installation was successful with the following command.

```sh
$ autoseq --help
     _         _       ____             
    / \  _   _| |_ ___/ ___|  ___  __ _ 
   / _ \| | | | __/ _ \___ \ / _ \/ _` |
  / ___ \ |_| | || (_) |__) |  __/ (_| |
 /_/   \_\__,_|\__\___/____/ \___|\__, |
                                     |_| 🐍
                         version: 0.1.0


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


Launching autoseq pipeline
--------------------------

To launch autoseq pipeline, you need to specify different config files for autoseq-genome and input samples. You also need to specify the path to input and output directories, path to singularity containers, number of cores and profile options. 


**Input config file**

The input config json file has to be in the following format.
```
{
    "sdid": "P-NA12877",
    "T": "",
    "N": ["LB-P-NA12877-N-03098850-TD1-C31"],
    "CFDNA": ["LB-P-NA12877-CFDNA-03098850-TD1-C31"]
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
