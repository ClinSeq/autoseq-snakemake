# 1 Quick start

## 1.1 Requirements

**System Requirements:**

The system requirements for NGS analysis vary depending on the project's scale, the type of analysis (e.g., whole-genome sequencing, RNA-seq, or targeted sequencing), and the specific tools employed. In general, the following features are highly recommended as a baseline for system specifications. 

|||
|----------------------|----------------|
| **Operating System** |  Any Linux Distrubution (Tested on Ubuntu 22.04)  |
| **RAM**              |  minimum of 64 GB for small dataset. 128G or more RAM is preferred.|
| **Disk Space**       |  500G is required for complete installation of this pipeline and supported files. But disk space for data is additional.|
| **CPU cores**        |  minimum of 8 CPU cores is required. |

<br>
**Note:** These are the minimum system requirements, and they may vary depending on the pipeline version and the size of the input files.

**Software Requirements:**

To run the Autoseq pipeline, several dependencies are needed, including Python, Conda, and optionally, Singularity. We highly recommend using Singularity containers, as they help avoid most dependency conflicts that can occur when working within a conda environment.   

**Dependencies:**

```
* python =3.8.12
* Singularity > 3.0 
* conda
* snakemake==6.2.1  
```

## 1.2 Installation

You can visit Singularity's [installation page](https://sylabs.io/guides/3.0/user-guide/installation.html) to download and install singularity which is compatible with your operating system.

To install conda, please refer to the [Conda installation page](https://docs.anaconda.com/miniconda/#quick-command-line-install)


**Autoseq Pipeline Installation**

Once the required dependencies are installed, you can proceed to install Autoseq-Snakemake in a separate Conda environment using the following commands. It is recommended to use Python version 3.8.12 to create the Conda environment, as other versions may lead to dependency issues. 

```sh
conda create --name autoseq python=3.8.12
conda activate autoseq

git clone https://github.com/Clinseq/autoseq-snakemake.git 
pip install -e autoseq-snakemake/

cd autoseq-snakemake/
conda env create -f env/base.yml
```

These commands will create a Conda environment with Python 3.8.12 and install all required tools as specified in the setup.py file. Once the installation is complete, you can verify its success by running the following command: `autoseq --help`.

```sh
$ autoseq --help
     _         _       ____             
    / \  _   _| |_ ___/ ___|  ___  __ _ 
   / _ \| | | | __/ _ \___ \ / _ \/ _` |
  / ___ \ |_| | || (_) |__) |  __/ (_| |
 /_/   \_\__,_|\__\___/____/ \___|\__, |
                                     |_| 🐍
                         version: 3.4.0


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

#### Downloading the Reference Genome:

Currently, we use the human reference genome GRCh37 for all analyses. The associated files (such as BED, GTF, interval lists, etc.) have been customized specifically for liquid biopsy samples. These files are available for download from our S3 cloud storage upon request.  

## 1.3 Launching autoseq pipeline

Once you have successfully installed all the dependencies and the Autoseq pipeline, you can begin launching your samples. Generally, the Autoseq pipeline can be initiated in two different ways: 

* Using singularity containers (highly recommended)
* Using conda environment (not recommended)

There are three steps involved in launching the Autoseq pipeline. The first two steps remain the same, regardless of whether you use Singularity or Conda. Therefore, we will first discuss these two steps before outlining the specific procedures for launching the pipeline with either Singularity or Conda.

**Step 1:** Preparing the Sample List File
First, ensure that the input directory name adheres to the [recommended](barcodes.md) format. Naming the directories correctly is crucial, as Autoseq uses this format to automatically select the appropriate BED, interval_list, or other files. Next, create a sample list file in `/path/to/project_name/sample_list/` according to the [recommended](barcodes.md) format. You can use the following command to accomplish this:

```
find /path/to/input_directory/ -maxdepth 1 -name "PROJECT*$(date '+%Y%m%d')" \
    | xargs -I {} basename {} | sort -V > \
    /path/to/project_name/sample_lists/clinseqBarcodes_`date "+%Y-%m-%d"`.txt 
```

**Step 2:** Creating hte Config File. 
Once you have created the sample list file as described above, use the `autoseq config` command to generate the config file:

```
mkdir -p /path/to/project_name/config/$(date "+%Y-%m-%d")

screen -S autoseq_run
conda activate autoseq
autoseq config --outdir /path/to/project_name/config/$(date "+%Y-%m-%d") \
    /path/to/project_name/sample_lists/clinseqBarcodes_$(date "+%Y-%m-%d").txt
```

The above command will create a directory with today's date and generate an input config file in `/path/to/project_name/config/YYYY-MM-DD/`. This config file will be used by the Autoseq pipeline in conjunction with the reference config file to launch the samples.

### Launching pipeline with Singularity

It is highly recommended to launch Autoseq using Singularity, as all required dependencies have been containerized, allowing the pipeline to run more consistently compared to using a Conda environment.

To use Singularity, you first need to build the containers for each tool. We have provided a separate script to build all the necessary Singularity images, which you can do using the following commands:

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

These commands will create all the necessary Singularity images, which can then be used to run the pipeline.

**Note:** This is a one-time process. You only need to create the Singularity images with the commands above when launching the pipeline for the first time.

**Step3:** Launching the Pipeline.
Once you have prepared your config files as described in Step 2, you can run the Autoseq pipeline with the following command. It is highly recommended to launch the pipeline within a screen session so you can monitor the process and address any errors that may arise:
```
screen -r autoseq_run

autoseq launch -r autoseq-genome/autoseq-genome.json \
    --samples /path/to/sample.json \
    --outdir /path/to/autoseq-output/ \
    --libdir /path/to/input_directory/ --use-singularity \
    --singularity /path/to/container_dir \
    --umi --cores 8 --smk-opt " --singularity-args \
    '--bind /path/to/autoseq-snakemake/:/path/to/autoseq-snakemake/'"
```

For more information about each parameter used in `autoseq launch`, please visit [pipeline parameters](quick_start.md/#pipeline-parameters)

The above command will initiate the Autoseq pipeline using Singularity images. If the pipeline completes successfully, you will see an `analysis_finished` file in the specified output directory (`/path/to/autoseq-output/sdid/*/`). If the pipeline does not complete successfully, the `analysis_finished` file will not be generated, and you should check the analysis log file to troubleshoot the error.

### Launching pipeline with conda:

If you prefer to launch the Autoseq pipeline using a Conda environment (which is not recommended), you will first need to create all the necessary Conda environments with the following commands.

**NOTE:** This is a one-time process. You only need to create the Conda environments with the commands below when launching the pipeline for the first time. You can check if these Conda environments already exist using the command `conda env list`.

```
conda env create -f env/base.yml
conda env create -f env/ensemblvep.yml
conda env create -f env/franken.yml
conda env create -f env/gatk_3.yml
conda env create -f env/liqbiocna-env.yml
conda env create -f env/purecn-env.yml
conda env create -f env/somaticseqenv.yml
conda env create -f env/svcallerenv.yml
```

These commands will install all the required tools specified in `/autoseq-snakemake/env/*.yml` files within each Conda environment. The pipeline will then utilize these newly created Conda environments to process your samples.

**Step3: Launching the Pipeline.**

Once you have prepared your config files as mentioned in Step 2, you can run the Autoseq pipeline with the following command. It is advisable to execute this command within a screen session so that your job continues to run even if you accidentally close the terminal, allowing you to monitor the progress of your analysis:
```
screen -r autoseq_run
conda activate autoseq

autoseq launch -r autoseq-genome/autoseq-genome.json \
        --samples /path/to/sample.json \
        --outdir /path/to/autoseq-output/ \
        --libdir /path/to/INBOX/ --umi --cores 8
```

To learn more about the description of each of these parameters, visit [pipeline parameters](quick_start.md#pipeline-parameters)

### Error Tracing:

If the pipeline fails for any reason, you can trace the error by reopening the screen session using the command `screen -r autoseq_run` or by checking the analysis log. The analysis log is typically located at `/path/to/autoseq-output/sdid/*/.snakemake/log/date_time.snakemake.log` After opening the log file, identify the rule where the error occurred and then check the corresponding log file, which can be found in`/path/to/autoseq-output/sdid/*/logs/`. If you are unable to resolve the error, please feel free to contact our bioinformatics team for assistance.

### Relaunching Sample:

Once you have identified and resolved the issue, you can relaunch the sample using the Snakemake option `--smk-opt "--rerun-incomplete"`. This option allows Snakemake to restart the analysis from the point where the error occurred, saving time. The following commands will help you relaunch the failed sample.

**If you are using conda to launch samples**
```
screen -r autoseq_run

autoseq launch -r autoseq-genome/autoseq-genome.json \
        --samples /path/to/sample.json \
        --outdir /path/to/autoseq-output/ \
        --libdir /path/to/INBOX/ --umi --cores 8 \
        --smk-opt " --rerun-incomplete "
```

**If you are using singularity to launch samples**
```
screen -r autoseq_run

autoseq launch -r autoseq-genome/autoseq-genome.json \
    --samples /path/to/sample.json \
    --outdir /path/to/autoseq-output/ \
    --libdir /path/to/INBOX/ --use-singularity \
    --singularity /path/to/container_dir \
    --umi --cores 8 --smk-opt " --rerun-incomplete --singularity-args \
    '--bind /path/to/project_name/:/path/to/project_name/'"
```

For more details on the parameters used in `autoseq launch`, please refer to [pipeline parameters](quick_start.md#pipeline-parameters)

### Results:

Once the pipeline successfully completes, you can view the results in `/path/to/autoseq-output/sdid/*/`. The result directory will be organized as follows. For detailed information on each tool used in the pipeline, visit the [autoseq pipeline](autoseq_pipeline.md) page

```
.
|-- bams
|-- cnv
|-- contamination
|-- IGVnav
|-- logs
|-- msings-PROJECT-SDID-TYPE-SAMPLEID-PREPID-CAPTUREID
|-- multiqc
|-- purecn
|-- qc
|-- svs
|-- variants
|-- analysis_finished
|-- config_PROJECT-SDID-T-SAMPLEID-PREPID-CAPTUREID\
        _PROJECT-SDID-N-SAMPLEID-PREPID-CAPTUREID.yml
|-- PROJECT-SDID-T-SAMPLEID-PREPID-CAPTUREID\
        _PROJECT-SDID-N-SAMPLEID-PREPID-CAPTUREID_jobdb.json
|-- PROJECT-SDID-T-SAMPLEID-PREPID-CAPTUREID\
        _PROJECT-SDID-N-SAMPLEID-PREPID-CAPTUREID-igvnav-input.txt
|-- PROJECT-SDID-N-SAMPLEID-PREPID-CAPTUREID-igvnav-input.txt
```

**Bams:** This directory contains all the generated BAM files, although intermediate BAM files will be removed.

**cnv:** This folder includes all the files generated by copy number variation (CNV) calling tools (JUMBLE).

**svs:** This directory contains the output files from structural variant callers (svcaller and GRIDSS).

**Variants:** This directory houses all the files generated during the variant calling process.

**QC:** This directory includes files from quality control tools such as Picard, FastQC, and Samtools, as well as QC overview plots and liqbio-CNA plots.

**Multiqc:** This folder contains the MultiQC results, which summarize all QC tools' outputs.

**Logs:** This directory stores the log files for each tool used in the pipeline. Log files for the alignment process will be saved directly in this directory, while logs for Skewer, svs, variants, Samtools, Picard, FastQC, contamination, and clustering will be organized into subdirectories.

**Analysis Finished File:** An "analysis_finished" file will be generated once the pipeline completes successfully. It records the date and time when the analysis finished. If an error occurs at any point in the pipeline, this file will not be created.

**Config File:** Based on the input parameters for each sample, a configuration file will be automatically generated by modifying `autoseq-snakemake/config.yml` according to the specific inputs. The `autoseq-snakemake/config.yml` file contains default information, such as additional parameters for each tool and the gnomAD database location. Additional details, such as the locations of the Singularity containers, input/output directories, reference config files, and sample-specific config files, will be appended based on the user's input. All of this information will be stored in a file named `config_PB-P-*-T-*-KHYYYYMMDD-CYYYYMMDD_PB-P-*-N-*-KHYYYYMMDD-CYYYYMMDD.yml`.

**JobDB File:** The file `PB-P-*-T-*-KHYYYYMMDD-CYYYYMMDD_PB-P-*-N-*-KHYYYYMMDD-CYYYYMMDD_jobdb.json` will contain the IDs of each submitted job along with their corresponding log files.

**IGVNavigator Input Files:** The files `PB-P-*-N-*-KH-WG-igvnav-input.txt` and `PB-P-*-T-*-KH-WG-PB-P-*-N-*-KH-WG-igvnav-input.txt` files contains information that can be visualized when using the tool IGVNavigator. IGVNavigator is an add-on tool for IGV Visualizer which can assist in analyzing variants during manual curation.

#### Pipeline Parameters

|Parameter|Description|
|------------------|-----------------------------|
|--loglevel|level of logging|
|--verbose/-v|Print verbose output to the console.|
|--outdir|output directory|
|--ref/-r|json file with reference files to use|
|--samples|json file contains list of samples|
|--outdir|output directory|
|--libdir|directory to search libraries|
|--configfile|configuration file for params|
|--cluster-config|configuration file for different HPC|
|--scratch|path to /tmp/scratch|
|--dryrun/--run|dryrun for testing snakemake workflow|
|--umi|To process the data with UMI- Unique Molecular Identifier|
|--profile|job schedulers eg. SLURM|
|--pipeline|Pipeline to be launched|
|--normal-bam/-n|Normal bam files dir, Applicable only to tumor only pipeline|
|--use-singularity|To use singularity|
|--singularity|Path to singularity image|
|--smk-opt|snakemake option|
|--cores|max number of cores|

