# 1 Quick start

## 1.1 Requirements

**System Requirements:**

Next Generation Sequencing Analysis, produce huge volume of data. Hence, to begin NGS data analysis, it is crutial for any research lab to have a highly powerful compute environment. The system requirement varies depending on the scale of the project, type of analysis (e.g., whole-genome sequencing, RNA-seq, or targeted sequencing), and specific set of tools used for each analysis. In general, it is highly recommended to have a system with the following features. 

|||
|----------------------|----------------|
| **Operating System** |  Any Linux Distrubution (Tested on Ubuntu 22.04)  |
| **RAM**              |  minimum of 64 GB for small dataset. 128G or more RAM is preferred.|
| **Disk Space**       |  500G is required for complete installation and supported files. But disk space for data is additional.|
| **CPU cores**        |  minimum of 8 CPU cores is required. |

<br>
**Note:** These are the minimum system requirements and may vary with pipeline version and size of input files.

**Software Requirements:**

To run autoseq pipeline, we need some dependencies which are python, conda, and singularity (optional). We highly recommend using singularity containers since it avoids most of the dependency conflicts which can arise if you are using conda environment.   

**Dependencies:**

```
* python =3.8.12
* Singularity > 3.0 
* conda
* snakemake==6.2.1  
```

## 1.2 Installation

If you prefer to run this pipeline using singularity (which is highly recommended), you can visit their [installation page](https://sylabs.io/guides/3.0/user-guide/installation.html) to download and install singularity which is suitable for your operating system.

To install conda, you can visit [this page](https://docs.anaconda.com/miniconda/#quick-command-line-install)


**Autoseq Pipeline Installation**

Once you have successfully installed the above two dependencies, you can download and install autoseq-snakemake in a seperate conda environment with the following command. It is recommeneded to use python version 3.8.12 to create the conda environment, since other python version may cause dependency conflicts. 

```sh
conda create --name autoseq python=3.8.12
conda activate autoseq

git clone https://github.com/Clinseq/autoseq-snakemake.git 
pip install -e autoseq-snakemake/
```

The above command will create a conda environment with python 3.8.12 and install all the required tools mentioned in setup.py. Once you have successfully installed the pipeline, you can check if the installation was successful by typing the command `autoseq --help`.

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

#### Downloading Reference Genome:

Currently, we are using human reference genome GRCh37 for all our analysis. And, we have customized all the associated files (such as bed, gtf, interval_list, etc) specifically for liquid biopsy samples. It is available in our S3 cloud for download upon request.  

## 1.3 Launching autoseq pipeline

Once you have successfully installed all the dependencies and autoseq pipeline, you can start launching your samples. In general, autoseq pipeline can be launched in two different ways. 

* Using singularity containers (highly recommended)
* Using conda environment (not recommended)

There are 3 steps in launching autoseq pipeline among which the first two steps are same irrespective of you use singularity or conda to launch pipeline. Hence, let's discuss about those 2 steps first and then we'll discuss about the steps to run the pipeline specifically for singularity or conda.

**Step 1:** Preparing sample list file.
First, make sure that the input directory name is consistent with the [recommended](barcodes.md) format. Naming the files in this format is crutial, because, autoseq uses this format to select the appropriate bed (or) interval_list (or) other files automatically. Then, you need to create a sample list file inside `/path/to/project_name/sample_list/` as per the [recommended](barcodes.md) format. You can use the following command to create the same.

```
find /path/to/input_directory/ -maxdepth 1 -name "PROJECT*$(date '+%Y%m%d')" \
    | xargs -I {} basename {} | sort -V > \
    /path/to/project_name/sample_lists/clinseqBarcodes_`date "+%Y-%m-%d"`.txt 
```

**Step 2:** Creating config file. 
Once you have created the sample list file as mentioned above, you can use `autoseq config` command to create the config file.

```
mkdir -p /path/to/config/$(date "+%Y-%m-%d")

screen -S autoseq_run
conda activate autoseq
autoseq config --outdir /path/to/config/$(date "+%Y-%m-%d") \
    /path/to/project_name/sample_lists/clinseqBarcodes_$(date "+%Y-%m-%d").txt
```

The above command will create a directory in today's date, and then create input config file in `/path/to/config/YYYY-MM-DD/`. This config file is then used by autoseq pipeline along with reference config file to launch samples.

### Launching pipeline with Singularity

It is highly recommended to launch autoseq using singularity because all the required dependencies has been containerized, so that the pipeline can run more consistently compared to running the pipeline using conda environment.

To use singularity, first you need to build the containers for each tool. We have written a seperate script to build all the required singularity images, for which, you can use the following code.

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

This should create all the required singularity images which can then be used to run the pipeline.

**Note:** This is a one time process. Hence, you need to create the singularity images with the above command only when you are launching the pipeline for the first time.

**Step3:** Launching Pipeline.
Once you have prepared your config files as mentioned in step 2, you can run autoseq pipeline with the following command. It is highly recommended to launch the pipeline in screen, so that you can monitor the pipeline if incase any error occurs.
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

To know more about each parameters used in `autoseq launch`, please visit [pipeline parameters](quick_start.md/#pipeline-parameters)

The above command will launch the autoseq pipeline using singularity images. If the pipeline gets completed successfully, you will be able to see `analysis_finished` file in the specified output directory (`/path/to/autoseq-output/sdid/*/`). If pipeline is not completed successfully, `analysis_finished` will not be generated, and you need to check the analysis log file to fix the error. 

### Launching pipeline with conda:

If you prefer to launch autoseq pipeline using conda environment (which is not recommended), first, you need to create all the conda environment with the following command.

**NOTE:** This is a one time process. Hence, you need to create the conda environment with the following command only when you are launching the pipeline for the first time. You can also check if these conda environments are already available using the command `conda env list`.

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

These command will install all the required tools mentioned in .yml file inside each conda environment. Now, our pipeline can use these newly created conda environment to run your samples.

**Step3: Launching Pipeline.**

Once you have prepared your config files as mentioned in step 2, you can run autoseq pipeline with the following command. It is highly recommended to run the below command inside screen, so that your job will continue to run even if you accidently close the terminal.
```
screen -r autoseq_run
conda activate autoseq

autoseq launch -r autoseq-genome/autoseq-genome.json \
        --samples /path/to/sample.json \
        --outdir /path/to/autoseq-output/ \
        --libdir /path/to/INBOX/ --umi --cores 8
```

To know about the description of each of these parameters, visit [pipeline parameters](quick_start.md#pipeline-parameters)

### Tracing Error:

If the pipeline failes for any reason, you can check the error either by opening screen using the command `screen -r autoseq_run` or you can check the analysis log. The analysis log will usually be saved under `/path/to/autoseq-output/sdid/*/.snakemake/log/date_time.snakemake.log` file. Once you open this file, you can check for the rule in which the error has occured and open the corresponding log file. The log file for each tool will be present in `/path/to/autoseq-output/sdid/*/logs/`. Check if you can resolve the error, else, please feel free to reachout to our bioinformatics team.

### Relaunching Sample:

Once you have identified and resolved the error, you can relaunch the sample with snakemake option `--smk-opt "--rerun-incomplete"`. This option will enable snakemake to re-start the analysis from the step where the error has occured, so that you can save time. Here is the command to re-launch the failed sample.

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

To know about the description of each of these parameters, visit [pipeline parameters](quick_start.md#pipeline-parameters)

### Results:

Once the pipeline gets completed successfully, you can view the results under `/path/to/autoseq-output/sdid/*/`. The result directory will appear in the following structure. To know more about each tool used in the pipeline, please visit [autoseq pipeline](autoseq_pipeline.md) page

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
|-- config_PROJECT-SDID-T-SAMPLEID-PREPID-CAPTUREID_PROJECT-SDID-N-SAMPLEID-PREPID-CAPTUREID.yml
|-- PROJECT-SDID-T-SAMPLEID-PREPID-CAPTUREID_PROJECT-SDID-N-SAMPLEID-PREPID-CAPTUREID_jobdb.json
|-- PROJECT-SDID-T-SAMPLEID-PREPID-CAPTUREID_PROJECT-SDID-N-SAMPLEID-PREPID-CAPTUREID-igvnav-input.txt
|-- PROJECT-SDID-N-SAMPLEID-PREPID-CAPTUREID-igvnav-input.txt
```

**Bams:** This directory contains all the BAM files that were generated during analysis; however, intermediate BAM files will be removed.

**cnv:** This directory contains all the files generated by the copy number calling tools (JUMBLE).

**svs:** This directory contains all the files generated by structural variant callers (svcaller and GRIDSS).

**Variants:** This directory contains all the files generated during variant calling stage.

**QC:** This directory contains all the files generated by quality check tools such as picard, fastqc, samtools, etc. It also contains QC overview plot and liqbio-cna plot.

**Multiqc:** This directory contains the results of multiqc tool.

**Logs:** This directory contains all the log file of each tool that were created while running the pipeline. Log file for alignment process will be stored directly in this directory; however, the log files of skewer, svs, variants, samtools, picard, fastqc, contamination, and cluster will be stored in seperate directory.

**Analysis finished file:** A file named "analysis_finished" will be generated once the pipeline is finished completely. It contains date and time at which the analysis got completed. If any error occurs during any stage of the pipeline, this file will not be generated.

**Config file:** Depending on the input parameters for each sample, a configuration file will be generated automatically by modifying autoseq-snakemake/config.yml as per the input parameters. The autoseq-snakemake/config.yml file contains default information such as additional parameters for each tool, gnomad location etc. On top of this default information additional information such as location of each singularity containers, input and output directory, reference config file, sample config file etc will be added based on the parameters provided by the users. All of these information will be stored in config_PB-P-*-T-*-KHYYYYMMDD-CYYYYMMDD_PB-P-*-N-*-KHYYYYMMDD-CYYYYMMDD.yml file.

**JobDB file:** PB-P-*-T-*-KHYYYYMMDD-CYYYYMMDD_PB-P-*-N-*-KHYYYYMMDD-CYYYYMMDD_jobdb.json file will contain ID of each submitted job and their corresponding log file.

**IGVNavigator input files:** PB-P-*-N-*-KH-WG-igvnav-input.txt and PB-P-*-T-*-KH-WG-PB-P-*-N-*-KH-WG-igvnav-input.txt files contains information that can be visualized when using the tool IGVNavigator. IGVNavigator is an add-on tool for IGV Visualizer which can assist in analyzing variants during manual review.

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

