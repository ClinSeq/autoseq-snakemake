## Our Server:
At Karolinska Institutet, we use internal High Performance Computing (HPC) under Infrastructure as a Service (IaaS) provided by the Karolinska Institutet's Central IT. Our setup of VMs has one head node and five compute nodes. The head node contains with 12 CPUs, 64 GB of RAM, and 500 GB of SSD storage. The compute node contains 48 CPUs each, 252 GB of RAM, and 250 TB of disk space. 

## Workload Manager:
To efficiently manage these compute nodes and the compute jobs on them, we use a specialized software called a workload manager, which allocates resources for each job, tracks progress, monitors performance, and reports the job status. SLURM is our workload manager of choice. Read more about Slurm on their documentation [page](https://slurm.schedmd.com/documentation.html). In autoseq pipeline, we use Slurm as the workload manager, responsible for resource allocation. When working with HPC, users only need to log in to the head node to submit jobs. Slurm automatically assigns jobs to compute nodes based on the requested resources. Once a compute node completes the analysis, it sends a status report back to the head node. This system ensures that all analyses are managed efficiently by Slurm. If you're working on a different server and Slurm is not installed, please contact your IT support team to facilitate its installation. 

## Data Organization in the HPC:
Organizing project data folder is essential for easy navigation and access. Below is an overview of how data is structured on our HPC.

All projects and the autoseq pipeline are stored on a scalable NFS with all nodes as hosts. The structure is as follows:

```
ls /nfs/
|-- PIPELINE
|-- Project1
|-- Project2
.
.
.
```

### Pipeline Structure:
The autoseq pipeline is located in the `/nfs/PIPELINE/` directory, which is organized as follows:
```
ls /nfs/PIPELINE/
|-- autoseq-genome
|-- autoseq-snakemake
|-- containers
|-- tools
```

* The `autoseq-genome` directory houses the reference genome (GRCh37) and associated files. We currently support only GRCh37 but we are working on upgrading to the GRCh38 reference genome.
* The `autoseq-snakemake` directory contains all pipeline scripts, Conda environments, documentation, and test samples for verifying pipeline functionality.
* The `containers` directory stores all Singularity images used by the pipeline.
* The `tool` directory contains additional tools required for the autoseq pipeline.

Within the `autoseq-snakemake` directory, there are subdirectories for each of the five pipelines we use. Each pipeline has a `rules` directory containing the specific rules, along with a `Snakefile`. The structure for these subdirectories is as follows:

```
|-- docs
|-- env
|-- pipeline
│   │-- autoseq
│   │   │--rules
│   │   │   │-- alignment.smk
│   │   │   │-- cnvcalling.smk
│   │   │   │-- germline.smk
│   │   │   │-- pre_processing.smk
│   │   │   │-- qc.smk
│   │   │   │-- somatic.smk
│   │   │   │-- split_targets.smk
│   │   │   │-- structuralvariants.smk
│   │   │   │-- umi_processing.smk
│   │   │   │-- vep.smk
│   │   │-- Snakefile
│   │-- autoseq-rerun
│   │   │--rules
│   │   │   │-- cnvcalling.smk
│   │   │   │-- germline.smk
│   │   │   │-- qc.smk
│   │   │   │-- somatic.smk
│   │   │   │-- structuralvariants.smk
│   │   │   │-- vep.smk
│   │   │-- Snakefile
│   │-- autoseq-sd
│   │   │--rules
│   │   │   │-- qc.smk
│   │   │   │-- split_targets.smk
│   │   │   │-- structuralvariants.smk
│   │   │   │-- umi_processing.smk
│   │   │   │-- variant_calling.smk
│   │   │-- Snakefile
│   │-- autoseq-wgs
│   │   │--rules
│   │   │   │-- alignment.smk
│   │   │   │-- cnvcalling.smk
│   │   │   │-- germline.smk
│   │   │   │-- pre_processing.smk
│   │   │   │-- qc.smk
│   │   │   │-- somatic.smk
│   │   │   │-- structuralvariants.smk
│   │   │   │-- vep.smk
│   │   │-- Snakefile
│   │-- tumor_only
│   │   │--rules
│   │   │   │-- cnvcalling.smk
│   │   │   │-- pre_processing.smk
│   │   │   │-- qc.smk
│   │   │   │-- structuralvariants.smk
│   │   │   │-- variant_calling.smk
│   │   │   │-- vep.smk
│   │   │-- Snakefile
|-- tests
```

### Project Structure:
We manage different projects, all organized under the `/nfs/` directory. Each project follows a standardized structure, as illustrated below:

```
|-- Project1
│   │-- INBOX
│   │-- autoseq-output
│   │-- config
│   │-- logs
│   │-- sample_lists
```

* The `INBOX` directory contains all sample FASTQ files, and the directory names must follow the pattern specified in the [General Description page](barcodes.md)
* The `autoseq-output` directory stores the analysis results.
* The `config` directory contains configuration files, including tool-specific parameters used in the pipeline.
* The `logs` directory holds the "no hungup" logs for each sample analysis.
* The `sample_lists` directory includes files named `clinseqBarcodes_YYYY-MM-DD.txt`, each listing the samples analyzed on a particular date.

**Note:** When starting a new project or creating a new pipeline on our server, it is crucial to maintain the above directory structure for consistency.

## Virtual Environment

To streamline the analysis process, we have set up a virtual environment on the Ravenclaw server with all necessary dependencies. You can activate this environment using the command `prod_up`.

## Launching Multiple Samples on the Server

When launching multiple samples on the server, it is essential to submit your jobs via Slurm. Slurm will manage resource allocation for each job. Our Ravenclaw server is already equipped with Slurm, and the virtual environment `prod_up` is pre-configured.

If you are using a different server and need to specify memory or time requirements for certain jobs, you will need to modify the cluster configuration file accordingly. The cluster configuration file is located at `/nfs/PIPELINE/autoseq-snakemake/pipeline/scheduler/cluster_config.json`

Below is an example structure of the cluster configuration file:

```
{
    "__default__": {
        "time": "48:00:00",
        "output": "logs/cluster/{rule}-%j.out",
        "error": "logs/cluster/{rule}.err",
        "partition": "core"
    },
    "rule_name_1": {
        "time": "100:00:00"
    },
    "rule_name_2": {
        "time": "200:00:00",
        "mem": "12000"
    },
}
```

In this configuration, the `__default__` settings apply to all rules, while any parameters specified under individual rules (such as `rule_name_1` and `rule_name_2`) will override the default settings. You can pass this configuration file to the autoseq pipeline using the `--cluster-config` parameter.

Once the cluster configuration file is prepared, you are ready to launch your samples. The pipeline launch involves three key steps: 

 * Naming the input directory according to the barcode format
 * Preparing the configuration file
 * Launching the samples

#### Renaming the Input Directory

When running batch samples, it is critical to name the input directory following the format  `PROJECT-SDID-TYPE-SAMPLEID-PREPID-CAPTUREID`. Autoseq relies on this structure, especially the prepID and captureID, to select appropriate reference files for analysis. Ensure that the prepID and captureID include the received date or launch date (in YYYYMMDD format) of the sample as this allows Autoseq to retrieve all samples for processing on that specific date. For more details on this format, please refer to the [General Description](barcodes.md) page. Below is an example of the directory structure:

```
.
|-- INBOX/
|   |-- PROJECT-P-*-N-*-KHYYYYMMDD-CYYYYMMDD/
|   |   |-- *_DNA-B-*-00_S11_L001_R1_001.fastq.gz
|   |   |-- *_DNA-B-*-00_S11_L001_R2_001.fastq.gz
|   |   |-- *_DNA-B-*-00_S11_L002_R1_001.fastq.gz
|   |   |-- *_DNA-B-*-00_S11_L002_R2_001.fastq.gz
|   |-- PROJECT-P-*-T-*-KHYYYYMMDD-CYYYYMMDD/
|   |   |-- *_DNA-T-*-00_S11_L001_R1_001.fastq.gz
|   |   |-- *_DNA-T-*-00_S11_L001_R2_001.fastq.gz
|   |   |-- *_DNA-T-*-00_S11_L002_R1_001.fastq.gz
|   |   |-- *_DNA-T-*-00_S11_L002_R2_001.fastq.gz
|   |-- PROJECT-P-*-N-*-KHYYYYMMDD-CYYYYMMDD/
|   |   |-- *_DNA-B-*-00_S11_L001_R1_001.fastq.gz
|   |   |-- *_DNA-B-*-00_S11_L001_R2_001.fastq.gz
|   |   |-- *_DNA-B-*-00_S11_L002_R1_001.fastq.gz
|   |   |-- *_DNA-B-*-00_S11_L002_R2_001.fastq.gz
|   |-- PROJECT-P-*-CFDNA-*-KHYYYYMMDD-CYYYYMMDD/
|   |   |-- *_DNA-T-*-00_S11_L001_R1_001.fastq.gz
|   |   |-- *_DNA-T-*-00_S11_L001_R2_001.fastq.gz
|   |   |-- *_DNA-T-*-00_S11_L002_R1_001.fastq.gz
|   |   |-- *_DNA-T-*-00_S11_L002_R2_001.fastq.gz
```

#### Preparing the Configuration File

After setting up the input directory as described, the next step is to create the configuration file. This file, in JSON format, includes details such as the SDID, Tumor ID/CFDNA ID, and Normal ID (matching the directory name). Autoseq uses this information to locate specific directories within `/nfs/project_name/INBOX/`. To generate the configuration file, use the `autoseq config` command, which creates the config file in `/nfs/project_name/config/YYYY-MM-DD/` with file name `SDID.json`

Here’s the command sequence for generating the configuration file:

```
screen -S autoseq_run # use 'screen -r autoseq_run' if the screen is already active. 
                      # You can check it using the command "screen -ls"
prod_up
find /nfs/project_name/INBOX/ -maxdepth 1 \
        -name "PROJECT*$(date '+%Y%m%d')" | \
        xargs -I {} basename {} | sort -V > \
        /nfs/project_name/sample_lists/clinseqBarcodes_`date "+%Y-%m-%d"`.txt
mkdir -p /nfs/project_name/config/$(date "+%Y-%m-%d")
autoseq config --outdir /nfs/project_name/config/$(date "+%Y-%m-%d") \
        /nfs/project_name/sample_lists/clinseqBarcodes_$(date "+%Y-%m-%d").txt
# You can use ctrl+a then d to exit the screen.
```

This script performs the following actions:

* Searches all directories in `/nfs/project_name/INBOX/` ending with today’s date.
* Sorts these directory names and writes them into `/nfs/project_name/sample_lists/` as `clinseqBarcodes_today's_date`.
* Creates a new directory within `/nfs/project_name/config/` with today’s date.
* Generates individual JSON configuration files for each sample, structured as shown below:

```
.
|-- config
|   |-- YYYY-MM-DD/
|   |   |-- P-*.json
|   |   |-- P-*.json
```

Each configuration file will have the following format:

```
# P-*.json
{
    "sdid": "P-*",
    "T": ["PROJECT-P-*-T-*-KHYYYYMMDD-CYYYYMMDD"],
    "N": ["PROJECT-P-*-N-*-KHYYYYMMDD-CYYYYMMDD"],
    "CFDNA": ""
}
```

Following these steps ensures accurate and organized setup for batch analysis using the autoseq pipeline.

#### Launching Samples

After configuring all necessary files, you can launch multiple samples in batch using the following shell script.

**Note:**  It is highly recommended to start your analysis within a `screen` session to prevent the process from terminating if the terminal is closed accidentally.

```
screen -r autoseq_run
prod_up                # run this only if production environment is not active.
ref=/nfs/PIPELINE/autoseq-genome/autoseq-genome.json
libdir=/nfs/project_name/INBOX/
outdir=/nfs/project_name/autoseq-output/
cores=8
configs=(`find /nfs/project_name/config/$(date "+%Y-%m-%d") \
        -name "P-*json" | sort -r`)
for config in ${configs[@]}; do
    echo $config
    sdid=`basename $config |cut -f 1 -d "."`;
    echo $sdid
    nohup autoseq launch -r $ref --samples $config \
        --outdir $outdir --libdir $libdir \
        --cluster-config /path/to/cluster_config_specific_server.json \
        --use-singularity --singularity /nfs/PIPELINE/containers/ \
        --scratch /path/to/tmp --umi --cores $cores --profile slurm \
        --smk-opt "--latency-wait 60 --singularity-args \
        '--bind /base-path/:/base-path/'" >> logs/$sdid.nohub.log &
    sleep 250
done
# You can use ctrl+a then d to exit the screen.
```

In this script:

* Paths to the reference file (ref), INBOX directory (libdir), output directory (outdir), and configuration files (for the specific date) are initialized.
* For each sample configuration JSON file, the sdid is extracted, and the `autoseq launch` command is executed to initiate each job on the Slurm cluster.
* `nohup` is used to prevent interruptions, and `sleep 250` ensures a 250-second delay before the next job submission.

To learn more about each parameter used in autoseq launch, refer to [pipeline parameters](quick_start.md/#pipeline-parameters)

After launching jobs, track their status with the following slurm command:

```
squeue -o "%.7i %.4P %a %.60j %.20u %.8T %.10M %.9l %.6D %.6C %.6m %R"
##   [OR]
squeue -o "%.7i %.4P %a %.60j %.20u %.8T %.10M %.9l %.6D %.6C %.6m %R" | grep "username"
```

If job statuses show anything other than `PENDING` or `RUNNING` (such as `DependencyNeverSatisfied`), you may need to troubleshoot manually or contact bioinformatics or IT support. Additionally, you can verify completed samples by checking for the `analysis_finished` file in `/nfs/project_name/autoseq-output/sdid/*/`. Use the following command to list completed samples:

```
ls /nfs/project_name/autoseq-output/*/*/analysis_finished | \
         grep -E "sdid1|sdid2|sdid3|....|sdidn"
```

## Relaunching Failed Samples
If any jobs fail, address the issue and then relaunch the failed samples. First, cancel all Slurm jobs for the failed samples with:

```
jids=($(squeue -o "%j %i" | grep -E "sdid1|sdid2|sdid3|....|sdidN" | \
         cut -f 2 -d " "))
echo ${jids[@]}
for jid in ${jids[@]}; do
  scancel $jid
done
```

Next, confirm that no background jobs are running:

```
ps aux | grep "username"
```

If any background jobs are detected, terminate them with:

```
kill -9 $(ps aux | grep -E "sdid1|sdid2|sdid3|....|sdidN" | \
        grep -v "grep" | awk '{print $2}')
```

In cases of abrupt pipeline termination (server crash or other issues), the output directory for a sample may remain locked. Unlock these directories before relaunching by noting the config file paths for failed jobs and using the following:

```
screen -r autoseq_run
prod_up                # run this only if production environment is not active.
ref=/nfs/PIPELINE/autoseq-genome/autoseq-genome.json
libdir=/nfs/project_name/INBOX/
outdir=/nfs/project_name/autoseq-output/
cores=8

## Provide the config files for failed samples below.
configs=(/nfs/project_name/config/date/sdid1.json \
        /nfs/project_name/config/date/sdid2.json \
        /nfs/project_name/config/date/sdid3.json ... \
        /nfs/project_name/config/date/sdidN.json)  

for config in ${configs[@]}; do
    echo $config
    sdid=`basename $config |cut -f 1 -d "."`;
    echo $sdid
    nohup autoseq launch -r $ref --samples $config \
        --outdir $outdir --libdir $libdir \
        --cluster-config /path/to/cluster_config.json \
        --use-singularity --singularity /nfs/PIPELINE/containers/ \
        --scratch /path/to/tmp --umi --cores $cores --profile slurm \
        --smk-opt "--latency-wait 5 --unlock \
        --singularity-args '--bind /base-path/:/base-path/'"
    sleep 10
done
```

This command will unlock the specified sample directories. You can now relaunch failed samples as follows: 

```
screen -r autoseq_run
prod_up
ref=/path/to/autoseq-genome/autoseq-genome.json
libdir=/path/to/INBOX/
outdir=/path/to/autoseq-output/
cores=8

## Provide the config files for failed samples below.
configs=(/path/to/config/date/sdid1.json \
        /path/to/config/date/sdid2.json  \
        /path/to/config/date/sdid3.json ... \
        /path/to/config/date/sdidN.json)  

for config in ${configs[@]}; do
    echo $config
    sdid=`basename $config |cut -f 1 -d "."`;
    echo $sdid
    nohup autoseq launch -r $ref --samples $config \
        --outdir $outdir --libdir $libdir \
        --cluster-config /path/to/cluster_config.anchorage.json \
        --use-singularity --singularity /nfs/PIPELINE/containers/ \
        --scratch /path/to/tmp --umi --cores $cores --profile slurm \
        --smk-opt "--latency-wait 60 --rerun-incomplete \
        --singularity-args '--bind /base-path/:/base-path/'"
    sleep 250
done
```
