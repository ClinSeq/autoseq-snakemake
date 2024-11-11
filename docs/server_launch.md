## Our Server:
At Karolinska Institutet, we operate a High Performance Computing (HPC) system named "Ravenclaw." The server comprises one head node and five compute nodes. The head node is equipped with 12 CPUs, each with 6 cores, 64 GB of RAM, and 500 GB of SSD storage. Each compute node has 48 CPUs, 252 GB of RAM, and 250 TB of disk space. 

## Workload Manager:
To efficiently manage these compute nodes, we use a specialized software called a workload manager, which allocates resources for each job, tracks progress, monitors performance, and reports the job status. In our autoseq pipeline, we use Slurm as the workload manager, responsible for resource allocation. When working with HPC, users only need to log in to the head node to submit jobs. Slurm automatically assigns jobs to compute nodes based on the requested resources. Once a compute node completes the analysis, it sends a status report back to the head node. This system ensures that all analyses are managed efficiently by Slurm. If you're working on a different server and Slurm is not installed, please contact your IT support team to facilitate its installation. 

## Data Organization in Ravenclaw:
Organizing project data on the Ravenclaw server is essential for easy navigation and access. Below is an overview of how data is structured within the Ravenclaw server.

All projects and the autoseq pipeline are stored on a scalable SSD under the `/nfs/` directory. The structure is as follows:

```
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
|-- autoseq-genome
|-- autoseq-snakemake
|-- containers
|-- tools
```

* The `autoseq-genome` directory houses the reference genome (GRCh37) and associated files.
* The `autoseq-snakemake` directory contains all pipeline scripts, Conda environments, documentation, and test samples for verifying pipeline functionality.
* The `containers` directory stores all Singularity images used by the pipeline.
* The `tool` directory contains additional tools required for the autoseq pipeline.

Within the `autoseq-snakemake` directory, there are subdirectories for each of the five pipelines we use. Each pipeline has a `rules` directory containing the specific rules, along with a `Snakefile`. The structure for these subdirectories is as follows:

```
|-- docs
|-- env
|-- pipeline
│   │-- autoseq
│   │   │--rules
│   │   │   │-- alignment.smk
│   │   │   │-- cnvcalling.smk
│   │   │   │-- germline.smk
│   │   │   │-- pre_processing.smk
│   │   │   │-- qc.smk
│   │   │   │-- somatic.smk
│   │   │   │-- split_targets.smk
│   │   │   │-- structuralvariants.smk
│   │   │   │-- umi_processing.smk
│   │   │   │-- vep.smk
│   │   │-- Snakefile
│   │-- autoseq-rerun
│   │   │--rules
│   │   │   │-- cnvcalling.smk
│   │   │   │-- germline.smk
│   │   │   │-- qc.smk
│   │   │   │-- somatic.smk
│   │   │   │-- structuralvariants.smk
│   │   │   │-- vep.smk
│   │   │-- Snakefile
│   │-- autoseq-sd
│   │   │--rules
│   │   │   │-- qc.smk
│   │   │   │-- split_targets.smk
│   │   │   │-- structuralvariants.smk
│   │   │   │-- umi_processing.smk
│   │   │   │-- variant_calling.smk
│   │   │-- Snakefile
│   │-- autoseq-wgs
│   │   │--rules
│   │   │   │-- alignment.smk
│   │   │   │-- cnvcalling.smk
│   │   │   │-- germline.smk
│   │   │   │-- pre_processing.smk
│   │   │   │-- qc.smk
│   │   │   │-- somatic.smk
│   │   │   │-- structuralvariants.smk
│   │   │   │-- vep.smk
│   │   │-- Snakefile
│   │-- tumor_only
│   │   │--rules
│   │   │   │-- cnvcalling.smk
│   │   │   │-- pre_processing.smk
│   │   │   │-- qc.smk
│   │   │   │-- structuralvariants.smk
│   │   │   │-- variant_calling.smk
│   │   │   │-- vep.smk
│   │   │-- Snakefile
|-- tests
```

### Project Structure:
We are currently managing five different projects, all organized under the `/nfs/` directory. Each project follows a standardized structure, as illustrated below:

```
|-- Project1
│   │-- INBOX
│   │-- autoseq-output
│   │-- config
│   │-- logs
│   │-- sample_lists
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

Once you have prepared your cluster config file, you can start launching your samples. There are 3 essential steps to launch the pipeline. 

 * Ensure the directory name is consistent with barcode format.
 * Preparing config file
 * Launching samples.

#### Rename input directory

While running samples in batch, it is essential to create the input directory name in the format `PROJECT-SDID-TYPE-SAMPLEID-PREPID-CAPTUREID`, because, autoseq will use all the information present in this format (especially prepID and captureID) to select appropriate reference files during analysis. Please make sure to mention the date (in YYYYMMDD format) on which you are launching the sample in prepID and captureID, because, autoseq will use this information to retrive all samples that has to be launched on a specific day. To know more about this format, kindly visit [General Description](barcodes.md) page. An example of the input directory format is shown below.

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

#### Preparing config file

Once you have prepared your input directory as mentioned above, you need to create config file. The config file contains information such as SDID, Tumor ID/CFDNA ID, and Normal ID (which is same as directory name) in json format. Autoseq will use this information to search for specific directory inside `/nfs/project_name/INBOX/`. You can automatically create the config file using the command `autoseq config` which will create the config file inside `/nfs/project_name/config/YYYY-MM-DD/` with file name as `SDID.json`

You can use the following command to create config file.

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
# You can use ctrl+d to exit the screen.
```

The above script will search all the directoriers inside `/nfs/project_name/INBOX/` which are ending with today's date, sort those directory names and write them into `/nfs/project_name/sample_lists/` with file name as `clinseqBarcodes_today's_date`. It will then create a new directory with today's date inside `/nfs/project_name/config/`. Finally, the `autoseq config` script will create a seperate configuration files in json format for each input sample, and the structure of the configuration file will be as follow:

```
.
|-- config
|   |-- YYYY-MM-DD/
|   |   |-- P-*.json
|   |   |-- P-*.json
```

The contents of each configuration file will resemble the following format:

```
# P-*.json
{
    "sdid": "P-*",
    "T": ["PROJECT-P-*-T-*-KHYYYYMMDD-CYYYYMMDD"],
    "N": ["PROJECT-P-*-N-*-KHYYYYMMDD-CYYYYMMDD"],
    "CFDNA": ""
}
```

#### Launching Samples

Once you have prepared all the configuration files successfully, you can launch multiple samples in batch using the following shell script.

**Note:** It is highly recommended to start your analysis using `screen`, so that the analysis will not terminate even if you accidently close the terminal.

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
# You can use ctrl+d to exit the screen.
```

In the above script, first we are initializing reference config file path (ref), path to INBOX (libdir), path to output (outdir), number of cores, and config files (present within specific date). Then for each sample json file present in config file, we are extracting the sdid and launching autoseq pipeline with the command `autoseq launch` which will launch each job on slurm cluster. We are using `nohup` to ensure that the system does not hungup while launching job. Finally `sleep 250` is used to wait for 250-second before submitting the next job. To know more about each parameters used in `autoseq launch`, please visit [pipeline parameters](quick_start.md/#pipeline-parameters)

Once you have launched the job, you can track the status of job using the following command which will list all the jobs and their status.

```
squeue -o "%.7i %.4P %a %.60j %.20u %.8T %.10M %.9l %.6D %.6C %.6m %R"
##   [OR]
squeue -o "%.7i %.4P %a %.60j %.20u %.8T %.10M %.9l %.6D %.6C %.6m %R" | grep "username"
```

If the status of job shows anything other than `PENDING` or `RUNNING` (for example: `DependencyNeverSatisfied`), you may need to inspect the error manually or reachout to bioinformatics or IT support team. Additionally, you can check for `analysis_finished` file under `/nfs/project_name/autoseq-output/sdid/*/`. This file will be generated only if the entire pipeline gets generated successfully. You can use the following unix command to check the list of sdid that have `analysis_finished` file.

```
ls /nfs/project_name/autoseq-output/*/*/analysis_finished | \
         grep -E "sdid1|sdid2|sdid3|....|sdidn"
```

## Relaunching failed samples:
While running multiple jobs on batch, if any of the job/jobs failed, we need to fix the issue and re-launch the failed samples. But before re-launching the samples, first we need to cancel all the slurm jobs for the failed samples. You can do that with the following command.

```
jids=($(squeue -o "%j %i" | grep -E "sdid1|sdid2|sdid3|....|sdidN" | \
         cut -f 2 -d " "))
echo ${jids[@]}
for jid in ${jids[@]}; do
  scancel $jid
done
```

Then we need to ensure that there are no background jobs running. You can use the following command to check the same.

```
ps aux | grep "username"
```

If there are any background jobs running, you can kill them with the following command.

```
kill -9 $(ps aux | grep -E "sdid1|sdid2|sdid3|....|sdidN" | \
        grep -v "grep" | awk '{print $2}')
```

When snakemake starts to run any sample, it usually keeps the output directory locked for that particular sample. If the pipeline failes abruptly (due to server crash or any other reason) the sample directory will remain locked. Hence, we need to unlock such samples directories before re-launching the sample. To unlock such directories, first note down the config file path for the failed jobs and provide them in the code below. 

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
        --cluster-config /path/to/cluster_config.anchorage.json \
        --use-singularity --singularity /nfs/PIPELINE/containers/ \
        --scratch /path/to/tmp --umi --cores $cores --profile slurm \
        --smk-opt "--latency-wait 5 --unlock \
        --singularity-args '--bind /base-path/:/base-path/'"
    sleep 10
done
```

The above command will unlock all the sample directories mentioned in `configs` and we can now re-launch the failed samples with the following command. 

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
The above command will re-launch the samples mentioned in config files.