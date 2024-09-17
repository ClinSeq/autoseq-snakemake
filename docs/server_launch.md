## Our Server:
In Karolinska Institutet, we have a High Performance Computer called ravenclaw. It consists of one head node and 5 compute nodes. The head node has 12 CPUs with 2 cores in each CPUs and 64GB of RAM. Each of the compute nodes have 48 CPUs with 252GB RAM and 250TB of disk space. 

## Workload Manager:
In order to manage all these compute nodes effectively, we need a special software called workload manager which can effectively allocate resources for each job for specific amount of time, track, monitor, and report status of job. In our autoseq pipeline, we are using a workload manager called Slurm which takes care of the resource allocation. While working with HPC, we have to login only to the head node to submit our jobs. Slurm will in turn take care of assigning the job to each of the compute nodes depending on the resources requested by the users. Once the analysis is completed by the compute node, it will send status report to the head node. Thus the entire analysis are effectively managed by Slurm. If you are using a different server, and if you don't have Slurm installed in it, kindly contact your IT support team to install Slurm. 

## Data Organization in Ravenclaw:
While working on server, it is crutial to organize all the projects in specific structure such that it is easy to understand and access each projects. Here we have discussed about how data is organization in Ravenclaw server. 

All the projects and autoseq pipeline are stored in an extensible SSD under `/nfs/`. The directory structure is as follow
```
|-- PIPELINE
|-- Project1
|-- Project2
.
.
.
```

### Pipeline Structure:
Autoseq pipeline is located under `/nfs/PIPELINE/` and it has the following structure.
```
|-- autoseq-genome
|-- autoseq-snakemake
|-- containers
|-- tools
```

The `autoseq-genome` directory contains the reference genome (GRCh37) and all the associated files. The `autoseq-snakemake` directory contains all the pipelines, their conda environment, documentation and some test samples to check the pipeline. The `containers` has all the singularity images used in the pipeline and the `tools` directory contains all the additional tools used in autoseq pipeline. The `autoseq-snakemake` directory contains further subdirectories for each of the 5 pipeline we use. Inside each pipeline we have `rules` directory containing all rules and a `Snakefile`. Their structure is shown below.

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
So far we have been working on x different projects and each of these projects are listed under `/nfs/`. All of the projects follow similar structure, one such example is shown below.
```
|-- Project1
│   │-- INBOX
│   │-- autoseq-output
│   │-- config
│   │-- logs
│   │-- sample_lists
```

The `INBOX` contains all the sample fastq files; their directory name has to follow pattern mentioned in [General Description page](barcodes.md). The `autoseq-output` directory contains the results of the analysis. `config` directory contains the configuration files (which contains tool specific parameters) used in the pipeline. The `logs` directory contains the no hungup analysis log for each sample. And, the `sample_lists` directory contains list of files with file name as `clinseqBarcodes_YYYY-MM-DD.txt`. Each of these files contains list of samples for which analysis were started on the specific date.

**Note:** If you are starting analysis for any new project in our server or creating a new pipeline, it is crutial to maintain the above directory structure. 

## Virtual Environment

In order to make the analysis easier, we have configured a virtual environment in our ravenclaw server with all the required dependencies. You can use the command `prod_up` to activate this virtual environment.

## Launching multiple samples in server:

If you wish to launch multiple samples in server, it essential to ensure that you are submitting your jobs through Slurm, so that Slurm can take care of allocating resources for each of your job. We have already installed Slurm and configured the virtual environment called `prod_up` in our ravenclaw server. If you are using a different server and if you wish to specify any memory requirement for any specific job, you need to modify the cluster configuration file as per your requirement. The cluster configuration file can be found in `/nfs/PIPELINE/autoseq-snakemake/pipeline/scheduler/cluster_config.json`

Here is an example of cluster configuration file structure.

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

Here, the `__default__` will be applied to all the rules, and the parameters which you have specified inside each rule will overwrite the parameters provided in `__default__`. You can pass this config file to autoseq pipeline with `--cluster-config` parameter.

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
|   |-- PROJECT-P-*-T-*-KHYYYYMMDD-CYYYYMMDD/
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
autoseq config --outdir /path/to/config/$(date "+%Y-%m-%d") \
        /path/to/sample_lists/clinseqBarcodes_$(date "+%Y-%m-%d").txt
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

In the above script, first we are initializing reference config file path (ref), path to INBOX (libdir), path to output (outdir), number of cores, and config files. Then for each sample json file present in config file, we are extracting the sdid and launching autoseq pipeline with the command `autoseq launch` which will launch each job on slurm cluster. We are using `nohup` to ensure that the system does not hungup while launching job. Finally `sleep 250` is used to wait for 250-second before submitting the next job. To know more about each parameters used in `autoseq launch`, please visit [pipeline parameters](quick_start.md/#pipeline-parameters)

Once you have launched the job, you can track the status of job using the following command which will list all the jobs and their status.

```
squeue -o "%.7i %.4P %a %.60j %.20u %.8T %.10M %.9l %.6D %.6C %.6m %R"
```

If the status of job shows anything other than `PENDING` or `RUNNING` (for example: `DependencyNeverSatisfied`), you may need to inspect the error manually or reachout to bioinformatician or IT support team. Additionally, you can check for `analysis_finished` file under `/nfs/project_name/autoseq-output/sdid/*/`. This file will be generated only if the entire pipeline gets generated successfully. You can use the following unix command to check the list of sdid that have `analysis_finished` file.

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