# Launching Samples in Server

## Launching multiple samples:

Autoseq also provides the capability to run samples in batches. To process multiple samples, first ensure that each sample's input directory is correctly formatted as described earlier. For instance, if your input directory follows this structure:

```
.
|-- INBOX/
|   |-- batch_number1/
|   |   |--DNA-B-*****-**/
|   |   |  |-- *********_DNA-B-*****-00_S11_L001_R1_001.fastq.gz
|   |   |  |-- *********_DNA-B-*****-00_S11_L001_R2_001.fastq.gz
|   |   |  |-- *********_DNA-B-*****-00_S11_L002_R1_001.fastq.gz
|   |   |  |-- *********_DNA-B-*****-00_S11_L002_R2_001.fastq.gz
|   |   |--DNA-T-*****-**/
|   |   |  |-- *********_DNA-T-*****-00_S11_L001_R1_001.fastq.gz
|   |   |  |-- *********_DNA-T-*****-00_S11_L001_R2_001.fastq.gz
|   |   |  |-- *********_DNA-T-*****-00_S11_L002_R1_001.fastq.gz
|   |   |  |-- *********_DNA-T-*****-00_S11_L002_R2_001.fastq.gz
|   |-- batch_number2/
|   |   |--DNA-B-*****-**/
|   |   |  |-- *********_DNA-B-*****-00_S11_L001_R1_001.fastq.gz
|   |   |  |-- *********_DNA-B-*****-00_S11_L001_R2_001.fastq.gz
|   |   |  |-- *********_DNA-B-*****-00_S11_L002_R1_001.fastq.gz
|   |   |  |-- *********_DNA-B-*****-00_S11_L002_R2_001.fastq.gz
|   |   |--DNA-T-*****-**/
|   |   |  |-- *********_DNA-T-*****-00_S11_L001_R1_001.fastq.gz
|   |   |  |-- *********_DNA-T-*****-00_S11_L001_R2_001.fastq.gz
|   |   |  |-- *********_DNA-T-*****-00_S11_L002_R1_001.fastq.gz
|   |   |  |-- *********_DNA-T-*****-00_S11_L002_R2_001.fastq.gz
```

You can use the following shell script to create symbolic links to all the input files in a directory with the correct naming convention:

```
for dpath in /path/to/INBOX/batch_number/DNA-*;do
    base=`basename $dpath`;
    sampletype=`echo $base | awk -F "-" '{if ($2 == "B") {print "N"} else {print $2}}'`
    sdid=`echo $base | awk -F "-" '{if (NF == 4) {print $3$4} else {print $4$5}}'
    barcode=`echo PB-P-$sdid-$sampletype-$sdid-KH$(date '+%Y%m%d')-C$(date '+%Y%m%d')`  ## added date, please check.
    mkdir /path/to/INBOX/$barcode
    ln -s $dpath/* /path/to/INBOX/$barcode/
    echo $base $barcode
done
```

**NOTE:** If your input files follow a different format, adjustments to the above script may be necessary.

The shell script above will generate symbolic links for your input files in the following structure:

```
.
|-- INBOX/
|   |-- PB-P-*****-N-*****-KH20241026-C20241026/
|   |   |-- *********_DNA-B-*****-00_S11_L001_R1_001.fastq.gz
|   |   |-- *********_DNA-B-*****-00_S11_L001_R2_001.fastq.gz
|   |   |-- *********_DNA-B-*****-00_S11_L002_R1_001.fastq.gz
|   |   |-- *********_DNA-B-*****-00_S11_L002_R2_001.fastq.gz
|   |-- PB-P-*****-T-*****-KH20241026-C20241026/
|   |   |-- *********_DNA-T-*****-00_S11_L001_R1_001.fastq.gz
|   |   |-- *********_DNA-T-*****-00_S11_L001_R2_001.fastq.gz
|   |   |-- *********_DNA-T-*****-00_S11_L002_R1_001.fastq.gz
|   |   |-- *********_DNA-T-*****-00_S11_L002_R2_001.fastq.gz
|   |-- PB-P-*****-N-*****-KH20241026-C20241026/
|   |   |-- *********_DNA-B-*****-00_S11_L001_R1_001.fastq.gz
|   |   |-- *********_DNA-B-*****-00_S11_L001_R2_001.fastq.gz
|   |   |-- *********_DNA-B-*****-00_S11_L002_R1_001.fastq.gz
|   |   |-- *********_DNA-B-*****-00_S11_L002_R2_001.fastq.gz
|   |-- PB-P-*****-T-*****-KH20241026-C20241026/
|   |   |-- *********_DNA-T-*****-00_S11_L001_R1_001.fastq.gz
|   |   |-- *********_DNA-T-*****-00_S11_L001_R2_001.fastq.gz
|   |   |-- *********_DNA-T-*****-00_S11_L002_R1_001.fastq.gz
|   |   |-- *********_DNA-T-*****-00_S11_L002_R2_001.fastq.gz
```


After creating the symbolic links for all your input samples using the script above, you can generate configuration files for each of these samples using the following shell script:

```
find /path/to/INBOX/ -maxdepth 1 -name "SARC*$(date '+%Y%m%d')" | xargs -I {} basename {} | sort -V > /path/to/sample_lists/clinseqBarcodes_`date "+%Y-%m-%d"`.txt   ## added date, please check.
mkdir -p /path/to/config/$(date "+%Y-%m-%d")
/path/to/autoseq-snakemake/autoseq config --outdir /path/to/config/$(date "+%Y-%m-%d") /path/to/sample_lists/clinseqBarcodes_$(date "+%Y-%m-%d").txt
```

This script will create configuration files for each input sample in the following structure:

```
.
|-- config
|   |-- 2024-10-26/
|   |   |-- P-*******.json
|   |   |-- P-*******.json
```

The contents of each configuration file will resemble the following format:

```
# P-*******.json
{
    "sdid": "P-*******",
    "T": ["PB-P-*****-T-*****-KH20241026-C20241026"],
    "N": ["PB-P-*****-N-*****-KH20241026-C20241026"],
    "CFDNA": ""
}

# P-*******.json
{
    "sdid": "P-*******",
    "T": ["PB-P-*****-T-*****-KH20241026-C20241026"],
    "N": ["PB-P-*****-N-*****-KH20241026-C20241026"],
    "CFDNA": ""
}
```

Once all configuration files have been successfully created, you can launch multiple samples in batches using the following shell script:

```
screen -S autoseq_run
prod_up
ref=/path/to/autoseq-genome/autoseq-genome.json
libdir=/path/to/INBOX/
outdir=/path/to/autoseq-output/
cores=8
configs=(`find /path/to/config/$(date "+%Y-%m-%d") -name "P-*json" | sort -r`)
for config in ${configs[@]}; do
    echo $config
    sdid=`basename $config |cut -f 1 -d "."`;
    echo $sdid
    nohup autoseq launch -r $ref --samples $config --outdir $outdir --libdir $libdir \
        --cluster-config /path/to/autoseq-snakemake/pipeline/scheduler/cluster_config_specific_server.json \
        --use-singularity --singularity /path/to/containers/ \
        --scratch /path/to/tmp --umi --cores $cores --profile slurm \
        --smk-opt "--latency-wait 60 " >> logs/$sdid.nohub.log &
    sleep 250
done
```

This script will submit your jobs to the Slurm cluster sequentially, with a 250-second interval between each job submission.

## Relaunching failed samples:
While running multiple jobs on batch, if any of the job/jobs failed, we can use the following set of command to re-launch the failed samples. Before re-launching the samples, first we need to ensure that there are no background jobs running. You can use the following command to check the same.

```
ps aux | grep "username"
```

If there are any background jobs running, you can kill them with the following command.

```
kill -9 $(ps aux | grep -E "sdid1|sdid2|sdid3|....|sdidN" | grep -v "grep" | awk '{print $2}')
```

Usually, when snakemake starts to run any sample and if it failes abruptly (due to server crash or any other reason), it usually keeps the output directory locked for the particular sample. To unlock the directory, you can use the following set of command.

```
screen -r autoseq_run
prod_up   # run this only if production environment is not active.
ref=/path/to/autoseq-genome/autoseq-genome.json
libdir=/path/to/INBOX/
outdir=/path/to/autoseq-output/
cores=8
configs=(/path/to/config/date/sdid1.json /path/to/config/date/sdid2.json /path/to/config/date/sdid3.json ... /path/to/config/date/sdidN.json)  ## Provide the config files for failed samples here.

for config in ${configs[@]}; do
    echo $config
    sdid=`basename $config |cut -f 1 -d "."`;
    echo $sdid
    nohup autoseq launch -r $ref --samples $config --outdir $outdir --libdir $libdir \
        --cluster-config /nfs/PIPELINE/autoseq-snakemake/pipeline/scheduler/cluster_config.anchorage.json \
        --use-singularity --singularity /nfs/PIPELINE/containers/ \
        --scratch /nfs/KODIAK2/PSFF/re-run/tmp --umi --cores $cores --profile slurm \
        --smk-opt "--latency-wait 5 --unlock"
    sleep 10
done
```

The above command will unlock all the sample directories mentioned in the configs and we can now re-launch the failed samples with the following command. 

```
screen -r autoseq_run
prod_up
ref=/path/to/autoseq-genome/autoseq-genome.json
libdir=/path/to/INBOX/
outdir=/path/to/autoseq-output/
cores=8
configs=(/path/to/config/date/sdid1.json /path/to/config/date/sdid2.json /path/to/config/date/sdid3.json ... /path/to/config/date/sdidN.json)  ## Provide the config files for failed samples here.

for config in ${configs[@]}; do
    echo $config
    sdid=`basename $config |cut -f 1 -d "."`;
    echo $sdid
    nohup autoseq launch -r $ref --samples $config --outdir $outdir --libdir $libdir \
        --cluster-config /nfs/PIPELINE/autoseq-snakemake/pipeline/scheduler/cluster_config.anchorage.json \
        --use-singularity --singularity /nfs/PIPELINE/containers/ \
        --scratch /nfs/KODIAK2/PSFF/re-run/tmp --umi --cores $cores --profile slurm \
        --smk-opt "--latency-wait 60 --rerun-incomplete"
    sleep 250
done
```
The above command will re-launch the samples mentioned in config files.