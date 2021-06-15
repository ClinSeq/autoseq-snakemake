# Quick start

Installation
-------------
Dependencies:
* python >3.6
* Singularity > 3.0 

To install singularity [follow this link](https://sylabs.io/guides/3.0/user-guide/installation.html)

Download and install autoseq-snakemake and the requirements using pip.
```sh
git clone https://github.com/Clinseq/autoseq-snakemake.git 
pip install autoseq-snakemake/
```

Pull singularity container 
```sh
singularity pull --arch amd64 library://imsarath/default/autoseq-smk:latest 
```

Command Line Usage
------------------

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

Prepare the reference files
-----------------------------


LiqBio Pipeline
---------------

```sh
autoseq --ref ref.json --outdir /path/to/outdir --jobdb jobdb.json --cores 5 --runner_name slurmrunner --libdir /path/to/libdir liqbio sample.json
```

