# Quick start

Installation
-------------

Dependencies:

```
* python =3.8.12
* Singularity > 3.0 
* conda (package manager)
```
To install singularity follow this [link](https://sylabs.io/guides/3.0/user-guide/installation.html)

Download and install autoseq-snakemake and the requirements using pip.

```sh
$ conda create --name autoseq_doc python=3.8.12

$ conda activate autoseq_doc

$ git clone https://github.com/Clinseq/autoseq-snakemake.git 

$ pip install -e autoseq-snakemake/

$ cd autoseq-snakemake

$ pip install "setuptools<58" --upgrade

$ pip install wheel

$ pip uninstall pyvcf

$ pip install pyvcf==0.6.8  # If any error occurs at this stage, you need to resolve them manually.

$ nano env/base.yml

# Comment out the first line "#name: base". Because we are installing all the packages in the same environment.

# To install all required tools in current environment, run the following command.
$ conda env update --file env/base.yml

# env Variable

$ GRIDSS_JAR=/path/to/autoseq-snakemake/pipeline/scripts/gridss-2.10.2-gridss-jar-with-dependencies.jar

$ GRIDSS_SCRIPT = /path/to/tools/gridss-2.13.2/

$ SAGE_JAR = /path/to/tools/sage_v3.2.3.jar
```

Pull singularity containers 

```sh
$ singularity pull --arch amd64 library://imsarath/default/autoseq-smk:latest 

$ singularity pull library://imsarath/default/gridss:latest

# There are additional singularity images present in Ravenclaw server. They are autoseq-ensemblvep.sif, autoseq-franken.sif, autoseq-gatk3.sif, autoseq-jumble.sif, autoseq-purecn.sif, autoseq-somaticseq.sif, autoseq-svcaller.sif.

# Make sure you have copied all these images before launching the pipeline.
```

Launch autoseq pipeline

```

$ autoseq launch -r autoseq-genome/autoseq-genome.json --samples /path/to/sample.json --outdir /path/to/autoseq-output/ --libdir /path/to/INBOX/ --use-singularity --singularity /path/to/container_dir --umi --cores 8 --profile slurm --smk-opt " --singularity-args '--bind /home/curator-analyst-5/projects/autoseq-snakemake/:/home/curator-analyst-5/projects/autoseq-snakemake/'"

```

While running the autoseq, if you get an error stating ModuleNotFoundError in line 4 of /path/to/autoseq-snakemake/pipeline/utils/utils.py:  No module named 'pipeline', add the following lines

```sh

# In autoseq-snakemake/pipeline/cli.py file, add the following line in line no 17 (just after importing all the packages)

os.environ["AUTOSEQ_BASE_PATH"] = os.paht.dirname(os.path.dirname(os.path.abspath(__file__)))

# Then in autoseq-snakemake/pipeline/utils/utils.py file, add the following line in line no 7 (just after importing all the packages)

sys.path.append(os.environ.get("AUTOSEQ_BASE_PATH"))
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



