# Quick start

Installation
-------------

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

This document covers only the installation procedure for required dependencies.

#### Installing Singularity

First, we need to download all the prerequisites using the following command.

```
sudo apt-get update
sudo apt-get install -y build-essential libssl-dev uuid-dev libgpgme11-dev \
    squashfs-tools libseccomp-dev wget pkg-config git cryptsetup debootstrap \
    libglib2.0-dev runc
```

The above command will download all the prerequisites to install singularity. We can install singularity only through GO. Hence, to install GO we need to first ensure that GO is not already installed in our system. We can delete the pre-existing GO with the following command.

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

#### Installing Miniconda

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


### Autoseq Pipeline Installation

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

$ pip3 install pulp==2.7.0

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

While running the autoseq pipeline, if you get an error stating `ModuleNotFoundError in line 4 of /path/to/autoseq-snakemake/pipeline/utils/utils.py:  No module named 'pipeline'`, add the following lines

```sh

# In autoseq-snakemake/pipeline/cli.py file, add the following line in line no 17 (just after importing all the packages)

os.environ["AUTOSEQ_BASE_PATH"] = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

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



