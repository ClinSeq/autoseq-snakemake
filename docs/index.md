# Autoseq: A clinical genomics workflow management platform to analyze next generation sequencing data from cancer samples.

Autoseq is a custom pipeline with supplementary support modules, designed primarily for the analysis of data from high-throughput sequencing of liquid biopsies. However, the pipeline performs equally well with tissue-derived data, with only minor adjustments to the recommended settings. The pipeline is continuously developed and supported by Johan Lindberg's Cancer Genomics Lab at Karolinska Institutet, Stockholm, Sweden.

This pipeline is used to analyze cell-free DNA from blood samples in the PROBIO study. There are two versions of the pipeline: one with standard data analysis steps and another that includes additional UMI processing. Currently, this pipeline supports only the GRCh37 reference genome. 

This pipeline reports different types of variants such as 

* Small Somatic Variants
* Germline Variants
* Copy Number Variants
* Structural Variants

To know more about different tools that we use in this pipeline, visit [autoseq pipeline page](autoseq_pipeline.md). 

|||
|------------------|-----------------------------|
|Source code | [GitHub](https://github.com/ClinSeq/autoseq-snakemake/)|
|License | [Apache License 2.0](http://www.apache.org/licenses/LICENSE-2.0)|
|Packages | [Python >3.6](https://www.python.org/downloads/), [Snakemake](https://snakemake.readthedocs.io/en/stable/), [Conda](https://docs.conda.io/en/latest/), [Singularity](https://cloud.sylabs.io/)|
|Q&A | [Questions & Answers](https://github.com/ClinSeq/autoseq-snakemake/issues)|
|Latest Version | [v3.5.2](https://github.com/ClinSeq/autoseq-snakemake/releases/tag/v3.5.2)|


**Authors**

* [Sarath Kumar Murugan](https://github.com/imsarath) (sarath.murugan@ki.se)
* [Venkatesh Chellappa Patel](https://github.com/drvenki) (venkatesh.chellappa@ki.se)
* [Rebecka Bergstrom](https://github.com/rebber) (rebecka.bergstrom@ki.se)
* [Markus Mayrhofer](https://github.com/mayrhofer) (markus.mayrhofer@ki.se)

**Index**

* [Quick Start](quick_start.md)
    * [Requirements](quick_start.md/#11-requirements)
    * [Installation](quick_start.md/#12-installation)
    * [Launching Autoseq Pipeline](quick_start.md/#13-launching-autoseq-pipeline)
* [Launching Samples in Server](launch_pipeline.md)
    * [Breif Description About Our Server](launch_pipeline.md#breif-description-about-our-server)
    * [Workload Manager: Slurm](launch_pipeline.md#workload-manager)
    * [Data Organization in Ravenclaw](launch_pipeline.md#data-organization-in-the-hpc)
    * [Virtual Environment](launch_pipeline.md#virtual-environment)
    * [Launching multiple samples in server](launch_pipeline.md#launching-multiple-samples-on-the-server)
    * [Relaunching failed samples](launch_pipeline.md#relaunching-failed-samples)
* [Autoseq Pipeline](autoseq_pipeline.md)
* [General Description](barcodes.md)
* [Autoseq scripts](scripts.md)
