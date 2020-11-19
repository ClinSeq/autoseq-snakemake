# autoseq-snakemake
snakemake version of autoseq pipeline

<h1>Autoseq framework</h1>

<b>Autoseq</b> is a custom-pipeline with additional support modules aimed primarily for the analysis of data from high-throughput sequencing of liquid biopsies. However the pipeline performs equally well with data from tissues with slight change in recommended settings. The pipeline is continuously developed and supported at Johan Lindberg's Cancer Genomics lab at Karolinska Institutet, Stockholm, Sweden. 

This is the latest port of the pipeline to snakemake framework with improvements in containerization, workflow management and upgrades in individual steps like variant calling, alignment and downstream processing of variants. Previous versions of autoseq pipeline have now been discontinued and can be found here (https://github.com/ClinSeq/autoseq). In this version, the focus has been on the forward compatibility of the tools, deployability to cloud environments and improvements in parallelization of compute-intensive and memory-heavy steps.

<h2>Structure</h2>
The repo is organized into four essential compartments -
<ul>
<li><code>rules</code> - all snakemake rule files for various steps and processes in the pipeline</li>
<li><code>config and settings</code> - configuration (*.yml) files for the pipeline and workflow settings</li>
<li><code>env</code> - all necessary conda and linux environments</li>
<li><code>scripts</code> - supporting scripts and utilities for additional functionalities for the pipeline</li>
</ul>

<h2>Pipeline</h2>
NGS-data analysis pipeline written in python for very deep targeted resequencing data from custom panels and works well with targeted or whole-exome data as well.
Contains the essential steps (QC, trimming, alignment, realignment, variant calling, prioritization, reporting). Variants from raw VCF file are annotated using VeP! and run through a set of semantic filters to eliminate irrelevant, invalid and non-significant calls. These filters (frequency, functional significance, relevance to type of disease, etc.) have been implemented after elaborate testing with different settings for both data and the callers. 
The data points in results files are suitable for manual curation in IGV after which they could be exported as a text file or HTML report or a PDF report.

<h2>Settings</h2>
Autoseq pipeline requires the input to be given in a .json format which contains the sample names following a strict naming convention. 
<p>The <code>sample.json</code> file has the format</p>
<code>
{
    "sdid": "NA12877",
    "panel": {
        "T": "NA12877-T-03098849-TD1-TT1",
        "N": "NA12877-N-03098121-TD1-TT1",
        "CFDNA": ["NA12877-CFDNA-03098850-TD1-TT1", "NA12877-CFDNA-03098850-TD2-TT1"]
    },
    "wgs": {
        "T": "NA12877-T-03098849-TD1-WGS",
        "N": "NA12877-N-03098121-TD1-WGS",
        "CFDNA": ["NA12877-CFDNA-03098850-TD1-WGS"]
    }
}
</code>






