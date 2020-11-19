# autoseq-snakemake
snakemake version of autoseq pipeline

<h1>Autoseq framework</h1>

<b>Autoseq</b> is a custom-pipeline with additional support modules aimed primarlity for the analysis of data from high-throughput sequencing of liquid biopsies. However the pipeline performs equally well with data from tissues with slight change in recommended settings. The pipeline is continuously developed and supported at Johan Lindberg's Cancer Genomics lab at Karolinska Institutet, Stockholm, Sweden. 

This is the latest port of the pipeline to snakemake framework with improvements in containerization, workflow management and upgrades in individual steps like variant calling, alignment and downstream processing of variants. Previous versions of autoseq pipeline have now been discontinued and can be found here (https://github.com/ClinSeq/autoseq). In this version, the focus has been on the forward compatibility of the tools, deployability to cloud environments and improvements in parallelization of compute-intensive and memory-heavy steps.

<h2>Structure</h2>
The repo is organized into four essential compartments -
<ul>
<li><code>rules</code> - <code>all snakemake rule files for various steps and processes in the pipeline</code></li>
<li><code>config and settings</code> - <code>configuration (*.yml) files for the pipeline and workflow settings</code></li>
<li><code>env</code> - <code>all necessary conda and linux environments</code></li>
<li><code>scripts</code> - <code>supporting scripts and utilities for additional functionalities for the pipeline</code></li>
</ul>










