# Autoseq - Pipeline

Autoseq pipeline is specifically designed to run liquid biopsy samples, however, it performs equally well with tissue biopsy samples with minor modification to the recommended settings. This pipeline requires both tumor and matched normal samples; and the input fastq file has to be either in the format of fastq.gz or fq.gz. Additionally, the input file name has to be specified in the format `PROJECT-SDID-TYPE-SAMPLEID-PREPID-CAPTUREID`. For example: PB-P-00462065-CFDNA-04055058-KH20221214-C420221214.fastq.gz. To know more about this format, visit [General Description](barcodes.md) page. Currently, we have 5 different pipeline in our workflow. They are 

1. Autoseq Pipeline
2. Autoseq-rerun Pipeline
3. Autoseq-wgs Pipeline
4. Autoseq-sd Pipeline
5. Tumor-only Pipeline

## 1. Autoseq Pipeline

Autoseq pipeline is specifically designed for targeted variant calling analysis. Different capture kits focuses on different target genes, for example, in probio, we are targeting 94 genes. If you are using any new capture kit, you can specify the corresponding reference files such as bed, interval_list, etc under `/nfs/PIPELINE/autoseq-genome/intervals/targets/` and specify the same in autoseq-genome.json file in the below format. Also, remember to update `capture_kit_loopkup` dictionary in `utils.py` script so that autoseq pipeline can automatically take the reference file during analysis. 

```
"targets": {
        "alascca_core": {
            "blacklist-bed": null, 
            "msings-baseline": null, 
            "msings-bed": null, 
            "msings-msi_intervals": null, 
            "msisites": "intervals/targets/alascca_core.slopped20.msisites.tsv", 
            "purecn_targets": null, 
            "targets-bed-slopped20": "intervals/targets/alascca_core.slopped20.bed", 
            "targets-bed-slopped20-gz": "intervals/targets/alascca_core.slopped20.bed.gz", 
            "targets-interval_list": "intervals/targets/alascca_core.interval_list", 
            "targets-interval_list-slopped20": "intervals/targets/alascca_core.slopped20.interval_list"
        },
        "your_new_capture_kit_name": {
            "blacklist-bed": null, 
            "msings-baseline": null, 
            "msings-bed": null, 
            "msings-msi_intervals": null, 
            "msisites": "/path/to/your_new_capture_kit.slopped20.msisites.tsv", 
            "purecn_targets": null, 
            "targets-bed-slopped20": "/path/to/your_new_capture_kit.slopped20.bed", 
            "targets-bed-slopped20-gz": "/path/to/your_new_capture_kit.slopped20.bed.gz", 
            "targets-interval_list": "/path/to/your_new_capture_kit.interval_list", 
            "targets-interval_list-slopped20": "/path/to/your_new_capture_kit.slopped20.interval_list"
        } 
}
```

While running the pipeline, you can run it either with or without umi parameter. If you specify umi parameter, additional UMI preprocessing step will be performed.

This pipeline can be broadly classified into 7 major steps. They are

* Preprocessing and UMI Processing
* Germline Variant Calling
* Somatic Variant Calling
* Copy Number Variant Calling
* Structural Variant Calling
* Microsatellite Instability
* QC steps 

![Autoseq workflow](img/autoseq_overall_diagram.png)
**Figure 1a:** This diagram shows the overall workflow for autoseq pipeline. Here different type of processes are highlighed with different colors. The preprocessing, UMI processing, microsatellite instability, quality check, somatic variant calling, germline variant calling, copy number analysis, structural variant calling and plots are mentioned in grey, orange, plum, blue, dark orange, green, light blue, dark yellow and light orange respectively. 


**Preprocessing and UMI Processing**

This step covers various sub-steps such as trimming, UMI processing, bwa alignment, indel realignment etc. For trimming the input fastq file, we are using a tool called skewer. It automatically identifies adapter sequence and removes them from each reads. Once trimming process is completed, the resulting trimmed fastq file will be passed to UMI processing step. For UMI pre-procssing we are using fgbio's groupreadsbyumi, callduplexconsensus, and filterconsensus. During UMI preprocessing stage, the input fastq files be converted to ubam file. Then the bwa alignment will be carried out for the trimmed fastq files and the resulting bam file will be merged with the previously generated ubam file, so that the UMI information will be preserved in the final bam file. Then an indel realignment step is performed in the bam file with the help of gatk's targetcreator and indelrealigner to improve the accuracy of alignment near indel regions.

Once the bam file with UMI information is generated and indels regions are re-aligned, we will group reads based on UMI using fgbio's groupreadsbyumi and consequtively call consensus reads with fgbio's callduplexconsensus. The resulting bam file will then be converetd to fastq file, bwa aligned, and indel regions will be re-aligned.    

Once we group reads based on consensus, we can observe that true variants are present in most of the reads that has same UMI, where as false variants that arise due to sequencing error will present only in a small fractions of reads that has same UMI. Fgbio's FilterConsensusReads uses this and few other concept to filter variants that are probably arising from sequencing artifacts. To know more about all the filters that are applied, you can visit [fgbio's FilterConsensusReads page](http://fulcrumgenomics.github.io/fgbio/tools/latest/FilterConsensusReads.html). Finally, fgbio's clipbam is used to clip N bases from the end of either R1 or R2 so that variants present in end of consensus reads are not counted twice (if not clipped, the variants present in both R1 and R2 will be counted as 2 variants; though technically only 1 variant is present in forward and reverse strand.)

Once the preprocessing and UMI processing steps are completed, the remaining steps can be run in parallel.

**Germline Variant Calling**

For normal sample, we use [GATK's haplotypecaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) to identify all germline variants (including SNPs and INDELs). This tool locates the positions where there is a probable mismatch in bam file and reassembles reads in that region, thus calling variants with increased accuracy. Once the variant call file is generated, it is annotated using [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html).

**Somatic Variant Calling**

To identify somatic variants, we use two different tools. They are [GATK's Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) and [SAGE somatic](https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md). Mutect2 uses a Bayesian somatic genotyping model to identify variants, while SAGE somatic algorithm uses 9 keys steps to call somatic variants. They are Alt Specific Base Quality Recalibration, Candidate Variants, Tumor Counts and Quality, Normal Counts and Quality, Soft Filters, Phasing, De-duplication and Gene Panel Coverage. The variants identified from these two tools will then be merged using [somaticseq](https://github.com/bioinform/somaticseq/blob/master/docs/Manual.pdf) which uses machine learning model to filter out false positive from these two variant callers.

**Copy Number Variant Calling**

For calling copy number variants, we use a tool called JUMBLE. Once copy number analysis is completed, autoseq pipeline will also generate two plots. They are liqbiocna plot and franken plot.  

**Structural Variant Calling**

For structural variant calling, we use [svcaller](https://github.com/tomwhi/svcaller) and [GRIDSS](https://github.com/PapenfussLab/gridss). svcaller uses an algorithm similar to Delly. GRIDSS calls variants based on alignment-guided positional de Bruijn graph genome-wide break-end assembly, split read, and read pair evidence.

**Microsatellite Instability**

Microsatellite Instability provides a glimps into the underlying issues with DNA mismatch repair mechanism. To detect the MSI status, we use two tools called [msiNGS](https://bitbucket.org/uwlabmed/msings/src/master/) and [MSIsensor](https://github.com/ding-lab/msisensor). msiNGS compares the number of differently sized repeats within the microsatellite locus of tumor and normal sample. Loci were considered unstable if the number of repeats in tumor was statistically greater than the number of repeats in normal. MSIsensor uses Pearson's Chi-Squared Test to compare the number of repeats in microsatellite locus between tumor and normal sample.

**QC steps**

Autoseq pipeline utilizes various tools to ensure the generated results in various stages are of good quality. In the initial fastq file, we use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to check the quality of input reads. Once the fastq files are aligned, we use [samtools flagstats](http://www.htslib.org/doc/samtools-flagstat.html) to check the overall quality of aligned BAM file. Further, we use various tools in Picard such as [CollectHsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036856051-CollectHsMetrics-Picard) to identify GC bias, [MarkDups](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard) to identify duplicate reads, [CollectInsertSizeMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard) and [collectOxoGMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037428231-CollectOxoGMetrics-Picard) to validate the library construction process. GATK's [ContEst](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tooldocs/3.6-0/org_broadinstitute_gatk_tools_walkers_cancer_contamination_ContEst.html) tool is used to estimate the level of cross contamination between sequencing data. And, finally, purecn is used to estimate tumor purity and ploidy, copy number and loss of heterozygosity.

### Pipeline Detailed Workflow

![Autoseq Complete Workflow](img/autoseq_detailed_workflow.png)
**Figure 1b:** Complete DAG diagram for the autoseq pipeline is shown in this diagram. Here, the postprocessing step is highlighted in light yellow. All other color code is similar to the above diagram.



## 2. Autoseq Re-run Pipeline

This pipeline was designed specifically to start analysis from bam file. If the pipeline gets aborted due to any reason after completing the alignment process, or if we do not have fastq file, rather have duplicates removed bam file available, we can use this pipeline. All the tools used here are same as autoseq pipeline.

![Autoseq rerun](img/autoseq_rerun.png)
**Figure 2:** This diagram shows the overall workflow for autoseq rerun pipeline. Here different type of processes are highlighed with different colors. The preprocessing, quality check, somatic variant calling, germline variant calling, copy number analysis, structural variant calling and plots are mentioned in grey, blue, dark orange, green, light blue, dark yellow and light orange respectively.

### Launching Rerun Pipeline:

Launching autoseq-rerun pipeline is similar to autoseq pipeline; however, here you need to specify the path to bam file in `--libdir` parameter. Additionally, you need to specify `--pipeline autoseq-rerun`. Below, you can find an example command to launch autoseq-rerun pipeline.

```
autoseq launch -r /path/to/autoseq-genome/autoseq-genome.json \
        --samples /path/to/sample.json \ 
        --outdir /path/to/autoseq-output/ \
        --libdir /path/to/autoseq-output/sdid/ \ 
        PROJECT-SDID-CFDNA-SAMPLEID-PREPID-CAPTUREID_\
        PROJECT-SDID-N-SAMPLEID-PREPID-CAPTUREID/bams/ \
        --use-singularity --singularity /path/to/container_dir \
        --cores 8 --pipeline autoseq-rerun --smk-opt " --singularity-args\
        '--bind /base-path/:/base-path/'"
```

## 3. Autoseq WGS Pipeline

This pipeline is used to analyze whole genome sequencing data. Unlike autoseq pipeline we are not using unique molecular identifier (UMI) in WGS pipeline, and we are running haplotype caller (rule gatk4_haplotypecaller_jc) on tumor sample as well to identify germline variants that were present in tumor sample. The results of this analysis will be used in copy number analysis for generating liqbio plot and franken plot. Additionally, in quality metrics, we are using WGSmetrics instead of CollectHsMetrics. This pipeline also requires tumor as well as matched normal sample for analysis. All the other steps used in this pipeline is similar to autoseq pipeline. 

![Autoseq WGS workflow](img/autoseq_wgs.png)
**Figure 3:** This diagram shows the overall workflow for autoseq whole genome sequencing pipeline. Here different type of processes are highlighed with different colors. The preprocessing, quality check, somatic variant calling, germline variant calling, copy number analysis, structural variant calling and plots are mentioned in grey, blue, dark orange, green, light blue, dark yellow and light orange respectively.

### Launching WGS Pipeline:

Launching WGS pipeline is similar to launching autoseq pipeline, you just have to specify an additional parameter which is `--pipeline autoseq-wgs`. Below, you can find an example command to launch WGS sample.

```
autoseq launch -r /path/to/autoseq-genome/autoseq-genome.json \
        --samples /path/to/sample.json \ 
        --outdir /path/to/autoseq-output/ \ 
        --libdir /path/to/input_directory/ --use-singularity \ 
        --singularity /path/to/container_dir \
        --cores 8 --pipeline autoseq-wgs --smk-opt " --singularity-args\
        '--bind /base-path/:/base-path/'"
```

### Launching WGS Pipeline in Ravenclaw Server:

In KI, we are receving WGS samples in batch wise and their directory name does not follow the [recommended](quick_start/barcodes.md) format. Hence, we need to create new directory inside input directory with correct format, and create symbolic link of all input fastq files inside this directory. You can use the following command to do the same.

```
for dpath in /path/to/INBOX/batch_5/DNA-*;do
	base=`basename $dpath`;
	sampletype=`echo $base | awk -F "-" '{if ($2 == "B") {print "N"} \
                else {print $2}}'`
	sdid=`echo $base | awk -F "-" '{if (NF == 4) {print $3$4} \
                else {print $4$5}}' | sed -e "s/WGS//g"`
	barcode=`echo SARC-P-$sdid-$sampletype-$sdid-\
                KH$(date '+%Y%m%d')-WG$(date '+%Y%m%d')`
	mkdir /path/to/INBOX/$barcode
	ln -s $dpath/* /path/to/INBOX/$barcode
	echo $base $barcode
done
```

Now, your input directories are in correct format, and all the fastq files were symbolic linked inside the input directory. You can now proceed further to prepare config files by following the procedure mentioned in [launching samples in server](server_launch.md/#preparing-config-file) page.

## 4. Autoseq Tumor Only Pipeline:

This pipeline can be used when we do not have a matched normal for your tumor sample. For such samples, we need to provide a Panel Of Normal (PON) bam file with the option `--normal /path/to/panel_of_normal.bam`. In Karolinska Institute, we have specifically designed panel of normal for each projects that we are handling. These data can be downloaded upon request.

![Autoseq tumor only workflow](img/tumor_only.png)
**Figure 4:** This diagram shows the overall workflow for autoseq whole genome sequencing pipeline. Here different type of processes are highlighed with different colors. The preprocessing, quality check, somatic variant calling, germline variant calling, copy number analysis, structural variant calling and plots are mentioned in grey, blue, dark orange, green, light blue, dark yellow and light orange respectively.

### Launching Tumor-Only Pipeline:

Launching tumor-only pipeline is similar to launching autoseq pipeline, you just have to specify two additional parameters which are`--pipeline tumor_only` and the location of panel of normal bam file with the parameter `--normal-bam /path/to/panel_of_normal/bams/`. Below, you can find an example command to launch tumor only sample.

```
autoseq launch -r /path/to/autoseq-genome/autoseq-genome.json \
        --samples /path/to/sample.json \ 
        --outdir /path/to/autoseq-output/ \ 
        --libdir /path/to/input_directory/ --use-singularity \ 
        --singularity /path/to/container_dir \
        --pipeline tumor_only --normal-bam /path/to/panel_of_normal/bams/ \
        --cores 8 --pipeline tumor_only --smk-opt " --singularity-args\
        '--bind /base-path/:/base-path/'"
```