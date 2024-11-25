# Autoseq internal scripts 


### create_contest_vcfs.py

Create the contamination test VCF input file from the designated "genotype" and "eval" target bed region files and population allele frequency VCF file.

```
create_contest_vcfs.py normal.slopped20.bed tumor.slopped20.bed \
        swegen_common.vcf.gz --tmpdir /path/to/temp-dir/ \
        --output-filename /path/to/output_file.vcf 2> /path/to/log_file.log
```

### vcf_add_sample.py

Add DP, RO and AO for a new sample from a BAM file. Filters variants that are not simple substitutions. Optionally filter homozygous variants from the output file (--filter_hom)

```
usage: vcf_add_sample.py --samplename tumor_id variants.vcf(.gz) \
         aligned.bam > new.vcf
```

### contest_to_contam_caveat.py

Generate a call on contamination from contest output, and output it as a JSON file.

Example output contents: `{"CALL": "OK"}`

```
contest_to_contam_caveat.py /path/to/contamination/id.contest.txt \
        > /path/to/outputdir//qc/id-contam-qc-call.json 2> \
        /path/to/log_file.log
```

### generate_allelic_fraction_bedGraph.py

For each variants in the VCF file, generates allelic fraction. 

```
generate_allelic_fraction_bedGraph.py --output \
        /path/to/variants/id-and-id.germline-variants-with-somatic-afs.vcf.gz \
        /path/to/variants/id-and-id.germline-variants-with-somatic-afs.vcf.gz 

```

### generateIGVnavInput.py

Script to create list of SNVs from merged somatic or germline vcfs with oncogenicity annoation for IGVNav.

```
generateIGVnavInput.py /path/to/variants/id-id-all.somatic.vep.vcf \
        OncoKB_6Mar19_v1.9.txt somatic \
        --cgc variant_annotation_for_curator_v1_2023-02-17.txt \
        --output id-id-igvnav-input.txt
```

### generateIGVnavInput_SV.py

Script to create mut files from different structural variant caller and combine all mut files with gene annotation.

```
generateIGVnavInput_SV.py --input /path/to/svs/svcaller/ \
        --sdid tumor_sdid --tool svcaller \
        --vcftype somatic --output /path/to/svs/igv/
```

### generate_symlinks.py

Generates symlinks required for IGVnav inputs - curations.

```
generate_symlinks.py --targets targets-bed-slopped20.bed \
        --outdir /path/to/outdir/
```

### liqbioCNA_Interactive_plots.R

To create static and interactive franken plot.

```
liqbioCNA_Interactive_plots.R  --tumor_cnr /path/to/cnv/tumor.cnr \
        --tumor_cns /path/to/cnv/tumor.cns \
        --normal_cnr /path/to/cnv/tumor.cnr \
        --normal_cns /path/to/cnv/tumor.cns \
        --het_snps_vcf id-and-id.germline-variants-with-somatic-afs.vcf.gz \
        --purecn_csv /path/to/purecn/tumor.csv \
        --purecn_genes_csv /path/to/purecn/tumor_purecn_genes.csv \
        --purecn_loh_csv /path/to/purecn/tumor_purecn_loh_csv \
        --purecn_variants_csv /path/to/purecn/tumor_purecn_variants_csv \
        --svcaller_T_DEL /path/to/svs/svcaller/id-tumor-DEL.gtf \
        --svcaller_T_DUP /path/to/svs/svcaller/id-tumor-DUP.gtf \
        --svcaller_T_INV /path/to/svs/svcaller/id-tumor-INV.gtf \
        --svcaller_T_TRA /path/to/svs/svcaller/id-tumor-TRA.gtf \
        --svcaller_N_DEL /path/to/svs/svcaller/id-normal-DEL.gtf \
        --svcaller_N_DUP /path/to/svs/svcaller/id-normal-DUP.gtf \
        --svcaller_N_INV /path/to/svs/svcaller/id-normal-INV.gtf \
        --svcaller_N_TRA /path/to/svs/svcaller/id-normal-TRA.gtf \
        --germline_mut_vcf /path/to/id-all.germline.vep.vcf \
        --somatic_mut_vcf /path/to/id-all.somatic.vep.vcf \
        --plot_png /path/to/qc/id-tumor-liqbio-cna.png \
        --plot_png_normal /path/to/qc/id-normal-liqbio-cna.png \
        --cna_json /path/to/id-liqbio-cna.json \
        --purity_json /path/to/qc/id-liqbio-purity.json \
        --gene_track gene_track_output.csv
```

### PureCN.R

Customized script to estimate tumor purity, copy number, and loss of heterozygosity (LOH), and classifies single nucleotide variants (SNVs) by somatic status and clonality. 

```
PureCN.R  --out /path/to/purecn/ \
        --sampleid tumor_id \
        --segfile /path/to/output/cnv/ID_dnacopy.seg \
        --tumor /path/to/output/cnv/ID.cnr \
        --vcf /path/to/output/variants/vardict-somatic-purecn.vcf.gz \
        --genome hg19  --funsegmentation none \
        --minpurity 0.05  --hzdev 0.1 --maxnonclonal 0.2 --minaf 0.01 \ 
        --error 0.0005 \
        --postoptimize
```

### QC_overview.R

This is a script creating plots showing the distribution of some QC values, and where the current normal and tumor sample lies within those distributions. QC metrics shown are read count, duplication rate, coverage, fold enrichment, on-bait rate, insert size & contamination

```
QC_overview.R -s normal_id:tumor_id -d /path/to/output_directory/ \
        -o /path/to/output/qc/ID.qc_overview.pdf
        -m /path/to/output/
```

### frankenscript.R

This script is used to generate franken plot for each chromosome.

```
frankenscript.R  --tumor_cnr tumor.cnr \
        --tumor_cns tumor.cns \
        -normal_cnr normal.cnr \
        --normal_cns normal.cns \
        --frankenplot_Rmd /path/to/scripts/frankenplot.Rmd \
        --het_snps_vcf normal-and-tumor.germline-variants-with-somatic-afs.vcf.gz \
        --purecn_csv /path/to/purecn/cancer.csv \
        --purecn_genes_csv /path/to/purecn/cancer_genes.csv \
        --purecn_loh_csv  /path/to/purecn/cancer_loh.csv \
        --purecn_variants_csv  /path/to/purecn/cancer_variants.csv \
        --svcaller_T_DEL /path/to/svs/svcaller/tumor_del.gtf \
        --svcaller_T_DUP /path/to/svs/svcaller/tumor_dup.gtf \
        --svcaller_T_INV /path/to/svs/svcaller/tumor_inv.gtf \
        --svcaller_T_TRA /path/to/svs/svcaller/tumor_tra.gtf \
        --svcaller_N_DEL /path/to/svs/svcaller/normal_del.gtf \
        --svcaller_N_DUP /path/to/svs/svcaller/normal_dup.gtf \
        --svcaller_N_INV /path/to/svs/svcaller/normal_inv.gtf \
        --svcaller_N_TRA /path/to/svs/svcaller/normal_tra.gtf \
        --germline_mut_vcf /path/to/variants/germline-all.germline.vep.vcf \
        --somatic_mut_vcf /path/to/variants/somatic-all.somatic.vep.vcf \
        --output /path/to/qc/normal-tumor-frankenplot.html
```

### generate_evidence_bam.py

This scripts takes gridss variant file and gridss bam file and applies supplementary read filter and SAME_START_READS filter (i.e checks if reads have same chromosome number, start and end position which likely indicates duplicate reads; hence, not a true variant). 

```
generate_evidence_bam.py \
        --vcf /path/to/svs/gridss/tumor(or)normal-gridss.filtered.svannotated.vcf \
        --bam /path/to/svs/gridss/tumor(or)normal-gridss.sv.bam --filter-vcf \
        --output /path/to/svs/gridss/tumor(or)normal-gridss.evidence.bam
```

### gridss_svannotate.R

This script adds simple event type annotation (TRA, INV, INS, DEL, DUP) based on the breakend position and orientation to the end of INFO column.

```
gridss_svannotate.R -v /path/to/svs/gridss/tumor(or)normal-gridss.filtered.vcf \
        -o /path/to/svs/gridss/tumor(or)normal-gridss.filtered.svannotated.vcf
```

### run_msings.sh

This script takes bed file, reference genome, msi_intervals file, msings baseline file and bam file as input 
In first step, it uses samtools mpileup to count reads that has depth greater than 6 for each location mentioned in the interval_list file. The results are then passed to varscan which filters reads that has base quality less than 10. It then uses analyser.py script to calculate the count, percentage and distribution of indels within each range mentioned in the bed file. Finally, it uses count_msi_samples to count the MSI locations.

```
run_msings.sh -b msings.bed \
        -f /path/to/genome/human_g1k_v37_decoy.fasta \
        -i msings.intervals \
        " -n msings.baseline \
        " -o /path/to/msings-id/  tumor.bam
```

### vcfsorter.pl

This code takes fasta file and vcf file as input and sorts the vcf file based on the chromosome order present in fasta file. 

```
vcfsorter.pl /path/to/genome/human_g1k_v37_decoy.dict input.vcf
```

### jumble-run.R

This code is used to perform copy number analysis on duplicates removed bam file.
``` 
jumble-run.R -r panel_specific_twist.bed.reference.RDS -b id_nodups.bam
```
