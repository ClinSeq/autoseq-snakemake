# Autoseq internal scripts 


### `create_contest_vcfs.py`

Create the contamination test VCF input file from the designated "genotype" and "eval" target bed region files and population allele frequency VCF file.

### `vcf_add_sample.py` 

Add DP, RO and AO for a new sample from a BAM file. Filters variants that are not simple substitutions. Optionally filter homozygous variants from the output file (--filter_hom)

```
usage: vcf_add_sample.py --samplename variants.vcf(.gz) aligned.bam > new.vcf
```
### `contest_to_contam_caveat.py`

Generate a call on contamination from contest output, and output it as a JSON file.

Example output contents: `{"CALL": "OK"}`


### `generate_allelic_fraction_bedGraph.py`

For each variants in the VCF file, generates allelic fraction. 

### `generateIGVnavInput.py`

Script to create list of SNVs from merged somatic or germline vcfs with oncogenicity annoation for IGVNav.

### `generateIGVnavInput_SV.py`

Script to create mut files from different structural variant caller and combine all mut files with gene annotation.

### `generate_symlinks.py`

Generates symlinks required for IGVnav inputs - curations.

### `liqbioCNA_Interactive_plots.R`

To create static and interactive franken plot.

### `PureCN.R`

Customized script to estimate tumor purity, copy number, and loss of heterozygosity (LOH), and classifies single nucleotide variants (SNVs) by somatic status and clonality. 

### `QC_overview.R`

This is a script creating plots showing the distribution of some QC values, and where the current normal and tumor sample lies within those distributions. QC metrics shown are read count, duplication rate, coverage, fold enrichment, on-bait rate, insert size & contamination

### `frankenscript.R`

### `generate_evidence_bam.py`

This scripts takes gridss variant file and gridss bam file and applies supplementary read filter and SAME_START_READS filter (i.e checks if reads have same chromosome number, start and end position which likely indicates duplicate reads; hence, not a true variant). 

### `gridss_svannotate.R`

This script adds simple event type annotation (TRA, INV, INS, DEL, DUP) based on the breakend position and orientation to the end of INFO column.

### `run_msings.sh`

This script takes bed file, reference genome, msi_intervals file, msings baseline file and bam file as input 
In first step, it uses samtools mpileup to count reads that has depth greater than 6 for each location mentioned in the interval_list file. The results are then passed to varscan which filters reads that has base quality less than 10. It then uses analyser.py script to calculate the count, percentage and distribution of indels within each range mentioned in the bed file. Finally, it uses count_msi_samples to count the MSI locations.

### `vcfsorter.pl`

This code takes fasta file and vcf file as input and sorts the vcf file based on the chromosome order present in fasta file. 

### `jumble-run.R`
