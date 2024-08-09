
##################################################################
### For split_target rule, NEED TO USE ALL CHROMOSOMES in small-design
### Since, SD doesn't have all chromosomes. which will affect the
### realignment spliting and merging process
###
##################################################################
rule split_target_1:
    input:
        target = sd_targets['targets-bed']
    output:
        expand(outdir + "/bams/split_targets/target.sd.{chr}.bed", chr = all_chromosomes)
    params: 
        outdir = outdir
    shell:
        "mkdir -p {params.outdir}/bams/split_targets/ && "
        "for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do "
        " grep -w $chr {input.target} > {params.outdir}/bams/split_targets/target.sd.$chr.bed || true; "
        "done"



rule split_target_2:
    input:
        target = sd_targets[sd_capture_snv]['targets-bed']
    output:
        expand(outdir + "/bams/split_targets/target.snv.{chr}.bed", chr = all_chromosomes)
    params: 
        outdir = outdir
    shell:
        "mkdir -p {params.outdir}/bams/split_targets/ && "
        "for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do "
        " grep -w $chr {input.target} > {params.outdir}/bams/split_targets/target.snv.$chr.bed || true; "
        "done"


rule splitbam_umimapped_1:
    input:
        mapped = outdir + "/bams/{sample}_umimapped.bam",
        nochr = reference["no_chr"]
    output:
        expand(outdir + "/bams/split_targets/bam/{{sample}}_umimapped.{chr}.bam", chr = all_chromosomes),
        outdir + "/bams/split_targets/bam/{sample}_umimapped.nochr.bam"
    params:
        output_dir = outdir + "/bams/split_targets/bam/",
        all_chrom = all_chromosomes
    threads: 8
    shell:
        """
        prefix=$(basename {input.mapped} .bam)
        no_chr={params.output_dir}/${{prefix}}.nochr.bam
        all_chrom=({params.all_chrom})
        samtools view -@ {threads} -L {input.nochr} -o ${{no_chr}} {input.mapped}
        for chr in ${{all_chrom[@]}}; do
            samtools view -@ {threads} -b {input.mapped} ${{chr}} > {params.output_dir}/${{prefix}}.${{chr}}.bam
            samtools index {params.output_dir}/${{prefix}}.${{chr}}.bam
        done
        """


rule splitbam_umimapped_2:
    input:
        mapped = outdir + "/bams/{sample}_umimapped-2.bam",
        nochr = reference["no_chr"]
    output:
        expand(outdir + "/bams/split_targets/bam/{{sample}}_umimapped-2.{chr}.bam", chr = all_chromosomes),
        outdir + "/bams/split_targets/bam/{sample}_umimapped-2.nochr.bam"
    params:
        output_dir = outdir + "/bams/split_targets/bam/",
        all_chrom = all_chromosomes
    threads: 8
    shell:
        """
        prefix=$(basename {input.mapped} .bam)
        no_chr={params.output_dir}/${{prefix}}.nochr.bam
        all_chrom=({params.all_chrom})
        samtools view -@ {threads} -L {input.nochr} -o ${{no_chr}} {input.mapped}
        for chr in ${{all_chrom[@]}}; do
            samtools view -@ {threads} -b {input.mapped} ${{chr}} > {params.output_dir}/${{prefix}}.${{chr}}.bam
            samtools index {params.output_dir}/${{prefix}}.${{chr}}.bam
        done
        """

            

rule samtools_merge_realign_1:
    input:
        expand(outdir + "/bams/split_targets/bam/{{sample}}_realigned-1.{chr}.bam", chr = all_chromosomes),
        outdir + "/bams/split_targets/bam/{sample}_umimapped.nochr.bam"
    output:
        bam = outdir + "/bams/{sample}_realigned-1.bam",
        bai = outdir + "/bams/{sample}_realigned-1.bam.bai"
    shell:
        """
        InputBams=({input})
        bamfiles=${{InputBams[*]}}
        samtools merge -c -p {output.bam} ${{bamfiles}}
        samtools index {output.bam}
        rm $bamfiles
        """ 


rule samtools_merge_realign_2:
    input:
        expand(outdir + "/bams/split_targets/bam/{{sample}}_realigned-2.{chr}.bam", chr = all_chromosomes),
        outdir + "/bams/split_targets/bam/{sample}_umimapped-2.nochr.bam"
    output:
        bam = outdir + "/bams/{sample}_realigned-2.bam",
        bai = outdir + "/bams/{sample}_realigned-2.bam.bai"
    shell:
        """
        InputBams=({input})
        bamfiles=${{InputBams[*]}}
        samtools merge -c -p {output.bam} ${{bamfiles}}
        samtools index {output.bam}
        rm ${{bamfiles}}
        """   
