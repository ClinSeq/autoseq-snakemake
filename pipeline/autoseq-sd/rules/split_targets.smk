

rule split_target_1:
    input:
        target = reference['small_design'][capture_p2]['targets-bed']
    output:
        expand(outdir + "/bams/split_targets/target.p2.{chr}.bed", chr = p2_chromo)
    params: 
        outdir = outdir
    shell:
        "mkdir -p {params.outdir}/bams/split_targets/ && "
        "for chr in `cut -f 1 {input.target} | sort | uniq`; do "
        " grep -w $chr {input.target} > {params.outdir}/bams/split_targets/target.$chr.bed; "
        "done"



rule split_target_2:
    input:
        target = reference['small_design'][capture_s2]['targets-bed']
    output:
        expand(outdir + "/bams/split_targets/target.s2.{chr}.bed", chr = s2_chromo)
    params: 
        outdir = outdir
    shell:
        "mkdir -p {params.outdir}/bams/split_targets/ && "
        "for chr in `cut -f 1 {input.target} | sort | uniq`; do "
        " grep -w $chr {input.target} > {params.outdir}/bams/split_targets/target.$chr.bed; "
        "done"


rule splitbam_umimapped_1:
    input:
        mapped = outdir + "/bams/{sample}_umimapped.bam",
        nochr = reference["no_chr"]
    output:
        expand(outdir + "/bams/split_targets/bam/{{sample}}_umimapped.{chr}.bam", chr = p2_chromo),
        outdir + "/bams/split_targets/bam/{sample}_umimapped.nochr.bam"
    threads: 8
    run:
        output_dir = outdir + "/bams/split_targets/bam/"
        bam = input.mapped
        prefix = os.path.basename(bam).split('.bam')[0]
        no_chr = output_dir + "/{}.nochr.bam".format(prefix)
        cmd = "samtools view  -L {} -o {} {} ".format(input.nochr, no_chr, bam)
        shell(cmd)
        for chr in all_chromosomes:
            run_cmd = "samtools view -b {} {} ".format(bam, chr) + \
                        " > {}/{}.{}.bam && ".format(output_dir, prefix, chr) + \
                        " samtools index {}/{}.{}.bam ".format(output_dir, prefix, chr)
            shell(run_cmd)


rule splitbam_umimapped_2:
    input:
        mapped = outdir + "/bams/{sample}_umimapped-2.bam",
        nochr = reference["no_chr"]
    output:
        expand(outdir + "/bams/split_targets/bam/{{sample}}_umimapped-2.{chr}.bam", chr = s2_chromo),
        outdir + "/bams/split_targets/bam/{sample}_umimapped-2.nochr.bam"
    threads: 8
    run:
        output_dir = outdir + "/bams/split_targets/bam/"
        bam = input.mapped
        prefix = os.path.basename(bam).split('.bam')[0]
        no_chr = output_dir + "/{}.nochr.bam".format(prefix)
        cmd = "samtools view  -L {} -o {} {} ".format(input.nochr, no_chr, bam)
        shell(cmd)
        for chr in all_chromosomes:
            run_cmd = "samtools view -b {} {} ".format(bam, chr) + \
                        " > {}/{}.{}.bam && ".format(output_dir, prefix, chr) + \
                        " samtools index {}/{}.{}.bam ".format(output_dir, prefix, chr)
            shell(run_cmd)

            

rule samtools_merge_realign_1:
    input:
        expand(outdir + "/bams/split_targets/bam/{{sample}}_realigned-1.{chr}.bam", chr = p2_chromo),
        outdir + "/bams/split_targets/bam/{sample}_umimapped.nochr.bam"
    output:
        outdir + "/bams/{sample}_realigned-1.bam"
    run:
        bamfiles = " ".join(input)
        shell("samtools merge -c -p {output} {bamfiles}")
        shell("samtools index {output} ")
        shell("rm {bamfiles}")        


rule samtools_merge_realign_2:
    input:
        expand(outdir + "/bams/split_targets/bam/{{sample}}_realigned-2.{chr}.bam", chr = s2_chromo),
        outdir + "/bams/split_targets/bam/{sample}_umimapped-2.nochr.bam"
    output:
        outdir + "/bams/{sample}_realigned-2.bam"
    run:
        bamfiles = " ".join(input)
        shell("samtools merge -c -p {output} {bamfiles}")
        shell("samtools index {output} ")
        shell("rm {bamfiles}")
