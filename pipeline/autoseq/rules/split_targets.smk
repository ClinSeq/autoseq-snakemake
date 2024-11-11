import os

capture_name = get_capture_name(CANCER_CAPTURE.capture_kit_id)

rule split_targets:
    input:
        target = reference['targets'][capture_name]['targets-bed-slopped20']
    output:
        expand(outdir + "/bams/split_targets/target.{chr}.bed", chr = all_chromosomes)
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
        expand(outdir + "/bams/split_targets/bam/{{sample}}_umimapped.{chr}.bam", chr = all_chromosomes),
        outdir + "/bams/split_targets/bam/{sample}_umimapped.nochr.bam"
    threads: 8
    run:
        output_dir = outdir + "/bams/split_targets/bam/"
        bam = input.mapped
        prefix = os.path.basename(bam).split('.bam')[0]
        no_chr = output_dir + "/{}.nochr.bam".format(prefix)
        cmd = "samtools view -@ {} -L {} -o {} {} ".format(threads, input.nochr, no_chr, bam)
        shell(cmd)
        for chr in all_chromosomes:
            run_cmd = "samtools view -@ {} -b {} {} ".format(threads, bam, chr) + \
                        " > {}/{}.{}.bam && ".format(output_dir, prefix, chr) + \
                        " samtools index {}/{}.{}.bam ".format(output_dir, prefix, chr)
            shell(run_cmd)


rule splitbam_umimapped_2:
    input:
        mapped = outdir + "/bams/{sample}_umimapped-2.bam",
        nochr = reference["no_chr"]
    output:
        expand(outdir + "/bams/split_targets/bam/{{sample}}_umimapped-2.{chr}.bam", chr = all_chromosomes),
        outdir + "/bams/split_targets/bam/{sample}_umimapped-2.nochr.bam"
    threads: 8
    run:
        output_dir = outdir + "/bams/split_targets/bam/"
        bam = input.mapped
        prefix = os.path.basename(bam).split('.bam')[0]
        no_chr = output_dir + "/{}.nochr.bam".format(prefix)
        cmd = "samtools view -@ {} -L {} -o {} {} ".format(threads, input.nochr, no_chr, bam)
        shell(cmd)
        for chr in all_chromosomes:
            run_cmd = "samtools view -@ {} -b {} {} ".format(threads, bam, chr) + \
                        " > {}/{}.{}.bam && ".format(output_dir, prefix, chr) + \
                        " samtools index {}/{}.{}.bam ".format(output_dir, prefix, chr)
            shell(run_cmd)

            

rule samtools_merge_realign_1:
    input:
        expand(outdir + "/bams/split_targets/bam/{{sample}}_realigned-1.{chr}.bam", chr = all_chromosomes),
        outdir + "/bams/split_targets/bam/{sample}_umimapped.nochr.bam"
    output:
        bam = outdir + "/bams/{sample}_realigned-1.bam",
        bai = outdir + "/bams/{sample}_realigned-1.bam.bai",
    run:
        bamfiles = " ".join(input)
        shell("samtools merge -c -p {output.bam} {bamfiles}")
        shell("samtools index {output.bam} ")
        shell("rm {bamfiles}")        


rule samtools_merge_realign_2:
    input:
        expand(outdir + "/bams/split_targets/bam/{{sample}}_realigned-2.{chr}.bam", chr = all_chromosomes),
        outdir + "/bams/split_targets/bam/{sample}_umimapped-2.nochr.bam"
    output:
        bam = outdir + "/bams/{sample}_realigned-2.bam",
        bai = outdir + "/bams/{sample}_realigned-2.bam.bai",
    run:
        bamfiles = " ".join(input)
        shell("samtools merge -c -p {output.bam} {bamfiles}")
        shell("samtools index {output.bam} ")
        shell("rm {bamfiles}")
