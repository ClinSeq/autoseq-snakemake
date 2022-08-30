

rule bwa_mem_alignment_normal:
    input:
        fq1 = os.path.join(nskewer_outdir, "{prefix}" + ns1),
        fq2 = os.path.join(nskewer_outdir, "{prefix}" + ns2),
        bwa_index = reference['bwaIndex']
    output:
        bamfile = outdir + "/bams/"+ normal_barcode + "/{prefix}.bam"
    wildcard_constraints:
        prefix = "|".join(nfq_prefix)
    params:
        readgroup = get_readgroup(normal_barcode),
        remove_duplicates = params['samblaster']['rm_dup'],
        tmpprefix = os.path.join(params['scratch'], 
                                "samtools-{}".format(str(uuid.uuid4())))
    threads: params['bwa']['threads']
    log:
        bwalog = outdir + "/logs/bwa_{prefix}.log",
        samblasterlog = outdir + "/logs/samblaster_{prefix}.log"
    shell:
        "bwa mem -M -v 1 -R  {params.readgroup} -t  {threads}"  
            " {input.bwa_index}  {input.fq1} {input.fq2}  2> {log.bwalog} "
            " | samblaster -M --addMateTags  {params.remove_duplicates} 2> {log.samblasterlog} "
            " | samtools view -Sb -u - | samtools sort  -T {params.tmpprefix} -@ {threads} "
            " -o  {output.bamfile}  -  && samtools index {output.bamfile} && "
            " rm {input.fq1} {input.fq2} "


rule bwa_mem_alignment_tumor:
    input:
        fq1 = os.path.join(tskewer_outdir, "{prefix}" + ts1),
        fq2 = os.path.join(tskewer_outdir, "{prefix}" + ts2),
        bwa_index = reference['bwaIndex']
    output:
        bamfile = outdir + "/bams/" + tumor_barcode + "/{prefix}.bam"
    wildcard_constraints:
        prefix = "|".join(tfq_prefix)
    params:
        readgroup = get_readgroup(tumor_barcode),
        remove_duplicates = params['samblaster']['rm_dup'],
        tmpprefix = os.path.join(params['scratch'], 
                                "samtools-{}".format(str(uuid.uuid4())))
    threads: params['bwa']['threads']
    log:
        bwalog = outdir + "/logs/bwa_{prefix}.log",
        samblasterlog = outdir + "/logs/samblaster_{prefix}.log"
    shell:
        "bwa mem -M -v 1 -R  {params.readgroup} -t  {threads}"  
            " {input.bwa_index}  {input.fq1} {input.fq2}  2> {log.bwalog} "
            " | samblaster -M --addMateTags  {params.remove_duplicates} 2> {log.samblasterlog} "
            " | samtools view -Sb -u - | samtools sort  -T {params.tmpprefix} -@ {threads} "
            " -o  {output.bamfile}  -  && samtools index {output.bamfile} && "
            " rm {input.fq1} {input.fq2} "


rule samtools_merge_normal:
    input:
        expand(outdir + "/bams/" + normal_barcode + "/{prefix}.bam", prefix = nfq_prefix)
    output:
        outdir + "/bams/{}.bam".format(normal_barcode)
    run:
        bamfiles = " ".join(input)
        shell("samtools merge -c -p {output} {bamfiles}")
        shell("samtools index {output} ")
        shell("rm {bamfiles}")


rule samtools_merge_tumor:
    input:
        expand(outdir + "/bams/" + tumor_barcode + "/{prefix}.bam", prefix = tfq_prefix)
    output:
        outdir + "/bams/{}.bam".format(tumor_barcode)
    run:
        bamfiles = " ".join(input)
        shell("samtools merge -c -p {output} {bamfiles}")
        shell("samtools index {output} ")
        shell("rm {bamfiles}")


rule samtools_splitbam:
    input:
        mapped = outdir + "/bams/{sample}.bam",
        nochr = reference["no_chr"]
    output:
        expand(outdir + "/bams/split_targets/bam/{{sample}}.{chr}.bam", chr = all_chromosomes),
        outdir + "/bams/split_targets/bam/{sample}.nochr.bam"
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


rule gatk3_targetcreator:
    input:
        bam = outdir + "/bams/split_targets/bam/{sample}.{chr}.bam",
        reference_genome = reference['reference_genome'],
        known_1kg = reference["1KG"],
        known_mills_gs = reference["Mills_and_1KG_gold_standard"],
    output:
        target_intervals = outdir + "/bams/split_targets/{sample}_{chr}.intervals"
    params:
        java_options = params['gatk3']['target_creator']['java_options'],
        extra = params['gatk3']['target_creator']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "realignerTC-{}".format(str(uuid.uuid4())))
    threads: params['gatk3']['target_creator']['threads']
    log:
        outdir + "/logs/gatk_realigner_targetcreator_{sample}_{chr}.log"
    shell:
        "source activate gatk_3 && "
        "gatk3 {params.java_options} -Djava.io.tmpdir={params.tmpdir} "
            " -T RealignerTargetCreator "
            " -R {input.reference_genome} "
            " -known {input.known_1kg} "
            " {params.extra} "
            " -known {input.known_mills_gs} "
            " -I {input.bam} "
            " -o {output.target_intervals} 2> {log} "


rule gatk3_indelrealigner:
    input:
        bam = outdir + "/bams/split_targets/bam/{sample}.{chr}.bam",
        reference_genome = reference['reference_genome'],
        known_1kg = reference["1KG"],
        known_mills_gs = reference["Mills_and_1KG_gold_standard"],
        target_intervals = outdir + "/bams/split_targets/{sample}_{chr}.intervals"
    output:
        bam = outdir + "/bams/split_targets/bam/{sample}_realigned.{chr}.bam",
    params:
        java_options = params['gatk3']['indel_realigner']['java_options'],
        extra = params['gatk3']['indel_realigner']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "indelrealigner-{}".format(str(uuid.uuid4())))
    threads: params['gatk3']['indel_realigner']['threads']
    log:
        outdir + "/logs/gatk_indel_realigner_{sample}_{chr}.log"
    shell:
        "source activate gatk_3 && "
        "gatk3 {params.java_options} -Djava.io.tmpdir={params.tmpdir} "
            " -T IndelRealigner  "
            " -R {input.reference_genome} "
            " -targetIntervals {input.target_intervals} "
            " -known {input.known_1kg} "
            " -known {input.known_mills_gs} "
            " {params.extra}"
            " -I {input.bam} "
            " -o {output.bam} 2> {log} && "
            " rm  {input.bam} "


rule samtools_merge_realign:
    input:
        expand(outdir + "/bams/split_targets/bam/{{sample}}_realigned.{chr}.bam", chr = all_chromosomes),
        outdir + "/bams/split_targets/bam/{sample}.nochr.bam"
    output:
        outdir + "/bams/{sample}_realigned.bam"
    run:
        bamfiles = " ".join(input)
        shell("samtools merge -c -p {output} {bamfiles}")
        shell("samtools index {output} ")
        shell("rm {bamfiles}") 


rule picard_markdups:
    input:
        bam = outdir + "/bams/{sample}_realigned.bam"
    output:
        bam = outdir + "/bams/{sample}_nodups.bam",
        metrics = outdir + "/bams/{sample}-picard-markdup.metrics.txt"
    params:
        rmdups = params['picard']['markdup']['rmdups'],
        java_options = params['picard']['markdup']['java_options'],
        extra = params['picard']['markdup']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "picard-markdups-{}".format(str(uuid.uuid4())))
    threads: params['picard']['markdup']['threads']
    log: outdir + "/logs/picard_markdups_{sample}.log"
    shell:
        "picard {params.java_options} -Djava.io.tmpdir={params.tmpdir} "
                " MarkDuplicates "
                " INPUT={input.bam} " 
                " METRICS_FILE={output.metrics} "
                " {params.extra} "
                " OUTPUT=/dev/stdout REMOVE_DUPLICATES={params.rmdups} "
                " | samtools sort -@ {threads} -T {params.tmpdir} -o {output.bam} 2> {log}"
                " && samtools index {output.bam} "


rule rm_interbamfiles:
    input:
        expand(outdir + "/bams/{sample}.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_realigned.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_nodups.bam", sample=all_clinseq_barcodes)
    output:
        outdir + "/bams/intermediate_bamfiles.removed"
    log:
        outdir + "/logs/remove_intermediate_{sample}.log".format(sample="_".join(all_clinseq_barcodes))
    run:
        del_bam = [bam for bam in input if 'nodups' not in bam]
        bamfiles = " ".join(del_bam)
        shell("rm {bamfiles} 2> {log} ")
        shell("touch {output} ")
