import os
import uuid 


rule skewer_trim_pe:
    input:
        libdir + "/{sample}/"
    output:
        directory(outdir + "/fastqs/skewer/{sample}/")
    params:
        scratch = params['scratch']
    threads: 8
    log:
        outdir + "logs/skewer/skewer_{sample}.log"
    run:
        fq1_files, fq2_files = find_fastqs(wildcards.sample, libdir)
        fq1_abs = [normpath(x) for x in fq1_files]
        fq2_abs = [normpath(x) for x in fq2_files]
        pairs = [(fq1_abs[k], fq2_abs[k]) for k in range(len(fq1_abs))]
        
        for fq1, fq2 in pairs:
            tmpdir = os.path.join(params.scratch, "skewer-{}".format(str(uuid.uuid4())))
            prefix = "{}/skewer".format(tmpdir)

            pre_fq1 = prefix + "-trimmed-pair1.fastq.gz"
            pre_fq2 = prefix + "-trimmed-pair2.fastq.gz"

            out_fq1 = os.path.join(output[0], os.path.basename(fq1))
            out_fq2 = os.path.join(output[0], os.path.basename(fq2))

            shell(
                " mkdir {tmpdir} && "
                " skewer -z -t {threads} --quiet "
                " -o {prefix} "
                " {fq1} {fq2} && "
                " cp {pre_fq1} {out_fq1} && "
                " cp {pre_fq2} {out_fq2} && "
                " rm -rf {tmpdir} "
            )


rule cat_fastq:
    input:
        outdir + "/fastqs/skewer/{sample}/"
    output:
        fq1 = outdir + "/fastqs/{sample}_concatenated_1.fastq.gz",
        fq2 = outdir + "/fastqs/{sample}_concatenated_2.fastq.gz"
    run:
        libirary = outdir + "/fastqs/skewer/
        fq1_files, fq2_files = find_fastqs(wildcards.sample, libirary)
        fq1_flist = " ".join(fq1_files)
        fq2_flist = " ".join(fq2_files)

        shell(
            " cat {fq1_list} > {output.fq1} && "
            " cat {fq2_list} > {output.fq2} "
        )


rule bwa_mem_alignment:
    input:
        fq1 = outdir + "/fastqs/{sample}_concatenated_1.fastq.gz",
        fq2 = outdir + "/fastqs/{sample}_concatenated_2.fastq.gz",
        bwa_index = reference['bwaIndex']
    output:
        bamfile = outdir + "/bams/{sample}.bam"
    params:
        readgroup = lambda wildcards: get_readgroup(wildcards),
        remove_duplicates = params['samblaster']['rm_dup'],
        tmpprefix = os.path.join(params['scratch'], 
                                "samtools-{}".format(str(uuid.uuid4())))
    threads: params['bwa']['threads']
    log:
        bwalog = outdir + "/logs/bwa_{sample}.log",
        samblasterlog = outdir + "/logs/samblaster_{sample}.log"
    shell:
        "bwa mem -M -v 1 -R  {params.readgroup} -t  {threads}"  
            " {input.bwa_index}  {input.fq1} {input.fq2}  2> {log.bwalog} "
            " | samblaster -M --addMateTags  {params.remove_duplicates} 2> {log.samblasterlog} "
            " | samtools view -Sb -u - | samtools sort  -T {params.tmpprefix} -@ {threads} "
            " -o  {output.bamfile}  -  && samtools index {output.bamfile} "


rule gatk3_targetcreator:
    input:
        bam = outdir + "/bams/{sample}.bam",
        reference_genome = reference['reference_genome'],
        target_region = lambda wildcards: get_targets(wildcards, reference),
        known_1kg = reference["1KG"],
        known_mills_gs = reference["Mills_and_1KG_gold_standard"],
    output:
        target_intervals = outdir + "/bams/{sample}.intervals"
    params:
        java_options = params['gatk3']['realigner_target']['java_options'],
        jarfile = params['gatk3']['jarfile'],
        extra = params['gatk3']['realigner_target']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "realignerTC-{}".format(str(uuid.uuid4())))
    threads: params['gatk3']['realigner_target']['threads']
    log:
        outdir + "/logs/gatk_realigner_targetcreator_{sample}.log"
    shell:
        "source activate gatk_3 && "
        "gatk3 {params.java_options} -Djava.io.tmpdir={params.tmpdir} "
            " -T RealignerTargetCreator "
            " -R {input.reference_genome} "
            " -known {input.known_1kg} "
            " {params.extra} "
            " -L {input.target_region} "
            " -known {input.known_mills_gs} "
            " -I {input.bam} "
            " -o {output.target_intervals} 2> {log} "


rule gatk3_indelrealigner:
    input:
        bam = outdir + "/bams/{sample}.bam",
        reference_genome = reference['reference_genome'],
        target_region = lambda wildcards: get_targets(wildcards, reference),
        known_1kg = reference["1KG"],
        known_mills_gs = reference["Mills_and_1KG_gold_standard"],
        target_intervals = outdir + "/bams/{sample}.intervals"
    output:
        bam = outdir + "/bams/{sample}_realigned.bam",
    params:
        jarfile = params['gatk3']['jarfile'],
        java_options = params['gatk3']['indel_realigner']['java_options'],
        extra = params['gatk3']['indel_realigner']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "indelrealigner-{}".format(str(uuid.uuid4())))
    threads: params['gatk3']['indel_realigner']['threads']
    log:
        outdir + "/logs/gatk_indel_realigner_{sample}.log"
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
            " -o {output.bam} 2> {log} "


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
                " | samtools sort -@ {threads} -T {params.tmpdir} -o {output.bam} "
                " && samtools index {output.bam} 2> {log}"