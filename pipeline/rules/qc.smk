import os
from functools import reduce


rule fastqc:
    input:
        libdir + "/{sample}/"
    output:
        directory(outdir + "/qc/fastqc/{sample}")
    params:
    threads: params['fastqc']['threads']
    log:
        outdir + "/logs/fastqc/fastqc_{sample}.log"
    run:
        fq_files = reduce(lambda r1, r2: r1 + r2, 
                          find_fastqs(wildcards.sample, libdir))
        
        shell("mkdir -p {output}")
        for fq in fq_files:
            shell("fastqc -o {output} --nogroup {fq}")


rule picard_collectinsertsize:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam"
    output:
        metrics = outdir + "/qc/picard/{sample}.picard-insertsize.txt"
    params:
        java_options = params['picard']['collectinsertsize']['java_options'],
        tmpdir = params['scratch'],
    threads: params['picard']['collectinsertsize']['threads']
    log:
        outdir + "/logs/picard/picard_insertsize_{sample}.log"
    shell:
        "picard  {params.java_options} -Djava.io.tmpdir={params.tmpdir} " 
            "CollectInsertSizeMetrics "
            "H=/dev/null "
            "I={input.bam} "
            "O={output.metrics} 2> {log}"


rule picard_collectoxog:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam",
        reference_genome = reference['reference_genome']
    output:
        metrics = outdir + "/qc/picard/{sample}.picard-oxog.txt"
    params:
        java_options = params['picard']['collectoxog']['java_options'],
        tmpdir = os.path.join(params['scratch'], 
                                "picard-oxog-{}".format(str(uuid.uuid4())))
    threads: params['picard']['collectoxog']['threads']
    log:
        outdir + "/logs/picard/picard_xoxg_{sample}.log"
    shell:
        "picard  {params.java_options} -Djava.io.tmpdir={params.tmpdir} " 
            "CollectOxoGMetrics "
            "I={input.bam} "
            "R={input.reference_genome} "
            "O={output.metrics} "


rule picard_collecthsmetrics:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam",
        reference_genome = reference['reference_genome'],
        target_region = lambda wildcards: get_targets(wildcards, reference,
                                                      'targets-interval_list-slopped20'),
        bait_regions = lambda wildcards: get_targets(wildcards, reference,
                                                      'targets-interval_list-slopped20')
    output:
        metrics = outdir + "/qc/picard/{sample}.picard-hsmetrics.txt"
    params:
        bait_name = lambda wildcards: get_target_name(wildcards),
        java_options = params['picard']['collecthsmetrics']['java_options'],
        tmpdir = os.path.join(params['scratch'], 
                                "picard-hsmetrics-{}".format(str(uuid.uuid4())))
    threads: params['picard']['collecthsmetrics']['threads']
    log:
        outdir + "/logs/picard/picard_hsmetrics_{sample}.log"
    shell:
        "picard {params.java_options}  -Djava.io.tmpdir={params.tmpdir} " 
            "CollectHsMetrics  "
            "I={input.bam}   "
            "R={input.reference_genome} "
            "O={output.metrics} "
            "TI={input.target_region}  "
            "BI={input.bait_regions} "
            "BAIT_SET_NAME={params.bait_name} "
            "METRIC_ACCUMULATION_LEVEL=LIBRARY "

