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
                          find_fastq(wildcards.sample, libdir))
        
        shell("mkdir -p {output}")
        for fq in fq_files:
            shell("fastqc -o {output} --nogroup {fq}")


rule picard_collectinsertsize:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam"
    output:
        metrics = outdir + "/qc/picard/{}.picard-insertsize.txt"
    params:
        java_options = params['picard']['collectinsertsize']['java_options'],
        tmpdir = os.path.join(params['scratch'], 
                                "picard-insertsize-{}".format(str(uuid.uuid4())))
    threads: params['picard']['collectinsersize']['threads']
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
        metrics = outdir + "/qc/picard/{}.picard-oxog.txt"
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


# rule picard_collecthsmetrics:
#     input:
#         bam = outdir + "/bams/{sample}_nodups.bam",
#         reference_genome = reference['reference_genome'],
#         target_region = 
#     output:
#     params:
#     threads:
#     log:
#     shell:
#         "picard {params.java_options}  -Djava.io.tmpdir={self.scratch} 
#             "CollectHsMetrics  "
#             "I={self.input}   "
#             "R={self.reference_sequence} "
#             "O={self.output_metrics} "
#             "TI={self.target_regions}  "
#             "BI={self.bait_regions} "
#             "BAIT_SET_NAME={self.bait_name} "
#             "METRIC_ACCUMULATION_LEVEL={self.accumulation_level} "