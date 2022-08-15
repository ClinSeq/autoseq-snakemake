

rule fastqc:
    input:
        libdir + "/{sample}/"
    output:
        directory(outdir + "/qc/fastqc/{sample}")
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
            "O={output.metrics} 2> {log} "


rule samtools_flagstat:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam"
    output:
        outdir + "/qc/samtools/{sample}-flagstats.json"
    threads: params['samtools']['flagstat']['threads']
    log:
        outdir + "/logs/samtools/samtools_flagstat_{sample}.log"
    shell:
        "samtools flagstat -@ {threads} -O json {input.bam} > {output} 2> {log} "


rule picard_wgsmetrics:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam",
        reference_genome = reference['reference_genome']
    output:
        metrics = outdir + "/qc/picard/{sample}.picard-wgsmetrics.txt"
    params:
        java_options = params['picard']['collectwgsmetrics']['java_options']
    threads: params['picard']['collectwgsmetrics']['threads']
    log:
        outdir + "/logs/picard/picard_collectwgsmetrics_{sample}.log"
    shell:
        "picard  {params.java_options}  " 
            "CollectWgsMetrics "
            "I={input.bam} "
            "R={input.reference_genome} "
            "O={output.metrics} 2> {log} "

