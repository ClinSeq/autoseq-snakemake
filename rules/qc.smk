import os

rule fastqc:
    input:
        libdir + "/{sample}/"
    output:
        directory(outdir + "/qc/fastqc/{sample}")
    params:
    threads: params['fastqc']['threads']
    log:
        outdir + "logs/fastqc/fastqc_{sample}.log"
    run:
        fq_files = reduce(lambda r1, r2: r1 + r2, 
                          find_fastq(wildcards.sample, libdir))
        
        for fq in fq_files:
            shell(
                "fastqc -o {output} --nogroup {fq}"
                )
