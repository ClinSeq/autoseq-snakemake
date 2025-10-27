

rule call_typeDPYD:
    input:
        bam = outdir + "/bams/{sample}_markdups.bam"
    output:
        csv = outdir + "/variants/{sample}.typeDPYD.csv",
        json = outdir + "/variants/{sample}.typeDPYD.json",
    params:
        extra = params['dpyd']['extra']
    threads: params['dpyd']['threads']
    container: containers['dpyd']
    log:
        outdir + "/logs/typeDPYD_{sample}.log"
    shell:
        """
        type_DPYD.R -b {input.bam} \\
            {params.extra}    \\ 
            -c {output.csv}   \\
            -j {output.json}  2> {log}

        """