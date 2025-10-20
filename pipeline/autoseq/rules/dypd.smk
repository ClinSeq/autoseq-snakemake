

rule call_typeDPYD:
    input:
        bam = outdir + "/bams/{sample}_clipoverlap.bam"
    output:
        csv = outdir + "/variants/{sample}.typeDPYD.csv",
        json = outdir + "/variants/{sample}.typeDPYD.json",
    params:
        extra = params['dypd']['typeDPYD']['extra']
    threads: params['dypd']['typeDPYD']['threads']
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