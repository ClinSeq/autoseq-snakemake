

rule svcaller_run:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam",
        reference = reference["reference_genome"]
    output:
        gtf = outdir + "/svs/{sample}-{events}.gtf",
        bam = outdir + "/svs/{sample}-{events}.bam",
    params: 
        tmpdir = params['scratch']
    threads: params['svcaller']['threads']
    log:
        outdir + "/logs/svs/svcaller-{sample}-{events}.log"
    shell:
        "source activate svcallerenv  && "
        "svcaller run-all --tmp-dir {params.tmpdir} --event-type {wildcards.events} "
        " --fasta-filename {input.reference}  "
        " --filter-event-overlap "
        " --events-gtf {output.gtf} "
        " --events-bam {output.bam} {input.bam} && "
        "source deactivate"



# rule sveffect_predict:
#     input:
#         gtf = outdir + "/svs/{sample}-{events}.gtf",
#     output:
#     params:
#     threads:
#     log:
#     shell:
#         "source activate svcallerenv  && "
#         "sveffect make-bed --del-gtf {DEL.gtf} "
#         " --dup-gtf {DUP.gtf} "
#         " --inv-gtf {INV.gtf} " 
#         " --tra-gtf {TRA.gtf} "
#         " {combined.bed} &&  "
#         "sveffect predict --ts-regions {autoseq-genome/intervals/ts_regions.bed} "
#         " --ar-regions {autoseq-genome/intervals/ar_regions.bed} "
#         " --fusion-regions {autoseq-genome/intervals/fusion_regions.bed} "
#         " --effects-filename {effects.json} {combined.bed} && "
#         "source deactivate"