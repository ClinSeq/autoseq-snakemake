
## read chrsize file for target regions in indelrealigner

chrsizes = dict()
fchrsizes = reference['chrsizes'] 

with open(fchrsizes, 'r') as fh:
    for line in fh.readlines():
        data = line.strip().split()
        chrsizes[data[0]] = "-".join(['1', str(data[1])])

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
    threads: 8
    run:
        bamfiles = " ".join(input)
        shell("samtools merge -@ {threads} -c -p {output} {bamfiles}")
        shell("samtools index {output} ")
        shell("rm {bamfiles}")


rule samtools_merge_tumor:
    input:
        expand(outdir + "/bams/" + tumor_barcode + "/{prefix}.bam", prefix = tfq_prefix)
    output:
        outdir + "/bams/{}.bam".format(tumor_barcode)
    threads: 8
    run:
        bamfiles = " ".join(input)
        shell("samtools merge -@ {threads} -c -p {output} {bamfiles}")
        shell("samtools index {output} ")
        shell("rm {bamfiles}")


rule picard_markdups:
    input:
        bam = outdir + "/bams/{sample}.bam"
    output:
        bam = outdir + "/bams/{sample}_nodups.bam",
        metrics = outdir + "/qc/picard/{sample}-picard-markdup.metrics.txt"
    params:
        rmdups = params['picard']['markdup']['rmdups'],
        java_options = params['picard']['markdup']['java_options'],
        extra = params['picard']['markdup']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "picard-markdups-{}".format(str(uuid.uuid4())))
    threads: 8
    log: outdir + "/logs/picard_markdups_{sample}.log"
    shell:
        "picard {params.java_options} -Djava.io.tmpdir={params.tmpdir} "
            " MarkDuplicates "
            " INPUT={input.bam} " 
            " METRICS_FILE={output.metrics} "
            " {params.extra} "
            " OUTPUT=/dev/stdout REMOVE_DUPLICATES={params.rmdups} "
            " | samtools sort -m 2G -@ {threads} -T {params.tmpdir} -o {output.bam} 2> {log}"
            " && samtools index {output.bam} "
            " && rm -rf {params.tmpdir} "


rule rm_interbamfiles:
    input:
        expand(outdir + "/bams/{sample}.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_nodups.bam", sample=all_clinseq_barcodes)
        # expand(outdir + "/bams/{sample}_realigned.bam", sample=all_clinseq_barcodes),
    output:
        outdir + "/bams/intermediate_bamfiles.removed"
    log:
        outdir + "/logs/remove_intermediate_{sample}.log".format(sample="_".join(all_clinseq_barcodes))
    run:
        del_bam = [bam for bam in input if 'nodups' not in bam]
        bamfiles = " ".join(del_bam)
        shell("rm {bamfiles} 2> {log} ")
        shell("touch {output} ")
