
## read chrsize file for target regions in indelrealigner

chrsizes = dict()
fchrsizes = reference['chrsizes'] 

with open(fchrsizes, 'r') as fh:
    for line in fh.readlines():
        data = line.strip().split()
        chrsizes[data[0]] = "-".join(['1', str(data[1])])

rule bwa_mem2_alignment_normal:
    input:
        fq1 = os.path.join(nskewer_outdir, "{n}.{prefix}" + ns1),
        fq2 = os.path.join(nskewer_outdir, "{n}.{prefix}" + ns2),
        bwa_index = reference['bwa-mem2-idx']
    output:
        bamfile = outdir + "/bams/"+ normal_barcode + "/{n}.{prefix}.bam"
    wildcard_constraints:
        prefix = "|".join(nfq_prefix),
        n = "|".join(fastp_split[0:nsplit])
    params:
        readgroup = get_readgroup(normal_barcode)
    threads: 12
    container: containers['mulled-v2']
    log:
        outdir + "/logs/bwa_{n}.{prefix}.log",
    shell:
        """
        bwa-mem2 mem \\
            -Y \\
            -K 100000000 \\
            -R {params.readgroup} \\
            -t {threads} \\
            {input.bwa_index} \\
            {input.fq1} \\
            {input.fq2} | \\
            sambamba view \\
                --sam-input \\
                --format bam \\
                --compression-level 0 \\
                --nthreads {threads} \\
                /dev/stdin | \\
            sambamba sort \\
                --nthreads {threads} \\
                --out {output.bamfile} \\
                /dev/stdin 2> {log}
        rm {input.fq1} {input.fq2}
        """


rule bwa_mem2_alignment_tumor:
    input:
        fq1 = os.path.join(tskewer_outdir, "{n}.{prefix}" + ts1),
        fq2 = os.path.join(tskewer_outdir, "{n}.{prefix}" + ts2),
        bwa_index = reference['bwa-mem2-idx']
    output:
        bamfile = outdir + "/bams/" + tumor_barcode + "/{n}.{prefix}.bam"
    wildcard_constraints:
        prefix = "|".join(tfq_prefix),
        n = "|".join(fastp_split)
    params:
        readgroup = get_readgroup(tumor_barcode),
    threads: 16
    container: containers['mulled-v2']
    log:
        outdir + "/logs/bwa_{n}.{prefix}.log"
    shell:
        """
        bwa-mem2 mem \\
            -Y \\
            -K 100000000 \\
            -R {params.readgroup} \\
            -t {threads} \\
            {input.bwa_index} \\
            {input.fq1} \\
            {input.fq2} | \\
            sambamba view \\
                --sam-input \\
                --format bam \\
                --compression-level 0 \\
                --nthreads {threads} \\
                /dev/stdin | \\
            sambamba sort \\
                -m 1G  \\
                --nthreads {threads} \\
                --out {output.bamfile} \\
                /dev/stdin 2> {log}
        rm {input.fq1} {input.fq2}
        """


rule sambamba_merge_normal:
    input:
        expand(outdir + "/bams/" + normal_barcode + "/{n}.{prefix}.bam", n = fastp_split[0:nsplit], prefix = nfq_prefix)
    output:
        bam = outdir + "/bams/{}.bam".format(normal_barcode)
    params:
        extra = '',
        bamdir = outdir + "/bams/{}".format(normal_barcode)
    threads: 8
    container: containers['mulled-v2']
    log:
        outdir + "/logs/sambamba_merge_{}.log".format(normal_barcode)
    shell:
        """
        InputBams=({input})
        bamfiles=${{InputBams[*]}}
        sambamba merge \\
            --nthreads {threads} \\
            {output.bam} \\
            ${{bamfiles}} 2> {log}
        sambamba index {output.bam}
        rm ${{bamfiles}} && rm -rf {params.bamdir}
        """


rule sambamba_merge_tumor:
    input:
        expand(outdir + "/bams/" + tumor_barcode + "/{n}.{prefix}.bam", n = fastp_split, prefix = tfq_prefix)
    output:
        bam = outdir + "/bams/{}.bam".format(tumor_barcode)
    params:
        extra = '',
        bamdir = outdir + "/bams/{}".format(tumor_barcode)
    threads: 8
    container: containers['mulled-v2']
    log: 
        outdir + "/logs/sambamba_merge_{}.log".format(tumor_barcode)
    shell:
        """
        InputBams=({input})
        bamfiles=${{InputBams[*]}}
        sambamba merge \\
            --nthreads {threads} \\
            {output.bam} \\
            ${{bamfiles}} 2> {log}
        sambamba index {output.bam}
        rm ${{bamfiles}} && rm -rf {params.bamdir}
        """


rule hmftools_redux:
    input:
        bam = outdir + "/bams/{sample}.bam",
        unmap_regions = reference['unmap_regions'],
        ref_genome = reference['reference_genome']
    output:
        bam = outdir + "/bams/{sample}_markdups.bam"
    params:
        sample_id = "{sample}",
        ref_genome_ver = "37",
        extra = ""
    threads: 8
    container: containers['hmftools-redux']
    log: outdir + "/logs/hmftools_markdups_{sample}.log"
    shell:
        """
        java_mem="-Xmx120g"
        if [[ "{params.sample_id}" == *-N-* ]]; then
            java_mem="-Xmx32g"
        fi
        redux \\
            $java_mem \\
            -bamtool $(type -p sambamba) \\
            -sample {params.sample_id} \\
            -input_bam {input.bam} \\
            -form_consensus \\
            -unmap_regions {input.unmap_regions} \\
            -ref_genome {input.ref_genome} \\
            -ref_genome_version {params.ref_genome_ver} \\
            -write_stats \\
            -threads {threads} \\
            -output_bam {output.bam} 2> {log}
        sambamba index {output.bam}
        """


rule rm_interbamfiles:
    input:
        expand(outdir + "/bams/{sample}.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_markdups.bam", sample=all_clinseq_barcodes)
    output:
        outdir + "/bams/intermediate_bamfiles.removed"
    log:
        outdir + "/logs/remove_intermediate_{sample}.log".format(sample="_".join(all_clinseq_barcodes))
    run:
        del_bam = [bam for bam in input if 'markdups' not in bam]
        bamfiles = " ".join(del_bam)
        shell("rm {bamfiles} 2> {log} ")
        shell("touch {output} ")
