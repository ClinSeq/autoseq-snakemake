

##### RNA sample barcode ################

rfq_prefix, rs1, rs2 = get_fqwildcards(rna_barcode, libdir)

rule star_alignment:
    input:
        fq1 = libdir + "/" + rna_barcode + "/{prefix}" + rs1,
        fq2 = libdir + "/" + rna_barcode + "/{prefix}" + rs2,
        genome_star_index = reference["star_index"],
    output:
        bam = outdir + "/bams/" + rna_barcode + "/{prefix}.bam" 
    wildcard_constraints:
        prefix = "|".join(rfq_prefix)
    params:
        output_prefix = outdir + "/bams/" + rna_barcode + "/{prefix}",
        sample_id = rna_barcode,
        args = ""
    threads: 16
    container: containers["autoseq-rnastar"]
    log:
        outdir + "/logs/star_alignment_{prefix}.log"
    shell:
        """
        STAR \\
            {params.args} \\
            --readFilesIn {input.fq1} {input.fq2} \\
            --genomeDir {input.genome_star_index} \\
            --runThreadN {threads} \\
            --readFilesCommand zcat \\
            --alignSJstitchMismatchNmax 5 -1 5 5 \\
            --alignSplicedMateMapLmin 35 \\
            --alignSplicedMateMapLminOverLmate 0.33 \\
            --chimJunctionOverhangMin 10 \\
            --chimOutType WithinBAM SoftClip \\
            --chimScoreDropMax 30 \\
            --chimScoreJunctionNonGTAG 0 \\
            --chimScoreMin 1 \\
            --chimScoreSeparation 1 \\
            --chimSegmentMin 10 \\
            --chimSegmentReadGapMax 3 \\
            --limitOutSJcollapsed 3000000 \\
            --outBAMcompression 0 \\
            --outFilterMatchNmin 35 \\
            --outFilterMatchNminOverLread 0.33 \\
            --outFilterMismatchNmax 3 \\
            --outFilterMultimapNmax 10 \\
            --outFilterScoreMinOverLread 0.33 \\
            --outSAMattributes All \\
            --outSAMattrRGline ID:{params.sample_id} SM:{params.sample_id} \\
            --outSAMtype BAM Unsorted \\
            --outSAMunmapped Within \\
            --outFileNamePrefix {params.output_prefix} \\
            --runRNGseed 0
        mv {params.output_prefix}Aligned.out.bam {output.bam}
        """


rule star_sort_bam:
    input:
        bam = outdir + "/bams/" + rna_barcode + "/{prefix}.bam"
    output:
        bam = outdir + "/bams/" + rna_barcode + "/{prefix}_sorted.bam"
    wildcard_constraints:
        prefix = "|".join(rfq_prefix)
    params:
        output_prefix = outdir + "/bams/" + rna_barcode + "/{prefix}"
    threads: 8
    log:
        outdir + "/logs/star_sort_bam_{prefix}.log"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.bam} 
        rm {input.bam}
        """


rule sambamba_merge:
    input:
        expand(outdir + "/bams/" + rna_barcode + "/{prefix}_sorted.bam", prefix = rfq_prefix)
    output:
        bam = outdir + "/bams/{}.bam".format(rna_barcode)
    params:
        extra = '',
        bamdir = outdir + "/bams/{}".format(rna_barcode)
    threads: 8
    container: containers['hmftools-redux']
    log:
        outdir + "/logs/sambamba_rna_merge_{}.log".format(rna_barcode)
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


rule gatk4_markduplicates:
    input:
        bam = outdir + "/bams/{}.bam".format(rna_barcode)
    output:
        bam = outdir + "/bams/{}.dedup.bam".format(rna_barcode),
        metrics = outdir + "/bams/{}.dedup.metrics".format(rna_barcode)
    params:
        bai = outdir + "/bams/{}.dedup.bai".format(rna_barcode),
        tmpdir = params['scratch'],
        java_options = '"-Xmx16G"',
    threads: 8
    log:
        outdir + "/logs/gatk4_markduplicates_{}.log".format(rna_barcode)
    shell:
        """
        gatk --java-options {params.java_options} MarkDuplicates \\
            --INPUT {input.bam} \\
            --OUTPUT {output.bam} \\
            --METRICS_FILE {output.metrics} \\
            --TMP_DIR {params.tmpdir} \\
            --CREATE_INDEX 
        mv {params.bai} {output.bam}.bai
        """