
rule fgbio_fastqtobam:
    input:
        fq1 = outdir + "/fastqs/{sample}_concatenated_1.fastq.gz",
        fq2 = outdir + "/fastqs/{sample}_concatenated_2.fastq.gz"
    output:
        bam = outdir + "/bams/{sample}_unmapped.bam"
    params:
        sample = lambda wildcards: compose_sample_str(extract_unique_capture(wildcards.sample)),
        library = lambda wildcards: parse_prep_id(wildcards.sample),
        java_options = params['fgbio']['fastqtobam']['java_options'],
        tmpdir = os.path.join(params['scratch'], 
                    "fgbio_fastqtobam-{}".format(str(uuid.uuid4())))
    threads: params['fgbio']['fastqtobam']['threads']
    log: outdir + "/logs/fgbio/fastqtobam_{sample}.log"
    shell:
        "fgbio {params.java_options} --tmp-dir {params.tmpdir} FastqToBam "
            " -i {input.fq1}  {input.fq2} " 
            " -o {output.bam} "
            " --sample {params.sample} "
            " --library {params.library} "
            " -r 3M2S+T 3M2S+T "
            " -s true 2> {log} "
            " && rm -rf {params.tmpdir}"


rule bwa_umialignment:
    input:
        bam = outdir + "/bams/{sample}_unmapped.bam",
        reference_genome = reference['bwaIndex']
    output:
        bam = outdir + "/bams/{sample}_umimapped.bam",
        bai = outdir + "/bams/{sample}_umimapped.bai"
    params:
        java_options = params['picard']['merge_bam']['java_options'],
        tmpdir = os.path.join(params['scratch'], 
                    "bwa-umialignment-{}".format(str(uuid.uuid4())))
    threads: params['bwa']['threads']
    log: outdir + "/logs/fgbio/bwa_umialignment_{sample}.log"
    shell:
        "picard SamToFastq I={input.bam} F=/dev/stdout INTERLEAVE=true TMP_DIR={params.tmpdir} "
            " | bwa mem -p -t {threads} {input.reference_genome} /dev/stdin " 
            " | picard -Djava.io.tmpdir={params.tmpdir}  {params.java_options} MergeBamAlignment "
            " UNMAPPED={input.bam} "
            " ALIGNED=/dev/stdin "
            " O={output.bam}"
            " R={input.reference_genome} "
            " SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true "
            " MAX_GAPS=-1 ORIENTATIONS=FR CREATE_INDEX=true "
            " TMP_DIR={params.tmpdir} && rm -rf {params.tmpdir} 2> {log} "


rule gatk3_targetcreator_umi_1:
    input:
        bam = outdir + "/bams/split_targets/bam/{sample}_umimapped.{chr}.bam",
        reference_genome = reference['reference_genome'],
        target_region = outdir + "/bams/split_targets/target.sd.{chr}.bed",
        known_1kg = reference["1KG"],
        known_mills_gs = reference["Mills_and_1KG_gold_standard"]
    output:
        target_intervals = outdir + "/bams/split_targets/{sample}_umi_{chr}.intervals"
    params:
        java_options = params['gatk3']['target_creator']['java_options'],
        extra = params['gatk3']['target_creator']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "realignerTC-{}".format(str(uuid.uuid4())))
    threads: params['gatk3']['target_creator']['threads']
    container: containers['gatk3']
    log:
        outdir + "/logs/gatk_realigner_targetcreator_umi_1_{sample}_{chr}.log"
    shell:
        "source activate gatk_3 && "
        "gatk3 {params.java_options} -Djava.io.tmpdir={params.tmpdir} "
            " -T RealignerTargetCreator "
            " -R {input.reference_genome} "
            " -known {input.known_1kg} "
            " {params.extra} "
            " -L {input.target_region} "
            " -known {input.known_mills_gs} "
            " -I {input.bam} "
            " -o {output.target_intervals} 2> {log} "
            " && rm -rf {params.tmpdir} "


rule gatk3_indelrealigner_umi_1:
    input:
        bam = outdir + "/bams/split_targets/bam/{sample}_umimapped.{chr}.bam",
        reference_genome = reference['reference_genome'],
        target_region = outdir + "/bams/split_targets/target.sd.{chr}.bed",
        known_1kg = reference["1KG"],
        known_mills_gs = reference["Mills_and_1KG_gold_standard"],
        target_intervals = outdir + "/bams/split_targets/{sample}_umi_{chr}.intervals"
    output:
        bam = outdir + "/bams/split_targets/bam/{sample}_realigned-1.{chr}.bam"
    params:
        java_options = params['gatk3']['indel_realigner']['java_options'],
        extra = params['gatk3']['indel_realigner']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "indelrealigner-{}".format(str(uuid.uuid4())))
    threads: params['gatk3']['indel_realigner']['threads']
    container: containers['gatk3']
    log:
        outdir + "/logs/gatk_indel_realigner_umi_1_{sample}_{chr}.log"
    shell:
        "source activate gatk_3 && "
        "gatk3 {params.java_options} -Djava.io.tmpdir={params.tmpdir} "
            " -T IndelRealigner  "
            " -R {input.reference_genome} "
            " -targetIntervals {input.target_intervals} "
            " -known {input.known_1kg} "
            " -known {input.known_mills_gs} "
            " {params.extra}"
            " -I {input.bam} "
            " -o {output.bam} 2> {log} && "
            "rm  {input.bam} "
            " && rm -rf {params.tmpdir} "


rule fgbio_groupreadsbyumi:
    input:
        bam = outdir + "/bams/{sample}_realigned-1.bam"
    output:
        bam = outdir + "/bams/{sample}_groupedbyumi.bam",
        histogram = outdir + "/bams/{sample}-groupedbyumi.bam.fs.txt"
    params:
        java_options = params['fgbio']['groupreadsbyumi']['java_options'],
        tmpdir = os.path.join(params['scratch'], 
                                "groupreadsbyumi-{}".format(str(uuid.uuid4())))
    threads: params['fgbio']['groupreadsbyumi']['threads']
    log: outdir + "/logs/fgbio_groupreadsbyumi_{sample}.log"
    shell:
        " fgbio {params.java_options} --tmp-dir {params.tmpdir} GroupReadsByUmi  "
              " -i {input.bam} "
              " -o {output.bam} "
              " --strategy paired "
              " --family-size-histogram {output.histogram} 2> {log} "
              " && rm -rf {params.tmpdir} "



rule fgbio_callduplexconsensus:
    input:
        bam = outdir + "/bams/{sample}_groupedbyumi.bam"
    output:
        bam = outdir + "/bams/{sample}_consensus.bam"
    params:
        java_options = params['fgbio']['callduplexconsensus']['java_options'],
        extra = params['fgbio']['callduplexconsensus']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "groupreadsbyumi-{}".format(str(uuid.uuid4())))
    threads: params['fgbio']['callduplexconsensus']['threads']
    log: outdir + "/logs/fgbio_callduplexconsensus_{sample}.log"
    shell:
        "fgbio {params.java_options} --tmp-dir {params.tmpdir} CallDuplexConsensusReads  "
             " -i {input.bam}  "
             " -o  {output.bam} "
             " --threads {threads} " 
             " {params.extra} "
             " && rm -rf {params.tmpdir}"


rule bwa_umialignment_2:
    input:
        bam = outdir + "/bams/{sample}_consensus.bam",
        reference_genome = reference['bwaIndex']
    output:
        bam = outdir + "/bams/{sample}_umimapped-2.bam",
        bai = outdir + "/bams/{sample}_umimapped-2.bai"
    params:
        java_options = params['picard']['merge_bam']['java_options'],
        tmpdir = os.path.join(params['scratch'], 
                    "bwa-umialignment-{}".format(str(uuid.uuid4())))
    threads: params['bwa']['threads']
    log: outdir + "/logs/fgbio/bwa_umialignment_2_{sample}.log"
    shell:
        "picard SamToFastq I={input.bam} F=/dev/stdout INTERLEAVE=true TMP_DIR={params.tmpdir} "
            " | bwa mem -p -t {threads} {input.reference_genome} /dev/stdin " 
            " | picard -Djava.io.tmpdir={params.tmpdir}  {params.java_options} MergeBamAlignment "
            " UNMAPPED={input.bam} "
            " ALIGNED=/dev/stdin "
            " O={output.bam}"
            " R={input.reference_genome} "
            " SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true "
            " MAX_GAPS=-1 ORIENTATIONS=FR CREATE_INDEX=true "
            " TMP_DIR={params.tmpdir} && rm -rf {params.tmpdir} "


rule gatk3_targetcreator_umi_2:
    input:
        bam = outdir + "/bams/split_targets/bam/{sample}_umimapped-2.{chr}.bam",
        reference_genome = reference['reference_genome'],
        target_region = outdir + "/bams/split_targets/target.snv.{chr}.bed",
        known_1kg = reference["1KG"],
        known_mills_gs = reference["Mills_and_1KG_gold_standard"]
    output:
        target_intervals = outdir + "/bams/split_targets/bam/{sample}_consensus_{chr}.intervals"
    params:
        java_options = params['gatk3']['target_creator']['java_options'],
        extra = params['gatk3']['target_creator']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "realignerTC-{}".format(str(uuid.uuid4())))
    threads: params['gatk3']['target_creator']['threads']
    container: containers['gatk3']
    log:
        outdir + "/logs/gatk_realigner_targetcreator_umi_2_{sample}.{chr}.log"
    shell:
        "source activate gatk_3 && "
        "gatk3 {params.java_options} -Djava.io.tmpdir={params.tmpdir} "
            " -T RealignerTargetCreator "
            " -R {input.reference_genome} "
            " -known {input.known_1kg} "
            " {params.extra} "
            " -L {input.target_region} "
            " -known {input.known_mills_gs} "
            " -I {input.bam} "
            " -o {output.target_intervals} 2> {log} "
            " && rm -rf {params.tmpdir} "


rule gatk3_indelrealigner_umi_2:
    input:
        bam = outdir + "/bams/split_targets/bam/{sample}_umimapped-2.{chr}.bam",
        reference_genome = reference['reference_genome'],
        target_region = outdir + "/bams/split_targets/target.snv.{chr}.bed",
        known_1kg = reference["1KG"],
        known_mills_gs = reference["Mills_and_1KG_gold_standard"],
        target_intervals = outdir + "/bams/split_targets/bam/{sample}_consensus_{chr}.intervals"
    output:
        bam = outdir + "/bams/split_targets/bam/{sample}_realigned-2.{chr}.bam"
    params:
        java_options = params['gatk3']['indel_realigner']['java_options'],
        extra = params['gatk3']['indel_realigner']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "indelrealigner-{}".format(str(uuid.uuid4())))
    threads: params['gatk3']['indel_realigner']['threads']
    container: containers['gatk3']
    log:
        outdir + "/logs/gatk_indel_realigner_umi_2_{sample}.{chr}.log"
    shell:
        "source activate gatk_3 && "
        "gatk3 {params.java_options} -Djava.io.tmpdir={params.tmpdir} "
            " -T IndelRealigner  "
            " -R {input.reference_genome} "
            " -targetIntervals {input.target_intervals} "
            " -known {input.known_1kg} "
            " -known {input.known_mills_gs} "
            " {params.extra}"
            " -I {input.bam} "
            " -o {output.bam} 2> {log} && "
            "rm  {input.bam} "


rule fgbio_filterconsensus:
    input:
        bam = outdir + "/bams/{sample}_realigned-2.bam",
        reference_genome = reference['reference_genome']
    output:
        bam = outdir + "/bams/{sample}_consensus_filtered.bam",
        bai = outdir + "/bams/{sample}_consensus_filtered.bai"
    params:
        java_options = params['fgbio']['filterconsensus']['java_options'],
        error_rate = params['fgbio']['filterconsensus']['error_rate'],
        base_quality = params['fgbio']['filterconsensus']['base_quality'],
        extra = " --min-reads 1 1 0 --reverse-per-base-tags true --require-single-strand-agreement true ",
        tmpdir = os.path.join(params['scratch'], 
                                "fgbio-filterconsensus-{}".format(str(uuid.uuid4())))
    threads: params['fgbio']['filterconsensus']['threads']
    log: 
        outdir + "/logs/fgbio_filter_consensus_{sample}.log"
    shell:
        "fgbio {params.java_options} --tmp-dir {params.tmpdir} FilterConsensusReads "
            " -i {input.bam}  "
            " -o {output.bam} "
            " --ref {input.reference_genome} "
            " {params.extra} "
            " {params.error_rate} "
            " {params.base_quality} "
            " && rm -rf {params.tmpdir} "


rule fgbio_clipbam:
    input:
        bam = outdir + "/bams/{sample}_consensus_filtered.bam",
        reference_genome = reference['reference_genome']
    output:
        bam = outdir + "/bams/{sample}_clipoverlap.bam",
        bai = outdir + "/bams/{sample}_clipoverlap.bai",
        metrics_txt = outdir + "/bams/{sample}_clipoverlap_metrix.txt"
    params:
        java_options = params['fgbio']['clipbam']['java_options'],
        tmpdir = os.path.join(params['scratch'], 
                                "fgbio-clipbam-{}".format(str(uuid.uuid4())))
    threads: params['fgbio']['clipbam']['threads']
    log:
        outdir + "/logs/fgbio_clipbam_{sample}.log"
    shell:
        "fgbio {params.java_options} --tmp-dir {params.tmpdir} ClipBam "
            " -i  {input.bam}  "
            " -o {output.bam}  "
            " -m  {output.metrics_txt} "
            " --ref {input.reference_genome} "
            " --clip-overlapping-reads true "
            " && rm -rf {params.tmpdir} "


rule picard_markdups:
    input:
        bam = outdir + "/bams/{sample}_realigned-1.bam"
    output:
        bam = outdir + "/bams/{sample}_nodups.bam",
        metrics = outdir + "/qc/picard/{sample}-picard-markdup.metrics.txt"
    params:
        rmdups = params['picard']['markdup']['rmdups'],
        java_options = params['picard']['markdup']['java_options'],
        extra = params['picard']['markdup']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "picard-markdups-{}".format(str(uuid.uuid4())))
    threads: params['picard']['markdup']['threads']
    log: outdir + "/logs/picard_markdups_{sample}.log"
    shell:
        "picard {params.java_options} -Djava.io.tmpdir={params.tmpdir} "
                " MarkDuplicates "
                " INPUT={input.bam} " 
                " METRICS_FILE={output.metrics} "
                " {params.extra} "
                " OUTPUT=/dev/stdout REMOVE_DUPLICATES={params.rmdups} "
                " | samtools sort -@ {threads} -T {params.tmpdir} -o {output.bam} "
                " && samtools index {output.bam} 2> {log}"
                " && rm -rf {params.tmpdir} "


rule rm_interbamfiles:
    input:
        expand(outdir + "/bams/{sample}_unmapped.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_umimapped.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_umimapped.bai", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_realigned-1.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_realigned-1.bam.bai", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_groupedbyumi.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_consensus.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_umimapped-2.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_umimapped-2.bai", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_realigned-2.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_realigned-2.bam.bai", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_consensus_filtered.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_consensus_filtered.bai", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_clipoverlap.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_nodups.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/fastqs/{sample}_concatenated_1.fastq.gz", sample=all_clinseq_barcodes),
        expand(outdir + "/fastqs/{sample}_concatenated_2.fastq.gz", sample=all_clinseq_barcodes)
    output:
        outdir + "/bams/intermediate_bamfiles.removed"
    params:
        split_targets = outdir + "/bams/split_targets/"
    log:
        outdir + "/logs/remove_intermediate_{sample}.log".format(sample="_".join(all_clinseq_barcodes))
    run:
        del_bam = [bam for bam in input if 'clipoverlap' not in bam and 'nodups' not in bam]
        bamfiles = " ".join(del_bam)
        shell("rm {bamfiles} 2> {log} ")
        shell("rm -rf {params.split_targets} 2>> {log} ")
        shell("touch {output} ")