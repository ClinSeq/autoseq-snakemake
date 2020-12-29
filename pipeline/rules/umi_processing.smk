import os
import uuid 


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
    log: outdir + "logs/fgbio/fastqtobam_{sample}.log"
    shell:
        "fgbio {params.java_options} --tmp-dir {params.tmpdir} FastqToBam "
            " -i {input.fq1}  {input.fq2} " 
            " -o {output.bam} "
            " --sample {params.sample} "
            " --library {params.library} "
            " -r 3M2S+T 3M2S+T "
            " -s true "
            " && rm -rf {params.tmpdir}"


rule bwa_umialignment:
    input:
        bam = outdir + "/bams/{sample}_unmapped.bam",
        reference_genome = reference['bwaIndex']
    output:
        bam = outdir + "/bams/{sample}_umimapped.bam"
    params:
        java_options = params['picard']['merge_bam']['java_options'],
        tmpdir = os.path.join(params['scratch'], 
                    "bwa-umialignment-{}".format(str(uuid.uuid4())))
    threads: params['bwa']['threads']
    log: outdir + "logs/fgbio/bwa_umialignment_{sample}.log"
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


rule gatk3_targetcreator_umi_1:
    input:
        bam = outdir + "/bams/{sample}_umimapped.bam",
        reference_genome = reference['reference_genome'],
        target_region = lambda wildcards: get_targets(wildcards, reference),
        known_1kg = reference["1KG"],
        known_mills_gs = reference["Mills_and_1KG_gold_standard"]
    output:
        target_intervals = outdir + "/bams/{sample}_umi.intervals"
    params:
        java_options = params['gatk3']['realigner_target']['java_options'],
        jarfile = params['gatk3']['jarfile'],
        extra = params['gatk3']['realigner_target']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "realignerTC-{}".format(str(uuid.uuid4())))
    threads: params['gatk3']['realigner_target']['threads']
    log:
        outdir + "/logs/gatk_realigner_targetcreator_{sample}.log"
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
            " -o {output.target_intervals} 2> {log} && "
        "source deactivate "


rule gatk3_indelrealigner_umi_1:
    input:
        bam = outdir + "/bams/{sample}_umimapped.bam",
        reference_genome = reference['reference_genome'],
        target_region = lambda wildcards: get_targets(wildcards, reference),
        known_1kg = reference["1KG"],
        known_mills_gs = reference["Mills_and_1KG_gold_standard"],
        target_intervals = outdir + "/bams/{sample}_umi.intervals"
    output:
        bam = outdir + "/bams/{sample}_realigned-1.bam"
    params:
        jarfile = params['gatk3']['jarfile'],
        java_options = params['gatk3']['indel_realigner']['java_options'],
        extra = params['gatk3']['indel_realigner']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "indelrealigner-{}".format(str(uuid.uuid4())))
    threads: params['gatk3']['indel_realigner']['threads']
    log:
        outdir + "/logs/gatk_indel_realigner_{sample}.log"
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
        "source deactivate "


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
              " --family-size-histogram {output.histogram} "
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
        bam = outdir + "/bams/{sample}_umimapped-2.bam"
    params:
        java_options = params['picard']['merge_bam']['java_options'],
        tmpdir = os.path.join(params['scratch'], 
                    "bwa-umialignment-{}".format(str(uuid.uuid4())))
    threads: params['bwa']['threads']
    log: outdir + "logs/fgbio/bwa_umialignment_{sample}.log"
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
        bam = outdir + "/bams/{sample}_umimapped-2.bam",
        reference_genome = reference['reference_genome'],
        target_region = lambda wildcards: get_targets(wildcards, reference),
        known_1kg = reference["1KG"],
        known_mills_gs = reference["Mills_and_1KG_gold_standard"]
    output:
        target_intervals = outdir + "/bams/{sample}_consensus.intervals"
    params:
        java_options = params['gatk3']['realigner_target']['java_options'],
        jarfile = params['gatk3']['jarfile'],
        extra = params['gatk3']['realigner_target']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "realignerTC-{}".format(str(uuid.uuid4())))
    threads: params['gatk3']['realigner_target']['threads']
    log:
        outdir + "/logs/gatk_realigner_targetcreator_{sample}.log"
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


rule gatk3_indelrealigner_umi_2:
    input:
        bam = outdir + "/bams/{sample}_umimapped-2.bam",
        reference_genome = reference['reference_genome'],
        target_region = lambda wildcards: get_targets(wildcards, reference),
        known_1kg = reference["1KG"],
        known_mills_gs = reference["Mills_and_1KG_gold_standard"],
        target_intervals = outdir + "/bams/{sample}_consensus.intervals"
    output:
        bam = outdir + "/bams/{sample}_realigned-2.bam"
    params:
        jarfile = params['gatk3']['jarfile'],
        java_options = params['gatk3']['indel_realigner']['java_options'],
        extra = params['gatk3']['indel_realigner']['extra'],
        tmpdir = os.path.join(params['scratch'], 
                                "indelrealigner-{}".format(str(uuid.uuid4())))
    threads: params['gatk3']['indel_realigner']['threads']
    log:
        outdir + "/logs/gatk_indel_realigner_{sample}.log"
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
            " -o {output.bam} 2> {log} "


rule fgbio_filterconsensus:
    input:
        bam = outdir + "/bams/{sample}_realigned-2.bam",
        reference_genome = reference['reference_genome']
    output:
        bam = outdir + "/bams/{sample}_consensus_filtered.bam"
    params:
        java_options = params['fgbio']['filterconsensus']['java_options'],
        error_rate = params['fgbio']['filterconsensus']['error_rate'],
        base_quality = params['fgbio']['filterconsensus']['base_quality'],
        extra = params['fgbio']['filterconsensus']['extra'],
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


rule fgbio_clipbam:
    input:
        bam = outdir + "/bams/{sample}_consensus_filtered.bam",
        reference_genome = reference['reference_genome']
    output:
        bam = outdir + "/bams/{sample}_clipoverlap.bam",
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
