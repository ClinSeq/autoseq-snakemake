

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
            shell("fastqc -o {output} -t {threads} --nogroup {fq}")


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
            " && rm -rf {params.tmpdir} "


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
            " && rm -rf {params.tmpdir} "


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



rule gatk3_contest_cancer:
    input:
        reference_genome = reference['reference_genome'],
        normal_bam = capture_to_results[NORMAL_CAPTURE].bamfile,
        cancer_bam = capture_to_results[CANCER_CAPTURE].bamfile,
        popvcf = reference["pop_vcf"]
    output:
        "{}/contamination/{}.contest.txt".format(outdir, CANCER_CAPTURE_STR)
    params:
        tmpdir = params['scratch'],
        min_genotype_ratio = params['contest_cancer']['min_genotype_ratio']
    threads: params['contest_cancer']['threads']
    container: containers['gatk3']
    log:
        outdir + "/logs/contamination/contest-{}.log".format(CANCER_CAPTURE_STR)
    shell:
        "source activate gatk_3 && "
        "gatk3 -Xmx15g -Djava.io.tmpdir={params.tmpdir} -T ContEst  "
            "-R {input.reference_genome}  "
            "-I:eval {input.cancer_bam}  "
            "-I:genotype {input.normal_bam} "
            "--popfile {input.popvcf}  "
            "--min_genotype_ratio {params.min_genotype_ratio}  "
            " -o {output} 2> {log} "
            " && rm -rf {params.tmpdir} "


rule gatk3_contest_normal:
    input:
        reference_genome = reference['reference_genome'],
        normal_bam = capture_to_results[NORMAL_CAPTURE].bamfile,
        cancer_bam = capture_to_results[CANCER_CAPTURE].bamfile,
        popvcf = reference["pop_vcf"]
    output:
        "{}/contamination/{}.contest.txt".format(outdir, NORMAL_CAPTURE_STR)
    params:
        tmpdir = params['scratch'],
        min_genotype_ratio = params['contest_cancer']['min_genotype_ratio']
    threads: params['contest_cancer']['threads']
    container: containers['gatk3']
    log:
        outdir + "/logs/contamination/contest-{}.log".format(NORMAL_CAPTURE_STR)
    shell:
        "source activate gatk_3 && "
        "gatk3 -Xmx15g -Djava.io.tmpdir={params.tmpdir} -T ContEst  "
            "-R {input.reference_genome}  "
            "-I:eval  {input.normal_bam} "
            "-I:genotype {input.cancer_bam} "
            "--popfile {input.popvcf}  "
            "--min_genotype_ratio {params.min_genotype_ratio}  "
            " -o {output} 2> {log} "
            " && rm -rf {params.tmpdir} "


rule contam_caveat:
    input:
        "{}/contamination/{}.contest.txt".format(outdir, CANCER_CAPTURE_STR)
    output:
        "{}/qc/{}-contam-qc-call.json".format(outdir, CANCER_CAPTURE_STR)
    threads: params['contam_caveat']['threads']
    log:
        "{}/logs/contamination/{}-contam-caveat.log".format(outdir, CANCER_CAPTURE_STR)
    shell:
        "contest_to_contam_caveat.py  "
        " {input}  > {output} 2> {log} "


rule multiqc:
    input:
        PICARD_QC,
        "{}/qc/{}-contam-qc-call.json".format(outdir, CANCER_CAPTURE_STR)
    output:
        directory("{}/multiqc".format(outdir))
    params:
        basefn = "{}-multiqc".format(CANCER_CAPTURE_STR),
        outdir = outdir,
    threads: params['multiqc']['threads']
    log:
        "{}/logs/multiqc-{}-{}.log".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "multiqc  {params.outdir} "
               "-o {output} "
               "-n {params.basefn} "
               "-k json " 
               " --data-dir --zip-data-dir -v -f 2> {log} "


rule overview_plot:
    input:
        PICARD_QC,
        expand(outdir + "/qc/samtools/{sample}-flagstats.json", sample = all_clinseq_barcodes),
        "{}/contamination/{}.contest.txt".format(outdir, CANCER_CAPTURE_STR),
        "{}/contamination/{}.contest.txt".format(outdir, NORMAL_CAPTURE_STR),
        "{}/qc/{}-contam-qc-call.json".format(outdir, CANCER_CAPTURE_STR)
    output:
        "{}/qc/{}.qc_overview.pdf".format(outdir, "_".join(samples_of_interest))
    params:
        samples = ":".join(samples_of_interest),
        mainpath = dirname(dirname(outdir)),
        outdir = outdir
    container: containers['purecn']
    log:
        "{}/logs/qc_overview-{}-{}.log".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "source activate purecn-env && "
        "QC_overview.R  -s {params.samples} "
                        "-d {params.outdir} "
                        "-o {output} --wgs "
                        "-m {params.mainpath} 2> {log} "
