
nodups_suffix = "-nodups.bam" if CANCER_CAPTURE.project == "AL" else "_nodups.bam"


rule picard_collectinsertsize:
    input:
        bam = outdir + "/bams/{sample}" + nodups_suffix
    output:
        metrics = outdir + "/qc/picard/{sample}.picard-insertsize.txt"
    params:
        java_options = params['picard']['collectinsertsize']['java_options'],
        tmpdir = os.path.join(params['scratch'], 
                                "picard-isize-{}".format(str(uuid.uuid4())))
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
        bam = outdir + "/bams/{sample}" + nodups_suffix,
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


rule picard_collecthsmetrics:
    input:
        bam = outdir + "/bams/{sample}" + nodups_suffix,
        reference_genome = reference['reference_genome'],
        target_region = lambda wildcards: get_targets(wildcards, reference,
                                                      'targets-interval_list-slopped20'),
        bait_regions = lambda wildcards: get_targets(wildcards, reference,
                                                      'targets-interval_list-slopped20')
    output:
        metrics = outdir + "/qc/picard/{sample}.picard-hsmetrics.txt"
    params:
        bait_name = lambda wildcards: get_target_name(wildcards),
        java_options = params['picard']['collecthsmetrics']['java_options'],
        tmpdir = os.path.join(params['scratch'], 
                                "picard-hsmetrics-{}".format(str(uuid.uuid4())))
    threads: params['picard']['collecthsmetrics']['threads']
    log:
        outdir + "/logs/picard/picard_hsmetrics_{sample}.log"
    shell:
        "picard {params.java_options}  -Djava.io.tmpdir={params.tmpdir} " 
            "CollectHsMetrics  "
            "I={input.bam}   "
            "R={input.reference_genome} "
            "O={output.metrics} "
            "TI={input.target_region}  "
            "BI={input.bait_regions} "
            "BAIT_SET_NAME={params.bait_name} "
            "METRIC_ACCUMULATION_LEVEL=LIBRARY 2> {log} "
            " && rm -rf {params.tmpdir} " 


rule samtools_flagstat:
    input:
        bam = outdir + "/bams/{sample}" + nodups_suffix
    output:
        outdir + "/qc/samtools/{sample}-flagstats.json"
    threads: params['samtools']['flagstat']['threads']
    log:
        outdir + "/logs/samtools/samtools_flagstat_{sample}.log"
    shell:
        "samtools flagstat -@ {threads} -O json {input.bam} > {output} 2> {log} "


rule create_popvcf:
    input:
        popvcf = reference["swegene_common"],
        normal_target = reference['targets'][get_capture_name(NORMAL_CAPTURE.capture_kit_id)]['targets-bed-slopped20'],
        cancer_target = reference['targets'][get_capture_name(CANCER_CAPTURE.capture_kit_id)]['targets-bed-slopped20']
    output:
        outdir + "/contamination/pop_vcf_{}-{}.vcf".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        tmpdir = params['scratch']
    threads: params['create_popvcf']['threads']
    log:
        outdir + "/logs/contamination/pop_vcf_{}-{}.log".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "create_contest_vcfs.py {input.normal_target} "
            " {input.cancer_target} "
            " {input.popvcf}  "
            " --tmpdir {params.tmpdir}  "
            " --output-filename {output} 2> {log} "


rule gatk3_contest_cancer:
    input:
        reference_genome = reference['reference_genome'],
        normal_bam = capture_to_results[NORMAL_CAPTURE].bamfile,
        cancer_bam = capture_to_results[CANCER_CAPTURE].bamfile,
        popvcf = outdir + "/contamination/pop_vcf_{}-{}.vcf".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    output:
        "{}/contamination/{}.contest.txt".format(outdir, CANCER_CAPTURE_STR)
    params:
        tmpdir = os.path.join(params['scratch'], 
                                "gatk3-contest-cancer-{}".format(str(uuid.uuid4()))),
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
        popvcf = outdir + "/contamination/pop_vcf_{}-{}.vcf".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    output:
        "{}/contamination/{}.contest.txt".format(outdir, NORMAL_CAPTURE_STR)
    params:
        tmpdir = os.path.join(params['scratch'], 
                                "gatk3-contest-normal-{}".format(str(uuid.uuid4()))),
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


rule purecn:
    input:
        cnr = capture_to_results[CANCER_CAPTURE].cnr,
        seg = capture_to_results[CANCER_CAPTURE].seg,
        vardict_vcf = "{}/variants/{}-{}.vardict-somatic-purecn.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        csv = "{}/purecn/{}.csv".format(outdir, CANCER_CAPTURE_STR),
        genes_csv = "{}/purecn/{}_genes.csv".format(outdir, CANCER_CAPTURE_STR),
        variants_csv = "{}/purecn/{}_variants.csv".format(outdir, CANCER_CAPTURE_STR),
        loh_csv = "{}/purecn/{}_loh.csv".format(outdir, CANCER_CAPTURE_STR)
    params:
        tumorid = CANCER_CAPTURE_STR,
        minaf = params['purecn']['minaf'],
        maxnonclonal = params['purecn']['maxnonclonal'],
        outdir = "{}/purecn".format(outdir)
    threads: params['purecn']['threads']
    container: containers['purecn']
    log:
        "{}/logs/{}-purecn.log".format(outdir, CANCER_CAPTURE_STR)
    shell:
        "source activate purecn-env && "
        "PureCN.R  --out {params.outdir} "
        " --sampleid {params.tumorid} "
        " --segfile {input.seg} "
        " --tumor {input.cnr} "
        " --vcf {input.vardict_vcf} "
        " --genome hg19  --funsegmentation none "
        " --minpurity 0.05  --hzdev 0.1  {params.maxnonclonal} "
        " {params.minaf}  --error 0.0005 "
        " --postoptimize &&  "
        " conda deactivate && "
        " touch {output.csv} "
        " {output.genes_csv} "
        " {output.variants_csv} "
        " {output.loh_csv} 2> {log} "


cancer_capture_name = get_capture_name(CANCER_CAPTURE.capture_kit_id)
msisensor_prefix = f"{params['scratch']}/msisensor-{uuid.uuid4()}" 

rule msisensor:
    input:
        normal_bam = capture_to_results[NORMAL_CAPTURE].bamfile,
        tumor_bam = capture_to_results[CANCER_CAPTURE].bamfile,
        msi_sites = reference['targets'][cancer_capture_name]['msisites']
    output:
        "{}/msisensor-{}-{}.tsv".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        prefix = msisensor_prefix,
        table = "{}".format(msisensor_prefix),
        dis = "{}_dis".format(msisensor_prefix),
        germline = "{}_germline".format(msisensor_prefix),
        somatic = "{}_somatic".format(msisensor_prefix)
    threads: params["msisensor"]["threads"]
    log:
        "{}/logs/msisensor-{}-{}.log".format(outdir,NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "msisensor msi " 
            "-d {input.msi_sites} "
            "-n {input.normal_bam} "
            "-t {input.tumor_bam} "
            "-o {params.prefix} "
            "-b {threads} 2> {log}&& "
            " cp {params.prefix} {output} && "
            " rm {params.table} {params.dis} {params.germline} {params.somatic} "


bam_name = os.path.basename(capture_to_results[CANCER_CAPTURE].bamfile).split('.bam')[0]
msings_outdir = "{}/msings-{}".format(outdir, CANCER_CAPTURE_STR)
msings_output = "{}/{}/{}.MSI_Analysis.txt".format(msings_outdir, bam_name, bam_name)

rule msings:
    input:
        bam = capture_to_results[CANCER_CAPTURE].bamfile,
        reference_genome = reference["reference_genome"],
        msings_baseline = reference['targets'][cancer_capture_name]['msings-baseline'],
        msings_bed = reference['targets'][cancer_capture_name]['msings-bed'],
        msings_intervals = reference['targets'][cancer_capture_name]['msings-msi_intervals']
    output:
        msings = msings_output
    params:
        prefix = msings_outdir
    container: containers['gatk3']
    log:
        "{}/logs/msings-{}-{}.log".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "source activate gatk_3 && "
        "run_msings.sh -b {input.msings_bed} "
        " -f {input.reference_genome} "
        " -i {input.msings_intervals} "
        " -n {input.msings_baseline} "
        " -o {params.prefix}  {input.bam} 2> {log} "
        

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
        expand(outdir + "/qc/samtools/{sample}-flagstats.json", sample = clinseq_barcodes),
        "{}/contamination/{}.contest.txt".format(outdir, CANCER_CAPTURE_STR),
        "{}/contamination/{}.contest.txt".format(outdir, NORMAL_CAPTURE_STR),
        "{}/qc/{}-contam-qc-call.json".format(outdir, CANCER_CAPTURE_STR),
        msings_output
    output:
        "{}/qc/{}.qc_overview.pdf".format(outdir, "_".join(samples_capture_str))
    params:
        samples = ":".join(samples_capture_str),
        mainpath = os.path.dirname(os.path.dirname(outdir)),
        outdir = outdir
    container: containers['purecn']
    log:
        "{}/logs/qc_overview-{}-{}.log".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "source activate purecn-env && "
        "QC_overview.R  -s {params.samples} "
                        "-d {params.outdir} "
                        "-o {output} "
                        "-m {params.mainpath} 2> {log} "