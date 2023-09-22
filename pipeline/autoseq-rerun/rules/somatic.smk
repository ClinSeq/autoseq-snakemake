

rule gatk4_mutect2:
    input:
        reference = reference['reference_genome'],
        normal_bam = normalBam,
        tumor_bam = cancerBam,
        interval_list = reference['targets'][capture_name]['targets-interval_list-slopped20']
    output:
        bam = "{}/variants/mutect/{}-{}-mutect.bam".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        filtered_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        normalized_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered-normalized.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        java_options = params['gatk4']['mutect2']['java_options'],
        normalid = compose_sample_str(NORMAL_CAPTURE),
        tumorid = compose_sample_str(CANCER_CAPTURE),
        tmpdir = params['scratch']
    threads: params['gatk4']['threads']
    log:
        "{}/logs/variants/{}-{}-gatk4-mutect-somatic.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        "gatk --java-options '{params.java_options} -Djava.io.tmpdir={params.tmpdir}' "
        " Mutect2  -R {input.reference} "
        " -I {input.tumor_bam} "
        " -I {input.normal_bam} "
        " -tumor {params.tumorid}  " 
        " -normal {params.normalid} " 
        " -L {input.interval_list} " 
        " --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter "
        " -bamout {output.bam} -O {output.vcf} 2> {log} && "
        " gatk --java-options '-Xmx10g -Djava.io.tmpdir={params.tmpdir}' "
        " FilterMutectCalls  -R {input.reference} "
        " --max-alt-allele-count 2 "
        " -V {output.vcf}  "
        " -O {output.filtered_vcf} 2>> {log} && "
        " vt decompose -s {output.filtered_vcf} "
        " | vt normalize  -r {input.reference} - "
        " | bgzip > {output.normalized_vcf} 2>> {log} "


somatic_vcf['mutect2'] = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered-normalized.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)


rule sage_somatic:
    input:
        normal_bam = normalBam,
        tumor_bam = cancerBam,
        reference = reference['reference_genome'],
        panel_bed = reference['targets'][capture_name]['targets-bed-slopped20'],
        known_hotspots = reference['wgs']['hartwig']['known-hotspots-somatic-vcf'],
        high_confi_bed = reference['wgs']['hartwig']['NA12878-highconf-bed'],
        ensembl_dir = reference['wgs']['hartwig']['ensembl-dir']
    output:
        "{}/variants/{}-{}-hartwig-sage-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        normalid = compose_sample_str(NORMAL_CAPTURE),
        tumorid = compose_sample_str(CANCER_CAPTURE),
        jarfile = os.environ.get('SAGE_JAR'),
        workdir = outdir + "/variants",
        minaf = 0.0002,
        hpmintq = 150,
        min_mapq = 20,
        min_baseq = 30,
        min_paneltq = 250
    threads: params['sage']['threads']
    container: containers['gridss']
    log:
        "{}/logs/variants/{}-{}-sage-somatic.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        " source activate gridss-env && "
        " bedtools merge -s -i {input.panel_bed} > {params.workdir}/targets_nonoverlap.bed 2> {log} && "
        " java -Xms4G -Xmx32G -cp {params.jarfile} "
        " com.hartwig.hmftools.sage.SageApplication -threads 16 "
        " -reference {params.normalid} -reference_bam {input.normal_bam}"
        " -tumor {params.tumorid} -tumor_bam {input.tumor_bam} "
        " -ref_genome_version 37  -ref_genome {input.reference} "
        " -hotspots {input.known_hotspots} "
        " -panel_bed {params.workdir}/targets_nonoverlap.bed "
        " -ensembl_data_dir {input.ensembl_dir} "
        " -hard_min_tumor_vaf {params.minaf} -hotspot_min_tumor_qual {params.hpmintq} "
        " -hotspot_min_tumor_vaf {params.minaf} -min_map_quality {params.min_mapq} "
        " -min_avg_base_qual {params.min_baseq} -panel_min_tumor_qual {params.min_paneltq} "
        " -panel_min_tumor_vaf {params.minaf} -panel_only -write_bqr_data "
        " -out {output} 2>> {log} && rm {params.workdir}/targets_nonoverlap.bed"
    

rule sage_splitvcf:
    input:
        "{}/variants/{}-{}-hartwig-sage-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        snv = "{}/variants/{}-{}-hartwig-sage-somatic.snvs.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        indel =  "{}/variants/{}-{}-hartwig-sage-somatic.indels.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: 1
    params:
        tmp_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR), 
        snv = "{}/variants/{}-{}-hartwig-sage-somatic.snvs.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        indel =  "{}/variants/{}-{}-hartwig-sage-somatic.indels.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    container: containers['somaticseq']
    log:
        "{}/logs/variants/{}-{}-sage-somatic-splitvcf.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "source activate somaticseqenv && "
        "bcftools filter -e 'FILTER!=\"PASS\"' {input} > {params.tmp_vcf} && "
        "splitVcf.py -infile {params.tmp_vcf} -snv {params.snv} -indel {params.indel} && "
        " bgzip {params.snv} && bgzip {params.indel} "


somatic_vcf['sage_snv'] = "{}/variants/{}-{}-hartwig-sage-somatic.snvs.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
somatic_vcf['sage_indel'] = "{}/variants/{}-{}-hartwig-sage-somatic.indels.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)


rule somaticseq_merge:
    input:
        **somatic_vcf,
        reference = reference['reference_genome'],
        normal_bam = normalBam,
        tumor_bam = cancerBam
    output:
        rundir = directory("{}/variants/{}-{}-somaticseq".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)),
        consensus_snv = "{}/variants/{}-{}-somaticseq/Consensus.sSNV.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        consensus_indel = "{}/variants/{}-{}-somaticseq/Consensus.sINDEL.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        tmpdir = params['scratch']
    threads: params['somaticseq']['threads']
    container: containers['somaticseq']
    log:
        "{}/logs/variants/{}-{}-somaticseq.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "source activate somaticseqenv && "
        "run_somaticseq.py  --output-directory {output.rundir} "
        " --genome-reference {input.reference} paired "
        " --tumor-bam-file {input.tumor_bam} " 
        " --normal-bam-file {input.normal_bam} " 
        " --mutect2-vcf {input.mutect2} "
        " --arbitrary-snvs {input.sage_snv} "
        " --arbitrary-indels {input.sage_indel} 2> {log} "


rule gatk3_combinevariants:
    input:
        reference = reference['reference_genome'],
        consensus_snv = "{}/variants/{}-{}-somaticseq/Consensus.sSNV.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        consensus_indel = "{}/variants/{}-{}-somaticseq/Consensus.sINDEL.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        "{}/variants/{}-{}-all.somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: params['somaticseq']['threads']
    container: containers['gatk3']
    log:
        "{}/logs/variants/{}-{}-combine_somaticvcf.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:    
        " source activate gatk_3 && "
        " gatk3 -T CombineVariants "
        " -R {input.reference} --variant {input.consensus_snv} " 
        " --variant {input.consensus_indel} " 
        " --assumeIdenticalSamples  | bgzip > {output} && "
        " tabix -p vcf {output} 2> {log} "