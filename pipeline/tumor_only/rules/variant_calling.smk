import os
import sys
import uuid

capture_name = get_capture_name(CANCER_CAPTURE.capture_kit_id)
somatic_vcf = dict()



rule gatk4_mutect2:
    input:
        reference = reference['reference_genome'],
        normal_bam = capture_to_results[NORMAL_CAPTURE].umibam,
        tumor_bam = capture_to_results[CANCER_CAPTURE].umibam,
        interval_list = reference['targets'][capture_name]['targets-interval_list-slopped20']
    output:
        bam = "{}/variants/mutect/{}-{}-mutect.bam".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        filtered_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
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
        " -O {output.filtered_vcf} 2>> {log} "


somatic_vcf['mutect2'] = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)


rule sage_somatic:
    input:
        normal_bam = capture_to_results[NORMAL_CAPTURE].umibam,
        tumor_bam = capture_to_results[CANCER_CAPTURE].umibam,
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
        " bedtools merge -i {input.panel_bed} > {params.workdir}/targets_nonoverlap.bed 2> {log} && "
        " java -Xms4G -Xmx32G -cp $SAGE_JAR "
        " com.hartwig.hmftools.sage.SageApplication -threads 16 "
        " -tumor {params.tumorid} -tumor_bam {input.tumor_bam} "
        " -reference {params.normalid} -reference_bam {input.normal_bam}"
        " -ref_genome_version 37  -ref_genome {input.reference} "
        " -hotspots {input.known_hotspots} "
        " -panel_bed {params.workdir}/targets_nonoverlap.bed "
        " -ensembl_data_dir {input.ensembl_dir} "
        " -hard_min_tumor_vaf {params.minaf} -hotspot_min_tumor_qual {params.hpmintq} "
        " -hotspot_min_tumor_vaf {params.minaf} -min_map_quality {params.min_mapq} "
        " -min_avg_base_qual {params.min_baseq} -panel_min_tumor_qual {params.min_paneltq} "
        " -panel_min_tumor_vaf {params.minaf} -panel_only -write_bqr_data "
        " -out {output} 2>> {log} && rm {params.workdir}/targets_nonoverlap.bed"


rule bcftools_filter:
    input:
        sage_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        mutect_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    output:
        sage_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        mutect_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.pass.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    threads: 4
    params:
        tmp_sage_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        tmp_mutect_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.pass.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    log:
        "{}/logs/variants/{}-{}-sage-somatic-filter.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "bcftools filter -e 'FILTER!=\"PASS\"' {input.sage_vcf} > {params.tmp_sage_vcf} && "
        "bcftools filter -e 'FILTER!=\"PASS\"' {input.mutect_vcf} > {params.tmp_mutect_vcf} && "
        "bgzip {params.tmp_sage_vcf} && tabix -p vcf {output.sage_vcf} && "
        "bgzip {params.tmp_mutect_vcf} && tabix -p vcf {output.mutect_vcf} "


rule bcftools_concat:
    input:
        sage_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        mutect_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.pass.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    output:
        "{}/variants/{}-{}-all.somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: 8
    params:
        ordered_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.reordered.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    log:
        "{}/logs/variants/{}-{}-bcftools-concat.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "bcftools query -l {input.sage_vcf} | sort > sample_names.txt && "
        "bcftools view -Oz -S sample_names.txt {input.sage_vcf} -o {params.ordered_vcf} && "
        "tabix -p vcf {params.ordered_vcf} &&  "
        "bcftools concat -a -D {input.mutect_vcf} {params.ordered_vcf} "
        " | bgzip > {output}  2> {log} && " 
        "tabix -p vcf {output} && rm sample_names.txt {params.ordered_vcf}* "


# rule sage_splitvcf:
#     input:
#         "{}/variants/{}-{}-hartwig-sage-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     output:
#         snv = "{}/variants/{}-{}-hartwig-sage-somatic.snvs.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
#         indel =  "{}/variants/{}-{}-hartwig-sage-somatic.indels.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     threads: 1
#     params:
#         tmp_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR), 
#         snv = "{}/variants/{}-{}-hartwig-sage-somatic.snvs.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
#         indel =  "{}/variants/{}-{}-hartwig-sage-somatic.indels.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     container: containers['somaticseq']
#     log:
#         "{}/logs/variants/{}-{}-sage-somatic-splitvcf.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     shell:
#         "source activate somaticseqenv && "
#         "bcftools filter -e 'FILTER!=\"PASS\"' {input} > {params.tmp_vcf} && "
#         "splitVcf.py -infile {params.tmp_vcf} -snv {params.snv} -indel {params.indel} && "
#         " bgzip {params.snv} && bgzip {params.indel} "


# somatic_vcf['sage_snv'] = "{}/variants/{}-{}-hartwig-sage-somatic.snvs.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
# somatic_vcf['sage_indel'] = "{}/variants/{}-{}-hartwig-sage-somatic.indels.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)


# rule somaticseq_merge:
#     input:
#         **somatic_vcf,
#         reference = reference['reference_genome'],
#         normal_bam = capture_to_results[NORMAL_CAPTURE].umibam,
#         tumor_bam = capture_to_results[CANCER_CAPTURE].umibam
#     output:
#         rundir = directory("{}/variants/{}-{}-somaticseq".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)),
#         consensus_snv = "{}/variants/{}-{}-somaticseq/Consensus.sSNV.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
#         consensus_indel = "{}/variants/{}-{}-somaticseq/Consensus.sINDEL.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     params:
#         tmpdir = params['scratch']
#     threads: params['somaticseq']['threads']
#     container: containers['somaticseq']
#     log:
#         "{}/logs/variants/{}-{}-somaticseq.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     shell:
#         "source activate somaticseqenv && "
#         "run_somaticseq.py  --output-directory {output.rundir} "
#         " --genome-reference {input.reference} paired "
#         " --tumor-bam-file {input.tumor_bam} " 
#         " --normal-bam-file {input.normal_bam} " 
#         " --mutect2-vcf {input.mutect2} "
#         " --arbitrary-snvs {input.sage_snv} "
#         " --arbitrary-indels {input.sage_indel} 2> {log}"


# rule gatk3_combinevariants:
#     input:
#         reference = reference['reference_genome'],
#         consensus_snv = "{}/variants/{}-{}-somaticseq/Consensus.sSNV.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
#         consensus_indel = "{}/variants/{}-{}-somaticseq/Consensus.sINDEL.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     output:
#         "{}/variants/{}-{}-all.somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     threads: params['somaticseq']['threads']
#     container: containers['gatk3']
#     log:
#         "{}/logs/variants/{}-{}-combine_somaticvcf.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     shell:    
#         " source activate gatk_3 && "
#         " gatk3 -T CombineVariants "
#         " -R {input.reference} --variant {input.consensus_snv} " 
#         " --variant {input.consensus_indel} " 
#         " --assumeIdenticalSamples  | bgzip > {output} && "
#         " tabix -p vcf {output} 2> {log} "


# germline variant callers on tumor sample
rule gatk4_haplotypecaller:
    input:
        bam = capture_to_results[CANCER_CAPTURE].umibam,
        reference = reference['reference_genome'],
        dbsnp = reference["dbSNP"],
        interval_list = reference['targets'][capture_name]['targets-interval_list-slopped20'],
    output:
        vcf = "{}/variants/haplotypecaller/{}.haplotypecaller-germline.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    params:
        java_options = params["gatk4"]["haplotypecaller"]["java_options"]
    threads: params["gatk4"]["threads"]
    log:
        "{}/logs/variants/haplotypecaller/{}.haplotypecaller-germline.log".format(outdir, CANCER_CAPTURE_STR)
    shell:
        "gatk --java-options '{params.java_options}' "
            " HaplotypeCaller   "
            " -R {input.reference}  "
            " -I {input.bam}  "
            " -L {input.interval_list} "
            " --dbsnp {input.dbsnp} "
            " -O {output.vcf} 2> {log} "


rule vt_decomp_norm_hc:
    input:
        reference = reference['reference_genome'],
        vcf = "{}/variants/haplotypecaller/{}.haplotypecaller-germline.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    output:
        "{}/variants/{}-merged.germline.split_norm.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    threads: 1
    log:
        "{}/logs/variants/haplotypecaller/{}.vt_decomp_norm_haplotypecaller.log".format(outdir, CANCER_CAPTURE_STR)
    shell:
        "vt decompose -s {input.vcf} | vt normalize -r {input.reference} - "
        " | bgzip > {output} 2> {log} && "
        " tabix -p vcf {output} "


### SNP BAFs ###
rule vcf_add_sample:
    input:
        vcf = "{}/variants/{}-merged.germline.split_norm.gnomADg.vep.SNPs.vcf.gz".format(outdir, CANCER_CAPTURE_STR),
        tumor_bam = capture_to_results[CANCER_CAPTURE].umibam
    output:
        "{}/variants/{}-merged.germline.split_norm.gnomADg.vep.SNPs.BAF.vcf.gz".format(outdir, CANCER_CAPTURE_STR),
    params:
        tumorid = CANCER_SAMPLE_STR,
        tmpdir = os.path.join(params['scratch'], 
                    "vcfaddsample-{}".format(str(uuid.uuid4())))
    threads: params['vcfaddsample']['threads']
    log:
        "{}/logs/variants/{}.vcf_add_sample.log".format(outdir, CANCER_CAPTURE_STR)
    shell:
        "vcf_filter.py --no-filtered {input.vcf} sq --site-quality 5 "
        " | bcftools view --exclude 'FORMAT/AD=\".\"' -Oz > {params.tmpdir}.vcf.gz  && "
        " vcf_add_sample.py --filter_hom --samplename {params.tumorid}-BAF "
        " {params.tmpdir}.vcf.gz {input.tumor_bam} | bgzip > {output} 2> {log} && "
        " tabix -p vcf {output} && rm -v {params.tmpdir}.vcf.gz "


rule bcftools_dbsnp_annotation:
    input:
        vcf = "{}/variants/{}-merged.germline.split_norm.gnomADg.vep.SNPs.vcf.gz".format(outdir, CANCER_CAPTURE_STR),
        dbsnp = reference["dbSNP"]
    output:
        "{}/variants/{}-merged.germline.split_norm.gnomADg.vep.SNPs.dbsnpids.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    log:
        "{}/logs/variants/{}.bcftools_dbsnp_annotation.log".format(outdir, CANCER_CAPTURE_STR)
    shell:
        "bcftools annotate --annotation {input.dbsnp} --columns ID --output-type u {input.vcf} "
        " | bcftools view --exclude 'FORMAT/AD=\".\"' --output-file {output} -Oz 2> {log} && "
        "tabix -p vcf {output} "


rule somatic_generateIGVnav:
    input:
        somatic = "{}/variants/{}-{}-all.somatic.gnomADg.noSNPs.brcaEx.vep.vcf.gz".format(outdir, 
                                                            CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        oncokb = reference['oncokb'],
        cgcann = reference['cgcann']
    output:
        "{}/{}-somatic-igvnav-input.txt".format(outdir, CANCER_CAPTURE_STR)
    params:
        tmp_vcf = "{}/variants/{}-{}-all.somatic.gnomADg.noSNPs.brcaEx.tmp.vep.vcf".format(outdir, 
                                                            CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        vcftype = "somatic"
    log:
        "{}/logs/variants/{}.somatic_generateIGVnav.log".format(outdir, CANCER_CAPTURE_STR)
    shell:
        "zcat {input.somatic} > {params.tmp_vcf} && "
        "generateIGVnavInput.py {params.tmp_vcf} {input.oncokb} {params.vcftype} "
        " --cgc {input.cgcann} --output {output} 2> {log} "


rule germline_generateIGVnav:
    input:
        vcf = "{}/variants/{}-merged.germline.split_norm.brcaEx.vep.vcf.gz".format(outdir, 
                                                                        CANCER_CAPTURE_STR),
        oncokb = reference['oncokb'],
        cgcann = reference['cgcann']
    output:
        "{}/{}-germline-igvnav-input.txt".format(outdir, CANCER_CAPTURE_STR)
    params:
        vcftype = "germline",
        tmp_vcf = "{}/variants/{}-merged.germline.tmp.vcf".format(outdir, 
                                                                     CANCER_CAPTURE_STR)
    log:
        "{}/logs/variants/{}.germline_generateIGVnav.log".format(outdir, CANCER_CAPTURE_STR)
    shell:
        "zcat {input.vcf} | awk -F '\\t' -v OFS='\\t' '{{if ($1 ~ /^#/) print $0; else if ($5 != \"*\") {{print $0}}}}' " 
        " | bcftools view --exclude 'FORMAT/AD=\".\"' -Ov > {params.tmp_vcf} && "
        "generateIGVnavInput.py {params.tmp_vcf} {input.oncokb} {params.vcftype} "
        " --cgc {input.cgcann} --output {output} 2> {log} && "
        "rm -v {params.tmp_vcf} "
        

