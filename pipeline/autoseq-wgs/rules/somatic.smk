

somatic_vcf = dict()
normalBam = capture_to_results[NORMAL_CAPTURE].bamfile
cancerBam = capture_to_results[CANCER_CAPTURE].bamfile

targets_dir = reference["wgs"]["targets"]["bed"]
suffix = list()
for fn in os.listdir(targets_dir):
    suffix.append(fn.split('.')[1])


intervals_dir = reference["wgs"]["targets"]["interval_list"]
mutect_vcf_prefix = "{}/variants/mutect/{}-{}-gatk-mutect-somatic".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
mutext_log_prefix = "{}/logs/variants/{}-{}-gatk4-mutect-somatic".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 

rule gatk4_mutect2:
    input:
        reference = reference['reference_genome'],
        normal_bam = normalBam,
        tumor_bam = cancerBam,
        interval_list = intervals_dir + "human_g1k_v37_decoy.{suf}.interval_list"
    output:
        bam = mutect_vcf_prefix + "{suf}.bam",
        vcf = mutect_vcf_prefix + "{suf}.vcf.gz",
        filtered_vcf = mutect_vcf_prefix + "{suf}-filtered.vcf.gz"
    wildcard_constraints:
        suf = "|".join(suffix)
    params:
        java_options = params['gatk4']['mutect2']['java_options'],
        normalid = compose_sample_str(NORMAL_CAPTURE),
        tumorid = compose_sample_str(CANCER_CAPTURE),
        tmpdir = params['scratch']
    threads: params['gatk4']['threads']
    log:
        mutext_log_prefix + "{suf}.log"        
    shell:
        "gatk --java-options '{params.java_options} -Djava.io.tmpdir={params.tmpdir}' "
        " Mutect2  -R {input.reference} "
        " -I {input.tumor_bam} "
        " -I {input.normal_bam} "
        " -tumor {params.tumorid}  " 
        " -normal {params.normalid} " 
        " -L {input.interval_list} "
        " -bamout {output.bam} -O {output.vcf} 2> {log} && "
        " gatk --java-options '-Xmx10g -Djava.io.tmpdir={params.tmpdir}' "
        " FilterMutectCalls  -R {input.reference} "
        " -V {output.vcf}  "
        " -O {output.filtered_vcf} 2>> {log} "


rule mutect2_vcfmerge:
    input:
        expand(mutect_vcf_prefix + "{suf}-filtered.vcf.gz", suf=suffix)
    output:
        mutect_vcf_prefix + "-filtered.vcf.gz"     
    params:
        vcf_prefix = mutect_vcf_prefix
    threads: params['bcftools']['threads']
    log:
        "{}/logs/variants/{}-{}-mutect-somatic-merge.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        "bcftools concat --threads {threads} -a -D "
        " -O z {input} > {output} 2> {log} && "
        " tabix -p vcf {output} && "
        " rm {params.vcf_prefix}000*.vcf.gz* "
        " {params.vcf_prefix}000*.ba* "


rule mutect2_normalize:
    input:
        vcf = mutect_vcf_prefix + "-filtered.vcf.gz",
        reference = reference['reference_genome']
    output:
        normalized_vcf = mutect_vcf_prefix + "-filtered-normalized.vcf.gz"
    params:
        tmpdir = params['scratch']
    threads: params['bcftools']['threads']
    log:
        "{}/logs/variants/{}-{}-mutect-somatic-merge.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        " vt decompose -s {input.vcf} "
        " | vt normalize  -r {input.reference} - "
        " | bgzip > {output.normalized_vcf} 2> {log} "


somatic_vcf['mutect2'] = mutect_vcf_prefix + "-filtered-normalized.vcf.gz"

rule sage_somatic:
    input:
        normal_bam = normalBam,
        tumor_bam = cancerBam,
        reference = reference['reference_genome'],
        panel_bed = reference['wgs']['hartwig']['actionable-somatic-panel-bed'],
        known_hotspots = reference['wgs']['hartwig']['known-hotspots-somatic-vcf'],
        high_confi_bed = reference['wgs']['hartwig']['NA12878-highconf-bed'],
        ensembl_dir = reference['wgs']['hartwig']['ensembl-dir']
    output:
        "{}/variants/{}-{}-hartwig-sage-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        normalid = compose_sample_str(NORMAL_CAPTURE),
        tumorid = compose_sample_str(CANCER_CAPTURE),
        jarfile = os.environ.get('SAGE_JAR')
    threads: params['sage']['threads']
    container: containers['gridss']
    log:
        "{}/logs/variants/{}-{}-sage-somatic.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        " source activate gridss-env && "
        " java -Xms4G -Xmx32G -cp $SAGE_JAR "
        " com.hartwig.hmftools.sage.SageApplication -threads 16 "
        " -reference {params.normalid} -reference_bam {input.normal_bam}"
        " -tumor {params.tumorid} -tumor_bam {input.tumor_bam} "
        " -ref_genome_version 37  -ref_genome {input.reference} "
        " -hotspots {input.known_hotspots} "
        " -panel_bed {input.panel_bed} "
        " -high_confidence_bed {input.high_confi_bed} "
        " -ensembl_data_dir {input.ensembl_dir} "
        " -out {output} 2> {log} "


rule bcftools_filter:
    input:
        sage_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        mutect_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered-normalized.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        sage_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        mutect_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.pass.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: 4
    params:
        tmp_sage_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        tmp_mutect_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.pass.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
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
        mutect_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.pass.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        "{}/variants/{}-{}-all.somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: 8
    params:
        ordered_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.reordered.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
    log:
        "{}/logs/variants/{}-{}-bcftools-concat.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "bcftools query -l {input.sage_vcf} | sort > sample_names.txt && "
        "bcftools view -Oz -S sample_names.txt {input.sage_vcf} -o {params.ordered_vcf} && "
        "tabix -p vcf {params.ordered_vcf} &&  "
        "bcftools concat -a -D {input.mutect_vcf} {params.ordered_vcf} "
        " | bgzip > {output}  2> {log} && " 
        "tabix -p vcf {output} && rm sample_names.txt {params.ordered_vcf}* "


rule somatic_generateIGVnav:
    input:
        somatic = "{}/variants/{}-{}-all.somatic.vep.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        oncokb = reference['oncokb'],
        cgcann = reference['cgcann'],
        curation_ann = reference['curation_ann'],
        cancer_hotspot_snv = reference['cancer_hotspot_snv'],
        cancer_hotspot_indel = reference['cancer_hotspot_indel']
    output:
        "{}/{}-{}-igvnav-input.txt".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        vcftype = "somatic"
    log:
        "{}/logs/{}-{}.somatic_generate_igvnav.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        "generateIGVnavInput.py {input.somatic} {input.oncokb} {params.vcftype} "
        " -c {input.curation_ann} "
        " --hotspot-snv {input.cancer_hotspot_snv} "
        " --hotspot-indel {input.cancer_hotspot_indel} "
        " --wgs --cgc {input.cgcann} --output {output} 2> {log} "


haplotype_jc_vcf_prefix = "{}/variants/haplotypecaller/{}-{}.haplotypecaller-joint-calling".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
haplotype_jc_log_prefix = "{}/logs/variants/haplotypecaller/{}-{}.haplotypecaller-joint-calling".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)


rule gatk4_haplotypecaller_jc:
    input:
        normal_bam = normalBam,
        tumor_bam = cancerBam,
        reference = reference['reference_genome'],
        dbsnp = reference["dbSNP"],
        interval_list = intervals_dir + "human_g1k_v37_decoy.{suf}.interval_list"
    output:
        vcf = haplotype_jc_vcf_prefix + ".{suf}.vcf.gz",
    wildcard_constraints:
        suf = "|".join(suffix)
    params:
        java_options = params["gatk4"]["haplotypecaller"]["java_options"]
    threads: 4
    log:
        haplotype_jc_log_prefix + ".{suf}.log"
    shell:
        "gatk --java-options '{params.java_options}' "
        " HaplotypeCaller   "
        " -R {input.reference}  "
        " -I {input.normal_bam}  "
        " -I {input.tumor_bam} "
        " -L {input.interval_list} "
        " --dbsnp {input.dbsnp} "
        " -O {output.vcf} 2> {log} "


rule picard_vcfmerge:
    input:
        expand(haplotype_jc_vcf_prefix + ".{suf}.vcf.gz", suf = suffix)
    output:
        vcf = haplotype_jc_vcf_prefix + ".vcf.gz"
    params:
        java_options = params["gatk4"]["haplotypecaller"]["java_options"]
    threads: params["gatk4"]["threads"]
    log:
        haplotype_jc_log_prefix + ".log"
    shell:
        """
        InputBams=({input})
        ibams=`for i in ${{InputBams[@]}}; do echo "I="$i; done | tr '\\n' '\\t' `
        picard MergeVcfs $ibams O={output.vcf} 2> {log}
        """


rule make_allelic_fraction_track:
    input:
        vcf = haplotype_jc_vcf_prefix + ".vcf.gz"
    output:
        "{}/variants/{}-and-{}.germline-variants-somatic-afs.bedGraph".format(
            outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    log:
       "{}/logs/{}-{}.make_allelic_fraction_track.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        "generate_allelic_fraction_bedGraph.py  --output {output}  {input.vcf} 2> {log} "


rule gatk4_haplotypecaller_tumor:
    input:
        reference = reference["reference_genome"],
        tumor_bam = cancerBam,
        germline_vcf = "{}/variants/{}-all.germline.vep.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
    output:
        vcf = "{}/variants/{}-haplotypecaller.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    params:
        java_options = params["gatk4"]["haplotypecaller"]["java_options"],
        tmpdir = params['scratch']
    threads: params['gatk4']['threads']
    log:
        "{}/logs/variants/{}-{}-gatk4-haplotypecaller-tumor.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        "gatk --java-options '{params.java_options}' "
        " HaplotypeCaller -R {input.reference} "
        " -I {input.tumor_bam} "
        " -O {output.vcf} "
        " -ERC GVCF "
        " -L {input.germline_vcf} 2> {log} "


rule gatk4_haplotypecaller_genotype_tumor:
    input:
        reference = reference["reference_genome"],
        vcf = "{}/variants/{}-haplotypecaller.vcf.gz".format(outdir, CANCER_CAPTURE_STR),
        germline_vcf = "{}/variants/{}-all.germline.vep.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
    output:
        vcf = "{}/variants/{}-haplotypecaller.genotyped.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    params:
        java_options = params["gatk4"]["haplotypecaller"]["java_options"],
        vcf = "{}/variants/{}-haplotypecaller.genotyped.vcf".format(outdir, CANCER_CAPTURE_STR)
    threads: params['gatk4']['threads']
    log: 
        "{}/logs/variants/{}-{}-haplotypecaller-tumor-genotyped.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        "gatk --java-options '{params.java_options}' "
        " GenotypeGVCFs -R {input.reference} "
        "  -V {input.vcf}  "
        "  -O {params.vcf} "
        "  -L {input.germline_vcf} 2> {log} && "
        " bgzip {params.vcf} && "
        " tabix -p vcf {output.vcf} 2>> {log} "
    

rule bcftools_merge:
    input:
        germline_vcf = "{}/variants/{}-all.germline.vep.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
        genotyped_tvcf = "{}/variants/{}-haplotypecaller.genotyped.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    output:
        vcf = "{}/variants/{}-{}.germline_variants_with_taf.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    threads: params['gatk4']['threads']
    log:
       "{}/logs/variants/{}-{}-bcftools-merge.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)  
    shell:
        "bcftools merge {input.germline_vcf} "
        " {input.genotyped_tvcf} "
        " -O z -o {output.vcf} 2> {log} "
