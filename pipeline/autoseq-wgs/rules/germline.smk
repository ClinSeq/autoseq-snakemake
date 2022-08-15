bamfile = capture_to_results[NORMAL_CAPTURE].umibam if umi else capture_to_results[NORMAL_CAPTURE].bamfile

rule gatk4_haplotypecaller:
    input:
        bam = bamfile,
        reference = reference['reference_genome'],
        dbsnp = reference["dbSNP"]
    output:
        vcf = "{}/variants/haplotypecaller/{}.haplotypecaller-germline.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
        normalized_vcf = "{}/variants/haplotypecaller/{}.haplotypecaller-germline-normalized.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
    params:
        java_options = params["gatk4"]["haplotypecaller"]["java_options"]
    threads: params["gatk4"]["threads"]
    log:
        "{}/logs/variants/haplotypecaller/{}.haplotypecaller-germline.log".format(outdir, NORMAL_CAPTURE_STR)
    shell:
        "gatk --java-options '{params.java_options}' "
            " HaplotypeCaller   "
            " -R {input.reference}  "
            " -I {input.bam}  "
            " --dbsnp {input.dbsnp} "
            " -O {output.vcf} && "
        " vt decompose -s {output.vcf} "
        " | vt normalize  -r {input.reference} - "
        " | bgzip > {output.normalized_vcf}  && "
        " tabix -p vcf {output.normalized_vcf} 2> {log} "


rule strelka_germline:
    input:
        bam = bamfile,
        reference = reference['reference_genome']
    output:
        vcf = "{}/variants/{}-strelka-germline.passed.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
    params:
        rundir = directory("{}/variants/{}-strelka-germline".format(outdir, NORMAL_CAPTURE_STR))
    threads: params["strelka"]["threads"]
    log:
        "{}/logs/variants/{}.strelka-germline.log".format(outdir, NORMAL_CAPTURE_STR)
    shell:
        "source activate gatk_3 && "
        "configureStrelkaGermlineWorkflow.py  --bam {input.bam} "
        " --ref {input.reference} --targeted "
        " --callRegions {input.call_region} "
        " --runDir {params.rundir} && "
        " {params.rundir}/runWorkflow.py -m local -j {threads} 2> {log} && "
        "zcat {params.rundir}/results/variants/variants.vcf.gz "
        " | awk 'BEGIN {{ OFS = \"\t\"}} /^#/ {{ print $0 }} {{if($7==\"PASS\") print $0 }}' "
        " | vt decompose -s - | vt normalize  -r {input.reference} - "
        " | bgzip > {output.vcf} && "
        " tabix -p vcf {output.vcf} && rm -rf {params.rundir} 2>> {log} "


rule gatk3_mergevcf:
    input:
        reference = reference['reference_genome'],
        haplotypecaller = "{}/variants/haplotypecaller/{}.haplotypecaller-germline-normalized.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
        strelka = "{}/variants/{}-strelka-germline.passed.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
    output:
        "{}/variants/{}-all.germline.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
    threads: params['gatk3']['threads']
    log:
        "{}/logs/variants/{}.combine-germline.log".format(outdir, NORMAL_CAPTURE_STR)
    shell:
        "source activate gatk_3 && "
        "gatk3 -T CombineVariants "
        " -R {input.reference} "
        " --variant:haplotypecaller {input.haplotypecaller} "
        " --variant:strelka {input.strelka} "
        " -genotypeMergeOptions PRIORITIZE "
        " -priority haplotypecaller,strelka "
        " | bgzip > {output}  && "
        " tabix -p vcf {output} 2> {log} "
    

rule germline_generateIGVnav:
    input:
        vcf = "{}/variants/{}-all.germline.vep.vcf".format(outdir, NORMAL_CAPTURE_STR),
        oncokb = reference['oncokb']
    output:
        "{}/{}-igvnav-input.txt".format(outdir, NORMAL_CAPTURE_STR)
    params:
        vcftype = "germline"
    shell:
        "generateIGVnavInput.py {input.vcf} {input.oncokb} {params.vcftype} --output {output} "
