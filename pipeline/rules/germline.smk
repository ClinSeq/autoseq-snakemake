

rule gatk4_haplotypecaller:
    input:
        bam = capture_to_results[NORMAL_CAPTURE].umibam,
        reference = reference['reference_genome'],
        dbsnp = reference["dbSNP"],
        interval_list = reference['targets'][get_capture_name(NORMAL_CAPTURE.capture_kit_id)]['targets-interval_list-slopped20'],
    output:
        vcf = "{}/variants/haplotypecaller/{}.haplotypecaller-germline.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
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
            " -L {input.interval_list} "
            " --dbsnp {input.dbsnp} "
            " -O {output.vcf} "


rule strelka_germline:
    input:
        bam = capture_to_results[NORMAL_CAPTURE].umibam,
        reference = reference['reference_genome'],
        call_region = reference['targets'][get_capture_name(NORMAL_CAPTURE.capture_kit_id)]['targets-bed-slopped20-gz'],
    output:
        rundir = directory("{}/variants/{}-strelka-germline"),
        vcf = "{}/variants/{}-strelka-germline/results/variants/{}-strelka.passed.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: params["strelka"]["threads"]
    log:
        "{}/logs/variants/{}.strelka-germline.log".format(outdir, NORMAL_CAPTURE_STR)
    shell:
        "configureStrelkaGermlineWorkflow.py  --bam {input.bam} "
        " --ref {input.reference} --targeted "
        " --callRegions {input.call_region} "
        " --runDir {output.rundir} && "
        " {output.rundir}/runWorkflow.py -m local -j 20 && "
        "zcat {output.rundir}/results/variants/variants.vcf.gz "
        " | awk 'BEGIN {{ OFS = \"\t\"}} /^#/ {{ print $0 }} {{if($7==\"PASS\") print $0 }}' "
        " | bgzip > {output.vcf} && "
        " tabix -p vcf {output.vcf} "


rule gatk3_mergevcf:
    input:
        reference = reference['reference_genome'],
        haplotypecaller = "{}/variants/haplotypecaller/{}.haplotypecaller-germline.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
        strelka = "{}/variants/{}-strelka-germline/results/variants/{}-strelka.passed.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        "{}/variants/{}-all.germline.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
    threads: params['gatk3']['threads']
    log:
        "{}/logs/variants/{}.combine-germline.log".format(outdir, NORMAL_CAPTURE_STR)
    shell:
        "source activate gatk_3"
        "gatk3 -T CombineVariants "
        " -R {input.reference} "
        " --variant:haplotypecaller {input.haplotypecaller} "
        " --variant:strelka {input.strelka} "
        " -genotypeMergeOptions PRIORITIZE "
        " -priority haplotypecaller,strelka "
        " | bgzip > {output}  && "
        " tabix -p vcf {output}"