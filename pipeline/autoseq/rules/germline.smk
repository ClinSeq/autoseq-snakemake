
bamfile = capture_to_results[NORMAL_CAPTURE].umibam if umi else capture_to_results[NORMAL_CAPTURE].bamfile

rule gatk4_haplotypecaller:
    input:
        bam = bamfile,
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
            " -O {output.vcf}  2> {log} "


# rule strelka_germline:
#     input:
#         bam = bamfile,
#         reference = reference['reference_genome'],
#         call_region = reference['targets'][get_capture_name(NORMAL_CAPTURE.capture_kit_id)]['targets-bed-slopped20-gz'],
#     output:
#         vcf = "{}/variants/{}-strelka-germline.passed.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
#     params:
#         rundir = directory("{}/variants/{}-strelka-germline".format(outdir, NORMAL_CAPTURE_STR))
#     threads: params["strelka"]["threads"]
#     log:
#         "{}/logs/variants/{}.strelka-germline.log".format(outdir, NORMAL_CAPTURE_STR)
#     shell:
#         "source activate gatk_3 && "
#         "configureStrelkaGermlineWorkflow.py  --bam {input.bam} "
#         " --ref {input.reference} --targeted "
#         " --callRegions {input.call_region} "
#         " --runDir {params.rundir} && "
#         " {params.rundir}/runWorkflow.py -m local -j {threads} && "
#         "zcat {params.rundir}/results/variants/variants.vcf.gz "
#         " | awk 'BEGIN {{ OFS = \"\t\"}} /^#/ {{ print $0 }} {{if($7==\"PASS\") print $0 }}' "
#         " | vt decompose -s - | vt normalize  -r {input.reference} - "
#         " | bgzip > {output.vcf} && "
#         " tabix -p vcf {output.vcf} && rm -rf {params.rundir} 2> {log} "


# rule gatk3_mergevcf:
#     input:
#         reference = reference['reference_genome'],
#         haplotypecaller = "{}/variants/haplotypecaller/{}.haplotypecaller-germline-normalized.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
#         strelka = "{}/variants/{}-strelka-germline.passed.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
#     output:
#         "{}/variants/{}-all.germline.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
#     threads: params['gatk3']['threads']
#     log:
#         "{}/logs/variants/{}.combine-germline.log".format(outdir, NORMAL_CAPTURE_STR)
#     shell:
#         "source activate gatk_3 && "
#         "gatk3 -T CombineVariants "
#         " -R {input.reference} "
#         " --variant:haplotypecaller {input.haplotypecaller} "
#         " --variant:strelka {input.strelka} "
#         " -genotypeMergeOptions PRIORITIZE "
#         " -priority haplotypecaller,strelka "
#         " | bgzip > {output}  && "
#         " tabix -p vcf {output} 2> {log} "
    

rule germline_generateIGVnav:
    input:
        vcf = "{}/variants/{}-all.germline.vep.vcf".format(outdir, NORMAL_CAPTURE_STR),
        oncokb = reference['oncokb'],
        cgcann = reference["cgcann"]
    output:
        "{}/{}-igvnav-input.txt".format(outdir, NORMAL_CAPTURE_STR)
    params:
        vcftype = "germline"
    shell:
        "generateIGVnavInput.py {input.vcf} {input.oncokb} {params.vcftype} "
        " --cgc {input.cgcann} --output {output} "
