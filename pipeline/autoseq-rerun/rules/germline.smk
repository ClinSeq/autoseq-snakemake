
normalBam = capture_to_results[NORMAL_CAPTURE].umibam if umi else capture_to_results[NORMAL_CAPTURE].bamfile

rule gatk4_haplotypecaller:
    input:
        bam = normalBam,
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


rule germline_generateIGVnav:
    input:
        vcf = "{}/variants/{}-{}.germline_variants_with_taf.vcf".format(outdir, 
                NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        oncokb = reference['oncokb'],
        cgcann = reference["cgcann"],
        cancer_hotspot_snv = reference['cancer_hotspot_snv'],
        cancer_hotspot_indel = reference['cancer_hotspot_indel']
    output:
        "{}/{}-igvnav-input.txt".format(outdir, NORMAL_CAPTURE_STR)
    params:
        vcftype = "germline"
    shell:
        "generateIGVnavInput.py {input.vcf} {input.oncokb} {params.vcftype} "
        " --hotspot-snv {input.cancer_hotspot_snv} "
        " --hotspot-indel {input.cancer_hotspot_indel} "
        " --cgc {input.cgcann} --output {output} "