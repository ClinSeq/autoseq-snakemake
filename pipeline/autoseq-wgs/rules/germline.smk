bamfile = capture_to_results[NORMAL_CAPTURE].bamfile

haplotype_vcf_prefix = "{}/variants/haplotypecaller/{}.haplotypecaller-germline".format(outdir, NORMAL_CAPTURE_STR)
haplotype_log_prefix = "{}/logs/variants/haplotypecaller/{}.haplotypecaller-germline".format(outdir, NORMAL_CAPTURE_STR)

rule gatk4_haplotypecaller:
    input:
        bam = bamfile,
        reference = reference['reference_genome'],
        dbsnp = reference["dbSNP"],
        interval_list = intervals_dir + "human_g1k_v37_decoy.{suf}.interval_list"
    output:
        vcf = haplotype_vcf_prefix + ".{suf}.vcf.gz",
    wildcard_constraints:
        suf = "|".join(suffix)
    params:
        java_options = params["gatk4"]["haplotypecaller"]["java_options"]
    threads: params["gatk4"]["threads"]
    log:
        haplotype_log_prefix + ".{suf}.log"
    shell:
        "gatk --java-options '{params.java_options}' "
        " HaplotypeCaller   "
        " -R {input.reference}  "
        " -I {input.bam}  "
        " -L {input.interval_list} "
        " --dbsnp {input.dbsnp} "
        " -O {output.vcf} 2> {log}  "


rule haplotypecaller_vcfmerge:
    input:
        expand(haplotype_vcf_prefix + ".{suf}.vcf.gz", suf=suffix)
    output:
        hp_vcf = haplotype_vcf_prefix + ".vcf.gz",
        vcf = "{}/variants/{}-all.germline.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
    threads: params['bcftools']['threads']
    log:
        "{}/logs/variants/{}-haplotypecaller-germline-merge.log".format(outdir, NORMAL_CAPTURE_STR)
    shell:
        "bcftools concat -D --threads {threads} -a "
        " -O z {input} > {output.hp_vcf} 2> {log} && "
        " cp {output.hp_vcf} {output.vcf} && "
        " tabix -p vcf {output.vcf} "
    

rule germline_generateIGVnav:
    input:
        vcf = "{}/variants/{}-all.germline.vep.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
        oncokb = reference['oncokb'],
        cgcann = reference['cgcann']
    output:
        "{}/{}-igvnav-input.txt".format(outdir, NORMAL_CAPTURE_STR)
    params:
        vcftype = "germline"
    shell:
        "generateIGVnavInput.py {input.vcf} {input.oncokb} {params.vcftype} --wgs --cgc {input.cgcann} --output {output} "
