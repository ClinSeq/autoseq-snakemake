

rule vep_annotation:
    input:
        reference = reference["reference_genome"],
        brca_exchange = reference["brca_exchange"],
        vep_dir = reference['vep_dir'],
        germline = "{}/variants/{}-all.germline.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
        somatic = "{}/variants/{}-{}-all.somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        germline = "{}/variants/{}-all.germline.vep.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
        somatic = "{}/variants/{}-{}-all.somatic.vep.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        germline = "{}/variants/{}-all.germline.vep.vcf".format(outdir, NORMAL_CAPTURE_STR),
        somatic = "{}/variants/{}-{}-all.somatic.vep.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: params['vep']['threads']
    container: containers['ensemblvep']
    log:
        outdir + "/logs/{}-{}-vep-annotation.log".format(CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "vep --vcf --output_file STDOUT " 
            " --pick --dir {input.vep_dir} "
            " --fasta {input.reference} --filter_common "
            " --check_existing  --total_length --allele_number "
            " --no_escape --no_stats --everything --offline "
            " --custom {input.brca_exchange},BrcaEx,vcf,exact,0,ClinicalSignificance "
            " --fork {threads} "
            " -i {input.germline} > {params.germline} 2> {log} && "
            " bgzip {params.germline} && tabix -p vcf {output.germline} && "
        "vep --vcf --output_file STDOUT " 
            " --pick --dir {input.vep_dir} "
            " --fasta {input.reference} --filter_common "
            " --check_existing  --total_length --allele_number "
            " --no_escape --no_stats --everything --offline "
            " --custom {input.brca_exchange},BrcaEx,vcf,exact,0,ClinicalSignificance "
            " --fork {threads} "
            " -i {input.somatic} > {params.somatic} 2>> {log} && "
            " bgzip {params.somatic} && tabix -p vcf {output.somatic} "

