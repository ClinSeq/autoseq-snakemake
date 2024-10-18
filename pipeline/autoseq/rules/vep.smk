

rule vep_annotation:
    input:
        reference = reference["reference_genome"],
        brca_exchange = reference["brca_exchange"],
        vep_dir = reference['vep_dir'],
        germline = "{}/variants/haplotypecaller/{}.haplotypecaller-germline.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
        somatic = "{}/variants/{}-{}-all.somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        germline = "{}/variants/{}-all.germline.vep.vcf".format(outdir, NORMAL_CAPTURE_STR),
        somatic = "{}/variants/{}-{}-all.somatic.vep.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: params['vep']['threads']
    log:
        vep_log = outdir + "/logs/vep_annotation_{}.log".format(CANCER_CAPTURE_STR)
    container: containers['ensemblvep_v112']
    shell:
        "vep --vcf --output_file STDOUT " 
            " --pick --dir {input.vep_dir} "
            " --fasta {input.reference} "
            " --check_existing  --total_length --allele_number "
            " --no_escape --no_stats --everything --offline "
            " --custom {input.brca_exchange},BrcaEx,vcf,exact,0,ClinicalSignificance "
            " --fork {threads} --filter_common  --format vcf "
            " -i {input.germline} > {output.germline} 2> {log.vep_log} &&"
        "vep --vcf --output_file STDOUT " 
            " --pick --dir {input.vep_dir} "
            " --fasta {input.reference} "
            " --check_existing  --total_length --allele_number "
            " --no_escape --no_stats --everything --offline "
            " --custom {input.brca_exchange},BrcaEx,vcf,exact,0,ClinicalSignificance "
            " --fork {threads} --format vcf "
            " -i {input.somatic} > {output.somatic} 2>> {log.vep_log}"

