

rule vep_annotation:
    input:
        reference = reference["reference_genome"],
        brca_exchange = reference["brca_exchange"],
        vep_dir = reference['vep_dir'],
        germline = "{}/variants/{}-all.germline.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
        somatic = "{}/variants/{}-{}-all.somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        germline = "{}/variants/{}-all.germline.vep.vcf".format(outdir, NORMAL_CAPTURE_STR),
        somatic = "{}/variants/{}-{}-all.somatic.vep.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: params['vep']['threads']
    run:
        ensembl_vep = "source activate ensembl-vep"
        vep_cmd = "vep --vcf --output_file STDOUT " + \
                " --pick --filter_common " + \
                " --dir {} ".format(input.vep_dir) + \
                " --fasta {} ".format(input.reference) + \
                " --check_existing  --total_length --allele_number " + \
                " --no_escape --no_stats --everything --offline " + \
                " --custom {},,vcf,exact,0,ClinicalSignificance ".format(input.brca_exchange) + \
                " --fork 16 "
        
        germline_cmd = vep_cmd + " -i {}  >  {}".format(input.germline, output.germline)
        somatic_cmd = vep_cmd + " -i {}  >  {}".format(input.somatic, output.somatic)
        
        shell(" && ".join([ensembl_vep, germline_cmd, somatic_cmd]))
