
rule vep_annotation_somatic:
    input:
        reference = reference["reference_genome"],
        gnomAD = config["gnomAD"],
        vep_dir = reference['vep_dir'],
        somatic = "{}/variants/{}-{}-all.somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        "{}/variants/{}-{}-all.somatic.gnomADg.vep.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: params['vep']['threads']
    shell:
        "source activate ensembl-vep && "
        "vep --vcf --output_file STDOUT " 
            " --pick --filter_common "
            " --dir {input.vep_dir} "
            " --fasta {input.reference} "
            " --check_existing  --total_length --allele_number "
            " --no_escape --no_stats --everything --offline "
            " --custom {input.gnomAD},gnomADg,vcf,exact,0,AF "
            " --fork {threads} "
            " -i {input.somatic} | bgzip > {output} && "
        " tabix -p vcf {output} "


# Create files of somatic calls with only SNPs and only non-SNPs, 
# with "SNPs" defined as variants with gnomADg_AF > 1%:
rule vep_filter_somatic:
    input:
        "{}/variants/{}-{}-all.somatic.gnomADg.vep.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        snps = "{}/variants/{}-{}-all.somatic.gnomADg.SNPs.vep.vcf.gz".format(outdir, 
                                            CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        nosnps = "{}/variants/{}-{}-all.somatic.gnomADg.noSNPs.vep.vcf.gz".format(outdir, 
                                            CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: 1
    shell:
        "source activate ensembl-vep && "
        "zcat {input} | filter_vep --filter \"gnomADg_AF > 0.01\" "
        " | bgzip > {output.snps} && "
        " tabix -p vcf {output.snps} && "
        "zcat {input} | filter_vep --filter \"gnomADg_AF <= 0.01 or not gnomADg_AF\" "
        " | bgzip > {output.nosnps} && "
        " tabix -p vcf {output.nosnps} "


rule vep_annotation_noSNPs:
    input:
        reference = reference["reference_genome"],
        brca_exchange = reference["brca_exchange"],
        vep_dir = reference['vep_dir'],
        vcf = "{}/variants/{}-{}-all.somatic.gnomADg.noSNPs.vep.vcf.gz".format(outdir, 
                                                CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        "{}/variants/{}-{}-all.somatic.gnomADg.noSNPs.brcaEx.vep.vcf.gz".format(outdir, 
                                                CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: params['vep']['threads']
    shell:
        "source activate ensembl-vep && "
        "vep --vcf --output_file STDOUT " 
            " --pick --filter_common "
            " --dir {input.vep_dir} "
            " --fasta {input.reference} "
            " --check_existing  --total_length --allele_number "
            " --no_escape --no_stats --everything --offline "
            " --custom {input.brca_exchange},BrcaEx,vcf,exact,0,ClinicalSignificance "
            " --fork {threads} "
            " -i {input.vcf} | bgzip > {output} && "
        " tabix -p vcf {output} "


rule vep_annotation_germline:
    input:
        reference = reference["reference_genome"],
        gnomAD = config["gnomAD"],
        vep_dir = reference['vep_dir'],
        germline = "{}/variants/{}-merged.germline.split_norm.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    output:
        "{}/variants/{}-merged.germline.split_norm.gnomADg.vep.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    threads: params['vep']['threads']
    shell:
        "source activate ensembl-vep && "
        "vep --vcf --output_file STDOUT " 
            " --pick --filter_common "
            " --dir {input.vep_dir} "
            " --fasta {input.reference} "
            " --check_existing  --total_length --allele_number "
            " --no_escape --no_stats --everything --offline "
            " --custom {input.gnomAD},gnomADg,vcf,exact,0,AF "
            " --fork {threads} "
            " -i {input.germline} | bgzip > {output} && "
        " tabix -p vcf {output} "


# Create files of germline calls with only SNPs and only non-SNPs, 
# with "SNPs" defined as variants with gnomADg_AF > 1%
rule vep_filter_germline:
    input:
        "{}/variants/{}-merged.germline.split_norm.gnomADg.vep.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    output:
        snps = "{}/variants/{}-merged.germline.split_norm.gnomADg.vep.SNPs.vcf.gz".format(outdir, CANCER_CAPTURE_STR),
        nosnps = "{}/variants/{}-merged.germline.split_norm.gnomADg.vep.noSNPs.vcf.gz".format(outdir, CANCER_CAPTURE_STR),
    threads: 1
    shell:
        "source activate ensembl-vep && "
        "zcat {input} | filter_vep --filter \"gnomADg_AF > 0.01\" "
        " | bgzip > {output.snps} && "
        " tabix -p vcf {output.snps} && "
        "zcat {input} | filter_vep --filter \"gnomADg_AF <= 0.01 or not gnomADg_AF\" "
        " | bgzip > {output.nosnps} && "
        " tabix -p vcf {output.nosnps} "


# Run vep with autoseq configurations, on unfiltered "germline" vcfs (needed for creating germline variant table)
rule vep_annotation_brcaex:
    input:
        reference = reference["reference_genome"],
        brca_exchange = reference["brca_exchange"],
        vep_dir = reference['vep_dir'],
        germline = "{}/variants/{}-merged.germline.split_norm.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    output:
        "{}/variants/{}-merged.germline.split_norm.brcaEx.vep.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    threads: params['vep']['threads']
    shell:
        "source activate ensembl-vep && "
        "vep --vcf --output_file STDOUT " 
            " --pick --filter_common "
            " --dir {input.vep_dir} "
            " --fasta {input.reference} "
            " --check_existing  --total_length --allele_number "
            " --no_escape --no_stats --everything --offline "
            " --custom {input.brca_exchange},BrcaEx,vcf,exact,0,ClinicalSignificance "
            " --fork {threads} "
            " -i {input.germline} | bgzip > {output} && "
        " tabix -p vcf {output} "


