

rule jumblerun_cnv:
    input:
        bam = outdir + "/bams/{sample}_markdups.bam",
        reference = reference['wgs']['jumble-ref']
    output:
        cns = outdir + "/cnv/{sample}.cns",
        cnr = outdir + "/cnv/{sample}.cnr",
        seg = outdir + "/cnv/{sample}_dnacopy.seg"
    params:
        outdir = outdir + "/cnv/",
        prefix = outdir + "/cnv/{sample}_markdups.bam"
    threads: params['jumble']['threads']
    container: containers['jumble']
    log:
        outdir + "/logs/variants/{sample}-jumblerun-cnv.log"
    shell:
        "jumble-run.R -r {input.reference} "
        " -b {input.bam} " 
        " -o {params.outdir} 2> {log} && "
        "mv {params.prefix}.cnr {output.cnr} && "
        "mv {params.prefix}.cns {output.cns} && "
        "mv {params.prefix}_dnacopy.seg {output.seg} "


rule cnv_annotation:
    input:
        cns = outdir + "/cnv/{sample}.cns",
        curation_ann = reference['curation_ann'],
    output:
        outdir + "/cnv/{sample}_ann.cns"
    threads: 2
    log: outdir + "/logs/variants/{sample}-cnv-annotation.log"
    shell:
        "annotate_cnvs.py -i {input.cns} -c {input.curation_ann} -o {output} 2> {log} "


rule cnvkit_tracks:
    input:
        cns = outdir + "/cnv/{sample}.cns",
        cnr = outdir + "/cnv/{sample}.cnr"        
    output:
        profile_bedgraph = outdir + "/cnv/{sample}_profile.bedGraph",
        segments_bedgraph = outdir + "/cnv/{sample}_segments.bedGraph",
    threads: params['cnv_tracks']['threads']
    shell:
        "awk -F$'\\t' -v OFS='\\t' '$1 != \"chromosome\" {{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$6}}' "
        " {input.cnr} > {output.profile_bedgraph} "
        " && awk -F$'\\t' -v OFS='\\t' '$1 != \"chromosome\" {{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5}}' "
        " {input.cns} > {output.segments_bedgraph} "


rule liqbiocna_plot:
    input:
        capture_to_results[NORMAL_CAPTURE].svs.values(),
        capture_to_results[CANCER_CAPTURE].svs.values(),
        germline_vcf = "{}/variants/{}-all.germline.vep.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
        somatic_vcf = "{}/variants/{}-{}-all.somatic.vep.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        tumor_cns = capture_to_results[CANCER_CAPTURE].cns,
        tumor_cnr = capture_to_results[CANCER_CAPTURE].cnr,
        normal_cns = capture_to_results[NORMAL_CAPTURE].cns,
        normal_cnr = capture_to_results[NORMAL_CAPTURE].cnr,
        vcf_add_sample = "{}/variants/haplotypecaller/{}-{}.haplotypecaller-joint-calling.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        gene_track = reference['gene_track']
    output:
        liqbiocna_png = "{}/qc/{}-{}-liqbio-cna.png".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        normal_liqbiocna = "{}/qc/{}-liqbio-cna.png".format(outdir, NORMAL_CAPTURE_STR),
        cna_json = "{}/variants/{}-{}-liqbio-cna.json".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        purity_json = "{}/qc/{}-{}-liqbio-purity.json".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        normal_del = capture_to_results[NORMAL_CAPTURE].svs['DEL'],
        normal_dup = capture_to_results[NORMAL_CAPTURE].svs['DUP'],
        normal_inv = capture_to_results[NORMAL_CAPTURE].svs['INV'],
        normal_tra = capture_to_results[NORMAL_CAPTURE].svs['TRA'],
        tumor_del = capture_to_results[CANCER_CAPTURE].svs['DEL'],
        tumor_dup = capture_to_results[CANCER_CAPTURE].svs['DUP'],
        tumor_inv = capture_to_results[CANCER_CAPTURE].svs['INV'],
        tumor_tra = capture_to_results[CANCER_CAPTURE].svs['TRA'],
        purecn_mock = outdir + "/qc/purecn_mock.csv"
    threads: params['liqbiocna']['threads']
    container: containers['franken']
    log:
        outdir + "/logs/{}-{}_liqbiocna_plot.log".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "source activate franken && " 
        "touch {params.purecn_mock} && "
        "liqbioCNA_Interactive_plots_WGS.R  --tumor_cnr {input.tumor_cnr} "
                    "  --tumor_cns {input.tumor_cns} "
                    "  --normal_cnr {input.normal_cnr} "
                    "  --normal_cns {input.normal_cns} "
                    "  --het_snps_vcf {input.vcf_add_sample} "
                    "  --purecn_csv {params.purecn_mock} "
                    "  --purecn_genes_csv {params.purecn_mock} "
                    "  --purecn_loh_csv {params.purecn_mock} "
                    "  --purecn_variants_csv {params.purecn_mock} "
                    "  --svcaller_T_DEL {params.tumor_del} "
                    "  --svcaller_T_DUP {params.tumor_dup} "
                    "  --svcaller_T_INV {params.tumor_inv} "
                    "  --svcaller_T_TRA {params.tumor_tra} "
                    "  --svcaller_N_DEL {params.normal_del} "
                    "  --svcaller_N_DUP {params.normal_dup} "
                    "  --svcaller_N_INV {params.normal_inv} "
                    "  --svcaller_N_TRA {params.normal_tra} "
                    "  --germline_mut_vcf {input.germline_vcf} "
                    "  --somatic_mut_vcf {input.somatic_vcf} "
                    "  --plot_png {output.liqbiocna_png} "
                    "  --plot_png_normal {output.normal_liqbiocna} "
                    "  --cna_json {output.cna_json} "
                    "  --purity_json {output.purity_json} "
                    "  --gene_track {input.gene_track} && "
        " source deactivate "


rule franken_plot:
    input:
        capture_to_results[NORMAL_CAPTURE].svs.values(),
        capture_to_results[CANCER_CAPTURE].svs.values(),
        germline_vcf = "{}/variants/{}-{}.germline_variants_with_taf.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        somatic_vcf = "{}/variants/{}-{}-all.somatic.vep.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        tumor_cns = capture_to_results[CANCER_CAPTURE].cns,
        tumor_cnr = capture_to_results[CANCER_CAPTURE].cnr,
        normal_cns = capture_to_results[NORMAL_CAPTURE].cns,
        normal_cnr = capture_to_results[NORMAL_CAPTURE].cnr,
        vcf_add_sample = "{}/variants/haplotypecaller/{}-{}.haplotypecaller-joint-calling.vcf.gz".format(
                          outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        gene_track = reference['gene_track'],
        dpyd_json_T = outdir + "/variants/{}.typeDPYD.json".format(tumor_barcode),
        dpyd_json_N = outdir + "/variants/{}.typeDPYD.json".format(normal_barcode),
        dpyd_csv_T = outdir + "/variants/{}.typeDPYD.csv".format(tumor_barcode),
        dpyd_csv_N = outdir + "/variants/{}.typeDPYD.csv".format(normal_barcode)
    output:
        frankenplot = "{}/qc/{}-{}-frankenplot.html".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        normal_del = capture_to_results[NORMAL_CAPTURE].svs['DEL'],
        normal_dup = capture_to_results[NORMAL_CAPTURE].svs['DUP'],
        normal_inv = capture_to_results[NORMAL_CAPTURE].svs['INV'],
        normal_tra = capture_to_results[NORMAL_CAPTURE].svs['TRA'],
        tumor_del = capture_to_results[CANCER_CAPTURE].svs['DEL'],
        tumor_dup = capture_to_results[CANCER_CAPTURE].svs['DUP'],
        tumor_inv = capture_to_results[CANCER_CAPTURE].svs['INV'],
        tumor_tra = capture_to_results[CANCER_CAPTURE].svs['TRA'],
        frankenplot_rmd = os.environ.get('FRANKEN_RMD')
    threads: params['liqbiocna']['threads']
    container: containers['jumble']
    shell:
        "frankenscript.R  --tumor_cnr {input.tumor_cnr} "
                    "  --tumor_cns {input.tumor_cns} "
                    "  --normal_cnr {input.normal_cnr} "
                    "  --normal_cns {input.normal_cns} "
                    "  --frankenplot_Rmd {params.frankenplot_rmd} "
                    "  --het_snps_vcf {input.vcf_add_sample} "
                    "  --svcaller_T_DEL {params.tumor_del} "
                    "  --svcaller_T_DUP {params.tumor_dup} "
                    "  --svcaller_T_INV {params.tumor_inv} "
                    "  --svcaller_T_TRA {params.tumor_tra} "
                    "  --svcaller_N_DEL {params.normal_del} "
                    "  --svcaller_N_DUP {params.normal_dup} "
                    "  --svcaller_N_INV {params.normal_inv} "
                    "  --svcaller_N_TRA {params.normal_tra} "
                    "  --germline_mut_vcf {input.germline_vcf} "
                    "  --somatic_mut_vcf {input.somatic_vcf} "
                    "  --dpyd_json_T  {input.dpyd_json_T} "
                    "  --dpyd_json_N  {input.dpyd_json_N} "
                    "  --dpyd_csv_T {input.dpyd_csv_T} "
                    "  --dpyd_csv_N {input.dpyd_csv_N} "
                    "  --output {output.frankenplot} || true  && "
        "touch  {output.frankenplot} "
