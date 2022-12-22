import os
import uuid 

capture_name = get_capture_name(CANCER_CAPTURE.capture_kit_id)

# rule cnvkit:
#     input:
#         bam = outdir + "/bams/{sample}_nodups.bam",
#         reference = lambda wildcards: get_cnvkitref(wildcards, reference)
#     output:
#         cns = outdir + "/cnv/{sample}.cns",
#         cnr = outdir + "/cnv/{sample}.cnr"
#     params:
#         prefix = os.path.basename(outdir + "/bams/{sample}_nodups"),
#         tmpdir = os.path.join(params['scratch'], 
#                     "cnvkit-{}".format(str(uuid.uuid4())))
#     threads: params['cnvkit']['threads']
#     shell:
#         "mkdir -p {params.tmpdir} && "
#         "cnvkit.py batch {input.bam}  -r {input.reference} "
#         " -d {params.tmpdir} "
#         " && cp {params.tmpdir}/{params.prefix}.cns {output.cns}  "
#         " && cp {params.tmpdir}/{params.prefix}.cnr {output.cnr}  "
#         " && rm -r {params.tmpdir}"


# rule cnvkit_cnstoseg:
#     input:
#         cns = outdir + "/cnv/{sample}.cns",
#     output:
#         seg = outdir + "/cnv/{sample}_dnacopy.seg",
#     threads: params['cnstoseg']['threads']
#     shell:
#         "cnvkit.py export seg  -o {output.seg}  {input.cns}"


rule jumblerun_cnv:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam",
        reference = reference['targets'][capture_name]['jumble-ref']
    output:
        cns = outdir + "/cnv/{sample}.cns",
        cnr = outdir + "/cnv/{sample}.cnr",
        seg = outdir + "/cnv/{sample}_dnacopy.seg"
    params:
        outdir = outdir + "/cnv/"
    threads: params['jumble']['threads']
    log:
        outdir + "/logs/variants/{sample}-jumblerun-cnv.log"
    shell:
        "source activate jumble-env && "
        "jumble-run.R -r {input.reference} "
        " -b {input.bam} " 
        " -o {params.outdir} 2> {log} "


rule cnv_tracks:
    input:
        cns = outdir + "/cnv/{sample}.cns",
        cnr = outdir + "/cnv/{sample}.cnr"        
    output:
        profile_bedgraph = outdir + "/cnv/{sample}_profile.bedGraph",
        segments_bedgraph = outdir + "/cnv/{sample}_segments.bedGraph",
    threads: params['cnv_tracks']['threads']
    shell:
        "awk -F$'\\t' -v OFS='\\t' '$1 != \"chromosome\" {{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5}}' "
        " {input.cnr} > {output.profile_bedgraph} "
        " && awk -F$'\\t' -v OFS='\\t' '$1 != \"chromosome\" {{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5}}' "
        " {input.cns} > {output.segments_bedgraph} "


rule liqbiocna_plot:
    input:
        capture_to_results[CANCER_CAPTURE].svs.values(),
        germline_vcf = "{}/variants/{}-merged.germline.split_norm.brcaEx.vep.vcf.gz".format(outdir, CANCER_CAPTURE_STR),
        somatic_vcf = "{}/variants/{}-{}-all.somatic.gnomADg.noSNPs.brcaEx.vep.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        cns = capture_to_results[CANCER_CAPTURE].cns,
        cnr = capture_to_results[CANCER_CAPTURE].cnr,
        vcf_add_sample = "{}/variants/{}-merged.germline.split_norm.gnomADg.vep.SNPs.BAF.vcf.gz".format(outdir, CANCER_CAPTURE_STR),
        purecn_csv = "{}/purecn/{}.csv".format(outdir, CANCER_CAPTURE_STR),
        purecn_genes_csv = "{}/purecn/{}_genes.csv".format(outdir, CANCER_CAPTURE_STR),
        purecn_variants_csv = "{}/purecn/{}_variants.csv".format(outdir, CANCER_CAPTURE_STR),
        purecn_loh_csv = "{}/purecn/{}_loh.csv".format(outdir, CANCER_CAPTURE_STR),
        gene_track = reference['gene_track']
    output:
        liqbiocna_png = "{}/qc/{}-liqbio-cna.png".format(outdir, CANCER_CAPTURE_STR),
        normal_liqbiocna = "{}/qc/{}-normal-liqbio-cna.png".format(outdir, CANCER_CAPTURE_STR),
        cna_json = "{}/variants/{}-{}-liqbio-cna.json".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        purity_json = "{}/qc/{}-{}-liqbio-purity.json".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        tumor_del = capture_to_results[CANCER_CAPTURE].svs['DEL'],
        tumor_dup = capture_to_results[CANCER_CAPTURE].svs['DUP'],
        tumor_inv = capture_to_results[CANCER_CAPTURE].svs['INV'],
        tumor_tra = capture_to_results[CANCER_CAPTURE].svs['TRA'],
        tmp_germvcf = "{}/variants/tmp-{}-merged.germline.split_norm.brcaEx.vep.vcf".format(outdir, CANCER_CAPTURE_STR),
        mock_txt = os.path.join(outdir, "mock.txt")
    threads: params['liqbiocna']['threads']
    run:
        shell("source activate franken && " 
        "touch {params.mock_txt} && "
        "zcat {input.germline_vcf} | awk -F '\\t' -v OFS='\\t' "
        " '{{if ($1 ~ /^#/) print $0; else if ($5 != \"*\") {{$3=\".\"; print $0}}}}' > {params.tmp_germvcf} && "
        "liqbioCNA_Interactive_plots.R  --tumor_cnr {input.cnr} "
                    "  --tumor_cns {input.cns} "
                    "  --normal_cnr {input.cnr} "
                    "  --normal_cns {input.cns} "
                    "  --het_snps_vcf {input.vcf_add_sample} "
                    "  --purecn_csv {input.purecn_csv} "
                    "  --purecn_genes_csv {input.purecn_genes_csv} "
                    "  --purecn_loh_csv {input.purecn_loh_csv} "
                    "  --purecn_variants_csv {input.purecn_variants_csv} "
                    "  --svcaller_T_DEL {params.tumor_del} "
                    "  --svcaller_T_DUP {params.tumor_dup} "
                    "  --svcaller_T_INV {params.tumor_inv} "
                    "  --svcaller_T_TRA {params.tumor_tra} "
                    "  --svcaller_N_DEL {params.mock_txt} "
                    "  --svcaller_N_DUP {params.mock_txt} "
                    "  --svcaller_N_INV {params.mock_txt} "
                    "  --svcaller_N_TRA {params.mock_txt} "
                    "  --germline_mut_vcf {params.tmp_germvcf} "
                    "  --somatic_mut_vcf {input.somatic_vcf} "
                    "  --plot_png {output.liqbiocna_png} "
                    "  --plot_png_normal {output.normal_liqbiocna} "
                    "  --cna_json {output.cna_json} "
                    "  --purity_json {output.purity_json} "
                    "  --gene_track {input.gene_track} && "
        "rm {params.tmp_germvcf} && "
        " source deactivate ")


rule franken_plot:
    input:
        capture_to_results[CANCER_CAPTURE].svs.values(),
        somatic_vcf = "{}/variants/{}-{}-all.somatic.gnomADg.noSNPs.brcaEx.vep.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        cns = capture_to_results[CANCER_CAPTURE].cns,
        cnr = capture_to_results[CANCER_CAPTURE].cnr,
        vcf_add_sample = "{}/variants/{}-merged.germline.split_norm.gnomADg.vep.SNPs.BAF.vcf.gz".format(outdir, CANCER_CAPTURE_STR),
        purecn_csv = "{}/purecn/{}.csv".format(outdir, CANCER_CAPTURE_STR),
        purecn_genes_csv = "{}/purecn/{}_genes.csv".format(outdir, CANCER_CAPTURE_STR),
        purecn_variants_csv = "{}/purecn/{}_variants.csv".format(outdir, CANCER_CAPTURE_STR),
        purecn_loh_csv = "{}/purecn/{}_loh.csv".format(outdir, CANCER_CAPTURE_STR)
    output:
        frankenplot = "{}/qc/{}-frankenplot.html".format(outdir, CANCER_CAPTURE_STR)
    params:
        tumor_del = capture_to_results[CANCER_CAPTURE].svs['DEL'],
        tumor_dup = capture_to_results[CANCER_CAPTURE].svs['DUP'],
        tumor_inv = capture_to_results[CANCER_CAPTURE].svs['INV'],
        tumor_tra = capture_to_results[CANCER_CAPTURE].svs['TRA'],
        frankenplot_rmd = os.environ.get('FRANKEN_RMD')
    threads: params['liqbiocna']['threads']
    run:
        shell("source activate jumble-env && " 
        "frankenscript.R  --tumor_cnr {input.cnr} "
                    "  --tumor_cns {input.cns} "
                    "  --frankenplot_Rmd {params.frankenplot_rmd} "
                    "  --het_snps_vcf {input.vcf_add_sample} "
                    "  --purecn_csv {input.purecn_csv} "
                    "  --purecn_genes_csv {input.purecn_genes_csv} "
                    "  --purecn_loh_csv {input.purecn_loh_csv} "
                    "  --purecn_variants_csv {input.purecn_variants_csv} "
                    "  --svcaller_T_DEL {params.tumor_del} "
                    "  --svcaller_T_DUP {params.tumor_dup} "
                    "  --svcaller_T_INV {params.tumor_inv} "
                    "  --svcaller_T_TRA {params.tumor_tra} "
                    "  --somatic_mut_vcf {input.somatic_vcf} "
                    "  --output {output.frankenplot} && "
        " source deactivate ")
