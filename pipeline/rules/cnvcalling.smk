import os
import uuid 


rule cnvkit:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam",
        reference = lambda wildcards: get_cnvkitref(wildcards, reference)
    output:
        cns = outdir + "/cnv/{sample}.cns",
        cnr = outdir + "/cnv/{sample}.cnr"
    params:
        prefix = os.path.basename(outdir + "/bams/{sample}_nodups"),
        tmpdir = os.path.join(params['scratch'], 
                    "cnvkit-{}".format(str(uuid.uuid4())))
    threads: params['cnvkit']['threads']
    shell:
        "mkdir -p {params.tmpdir} && "
        "cnvkit.py batch {input.bam}  -r {input.reference} "
        " -d {params.tmpdir} "
        " && cp {params.tmpdir}/{params.prefix}.cns {output.cns}  "
        " && cp {params.tmpdir}/{params.prefix}.cnr {output.cnr}  "
        " && rm -r {params.tmpdir}"


# need to discuss
# rule cnvkit_fix:
#     input:
#         cns = outdir + "/cnv/{sample}.cns",
#         cnr = outdir + "/cnv/{sample}.cnr",
#         ref = lambda wildcards: get_cnvkitref(wildcards, reference)
#     output:
#         cns = outdir + "/cnv/{sample}-fixed.cns",
#         cnr = outdir + "/cnv/{sample}-fixed.cnr",
#     threads: params['cnvkit-fix']['threads']
#     shell:
#         "fix_cnvkit.py --input-cnr {input.cnr} "
#         " --input-cns {input.cns} "
#         " --input-reference {input.ref} "
#         " --output-cnr {output.cnr} "
#         " --output-cns {output.cns} "


rule cnvkit_cnstoseg:
    input:
        cns = outdir + "/cnv/{sample}.cns",
    output:
        seg = outdir + "/cnv/{sample}.seg",
    threads: params['cnstoseg']['threads']
    shell:
        "cnvkit.py export seg  -o {output.seg}  {input.cns}"


rule cnvkit_tracks:
    input:
        cns = outdir + "/cnv/{sample}.cns",
        cnr = outdir + "/cnv/{sample}.cnr"        
    output:
        profile_bedgraph = outdir + "/cnv/{sample}_profile.bedGraph",
        segments_bedgraph = outdir + "/cnv/{sample}_segments.bedGraph",
    threads: params['cnvkit_tracks']['threads']
    shell:
        "awk '$1 != \"chromosome\" {{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5}}' "
        " {input.cnr} > {output.profile_bedgraph} "
        " && awk '$1 != \"chromosome\" {{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5}}' "
        " {input.cns} > {output.segments_bedgraph} "


rule liqbiocna_plot:
    input:
        capture_to_results[NORMAL_CAPTURE].svs.values(),
        capture_to_results[CANCER_CAPTURE].svs.values(),
        germline_vcf = "{}/variants/{}-all.germline.vep.vcf".format(outdir, NORMAL_CAPTURE_STR),
        somatic_vcf = "{}/variants/{}-{}-all.somatic.vep.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        tumor_cns = capture_to_results[CANCER_CAPTURE].cns,
        tumor_cnr = capture_to_results[CANCER_CAPTURE].cnr,
        normal_cns = capture_to_results[NORMAL_CAPTURE].cns,
        normal_cnr = capture_to_results[NORMAL_CAPTURE].cnr,
        vcf_add_sample = "{}/variants/{}-and-{}.germline-variants-with-somatic-afs.vcf.gz".format(
            outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        purecn_csv = "{}/purecn/{}.csv".format(outdir, CANCER_CAPTURE_STR),
        purecn_genes_csv = "{}/purecn/{}_genes.csv".format(outdir, CANCER_CAPTURE_STR),
        purecn_variants_csv = "{}/purecn/{}_variants.csv".format(outdir, CANCER_CAPTURE_STR),
        purecn_loh_csv = "{}/purecn/{}_loh.csv".format(outdir, CANCER_CAPTURE_STR)
    output:
        liqbiocna_png = "{}/qc/{}-{}-liqbio-cna.png".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        normal_liqbiocna = "{}/qc/{}-liqbio-cna.png".format(outdir, NORMAL_CAPTURE_STR),
        cna_json = "{}/variants/{}-{}-liqbio-cna.json".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        purity_json = "{}/qc/{}-{}-liqbio-purity.json".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        normal = capture_to_results[NORMAL_CAPTURE].svs,
        tumor = capture_to_results[CANCER_CAPTURE].svs
    threads: params['liqbiocna']['threads']
    run:
        activate_liqbiocna = "source activate liqbiocna-env " 
        activate_flanken = "source activate flanken"

        running_cmd = "  --tumor_cnr {} ".format(input.tumor_cnr) + \
                    "  --tumor_cns {} ".format(input.tumor_cns) + \
                    "  --normal_cnr {} ".format(input.normal_cnr) + \
                    "  --normal_cns {} ".format(input.normal_cns)  + \
                    "  --het_snps_vcf {} ".format(input.vcf_add_sample) + \
                    "  --purecn_csv {} ".format(input.purecn_csv) + \
                    "  --purecn_genes_csv {} ".format(input.purecn_genes_csv) + \
                    "  --purecn_loh_csv {} ".format(input.purecn_loh_csv) + \
                    "  --purecn_variants_csv {} ".format(input.purecn_variants_csv) + \
                    "  --svcaller_T_DEL {} ".format(params.tumor['DEL']) + \
                    "  --svcaller_T_DUP {} ".format(params.tumor['DUP']) + \
                    "  --svcaller_T_INV {} ".format(params.tumor['INV']) + \
                    "  --svcaller_T_TRA {} ".format(params.tumor['TRA']) + \
                    "  --svcaller_N_DEL {} ".format(params.normal['DEL']) + \
                    "  --svcaller_N_DUP {} ".format(params.normal['DUP']) + \
                    "  --svcaller_N_INV {} ".format(params.normal['INV']) + \
                    "  --svcaller_N_TRA {} ".format(params.normal['TRA']) + \
                    "  --germline_mut_vcf {} ".format(input.germline_vcf) + \
                    "  --somatic_mut_vcf {} ".format(input.somatic_vcf) + \
                    "  --plot_png {} ".format(output.liqbiocna_png) + \
                    "  --plot_png_normal {} ".format(output.normal_liqbiocna) + \
                    "  --cna_json {} ".format(output.cna_json) + \
                    "  --purity_json {} ".format(output.purity_json) 
        
        deactivate_cmd = "source deactivate"

        run_static_plot = "liqbioCNA.R" + running_cmd
        run_interactive_plot = "liqbioCNA_Interactive_plots.R" + running_cmd

        cmd = " && ".join([activate_liqbiocna, run_static_plot, deactivate_cmd, activate_flanken, run_interactive_plot, deactivate_cmd])

        shell(cmd)
