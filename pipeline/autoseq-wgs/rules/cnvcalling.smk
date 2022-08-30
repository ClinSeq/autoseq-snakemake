

rule cnvkit:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam",
        reference_genome = reference['reference_genome']
    output:
        cns = outdir + "/cnv/{sample}.cns",
        cnr = outdir + "/cnv/{sample}.cnr"
    params:
        prefix = os.path.basename(outdir + "/bams/{sample}_nodups"),
        tmpdir = os.path.join(params['scratch'], 
                    "cnvkit-{}".format(str(uuid.uuid4())))
    threads: params['cnvkit']['threads']
    log:
        outdir + "/logs/{sample}_cnvkit.log"
    shell:
        "mkdir -p {params.tmpdir} && "
        "cnvkit.py batch {input.bam}  -p {threads} --fasta {input.reference_genome} "
        " -d {params.tmpdir} 2> {log} "
        " && cp {params.tmpdir}/{params.prefix}.cns {output.cns}  "
        " && cp {params.tmpdir}/{params.prefix}.cnr {output.cnr}  "
        " && rm -r {params.tmpdir}"


rule cnvkit_cnstoseg:
    input:
        cns = outdir + "/cnv/{sample}.cns",
    output:
        seg = outdir + "/cnv/{sample}.seg",
    threads: params['cnstoseg']['threads']
    log:
        outdir + "/logs/{sample}_cnvkit_cnstoseg.log"
    shell:
        "cnvkit.py export seg  -o {output.seg}  {input.cns} 2> {log} "


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
        purecn_loh_csv = "{}/purecn/{}_loh.csv".format(outdir, CANCER_CAPTURE_STR),
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
    threads: params['liqbiocna']['threads']
    log:
        outdir + "/logs/{}-{}_liqbiocna_plot.log".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    run:
        shell("source activate franken && " 
        "liqbioCNA_Interactive_plots.R  --tumor_cnr {input.tumor_cnr} "
                    "  --tumor_cns {input.tumor_cns} "
                    "  --normal_cnr {input.normal_cnr} "
                    "  --normal_cns {input.normal_cns} "
                    "  --het_snps_vcf {input.vcf_add_sample} "
                    "  --purecn_csv {input.purecn_csv} "
                    "  --purecn_genes_csv {input.purecn_genes_csv} "
                    "  --purecn_loh_csv {input.purecn_loh_csv} "
                    "  --purecn_variants_csv {input.purecn_variants_csv} "
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
        " source deactivate ")
