

sdid = samples["sdid"]
nf_input = {
    "normal_bam": outdir + "/bams/{}_markdups.bam".format(normal_barcode),
    "tumor_bam": outdir + "/bams/{}_markdups.bam".format(tumor_barcode),
}

if rna_barcode:
    include: "rna-alignment.smk"
    nf_input["rna_bam"] = outdir + "/bams/{}.dedup.bam".format(rna_barcode)


rule prepare_samplesheet:
    input:
        **nf_input
    output:
        outdir + "/nf-oncoanalyser/samplesheet_{}.csv".format(sdid)
    params:
        group_id = group_id,
        subject_id = sdid,
        normal_id = normal_barcode,
        tumor_id = tumor_barcode,
        rna_id = rna_barcode,
    log:
        outdir + "/logs/prepare_samplesheet_{}.log".format(sdid)
    run:
        header = "group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath"
        n_entry = f"{params.group_id},{params.subject_id},{params.normal_id},normal,dna,bam_markdups,{input.normal_bam}"
        t_entry = f"{params.group_id},{params.subject_id},{params.tumor_id},tumor,dna,bam_markdups,{input.tumor_bam}"
        shell(f"echo {header} > {output}")
        shell(f"echo {n_entry} >> {output}")
        shell(f"echo {t_entry} >> {output}")
        if "rna_bam" in input:
            r_entry = f"{params.group_id},{params.subject_id},{params.rna_id},tumor,rna,bam_markdups,{input.rna_bam}"
            shell(f"echo {r_entry} >> {output}")


rule run_oncoanalyser:
    input:
        outdir + "/nf-oncoanalyser/samplesheet_{}.csv".format(sdid)
    output:
        outdir + "/nf-oncoanalyser/{group}/orange/{tumor}.orange.pdf".format(group=group_id, tumor=tumor_barcode)
    params:
        outdir = outdir + "/nf-oncoanalyser/",
        oncoanalyser_base = "/nfs/PIPELINE/oncoanalyser",
        genome_gtf = reference['gencode_gtf'],
        nf_reference = config['nf_reference'],
    threads: 22
    log:
        outdir + "/logs/run_oncoanalyser_{}.log".format(sdid)
    run:
        conda_cmd = "source activate oncoanalyser && "
        cmd = "nextflow run {}/main.nf ".format(params.oncoanalyser_base)
        cmd += " -profile singularity  "
        cmd += " -config  {} ".format(params.nf_reference)
        cmd += "  --mode wgts --genome CustomGenome --genome_type no_alt --genome_version 37 --force_genome "
        cmd += "  --ref_data_genome_gtf {} ".format(params.genome_gtf)
        cmd += " --input {} ".format(input)
        cmd += " --outdir {} ".format(params.outdir)

        shell(f"{conda_cmd} {cmd} 2> {log}")