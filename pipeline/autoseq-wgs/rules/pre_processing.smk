

nfq_prefix, ns1, ns2 = get_fqwildcards(normal_barcode, libdir)
tfq_prefix, ts1, ts2 = get_fqwildcards(tumor_barcode, libdir)


nskewer_outdir = outdir + "/fastqs/skewer/" + normal_barcode
tskewer_outdir = outdir + "/fastqs/skewer/" + tumor_barcode


rule skewer_trim_pe_normal:
    input:
        fq1 = libdir + "/" + normal_barcode + "/{prefix}" + ns1,
        fq2 = libdir + "/" + normal_barcode + "/{prefix}" + ns2
    output:
        out_fq1 = os.path.join(nskewer_outdir, "{prefix}" + ns1),
        out_fq2 = os.path.join(nskewer_outdir, "{prefix}" + ns2)
    wildcard_constraints:
        prefix = "|".join(nfq_prefix)
    params:
        tmpdir = os.path.join(params['scratch'], "skewer-{}".format(str(uuid.uuid4()))),
        nout = outdir + "/fastqs/skewer/" + normal_barcode
    threads: 4
    log:
        outdir + "/logs/skewer/{prefix}-skewer.log"
    shell:
        " mkdir {params.tmpdir} && "
        " mkdir -p {params.nout} && "
        " skewer -z -t {threads} --quiet "
        " -o {params.tmpdir}/skewer "
        " {input.fq1} {input.fq2} && "
        " cp {params.tmpdir}/skewer-trimmed-pair1.fastq.gz {output.out_fq1} && "
        " cp {params.tmpdir}/skewer-trimmed-pair2.fastq.gz {output.out_fq2} && "
        " rm -rf {params.tmpdir} 2> {log} "


rule skewer_trim_pe_tumor:
    input:
        fq1 = libdir + "/" + tumor_barcode + "/{prefix}" + ts1,
        fq2 = libdir + "/" + tumor_barcode + "/{prefix}" + ts2,
    output:
        out_fq1 = os.path.join(tskewer_outdir, "{prefix}" + ts1),
        out_fq2 = os.path.join(tskewer_outdir, "{prefix}" + ts2)
    wildcard_constraints:
        prefix = "|".join(tfq_prefix)
    params:
        tmpdir = os.path.join(params['scratch'], "skewer-{}".format(str(uuid.uuid4()))),
        nout = outdir + "/fastqs/skewer/" + tumor_barcode
    threads: 4
    log:
        outdir + "/logs/skewer/{prefix}-skewer.log"
    shell:
        " mkdir {params.tmpdir} && "
        " mkdir -p {params.nout} && "
        " skewer -z -t {threads} --quiet "
        " -o {params.tmpdir}/skewer "
        " {input.fq1} {input.fq2} 2> {log} && "
        " cp {params.tmpdir}/skewer-trimmed-pair1.fastq.gz {output.out_fq1} && "
        " cp {params.tmpdir}/skewer-trimmed-pair2.fastq.gz {output.out_fq2} && "
        " rm -rf {params.tmpdir} 2>> {log} "
