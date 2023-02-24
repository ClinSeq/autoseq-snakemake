


for sample in all_clinseq_barcodes:
    if '-N-' in sample:
        normal_barcode = sample
    if '-CFDNA-' or '-T-' in sample:
        tumor_barcode = sample


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


# rule cat_normal_fastq:
#     input:
#         expand(nskewer_outdir + "/{prefix}" + ns1, prefix=nfq_prefix),
#         expand(nskewer_outdir + "/{prefix}" + ns2, prefix=nfq_prefix)
#     output:
#         fq1 = outdir + "/fastqs/{}_concatenated_1.fastq.gz".format(normal_barcode),
#         fq2 = outdir + "/fastqs/{}_concatenated_2.fastq.gz".format(normal_barcode)
#     run:
#         libirary = outdir + "/fastqs/skewer/"
#         fq1_files, fq2_files = find_fastqs(normal_barcode, libirary)
#         fq1_flist = " ".join(fq1_files)
#         fq2_flist = " ".join(fq2_files)

#         shell(
#             " cat {fq1_flist} > {output.fq1} && "
#             " cat {fq2_flist} > {output.fq2} && "
#             " rm {fq1_flist} {fq2_flist} "
#         )


# rule cat_tumor_fastq:
#     input:
#         expand(tskewer_outdir + "/{prefix}" + ns1, prefix=tfq_prefix),
#         expand(tskewer_outdir + "/{prefix}" + ns2, prefix=tfq_prefix)
#     output:
#         fq1 = outdir + "/fastqs/{}_concatenated_1.fastq.gz".format(tumor_barcode),
#         fq2 = outdir + "/fastqs/{}_concatenated_2.fastq.gz".format(tumor_barcode)
#     run:
#         libirary = outdir + "/fastqs/skewer/"
#         fq1_files, fq2_files = find_fastqs(tumor_barcode, libirary)
#         fq1_flist = " ".join(fq1_files)
#         fq2_flist = " ".join(fq2_files)

#         shell(
#             " cat {fq1_flist} > {output.fq1} && "
#             " cat {fq2_flist} > {output.fq2} && "
#             " rm {fq1_flist} {fq2_flist} "
#         )
