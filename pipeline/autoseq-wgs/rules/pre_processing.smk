

nfq_prefix, ns1, ns2 = get_fqwildcards(normal_barcode, libdir)
tfq_prefix, ts1, ts2 = get_fqwildcards(tumor_barcode, libdir)


nskewer_outdir = outdir + "/fastqs/skewer/" + normal_barcode
tskewer_outdir = outdir + "/fastqs/skewer/" + tumor_barcode
fastp_split = ["0001", "0002", "0003", "0004", "0005", "0006", "0007", "0008"]
nsplit = 2

rule fastp_splitfq_normal:
    input:
        fq1 = libdir + "/" + normal_barcode + "/{prefix}" + ns1,
        fq2 = libdir + "/" + normal_barcode + "/{prefix}" + ns2
    output:
        expand(outdir + "/fastqs/" + normal_barcode + "/{n}.{{prefix}}" + ns1, n = fastp_split[0:nsplit]),
        expand(outdir + "/fastqs/" + normal_barcode + "/{n}.{{prefix}}" + ns2, n = fastp_split[0:nsplit])
    wildcard_constraints:
        prefix = "|".join(nfq_prefix)
    params:
        nout = outdir + "/fastqs/" + normal_barcode + "/{prefix}",
        out_fq1 = os.path.join(outdir, "fastqs", normal_barcode, "{prefix}" + ns1),
        out_fq2 = os.path.join(outdir, "fastqs", normal_barcode, "{prefix}" + ns2),
        n_split = nsplit
    threads: nsplit
    log:
        outdir + "/logs/fastp/{prefix}-fastp.log"
    shell:
        """
        fastp --in1 {input.fq1} \\
            --in2  {input.fq2}    \\ 
            --disable_quality_filtering \\
            --disable_length_filtering   \\
            --disable_adapter_trimming   \\
            --disable_trim_poly_g  \\
            --split_prefix_digits=4 \\
            --split  {params.n_split}  \\  
            --thread {threads}      \\
            --out1 {params.out_fq1}    \\
            --out2 {params.out_fq2}  2> {log} 
        """


rule fastp_splitfq_tumor:
    input:
        fq1 = libdir + "/" + tumor_barcode + "/{prefix}" + ts1,
        fq2 = libdir + "/" + tumor_barcode + "/{prefix}" + ts2
    output:
        expand(outdir + "/fastqs/" + tumor_barcode + "/{n}.{{prefix}}" + ts1, n = fastp_split),
        expand(outdir + "/fastqs/" + tumor_barcode + "/{n}.{{prefix}}" + ts2, n = fastp_split)
    wildcard_constraints:
        prefix = "|".join(tfq_prefix)
    params:
        nout = outdir + "/fastqs/" + tumor_barcode + "/{prefix}",
        out_fq1 = os.path.join(outdir, "fastqs", tumor_barcode, "{prefix}" + ts1),
        out_fq2 = os.path.join(outdir, "fastqs", tumor_barcode, "{prefix}" + ts2),
        n_split = 8
    threads: 4
    log:
        outdir + "/logs/fastp/{prefix}-fastp.log"
    shell:
        """
        fastp --in1 {input.fq1} \\
            --in2  {input.fq2}    \\ 
            --disable_quality_filtering \\
            --disable_length_filtering   \\
            --disable_adapter_trimming   \\
            --disable_trim_poly_g  \\
            --split_prefix_digits=4 \\
            --split  {params.n_split}  \\  
            --thread {threads}      \\
            --out1 {params.out_fq1}    \\
            --out2 {params.out_fq2}  2> {log} 
        """


rule skewer_trim_pe_normal:
    input:
        fq1 = outdir + "/fastqs/" + normal_barcode + "/{n}.{prefix}" + ns1,
        fq2 = outdir + "/fastqs/" + normal_barcode + "/{n}.{prefix}" + ns2,
    output:
        out_fq1 = os.path.join(nskewer_outdir, "{n}.{prefix}" + ns1),
        out_fq2 = os.path.join(nskewer_outdir, "{n}.{prefix}" + ns2)
    wildcard_constraints:
        prefix = "|".join(nfq_prefix),
        n = "|".join(fastp_split[0:nsplit])
    params:
        tmpdir = os.path.join(params['scratch'], "skewer-{}".format(str(uuid.uuid4()))),
        nout = outdir + "/fastqs/skewer/" + normal_barcode
    threads: 4
    log:
        outdir + "/logs/skewer/{n}.{prefix}-skewer.log"
    shell:
        """
        mkdir {params.tmpdir}
        mkdir -p {params.nout}
        skewer -z -t {threads} --quiet \\
            -o {params.tmpdir}/skewer \\
             {input.fq1} {input.fq2} 
        # move files to output directory
        cp {params.tmpdir}/skewer-trimmed-pair1.fastq.gz {output.out_fq1} 
        cp {params.tmpdir}/skewer-trimmed-pair2.fastq.gz {output.out_fq2} 
        rm -rf {params.tmpdir} 2> {log} 
        """


rule skewer_trim_pe_tumor:
    input:
        fq1 = outdir + "/fastqs/" + tumor_barcode + "/{n}.{prefix}" + ts1,
        fq2 = outdir + "/fastqs/" + tumor_barcode + "/{n}.{prefix}" + ts2
    output:
        out_fq1 = os.path.join(tskewer_outdir, "{n}.{prefix}" + ts1),
        out_fq2 = os.path.join(tskewer_outdir, "{n}.{prefix}" + ts2)
    wildcard_constraints:
        prefix = "|".join(tfq_prefix),
        n = "|".join(fastp_split)
    params:
        tmpdir = os.path.join(params['scratch'], "skewer-{}".format(str(uuid.uuid4()))),
        nout = outdir + "/fastqs/skewer/" + tumor_barcode
    threads: 4
    log:
        outdir + "/logs/skewer/{n}.{prefix}-skewer.log"
    shell:
        """
        mkdir {params.tmpdir}
        mkdir -p {params.nout}
        skewer -z -t {threads} --quiet \\
            -o {params.tmpdir}/skewer \\
            {input.fq1} {input.fq2} 2> {log} 
        cp {params.tmpdir}/skewer-trimmed-pair1.fastq.gz {output.out_fq1}
        cp {params.tmpdir}/skewer-trimmed-pair2.fastq.gz {output.out_fq2}
        rm -rf {params.tmpdir} 2>> {log} 
        """
