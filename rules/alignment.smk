import os 


rule skewer_trim_pe:
    input:
        libdir + "/{sample}/"
    output:
        directory(outdir + "/fastqs/skewer/{sample}/")
    threads: 8
    log:
        outdir + "logs/skewer/skewer_{sample}.log"
    run:
        fq1_files, fq2_files = find_fastqs(wildcards.sample, libdir)
        fq1_abs = [normpath(x) for x in fq1_files]
        fq2_abs = [normpath(x) for x in fq2_files]
        pairs = [(fq1_abs[k], fq2_abs[k]) for k in range(len(fq1_abs))]
        
        for fq1, fq2 in pairs:
            tmpdir = os.path.join(params.scratch, "skewer-{uuid}".format(str(uuid.uuid4())))
            prefix = "{}/skewer".format(tmpdir)

            pre_fq1 = prefix + "-trimmed-pair1.fastq.gz"
            pre_fq2 = prefix + "-trimmed-pair2.fastq.gz"

            out_fq1 = os.path.join(output, os.path.basename(fq1))
            out_fq2 = os.path.join(output, os.path.basename(fq2))

            shell(
                " mkdir {tmpdir} && "
                " skewer -z -t {threads} --quiet "
                " -o {prefix} "
                " {fq1} {fq2} && "
                " cp {pre_fq1} {out_fq1} && "
                " cp {pre_fq2} {out_fq2} && "
                " rm -rf {tmpdir} "
            )


rule cat_fastq:
    input:
        outdir + "/fastqs/skewer/{sample}/"
    output:
        fq1 = outdir + "/fastqs/{sample}_concatenated_1.fastq.gz",
        fq2 = outdir + "/fastqs/{sample}_concatenated_2.fastq.gz"
    run:
        fq1_files, fq2_files = find_fastqs(wildcards.sample, input)
        fq1_flist = " ".join(fq1_files)
        fq2_flist = " ".join(fq2_files)

        shell(
            " cat {fq1_list} > {output.fq1} && "
            " cat {fq2_list} > {output.fq2} "
        )


rule bwa_mem_alignment:
    input:
