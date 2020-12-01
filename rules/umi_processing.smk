import os
import uuid 


rule fgbio_fastqtobam:
    input:
        fq1 = outdir + "/fastqs/{sample}_concatenated_1.fastq.gz",
        fq2 = outdir + "/fastqs/{sample}_concatenated_2.fastq.gz"
    output:
        bam = outdir + "/bams/{sample}_unmapped.bam"
    params:
        sample = lambda wildcards: compose_sample_str(extract_unique_capture(wildcards.sample)),
        library = lambda wildcards: parse_prep_id(wildcards.sample),
        java_options = params['fgbio']['fastqtobam']['java_options'],
        tmpdir = os.path.join(params['scratch'], 
                    "fgbio_fastqtobam-{}".format(str(uuid.uuid4())))
    threads: params['fgbio']['fastqtobam']['threads']
    log: outdir + "logs/fgbio/fastqtobam_{sample}.log"
    shell:
        """
        fgbio {params.java_options} --tmp-dir {params.tmpdir} FastqToBam 
            -i {input.fq1}  {input.fq2} 
            -o {output.bam} 
            --sample {params.sample} 
            --library {params.library} 
            -r 3M2S+T 3M2S+T 
            -s true 
            && rm -rf {params.tmpdir}
        """


rule bwa_umialignment_1:
    input:
        bam = outdir + "/bams/{sample}_unmapped.bam",
        reference_genome = reference['bwaIndex']
    output:
        bam = outdir + "/bams/{sample}_mapped-1.bam"
    params:
        java_options = params['picard']['merge_bam']['java_options'],
        tmpdir = params['scratch']
    threads: params['bwa']['threads']
    log: outdir + "logs/fgbio/bwa_umialignment_{sample}.log"
    shell:
        """
        picard SamToFastq I={input.bam} F=/dev/stdout INTERLEAVE=true TMP_DIR={params.tmpdir} 
            | bwa mem -p -t {threads} {input.reference_genome} /dev/stdin 
            | picard -Djava.io.tmpdir={params.tmpdir}  {params.java_options} MergeBamAlignment 
            UNMAPPED={input.bam} ALIGNED=/dev/stdin O={output.bam}
            R={input.reference_genome} 
            SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true 
            MAX_GAPS=-1 ORIENTATIONS=FR CREATE_INDEX=true 
            TMP_DIR={params.tmpdir} && rm -rf {params.tmpdir}"
        """