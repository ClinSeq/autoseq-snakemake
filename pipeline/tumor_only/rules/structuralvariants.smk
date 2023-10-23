import os
import uuid 


rule svcaller_run:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam",
        reference = reference["reference_genome"]
    output:
        gtf = outdir + "/svs/svcaller/{sample}-{events}.gtf",
        bam = outdir + "/svs/svcaller/{sample}-{events}.bam",
    params: 
        tmpdir = params['scratch']
    threads: params['svcaller']['threads']
    container: containers['svcaller']
    log:
        outdir + "/logs/svs/svcaller-{sample}-{events}.log"
    shell:
        "source activate svcallerenv  && "
        "svcaller run-all --tmp-dir {params.tmpdir} --event-type {wildcards.events} "
        " --fasta-filename {input.reference}  "
        " --filter-event-overlap "
        " --events-gtf {output.gtf} "
        " --events-bam {output.bam} {input.bam} 2> {log} && "
        "source deactivate"


rule sveffect_predict:
    input:
        unpack(lambda wildcards: get_capture_svs(wildcards, outdir)),
        ts_regions = reference["ts_regions"],
        ar_regions = reference["ar_regions"],
        fusion_regions = reference["fusion_regions"]
    output:
        combined_bed = outdir + "/svs/svcaller/{sample}_combined.bed",
        effects_json = outdir + "/svs/svcaller/{sample}_effects.json"
    threads: params['svcaller']['threads']
    container: containers['svcaller']
    log:
        outdir + "/logs/svs/sveffect-{sample}.log"
    shell:
        "source activate svcallerenv  && "
        "sveffect make-bed --del-gtf {input.DEL} "
        " --dup-gtf {input.DUP} "
        " --inv-gtf {input.INV} " 
        " --tra-gtf {input.TRA} "
        " {output.combined_bed} 2> {log} &&  "
        "sveffect predict --ts-regions {input.ts_regions} "
        " --ar-regions {input.ar_regions} "
        " --fusion-regions {input.fusion_regions} "
        " --effects-filename {output.effects_json} {output.combined_bed} 2>> {log}&& "
        "source deactivate"


rule svcaller_merge:
    input:
        unpack(lambda wildcards: get_capture_svs(wildcards, outdir)),
        DEL_bam = outdir + "/svs/svcaller/{sample}-DEL.bam",
        DUP_bam = outdir + "/svs/svcaller/{sample}-DUP.bam",
        INV_bam = outdir + "/svs/svcaller/{sample}-INV.bam",
        TRA_bam = outdir + "/svs/svcaller/{sample}-TRA.bam"
    output:
        svs_bam = outdir + "/svs/{sample}-svs.bam",
        svs_gtf = outdir + "/svs/{sample}-svs.gtf"
    threads: 8
    log:
        outdir + "/logs/svs/svcaller-merge-{sample}.log"
    shell:
        "samtools merge -c -p {output.svs_bam} {input.DEL_bam} "
        " {input.DUP_bam} {input.INV_bam} {input.TRA_bam} && "
        "samtools index {output.svs_bam} && "
        "cat {input.DEL} {input.DUP} {input.INV} {input.TRA} "
        " > {output.svs_gtf} "


rule generateIGVnavInput_svcaller:
    input:
        capture_to_results[CANCER_CAPTURE].svs.values()
    output:
        svcaller_tumor = "{}/svs/igv/{}_svcaller.mut".format(outdir, CANCER_CAPTURE_STR)
    params:
        cancer_str = CANCER_CAPTURE_STR,
        svs_dir = outdir + "/svs/svcaller/",
        igvout = outdir + "/svs/igv/"
    log:
        outdir + "/logs/generateIGVnavInput_svcaller-{}.log".format(CANCER_CAPTURE_STR)
    shell:
        "generateIGVnavInput_SV.py --input {params.svs_dir} "
                        " --sdid {params.cancer_str} "
                        " --tool svcaller " 
                        " --vcftype somatic " 
                        " --output {params.igvout} 2> {log} "


rule annotate_generateIGVnavInput:
    input:
        svcaller_tumor = "{}/svs/igv/{}_svcaller.mut".format(outdir, CANCER_CAPTURE_STR),
        genes = reference["genes_bed"],
        targets = reference['sv_filter'],
        cgcann = reference['cgcann'],
        exons = reference['exons_gtf']
    output:
        "{}/svs/igv/{}-sv-annotated.txt".format(outdir, CANCER_CAPTURE_STR)
    params:
        svs_dir = "{}/svs/igv/".format(outdir),
        capture_kit_id = CANCER_CAPTURE.capture_kit_id
    log:
        outdir + "/logs/generateIGVnavInput_annotate-{}.log".format(CANCER_CAPTURE_STR)
    shell:        
        "generateIGVnavInput_SV.py --input {params.svs_dir} --cgc {input.cgcann} "
            " --annotBed {input.genes} --target {params.capture_kit_id} {input.targets} "
            " --exons {input.exons} --output {output} 2> {log} "