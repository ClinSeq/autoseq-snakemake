import os
import uuid 


rule svcaller_run:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam",
        reference = reference["reference_genome"]
    output:
        gtf = outdir + "/svs/{sample}-{events}.gtf",
        bam = outdir + "/svs/{sample}-{events}.bam",
    params: 
        tmpdir = params['scratch']
    threads: params['svcaller']['threads']
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
        combined_bed = outdir + "/svs/{sample}_combined.bed",
        effects_json = outdir + "/svs/{sample}_effects.json"
    threads: params['svcaller']['threads']
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


rule generateIGVnavInput_svcaller:
    input:
        capture_to_results[CANCER_CAPTURE].svs.values()
    output:
        svcaller_tumor = "{}/svs/igv/{}_svcaller.mut".format(outdir, CANCER_CAPTURE_STR)
    params:
        cancer_str = CANCER_CAPTURE_STR,
        svs_dir = outdir + "/svs/",
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
        targets = reference['sv_filter']
    output:
        "{}/svs/igv/{}-sv-annotated.txt".format(outdir, CANCER_CAPTURE_STR)
    params:
        svs_dir = "{}/svs/igv/".format(outdir),
        capture_kit_id = CANCER_CAPTURE.capture_kit_id
    log:
        outdir + "/logs/generateIGVnavInput_annotate-{}.log".format(CANCER_CAPTURE_STR)
    shell:        
        "generateIGVnavInput_SV.py --input {params.svs_dir} "
                " --annotBed {input.genes} --target {params.capture_kit_id} {input.target_bed} "
                " --output {output} 2> {log} "