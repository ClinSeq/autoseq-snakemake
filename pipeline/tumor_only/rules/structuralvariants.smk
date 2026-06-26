import os
import uuid 


rule svcaller_run:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam",
        reference = reference["reference_genome"]
    output:
        gtf = outdir + "/svs/svcaller/{sample}-{events}.gtf"
    params: 
        tmpdir = os.path.join(params['scratch'], 
                                "svcaller-run-{}".format(str(uuid.uuid4())))
    threads: params['svcaller']['threads']
    container: containers['svcaller']
    log:
        outdir + "/logs/svs/svcaller-{sample}-{events}.log"
    shell:
        "touch {output.gtf}"


rule gridss_extract_overlapping_fragments:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam",
        target_bed = reference['targets'][capture_name]['targets-bed-slopped20']
    output:
        bam = outdir + "/svs/gridss/{sample}-gridss-targeted.bam"
    params:
        workdir = directory("{}/svs/gridss/".format(outdir))
    threads: params['gridss']['threads']
    container: containers['gridss']
    log:
        outdir + "/logs/svs/{sample}-gridss.log"
    shell:
        "source activate gridss-env && "
        "gridss_extract_overlapping_fragments -w {params.workdir} "
        " --targetbed  {input.target_bed} -j $GRIDSS_JAR "
        " -o {output.bam} {input.bam} && "
        "samtools index {output.bam} && "
        "rm -rf  {output.bam}.gridss.working/ "


rule gridss_svcalling_tumor:
    input:
        bam = outdir + "/svs/gridss/{}-gridss-targeted.bam".format(cancer_sample[0]),
        reference = reference["bwaIndex"]
    output:
        assembly_bam = "{}/svs/gridss/{}-assembly.bam".format(outdir, CANCER_CAPTURE_STR),
        vcf = "{}/svs/gridss/{}-gridss.vcf".format(outdir, CANCER_CAPTURE_STR),
        svbam = "{}/svs/gridss/{}-gridss.sv.bam".format(outdir, CANCER_CAPTURE_STR)
    params:
        jvmheap = '10g',
        basename = "{}-gridss-targeted.bam".format(tumor_barcode),
        workdir = directory("{}/svs/gridss/{}/".format(outdir, CANCER_CAPTURE_STR))
    threads: params['gridss']['threads']
    container: containers['gridss']
    log:
        outdir + "/logs/svs/gridss-{}.log".format(NORMAL_CAPTURE_STR)
    shell:
        "source activate gridss-env && "
        "gridss --reference {input.reference} "
        " --jvmheap {params.jvmheap} "
        " --jar $GRIDSS_JAR "
        " -c $GRIDSS_SCRIPT/gridss.properties "
        " --assembly {output.assembly_bam} "
        " --threads {threads} --steps  ALL "
        " --workingdir {params.workdir} "
        " --output {output.vcf}.gz {input.bam} 2> {log} && "
        "mv {params.workdir}{params.basename}.gridss.working/*sv.bam {output.svbam} && "
        " samtools index {output.svbam} && "
        " gzip -d {output.vcf}.gz  && "
        "rm -rf {params.workdir} "


rule gridss_svannotation:
    input:
        vcf = "{}/svs/gridss/{}-gridss.vcf".format(outdir, CANCER_CAPTURE_STR)
    output:
        vcf = "{}/svs/gridss/{}-gridss.svannotated.vcf".format(outdir, CANCER_CAPTURE_STR)
    threads: params["gridss_filter"]["threads"]
    container: containers['gridss']
    log:
        outdir + "/logs/svs/gridss-svannotation-{}.log".format(CANCER_CAPTURE_STR)
    shell:
        "source activate gridss-env && "
        "gridss_svannotate.R -v {input.vcf} -o {output.vcf} 2> {log} "


rule gridss_evidence_bam:
    input:
        vcf = "{}/svs/gridss/{}-gridss.svannotated.vcf".format(outdir, CANCER_CAPTURE_STR),
        tumor_svbam = "{}/svs/gridss/{}-gridss.sv.bam".format(outdir, CANCER_CAPTURE_STR),
    output:
        tumor_bam = "{}/svs/gridss/{}-gridss.evidence.bam".format(outdir, CANCER_CAPTURE_STR)
    threads: params["gridss_filter"]["threads"]
    log:
        outdir + "/logs/svs/gridss-evidence-{}-{}.log".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "generate_evidence_bam.py --vcf {input.vcf}"
        " --bam {input.tumor_svbam} --filter-vcf "
        " --output {output.tumor_bam} 2> {log} "


rule generateIGVnavInput_gridss:
    input:
        vcf = "{}/svs/gridss/{}-gridss.svannotated.vcf".format(outdir, CANCER_CAPTURE_STR),
        tumor_bam = "{}/svs/gridss/{}-gridss.evidence.bam".format(outdir, CANCER_CAPTURE_STR)
    output:
        tumor_mut = "{}/svs/igv/{}_tumor_pass_gridss.mut".format(outdir, CANCER_CAPTURE_STR)
    params:
        tprefix = outdir + "/svs/igv/{}".format(CANCER_CAPTURE_STR),
        sdid = "-".join(CANCER_CAPTURE_STR.split("-")[1:3])
    shell:
        "generateIGVnavInput_SV.py --input {input.vcf} "
                " --sdid {params.sdid} --tool gridss " 
                " --vcftype tumor --output {params.tprefix} "


rule annotate_generateIGVnavInput:
    input:
        gridss_tumor = "{}/svs/igv/{}_tumor_pass_gridss.mut".format(outdir, CANCER_CAPTURE_STR),
        genes = reference["genes_bed"],
        targets = reference['sv_filter'],
        cgcann = reference['cgcann'],
        curation_ann = reference['curation_ann'],
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
            " -c {input.curation_ann} --annotBed {input.genes} "
            " --target {params.capture_kit_id} {input.targets} "
            " --exons {input.exons} --output {output} 2> {log} "