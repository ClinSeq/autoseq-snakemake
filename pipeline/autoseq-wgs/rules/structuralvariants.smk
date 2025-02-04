

rule svcaller_run:
    input:
        bam = outdir + "/bams/{sample}_markdups.bam",
        reference = reference["reference_genome"]
    output:
        gtf = outdir + "/svs/svcaller/{sample}-{events}.gtf"
    params: 
        tmpdir = params['scratch']
    threads: params['svcaller']['threads']
    log:
        outdir + "/logs/svs/svcaller-{sample}-{events}.log"
    shell:
        "touch {output.gtf}"


rule gridss_svcalling_normal:
    input:
        normal_bam = capture_to_results[NORMAL_CAPTURE].bamfile,
        reference = reference["bwaIndex"]
    output:
        assembly_bam = "{}/svs/gridss/{}-assembly.bam".format(outdir, NORMAL_CAPTURE_STR),
        vcf = "{}/svs/gridss/{}-gridss.vcf".format(outdir, NORMAL_CAPTURE_STR),
        svbam = "{}/svs/gridss/{}-gridss.sv.bam".format(outdir, NORMAL_CAPTURE_STR)
    params:
        jvmheap = '10g',
        basename = os.path.basename(capture_to_results[NORMAL_CAPTURE].bamfile),
        workdir = directory("{}/svs/gridss/{}/".format(outdir, NORMAL_CAPTURE_STR))
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
        " --output {output.vcf}.gz {input.normal_bam} 2> {log} && "
        "mv {params.workdir}{params.basename}.gridss.working/*sv.bam {output.svbam} && "
        " samtools index {output.svbam} && "
        " gzip -d {output.vcf}.gz  && "
        "rm -rf {params.workdir} "


rule gridss_svcalling_somatic:
    input:
        normal_bam = capture_to_results[NORMAL_CAPTURE].bamfile,
        tumor_bam = capture_to_results[CANCER_CAPTURE].bamfile,
        reference = reference["bwaIndex"]
    output:
        assembly_bam = "{}/svs/gridss/{}-{}-assembly.bam".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        vcf = "{}/svs/gridss/{}-{}-gridss.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        svbam = "{}/svs/gridss/{}-gridss.sv.bam".format(outdir, CANCER_CAPTURE_STR)
    params:
        jvmheap = '10g',
        basename = os.path.basename(capture_to_results[CANCER_CAPTURE].bamfile),
        workdir = directory("{}/svs/gridss/{}/".format(outdir, CANCER_CAPTURE_STR))
    threads: params['gridss']['threads']
    container: containers['gridss']
    log:
        outdir + "/logs/svs/gridss-{}-{}.log".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "source activate gridss-env && "
        "gridss --reference {input.reference} "
        " --jvmheap {params.jvmheap} "
        " --jar $GRIDSS_JAR "
        " -c $GRIDSS_SCRIPT/gridss.properties "
        " --assembly {output.assembly_bam} "
        " --threads {threads} --steps  ALL "
        " --workingdir {params.workdir} "
        " --output {output.vcf} {input.normal_bam} {input.tumor_bam} 2> {log} && "
        "mv {params.workdir}{params.basename}.gridss.working/*sv.bam {output.svbam} && "
        " samtools index {output.svbam} && "
        "rm -rf {params.workdir} "


rule gridss_somatic_filter:
    input:
        vcf = "{}/svs/gridss/{}-{}-gridss.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    output:
        vcf = "{}/svs/gridss/{}-{}-gridss.filtered.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        pondir = reference["pondir"],
        script_dir = os.environ.get('GRIDSS_SCRIPT')
    threads: params["gridss_filter"]["threads"]
    container: containers['gridss']
    log:
        outdir + "/logs/svs/gridss-somatic-{}-{}.log".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "source activate gridss-env && "
        "Rscript $GRIDSS_SCRIPT/gridss_somatic_filter -p {params.pondir} "
        " -i {input.vcf} "
        " -o {output.vcf} "
        " -s $GRIDSS_SCRIPT 2> {log} && "
        " bgzip -d {output.vcf}.bgz "


rule gridss_svannotation:
    input:
        somatic_vcf = "{}/svs/gridss/{}-{}-gridss.filtered.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        normal_vcf = "{}/svs/gridss/{}-gridss.vcf".format(outdir, NORMAL_CAPTURE_STR)
    output:
        somatic_vcf = "{}/svs/gridss/{}-{}-gridss.filtered.svannotated.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        normal_vcf = "{}/svs/gridss/{}-gridss.svannotated.vcf".format(outdir, NORMAL_CAPTURE_STR),
    threads: params["gridss_filter"]["threads"]
    container: containers['gridss']
    log:
        outdir + "/logs/svs/gridss-svannotation-{}-{}.log".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "source activate gridss-env && "
        "gridss_svannotate.R -v {input.somatic_vcf} -o {output.somatic_vcf} 2> {log} && "
        "gridss_svannotate.R -v {input.normal_vcf} -o {output.normal_vcf} 2>> {log} "


rule gridss_extract_readnames:
    input:
        somatic_vcf = "{}/svs/gridss/{}-{}-gridss.filtered.svannotated.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        tumor_svbam = "{}/svs/gridss/{}-gridss.sv.bam".format(outdir, CANCER_CAPTURE_STR),
        normal_vcf = "{}/svs/gridss/{}-gridss.svannotated.vcf".format(outdir, NORMAL_CAPTURE_STR),
        normal_svbam = "{}/svs/gridss/{}-gridss.sv.bam".format(outdir, NORMAL_CAPTURE_STR)
    output:
        t_readnames = "{}/svs/gridss/{}-readnames.txt".format(outdir, CANCER_CAPTURE_STR),
        n_readnames = "{}/svs/gridss/{}-readnames.txt".format(outdir, NORMAL_CAPTURE_STR)
    threads: params["gridss_filter"]["threads"]
    log:
        outdir + "/logs/svs/gridss-extract-evidence-readnames-{}-{}.log".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "generate_evidence_bam.py --vcf {input.somatic_vcf}"
        "   --bam {input.tumor_svbam}  "
        "   --readnames {output.t_readnames} 2> {log} && "
        "generate_evidence_bam.py --vcf {input.normal_vcf}"
        "   --bam {input.normal_svbam}  "
        "   --readnames {output.n_readnames} 2>> {log} "


rule gridss_evidence_bam_tumor:
    input:
        somatic_vcf = "{}/svs/gridss/{}-{}-gridss.filtered.svannotated.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        tumor_svbam = "{}/svs/gridss/{}-gridss.sv.bam".format(outdir, CANCER_CAPTURE_STR),
        readnames = "{}/svs/gridss/{}-readnames.txt".format(outdir, CANCER_CAPTURE_STR),
    output:
        tumor_bam = "{}/svs/gridss/{}-gridss.evidence.bam".format(outdir, CANCER_CAPTURE_STR)
    threads: 8
    container: containers['gridss']
    log:
        outdir + "/logs/svs/gridss-evidence-{}.log".format(CANCER_CAPTURE_STR)
    shell:
        "source activate gridss-env && "
        " samtools view -@ {threads} -hb -N {input.readnames} "
        "   -o {output.tumor_bam} {input.tumor_svbam} 2> {log} &&    "
        " samtools index -@ {threads} {output.tumor_bam} 2>> {log}"


rule gridss_evidence_bam_normal:
    input:
        normal_vcf = "{}/svs/gridss/{}-gridss.svannotated.vcf".format(outdir, NORMAL_CAPTURE_STR),
        normal_svbam = "{}/svs/gridss/{}-gridss.sv.bam".format(outdir, NORMAL_CAPTURE_STR),
        readnames = "{}/svs/gridss/{}-readnames.txt".format(outdir, NORMAL_CAPTURE_STR),
    output:
        normal_bam = "{}/svs/gridss/{}-gridss.evidence.bam".format(outdir, NORMAL_CAPTURE_STR)
    threads: 8
    container: containers['gridss']
    log:
        outdir + "/logs/svs/gridss-evidence-{}.log".format(NORMAL_CAPTURE_STR)
    shell:
        "source activate gridss-env && "
        " samtools view -@ {threads} -hb -N {input.readnames} "
        "   -o {output.normal_bam} {input.normal_svbam} 2> {log} &&    "
        " samtools index -@ {threads} {output.normal_bam} 2>> {log}"


rule generateIGVnavInput_gridss:
    input:
        somatic_vcf = "{}/svs/gridss/{}-{}-gridss.filtered.svannotated.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        normal_vcf = "{}/svs/gridss/{}-gridss.svannotated.vcf".format(outdir, NORMAL_CAPTURE_STR),
        tumor_bam = "{}/svs/gridss/{}-gridss.evidence.bam".format(outdir, CANCER_CAPTURE_STR),
        normal_bam = "{}/svs/gridss/{}-gridss.evidence.bam".format(outdir, NORMAL_CAPTURE_STR)
    output:
        somatic_mut = "{}/svs/igv/{}_somatic_pass_gridss.mut".format(outdir, CANCER_CAPTURE_STR),
        normal_mut = "{}/svs/igv/{}_normal_pass_gridss.mut".format(outdir, NORMAL_CAPTURE_STR)
    params:
        nprefix = outdir + "/svs/igv/{}".format(NORMAL_CAPTURE_STR),
        tprefix = outdir + "/svs/igv/{}".format(CANCER_CAPTURE_STR),
        sdid = "-".join(CANCER_CAPTURE_STR.split("-")[1:3])
    shell:
        "generateIGVnavInput_SV.py --input {input.normal_vcf} "
                " --sdid {params.sdid} --tool gridss " 
                " --vcftype normal --output {params.nprefix} && "
        "generateIGVnavInput_SV.py --input {input.somatic_vcf} "
                " --sdid {params.sdid} --tool gridss " 
                " --vcftype somatic --output {params.tprefix} "


rule annotate_generateIGVnavInput:
    input:
        genes = reference["genes_bed"],
        targets = reference['sv_filter'],
        gridss_somatic = "{}/svs/igv/{}_somatic_pass_gridss.mut".format(outdir, CANCER_CAPTURE_STR),
        gridss_normal = "{}/svs/igv/{}_normal_pass_gridss.mut".format(outdir, NORMAL_CAPTURE_STR),
        cgcann = reference['cgcann'],
        exons = reference['exons_gtf']
    output:
        "{}/svs/igv/{}-{}-sv-annotated.txt".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        svs_dir = "{}/svs/igv/".format(outdir),
        capture_kit_id = CANCER_CAPTURE.capture_kit_id
    shell:        
        "generateIGVnavInput_SV.py --input {params.svs_dir}  --cgc {input.cgcann} "
            " --annotBed {input.genes} --target {params.capture_kit_id} {input.targets} "
            " --exons {input.exons} --output {output} "
