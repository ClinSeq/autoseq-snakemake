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
    log:
        outdir + "/logs/svs/svcaller-{sample}-{events}.log"
    shell:
        "source activate svcallerenv  && "
        "svcaller run-all --tmp-dir {params.tmpdir} --event-type {wildcards.events} "
        " --fasta-filename {input.reference}  "
        " --filter-event-overlap "
        " --events-gtf {output.gtf} "
        " --events-bam {output.bam} {input.bam} && "
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
    log:
        outdir + "/logs/svs/sveffect-{sample}.log"
    shell:
        "source activate svcallerenv  && "
        "sveffect make-bed --del-gtf {input.DEL} "
        " --dup-gtf {input.DUP} "
        " --inv-gtf {input.INV} " 
        " --tra-gtf {input.TRA} "
        " {output.combined_bed} &&  "
        "sveffect predict --ts-regions {input.ts_regions} "
        " --ar-regions {input.ar_regions} "
        " --fusion-regions {input.fusion_regions} "
        " --effects-filename {output.effects_json} {output.combined_bed} && "
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
        capture_to_results[NORMAL_CAPTURE].svs.values(),
        capture_to_results[CANCER_CAPTURE].svs.values()
    output:
        svcaller_normal = "{}/svs/igv/{}_svcaller.mut".format(outdir, NORMAL_CAPTURE_STR),
        svcaller_tumor = "{}/svs/igv/{}_svcaller.mut".format(outdir, CANCER_CAPTURE_STR)
    params:
        cancer_str = CANCER_CAPTURE_STR,
        normal_str = NORMAL_CAPTURE_STR,
        svs_dir = outdir + "/svs/svcaller/",
        igvout = outdir + "/svs/igv/"
    shell:
        "generateIGVnavInput_SV.py --input {params.svs_dir} "
                        " --sdid {params.cancer_str} "
                        " --tool svcaller " 
                        " --vcftype somatic " 
                        " --output {params.igvout} && "
        "generateIGVnavInput_SV.py --input {params.svs_dir} "
                        " --sdid {params.normal_str} "
                        " --tool svcaller " 
                        " --vcftype normal "
                        " --output {params.igvout} "


rule svaba_svcalling:
    input:
        normal_bam = capture_to_results[NORMAL_CAPTURE].bamfile,
        tumor_bam = capture_to_results[CANCER_CAPTURE].bamfile,
        reference = reference["bwaIndex"]
    output:
        somatic = "{}/svs/svaba/{}-{}.svaba.somatic.sv.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        germline = "{}/svs/svaba/{}-{}.svaba.germline.sv.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        prefix = "{}/svs/svaba/{}-{}".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        tmpdir = params['scratch']
    threads: params['svaba']['threads']
    log:
        outdir + "/logs/svs/svaba-{}-{}.log".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "svaba run -t {input.tumor_bam} "
        " -n {input.normal_bam} "
        " -G {input.reference} "
        " -p {threads} -a {params.prefix} "
        " && samtools sort -T {params.tmpdir} "
        " {params.prefix}.contigs.bam "
        " -o {params.prefix}.contigs.sort.bam && "
        "samtools index {params.prefix}.contigs.sort.bam"


rule svaba_annotate:
    input:
        somatic = "{}/svs/svaba/{}-{}.svaba.somatic.sv.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        germline = "{}/svs/svaba/{}-{}.svaba.germline.sv.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    output:
        somatic = "{}/svs/svaba/{}-{}.svaba.somatic.annotated.sv.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        germline = "{}/svs/svaba/{}-{}.svaba.germline.annotated.sv.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "annotate_svaba.py {input.somatic} > {output.somatic} && "
        "annotate_svaba.py {input.germline} > {output.germline} "


rule generateIGVnavInput_svaba:
    input:
        somatic = "{}/svs/svaba/{}-{}.svaba.somatic.annotated.sv.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        germline = "{}/svs/svaba/{}-{}.svaba.germline.annotated.sv.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    output:
        somatic = "{}/svs/igv/{}-{}_somatic_svaba.mut".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        germline = "{}/svs/igv/{}-{}_germline_svaba.mut".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        cancer_str = CANCER_CAPTURE_STR,
        normal_str = NORMAL_CAPTURE_STR,
        prefix = outdir + "/svs/igv/{}-{}".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        sdid = "-".join(NORMAL_CAPTURE_STR.split("-")[1:3])
    shell:
        "generateIGVnavInput_SV.py --input {input.somatic} "
                " --sdid {params.sdid}  --tool svaba "
                " --vcftype somatic  --output {params.prefix}  && "        
        "generateIGVnavInput_SV.py --input {input.germline} "
                " --sdid {params.sdid}  --tool svaba " 
                " --vcftype germline --output {params.prefix} "



rule lumpy_svcalling:
    input:
        normal_bam = capture_to_results[NORMAL_CAPTURE].bamfile,
        tumor_bam = capture_to_results[CANCER_CAPTURE].bamfile
    output:
        normal_discordants_bam = "{}/svs/lumpy/{}-discordants.bam".format(outdir, NORMAL_CAPTURE_STR),
        tumor_discordants_bam = "{}/svs/lumpy/{}-discordants.bam".format(outdir, CANCER_CAPTURE_STR),
        normal_splitters_bam = "{}/svs/lumpy/{}-splitters.bam".format(outdir, NORMAL_CAPTURE_STR),
        tumor_splitters_bam = "{}/svs/lumpy/{}-splitters.bam".format(outdir, CANCER_CAPTURE_STR),
        vcf = "{}/svs/lumpy/{}-{}-lumpy.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        tmpdir = os.path.join(params['scratch'], 
                    "lumpy-{}".format(str(uuid.uuid4())))
    threads: params['lumpy']['threads']
    log:
        outdir + "/logs/svs/lumpy-{}-{}.log".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "source activate gatk_3 && "
        "samtools view -@ {threads} -b -F 1294 {input.normal_bam}  "
        "  > {output.normal_discordants_bam}  && "
        "samtools view -@ {threads} -b -F 1294 {input.tumor_bam}  "
        "  > {output.tumor_discordants_bam} &&  "
        "samtools view -@ {threads} -h {input.normal_bam} "
        " | extractSplitReads_BwaMem -i stdin | "
        " samtools view -@ {threads} -Sb -  > {output.normal_splitters_bam} && "
        "samtools view -@ {threads} -h {input.tumor_bam} "
        " | extractSplitReads_BwaMem -i stdin | "
        " samtools view -@ {threads} -Sb -  > {output.tumor_splitters_bam}  && "
        "lumpyexpress -T {params.tmpdir} -B {input.tumor_bam},{input.normal_bam} "
        " -S {output.tumor_splitters_bam},{output.normal_splitters_bam}  " 
        " -D {output.tumor_discordants_bam},{output.normal_discordants_bam} "
        " -o {output.vcf} && "
        "samtools index {output.normal_discordants_bam} && "
        "samtools index {output.normal_splitters_bam} && "
        "samtools index {output.tumor_discordants_bam} && "
        "samtools index {output.tumor_splitters_bam} "


rule generateIGVnavInput_lumpy:
    input:
        vcf = "{}/svs/lumpy/{}-{}-lumpy.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    output:
        "{}/svs/igv/{}-{}_lumpy_len500_SU24.mut".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        "{}/svs/igv/{}-{}_lumpy_len1k_SU50.mut".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        prefix = outdir + "/svs/igv/{}-{}".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        sdid = "-".join(NORMAL_CAPTURE_STR.split("-")[1:3])
    shell:
        "generateIGVnavInput_SV.py --input {input.vcf} "
                " --sdid {params.sdid}  --tool lumpy " 
                " --vcftype somatic --output {params.prefix} "

envvars:
    "GRIDSS_JAR"

# rule gridss_svcalling_normal:
#     input:
#         normal_bam = capture_to_results[NORMAL_CAPTURE].bamfile,
#         reference = reference["bwaIndex"]
#     output:
#         assembly_bam = "{}/svs/gridss/{}-assembly.bam".format(outdir, NORMAL_CAPTURE_STR),
#         vcf = "{}/svs/gridss/{}-gridss.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
#     params:
#         gridss_jar = os.environ.get('GRIDSS_JAR'),
#         jvmheap = '10g',
#         workdir = directory("{}/svs/gridss/".format(outdir))
#     threads: params['gridss']['threads']
#     log:
#         outdir + "/logs/svs/gridss-{}.log".format(NORMAL_CAPTURE_STR)
#     shell:
#         "gridss.sh --reference {input.reference} "
#         " --jvmheap {params.jvmheap} "
#         " --jar {params.gridss_jar} "
#         " --assembly {output.assembly_bam} "
#         " --threads {threads} --steps  ALL "
#         " --workingdir {params.workdir} "
#         " --output {output.vcf} {input.normal_bam} "


# rule gridss_svcalling_tumor:
#     input:
#         tumor_bam = capture_to_results[CANCER_CAPTURE].bamfile,
#         reference = reference["bwaIndex"]
#     output:
#         assembly_bam = "{}/svs/gridss/{}-assembly.bam".format(outdir, CANCER_CAPTURE_STR),
#         vcf = "{}/svs/gridss/{}-gridss.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
#     params:
#         gridss_jar = os.environ.get('GRIDSS_JAR'),
#         jvmheap = '10g',
#         workdir = directory("{}/svs/gridss/".format(outdir))
#     threads: params['gridss']['threads']
#     log:
#         outdir + "/logs/svs/gridss-{}.log".format(CANCER_CAPTURE_STR)
#     shell:
#         "gridss.sh --reference {input.reference} "
#         " --jvmheap {params.jvmheap} "
#         " --jar {params.gridss_jar} "
#         " --assembly {output.assembly_bam} "
#         " --threads {threads} --steps  ALL "
#         " --workingdir {params.workdir} "
#         " --output {output.vcf} {input.tumor_bam}"


# rule gridss_svcalling_somatic:
#     input:
#         normal_bam = capture_to_results[NORMAL_CAPTURE].bamfile,
#         tumor_bam = capture_to_results[CANCER_CAPTURE].bamfile,
#         reference = reference["bwaIndex"]
#     output:
#         assembly_bam = "{}/svs/gridss/{}-{}-assembly.bam".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
#         vcf = "{}/svs/gridss/{}-{}-gridss.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
#     params:
#         gridss_jar = os.environ.get('GRIDSS_JAR'),
#         jvmheap = '10g',
#         workdir = directory("{}/svs/gridss/".format(outdir))
#     threads: params['gridss']['threads']
#     log:
#         outdir + "/logs/svs/gridss-{}-{}.log".format(NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
#     shell:
#         "gridss.sh --reference {input.reference} "
#         " --jvmheap {params.jvmheap} "
#         " --jar {params.gridss_jar} "
#         " --assembly {output.assembly_bam} "
#         " --threads {threads} --steps  ALL "
#         " --workingdir {params.workdir} "
#         " --output {output.vcf} {input.normal_bam} {input.tumor_bam}"


# rule gridss_somatic_filter:
#     input:
#         vcf = "{}/svs/gridss/{}-{}-gridss.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
#     output:
#         vcf = "{}/svs/gridss/{}-{}-gridss.filtered.vcf.bgz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
#     params:
#         pondir = reference["pondir"],
#         script_dir = " -s /opt/gridss -c /opt/gridss ",
#         plotdir = "{}/svs/gridss/".format(outdir),
#         outvcf = "{}/svs/gridss/{}-{}-gridss.filtered.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
#     threads: params["gridss_filter"]["threads"]
#     container: config['container']["gridss"]
#     shell:
#         "gridss_somatic_filter -p {params.pondir} "
#         " -i {input.vcf} "
#         " -o {params.outvcf} "
#         " --plotdir {params.plotdir} {params.script_dir}"


# rule generateIGVnavInput_gridss:
#     input:
#         normal_vcf = "{}/svs/gridss/{}-gridss.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
#         tumor_vcf = "{}/svs/gridss/{}-gridss.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
#     output:
#         normal_mut = "{}/svs/igv/{}_normal_pass_gridss.mut".format(outdir, NORMAL_CAPTURE_STR),
#         somatic_mut = "{}/svs/igv/{}_somatic_pass_gridss.mut".format(outdir, CANCER_CAPTURE_STR)
#     params:
#         nprefix = outdir + "/svs/igv/{}".format(NORMAL_CAPTURE_STR),
#         tprefix = outdir + "/svs/igv/{}".format(CANCER_CAPTURE_STR),
#         sdid = "-".join(NORMAL_CAPTURE_STR.split("-")[1:3])
#     shell:
#         "generateIGVnavInput_SV.py --input {input.normal_vcf} "
#                 " --sdid {params.sdid} --tool gridss " 
#                 " --vcftype normal --output {params.nprefix} && "
#         "generateIGVnavInput_SV.py --input {input.tumor_vcf} "
#                 " --sdid {params.sdid} --tool gridss " 
#                 " --vcftype somatic --output {params.tprefix} "


rule annotate_generateIGVnavInput:
    input:
        svcaller_normal = "{}/svs/igv/{}_svcaller.mut".format(outdir, NORMAL_CAPTURE_STR),
        svcaller_tumor = "{}/svs/igv/{}_svcaller.mut".format(outdir, CANCER_CAPTURE_STR),
        svaba_somatic = "{}/svs/igv/{}-{}_somatic_svaba.mut".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        svaba_germline = "{}/svs/igv/{}-{}_germline_svaba.mut".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        lumpy_len500 = "{}/svs/igv/{}-{}_lumpy_len500_SU24.mut".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        lumpy_1k = "{}/svs/igv/{}-{}_lumpy_len1k_SU50.mut".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        genes = reference["genes_bed"],
        targets = reference['sv_filter'],
        cgcann = reference['cgcann']
        # gridss_somatic = "{}/svs/igv/{}_somatic_pass_gridss.mut".format(outdir, CANCER_CAPTURE_STR),
    output:
        "{}/svs/igv/{}-{}-sv-annotated.txt".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        svs_dir = "{}/svs/igv/".format(outdir),
        capture_kit_id = CANCER_CAPTURE.capture_kit_id
    shell:        
        "generateIGVnavInput_SV.py --input {params.svs_dir} --cgc {input.cgcann} "
            " --annotBed {input.genes} --target {params.capture_kit_id} {input.targets} "
            " --output {output} "
