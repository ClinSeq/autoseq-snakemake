import sys

capture_name = get_capture_name(CANCER_CAPTURE.capture_kit_id)
somatic_vcf = dict()

rule vardict_somatic:
    input:
        reference = reference['reference_genome'],
        reference_dict = reference["reference_dict"],
        normal_bam = capture_to_results[NORMAL_CAPTURE].umibam,
        cancer_bam = capture_to_results[CANCER_CAPTURE].umibam,
        target_bed = reference['targets'][capture_name]['targets-bed-slopped20-gz']
    output:
        "{}/variants/vardict/{}-{}.vardict-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        normalid = compose_sample_str(NORMAL_CAPTURE),
        tumorid = compose_sample_str(CANCER_CAPTURE),
        min_alt_frac = params['vardict']['min_alt_frac'],
        min_num_reads = params['vardict']['min_num_reads'],
        target_bed = reference["targets"][capture_name]["blacklist-bed"]
    threads: params['vardict']['threads']
    log:
        "{}/logs/variants/{}-{}.vardict-somatic.vcf.gz.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    run:
        freq_filter = (" bcftools filter -e 'STATUS !~ \".*Somatic\"' 2> /dev/null "
                       "| %s -c 'from pipeline.utils.bcbio import depth_freq_filter_input_stream; import sys; print depth_freq_filter_input_stream(sys.stdin, %s, \"%s\")' " %
                       (sys.executable, 0, 'bwa'))

        somatic_filter = (" sed 's/\\.*Somatic\\\"/Somatic/' "  # changes \".*Somatic\" to Somatic
                          "| sed 's/REJECT,Description=\".*\">/REJECT,Description=\"Not Somatic via VarDict\">/' "
                          "| %s -c 'from pipeline.utils.bcbio import call_somatic; import sys; print call_somatic(sys.stdin.read())' " % sys.executable)
        
        blacklist_filter = ""

        if params.blacklist_bed:
            blacklist_filter = " | intersectBed -a . -b {} | ".format(params.blacklist_bed)
        
        if params.min_num_reads:
            min_num_reads = "-r {} ".format(params.min_num_reads)

        cmd = "vardict-java -G {} ".format(input.reference) + \
              "-f {} ".format(params.min_alt_frac) + \
              "-N {} ".format(params.tumorid) + \
              " {} ".format(min_num_reads) + \
              " -b \"{}|{}\" ".format(input.cancer_bam, input.normal_bam) + \
              " -c 1 -S 2 -E 3 -g 4 -Q 10 " + " {} ".format(input.target_bed) + \
              " | testsomatic.R " + \
              " | var2vcf_paired.pl -P 0.05 -m 4.25 -M " + "-f {} ".format(params.min_alt_frac) + \
              " -N \"{}|{}\" ".format(params.tumorid, params.normalid) + \
              " | " + freq_filter + " | " + somatic_filter + " | " + \
              " awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}' " + \
              " | awk -F$'\t' -v OFS='\t' '$1!~/^#/ && $4 == $5 \{next\} \{print\}' " + \
              " | vcfstreamsort -w 1000 " + \
              " | vt decompose -s - |vt normalize  -r {} - ".format(input.reference) + \
              " | bcftools view --apply-filters .,PASS " + \
              " | vcfsorter.pl {} /dev/stdin ".format(input.reference_dict) + \
              " {} ".format(blacklist_filter) + \
              " | bgzip > {output} && tabix -p vcf {output}".format(output=output)
        
        shell(cmd)


somatic_vcf['vardict'] = "{}/variants/vardict/{}-{}.vardict-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)


rule strelka_somatic:
    input:
        reference = reference['reference_genome'],
        normal_bam = capture_to_results[NORMAL_CAPTURE].umibam,
        tumor_bam = capture_to_results[CANCER_CAPTURE].umibam,
        target_bed = reference['targets'][capture_name]['targets-bed-slopped20-gz']
    output:
        rundir = directory("{}/variants/{}-{}-strelka-somatic".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)),
        snvs_vcf= "{}/variants/{}-{}-strelka-somatic/results/variants/somatic.passed.snvs.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        indels_vcf= "{}/variants/{}-{}-strelka-somatic/results/variants/somatic.passed.indels.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: params['strelka']['threads']
    log:
        "{}/logs/variants/{}-{}-strelka-somatic.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "configureStrelkaSomaticWorkflow.py --targeted "
        " --normalBam {input.normal_bam} "
        " --tumorBam {input.tumor_bam} "
        " --ref {input.reference} "
        " --callRegions {input.target_bed} "
        " --runDir {output.rundir} && "
        " {output.rundir}/runWorkflow.py -m local -j 20 && "
        " zcat {output.rundir}/results/variants/somatic.snvs.vcf.gz | "
        " awk 'BEGIN {{ OFS = \"\t\"}} /^#/ {{ print $0 }} {{if($7==\"PASS\") print $0 }}' "
        " | bgzip > {output.snvs_vcf} && "
        " tabix -p vcf {output.snvs_vcf} && "
        " zcat {output.rundir}/results/variants/somatic.indels.vcf.gz | "
        " awk 'BEGIN {{ OFS = \"\t\"}} /^#/ {{ print $0 }} {{if($7==\"PASS\") print $0 }}' "
        " | bgzip > {output.indels_vcf} && "
        " tabix -p vcf {output.indels_vcf}"


somatic_vcf['strelka_snvs'] = "{}/variants/{}-{}-strelka-somatic/results/variants/somatic.passed.snvs.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
somatic_vcf['strelka_indels'] = "{}/variants/{}-{}-strelka-somatic/results/variants/somatic.passed.indels.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)


rule gatk4_mutect2:
    input:
        reference = reference['reference_genome'],
        normal_bam = capture_to_results[NORMAL_CAPTURE].umibam,
        tumor_bam = capture_to_results[CANCER_CAPTURE].umibam,
        interval_list = reference['targets'][capture_name]['targets-interval_list-slopped20']
    output:
        bam = "{}/variants/mutect/{}-{}-mutect.bam".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        filtered_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
    params:
        java_options = params['gatk4']['mutect2']['java_options'],
        normalid = compose_sample_str(NORMAL_CAPTURE),
        tumorid = compose_sample_str(CANCER_CAPTURE),
        tmpdir = params['scratch']
    threads: params['gatk4']['threads']
    log:
        "{}/logs/variants/{}-{}-gatk4-mutect-somatic.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        "gatk --java-options '{params.java_options} -Djava.io.tmpdir={params.tmpdir}' "
        " Mutect2  -R {input.reference} "
        " -I {input.tumor_bam} "
        " -I {input.normal_bam} "
        " -tumor {params.tumorid}  " 
        " -normal {params.normalid} " 
        " -L {input.interval_list} " 
        " --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter "
        " -bamout {output.bam} -O {output.vcf} && "
        " gatk --java-options '-Xmx10g -Djava.io.tmpdir={params.tmpdir}' "
        " FilterMutectCalls  -R {input.reference} "
        " -V {output.vcf}  "
        " -O {output.filtered_vcf} "


somatic_vcf['mutect2'] = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)


rule varscan_somatic:
    input:
        reference = reference['reference_genome'],
        normal_bam = capture_to_results[NORMAL_CAPTURE].umibam,
        tumor_bam = capture_to_results[CANCER_CAPTURE].umibam,
        target_bed = reference['targets'][capture_name]['targets-bed-slopped20-gz']
    output:
        normal_pileup = "{}/variants/varscan/{}.pileup".format(outdir, NORMAL_CAPTURE_STR),
        tumor_pileup = "{}/variants/varscan/{}.pileup".format(outdir, CANCER_CAPTURE_STR),
        snp_vcf = "{}/variants/varscan/{}-{}-varscan.snp.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        indel_vcf = "{}/variants/varscan/{}-{}-varscan.indel.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        somatic_snp = "{}/variants/varscan/{}-{}-varscan.snp.Somatic.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        somatic_indel = "{}/variants/varscan/{}-{}-varscan.indel.Somatic.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
    params:
        java_options = params['varscan']['java_options'],
        extra = params['varscan']['extra'],
        tmpdir = params['scratch']
    threads: params['varscan']['threads']
    log:
        "{}/logs/variants/{}-{}-varscan-somatic.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)        
    shell:
        "samtools mpileup -C50 -f {input.reference} -l {input.target_bed} "
        " {input.normal_bam} > {output.normal_pileup} && "
        " samtools mpileup -C50 -f {input.reference} -l {input.target_bed} "
        " {input.tumor_bam} > {output.tumor_pileup} && "
        " varscan {params.java_options} -Djava.io.tmpdir={params.tmpdir} somatic "
        " {output.normal_pileup} {output.tumor_pileup} "
        " --output-snp {output.snp_vcf} "
        " --output-indel {output.indel_vcf} "
        " {params.extra} --output-vcf 1 && "
        " varscan {params.java_options} -Djava.io.tmpdir={params.tmpdir} processSomatic {output.snp_vcf} && "
        " varscan {params.java_options} -Djava.io.tmpdir={params.tmpdir} processSomatic {output.indel_vcf} "


somatic_vcf['varscan_snvs'] = "{}/variants/varscan/{}-{}-varscan.snp.Somatic.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
somatic_vcf['varscan_indels'] = "{}/variants/varscan/{}-{}-varscan.indel.Somatic.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)


rule somaticseq_merge:
    input:
        **somatic_vcf,
        reference = reference['reference_genome'],
        normal_bam = capture_to_results[NORMAL_CAPTURE].umibam,
        tumor_bam = capture_to_results[CANCER_CAPTURE].umibam
    output:
        rundir = directory("{}/variants/{}-{}-somaticseq".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)),
        consensus_snv = "{}/variants/{}-{}-somatic-seq/Consensus.sSNV.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        consensus_indel = "{}/variants/{}-{}-somatic-seq/Consensus.sINDEL.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        all_somatic = "{}/variants/{}-{}-all.somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        tmpdir = params['scratch']
    threads: params['somaticseq']['threads']
    log:
        "{}/logs/variants/{}-{}-somaticseq.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "source activate somaticseqenv && "
        "run_somaticseq.py  --output-directory {output.rundir} "
        " --genome-reference {input.reference} paired "
        " --tumor-bam-file {input.tumor_bam} " 
        " --normal-bam-file {input.normal_bam} " 
        " --mutect2-vcf {input.mutect2} " 
        " --varscan-snv {input.varscan_snvs} "
        " --varscan-indel {input.varscan_indels} "
        " --vardict-vcf {input.vardict} " 
        " --strelka-snv {input.strelka_snvs} "
        " --strelka-indel {input.strelka_indels} && "
        " source deactivate && source activate gatk_3 && "
        " gatk3 --java-options '-Djava.io.tmpdir={params.tmpdir}' -T CombineVariants "
        " -R {input.reference} --variant {output.consensus_snv} " 
        " --variant {output.consensus_indel} " 
        " --assumeIdenticalSamples  | bgzip > {output.all_somatic} && "
        " source deactivate && "
        " && tabix -p vcf {output.all_somatic} "
    