import os
import sys
import uuid

capture_name = get_capture_name(CANCER_CAPTURE.capture_kit_id)
somatic_vcf = dict()
normalBam = capture_to_results[NORMAL_CAPTURE].umibam if umi else capture_to_results[NORMAL_CAPTURE].bamfile
cancerBam = capture_to_results[CANCER_CAPTURE].umibam if umi else capture_to_results[CANCER_CAPTURE].bamfile


rule vardict_somatic:
    input:
        reference = reference['reference_genome'],
        reference_dict = reference["reference_dict"],
        normal_bam = normalBam,
        cancer_bam = cancerBam,
        target_bed = reference['targets'][capture_name]['targets-bed-slopped20']
    output:
        "{}/variants/vardict/{}-{}.vardict-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        normalid = compose_sample_str(NORMAL_CAPTURE),
        tumorid = compose_sample_str(CANCER_CAPTURE),
        min_alt_frac = params['vardict']['min_alt_frac'],
        min_num_reads = params['vardict']['min_num_reads'],
        blacklist_bed = reference["targets"][capture_name]["blacklist-bed"],
        pyexe = "/opt/conda/bin/python" if config['container']['base'] != ' ' else sys.executable
    threads: params['vardict']['threads']
    log:
        "{}/logs/variants/{}-{}.vardict-somatic.vcf.gz.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "VarDict -G {input.reference} "
            " -th {threads} "
            "-f {params.min_alt_frac} "
            "-N {params.tumorid} "
            " -b \"{input.cancer_bam}|{input.normal_bam}\" "
            " -c 1 -S 2 -E 3 -g 4 -Q 10  {input.target_bed} "
            " | testsomatic.R " 
            " | var2vcf_paired.pl -P 0.05 -m 4.25 -M -f {params.min_alt_frac} "
            " -N \"{params.tumorid}|{params.normalid}\" "
            " | bcftools filter -e 'STATUS !~ \".*Somatic\"' 2> /dev/null "
            " | {params.pyexe} -c 'from pipeline.utils.bcbio import depth_freq_filter_input_stream; import sys; print(depth_freq_filter_input_stream(sys.stdin, 0, \"bwa\"))' "
            " | sed 's/\\.*Somatic\\\"/Somatic/' " 
            " | sed 's/REJECT,Description=\".*\">/REJECT,Description=\"Not Somatic via VarDict\">/' "
            " | {params.pyexe} -c 'from pipeline.utils.bcbio import call_somatic; import sys; print(call_somatic(sys.stdin.read()))' "
            " | awk -F$'\\t' -v OFS='\\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) }} {{print}}' " 
            " | awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {{next}} {{print}}' " 
            " | vcfstreamsort -w 1000 " 
            " | vt decompose -s - |vt normalize  -r {input.reference} - "
            " | bcftools view --apply-filters .,PASS "
            " | vcfsorter.pl {input.reference_dict} /dev/stdin "
            " | bgzip > {output} 2> {log} && tabix -p vcf {output} "
        

somatic_vcf['vardict'] = "{}/variants/vardict/{}-{}.vardict-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)


rule strelka_somatic:
    input:
        reference = reference['reference_genome'],
        normal_bam = normalBam,
        tumor_bam = cancerBam,
        target_bed = reference['targets'][capture_name]['targets-bed-slopped20-gz']
    output:
        snvs_vcf = "{}/variants/{}-{}-strelka-somatic.passed.snvs.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        indels_vcf = "{}/variants/{}-{}-strelka-somatic.passed.indels.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        rundir = directory("{}/variants/{}-{}-strelka-somatic".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)),
    threads: params['strelka']['threads']
    log:
        "{}/logs/variants/{}-{}-strelka-somatic.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "source activate gatk_3 && "
        "configureStrelkaSomaticWorkflow.py --targeted "
        " --normalBam {input.normal_bam} "
        " --tumorBam {input.tumor_bam} "
        " --ref {input.reference} "
        " --callRegions {input.target_bed} "
        " --runDir {params.rundir} && "
        " {params.rundir}/runWorkflow.py -m local -j {threads} && "
        " zcat {params.rundir}/results/variants/somatic.snvs.vcf.gz | "
        " awk 'BEGIN {{ OFS = \"\t\"}} /^#/ {{ print $0 }} {{if($7==\"PASS\") print $0 }}' "
        " | vt decompose -s - | vt normalize  -r {input.reference} - "
        " | bgzip > {output.snvs_vcf} 2> {log} && "
        " tabix -p vcf {output.snvs_vcf} && "
        " zcat {params.rundir}/results/variants/somatic.indels.vcf.gz | "
        " awk 'BEGIN {{ OFS = \"\t\"}} /^#/ {{ print $0 }} {{if($7==\"PASS\") print $0 }}' "
        " | vt decompose -s - | vt normalize  -r {input.reference} - "
        " | bgzip > {output.indels_vcf} && "
        " tabix -p vcf {output.indels_vcf} && "
        " rm -rf {params.rundir}  "


somatic_vcf['strelka_snvs'] = "{}/variants/{}-{}-strelka-somatic.passed.snvs.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
somatic_vcf['strelka_indels'] = "{}/variants/{}-{}-strelka-somatic.passed.indels.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)


rule gatk4_mutect2:
    input:
        reference = reference['reference_genome'],
        normal_bam = normalBam,
        tumor_bam = cancerBam,
        interval_list = reference['targets'][capture_name]['targets-interval_list-slopped20']
    output:
        bam = "{}/variants/mutect/{}-{}-mutect.bam".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        filtered_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        normalized_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered-normalized.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
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
        " -bamout {output.bam} -O {output.vcf} 2> {log} && "
        " gatk --java-options '-Xmx10g -Djava.io.tmpdir={params.tmpdir}' "
        " FilterMutectCalls  -R {input.reference} "
        " -V {output.vcf}  "
        " -O {output.filtered_vcf} 2>> {log} && "
        " vt decompose -s {output.filtered_vcf} "
        " | vt normalize  -r {input.reference} - "
        " | bgzip > {output.normalized_vcf} 2>> {log} "


somatic_vcf['mutect2'] = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered-normalized.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)


rule varscan_somatic:
    input:
        reference = reference['reference_genome'],
        normal_bam = normalBam,
        tumor_bam = cancerBam,
        target_bed = reference['targets'][capture_name]['targets-bed-slopped20-gz']
    output:
        normal_pileup = "{}/variants/varscan/{}.pileup".format(outdir, NORMAL_CAPTURE_STR),
        tumor_pileup = "{}/variants/varscan/{}.pileup".format(outdir, CANCER_CAPTURE_STR),
        snp_vcf = "{}/variants/varscan/{}-{}-varscan.snp.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        indel_vcf = "{}/variants/varscan/{}-{}-varscan.indel.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        somatic_snp = "{}/variants/varscan/{}-{}-varscan.snp.Somatic.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        somatic_indel = "{}/variants/varscan/{}-{}-varscan.indel.Somatic.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        normalized_snp = "{}/variants/varscan/{}-{}-varscan.snp.Somatic.normalized.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
        normalized_indel = "{}/variants/varscan/{}-{}-varscan.indel.Somatic.normalized.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR),
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
        " {params.extra} --output-vcf 1 2> {log} && "
        " varscan {params.java_options} -Djava.io.tmpdir={params.tmpdir} processSomatic {output.snp_vcf} && "
        " varscan {params.java_options} -Djava.io.tmpdir={params.tmpdir} processSomatic {output.indel_vcf} && "
        " vt decompose -s {output.somatic_snp} | vt normalize  -r {input.reference} -  > {output.normalized_snp} && "
        " vt decompose -s {output.somatic_indel} | vt normalize  -r {input.reference} -  > {output.normalized_indel}  "


somatic_vcf['varscan_snvs'] = "{}/variants/varscan/{}-{}-varscan.snp.Somatic.normalized.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
somatic_vcf['varscan_indels'] = "{}/variants/varscan/{}-{}-varscan.indel.Somatic.normalized.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)


rule somaticseq_merge:
    input:
        **somatic_vcf,
        reference = reference['reference_genome'],
        normal_bam = normalBam,
        tumor_bam = cancerBam
    output:
        rundir = directory("{}/variants/{}-{}-somaticseq".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)),
        consensus_snv = "{}/variants/{}-{}-somaticseq/Consensus.sSNV.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        consensus_indel = "{}/variants/{}-{}-somaticseq/Consensus.sINDEL.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
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
        " source activate gatk_3 && "
        " gatk3 -T CombineVariants "
        " -R {input.reference} --variant {output.consensus_snv} " 
        " --variant {output.consensus_indel} " 
        " --assumeIdenticalSamples  | bgzip > {output.all_somatic} && "
        " source deactivate && "
        " tabix -p vcf {output.all_somatic} 2> {log} "


rule vardict_purecn:
    input:
        reference = reference['reference_genome'],
        reference_dict = reference["reference_dict"],
        dbsnp = reference["dbSNP"],
        normal_bam = normalBam,
        cancer_bam = cancerBam,
        target_bed = reference['targets'][capture_name]['targets-bed-slopped20']
    output:
        "{}/variants/{}-{}.vardict-somatic-purecn.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        tmpdir = params['scratch'],
        normalid = compose_sample_str(NORMAL_CAPTURE),
        tumorid = compose_sample_str(CANCER_CAPTURE),
        min_alt_frac = params['vardict']['min_alt_frac'],
        min_num_reads = 6,
        tmpvcf = f"{params['scratch']}/{uuid.uuid4()}.vcf.gz"
    threads: params['vardict']['threads']
    shell:
        "VarDict -G {input.reference} "
            " -th {threads} "
            "-f {params.min_alt_frac} "
            "-N {params.tumorid} "
            "-r {params.min_num_reads} "
            " -b \"{input.cancer_bam}|{input.normal_bam}\" "
            " -c 1 -S 2 -E 3 -g 4 -Q 10 {input.target_bed} "
            " | testsomatic.R " 
            " | var2vcf_paired.pl -P 0.9 -m 4.25 -f {params.min_alt_frac} "
            " -N \"{params.tumorid}|{params.normalid}\" "
            " | awk -F$'\\t' -v OFS='\\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) }} {{print}}' " 
            " |  awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {{next}} {{print}}'" 
            " | sed 's/Somatic;/Somatic;SOMATIC;/g' " 
            " | sed '/^#CHROM/i ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic event\">' " 
            " | vcfstreamsort -w 1000 " 
            " | bcftools view --apply-filters .,PASS " 
            " | vcfsorter.pl {input.reference_dict} /dev/stdin "
            " | bgzip > {params.tmpvcf} && tabix -p vcf {params.tmpvcf} && "
            " bcftools annotate --annotation {input.dbsnp} --columns ID "
            " --output-type z --output {output} {params.tmpvcf} "
            " && tabix -p vcf {output} && "
            "rm {params.tmpvcf} "


rule vcf_add_sample:
    input:
        germline_vcf = "{}/variants/{}-all.germline.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
        tumor_bam = cancerBam
    output:
        vcf = "{}/variants/{}-and-{}.germline-variants-with-somatic-afs.vcf.gz".format(
            outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        tumorid = CANCER_SAMPLE_STR,
        tmpdir = os.path.join(params['scratch'], 
                    "vcfaddsample-{}".format(str(uuid.uuid4())))
    threads: params['vcfaddsample']['threads']
    shell:
        "vcf_filter.py --no-filtered  {input.germline_vcf} "
        " sq --site-quality 5 | bgzip > {params.tmpdir}.vcf.gz && "
        " vcf_add_sample.py --filter_hom --samplename {params.tumorid}  " 
        " {params.tmpdir}.vcf.gz  {input.tumor_bam} "
        " | bgzip > {output.vcf} && "
        " tabix -p vcf {output.vcf} && " 
        " rm {params.tmpdir}.vcf.gz "


rule make_allelic_fraction_track:
    input:
        vcf = "{}/variants/{}-and-{}.germline-variants-with-somatic-afs.vcf.gz".format(
            outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    output:
        "{}/variants/{}-and-{}.germline-variants-somatic-afs.bedGraph".format(
            outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    shell:
        "generate_allelic_fraction_bedGraph.py  --output {output}  {input.vcf} "


rule somatic_generateIGVnav:
    input:
        somatic = "{}/variants/{}-{}-all.somatic.vep.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        vardict_vcf = "{}/variants/vardict/{}-{}.vardict-somatic.vep.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        oncokb = reference['oncokb'],
        cgcann = reference['cgcann']
    output:
        "{}/{}-{}-igvnav-input.txt".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        vcftype = "somatic"
    shell:
        "generateIGVnavInput.py {input.somatic} {input.oncokb} "
        " {params.vcftype} --cgc {input.cgcann} --output {output} "
