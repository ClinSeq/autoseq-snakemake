

somatic_vcf = dict()
normalBam = capture_to_results[NORMAL_CAPTURE].bamfile
cancerBam = capture_to_results[CANCER_CAPTURE].bamfile

targets_dir = reference["wgs"]["targets"]["bed"]
suffix = list()
for fn in os.listdir(targets_dir):
    suffix.append(fn.split('.')[1])

vardict_vcf_prefix = "{}/variants/vardict/{}-{}.vardict-somatic".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
vardict_log_prefix = "{}/logs/variants/{}-{}-vardict-somatic".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 

rule vardict_somatic:
    input:
        reference = reference['reference_genome'],
        reference_dict = reference["reference_dict"],
        normal_bam = normalBam,
        cancer_bam = cancerBam,
        target_bed = targets_dir + "human_g1k_v37_decoy.{suf}.bed"
    output:
        vardict_vcf_prefix + "{suf}.vcf.gz"
    wildcard_constraints:
        suf = "|".join(suffix)
    params:
        normalid = compose_sample_str(NORMAL_CAPTURE),
        tumorid = compose_sample_str(CANCER_CAPTURE),
        min_alt_frac = params['vardict']['min_alt_frac'],
        min_num_reads = params['vardict']['min_num_reads']
    threads: params['vardict']['threads']
    log:
        vardict_log_prefix + "{suf}.vcf.gz.log"
    shell:
        "VarDict -G {input.reference} "
            " -th {threads} --fisher "
            "-f {params.min_alt_frac} "
            "-N {params.tumorid} "
            " -b \"{input.cancer_bam}|{input.normal_bam}\" "
            " -c 1 -S 2 -E 3 {input.target_bed} "
            " | var2vcf_paired.pl -P 0.05 -M -f {params.min_alt_frac} "
            " -N \"{params.tumorid}|{params.normalid}\" "
            " | bgzip > {output} 2> {log} && tabix -p vcf {output} "
        

rule vardict_vcfmerge:
    input:
        expand(vardict_vcf_prefix + "{suf}.vcf.gz", suf=suffix)
    output:
        vardict_vcf_prefix + ".vcf.gz"
    threads: params['bcftools']['threads']
    log:
        "{}/logs/variants/{}-{}-vardict-somatic-merge.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        "bcftools concat --threads {threads} -a "
        " -O z {input} > {output} 2> {log} "


rule vardict_filter:
    input:
        vcf = vardict_vcf_prefix + ".vcf.gz",
        reference_dict = reference["reference_dict"],
        reference = reference['reference_genome']
    output:
        vardict_vcf_prefix + "-filtered.vcf.gz"
    params:
        pyexe = "/opt/conda/bin/python" if config['container']['base'] != ' ' else sys.executable
    threads: params['bcftools']['threads']
    log:
        "{}/logs/variants/{}-{}-vardict-somatic-filter.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        " bcftools filter -e 'STATUS !~ \".*Somatic\"' {input.vcf} 2> /dev/null "
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


rule vardict_purecn:
    input:
        vcf = vardict_vcf_prefix + ".vcf.gz",
        reference_dict = reference["reference_dict"],
        reference = reference['reference_genome'],
        dbsnp = reference["dbSNP"]
    output:
        vardict_vcf_prefix + "-purecn.vcf.gz"
    params:
        tmpvcf = f"{params['scratch']}/{uuid.uuid4()}.vcf.gz"
    threads: params['bcftools']['threads']
    log:
        "{}/logs/variants/{}-{}-vardict-somatic-purecn.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        " zcat {input.vcf} "
        " | awk -F$'\\t' -v OFS='\\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) }} {{print}}' " 
        " | awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {{next}} {{print}}'" 
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


somatic_vcf['vardict'] = "{}/variants/vardict/{}-{}.vardict-somatic-filtered.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)

intervals_dir = reference["wgs"]["targets"]["interval_list"]
mutect_vcf_prefix = "{}/variants/mutect/{}-{}-gatk-mutect-somatic".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
mutext_log_prefix = "{}/logs/variants/{}-{}-gatk4-mutect-somatic".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 

rule gatk4_mutect2:
    input:
        reference = reference['reference_genome'],
        normal_bam = normalBam,
        tumor_bam = cancerBam,
        interval_list = intervals_dir + "human_g1k_v37_decoy.{suf}.interval_list"
    output:
        bam = mutect_vcf_prefix + "{suf}.bam",
        vcf = mutect_vcf_prefix + "{suf}.vcf.gz",
        filtered_vcf = mutect_vcf_prefix + "{suf}-filtered.vcf.gz"
    wildcard_constraints:
        suf = "|".join(suffix)
    params:
        java_options = params['gatk4']['mutect2']['java_options'],
        normalid = compose_sample_str(NORMAL_CAPTURE),
        tumorid = compose_sample_str(CANCER_CAPTURE),
        tmpdir = params['scratch']
    threads: params['gatk4']['threads']
    log:
        mutext_log_prefix + "{suf}.log"        
    shell:
        "gatk --java-options '{params.java_options} -Djava.io.tmpdir={params.tmpdir}' "
        " Mutect2  -R {input.reference} "
        " -I {input.tumor_bam} "
        " -I {input.normal_bam} "
        " -tumor {params.tumorid}  " 
        " -normal {params.normalid} " 
        " -L {input.interval_list} "
        " -bamout {output.bam} -O {output.vcf} 2> {log} && "
        " gatk --java-options '-Xmx10g -Djava.io.tmpdir={params.tmpdir}' "
        " FilterMutectCalls  -R {input.reference} "
        " -V {output.vcf}  "
        " -O {output.filtered_vcf} 2>> {log} "
      

rule mutect2_vcfmerge:
    input:
        expand(mutect_vcf_prefix + "{suf}-filtered.vcf.gz", suf=suffix)
    output:
        mutect_vcf_prefix + "-filtered.vcf.gz"     
    threads: params['bcftools']['threads']
    log:
        "{}/logs/variants/{}-{}-mutect-somatic-merge.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        "bcftools concat --threads {threads} -a "
        " -O z {input} > {output} 2> {log} && "
        " tabix -p vcf {output} "


rule mutect2_normalize:
    input:
        vcf = mutect_vcf_prefix + "-filtered.vcf.gz",
        reference = reference['reference_genome']
    output:
        normalized_vcf = mutect_vcf_prefix + "-filtered-normalized.vcf.gz"
    params:
        tmpdir = params['scratch']
    threads: params['bcftools']['threads']
    log:
        "{}/logs/variants/{}-{}-mutect-somatic-merge.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        " vt decompose -s {input.vcf} "
        " | vt normalize  -r {input.reference} - "
        " | bgzip > {output.normalized_vcf} 2> {log} "


somatic_vcf['mutect2'] = mutect_vcf_prefix + "-filtered-normalized.vcf.gz"

rule sage_somatic:
    input:
        normal_bam = normalBam,
        tumor_bam = cancerBam,
        reference = reference['reference_genome'],
        panel_bed = reference['wgs']['hartwig']['actionable-somatic-panel-bed'],
        known_hotspots = reference['wgs']['hartwig']['known-hotspots-somatic-vcf'],
        high_confi_bed = reference['wgs']['hartwig']['NA12878-highconf-bed'],
        ensembl_dir = reference['wgs']['hartwig']['ensembl-dir']
    output:
        "{}/variants/{}-{}-hartwig-sage-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        normalid = compose_sample_str(NORMAL_CAPTURE),
        tumorid = compose_sample_str(CANCER_CAPTURE),
        jarfile = os.environ.get('SAGE_JAR')
    threads: params['sage']['threads']
    container: containers['gridss']
    log:
        "{}/logs/variants/{}-{}-sage-somatic.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        " source activate gridss-env && "
        " java -Xms4G -Xmx32G -cp {params.jarfile} "
        " com.hartwig.hmftools.sage.SageApplication -threads 16 "
        " -reference {params.normalid} -reference_bam {input.normal_bam}"
        " -tumor {params.tumorid} -tumor_bam {input.tumor_bam} "
        " -ref_genome_version 37  -ref_genome {input.reference} "
        " -hotspots {input.known_hotspots} "
        " -panel_bed {input.panel_bed} "
        " -high_confidence_bed {input.high_confi_bed} "
        " -ensembl_data_dir {input.ensembl_dir} "
        " -out {output} 2> {log} "


rule sage_splitvcf:
    input:
        "{}/variants/{}-{}-hartwig-sage-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        snv = "{}/variants/{}-{}-hartwig-sage-somatic.snvs.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        indel =  "{}/variants/{}-{}-hartwig-sage-somatic.indels.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: 1
    container: containers['somaticseq']
    log:
        "{}/logs/variants/{}-{}-sage-somatic-splitvcf.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "source activate somaticseqenv && "
        "splitVcf.py -infile {input} -snv {output.snv} -indel {output.indel}"


somatic_vcf['sage_snv'] = "{}/variants/{}-{}-hartwig-sage-somatic.snvs.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
somatic_vcf['sage_indel'] = "{}/variants/{}-{}-hartwig-sage-somatic.indels.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)


rule somaticseq_merge:
    input:
        **somatic_vcf,
        reference = reference['reference_genome'],
        normal_bam = normalBam,
        tumor_bam = cancerBam
    output:
        rundir = directory("{}/variants/{}-{}-somaticseq".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)),
        consensus_snv = "{}/variants/{}-{}-somaticseq/Consensus.sSNV.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        consensus_indel = "{}/variants/{}-{}-somaticseq/Consensus.sINDEL.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        tmpdir = params['scratch']
    threads: 16
    container: containers['somaticseq']
    log:
        "{}/logs/variants/{}-{}-somaticseq.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "source activate somaticseqenv && "
        "somaticseq_parallel.py  --output-directory {output.rundir} "
        " --genome-reference {input.reference} "
        " --threads {threads} paired "
        " --tumor-bam-file {input.tumor_bam} " 
        " --normal-bam-file {input.normal_bam} " 
        " --mutect2-vcf {input.mutect2} " 
        " --vardict-vcf {input.vardict} " 
        " --arbitrary-snvs {input.sage_snv} "
        " --arbitrary-indels {input.sage_indel} "


rule gatk3_combinevariants:
    input:
        reference = reference['reference_genome'],
        consensus_snv = "{}/variants/{}-{}-somaticseq/Consensus.sSNV.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        consensus_indel = "{}/variants/{}-{}-somaticseq/Consensus.sINDEL.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    output:
        "{}/variants/{}-{}-all.somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: 8
    container: containers['gatk3']
    log:
        "{}/logs/variants/{}-{}-somaticseq-vcfmerge.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        " source activate gatk_3 && "
        " gatk3 -T CombineVariants "
        " -R {input.reference} --variant {output.consensus_snv} " 
        " --variant {output.consensus_indel} " 
        " --assumeIdenticalSamples  | bgzip > {output} 2> {log} && "
        " tabix -p vcf {output} 2>> {log} "


rule somatic_generateIGVnav:
    input:
        somatic = "{}/variants/{}-{}-all.somatic.vep.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        oncokb = reference['oncokb'],
        cgcann = reference['cgcann']
    output:
        "{}/{}-{}-igvnav-input.txt".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        vcftype = "somatic"
    log:
        "{}/logs/{}-{}.somatic_generate_igvnav.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        "generateIGVnavInput.py {input.somatic} {input.oncokb} {params.vcftype} --wgs --cgc {input.cgcann} --output {output} 2> {log} "


haplotype_jc_vcf_prefix = "{}/variants/haplotypecaller/{}-{}.haplotypecaller-joint-calling".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
haplotype_jc_log_prefix = "{}/logs/variants/haplotypecaller/{}-{}.haplotypecaller-joint-calling".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)


rule gatk4_haplotypecaller_jc:
    input:
        normal_bam = normalBam,
        tumor_bam = cancerBam,
        reference = reference['reference_genome'],
        dbsnp = reference["dbSNP"],
        interval_list = intervals_dir + "human_g1k_v37_decoy.{suf}.interval_list"
    output:
        vcf = haplotype_jc_vcf_prefix + ".{suf}.vcf.gz",
    wildcard_constraints:
        suf = "|".join(suffix)
    params:
        java_options = params["gatk4"]["haplotypecaller"]["java_options"]
    threads: 4
    log:
        haplotype_jc_log_prefix + ".{suf}.log"
    shell:
        "gatk --java-options '{params.java_options}' "
        " HaplotypeCaller   "
        " -R {input.reference}  "
        " -I {input.normal_bam}  "
        " -I {input.tumor_bam} "
        " -L {input.interval_list} "
        " --dbsnp {input.dbsnp} "
        " -O {output.vcf} 2> {log} "


rule picard_vcfmerge:
    input:
        expand(haplotype_jc_vcf_prefix + ".{suf}.vcf.gz", suf = suffix)
    output:
        vcf = haplotype_jc_vcf_prefix + ".vcf.gz"
    params:
        java_options = params["gatk4"]["haplotypecaller"]["java_options"]
    threads: params["gatk4"]["threads"]
    log:
        haplotype_jc_log_prefix + ".log"
    run:
        input_bams = ["I="+i for i in input]
        ibams = " ".join(input_bams)
        shell("picard MergeVcfs {ibams} O={output.vcf} 2> {log}")


rule make_allelic_fraction_track:
    input:
        vcf = haplotype_jc_vcf_prefix + ".vcf.gz"
    output:
        "{}/variants/{}-and-{}.germline-variants-somatic-afs.bedGraph".format(
            outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    log:
       "{}/logs/{}-{}.make_allelic_fraction_track.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        "generate_allelic_fraction_bedGraph.py  --output {output}  {input.vcf} 2> {log} "

