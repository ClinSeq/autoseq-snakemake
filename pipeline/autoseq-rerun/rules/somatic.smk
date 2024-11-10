
somatic_vcf = dict()
cancerBam = capture_to_results[CANCER_CAPTURE].umibam if umi else capture_to_results[CANCER_CAPTURE].bamfile

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
        "normalId=`samtools view -H {input.normal_bam} | grep \"^@RG\" | tr \'\\t\' \'\\n\' | grep \'^SM\' | cut -d\':\' -f2 ` && "
        "tumorId=`samtools view -H {input.tumor_bam} | grep \"^@RG\" |  tr \'\\t\' \'\\n\' | grep \'^SM\' | cut -d\':\' -f2 ` && "
        "gatk --java-options '{params.java_options} -Djava.io.tmpdir={params.tmpdir}' "
        " Mutect2  -R {input.reference} "
        " -I {input.tumor_bam} "
        " -I {input.normal_bam} "
        " -tumor $tumorId  " 
        " -normal $normalId " 
        " -L {input.interval_list} " 
        " --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter "
        " -bamout {output.bam} -O {output.vcf} 2> {log} && "
        " gatk --java-options '-Xmx10g -Djava.io.tmpdir={params.tmpdir}' "
        " FilterMutectCalls  -R {input.reference} "
        " --max-alt-allele-count 2 "
        " -V {output.vcf}  "
        " -O {output.filtered_vcf} 2>> {log} && "
        " vt decompose -s {output.filtered_vcf} "
        " | vt normalize  -r {input.reference} - "
        " | bgzip > {output.normalized_vcf} 2>> {log} "


rule sage_somatic:
    input:
        normal_bam = normalBam,
        tumor_bam = cancerBam,
        reference = reference['reference_genome'],
        panel_bed = reference['targets'][capture_name]['targets-bed-slopped20'],
        known_hotspots = reference['wgs']['hartwig']['known-hotspots-somatic-vcf'],
        high_confi_bed = reference['wgs']['hartwig']['NA12878-highconf-bed'],
        ensembl_dir = reference['wgs']['hartwig']['ensembl-dir']
    output:
        "{}/variants/{}-{}-hartwig-sage-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        normalid = compose_sample_str(NORMAL_CAPTURE),
        tumorid = compose_sample_str(CANCER_CAPTURE),
        jarfile = os.environ.get('SAGE_JAR'),
        workdir = outdir + "/variants",
        minaf = 0.0002,
        hpmintq = 150,
        min_mapq = 20,
        min_baseq = 30,
        min_paneltq = 250
    threads: params['sage']['threads']
    container: containers['gridss']
    log:
        "{}/logs/variants/{}-{}-sage-somatic.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        " source activate gridss-env && "
        " bedtools merge -i {input.panel_bed} > {params.workdir}/targets_nonoverlap.bed 2> {log} && "
        " java -Xms4G -Xmx32G -cp {params.jarfile} "
        " com.hartwig.hmftools.sage.SageApplication -threads 16 "
        " -reference {params.normalid} -reference_bam {input.normal_bam}"
        " -tumor {params.tumorid} -tumor_bam {input.tumor_bam} "
        " -ref_genome_version 37  -ref_genome {input.reference} "
        " -hotspots {input.known_hotspots} "
        " -panel_bed {params.workdir}/targets_nonoverlap.bed "
        " -ensembl_data_dir {input.ensembl_dir} "
        " -hard_min_tumor_vaf {params.minaf} -hotspot_min_tumor_qual {params.hpmintq} "
        " -hotspot_min_tumor_vaf {params.minaf} -min_map_quality {params.min_mapq} "
        " -min_avg_base_qual {params.min_baseq} -panel_min_tumor_qual {params.min_paneltq} "
        " -panel_min_tumor_vaf {params.minaf} -panel_only -write_bqr_data "
        " -out {output} 2>> {log} && rm {params.workdir}/targets_nonoverlap.bed"


rule bcftools_filter:
    input:
        sage_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        mutect_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered-normalized.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    output:
        sage_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        mutect_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.pass.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    threads: 4
    params:
        tmp_sage_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        tmp_mutect_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.pass.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    log:
        "{}/logs/variants/{}-{}-sage-somatic-filter.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "bcftools filter -e 'FILTER!=\"PASS\"' {input.sage_vcf} > {params.tmp_sage_vcf} && "
        "bcftools filter -e 'FILTER!=\"PASS\"' {input.mutect_vcf} > {params.tmp_mutect_vcf} && "
        "bgzip {params.tmp_sage_vcf} && tabix -p vcf {output.sage_vcf} && "
        "bgzip {params.tmp_mutect_vcf} && tabix -p vcf {output.mutect_vcf} "


rule bcftools_concat:
    input:
        sage_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        mutect_vcf = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.pass.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    output:
        "{}/variants/{}-{}-all.somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    threads: 8
    params:
        ordered_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.reordered.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
    log:
        "{}/logs/variants/{}-{}-bcftools-concat.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    shell:
        "bcftools query -l {input.sage_vcf} | sort > sample_names.txt && "
        "bcftools view -Oz -S sample_names.txt {input.sage_vcf} -o {params.ordered_vcf} && "
        "tabix -p vcf {params.ordered_vcf} &&  "
        "bcftools concat -a -D {input.mutect_vcf} {params.ordered_vcf} "
        " | bgzip > {output}  2> {log} && " 
        "tabix -p vcf {output} && rm sample_names.txt "


# rule sage_splitvcf:
#     input:
#         "{}/variants/{}-{}-hartwig-sage-somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     output:
#         snv = "{}/variants/{}-{}-hartwig-sage-somatic.snvs.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
#         indel =  "{}/variants/{}-{}-hartwig-sage-somatic.indels.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     threads: 1
#     params:
#         tmp_vcf = "{}/variants/{}-{}-hartwig-sage-somatic.pass.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR), 
#         snv = "{}/variants/{}-{}-hartwig-sage-somatic.snvs.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
#         indel =  "{}/variants/{}-{}-hartwig-sage-somatic.indels.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     container: containers['somaticseq']
#     log:
#         "{}/logs/variants/{}-{}-sage-somatic-splitvcf.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     shell:
#         "source activate somaticseqenv && "
#         "bcftools filter -e 'FILTER!=\"PASS\"' {input} > {params.tmp_vcf} && "
#         "splitVcf.py -infile {params.tmp_vcf} -snv {params.snv} -indel {params.indel} && "
#         " bgzip {params.snv} && bgzip {params.indel} "


# somatic_vcf['mutect2'] = "{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered-normalized.vcf.gz".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
# somatic_vcf['sage_snv'] = "{}/variants/{}-{}-hartwig-sage-somatic.snvs.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
# somatic_vcf['sage_indel'] = "{}/variants/{}-{}-hartwig-sage-somatic.indels.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)


# rule somaticseq_merge:
#     input:
#         **somatic_vcf,
#         reference = reference['reference_genome'],
#         normal_bam = normalBam,
#         tumor_bam = cancerBam
#     output:
#         rundir = directory("{}/variants/{}-{}-somaticseq".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)),
#         consensus_snv = "{}/variants/{}-{}-somaticseq/Consensus.sSNV.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
#         consensus_indel = "{}/variants/{}-{}-somaticseq/Consensus.sINDEL.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     params:
#         tmpdir = params['scratch']
#     threads: params['somaticseq']['threads']
#     container: containers['somaticseq']
#     log:
#         "{}/logs/variants/{}-{}-somaticseq.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     shell:
#         "source activate somaticseqenv && "
#         "run_somaticseq.py  --output-directory {output.rundir} "
#         " --genome-reference {input.reference} paired "
#         " --tumor-bam-file {input.tumor_bam} " 
#         " --normal-bam-file {input.normal_bam} " 
#         " --mutect2-vcf {input.mutect2} "
#         " --arbitrary-snvs {input.sage_snv} "
#         " --arbitrary-indels {input.sage_indel} 2> {log} "


# rule gatk3_combinevariants:
#     input:
#         reference = reference['reference_genome'],
#         consensus_snv = "{}/variants/{}-{}-somaticseq/Consensus.sSNV.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
#         consensus_indel = "{}/variants/{}-{}-somaticseq/Consensus.sINDEL.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     output:
#         "{}/variants/{}-{}-all.somatic.vcf.gz".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     threads: params['somaticseq']['threads']
#     container: containers['gatk3']
#     log:
#         "{}/logs/variants/{}-{}-combine_somaticvcf.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
#     shell:    
#         " source activate gatk_3 && "
#         " gatk3 -T CombineVariants "
#         " -R {input.reference} --variant {input.consensus_snv} " 
#         " --variant {input.consensus_indel} " 
#         " --assumeIdenticalSamples  | bgzip > {output} && "
#         " tabix -p vcf {output} 2> {log} "
    

rule somatic_generateIGVnav:
    input:
        somatic = "{}/variants/{}-{}-all.somatic.vep.vcf".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR),
        oncokb = reference['oncokb'],
        cgcann = reference['cgcann']
    output:
        "{}/{}-{}-igvnav-input.txt".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)
    params:
        vcftype = "somatic"
    shell:
        "generateIGVnavInput.py {input.somatic} {input.oncokb} "
        " {params.vcftype} --cgc {input.cgcann} --output {output} "


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
        germline_vcf = "{}/variants/haplotypecaller/{}.haplotypecaller-germline.vcf.gz".format(outdir, NORMAL_CAPTURE_STR),
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


rule gatk4_haplotypecaller_tumor:
    input:
        reference = reference["reference_genome"],
        tumor_bam = cancerBam,
        germline_vcf = "{}/variants/{}-all.germline.vep.vcf".format(outdir, NORMAL_CAPTURE_STR)
    output:
        vcf = "{}/variants/{}-haplotypecaller.vcf".format(outdir, CANCER_CAPTURE_STR)
    params:
        tmpdir = params['scratch']
    threads: params['gatk4']['threads']
    log:
        "{}/logs/variants/{}-{}-gatk4-haplotypecaller-tumor.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        "gatk HaplotypeCaller -R {input.reference} "
        " -I {input.tumor_bam} "
        " -O {output.vcf} "
        " -ERC GVCF "
        " -L {input.germline_vcf} 2> {log} "


rule gatk4_haplotypecaller_genotype_tumor:
    input:
        reference = reference["reference_genome"],
        vcf = "{}/variants/{}-haplotypecaller.vcf".format(outdir, CANCER_CAPTURE_STR),
        germline_vcf = "{}/variants/{}-all.germline.vep.vcf".format(outdir, NORMAL_CAPTURE_STR)
    output:
        vcf = "{}/variants/{}-haplotypecaller.genotyped.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    params:
        vcf = "{}/variants/{}-haplotypecaller.genotyped.vcf".format(outdir, CANCER_CAPTURE_STR)
    threads: params['gatk4']['threads']
    log: 
        "{}/logs/variants/{}-{}-haplotypecaller-tumor-genotyped.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR) 
    shell:
        "gatk GenotypeGVCFs -R {input.reference} "
        "  -V {input.vcf}  "
        "  -O {params.vcf} "
        "  -L {input.germline_vcf} 2> {log} && "
        " bgzip {params.vcf} && "
        " tabix -p vcf {output.vcf} 2>> {log} "
    

rule bcftools_merge:
    input:
        germline_vcf = "{}/variants/{}-all.germline.vep.vcf".format(outdir, NORMAL_CAPTURE_STR),
        genotyped_tvcf = "{}/variants/{}-haplotypecaller.genotyped.vcf.gz".format(outdir, CANCER_CAPTURE_STR)
    output:
        vcf = "{}/variants/{}-{}.germline_variants_with_taf.vcf".format(outdir, NORMAL_CAPTURE_STR, CANCER_CAPTURE_STR)
    params:
        gvcf = "{}/variants/{}-all.germline.vep.vcf.gz".format(outdir, NORMAL_CAPTURE_STR)
    threads: params['gatk4']['threads']
    log:
       "{}/logs/variants/{}-{}-bcftools-merge.log".format(outdir, CANCER_CAPTURE_STR, NORMAL_CAPTURE_STR)  
    shell:
        "bgzip -c {input.germline_vcf} > {params.gvcf} && "
        "tabix -p vcf {params.gvcf} && "
        "bcftools merge {params.gvcf} {input.genotyped_tvcf} "
        " -O v -o {output.vcf} 2> {log} && "
        " rm {params.gvcf} "
