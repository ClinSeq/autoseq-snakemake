import os
import uuid 


rule cnvkit:
    input:
        bam = outdir + "/bams/{sample}_nodups.bam",
        reference = lambda wildcards: get_cnvkitref(wildcards, reference)
    output:
        cns = outdir + "/cnv/{sample}.cns",
        cnr = outdir + "/cnv/{sample}.cnr"
    params:
        prefix = os.path.basename(outdir + "/bams/{sample}_nodups"),
        tmpdir = os.path.join(params['scratch'], 
                    "cnvkit-{}".format(str(uuid.uuid4())))
    threads: params['cnvkit']['threads']
    shell:
        "mkdir -p {params.tmpdir} && "
        "cnvkit.py batch {input.bam}  -r {input.reference} "
        " -d {params.tmpdir} "
        " && cp {params.tmpdir}/{params.prefix}.cns {output.cns}  "
        " && cp {params.tmpdir}/{params.prefix}.cnr {output.cnr}  "
        " && rm -r {params.tmpdir}"


# need to discuss
# rule cnvkit_fix:
#     input:
#         cns = outdir + "/cnv/{sample}.cns",
#         cnr = outdir + "/cnv/{sample}.cnr",
#         ref = lambda wildcards: get_cnvkitref(wildcards, reference)
#     output:
#         cns = outdir + "/cnv/{sample}-fixed.cns",
#         cnr = outdir + "/cnv/{sample}-fixed.cnr",
#     threads: params['cnvkit-fix']['threads']
#     shell:
#         "fix_cnvkit.py --input-cnr {input.cnr} "
#         " --input-cns {input.cns} "
#         " --input-reference {input.ref} "
#         " --output-cnr {output.cnr} "
#         " --output-cns {output.cns} "


rule cnvkit_cnstoseg:
    input:
        cns = outdir + "/cnv/{sample}.cns",
    output:
        seg = outdir + "/cnv/{sample}.seg",
    threads: params['cnstoseg']['threads']
    shell:
        "cnvkit.py export seg  -o {output.seg}  {input.cns}"


rule cnvkit_tracks:
    input:
        cns = outdir + "/cnv/{sample}.cns",
        cnr = outdir + "/cnv/{sample}.cnr"        
    output:
        profile_bedgraph = outdir + "/cnv/{sample}_profile.bedGraph",
        segments_bedgraph = outdir + "/cnv/{sample}_segments.bedGraph",
    threads: params['cnvkit_tracks']['threads']
    shell:
        "awk '$1 != \"chromosome\" {{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5}}' "
        " {input.cnr} > {output.profile_bedgraph} "
        " && awk '$1 != \"chromosome\" {{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5}}' "
        " {input.cns} > {output.segments_bedgraph} "


