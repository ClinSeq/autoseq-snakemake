import json

include: "rules/utils.smk"
include: "rules/clinseq_barcodes.smk"
configfile: "config.yml"

samples = json.load(open(config['samples']))
params = config['params']
refjson = json.load(open(config['reference']))

# loading reference files
basepath = os.path.dirname(config['reference'])
reference = make_paths_absolute(refjson, basepath)

print(basepath)
# libdir
libdir = 'tests/libraries'

# Need to check clinseq_barcodes and data available for those valid barcodes
sampledata, all_clinseq_barcodes = check_sampledata(libdir, samples)

outdir = 'tests/' + sampledata['sdid']

wildcard_constraints:
    sample = "|".join(all_clinseq_barcodes)


rule all:
    input: 
        expand(outdir + "/bams/{sample}-nodups.bam", sample=all_clinseq_barcodes),
        expand(outdir + "/bams/{sample}_mapped-1.bam", sample=all_clinseq_barcodes)
                                    

include: "rules/alignment.smk"
include: "rules/umi_processing.smk"
