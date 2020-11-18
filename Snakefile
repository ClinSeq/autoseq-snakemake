import json

include: "rules/utils.smk"
include: "rules/clinseq_barcodes.smk"
configfile: "config.yml"

samples = json.load(open(config['samples']))


# libdir
libdir = 'tests/libraries'

# Need to check clinseq_barcodes and data available for those valid barcodes
sampledata, all_clinseq_barcodes = check_sampledata(libdir, samples)

outdir = 'tests/' + sampledata['sdid']

wildcard_constraints:
    sample = "|".join(all_clinseq_barcodes)


rule all:
    input: 
        expand(outdir + "/fastqs/{sample}_concatenated_{read}.fastq.gz", 
                                        sample=all_clinseq_barcodes,
                                        read=["1", "2"])

include: "rules/alignment.smk"
