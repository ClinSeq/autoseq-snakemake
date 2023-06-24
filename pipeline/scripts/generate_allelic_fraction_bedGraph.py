#!/usr/bin/env python

import sys
import vcf
import click
import logging

class InvalidInputVcfException(Exception):
    pass


@click.command()
@click.option('--loglevel', default='INFO', help='level of logging')
@click.option('--output', default="output.bedGraph", help="output file")
@click.argument('vcf-filename')
def extract_somatic_afs(vcf_filename, output, loglevel):
    logging.basicConfig(level=logging.DEBUG)

    vcf_reader = vcf.Reader(filename=vcf_filename)

    # Extract the names of all tumor (i.e. non-normal) samples:
    non_normal_samples = [sample for sample in vcf_reader.samples if sample.split("-")[3] != "N"]

    # Currently just cope with input files containing exactly one tumor sample:
    if len(non_normal_samples) != 1:
        raise InvalidInputVcfException("Input VCF must include exactly one tumor (i.e. non-normal) sample.")

    non_normal_sample = non_normal_samples[0]

    with open(output, 'w') as output_file:
        for record in vcf_reader:
            tumor_genotype = record.genotype(non_normal_sample)
            if tumor_genotype.data.DP is not None:
                if tumor_genotype.data.DP > 0:
                    allelic_fraction = float(tumor_genotype.data.AD[1])/float(tumor_genotype.data.DP)
                    print("{}\t{}\t{}\t{}".format(record.CHROM, record.POS, record.POS, round(allelic_fraction, 2)), file=output_file)


if __name__ == "__main__":
    sys.exit(extract_somatic_afs())
