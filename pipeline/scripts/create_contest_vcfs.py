#!/usr/bin/env python
import click, logging, gzip, re, sys
import pandas as pd
import pybedtools


def extract_allele_counts(vcf_fields):
    allele_freq_info_str = vcf_fields[7].split("|")[0]
    allele_freq_info_tups = [tuple(elem.split("=")) for elem in allele_freq_info_str.split(";")]
    allele_freq_info = dict(filter(lambda tup: len(tup) == 2, allele_freq_info_tups))
    alt_allele_count = float(re.match("[0-9]+", allele_freq_info["AC"]).group())
    total_allele_number = float(re.match("[0-9]+", allele_freq_info["AN"]).group())
    alt_allele_freq = alt_allele_count / total_allele_number
    return (alt_allele_count, total_allele_number, alt_allele_freq)


def write_contest_vcf(vcf_as_df, vcf_header, output_file):
    vcf_idx_reset = vcf_as_df.reset_index(drop=True)
    print >> output_file, vcf_header.strip()
    allele_count_tups = pd.DataFrame(list(vcf_idx_reset.apply(extract_allele_counts, axis = 1)))
    vcf_with_allele_freq = vcf_idx_reset.iloc[:,:7]
    vcf_with_allele_freq["alt_count"] = allele_count_tups[0]
    vcf_with_allele_freq["total_count"] = allele_count_tups[1]
    vcf_with_allele_freq["alt_freq"] = allele_count_tups[2]

    filter = (vcf_with_allele_freq.iloc[:, 9] > 0.1) & (vcf_with_allele_freq.iloc[:, 9] < 0.9)
    vcf_filtered = vcf_with_allele_freq.loc[filter,:]

    for index, row in vcf_filtered.iterrows():
        vcf_fields = list(row)
        last_elem = "AC=%d;AF=%1.3f;AN=%d;CEU={%s*=%1.3f,%s=%1.3f};set=CEU" % \
            (vcf_fields[7], vcf_fields[9], vcf_fields[8],
             vcf_fields[3], vcf_fields[9], vcf_fields[4], 1 - vcf_fields[9])
        print >> output_file, "%s\t%d\t%s\t%s\t%s\t%1.2f\t%s\t%s" % tuple(vcf_fields[:7] + [last_elem])


def extract_vcf_header(vcf_file):
    '''Extract the header from the vcf file. Currently, just assumes that the header is
    within the first 10000 lines of the specified file.'''

    header_str = ""
    line_idx = 0
    curr_line = vcf_file.readline()
    while curr_line != "" and line_idx < 10000:
        if curr_line[0] == "#":
            header_str += curr_line
        curr_line = vcf_file.readline()
        line_idx += 1

    return header_str


def extract_header_unknown_type(vcf_filename):
    file_extension = vcf_filename.split(".")[-1]
    if file_extension == "gz":
        with gzip.open(vcf_filename, 'rb') as vcf_file:
            return extract_vcf_header(vcf_file)
    else:
        if not file_extension == "vcf":
            raise ValueError("Invalid vcf file extension: " + vcf_filename)
        with open(vcf_filename) as vcf_file:
            return extract_vcf_header(vcf_file)


@click.command()
@click.option('--loglevel', default='INFO', help='level of logging')
@click.option('--tmpdir', default='/tmp', help='Temp directory for intermediate files')
@click.option('--max-maf', default=0.9, help='Maximum VCF record MAF')
@click.option('--min-maf', default=0.1, help='Minimum VCF record MAF')
@click.option('--output-filename', default='intersect_snps.vcf', help='Output filename')
@click.argument('genotype-target-bedfile')
@click.argument('eval-target-bedfile')
@click.argument('population-vcf')
def main(genotype_target_bedfile, eval_target_bedfile, population_vcf, output_filename, min_maf, max_maf, tmpdir, loglevel):
    """
Create the contest VCF input file from the designated "genotype" and "eval" target bed region
files and population allele frequency VCF file.
"""

    numeric_level = getattr(logging, loglevel, None)
    logging.basicConfig(level=numeric_level)

    # Get the intersection of the two bed files:
    logging.info("Getting genotype and eval bed intersections...")
    genotype_target_bed = pybedtools.BedTool(genotype_target_bedfile)
    eval_target_bed = pybedtools.BedTool(eval_target_bedfile)
    bed_intersection = genotype_target_bed.intersect(eval_target_bed)
    logging.info("Done.")

    # Get the intersection of the resulting BED file with the population VCF file:
    logging.info("Getting population vcf intersection...")
    pop_vcf = pybedtools.BedTool(population_vcf)
    vcf_bed_intersect = pop_vcf.intersect(bed_intersection)
    logging.info("Done.")

    # Write the VCF records to the output file, the format contest requires:
    logging.info("Writing vcf file output...")
    vcf_header = extract_header_unknown_type(population_vcf)
    vcf_bed_intersect_df = vcf_bed_intersect.to_dataframe()

    # Filter the VCF DF to exclude unwanted SNP types:
    snp_filt = vcf_bed_intersect_df.iloc[:, 3].str.contains("[ACGT]") & \
               vcf_bed_intersect_df.iloc[:, 4].str.contains("[ACGT]")
    vcf_filtered = vcf_bed_intersect_df.loc[snp_filt]

    with open(output_filename, 'w') as output_file:
        write_contest_vcf(vcf_filtered, vcf_header, output_file)
    logging.info("Done.")


if __name__ == "__main__":
    sys.exit(main())
