#!/usr/bin/env python

import os
import argparse
import shutil
import logging

import pysam

def extract_reads(samfile, name_indexed, read_names, output):
    """
    Extract reads from nodups bam file
    """

    outdir = os.path.dirname(output)
    bam_prefix = os.path.basename(output).split('.bam')[0]
    bam_sorted = str(os.path.join(outdir, bam_prefix + ".sorted.bam"))
    
    outfile = pysam.AlignmentFile(output, "w", template = samfile)
    
    for name in set(read_names):
        try:
            iterator = name_indexed.find(name)
            rmdups = set()
            for x in iterator:
                if x.is_supplementary:
                    continue

                mpos = ":".join(map(str, [x.reference_name, 
                                          x.reference_start, 
                                          x.reference_end]))
                if mpos not in rmdups:
                    # print(x, mpos)  
                    outfile.write(x)
                    rmdups.update(mpos)
        except KeyError:
            print(f"Reads could not find: {name}")
            pass
    
    outfile.close()
    pysam.sort(output, '-o', bam_sorted)
    shutil.move(bam_sorted, output)
    pysam.index(output)

    return 


def extract_readnames(vcf):
    """
    extract readnames from gridss vcf file
    """
    vcffile = pysam.VariantFile(vcf)
    read_names = list()
    for record in vcffile.fetch():
        read_names.extend(record.info['BPNAMES'])

    return set(read_names)


def apply_filters(vcf, name_indexed):
    """
    Apply SAME_START_READS filter to gridss variants

    Check if all supporting reads have same start position
    If so, Its most likely not a true variant.
    """
    vcffile = pysam.VariantFile(vcf)
    vcffile.header.add_meta('FILTER', 
                            items=[('ID', 'SAME_START_READS'), 
                                   ('Description', 
                                    'All supporting reads are having same start positions')])
    
    basedir = os.path.dirname(vcf)
    tmppath = basedir + "temp_filters.vcf"
    tmpvcf = pysam.VariantFile(tmppath, "w", header=vcffile.header)

    read_names = list()
    for record in vcffile.fetch():
        try:
            rn = set(record.info['BPNAMES'])
        except KeyError:
            logging.error("Need to add BPNAMES to INFO column in vcf file")
            exit 

        read_names.extend(rn)
        mpos = set()
        for name in rn:
            try:
                iterator = name_indexed.find(name)
                mpos.update(set([":".join(map(str, [x.reference_name, x.reference_start, x.reference_end]))
                                  for x in iterator if not x.is_supplementary]))
            except KeyError:
                print(f"Reads could not find: {name}")
                pass

        if len(mpos) == 1:
            record.filter.add("SAME_START_READS")
        
        tmpvcf.write(record)
    
    shutil.move(tmppath, vcf)

    return read_names


def setup_logging(loglevel="INFO"):
    """
    Set up logging
    :param loglevel: loglevel to use, one of ERROR, WARNING, DEBUG, INFO (default INFO)
    :return:
    """
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level,
            format='%(levelname)s %(asctime)s %(funcName)s - %(message)s')
    logging.info("Started log with loglevel %(loglevel)s" % {"loglevel": loglevel})


def main():
    parser = argparse.ArgumentParser(description=
        'Generate evidence bam for structural variants ')
    parser.add_argument('--bam', required=True, help="Input bam file ")
    parser.add_argument('--vcf', required=True, help="structural variants vcf as input")
    parser.add_argument('--filter-vcf', action='store_true', help="Apply SAME_START_READS filter on input vcf file")
    parser.add_argument('--output', help="evidence bam for gridss variants")
    args = parser.parse_args()
    
    setup_logging()

    logging.info(f"Loading {args.bam} ..")
    samfile = pysam.AlignmentFile(args.bam, "rb")
    logging.info(f"Indexing by read names {args.bam} ")
    name_indexed = pysam.IndexedReads(samfile)
    name_indexed.build()

    if args.filter_vcf:
        logging.info(f"Applying filters on {args.vcf} and extracting read names")
        read_names = apply_filters(args.vcf, name_indexed)
    else:
        logging.info(f"Extracting variants supporting reads - {args.vcf} ")
        read_names = extract_readnames(args.vcf)

    logging.info(f"Generating evidence Bam {args.output}")
    extract_reads(samfile, name_indexed, read_names, args.output)
    logging.info(f"Done - Generating evidence Bam {args.output}")


if __name__ == "__main__":
    main()

    