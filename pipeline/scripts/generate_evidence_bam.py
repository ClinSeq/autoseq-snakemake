#!/usr/bin/env python

import argparse
import glob
import gzip
import pysam


def parse_mut(mutfile):
    variants = list()
    with open(mutfile, 'r') as fh:
        header = fh.readline()
        for line in fh:
            variants.append(line.strip().split('\t'))

    return variants


def extract_asm(variants):
    """
    Extract assembly ids from vcf file
    """
    assm_ids = set()
    with open(variants, 'r') as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            
            data = line.strip().split('\t')
            info = data[7]
            BEID = info.split(';')[13]
            _assm = BEID.replace('BEID=', '')
            
            for i in _assm.split(','):
                assm_ids.add(i)
    
    return assm_ids


def extract_rp(assm_ids, assm_bam):
    """
    Extract read names from assemblies
    """
    samfile = pysam.AlignmentFile(assm_bam, "rb")
    read_names = list()
    for read in samfile.fetch():
        if read.query_name in assm_ids:
            for _tag in read.tags:
                if _tag[0] == 'ef':                
                    read_names.extend(_tag[1].split())
            
    return read_names


def extract_reads(bam, read_names):
    """
    Extract reads from tumor nodups bam file
    """
    samfile = pysam.AlignmentFile(bam, "rb")
    name_indexed = pysam.IndexedReads(samfile)
    name_indexed.build()

    outfile = pysam.AlignmentFile("output.bam", "wb", template = samfile)
    
    for name in read_names:
        try:
            iterator = name_indexed.find(name)
            for x in iterator:
                outfile.write(x)
        except KeyError:
            print(f"Reads could not find: {name}")
            pass
    return 


def extract_contigs(vcf, align):
    """
    Extract contigs and corresponding read names
    """
    contigs = set()
    with open(vcf, 'r') as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            data = line.strip().split('\t')
            info = ''
            for _info in data[7].split(';'):
                if 'SCTG=c_' in _info:
                    info = _info.replace("SCTG=", '')
            contigs.add(info)
    
    aln_fh = gzip.open(align, 'rb')
    read_names = list()
    for line in aln_fh:
        line = line.strip().decode()
        contig = [i for i in contigs if i in line]
        data = line.split()
        if contig and len(data) == 6:
            _tmp_rn = data[1].split('_')[2]
            read_names.append(_tmp_rn.split('--')[0])
    
    return read_names


def main():
    parser = argparse.ArgumentParser(description=
        'Generate evidence bam for structural variants ')
    parser.add_argument('--bam', required=True, help="Input bam file ")
    parser.add_argument('--svs', required=True, help="structural variants dir as input")
    parser.add_argument('--assembly', help="SV Assembly bam/ alignments file ")
    parser.add_argument('--vcf', help="List of variants as vcf file")
    parser.add_argument('--tool', help="Tool name - Variant callers (gridss, svaba, svcallers)")
    parser.add_argument('--output', help="output tab delimited file for IGVNav, format=output.bam")
    args = parser.parse_args()
    
    # mutfile for lumpy
    if args.svs:
        svaba_aln = glob.glob(args.svs + "/svaba/" + "*.alignments.txt.gz")
        gridss_vcf = glob.glob(args.svs + "/gridss/" + "*gridss.filtered.vcf.bgz")
        gridss_asm = glob.glob(args.svs + "/gridss/" + "*assembly.bam")

    print(svaba_aln)
    print(gridss_vcf)
        
    if args.mut:
        variants = parse_mut(args.mut)
    
    if args.tool == 'gridss':
        assm_ids = extract_asm(args.vcf)
        read_names = extract_rp(assm_ids, args.assembly_bam)
        
    if args.tool == 'svaba':
        read_names = extract_contigs(args.vcf, args.assembly)

    # extract_reads(args.bam, read_names)



if __name__ == "__main__":
    main()

    