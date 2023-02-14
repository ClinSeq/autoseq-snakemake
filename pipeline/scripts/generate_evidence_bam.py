#!/usr/bin/env python

import os
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


def extract_reads(bam, read_names, output):
    """
    Extract reads from tumor nodups bam file
    """
    samfile = pysam.AlignmentFile(bam, "rb")
    name_indexed = pysam.IndexedReads(samfile)
    name_indexed.build()

    outdir = os.path.dirname(output)
    bam_prefix = os.path.basename(output).split('.bam')[0]
    bam_sorted = str(os.path.join(outdir, bam_prefix + ".sorted.bam"))
    
    outfile = pysam.AlignmentFile(output, "wb", template = samfile)
    
    for name in read_names:
        try:
            iterator = name_indexed.find(name)
            for x in iterator:
                outfile.write(x)
        except KeyError:
            print(f"Reads could not find: {name}")
            pass
    
    outfile.close()
    pysam.sort(output, '-o', bam_sorted)
    pysam.index(bam_sorted)
    
    os.remove(output)

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


def extract_readnames(bam):
    """
    extract readnames from svcaller evidance bam files
    """
    samfile = pysam.AlignmentFile(bam, "rb")
    read_names = list()
    for read in samfile.fetch():
        read_names.append(read.query_name)
    
    return read_names


def main():
    parser = argparse.ArgumentParser(description=
        'Generate evidence bam for structural variants ')
    parser.add_argument('--bam', required=True, help="Input bam file ")
    parser.add_argument('--svs', required=True, help="structural variants dir as input")
    parser.add_argument('--output', help="output tab delimited file for IGVNav, format=output.bam")
    args = parser.parse_args()
    
    # extract files from svs dir
    if args.svs:
        svaba_vcf = glob.glob(args.svs + "/svaba/" + "*.svaba.somatic.annotated.sv.vcf")
        svaba_aln = glob.glob(args.svs + "/svaba/" + "*.alignments.txt.gz")
        gridss_vcf = glob.glob(args.svs + "/gridss/" + "*gridss.filtered.vcf")
        gridss_asm = glob.glob(args.svs + "/gridss/" + "*assembly.bam")
        svcaller_bams = glob.glob(args.svs + "/*-CFDNA-*.bam")
    
    read_names = set()

    print("Processing gridss output ...")
    if len(gridss_vcf) == 0:
        print("No gridss output !!")
    else:
        assm_ids = extract_asm(gridss_vcf[0])
        read_names.update(set(extract_rp(assm_ids, gridss_asm[0])))
    print("GRIDSS - Done")
        
    print("Processing svaba output ...")
    if len(svaba_aln) == 0:
        print("No svaba output !!")
    else:
        read_names.update(set(extract_contigs(svaba_vcf[0], svaba_aln[0])))
    print("SvABA - Done")

    print("Processing svcaller output ...")
    for event in svcaller_bams:
        read_names.update(set(extract_readnames(event)))
    print("SVCALLER - Done")

    extract_reads(args.bam, read_names, args.output)



if __name__ == "__main__":
    main()

    