#!/usr/bin/env python
# 
# Script to Generate IGVNav input file from combined vcf file with Oncogenicity Annoation
# Written for Liqbio pipeline on 5 March 2019
# Required Input files - vep annoated vcf file, OncoKB allvariants text file, somatic (or) germline analysis info

###################################################################
#modified on 26-6-19 : creating symlinks for IGVnav related files##
###################################################################


import argparse
import vcf
import os
import re
import pandas as pd


def csq_parsing(csq, csq_keys):
    # parsing CSQ taq from VeP annotation 
    # return canonical transcript csq taq as dict with corresponding keys

    csq_dicts = []
    can_trans = {}

    for transcript in csq:
        tmp = {csq_keys[idx]: ann for idx,ann in enumerate(transcript.split('|')) }
        csq_dicts.append(tmp)

    for trans in csq_dicts:
        if trans['CANONICAL'] == 'YES':
            can_trans = trans
    
    if not can_trans:
        can_trans = csq_dicts[0]
    
    return can_trans


def loadOncoKB(filepath):
    # Load the OncoKB database and converting it into json
    OncoKB = {}

    with open(filepath,'r', encoding='utf-8', errors='ignore') as oncokb:
        #header = oncokb.readline()
        for line in oncokb.readlines():   
            data = line.rstrip().split('\t')
            gene = data[3]
            protein_change = data[5]
            oncogenicity = data[6]
            mutation_effect = data[7]
            pmids = data[8] if len(data) > 8 else '' 
            if gene in OncoKB:
                OncoKB[gene].update({protein_change:[oncogenicity, mutation_effect, pmids]})
            else:
                OncoKB.update({gene:{protein_change:[oncogenicity, mutation_effect, pmids]}})
    return OncoKB


def loadCGC(filepath):
    """
    load Cancer Gene Census file
    """ 
    cgc_ann = dict() 
    with open(filepath, 'r') as fh:
        header = fh.readline()
        for line in fh.readlines():
            data = line.strip('\n').split('\t')
            ensemblID = data[1]
            gene_symbol = data[0]
            ann = data[4]
            cgc_ann[ensemblID] = (ann, gene_symbol)
    
    return cgc_ann


def create_symlink(travers_dir_name, src_dir, igvnav_dirname_dst, suffix):
    "Recursively Traverse through the directory and create symlink"
    for root, dirs, files in os.walk(travers_dir_name):
        for each_file in files:
            if each_file.endswith(suffix) and not os.path.exists(os.path.join(igvnav_dirname_dst,each_file)):
                os.symlink(os.path.join(root,each_file), os.path.join(igvnav_dirname_dst,each_file))
    return 


def load_hotspot(hs_fpath, vtype):
    """
    Load cancer hotspot files and convert it into pyranges
    for downstream processes.
    """
    hs_lookup = pd.read_csv(hs_fpath)

    if vtype == "snv":
        hs_lookup['start_aa'] = hs_lookup['amino_acid_position'] - 1
        hs_lookup['end_aa'] = hs_lookup['amino_acid_position'] + 1
    
    return hs_lookup


def tri_to_single_aa(hgvsp):
    """
    """
    aa_lookup = {
        'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
        'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
    }
    res = ""
    # split the code based on digit and special character
    tricode = re.split('([0-9]+|[a-zA-Z \s\n\.]+)', hgvsp)
    tricode = [i for i in tricode if i]

    for c in tricode:
        tmp = ''
        # check if any special character occurs 
        if(not any(not _.isalnum() for _ in c)):
            # frameshift
            if 'fs' in c:
                tmp = c.split('fs')
                c = tmp[0]
            if not c.isdigit() and c.upper() in aa_lookup:
                aa = aa_lookup[c.upper()]
                c = aa
        res += c
        if tmp:
            res += 'fs' + tmp[1]

    return res


def hotspot_aapos(gene, one_aa, aa_range, hs_lookup):
    """
    """
    ann = ''
    aa_pos = int(aa_range[0])
    res = hs_lookup.query('gene == @gene and hgvsp == @one_aa')
    if len(res) > 0:
        ann = "Hotspot"
    else:
        res = hs_lookup.query('gene == @gene and \
                               (start_aa <= @aa_pos <= end_aa)')
        ann = "Warmspot" if len(res) > 0 else ''
    
    return ann


def hotspot_aarange(gene, one_aa, aa_range, hs_lookup):
    """
    """
    ann = ''
    start, end = int(aa_range[0]), int(aa_range[1])
    res = hs_lookup.query('gene == @gene and \
                            (start_aa <= @start <= end_aa) and \
                            (start_aa <= @end <= end_aa) ')
    if len(res) > 0:
        ann = "Hotspot"
    else:
        res = hs_lookup.query('gene == @gene and \
                            ((start_aa <= @start <= end_aa) or \
                            (start_aa <= @end <= end_aa)) ')
        ann = "Warmspot" if len(res) > 0 else ''
    
    return ann


def hotspot_splice(gene, one_aa, aa_range, gpos, hs_lookup):
    """
    """
    ann = ''
    aa_pos = int(aa_range[0])
    res = hs_lookup.query('gene == @gene and amino_acid_position == @aa_pos')
    if len(res) > 0:
        genomic_pos = res['genomic_position'].iloc[0]
        ann = "Hotspot" if gpos in genomic_pos else ''
    
    return ann


def check_noncoding_hotspot(chr, pos):
    """
    """
    gpos = ":".join(map(str, [chr, pos]))
    lookup_loci = ['5:1295228', '5:1295250']
    lookup_interval = {"7" : (116411709, 116411902)}
    
    if gpos in lookup_loci:
        return True

    if chr in lookup_interval:
        start, end = lookup_interval[chr]
        if start <= pos <= end:
            return True

    return False


def annotate_hotspot(query, vartype, conseq, hotspot_lookup, gpos = ''):
    """
    Annotate each variant with cancer hotspot information
    """
    gene, hgvsp = query

    ann = ''
    one_aa = tri_to_single_aa(hgvsp.split('p.')[1])
    aa_range = re.findall(r'\d+', one_aa)

    if vartype == "snv":
        ann = hotspot_aapos(gene, one_aa, aa_range, hotspot_lookup)

    if 'inframe' in conseq:
        if '_' in one_aa and len(aa_range) > 1:
            ann = hotspot_aarange(gene, one_aa, aa_range, hotspot_lookup)
        else:
            ann = hotspot_aapos(gene, one_aa, aa_range, hotspot_lookup)
    
    if 'splice' in conseq and vartype == "snv":
        ann = hotspot_splice(gene, one_aa, aa_range, gpos, hotspot_lookup)

    return ann


if __name__ == "__main__": 
    # Parsing Commandline Arguments using argparse
    #

    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', help="Input VCF file for annotation")
    parser.add_argument('oncokb', help="OncoKB - all variants tab demilited file")
    parser.add_argument('vcftype', help="somatic (or) germline vcf")
    parser.add_argument('--wgs', action='store_true', default=False, help="tag to use WGS filter")
    parser.add_argument('--cgc', default=False, help="Cancer Gene Census Annotation")
    parser.add_argument('--hotspot-snv', help="Cancer Hotspot SNV list, as csv file")
    parser.add_argument('--hotspot-indel', help="Cancer Hotspot Indel list, as csv file")
    parser.add_argument('-v', '--vardict', help="Adding vardict long indels into IGVNav")
    parser.add_argument('--output', help="output tab demilited file for IGVNav", default='output.txt')
    args = parser.parse_args()

    #OncoKB_lookup = loadOncoKB("/home/chimera/genome-files/allAnnotatedVariants.txt")
    OncoKB_lookup = loadOncoKB(args.oncokb)
    wgs = args.wgs

    #vcf_reader = vcf.Reader(open("/home/chimera/Downloads/new_vcf_format.vcf", 'r'))
    if args.vcf.endswith('.gz'):
        vcf_reader = vcf.Reader(filename=args.vcf)
    else:
        vcf_reader = vcf.Reader(open(args.vcf, 'r'))

    vcftype = args.vcftype
    sdid = "-".join(os.path.basename(args.vcf).split('-')[1:3])
    sid = "-".join(os.path.basename(args.vcf).split('-')[0:5])
    output_file = open(args.output, 'w') 
    variants = list()
    cancer_genes = dict()

    if args.vardict:
        vardict_vcf = vcf.Reader(open(args.vardict, 'r'))

    if args.cgc:
        cgc_ann = loadCGC(args.cgc)


    igvnav_dirname_dst = os.path.join(os.path.dirname(os.path.abspath(args.output)), 'IGVnav')
    src_dir = os.path.abspath(os.path.dirname(os.path.abspath(args.output)))

    try:
        if not os.path.exists(igvnav_dirname_dst): os.mkdir(igvnav_dirname_dst)
        for each_input in [('variants','.vep.vcf')]:
            dir_name = os.path.join(src_dir,each_input[0])
            create_symlink(dir_name, src_dir, igvnav_dirname_dst, each_input[1])
    except Exception as e:
        print(e)

    #######################################################################################################################

    # output file headers 

    if vcftype == "somatic":
        output_file.write('\t'.join(['CHROM','START','END','REF','ALT', 'CALL', 'TAG', 'NOTES', 'GENE',  'ENSEMBLID', 'IMPACT', 'CONSEQUENCE', 'TRANSCRIPT', 'HGVSc', 'HGVSp', 'HOTSPOT', 'T_DP', 'T_ALT', 'T_VAF', 'N_DP', 'N_ALT', 'N_VAF', 'CLIN_SIG', 'RSID', 'gnomAD', 'BRCAEx', 'OncoKB', 'CGC_ANN', 'NUM_TOOLS']) + "\n")
    elif vcftype == "germline":
        output_file.write('\t'.join(['CHROM','START','END','REF','ALT', 'CALL', 'TAG', 'NOTES', 'GENE', 'ENSEMBLID', 'IMPACT', 'CONSEQUENCE', 'TRANSCRIPT', 'HGVSc','HGVSp',  'HOTSPOT', 'N_DP', 'N_ALT', 'N_VAF', 'T_VAF', 'CLIN_SIG', 'RSID', 'gnomAD', 'BRCAEx', 'OncoKB', 'CGC_ANN']) + "\n")


    if "CSQ" in vcf_reader.infos:
        csq_description = vcf_reader.infos["CSQ"].desc
        # Extract header keys (typically found in the description as 'Allele|Consequence|...|SIFT|PolyPhen')
        csq_headers = csq_description.split(": ")[1].split("|")

    
    hotspot_snv = load_hotspot(args.hotspot_snv, "snv") if args.hotspot_snv else None
    hotspot_indel = load_hotspot(args.hotspot_indel, "indel") if args.hotspot_indel else None

    for record in vcf_reader:
        try:
            # skip germline variants which doesn't have CSQ annotations
            if "CSQ" not in record.INFO:
                continue 

            canonical_trans = csq_parsing(record.INFO['CSQ'], csq_headers)
            gene = canonical_trans['SYMBOL']
            ensembl_id = canonical_trans['Gene']
            aa = canonical_trans['Amino_acids'].split('/')
            protein_position = canonical_trans['Protein_position'].split('/')
            clinsig =  canonical_trans['CLIN_SIG']
            impact = canonical_trans['IMPACT']
            brcaEx = canonical_trans['BrcaEx_ClinicalSignificance']
            hgvsp = canonical_trans['HGVSp']
        except Exception as e:
            raise e
        
        if 'gnomAD_AF' in canonical_trans:
            gnomAD = canonical_trans['gnomAD_AF']
        elif not wgs:
            gnomAD = canonical_trans['gnomADe_AF']
        else:
            gnomAD = canonical_trans['gnomADg_AF']

        oncogenicity = ''
        filter_col = ''
        is_CGC = False

        consequence = canonical_trans['Consequence']
        is_splice_variant = True if 'splice_region_variant' in consequence else False
        
        # Filter variants 
        if not record.FILTER: 
            filter_col = "PASS"
        else: 
            filter_col = record.FILTER[0]

        # Oncogenicity annotation from OncoKB
        if gene in OncoKB_lookup:
            if len(aa) > 1 and len(protein_position) > 0:
                protein_change = aa[0] +  protein_position[0] + aa[1]
                if protein_change in OncoKB_lookup[gene]:
                    oncogenicity = OncoKB_lookup[gene][protein_change]

        cgcann = ''
        if ensembl_id in cgc_ann:
            cgcann = cgc_ann[ensembl_id][0]
        
        ann_hotspot = 'Hotspot' if check_noncoding_hotspot(record.CHROM, record.POS) else ''
        if hgvsp and not (hotspot_snv.empty and hotspot_indel.empty):
            if record.is_snp:
                gpos = ":".join(map(str, [record.CHROM, record.POS]))
                ann_hotspot = annotate_hotspot((gene, hgvsp), "snv", consequence, hotspot_snv, gpos)
            else:
                ann_hotspot = annotate_hotspot((gene, hgvsp), "indel", consequence, hotspot_indel)

        # processing somatic vcf file
        if vcftype == "somatic":
            try:
                n_sample = [sample.sample for sample in record.samples if "-N-" in sample.sample ][0]
                t_sample = [sample.sample for sample in record.samples if "-N-" not in sample.sample ][0]
            except Exception as e:
                n_sample = "NORMAL"
                t_sample = "TUMOR"
            
            normal = record.genotype(n_sample)
            tumor = record.genotype(t_sample)
            
            ## for backward compatability
            normal_dp = normal['DP'] if hasattr(normal.data, 'DP') else sum(normal['DP4'])
            normal_alt = normal['AD'][1] if hasattr(normal.data, 'AD') else sum(normal['DP4'][2:])
            normal_vaf = normal['AF'] if hasattr(normal.data, 'AF') else normal['VAF']

            tumor_dp = tumor['DP'] if hasattr(tumor.data, 'DP') else sum(tumor['DP4'])
            tumor_alt = tumor['AD'][1] if hasattr(tumor.data, 'AD') else sum(tumor['DP4'][2:])
            tumor_vaf = tumor['AF'] if hasattr(tumor.data, 'AF') else tumor['VAF']

            num_tools = int(record.INFO['NUM_TOOLS']) if "NUM_TOOLS" in record.INFO else int(1)
            rsid = canonical_trans['Existing_variation']

            if (filter_col == 'PASS' or filter_col == 'LowQual') and \
                (impact == 'HIGH' or impact == 'MODERATE' or is_splice_variant) and not wgs:
                # forming variant string to remove duplicates
                # eg: 3-113275658-G-TTTTTTT
                tmp_str = "-".join(map(str, [record.CHROM, record.POS, record.REF, record.ALT[0]]))
                variants.append(tmp_str)
                output_file.write('\t'.join(map(str, [record.CHROM, record.POS-1, record.POS,
                                                    record.REF, record.ALT, '', '', '', gene, ensembl_id,
                                                    impact, canonical_trans['Consequence'], canonical_trans['Feature'],
                                                    canonical_trans['HGVSc'], canonical_trans['HGVSp'], ann_hotspot, tumor_dp, tumor_alt, 
                                                    tumor_vaf, normal_dp, normal_alt, normal_vaf, 
                                                    clinsig, rsid, gnomAD, brcaEx, oncogenicity, cgcann, num_tools])) + "\n")

            # filter for WGS samples
            if wgs and tumor_alt >= 5 and (impact == 'HIGH' or impact == 'MODERATE'):
                
                output_file.write('\t'.join(map(str, [record.CHROM, record.POS-1, record.POS,
                                                    record.REF, record.ALT, '', '', '', gene, ensembl_id,
                                                    impact, canonical_trans['Consequence'], canonical_trans['Feature'],
                                                    canonical_trans['HGVSc'], canonical_trans['HGVSp'], ann_hotspot, tumor_dp, tumor_alt, 
                                                    tumor_vaf, normal_dp, normal_alt, normal_vaf, 
                                                    clinsig, rsid, gnomAD, brcaEx, oncogenicity, cgcann, num_tools])) + "\n")
        
        elif vcftype == "germline":
            normal = record.samples[0]
            tumor_vaf = 0
            
            if len(record.samples) > 1 and len(record.ALT) == 1:
                tumor = record.samples[1]
                if tumor['AD'] and tumor['DP']:
                    tumor_vaf = float(tumor['AD'][1]) / float(tumor['DP']) 
                                

            if len(record.ALT) == 1 and filter_col == 'PASS' and \
                (impact == 'HIGH' or impact == 'MODERATE' or is_splice_variant) and not wgs:
                
                if "missense_variant" in canonical_trans['Consequence'] and 'pathogenic' not in clinsig:
                    continue
                
                if is_splice_variant and 'pathogenic' not in clinsig:
                    continue
                
                if normal['DP'] and normal['AD']:
                    normal_dp = normal['DP']
                    normal_alt = normal['AD'][1]
                    normal_vaf = float(normal_alt) / float(normal_dp)

                    output_file.write('\t'.join(map(str, [record.CHROM, record.POS-1, record.POS,
                                                            record.REF, record.ALT, '', '', '', gene, ensembl_id, 
                                                            impact, canonical_trans['Consequence'], 
                                                            canonical_trans['Feature'], canonical_trans['HGVSc'],
                                                            canonical_trans['HGVSp'], ann_hotspot, normal_dp , normal_alt,
                                                            round(normal_vaf, 2), round(tumor_vaf, 2), clinsig, 
                                                            record.ID, gnomAD, brcaEx, oncogenicity, cgcann])) + "\n")
            
            # WGS filter
            if (wgs and len(record.ALT) == 1 and impact == 'HIGH' and not gene.startswith("HLA")) \
                or (wgs and 'pathogenic' in clinsig):
                
                if cgcann == '':
                    continue
                
                if normal['DP'] and normal['AD']:
                    normal_dp = normal['DP']
                    normal_alt = normal['AD'][1]
                    normal_vaf = float(normal_alt)/float(normal_dp)

                    output_file.write('\t'.join(map(str, [record.CHROM, record.POS-1, record.POS,
                                                            record.REF, record.ALT, '', '', '', gene, ensembl_id,
                                                            impact, canonical_trans['Consequence'], 
                                                            canonical_trans['Feature'], canonical_trans['HGVSc'],
                                                            canonical_trans['HGVSp'], ann_hotspot, normal_dp , normal_alt,
                                                            round(normal_vaf, 2), round(tumor_vaf, 2), clinsig, record.ID, 
                                                            gnomAD, brcaEx, oncogenicity, cgcann])) + "\n")

