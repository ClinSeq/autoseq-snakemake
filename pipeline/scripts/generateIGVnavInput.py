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
import shutil
import json


def csq_parsing(csq, vcftype):
    # parsing CSQ taq from VeP annotation 
    # return canonical transcript csq taq as dict with corresponding keys

    csq_dicts = []
    can_trans = {}
    csq_keys_69 = ['Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','ALLELE_NUM','DISTANCE','STRAND','FLAGS','VARIANT_CLASS','SYMBOL_SOURCE','HGNC_ID','CANONICAL','TSL','APPRIS','CCDS','ENSP','SWISSPROT','TREMBL','UNIPARC','SOURCE','GENE_PHENO','SIFT','PolyPhen','DOMAINS','miRNA','HGVS_OFFSET','AF','AFR_AF','AMR_AF','EAS_AF','EUR_AF','SAS_AF','AA_AF','EA_AF','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF','MAX_AF','MAX_AF_POPS','CLIN_SIG','SOMATIC','PHENO','PUBMED','MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE','BrcaEx','BrcaEx_ClinicalSignificance']
    
    csq_keys_70 = ['Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','ALLELE_NUM','DISTANCE','STRAND','FLAGS','VARIANT_CLASS','SYMBOL_SOURCE','HGNC_ID','CANONICAL','TSL','APPRIS','CCDS','ENSP','SWISSPROT','TREMBL','UNIPARC','SOURCE','GENE_PHENO','SIFT','PolyPhen','DOMAINS','miRNA','HGVS_OFFSET','AF','AFR_AF','AMR_AF','EAS_AF','EUR_AF','SAS_AF','AA_AF','EA_AF','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF','MAX_AF','MAX_AF_POPS', 'FREQS','CLIN_SIG','SOMATIC','PHENO','PUBMED','MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE','BrcaEx','BrcaEx_ClinicalSignificance']
    
    csq_keys_85 = ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID', 'CANONICAL', 'MANE_SELECT', 'MANE_PLUS_CLINICAL', 'TSL', 'APPRIS', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'UNIPROT_ISOFORM', 'GIVEN_REF', 'USED_REF', 'BAM_EDIT', 'SOURCE', 'GENE_PHENO', 'SIFT', 'PolyPhen', 'DOMAINS', 'miRNA', 'HGVS_OFFSET', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'gnomADe_AF', 'gnomADe_AFR_AF', 'gnomADe_AMR_AF', 'gnomADe_ASJ_AF', 'gnomADe_EAS_AF', 'gnomADe_FIN_AF', 'gnomADe_NFE_AF', 'gnomADe_OTH_AF', 'gnomADe_SAS_AF', 'gnomADg_AF', 'gnomADg_AFR_AF', 'gnomADg_AMI_AF', 'gnomADg_AMR_AF', 'gnomADg_ASJ_AF', 'gnomADg_EAS_AF', 'gnomADg_FIN_AF', 'gnomADg_MID_AF', 'gnomADg_NFE_AF', 'gnomADg_OTH_AF', 'gnomADg_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'CLIN_SIG', 'SOMATIC', 'PHENO', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'TRANSCRIPTION_FACTORS', 'BrcaEx', 'BrcaEx_ClinicalSignificance']
     
    csq_keys_86 = ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID', 'CANONICAL', 'MANE_SELECT', 'MANE_PLUS_CLINICAL', 'TSL', 'APPRIS', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'UNIPROT_ISOFORM', 'GIVEN_REF', 'USED_REF', 'BAM_EDIT', 'SOURCE', 'GENE_PHENO', 'SIFT', 'PolyPhen', 'DOMAINS', 'miRNA', 'HGVS_OFFSET', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'gnomADe_AF', 'gnomADe_AFR_AF', 'gnomADe_AMR_AF', 'gnomADe_ASJ_AF', 'gnomADe_EAS_AF', 'gnomADe_FIN_AF', 'gnomADe_NFE_AF', 'gnomADe_OTH_AF', 'gnomADe_SAS_AF', 'gnomADg_AF', 'gnomADg_AFR_AF', 'gnomADg_AMI_AF', 'gnomADg_AMR_AF', 'gnomADg_ASJ_AF', 'gnomADg_EAS_AF', 'gnomADg_FIN_AF', 'gnomADg_MID_AF', 'gnomADg_NFE_AF', 'gnomADg_OTH_AF', 'gnomADg_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'FREQS', 'CLIN_SIG', 'SOMATIC', 'PHENO', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'TRANSCRIPTION_FACTORS', 'BrcaEx', 'BrcaEx_ClinicalSignificance']

    csq_keys_map = { 
        69 : csq_keys_69,
        70 : csq_keys_70,
        85 : csq_keys_85,
        86 : csq_keys_86
    }

    for transcript in csq:
        csq_len = len(transcript.split('|'))
        csq_keys = csq_keys_map[csq_len] 

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

 
# Parsing Commandline Arguments using argparse
#

parser = argparse.ArgumentParser()
parser.add_argument('vcf', help="Input VCF file for annotation")
parser.add_argument('oncokb', help="OncoKB - all variants tab demilited file")
parser.add_argument('vcftype', help="somatic (or) germline vcf")
parser.add_argument('--wgs', action='store_true', default=False, help="tag to use WGS filter")
parser.add_argument('--cgc', default=False, help="Cancer Gene Census Annotation")
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

###############generate IGVNAV symblins################################################################

def create_symlink(travers_dir_name, src_dir, igvnav_dirname_dst, suffix):
    "Recursively Traverse through the directory and create symlink"
    for root, dirs, files in os.walk(travers_dir_name):
        for each_file in files:
            if each_file.endswith(suffix) and not os.path.exists(os.path.join(igvnav_dirname_dst,each_file)):
                os.symlink(os.path.join(root,each_file), os.path.join(igvnav_dirname_dst,each_file))
    return 


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
    output_file.write('\t'.join(['CHROM','START','END','REF','ALT', 'CALL', 'TAG', 'NOTES', 'GENE',  'ENSEMBLID', 'IMPACT', 'CONSEQUENCE', 'TRANSCRIPT', 'HGVSc', 'HGVSp', 'T_DP', 'T_ALT', 'T_VAF', 'N_DP', 'N_ALT', 'N_VAF', 'CLIN_SIG', 'RSID', 'gnomAD', 'BRCAEx', 'OncoKB', 'CGC_ANN', 'NUM_TOOLS']) + "\n")
elif vcftype == "germline":
    output_file.write('\t'.join(['CHROM','START','END','REF','ALT', 'CALL', 'TAG', 'NOTES', 'GENE', 'ENSEMBLID', 'IMPACT', 'CONSEQUENCE', 'TRANSCRIPT', 'HGVSc','HGVSp', 'N_DP', 'N_ALT', 'N_VAF', 'T_VAF', 'CLIN_SIG', 'RSID', 'gnomAD', 'BRCAEx', 'OncoKB', 'CGC_ANN']) + "\n")

for record in vcf_reader:
    try:
        # skip germline variants which doesn't have CSQ annotations
        if "CSQ" not in record.INFO:
            continue 

        canonical_trans = csq_parsing(record.INFO['CSQ'], vcftype)
        gene = canonical_trans['SYMBOL']
        ensembl_id = canonical_trans['Gene']
        aa = canonical_trans['Amino_acids'].split('/')
        protein_position = canonical_trans['Protein_position'].split('/')
        clinsig =  canonical_trans['CLIN_SIG']
        impact = canonical_trans['IMPACT']
        brcaEx = canonical_trans['BrcaEx_ClinicalSignificance']
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

    # processing somatic vcf file
    if vcftype == "somatic":
        normal = record.genotype('NORMAL')
        normal_dp = sum(normal['DP4'])
        normal_alt = sum(normal['DP4'][2:])
        normal_vaf = normal['VAF']

        tumor = record.genotype('TUMOR')
        tumor_dp = sum(tumor['DP4'])
        tumor_alt = sum(tumor['DP4'][2:])
        tumor_vaf = tumor['VAF']

        num_tools = int(record.INFO['NUM_TOOLS'])
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
                                                  canonical_trans['HGVSc'], canonical_trans['HGVSp'], tumor_dp, tumor_alt, 
                                                  tumor_vaf, normal_dp, normal_alt, normal_vaf, 
                                                  clinsig, rsid, gnomAD, brcaEx, oncogenicity, cgcann, num_tools])) + "\n")

        # filter for WGS samples
        if wgs and num_tools >= 2 and tumor_alt >= 5 and (impact == 'HIGH' or impact == 'MODERATE'):
            
            output_file.write('\t'.join(map(str, [record.CHROM, record.POS-1, record.POS,
                                                  record.REF, record.ALT, '', '', '', gene, ensembl_id,
                                                  impact, canonical_trans['Consequence'], canonical_trans['Feature'],
                                                  canonical_trans['HGVSc'], canonical_trans['HGVSp'], tumor_dp, tumor_alt, 
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
                                                        canonical_trans['HGVSp'], normal_dp , normal_alt,
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
                                                        canonical_trans['HGVSp'], normal_dp , normal_alt,
                                                        round(normal_vaf, 2), round(tumor_vaf, 2), clinsig, record.ID, 
                                                        gnomAD, brcaEx, oncogenicity, cgcann])) + "\n")
                

## adding vardict indels into IGVNav 

if args.vardict and vcftype == "somatic":
    
    for record in vardict_vcf:
        # filter snvs 
        if record.INFO['TYPE'] == "SNV":
            continue

        ref = record.REF
        alt = record.ALT[0]

        try:
            # filter small indels - length less than 5 
            if len(ref) > 5 or len(alt) > 5:
                # to avoid duplicates    
                tmp_str = "-".join(map(str, [record.CHROM, record.POS, ref, alt]))
                if tmp_str in variants:
                    continue

                canonical_trans = csq_parsing(record.INFO['CSQ'], vcftype)

                gene = canonical_trans['SYMBOL']
                aa = canonical_trans['Amino_acids'].split('/')
                protein_position = canonical_trans['Protein_position'].split('/')
                clinsig =  canonical_trans['CLIN_SIG']
                gnomAD = canonical_trans['gnomAD_AF']
                brcaEx = canonical_trans['BrcaEx_ClinicalSignificance']
                impact = canonical_trans['IMPACT']
                oncogenicity = ''

                # Oncogenicity annotation from OncoKB
                if gene in OncoKB_lookup:
                    if len(aa) > 1 and len(protein_position) > 1:
                        protein_change = aa[0] +  protein_position[0] + aa[1]
                        if protein_change in OncoKB_lookup[gene]:
                            oncogenicity = OncoKB_lookup[gene][protein_change]
                
                normal = [sam for sam in record.samples if '-N-' in sam.sample][0]
                tumor = [sam for sam in record.samples if '-CFDNA-' in sam.sample or '-T-' in sam.sample][0]
            
                normal_dp = sum(normal['AD'])
                normal_alt = normal['AD'][1]
                normal_vaf = normal['AF']
                tumor_dp = sum(tumor['AD'])
                tumor_alt = tumor['AD'][1]
                tumor_vaf = tumor['AF']

                output_file.write('\t'.join(map(str, [record.CHROM, record.POS-1, record.POS,
                                                    record.REF, record.ALT, '', '', '', gene, 
                                                    impact, canonical_trans['Consequence'], 
                                                    canonical_trans['HGVSp'], tumor_dp, tumor_alt, 
                                                    tumor_vaf, normal_dp, normal_alt, normal_vaf, 
                                                    clinsig, gnomAD, brcaEx, oncogenicity, 1])) + "\n")
        except TypeError:
            print(record)
            pass
