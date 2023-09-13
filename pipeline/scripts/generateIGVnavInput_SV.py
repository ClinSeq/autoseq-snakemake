#!/usr/bin/env python

import subprocess
import argparse
import os
import glob
import re
import json
import pandas as pd
import pyranges as pr
import vcf



valid_chromo = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]


# structural variant color map for IGV SV representation
igv_color_map = {"DEL": "Non-coding_Transcript",
                 "DUP": "Truncating",
                 "INS": "Truncating",
                 "INV": "Indel",
                 "TRA": "Nonsense",
                 "CTX": "Nonsense",
                 "BND": "Nonsense"
                 }


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
            cgc_ann[gene_symbol] = (ann, ensemblID)
    
    return cgc_ann


def load_exons(gtfpath):
    """
    function to load gtf file
    """
    colnames = ["Chromosome", "source", "type", "Start", "End", "score", "Strand", "name", "feature"]
    exons_df = pd.read_csv(gtfpath, header = None, sep = "\t", names = colnames)

    return pr.PyRanges(exons_df)


def get_igvcolortype(mutfile, tool):
    """
    function to add type for IGV representation
    """
    outfile = open(mutfile.replace(".tmp", ""), "w")
    
    if tool == "svaba":
        outfile.write("\t".join(['CHROM','START','END','SDID','TYPE','SVTYPE','ALT','SUPPORT_normal', 'SUPPORT_tumor']) + '\n')
    else:
        outfile.write("\t".join(['CHROM','START','END','SDID','TYPE','SVTYPE','ALT','SUPPORT_READS']) + '\n')

    with open(mutfile, 'r') as fh:
        for line in fh:
            if line.startswith("CHROM"):
                continue
            data = line.strip().split('\t')
            chrom = data[0]
            start = data[1]
            end = data[2]
            sdid = data[3]
            svtype = data[4]
            if svtype == '':
                svtype = 'BND'
            igvtype = igv_color_map[svtype]
            alt = data[5]
            if tool == "svaba":
                support_reads = "\t".join([data[6], data[7]])
            else:
                support_reads = data[6]
            outfile.write("\t".join([chrom, start, end, sdid, igvtype, svtype, alt, support_reads]) + "\n")
    
    subprocess.call("rm {}".format(mutfile), shell=True)


def parse_svaba(input_vcf, SDID, output, vcftype):
    """
    vawk '{print $1, $2, $2+1, "P-00356971_svaba", "BND", $5, S$*$AD, S$*$DP}'
    """
    header = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSUPPORT_normal\tSUPPORT_tumor \
             \tDPnormal\tDPtumor\tGENES\"" + " > " + output + "_" + vcftype + "_svaba.mut"
    svaba_cmd = "vawk '{print $1, $2, $2+1, \"" + SDID + '_svaba_' + vcftype +"\",I$SVTYPE, $5, S$*$AD,S$*$DP}'" + \
                " " + input_vcf + " >> " + output + "_" + vcftype + "_svaba.mut.tmp"
    
    tmp_mut = output + "_" + vcftype + "_svaba.mut.tmp"
    subprocess.call(" && ".join([header, svaba_cmd]), shell=True)
    get_igvcolortype(tmp_mut, "svaba")


def parse_lumpy(input_vcf, SDID, output, vcftype):
    """
    ##CHROM START   END SDID    SVTYPE  ALT SVLENGTH    SUPPORT_READS    NOTES   GENES
    vawk '{if ((I$SVLEN>1000 || I$SVLEN<-1000) && $1 != "hs37d5" && $1 !~ "GL" && $5 !~ "hs37d5" &&
    I$SU>50 && I$SVTYPE ~ "BND") print $1, $2, $2+1, "P-00356971_lumpy", I$SVTYPE, $5, "NA", I$SU
    else if ((I$SVLEN>1000 || I$SVLEN<-1000) && $1 != "hs37d5" && $1 !~ "GL" && $5 !~ "hs37d5"
    && I$SU>50 && I$SVTYPE !~ "BND") print $1, $2, I$END, "P-00356971_lumpy", I$SVTYPE, $5, I$SVLEN, I$SU}
    """
    header1k_sup_50 = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSUPPORT_READS\"" + \
                      " > " + output + "_lumpy_len1k_SU50.mut.tmp"
    header500_sup_24 = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSUPPORT_READS\"" + \
                       " > " + output + "_lumpy_len500_SU24.mut.tmp"

    len1k_sup_50 = "vawk '{if ((I$SVLEN>1000 || I$SVLEN<-1000) && $1 != \"hs37d5\" && $1 !~ \"GL\" && $5 !~ \"hs37d5\" && I$SU>50 && I$SVTYPE ~ \"BND\") print $1, $2, $2+1, \""+ SDID + '_lumpy_' + vcftype +"\", I$SVTYPE, $5, \"NA\", I$SU ;"  + \
                " else if ((I$SVLEN>1000 || I$SVLEN<-1000) && $1 != \"hs37d5\" && $1 !~ \"GL\" && $5 !~ \"hs37d5\" && I$SU>50 && I$SVTYPE !~ \"BND\") print $1, $2, I$END, \""+ SDID + '_lumpy_' + vcftype +"\", I$SVTYPE, $5, I$SU}' " + input_vcf + \
                " >> " + output + "_lumpy_len1k_SU50.mut.tmp"

    len500_sup_24 = "vawk '{if ((I$SVLEN>500 || I$SVLEN<-500) && $1 != \"hs37d5\" && $1 !~ \"GL\" && $5 !~ \"hs37d5\" && I$SU>24 && I$SVTYPE ~ \"BND\") print $1, $2, $2+1, \"" + SDID + '_lumpy_' + vcftype + "\", I$SVTYPE, $5, I$SU ;" + \
                " else if ((I$SVLEN>500 || I$SVLEN<-500) && $1 != \"hs37d5\" && $1 !~ \"GL\" && $5 !~ \"hs37d5\" && I$SU>24 && I$SVTYPE !~ \"BND\") print $1, $2, I$END, \"" + SDID + '_lumpy_' + vcftype + "\", I$SVTYPE, $5, I$SU }' " + input_vcf + \
                " >> " + output + "_lumpy_len500_SU24.mut.tmp"
  
    # cmd = "awk 'NR>1 {OFS=\"\\t\";print $1, $2, $3, $5,\"lumpy\", $4, \"" + vcftype + "\", $6, $7}' " + output + "_lumpy_len500_SU24.mut" +  " >> " + output_dir + "/annotate_combined_sv.txt"

    tmp_len500_SU24 = output + "_lumpy_len500_SU24.mut.tmp"
    tmp_len1k_SU50 = output + "_lumpy_len1k_SU50.mut.tmp"
    subprocess.call(" && ".join([header1k_sup_50, header500_sup_24, len1k_sup_50, len500_sup_24]), shell=True)
    get_igvcolortype(tmp_len500_SU24, "lumpy")
    get_igvcolortype(tmp_len1k_SU50, "lumpy")


def parse_gridss(input_vcf, SDID, output, vcftype):
    """
    GRIDSS - vcf parsing 
    """
    # header = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSUPPORT_READS\"" + \
    #                   " > " + output + "_" + vcftype + "_pass_gridss.mut.tmp"

    # gridss_cmd = "less " + input_vcf  + " | vawk '{ if($7 == \"PASS\")  print $1, $2, $2+1, \""+ SDID + '_gridss_' + vcftype +"\", I$SIMPLE_TYPE, $5, I$VF}' " \
    #              " >> " + output + "_" + vcftype +"_pass_gridss.mut.tmp"
    
    # tmp_mut = output + "_" + vcftype + "_pass_gridss.mut.tmp"
    # subprocess.call(" && ".join([header, gridss_cmd]), shell=True)
    # get_igvcolortype(tmp_mut, "gridss")
    output_fname = "_".join([output, vcftype, "pass_gridss.mut"])
    outfile = open(output_fname, "w")

    outfile.write("\t".join(['CHROM','START','END','SDID','TYPE','SVTYPE','ALT','SUPPORT_READS']) + '\n')

    vcf_reader = vcf.Reader(filename=input_vcf)
    events = set()
    sdid = SDID + '_gridss_' + vcftype
    for record in vcf_reader:
        if record.FILTER != []:
            continue
        curr_event = record.INFO['EVENT']
        
        if curr_event not in events \
            and 'SIMPLE_TYPE' in record.INFO:

            svtype = record.INFO['SIMPLE_TYPE']
            chrom = record.CHROM
            start = record.POS
            end = record.POS + 1
            alt = str(record.ALT[0])
            
            if svtype != 'TRA':
                events.add(curr_event)
                end = ''.join(list(filter(str.isdigit, alt.split(':')[1])))

            support_reads = record.INFO['VF']

            outfile.write("\t".join(map(str, [chrom, start, end, sdid, 
                                     igv_color_map[svtype], svtype, alt, support_reads])) + "\n")
    

def parse_gtf(gtf, sdid, vcftype):
    """
    Parsing svcaller output - gtf
    """
    if 'DEL' in gtf:
        svtype = 'DEL'
    elif 'DUP' in gtf:
        svtype = 'DUP'
    elif 'TRA' in gtf:
        svtype = 'TRA'
    elif 'INV' in gtf:
        svtype = 'INV'

    sdid = sdid + '_svcaller_' + vcftype
    events_list = []
    gtf_set = set()
    with open(gtf, 'r') as gtf_fh:
        for line in gtf_fh.readlines():
            if line.strip():
                feature = line.strip().split('\t')[2]
                gene_id = line.strip().split('\t')[8]
                support_reads = line.strip().split('\t')[5]
                if feature == 'exon' and gene_id not in gtf_set:
                    gtf_set.add(gene_id)
                    genes_coords = re.search('gene_id "(.*)"; transcript_id', gene_id).group(1)
                    gene_a = genes_coords.split(',')[0]
                    chrom_a = gene_a.split(':')[0]
                    start_a = re.search(':(\d+)', gene_a).group(1)
                    end_a = re.search('-(\d+)', gene_a).group(1)
                    gene_b = genes_coords.split(',')[1]
                    chrom_b = gene_b.split(':')[0]
                    start_b = re.search(':(\d+)', gene_b).group(1)
                    end_b = re.search('-(\d+)', gene_b).group(1)
                    sv_length = int(end_b) - int(start_a)
                    alt = ",".join([gene_a, gene_b])
                    
                    if chrom_a != chrom_b:
                        sv_length = svtype
                    
                    if svtype == 'TRA':
                        event_1 = [chrom_a, start_a, end_a, sdid, igv_color_map[svtype], svtype, gene_b, support_reads]
                        event_2 = [chrom_b, start_b, end_b, sdid, igv_color_map[svtype], svtype, gene_a, support_reads]
                        events_list.extend([event_1, event_2])
                    else:
                        event = [chrom_a, start_a, end_b, sdid, igv_color_map[svtype], svtype, alt, support_reads]
                        events_list.append(event)
                    
    return events_list


def parse_svcaller(input_dir, SDID, output, vcftype):
    """
    parse sdid and extract corresponding gtf files to process
    """
    sdid = "-".join(SDID.split("-")[0:4])
    gtf_files = glob.glob(input_dir + "/" + sdid + "-*.gtf")
    mut_file = output + "/" + SDID + "_svcaller.mut"
    sdid = re.search("P-[A-Za-z0-9]*", SDID).group()
    events = []
    for gtf in gtf_files:
        events.extend(parse_gtf(gtf, sdid, vcftype))

    with open(mut_file, 'w') as mut_fh:
        mut_fh.write("\t".join(['CHROM','START','END','SDID','TYPE','SVTYPE','ALT', 'SUPPORT_READS']) + '\n')
        for event in events:
            mut_fh.write('\t'.join(map(str, event)) + '\n')
 

def combine_mut(input_dir, output_dir):
    """
    Function to combine all mut files prepared for IGVNav
    """
    files = glob.glob(input_dir + "/*.mut")
    cmd = []

    if os.path.exists(output_dir + "/annotate_combined_sv.txt"):
        print("annotate_combined_sv.txt file already exists!")
        return output_dir + "/annotate_combined_sv.txt"
    
    header2 = "echo \"CHROM\tSTART\tEND\tSVTYPE\tTOOL\tSDID\tSAMPLE\tALT\tSUPPORT_READS\"" + " >> " + output_dir + "/annotate_combined_sv.txt"
    subprocess.call(header2, shell=True)

    for file in files:
        vcftype = ''
        sup_reads = ''
        filebase = os.path.basename(file)
        if 'lumpy_len500_SU24' in file:
            cmd.append("awk -F'\\t' 'NR>1 {OFS=\"\\t\";print $1, $2, $3, $6,\"lumpy\", $4, \"somatic\", $7, $8}' " \
                        + file +  " >> " + output_dir + "/annotate_combined_sv.txt")
        elif 'svaba.mut' in file:
            vcftype = 'somatic' if 'somatic' in filebase else 'germline'
            sup_reads = '$9' if vcftype == 'SOMATIC' else '$8'
            cmd.append("awk -F'\\t' 'NR>1 {OFS=\"\\t\";print $1, $2, $3, $6,\"svaba\", $4, \"" + vcftype + "\", $7, " + sup_reads + "}' " +\
                file + " >> " + output_dir + "/annotate_combined_sv.txt")
        elif 'svcaller.mut' in file:
            vcftype = 'cfdna' if '-CFDNA-' in filebase else 'tumor' if '-T-' in filebase else 'germline'
            cmd.append("awk -F'\\t' 'NR>1 {OFS=\"\\t\";print $1, $2, $3, $6,\"svcaller\", $4, \"" + vcftype + "\", $7, $8}' " \
                        + file +  " >> " + output_dir + "/annotate_combined_sv.txt")
        elif 'gridss.mut' in file:
            vcftype = 'somatic' if 'somatic' in filebase else 'germline'            
            cmd.append("awk -F'\\t' 'NR>1 {OFS=\"\\t\";print $1, $2, $3, $6,\"gridss\", $4, \"" + vcftype + "\", $7, $8}' " \
                        + file +  " >> " + output_dir + "/annotate_combined_sv.txt")

    subprocess.call(" && ".join(cmd), shell=True)

    return output_dir + "/annotate_combined_sv.txt"


def load_bed(bed_file):
    """
    Loading genes from genes.bed file for SV annotations
    """
    genes = {}
    with open(bed_file, 'r') as genes_fh:
        genes_db = genes_fh.readlines()
        for each_entry in genes_db:
            data = each_entry.strip().split('\t')
            chrom = data[0]
            start = data[1]
            end = data[2]
            gene = data[3]

            if chrom in genes:
                genes[chrom].update({(start, end): gene})
            else:
                genes[chrom] = {(start, end): gene}
    return genes


def gene_annotation(chrom, start, end, genes):
    """
    Return gene name for given  SV event
    """
    gene = ''

    try:
        # annotate gene name
        for ranges, gene_name in genes[chrom].items():
            if int(ranges[0]) - 20 <= int(start) <= int(ranges[1]) + 20 \
                or int(ranges[0]) - 20 <= int(end) <= int(ranges[1]) + 20:
                gene = gene_name
                break
        
        if not gene:
            gene = 'None'

        return gene
    except KeyError:
        print("Warning! chromosome {chrom} is not valid".format(chrom=chrom))
        return 'NA'
        

def check_targets(chrom, start, end, targets):
    """
    function to filter out SVs using target intervals
    """
    if chrom not in targets:
        return False

    for i in targets[chrom]:
        if int(i["START"]) - 150 <= int(start) <= int(i["END"]) + 150 \
            or int(i["START"]) - 150 <= int(end) <= int(i["END"]) + 150:
            return True

    return False



def annotate_combined_sv(combined_file, genes, targets, capture, cgc_ann, output, exons):
    """
    Parsing combined sv list and apply gene annotation for each SV
    """
    # output_file = open(output, 'w')
    summary_columns = ['CHROM_A', 'START_A', 'END_A', 'CHROM_B', 'START_B', 'END_B',
                       'IGV_COORD', 'SVTYPE', 'SV_LENGTH', 'SUPPORT_READS', 'TOOL', 'SDID', 'SAMPLE',
                       'GENE_A', 'GENE_B', "GENE_A-GENE_B-sorted", "CGC_ANN", "CURATOR"]
    summary_sv = list()
    with open(combined_file, 'r') as fh:
        header = fh.readline()
        for line in fh.readlines():
            data = line.strip().split('\t')
            chrom_a = data[0]
            start_a = data[1]
            end_a = data[2]
            igv_coord_a = chrom_a + ':' + str(start_a)
            igv_coord_b = ''
            svtype = data[3]
            tool = data[4]
            sdid = data[5].split('_')[0]
            sample = data[6]
            alt = data[7]
            sup_reads = data[8] if len(data) == 9 else '.'
            svlength = 'NA'
            chrom_b = 'NA'

            if ':' in alt:
                chrom_b = ''.join(list(filter(str.isdigit, alt.split(':')[0])))
                end_b = ''.join(list(filter(str.isdigit, alt.split(':')[1])))
                start_b = int(end_b) - 1

                if 'X' in alt:
                    chrom_b = 'X'
                elif 'Y' in alt:
                    chrom_b = 'Y'
                
                if tool == "svcaller":
                    bps = data[7].split(",")
                    if len(bps) == 2:
                        bp_a, bp_b = bps
                        start_a = bp_a.split(':')[1].split('-')[0]
                        end_a = bp_a.split(':')[1].split('-')[1]
                        chrom_b = bp_b.split(':')[0]
                        start_b = bp_b.split(':')[1].split('-')[0]
                        end_b = bp_b.split(':')[1].split('-')[1]
                    else:
                        chrom_b = bps[0].split(':')[0]
                        start_b = bps[0].split(':')[1].split('-')[0]
                        end_b = bps[0].split(':')[1].split('-')[1]

                igv_coord_b = chrom_b + ':' + str(start_b)
            
            # Filtered invalid chromosome and decoy events 
            if chrom_a not in valid_chromo or chrom_b not in valid_chromo:
                continue
            
            # this section for backward compatability
            if svtype != 'TRA' and tool != 'svcaller' and tool != 'gridss':
                chrom_b = 'NA'
                start_b = 'NA'
                end_b = 'NA'
            
            # GRIDSS coord split 
            if svtype != 'TRA' and tool == 'gridss':
                end_a = int(start_a)
                start_a = int(start_a) - 1

            igv_coord = ' '.join([igv_coord_a, igv_coord_b])
            gene_a = gene_annotation(chrom_a, start_a, end_a, genes)

            if chrom_b != 'NA':
                svlength = abs(int(end_b)-int(start_a)) if chrom_a == chrom_b else 'NA'
                gene_b = gene_annotation(chrom_b, start_b, end_b, genes)
            else:
                gene_b = 'NA'
            
            cgcann = set()
            cgcann.add(cgc_ann[gene_a][0] if gene_a in cgc_ann else None) 
            cgcann.add(cgc_ann[gene_b][0] if gene_b in cgc_ann else None)

            if capture == "WG":
                curator = "NO"
                if list(filter(None, cgcann)) != []:
                    curator = "YES"
            else:
                if check_targets(chrom_a, start_a, end_a, targets) or \
                        check_targets(chrom_b, start_b, end_b, targets):
                    curator = "YES"
                else:
                    curator = "NO"
            
                

            if tool == 'gridss' and chrom_b == 'NA':
                svlength = abs(int(end_a)-int(start_a))
                # calculation for gridss INS svlength
                if svtype == "INS":
                    alt_seq = ''.join(list(filter(str.isalpha, alt)))
                    svlength = len(alt_seq)

            gene_a_b = [gene_a, gene_b]
            gene_a_b.sort()
            gene_a_b_sorted = ",".join(gene_a_b)

            summary_sv.append([chrom_a, start_a, end_a, chrom_b, start_b, end_b, igv_coord, svtype,
                                svlength, sup_reads, tool, sdid, sample, gene_a, 
                                gene_b, gene_a_b_sorted, ",".join(list(filter(None, cgcann))), curator])
        

        svs_df = pd.DataFrame(summary_sv, columns = summary_columns)
        # hard filter for gridss germline svs
        gf_idx = svs_df[(svs_df['SAMPLE'] == "germline") &  (svs_df['SUPPORT_READS'].astype('int') < 40)].index
        svs_df.loc[list(gf_idx), "CURATOR"] = "NO"

        # checking exons overlaps
        svs_df['idx'] = svs_df.index
        # fetching indices except TRA and off target svs
        t_idx = set(svs_df[(svs_df['SVTYPE'] != "TRA") & (svs_df['CURATOR'] == "YES")].index)
        
        # generating pyranges for gridss and svcaller
        svs_gridss_df = svs_df[["CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "SVTYPE", "idx", "CURATOR", "TOOL"]].rename(columns = {"CHROM_A": "Chromosome", "START_A": "Start", "END_B": "End"})
        svs_gridss_pr = pr.PyRanges(svs_gridss_df.loc[(svs_gridss_df['SVTYPE'] != "TRA") & (svs_gridss_df['CURATOR'] == "YES") & (svs_gridss_df['TOOL'] == "gridss") ])
        svs_svcaller_df = svs_df[["CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "SVTYPE", "idx", "CURATOR", "TOOL"]].rename(columns = {"CHROM_A": "Chromosome", "START_A": "Start", "END_B": "End"})
        svs_svcaller_pr = pr.PyRanges(svs_svcaller_df.loc[(svs_svcaller_df['SVTYPE'] != "TRA") & (svs_svcaller_df['CURATOR'] == "YES") & (svs_svcaller_df['TOOL'] == "svcaller")])

        # find exons overlapping svs
        hits_idx = set(svs_gridss_pr.intersect(exons).idx)
        hits_svc = svs_svcaller_pr.intersect(exons)
        _idx = set(hits_svc.idx) if hits_svc else set()
        hits_idx.update(_idx)
        
        # filter non overlapping svs by index
        filter_idx = t_idx - hits_idx
        svs_df.loc[list(filter_idx), "CURATOR"] = "NO"
        del svs_df['idx']

        summary_sv_df_sorted = svs_df.sort_values(["GENE_A-GENE_B-sorted", "CHROM_A", "START_A", "CHROM_B", "START_B", "TOOL"], 
                                                        ascending=[True, True, True, True, True, True])
        summary_sv_df_sorted.to_csv(output, sep = "\t", encoding = 'utf-8', index = False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=
        'A MUT file (.mut) is a tab-delimited text file that lists mutations. \
        The first row contains column headings and each subsequent row identifies a mutation. \
        IGV ignores the column headings.It reads the first five columns as shown below and \
        ignores all subsequent columns:  \
        1. chromosome \
        2. start location (location of the first base pair in the mutated region) \
        3. end location (location of the last base pair in the mutated region) \
        4. sample or patient ID \
        5. mutation type (for example, Synonymous, Missense, Nonsense, Indel, etc.)')
    parser.add_argument('--input', required=True, help="Input VCF or tab-delimited file")
    parser.add_argument('--annotBed', help="UCSC hg19 genes bed file with chrom, start, \
                        end and genesymbol")
    parser.add_argument('--target', help="capture kit ID and json file contains list of target genes interval",
                        nargs=2)
    parser.add_argument('--sdid', help="SDID from analysis")
    parser.add_argument('--vcftype', help="somatic (or) germline vcf (only for svaba)")
    parser.add_argument('--tool', help="Tool name - Variant callers")
    parser.add_argument('--cgc', help="Cancer Gene Census Annotation ")
    parser.add_argument('--exons', help="human exons coordinates as gtf file")
    parser.add_argument('--output', required=True,
                        help="output tab delimited file for IGVNav, format=output.mut")
    args = parser.parse_args()

    vcftype = args.vcftype
    input_file = args.input
    annotBed = args.annotBed
    sv_caller = args.tool
    sdid = args.sdid
    output = args.output
 
    if args.target:
        capture_kit, target_json = args.target

    if args.cgc:
        cgc_ann = loadCGC(args.cgc)
    
    if args.exons:
        exons = load_exons(args.exons)
    
    output_dir = os.path.dirname(output)

    if sv_caller == 'lumpy':
        parse_lumpy(input_file, sdid, output, vcftype)
    elif sv_caller == 'svaba':
        parse_svaba(input_file, sdid, output, vcftype)
    elif sv_caller == 'svcaller':
        parse_svcaller(input_file, sdid, output, vcftype)
    elif sv_caller == 'gridss':
        parse_gridss(input_file, sdid, output, vcftype)

    if annotBed:
        combined_input = combine_mut(input_file, output_dir)
        genes = load_bed(annotBed)
        fh = open(target_json, 'r')
        targets = json.load(fh)
        if capture_kit in targets:
            targets = targets[capture_kit]
        annotate_combined_sv(combined_input, genes, targets, 
                             capture_kit, cgc_ann, output, exons)

