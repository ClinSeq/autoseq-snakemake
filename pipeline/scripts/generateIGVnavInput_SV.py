#!/usr/bin/env python

import subprocess
import argparse
import os
import glob
import re
import pandas as pd


# structural variant color map for IGV SV representation
igv_color_map = {"DEL": "Non-coding_Transcript",
                 "DUP": "Truncating",
                 "INV": "Indel",
                 "TRA": "Nonsense",
                 "BND": "Nonsense"
                 }


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
            data = line.strip().split()
            chrom = data[0]
            start = data[1]
            end = data[2]
            sdid = data[3]
            svtype = data[4]
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
    header = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSUPPORT_READS\"" + \
                      " > " + output + "_" + vcftype + "_pass_gridss.mut.tmp"

    gridss_cmd = "zless " + input_vcf  + " | vawk '{ if($7 == \"PASS\")  print $1, $2, $2+1, \""+ SDID + '_gridss_' + vcftype +"\", I$SVTYPE, $5, I$VF}' " \
                 " >> " + output + "_" + vcftype +"_pass_gridss.mut.tmp"
    
    tmp_mut = output + "_" + vcftype + "_pass_gridss.mut.tmp"
    subprocess.call(" && ".join([header, gridss_cmd]), shell=True)
    get_igvcolortype(tmp_mut, "gridss")    


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
                    alt = gene_b
                    
                    if chrom_a != chrom_b:
                        sv_length = svtype
                    
                    if svtype == 'TRA':
                        event_1 = [chrom_a, start_a, end_a, sdid, igv_color_map[svtype], svtype, alt, support_reads]
                        event_2 = [chrom_b, start_b, end_b, sdid, igv_color_map[svtype], svtype, gene_a, support_reads]
                        events_list.extend([event_1, event_2])
                    else:
                        event = [chrom_a, start_a, end_a, sdid, igv_color_map[svtype], svtype, alt, support_reads]
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

    if not os.path.exists(output_dir + "/annotate_combined_sv.txt"):
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


def gene_annotation(chrom, start, end, genes, design):
    """
    Return gene name for given  SV event
    """
    gene = ''
    in_gene = ''

    try:
        # annotate gene name
        for ranges, gene_name in genes[chrom].items():
            if int(ranges[0]) - 20 <= int(start) <= int(ranges[1]) + 20 or int(ranges[0]) - 20 <= int(end) <= int(ranges[1]) + 20:
                gene = gene_name
                break
        if not gene:
            gene = 'None'

        # annotate pancancer info
        for pan_ranges, gene_name in design[chrom].items():
            if int(pan_ranges[0]) <= int(start) <= int(pan_ranges[1]) or int(pan_ranges[0]) <= int(end) <= int(pan_ranges[1]):
                in_gene = 'YES'
                break
        if not in_gene:
            in_gene = 'NO'

        return (gene, in_gene)
    except KeyError:
        return ('NA', 'NA')
        print("Warning! chromosome {chrom} is not valid".format(chrom=chrom))


def annotate_combined_sv(combined_file, genes, pancancer, output):
    """
    Parsing combined sv list and apply gene annotation for each SV
    """
    # output_file = open(output, 'w')
    summary_columns = ['CHROM_A', 'START_A', 'END_A', 'CHROM_B', 'START_B', 'END_B',
                       'IGV_COORD', 'SVTYPE', 'SV_LENGTH', 'SUPPORT_READS', 'TOOL', 'SDID', 'SAMPLE',
                       'GENE_A', 'IN_DESIGN_A', 'GENE_B', 'IN_DESIGN_B', "GENE_A-GENE_B-sorted"]
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
            sup_reads = data[8] if len(data) == 9 else '.'

            if ':' in data[7]:
                chrom_b = ''.join(list(filter(str.isdigit, data[7].split(':')[0])))
                start_b = ''.join(list(filter(str.isdigit, data[7].split(':')[1])))
                end_b = int(start_b) + 1
                
                if tool == "svcaller":
                    chrom_b = data[7].split(':')[0]
                    start_b = data[7].split(':')[1].split('-')[0]
                    end_b = data[7].split(':')[1].split('-')[1]
                
                igv_coord_b = chrom_b + ':' + str(start_b)
            else:
                chrom_b = 'NA'
                start_b = 'NA'
                end_b = 'NA'

            igv_coord = ' '.join([igv_coord_a, igv_coord_b])
            gene_a, in_gene_a = gene_annotation(chrom_a, start_a, end_a, genes, pancancer)

            if chrom_b != 'NA':
                gene_b, in_gene_b = gene_annotation(chrom_b, start_b, end_b, genes, pancancer)
            else:
                gene_b, in_gene_b = 'NA', 'NA'

            gene_a_b = [gene_a, gene_b]
            gene_a_b.sort()
            gene_a_b_sorted = ",".join(gene_a_b)
            summary_sv.append([chrom_a, start_a, end_a, chrom_b, start_b, end_b, igv_coord, svtype, int(end_a)-int(start_a), sup_reads, tool, sdid, sample, gene_a, in_gene_a, gene_b, in_gene_b, gene_a_b_sorted])
        summary_sv_df = pd.DataFrame(summary_sv, columns = summary_columns)
        summary_sv_df_sorted = summary_sv_df.sort_values(["IN_DESIGN_A", "IN_DESIGN_B", "GENE_A-GENE_B-sorted",
                                                        "CHROM_A", "START_A", "CHROM_B", "START_B", "TOOL"], 
                                                        ascending=[False, False, True, True, True, True, True, True])
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
    parser.add_argument('--targetBed', help="Bed file to annotate targets eg: pancancer")
    parser.add_argument('--sdid', help="SDID from analysis")
    parser.add_argument('--vcftype', help="somatic (or) germline vcf (only for svaba)")
    parser.add_argument('--tool', help="Tool name - Variant callers")
    parser.add_argument('--output', required=True,
                        help="output tab delimited file for IGVNav, format=output.mut")
    args = parser.parse_args()

    vcftype = args.vcftype
    input_file = args.input
    annotBed = args.annotBed
    sv_caller = args.tool
    sdid = args.sdid
    output = args.output
    target_bed = args.targetBed

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
        targets = load_bed(target_bed)
        annotate_combined_sv(combined_input, genes, targets, output)

