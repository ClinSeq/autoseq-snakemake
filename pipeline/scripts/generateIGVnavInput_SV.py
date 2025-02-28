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


def parse_gridss(input_vcf, SDID, output, vcftype):
    """
    GRIDSS - vcf parsing 
    """
    output_fname = "_".join([output, vcftype, "pass_gridss.mut"])
    outfile = open(output_fname, "w")

    outfile.write("\t".join(['CHROM','START','END','SDID','TYPE','SVTYPE','ALT','SUPPORT_READS', 'VAF']) + '\n')

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
            if vcftype == "somatic":
                vaf = ",".join(map(str, record.INFO['TAF']))
            else:
                vaf = str(record.samples[0]['AF'])

            outfile.write("\t".join(map(str, [chrom, start, end, sdid, 
                                     igv_color_map[svtype], svtype, alt, support_reads, vaf])) + "\n")
    

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
    
    header2 = "echo \"CHROM\tSTART\tEND\tSVTYPE\tTOOL\tSDID\tSAMPLE\tALT\tSUPPORT_READS\tVAF\"" + " >> " + output_dir + "/annotate_combined_sv.txt"
    subprocess.call(header2, shell=True)

    for file in files:
        vcftype = ''
        sup_reads = ''
        filebase = os.path.basename(file)
        if 'svcaller.mut' in file:
            vcftype = 'cfdna' if '-CFDNA-' in filebase else 'tumor' if '-T-' in filebase else 'germline'
            cmd.append("awk -F'\\t' 'NR>1 {OFS=\"\\t\";print $1, $2, $3, $6,\"svcaller\", $4, \"" + vcftype + "\", $7, $8, \"NA\"}' " \
                        + file +  " >> " + output_dir + "/annotate_combined_sv.txt")
        elif 'gridss.mut' in file:
            vcftype = 'somatic' if 'somatic' in filebase else 'tumor' if 'tumor' in filebase else 'germline'            
            cmd.append("awk -F'\\t' 'NR>1 {OFS=\"\\t\";print $1, $2, $3, $6,\"gridss\", $4, \"" + vcftype + "\", $7, $8, $9}' " \
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

    query_pr = pr.PyRanges(chromosomes =  chrom,
                           starts = [start],
                           ends = [end])
    
    hits = genes.intersect(query_pr)
    gene_names = ",".join(hits.Name) if hits else 'NA'

    return gene_names
        

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


def add_cgcann(gene, cgc):
    """
    Add CGC annotation
    """
    if ',' in gene:
        genes = gene.split(',') 
        cgc_ann = [cgc[i][0] for i in genes if i in cgc]
    else:
        cgc_ann = [cgc[gene][0] if gene in cgc else None]
            
    if not cgc_ann:
        return [None]

    return cgc_ann


def annotate_curate_info(curation_ann_df, gene, type, project):
    """
    function to annotate cancer specific curation information
    """
    curation_ann_df = curation_ann_df[(curation_ann_df['prefix'] == project) & (curation_ann_df['type'] == type)]
    curation_ann_df = curation_ann_df.set_index('gene')['comment'].to_dict()
    comments = []
    if ',' in gene:
        genes = gene.split(',')
        comments = [" ".join([i, curation_ann_df[i]]) for i in genes if i in curation_ann_df]
    else:
        comments = [" ".join([gene, curation_ann_df[gene]]) if gene in curation_ann_df else None]
    
    return ",".join(comments)


def annotate_combined_sv(combined_file, genes, targets, capture, cgc_ann, output, exons, curation_ann):
    """
    Parsing combined sv list and apply gene annotation for each SV
    """
    # output_file = open(output, 'w')
    summary_columns = ['CHROM_A', 'START_A', 'END_A', 'CHROM_B', 'START_B', 'END_B',
                       'IGV_COORD', 'SVTYPE', 'SV_LENGTH', 'SUPPORT_READS', 'VAF', 'TOOL', 'SDID',
                       'SAMPLE', 'GENE_A', 'GENE_B', "GENE_A-GENE_B-sorted", "CGC_ANN", "CURATOR", "CURATE"]
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
            sup_reads = data[8]
            vaf = data[9]
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
            
            cgcann = list()
            cgcann.extend(add_cgcann(gene_a, cgc_ann)) 
            cgcann.extend(add_cgcann(gene_b, cgc_ann))

            if capture == "WG":
                curator = "NO"
                if list(filter(None, cgcann)) != []:
                    curator = "YES"
            else:
                curator = "NO"
                if svtype == 'TRA': 
                    if check_targets(chrom_a, start_a, end_a, targets) or \
                        check_targets(chrom_b, start_b, end_b, targets):
                        curator = "YES"
                else:
                    if check_targets(chrom_a, start_a, end_b, targets):
                        curator = "YES"

            if tool == 'gridss' and svtype == "INS":
                alt_seq = ''.join(list(filter(str.isalpha, alt)))
                svlength = len(alt_seq)
    
            gene_a_b = [gene_a, gene_b]
            gene_a_b.sort()
            gene_a_b_sorted = ",".join(gene_a_b)
            
            ## annotate cancer specific curation information
            project = os.path.basename(output).split('-')[0]
            curate = ''
            if curation_ann is not None:
                curate = annotate_curate_info(curation_ann, gene_a_b_sorted, sample, project)

            summary_sv.append([chrom_a, start_a, end_a, chrom_b, start_b, end_b, igv_coord, svtype,
                                svlength, sup_reads, vaf, tool, sdid, sample, gene_a, 
                                gene_b, gene_a_b_sorted, ",".join(list(filter(None, set(cgcann)))), curator, curate])
        
        svs_df = pd.DataFrame(summary_sv, columns = summary_columns)
        # hard filter for gridss germline svs
        gf_idx = svs_df[(svs_df['TOOL'] == "gridss") & \
                        (svs_df['SAMPLE'] == "germline") & \
                        (svs_df['SUPPORT_READS'].astype('int') < 40)].index
        svs_df.loc[list(gf_idx), "CURATOR"] = "NO"

        # checking exons overlaps
        svs_df['idx'] = svs_df.index
        # fetching indices except TRA and off target svs
        t_idx = set(svs_df[(svs_df['SVTYPE'] != "TRA") & (svs_df['CURATOR'] == "YES")].index)
        
        # generating pyranges for gridss and svcaller
        svs_gridss_df = svs_df[["CHROM_A", "START_A", "END_A", 
                                "CHROM_B", "START_B", "END_B", 
                                "SVTYPE", "idx", "CURATOR", 
                                "TOOL"]].rename(columns = {
                                    "CHROM_A": "Chromosome", 
                                    "START_A": "Start", 
                                    "END_B": "End"
                                })
        svs_gridss_pr = pr.PyRanges(svs_gridss_df.loc[(svs_gridss_df['SVTYPE'] != "TRA") & \
                                                      (svs_gridss_df['CURATOR'] == "YES") & \
                                                      (svs_gridss_df['TOOL'] == "gridss") ])
        svs_svcaller_df = svs_df[["CHROM_A", "START_A", "END_A", 
                                  "CHROM_B", "START_B", "END_B", 
                                  "SVTYPE", "idx", "CURATOR", 
                                  "TOOL"]].rename(columns = {
                                      "CHROM_A": "Chromosome", 
                                      "START_A": "Start", 
                                      "END_B": "End"})
        svs_svcaller_pr = pr.PyRanges(svs_svcaller_df.loc[(svs_svcaller_df['SVTYPE'] != "TRA") & \
                                                          (svs_svcaller_df['CURATOR'] == "YES") & \
                                                          (svs_svcaller_df['TOOL'] == "svcaller")])

        # find exons overlapping svs
        hits_gridss = svs_gridss_pr.intersect(exons)
        hits_idx = set(hits_gridss.idx) if hits_gridss else set()
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
    parser.add_argument('-c', '--curation-ann', help="Cancer gene specific curation annotation")
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
    
    curation_ann = pd.read_csv(args.curation_ann) \
                        if args.curation_ann else None

    if sv_caller == 'svcaller':
        parse_svcaller(input_file, sdid, output, vcftype)
    elif sv_caller == 'gridss':
        parse_gridss(input_file, sdid, output, vcftype)

    if annotBed:
        combined_input = combine_mut(input_file, output_dir)
        genes = pr.read_bed(annotBed)
        fh = open(target_json, 'r')
        targets = json.load(fh)
        if capture_kit in targets:
            targets = targets[capture_kit]
        annotate_combined_sv(combined_input, genes, targets, 
                             capture_kit, cgc_ann, output, exons, curation_ann)

