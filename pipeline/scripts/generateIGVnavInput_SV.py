#!/usr/bin/env python

import subprocess
import argparse
import os
import glob
import re
import json
import pandas as pd
import vcf


CGC_genes = {'CXCR4', 'FOXO1', 'PDE4D', 'EXT1', 'PWWP2A', 'NUTM2B', 'ERC1', 'MLLT1', 'ARID1A', 'IDH2', 'HOXD13', 'CD274', 'ZMYM3', 'ERRFI1', 'MUTYH', 'PRRX1', 'CALR', 'SKP2', 'CREB3L1', 'SMARCA1', 'FLNA', 'MYD88', 'STAG2', 'ANK1', 'EWSR1', 'MGA', 'CDK12', 'CCND1', 'BTG2', 'MAP2K4', 'CCNB1IP1', 'RMI2', 'GPC5', 'ACSL3', 'BCL10', 'LRP1B', 'CCND3', 'BAZ1A', 'SMARCB1', 'PPP4R2', 'BRAF', 'GOLGA5', 'PDGFRB', 'NDRG1', 'FGFR2', 'ARID1B', 'RGPD3', 'ZNF429', 'MLLT10', 'DPYD', 'PMS2', 'CASP3', 'TNFAIP3', 'DNM2', 'HRAS', 'KDM6A', 'PPFIBP1', 'LEF1', 'PRSS1', 'SMO', 'FCRL4', 'GATA1', 'LATS1', 'AFF4', 'MYH11', 'PDE4DIP', 'IL6ST', 'MEF2B', 'TCL1B', 'TMSB4X', 'DDB2', 'HOXB13', 'COX6C', 'AR', 'PIK3CG', 'TENT5C', 'NT5C2', 'A1CF', 'COMT', 'RPL22', 'SKI', 'SS18L1', 'EPHB1', 'JAK3', 'GPS2', 'PDCD1LG2', 'MST1', 'PAX7', 'OLIG2', 'N4BP2', 'ERCC3', 'BCL2', 'TPR', 'GNA11', 'RPL5', 'COL5A1', 'PRDM16', 'SMAD4', 'LCK', 'TSC1', 'CHCHD7', 'ABL1', 'PER1', 'SND1', 'WWTR1', 'FAT3', 'CREB3L2', 'NOTCH1', 'CLP1', 'ERBB3', 'TFE3', 'EPS15', 'FANCE', 'KLK2', 'FUS', 'SIX2', 'MGMT', 'AKT3', 'LEPROTL1', 'RAD17', 'CDK4', 'PSIP1', 'CBFB', 'CLTC', 'SLC34A2', 'IGH', 'ZFHX3', 'NTRK1', 'HOXC13', 'KEAP1', 'NUP214', 'SDHD', 'GTF2I', 'SOX2', 'APC', 'PDPK1', 'HLA-B', 'TRAF7', 'UGT1A1', 'APOBEC3B', 'TRAF3', 'CSF3R', 'HIF1A', 'EPHA7', 'MAFB', 'MAP2K2', 'ZEB1', 'ABI1', 'LRIG3', 'AXL', 'USP9X', 'MRTFA', 'FEV', 'POLG', 'SUZ12', 'TFRC', 'FOXR1', 'EIF4A2', 'BCL11B', 'WDCP', 'PAX5', 'WAS', 'BRCA2', 'OMD', 'VAV1', 'SSX2', 'MUC4', 'FAT1', 'PRKACA', 'PHOX2B', 'CDH11', 'KAT6A', 'PATZ1', 'CSMD1', 'SPRED1', 'CCR4', 'FGFR3', 'SETBP1', 'RELN', 'TCF3', 'DNMT3A', 'ARHGEF12', 'DDX10', 'PPP6C', 'SFRP4', 'ERBB4', 'DDIT3', 'CD209', 'CASP8', 'RALGDS', 'CIC', 'CSMD3', 'DDR2', 'PLCG1', 'PDCD1', 'CYP2D6', 'RBM10', 'PTEN', 'HSP90AA1', 'BCL2L12', 'RPL10', 'PIK3CB', 'NCOR2', 'TRB', 'KLF4', 'PTPN13', 'TRIP11', 'XPA', 'CBLC', 'ERCC2', 'AMER1', 'ATR', 'FIP1L1', 'INPPL1', 'GATA3', 'LYL1', 'NFKB2', 'FLI1', 'AFF3', 'USP8', 'PDGFRA', 'CNTNAP2', 'RAD51C', 'ETNK1', 'CYP2C8', 'KCNJ5', 'EZR', 'BTG1', 'ATM', 'STAT6', 'PRKD1', 'SMC1A', 'HERPUD1', 'HLF', 'SOX9', 'EED', 'PAFAH1B2', 'FLCN', 'CDH17', 'PTPRK', 'AKAP9', 'FANCD2', 'MDC1', 'TLX3', 'TPM4', 'RECQL4', 'EPHA3', 'KDM5C', 'COL1A1', 'VTI1A', 'MALT1', 'PTPRB', 'DCC', 'KAT7', 'MALAT1', 'EP300', 'DDX3X', 'CBLB', 'ERCC5', 'SSX1', 'MSH6', 'PPP2R1A', 'NBN', 'RUNX1T1', 'PIK3R1', 'KEL', 'POLQ', 'IRF4', 'RICTOR', 'CREB1', 'MRE11', 'MSH2', 'ELL', 'RBM15', 'FANCL', 'SMAD3', 'PTPN11', 'JAZF1', 'TBX3', 'BCL9L', 'NIN', 'S100A7', 'RAP1GDS1', 'CREBBP', 'TNFRSF14', 'GNAQ', 'NFATC2', 'IKZF1', 'STAG1', 'DROSHA', 'CEP89', 'PIK3CA', 'MSI2', 'SGK1', 'TFEB', 'PAX8', 'TFPT', 'PTPRD', 'LSM14A', 'SDC4', 'BCL2L11', 'MSN', 'MLF1', 'ROS1', 'ELN', 'CRTC1', 'MLH3', 'TCL1A', 'FAM135B', 'ELK4', 'NUMA1', 'PRPF40B', 'GLI1', 'TAL1', 'HGF', 'NUP98', 'NUTM2D', 'USP6', 'GNAS', 'CBL', 'PCBP1', 'WWOX', 'XPC', 'IGL', 'ARHGAP26', 'BARD1', 'EZH2', 'ETV4', 'PLCB4', 'PICALM', 'BUB1B', 'PRCC', 'CUL3', 'DUSP22', 'CTNNA1', 'HIST1H3B', 'SRC', 'CHEK1', 'TRIM27', 'NUTM1', 'CRTC3', 'MACC1', 'SMARCE1', 'ALB', 'SNX29', 'AFF1', 'ATRX', 'KMT2A', 'CHD1', 'RPTOR', 'MYC', 'BCL3', 'GNA13', 'FANCM', 'AKT1', 'TP63', 'SETD1B', 'COL3A1', 'SOS1', 'INTS4', 'STK11', 'IGK', 'NCKIPSD', 'RHOH', 'SDHAF2', 'MTCP1', 'IRS4', 'HMGN2P46', 'AXIN2', 'FHIT', 'BRD4', 'ATP2B3', 'SEPT6', 'ARHGEF10', 'CDK6', 'BCL9', 'GMPS', 'PRKCB', 'CD28', 'MN1', 'KNSTRN', 'BCL7A', 'NOTCH2', 'SRSF2', 'CAMTA1', 'Other Biomarkers', 'PRDM14', 'TLX1', 'SPTA1', 'NF1', 'PBX1', 'MB21D2', 'CCR7', 'CDKN2C', 'TSHR', 'RAD21', 'SIRPA', 'PMS1', 'PLK2', 'CDH10', 'SHTN1', 'AFDN', 'TEC', 'SET', 'TPMT', 'BAP1', 'FOXO4', 'NCOA4', 'STRN', 'YWHAE', 'KMT2B', 'PIK3CD', 'PRKAR1A', 'SBDS', 'SETD2', 'ZRSR2', 'EIF3E', 'WRN', 'THRAP3', 'RAD51D', 'DICER1', 'SOX21', 'QKI', 'MDM2', 'STK19', 'HOXA9', 'RHOA', 'RAD50', 'BIRC6', 'BMP5', 'TRA', 'NPM1', 'ZNF479', 'CCND2', 'HMGA1', 'ZNF384', 'RXRA', 'SFPQ', 'DEK', 'SH2B3', 'TYK2', 'MDS2', 'ELF4', 'IRS2', 'TGFBR1', 'GRM3', 'SPEN', 'CARS', 'CD70', 'JAK2', 'ETV5', 'NFKBIE', 'NCOR1', 'NKX2-1', 'EML4', 'ERBB2', 'NCOA1', 'XPO1', 'SH3GL1', 'AXIN1', 'BIRC3', 'MLLT11', 'RSPO2', 'KLF6', 'MYB', 'NF2', 'MAP2K1', 'ACVR1B', 'KNL1', 'SESN2', 'GREM1', 'TAF15', 'HNRNPA2B1', 'CDC27', 'H3F3A', 'ANKRD11', 'IL7R', 'BCR', 'IKBKB', 'NSD2', 'PTPRC', 'CCDC6', 'SOX17', 'U2AF1', 'RNF213', 'MNX1', 'IL3', 'IGF2BP2', 'TNFRSF17', 'MAF', 'LZTR1', 'RASA1', 'SEPT9', 'POU5F1', 'NBEA', 'ARHGEF10L', 'BCL11A', 'ACVR1', 'ESR1', 'ELF3', 'FANCA', 'MAML2', 'SOCS1', 'IDH1', 'KDSR', 'VHL', 'CIITA', 'NR4A3', 'SMARCD1', 'MCL1', 'MUC1', 'RFWD3', 'PCM1', 'POU2AF1', 'EBF1', 'DNMT3B', 'FANCF', 'FANCG', 'STAT3', 'GPC3', 'RHEB', 'ERG', 'TBL1XR1', 'SDHA', 'SDHC', 'KAT6B', 'KTN1', 'NTRK3', 'KIF5B', 'STIL', 'PBRM1', 'CDC73', 'TRRAP', 'NCOA2', 'B2M', 'LMO2', 'ARNT', 'CDKN1A', 'CDKN2B', 'TSC2', 'CBFA2T3', 'TERT', 'ZCCHC8', 'P2RY8', 'FAS', 'IKZF3', 'FGFR1OP', 'HOXD11', 'NRG1', 'EPCAM', 'SDHB', 'DIS3', 'TMPRSS2', 'ABL2', 'JAK1', 'ALDH2', 'TCF12', 'TFG', 'CTNND2', 'MDM4', 'ARID2', 'PRDM1', 'GRIN2A', 'BMPR1A', 'WT1', 'RNF43', 'RAC1', 'HSP90AB1', 'ISX', 'GOPC', 'PTCH1', 'STAT5B', 'ZMYM2', 'LPP', 'SRSF3', 'PTK6', 'CDX2', 'UBR5', 'BRIP1', 'ATF1', 'KMT2C', 'TAF1', 'GPHN', 'RAF1', 'TLR4', 'MLLT3', 'GATA2', 'SPECC1', 'NSD1', 'ITGAV', 'CTNNB1', 'MLLT6', 'SRGAP3', 'DDX6', 'CD74', 'RIT1', 'RPS6KA4', 'FAT4', 'CDH1', 'REL', 'H3F3B', 'BCL6', 'PHF6', 'DNAJB1', 'CD79B', 'CDKN2A', 'CYSLTR2', 'NKX3-1', 'RB1', 'ETV1', 'ID3', 'IRF2', 'MLH1', 'SETDB1', 'PTPRT', 'TGFBR2', 'PRDM2', 'HOOK3', 'NACA', 'FAM47C', 'RANBP2', 'ASXL2', 'WIF1', 'BAX', 'LASP1', 'MAX', 'PTPN6', 'IL2', 'HIP1', 'CNTRL', 'ROBO2', 'MYCN', 'POLE', 'CNBD1', 'CRNKL1', 'SALL4', 'RABEP1', 'TCF7L2', 'ARHGAP35', 'CYLD', 'KDR', 'DCAF12L2', 'CTCF', 'TP53', 'RSPO3', 'RARA', 'SIX1', 'KDM5A', 'ARHGAP5', 'PGR', 'PIM1', 'BRCA1', 'ASXL1', 'DDX5', 'ACVR2A', 'CNOT3', 'SS18', 'CCNE1', 'FH', 'CHEK2', 'MUC16', 'EPOR', 'LYN', 'CASP9', 'TRIM33', 'FOXA1', 'FLT3', 'MSH3', 'GABRA6', 'FOXL2', 'BCORL1', 'ALK', 'TPM3', 'FES', 'SEPT5', 'MYCL', 'RUNX1', 'CHD4', 'PREX2', 'DUX4L1', 'CRLF2', 'FSTL3', 'CTNNA2', 'BTK', 'NTHL1', 'HSD3B1', 'CHD3', 'FUBP1', 'FEN1', 'ATIC', 'RAD51', 'PABPC1', 'MPL', 'RAD51B', 'CCNC', 'NTRK2', 'LMNA', 'CHST11', 'EIF1AX', 'NUP93', 'SMC3', 'ZNF521', 'SLC45A3', 'ARID5B', 'TET1', 'ZNF331', 'NQO1', 'WNK2', 'LMO1', 'MTOR', 'TAL2', 'LATS2', 'FGFR1', 'PIK3R2', 'CDKN1B', 'TNC', 'RPN1', 'TOP1', 'FOXP1', 'NRAS', 'GSK3B', 'DCTN1', 'FNBP1', 'NSD3', 'RAD54L', 'CHD2', 'PAX3', 'FBXW7', 'BLM', 'CSDE1', 'LIFR', 'BCOR', 'PARP1', 'SF3B1', 'IL21R', 'MEN1', 'TCEA1', 'POLD1', 'PRF1', 'FBLN2', 'PDGFB', 'COL2A1', 'DGCR8', 'SYK', 'DAXX', 'CANT1', 'CTLA4', 'CSF1R', 'MITF', 'SSX4', 'CUX1', 'PML', 'CPEB3', 'CHIC2', 'POLR2A', 'EGFR', 'MET', 'KIT', 'USP44', 'MAP3K1', 'FGFR4', 'TET2', 'IGF1R', 'NCOA3', 'HOXC11', 'GAS7', 'CARD11', 'MYO5A', 'HMGA2', 'FBXO11', 'NAB2', 'RGS7', 'JUN', 'EXT2', 'CLTCL1', 'TRD', 'MYH9', 'ACSL6', 'FKBP9', 'HNF1A', 'HOXA11', 'AKT2', 'CEBPA', 'TMEM127', 'RET', 'HOXA13', 'CD79A', 'FCGR2B', 'DNMT1', 'PPM1D', 'MECOM', 'ECT2L', 'ERCC4', 'SMARCA4', 'KRAS', 'HLA-A', 'CLIP1', 'RRAS2', 'PRKCI', 'ASPSCR1', 'SUFU', 'TRIM24', 'ZBTB16', 'ARAF', 'NONO', 'ZNRF3', 'CNBP', 'LHFPL6', 'ACKR3', 'HEY1', 'CTNND1', 'RAD52', 'BCLAF1', 'MAPK1', 'HLA-C', 'MED12', 'MYOD1', 'FAM131B', 'ATP1A1', 'RECQL', 'PLAG1', 'EPAS1', 'LARP4B', 'MAP3K13', 'PALB2', 'SMAD2', 'SPOP', 'ITK', 'PTPRS', 'POT1', 'NFE2L2', 'KMT2D', 'CACNA1D', 'ETV6', 'LCP1', 'FANCC', 'C15orf65', 'FOXO3', 'BRD3', 'NFIB', 'KIAA1549', 'FLT4', 'PPARG', 'KLF2'}


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

    vcf_reader = vcf.Reader(open(input_vcf, 'r'))
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


def annotate_combined_sv(combined_file, genes, targets, capture, varann, output):
    """
    Parsing combined sv list and apply gene annotation for each SV
    """
    # output_file = open(output, 'w')
    summary_columns = ['CHROM_A', 'START_A', 'END_A', 'CHROM_B', 'START_B', 'END_B',
                       'IGV_COORD', 'SVTYPE', 'SV_LENGTH', 'SUPPORT_READS', 'TOOL', 'SDID', 'SAMPLE',
                       'GENE_A', 'GENE_B', "GENE_A-GENE_B-sorted", "SOURCES", "CURATOR"]
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

            if ':' in alt:
                chrom_b = ''.join(list(filter(str.isdigit, alt.split(':')[0])))
                start_b = ''.join(list(filter(str.isdigit, alt.split(':')[1])))
                end_b = int(start_b) + 1

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
            
            if svtype != 'TRA' and not tool == 'svcaller':
                chrom_b = 'NA'
                start_b = 'NA'
                end_b = 'NA'

            igv_coord = ' '.join([igv_coord_a, igv_coord_b])
            gene_a = gene_annotation(chrom_a, start_a, end_a, genes)

            if chrom_b != 'NA':
                svlength = abs(int(end_b)-int(start_a)) if chrom_a == chrom_b else 'NA'
                gene_b = gene_annotation(chrom_b, start_b, end_b, genes)
            else:
                gene_b = 'NA'
            
            sources = ''
            if capture == "WG":
                curator = "NO"
                if gene_a in CGC_genes or gene_b in CGC_genes:
                    sources = 'CGC'

                if sources != '':
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
                                gene_b, gene_a_b_sorted, sources, curator])
        summary_sv_df = pd.DataFrame(summary_sv, columns = summary_columns)
        summary_sv_df_sorted = summary_sv_df.sort_values(["GENE_A-GENE_B-sorted", "CHROM_A", "START_A", "CHROM_B", "START_B", "TOOL"], 
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
    parser.add_argument('--varann', help="variant annotation from multiple resources like CGC, Oncokb")
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
 
    if args.target:
        capture_kit, target_json = args.target

    varann = ''
    if args.varann:
        varann = json.load(open(args.varann, 'r'))

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
        annotate_combined_sv(combined_input, genes, targets, capture_kit, varann, output)

