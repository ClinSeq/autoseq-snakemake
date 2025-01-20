#!/usr/bin/env Rscript

#### Dependencies and arguments ###########################
library(ggplot2)
library(plotly)
library(getopt)
library(RJSONIO)
library(VariantAnnotation)

library(data.table)
library(GenomicRanges)
library(stringr)

#long, short(NA), argmask, datatype, desc
#argmask 0=no arg, 1=req, 2=optional
args <- rbind(
    c("tumor_cnr", NA, 1, "character", "tumor bin file from CNVkit"),
    c("tumor_cns", NA, 1, "character", "tumor segment file from CNVkit"),
    c("normal_cnr", NA, 1, "character", "normal bin file from CNVkit"),
    c("normal_cns", NA, 1, "character", "normal segment file from CNVkit"),
    c("het_snps_vcf", NA, 1, "character", "heterozygous SNPs .vcf file"),
    c("purecn_csv", NA, 1, "character", "PureCN result .csv file"),
    c("purecn_genes_csv", NA, 1, "character", "PureCN result _genes.csv"),
    c("purecn_loh_csv", NA, 1, "character", "PureCN result _loh.csv"),
    c("purecn_variants_csv", NA, 1, "character", "PureCN result _variants.csv"),
    c("svcaller_T_DEL", NA, 1, "character", "Tumor SV caller DEL-events.gtf"),
    c("svcaller_T_DUP", NA, 1, "character", "Tumor SV caller DUP-events.gtf"),
    c("svcaller_T_INV", NA, 1, "character", "Tumor SV caller INV-events.gtf"),
    c("svcaller_T_TRA", NA, 1, "character", "Tumor SV caller TRA-events.gtf"),
    c("svcaller_N_DEL", NA, 1, "character", "Normal SV caller DEL-events.gtf"),
    c("svcaller_N_DUP", NA, 1, "character", "Normal SV caller DUP-events.gtf"),
    c("svcaller_N_INV", NA, 1, "character", "Normal SV caller INV-events.gtf"),
    c("svcaller_N_TRA", NA, 1, "character", "Normal SV caller TRA-events.gtf"),
    c("germline_mut_vcf", NA, 1, "character", "germline mutation vcf file"),
    c("somatic_mut_vcf", NA, 1, "character", "somatic mutation vcf file"),
    c("plot_png", NA, 1, "character", "plot .png file name"),
    c("plot_png_normal", NA, 1, "character", "normal CNV plot .png file name"),
    c("cna_json", NA, 1, "character", "CNA output json file name"),
    c("purity_json", NA, 1, "character", "purity output json file name"),
    c("gene_track", NA, 1, "character", "gene_track csv file name")
)


opts <- getopt(args)


### Chromosome size ####
chrsizes=structure(list(
    chr = c("1", "2", "3", "4", "5", "6", "7", "8",
            "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
            "20", "21", "22", "X", "Y", "MT"),
    size = c(249250621L, 243199373L,
             198022430L, 191154276L, 180915260L, 171115067L, 159138663L, 146364022L,
             141213431L, 135534747L, 135006516L, 133851895L, 115169878L, 107349540L,
             102531392L, 90354753L, 81195210L, 78077248L, 59128983L, 63025520L,
             48129895L, 51304566L, 155270560L, 59373566L, 16569L),
    cumstart = c(0,
                 249250621, 492449994, 690472424, 881626700, 1062541960, 1233657027,
                 1392795690, 1539159712, 1680373143, 1815907890, 1950914406, 2084766301,
                 2199936179, 2307285719, 2409817111, 2500171864, 2581367074, 2659444322,
                 2718573305, 2781598825, 2829728720, 2881033286, 3036303846, 3095677412
    )),
    .Names = c("chr", "size", "cumstart"),
    row.names = c(NA,-25L),
    class = "data.frame")


## Genes to plot  #######
genes=data.frame(label = c("AKT1", "APC", "AR", "ARID1A", "ARID2", "ATM", "ATR", "BARD1",
                           "BRAF", "BRCA1", "BRCA2", "BRIP1", "CCND1", "CDH1", "CDK12", "CDK4",
                           "CDK6", "CDKN1A", "CDKN1B", "CDKN2A", "CHD1", "CHEK2", "CTNNB1",
                           "CUL3", "DICER1", "DNMT3A", "FANCA", "FOXA1", "FOXO1", "HRAS", "IDH1",
                           "JAK1", "KDM6A", "KEAP1", "KMT2A", "KMT2C", "KMT2D", "KRAS", "MED12",
                           "MET", "MGA", "MLH1", "MLH3", "MRE11A", "MSH2", "MSH3", "MSH6", "MYC",
                           "NBN", "NCOR1", "NKX3-1", "NRAS", "PALB2", "PIK3CA", "PIK3CB", "PIK3CD",
                           "PIK3R1", "PIK3R2", "PMS1", "PMS2", "POLD1", "POLE", "PTEN", "RAD50",
                           "RAD51", "RAD51B", "RAD51C", "RAD51D", "RB1", "RNF43", "SETD2", "SF3B1",
                           "SPEN", "SPOP", "TMPRSS2", "TP53", "U2AF1", "XPO1", "ZBTB16", "ZFHX3",
                           "ZMYM3"),
                 chromosome = c("14", "5", "X", "1", "12", "11", "3", "2", "7", "17", "13", "17",
                                "11", "16", "17", "12", "7", "6", "12", "9", "5", "22", "3", "2",
                                "14", "2", "16", "14", "13", "11", "2", "1", "X", "19", "11", "7",
                                "12", "12", "X", "7", "15", "3", "14", "11", "2", "5", "2", "8",
                                "8", "17", "8", "1", "16", "3", "3", "1", "5", "19", "2", "7", "19",
                                "12", "10", "5", "15", "14", "17", "17", "13", "17", "3", "2", "1",
                                "17", "21", "17", "21", "2", "11", "16", "X"),
                 start = c(105235686, 112043195, 66764465, 27022524, 46123448, 108093211, 142168077,
                           215590370, 140419127, 41196312, 32889611, 59758627, 69455855, 68771128,
                           37617764, 58141510, 92234235, 36644305, 12867992, 21967751, 98190908,
                           29083731, 41236328, 225334867, 95552565, 25455845, 89803957, 38059189,
                           41129804, 532242, 209100951, 65298912, 44732757, 10596796, 118307205,
                           151832010, 49412758, 25357723, 70338406, 116312444, 41913422, 37034823,
                           75480467, 94152895, 47630108, 79950467, 47922669, 128747680, 90945564,
                           15932471, 23536206, 115247090, 23614488, 178865902, 138372860, 9711790,
                           67511548, 18263928, 190649107, 6012870, 50887461, 133200348, 89622870,
                           131891711, 40986972, 68286496, 56769934, 33426811, 48877887, 56429861,
                           47057919, 198254508, 16174359, 47676246, 42836478, 7565097, 44513066,
                           61704984, 113930315, 72816784, 70459474),
                 end = c(105262088, 112181936, 66950461, 27108595, 46301823, 108239829, 142297668,
                         215674428, 140624564, 41277500, 32973805, 59940882, 69469242, 68869451,
                         37721160, 58149796, 92465908, 36655116, 12875305, 21995300, 98262240,
                         29138410, 41301587, 225450110, 95624347, 25565459, 89883065, 38069245,
                         41240734, 537287, 209130798, 65432187, 44971847, 10614417, 118397539,
                         152133090, 49453557, 25403870, 70362303, 116438440, 42062141, 37107380,
                         75518235, 94227074, 47789450, 80172279, 48037240, 128753674, 91015456,
                         16121499, 23540440, 115259515, 23652631, 178957881, 138553780, 9789172,
                         67597649, 18281350, 190742355, 6048756, 50921273, 133263951, 89731687,
                         131980313, 41024354, 69196935, 56811703, 33448541, 49056122, 56494956,
                         47205457, 198299815, 16266955, 47755596, 42903043, 7590856, 44527697,
                         61765761, 114121398, 73093597, 70474996),
                 stringsAsFactors = F)
genes$cumstart <- genes$start + chrsizes$cumstart[match(genes$chromosome,chrsizes$chr)]
genes$cumend <- genes$end + chrsizes$cumstart[match(genes$chromosome,chrsizes$chr)]

### Read gene track ######################
data2 <- read.table(opts$gene_track, sep=",", na.strings=".", stringsAsFactors=FALSE, quote="", fill=FALSE , comment.char="#", header=TRUE)
data2$cumstart <- data2$Start + chrsizes$cumstart[match(data2$Chromosome,chrsizes$chr)]
data2$cumend <- data2$End + chrsizes$cumstart[match(data2$Chromosome,chrsizes$chr)]



data2_mRNA <- data2[data2$Feature == 'transcript',]
data2_3p <- data2[data2$Feature == 'three_prime_utr',]
data2_5p <- data2[data2$Feature == 'five_prime_utr',]
data2_exon <- data2[data2$Feature == 'exon',]




#save.image('ws.Rdata')
#stop('as desired')

### Read SNP allele ratio  #######

{
    vcf <- readVcf(opts$het_snps_vcf,genome = "GRCh37")
    vcf <- vcf[seqnames(vcf) %in% c(1:22,'X','Y')]
    
    g <- geno(vcf)
    r <- rowRanges(vcf)
    chr <- as.character(seqnames(r))
    pos <- data.frame(ranges(r))$start
    alf=NULL
    smoothedAi=NA
    if (length(pos)>0) {
        pos=as.numeric(pos)
        alf <- data.table(chromosome=chr, start=pos, end=pos, stringsAsFactors = F, cumstart=NA, cumend=NA)
    }
    
    
    
    if (!is.null(alf)) if( ! all(is.na(alf$chromosome)) ) {
        
        for(chr in chrsizes$chr){
            ix <- which(alf$chromosome == chr)
            alf$cumstart[ix] <- alf$start[ix] + chrsizes$cumstart[chrsizes$chr==chr]
            alf$cumend[ix] <- alf$end[ix] + chrsizes$cumstart[chrsizes$chr==chr]
        }
        
        alf$t <- as.numeric(sapply(g$AD[,2], "[", 2))/as.numeric(g$DP[,2])
        alf$t[is.nan(alf$t)]=NA # allele freq becomes NaN if cov=0. Then set to NA
        alf$n <- as.numeric(sapply(g$AD[,1], "[", 2))/as.numeric(g$DP[,1])
        alf$n[is.nan(alf$n)]=NA # allele freq becomes NaN if cov=0. Then set to NA
        alf$td <- as.numeric(g$DP[,2])
        alf$nd <- as.numeric(g$DP[,1])
        
        alf$ai=2*abs(alf$t-0.5)
        alf$ai_n=2*abs(alf$n-0.5)
        
        # require presence and heteroz in normal
        alf <- alf[nd > 15 & td > 50 & n > 0.3 & n < 0.7]
        
    }
}

#save.image('ws.1.Rdata')

### Read somatic point mutations #####
{
    vcf <- readVcf(opts$somatic_mut_vcf,genome = "GRCh37")
    vcf <- vcf[rowRanges(vcf)$FILTER == 'PASS']
    
    
    g <- geno(vcf)
    r=rowRanges(vcf)
    if (length(g)>0) { # if there are any somatic mutations...
        chr=as.character(seqnames(r))
        pos=data.frame(ranges(r))$start
        salf <- data.table(N=1:length(chr),chromosome=chr,pos=pos,stringsAsFactors = F)
        rownames(salf)=names(r)
        salf$REF=as.data.frame(r$REF)[,1]
        salf$ALT=as.data.frame(r$ALT)[,3]
        salf$cumpos <- NA
        for(chr in chrsizes$chr){
            ix <- which(salf$chromosome == chr)
            salf$cumpos[ix] <- salf$pos[ix] + chrsizes$cumstart[chrsizes$chr==chr]
        }
        
        salf$FILTER <- rowRanges(vcf)$FILTER
        
        salf$AF.T <- as.numeric(g$VAF[,2])
        salf$AO.T <- as.numeric(apply(g$DP4[,2,3:4],1,sum))  # sum alt forward and alt reverse
        salf$DP.T <- as.numeric(apply(g$DP4[,2,],1,sum))  # sum ref forward, ref reverse, alt forward and alt reverse
        salf$AO.N <- as.numeric(apply(g$DP4[,1,3:4],1,sum))  # sum alt forward and alt reverse
        salf$DP.N <- as.numeric(apply(g$DP4[,1,],1,sum))  # sum ref forward, ref reverse, alt forward and alt reverse
        
        salf$type='other'
        salf$type[isSNV(vcf)]='snv'
        salf$type[isDeletion(vcf)]='del'
        salf$type[isInsertion(vcf)]='ins'
        
        
        header=info(header(vcf))$Description
        ix=grep('Consequence annotations from Ensembl',header)
        header=strsplit(header[ix],'\\|')[[1]]
        header[1]='Allele'
        
        vep=info(vcf)$CSQ
        
        ## This loop creates a new "table" with all mutation effects.
        rowspermut=unlist(lapply(vep,length))
        table=#data.frame(
            matrix('',nrow = sum(rowspermut),ncol = length(header)+1)#,stringsAsFactors = F)
        colnames(table)=c('N',header)  # blir detta fel ibland?????
        
        
        for (i in 1:length(vep)) {
            #if (i %% 100 ==0) cat(i,'..')
            for (j in 1:length(vep[[i]])) { # for each effect
                t2=vep[[i]][j]
                t2=strsplit(t2,'[|]')[[1]] # pipe separated line with one effect of the mutation
                t2=c((i),t2)
                thisrow=sum(rowspermut[1:i])-rowspermut[i]+j
                table[thisrow,1:length(t2)]=t2
            }
        }
        
        # temp <- lapply(vep,function(vector) {
        #     as.data.table(t(strsplit(vector,'[|]')[[1]]))
        # })
        # table <- rbindlist(temp,fill=T)
        # table[,N:=1:.N]
        # colnames(table)=c(header,'N')
        
        table <- as.data.table(table)
        table[,N:=as.integer(N)]
        salf=merge(salf,table,by='N',all=T)
        
    } #end somatic mutations
    #salf=(salf[,-1])
    
    # mark the type
    salf$pch=rep(22,nrow(salf))
    salf$pch[salf$type=='snv']=21
    salf$pch[salf$type=='del']=24
    salf$pch[salf$type=='ins']=25
    
    
    # Icke-NA/intron på konsekvens
    salf$hasConsequence=!salf$Consequence %in% c("intron_variant", "synonymous_variant",
                                                 "splice_region_variant&intron_variant",
                                                 "3_prime_UTR_variant", "intergenic_variant",
                                                 "regulatory_region_variant", "upstream_gene_variant",
                                                 "downstream_gene_variant",
                                                 "intron_variant&non_coding_transcript_variant",
                                                 "5_prime_UTR_variant", "splice_region_variant&synonymous_variant",
                                                 "non_coding_transcript_exon_variant&non_coding_transcript_variant",
                                                 "intron_variant&NMD_transcript_variant", "TF_binding_site_variant",
                                                 "splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant",
                                                 "mature_miRNA_variant", "coding_sequence_variant&5_prime_UTR_variant")
    salf$hasConsequence[is.na(salf$hasConsequence)]=F
    
    
    ## Vilka är tydligt deleterious? (inkl inframe)
    salf$deleterious=rep(F,nrow(salf))
    salf$deleterious[unique(c(
        grep('start',salf$Consequence),
        grep('stop',salf$Consequence),
        grep('frame',salf$Consequence),
        grep('acceptor',salf$Consequence),
        grep('donor',salf$Consequence)
    ))]=T
    salf$deleterious[!salf$hasConsequence]=F
    
    # remove the rsids (SNPs missed in normal)
    salf <- salf[!str_detect(Existing_variation,'^rs')]
    
    # remove too few obs
    salf <- salf[AO.T > 3]
}
#save.image('ws.2.Rdata')

### Read germline point mutations ####
{
    vcf <- readVcf(opts$germline_mut_vcf,genome = "GRCh37")
    vcf <- expand(vcf)
    
    
    vep=info(vcf)$CSQ
    
    # convert to data.table
    vep <- as.data.table(vep)
    # retain only HIGH and pathogenic
    vep <- vep[str_detect(value,'pathogenic') | str_detect(value,'HIGH')]
    
    # subset vcf
    vcf <- vcf[unique(vep$group)]
    
    g <- geno(vcf)
    r=rowRanges(vcf)
    f <- rowRanges(vcf)
    if (length(g)>0) { # if there are any somatic mutations...
        chr=as.character(seqnames(r))
        pos=data.frame(ranges(r))$start
        galf <- data.table(N=1:length(chr),chromosome=chr,pos=pos,stringsAsFactors = F)
        rownames(galf)=names(r)
        galf$cumpos <- NA
        galf$REF <- as.character(f$REF)
        galf$ALT <- as.character(f$ALT)
        for(chr in chrsizes$chr){
            ix <- which(galf$chromosome == chr)
            galf$cumpos[ix] <- galf$pos[ix] + chrsizes$cumstart[chrsizes$chr==chr]
        }
        
        galf$AF <- g$AD[,1,2]/g$DP[,1]
        #galf$AF_n <- g$AD[,2,2]/g$DP[,2]    # <----- this should be the normal sample allele ratio.
        galf$AO <- g$AD[,1,2]
        galf$DP <- g$DP[,1]
        
        galf$type='other'
        galf$type[isSNV(vcf)]='snv'
        galf$type[isDeletion(vcf)]='del'
        galf$type[isInsertion(vcf)]='ins'
        
        
        header=info(header(vcf))$Description
        ix=grep('Consequence annotations from Ensembl',header)
        header=strsplit(header[ix],'\\|')[[1]]
        header[1]='Allele'
        
        
        vep=info(vcf)$CSQ
        ## This loop creates a new "table" with all mutation effects.
        rowspermut=unlist(lapply(vep,length))
        table=#data.frame(
            matrix('',nrow = sum(rowspermut),ncol = length(header)+1)#,stringsAsFactors = F)
        colnames(table)=c('N',header)  # blir detta fel ibland?????
        
        for (i in 1:length(vep)) {
            for (j in 1:length(vep[[i]])) { # for each effect
                t2=vep[[i]][j]
                t2=strsplit(t2,'[|]')[[1]] # pipe separated line with one effect of the mutation
                t2=c(i,t2)
                thisrow=sum(rowspermut[1:i])-rowspermut[i]+j
                table[thisrow,1:length(t2)]=t2
            }
        }
        
        # temp <- lapply(vep[rowspermut==1],function(vector) {
        #     t <- as.data.table(t(strsplit(vector,'[|]')[[1]]))
        #     # colnames(t)=c(header,'N')
        # })
        # table <- rbindlist(temp,fill=T)
        # 
        # colnames(table)=header
        
        table <- as.data.table(table)
        table$AF <- NULL  # remove AF col from vep data to not confuse with AF from galf
        
        # Add alleles as vep names them, to use for merge
        galf$Allele = ifelse(galf$type=="snv", yes = galf$ALT, no = substr(galf$ALT, start=2, stop=nchar(galf$ALT)))
        galf$Allele[which(galf$Allele== "")] = "-"
        # Merge with vep table
        table[,N:=as.integer(N)]
        galf=merge(galf,table,by=c('N','Allele'),all.x=T)
    } #end germline mutations
    galf$gnomAD_AF=as.numeric(galf$gnomAD_AF)  # Make gnom_AD numerical so it can be used
    galf=(galf[which(galf$AO>=12 & galf$AF>=0.2 & (galf$gnomAD_AF < 0.01 | is.na(galf$gnomAD_AF))),-1])
    # mark the type
    galf$pch=rep(22,nrow(galf))
    galf$pch[galf$type=='snv']=21
    galf$pch[galf$type=='del']=24
    galf$pch[galf$type=='ins']=25
    
    # Icke-NA/intron på konsekvens
    galf$hasConsequence=!galf$Consequence %in% c("intron_variant", "synonymous_variant",
                                                 "splice_region_variant&intron_variant",
                                                 "3_prime_UTR_variant", "intergenic_variant",
                                                 "regulatory_region_variant", "upstream_gene_variant",
                                                 "downstream_gene_variant",
                                                 "intron_variant&non_coding_transcript_variant",
                                                 "5_prime_UTR_variant", "splice_region_variant&synonymous_variant",
                                                 "non_coding_transcript_exon_variant", "non_coding_transcript_variant",
                                                 "non_coding_transcript_exon_variant&non_coding_transcript_variant",
                                                 "intron_variant&NMD_transcript_variant", "TF_binding_site_variant",
                                                 "splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant",
                                                 "mature_miRNA_variant", "coding_sequence_variant&5_prime_UTR_variant")
    galf$hasConsequence[is.na(galf$hasConsequence)]=F
    
    
    ## Vilka är tydligt deleterious? (inkl inframe)
    galf$deleterious=rep(F,nrow(galf))
    galf$deleterious[unique(c(
        grep('start',galf$Consequence),
        grep('stop',galf$Consequence),
        grep('frame',galf$Consequence),
        grep('acceptor',galf$Consequence),
        grep('donor',galf$Consequence)
    ))]=T
    galf$deleterious[!galf$hasConsequence]=F
    
    
    
}
#save.image('ws.3.Rdata')

### Read CNVkit copy number data   #####
{ #tumor
    segments <- fread(opts$tumor_cns,stringsAsFactors = F)
    bins <- fread(opts$tumor_cnr,stringsAsFactors = F)
    
    ## Segment start and end pos
    segments$cumstart <- NA
    segments$cumend   <- NA
    for(chr in chrsizes$chr){
        ix <- which(segments$chromosome == chr)
        segments$cumstart[ix] <- segments$start[ix] + chrsizes$cumstart[chrsizes$chr==chr]
        segments$cumend[ix] <- segments$end[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    }
    segments$centerpos <- segments$cumstart+(segments$cumend-segments$cumstart)/2
    
    ## Bin start and end pos
    bins$cumstart <- NA
    bins$cumend   <- NA
    for(chr in chrsizes$chr){
        ix <- which(bins$chromosome == chr)
        bins$cumstart[ix] <- bins$start[ix] + chrsizes$cumstart[chrsizes$chr==chr]
        bins$cumend[ix] <- bins$end[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    }
    bins$centerpos <- bins$cumstart+(bins$cumend-bins$cumstart)/2
    bins=bins[order(bins$cumstart),]
}

{ # normal
    segments_n <- fread(opts$normal_cns,stringsAsFactors = F)
    bins_n <- fread(opts$normal_cnr,stringsAsFactors = F)
    
    ## Segment start and end pos
    segments_n$cumstart <- NA
    segments_n$cumend   <- NA
    for(chr in chrsizes$chr){
        ix <- which(segments_n$chromosome == chr)
        segments_n$cumstart[ix] <- segments_n$start[ix] + chrsizes$cumstart[chrsizes$chr==chr]
        segments_n$cumend[ix] <- segments_n$end[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    }
    segments_n$centerpos <- segments_n$cumstart+(segments_n$cumend-segments_n$cumstart)/2
    
    ## Bin start and end pos
    bins_n$cumstart <- NA
    bins_n$cumend   <- NA
    for(chr in chrsizes$chr){
        ix <- which(bins_n$chromosome == chr)
        bins_n$cumstart[ix] <- bins_n$start[ix] + chrsizes$cumstart[chrsizes$chr==chr]
        bins_n$cumend[ix] <- bins_n$end[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    }
    bins_n$centerpos <- bins_n$cumstart+(bins_n$cumend-bins_n$cumstart)/2
    bins_n=bins_n[order(bins_n$cumstart),]
}

### Read structural variant files
{ # tumor
    t_strvs=NULL
    try( {
        sv <- read.delim(opts$svcaller_T_DEL,header=F,stringsAsFactors = F)
        colnames(sv)[c(1,4,5)]=c('chr','start','end')
        sv$type='DEL'
        t_strvs=rbind(t_strvs,sv)
    }, silent=T)
    try( {
        sv <- read.delim(opts$svcaller_T_DUP,header=F,stringsAsFactors = F)
        colnames(sv)[c(1,4,5)]=c('chr','start','end')
        sv$type='DUP'
        t_strvs=rbind(t_strvs,sv)
    }, silent=T)
    try( {
        sv <- read.delim(opts$svcaller_T_INV,header=F,stringsAsFactors = F)
        colnames(sv)[c(1,4,5)]=c('chr','start','end')
        sv$type='INV'
        t_strvs=rbind(t_strvs,sv)
    }, silent=T)
    try( {
        sv <- read.delim(opts$svcaller_T_TRA,header=F,stringsAsFactors = F)
        colnames(sv)[c(1,4,5)]=c('chr','start','end')
        sv$type='TRA'
        t_strvs=rbind(t_strvs,sv)
    }, silent=T)
    t_strvs$cumstart=NA
    t_strvs$cumend=NA
    for(i in 1:nrow(chrsizes)){
        ix <- which(t_strvs$chr == chrsizes$chr[i])
        t_strvs$cumstart[ix] <- t_strvs$start[ix] + chrsizes$cumstart[i]
        t_strvs$cumend[ix] <- t_strvs$end[ix] + chrsizes$cumstart[i]
    }
}
{ # normal
    n_strvs=NULL
    try( {
        sv <- read.delim(opts$svcaller_N_DEL,header=F,stringsAsFactors = F)
        colnames(sv)[c(1,4,5)]=c('chr','start','end')
        sv$type='DEL'
        n_strvs=rbind(n_strvs,sv)
    }, silent=T)
    try( {
        sv <- read.delim(opts$svcaller_N_DUP,header=F,stringsAsFactors = F)
        colnames(sv)[c(1,4,5)]=c('chr','start','end')
        sv$type='DUP'
        n_strvs=rbind(n_strvs,sv)
    }, silent=T)
    try( {
        sv <- read.delim(opts$svcaller_N_INV,header=F,stringsAsFactors = F)
        colnames(sv)[c(1,4,5)]=c('chr','start','end')
        sv$type='INV'
        n_strvs=rbind(n_strvs,sv)
    }, silent=T)
    try( {
        sv <- read.delim(opts$svcaller_N_TRA,header=F,stringsAsFactors = F)
        colnames(sv)[c(1,4,5)]=c('chr','start','end')
        sv$type='TRA'
        n_strvs=rbind(n_strvs,sv)
    }, silent=T)
    n_strvs$cumstart=NA
    n_strvs$cumend=NA
    for(i in 1:nrow(chrsizes)){
        ix <- which(n_strvs$chr == chrsizes$chr[i])
        n_strvs$cumstart[ix] <- n_strvs$start[ix] + chrsizes$cumstart[i]
        n_strvs$cumend[ix] <- n_strvs$end[ix] + chrsizes$cumstart[i]
    }
}

#save.image('ws.4.Rdata')

### Read PureCN files #####
{
    # files to read (purity/ploidy, mutations and snps, gene copy number and LOH, segmented copy number and LOH)
    purecn_files = c(opts$purecn_csv, opts$purecn_variants_csv, opts$purecn_genes_csv, opts$purecn_loh_csv)
    # variables to assign to
    purecn_variables = c("purecn_stat", "purecn_vars", "purecn_genes", "purecn_loh")
    # colnames to use if no data avialable and mock df created
    purecn_colnames = list(
        c("Sampleid","Purity","Ploidy","Sex","Contamination","Flagged","Failed","Curated","Comment"),
        c("Sampleid", "chr", "start", "end", "ID", "REF", "ALT", "SOMATIC.M0", "SOMATIC.M1", "SOMATIC.M2",
          "SOMATIC.M3", "SOMATIC.M4", "SOMATIC.M5", "SOMATIC.M6", "SOMATIC.M7", "GERMLINE.M0", "GERMLINE.M1",
          "GERMLINE.M2", "GERMLINE.M3", "GERMLINE.M4", "GERMLINE.M5", "GERMLINE.M6", "GERMLINE.M7", "GERMLINE.CONTHIGH",
          "GERMLINE.CONTLOW", "GERMLINE.HOMOZYGOUS", "ML.SOMATIC", "POSTERIOR.SOMATIC", "ML.M", "ML.C", "ML.M.SEGMENT",
          "M.SEGMENT.POSTERIOR", "M.SEGMENT.FLAGGED", "ML.AR", "AR", "AR.ADJUSTED", "MAPPING.BIAS", "ML.LOH",
          "CN.SUBCLONAL", "CELLFRACTION", "FLAGGED", "log.ratio", "depth", "prior.somatic", "prior.contamination",
          "on.target", "seg.id", "gene.symbol"),
        c("Sampleid", "gene.symbol", "chr", "start", "end", "C", "seg.mean", "seg.id", "number.targets", "gene.mean",
          "gene.min", "gene.max", "focal", "breakpoints", "type", "num.snps.segment", "M", "M.flagged", "loh"),
        c("Sampleid", "chr", "start", "end", "arm", "C", "M", "type")
    )
    # sample id to use if no data avialable and mock df created
    sampid = sub(".csv","", basename(purecn_files[1]))  # NB: is this correct, or should there also be the "-nodups" ending to conform with samples with PureCN data?
    for (i in 1:length(purecn_files)) {
        if (file.size(purecn_files[i]) == 0) {  # if no valid PureCN output produced, create mock df
            assign(purecn_variables[i],
                   data.frame(matrix(c(sampid, rep(NA, length(purecn_colnames[[i]])-1)), nrow=1,
                                     dimnames=list(1, purecn_colnames[[i]])), stringsAsFactors = F))
        } else {  # if data available in the PureCN output, read it
            assign(purecn_variables[i], fread(purecn_files[i], sep = ',', stringsAsFactors = F))
        }
    }
    purecn_loh$cumstart=NA
    purecn_loh$cumend=NA
    for(i in 1:nrow(chrsizes)){
        ix <- which(purecn_loh$chr == chrsizes$chr[i])
        if (length(ix) == 0) next()  # skip if the chromosome is not present in the loh table
        purecn_loh$cumstart[ix] <- purecn_loh$start[ix] + chrsizes$cumstart[i]
        purecn_loh$cumend[ix] <- purecn_loh$end[ix] + chrsizes$cumstart[i]
    }
}

#save.image('ws.5.Rdata')


### Purity  ####
{
    # select high-confidence mutations for purity estimate (not AR)
    ix=salf$AO.T>=12 & salf$chromosome != "X" & salf$chromosome != "Y"
    median_af=median(salf$AF.T[ix],na.rm = T)
    # consider the mutation(s) being present @ 1 copy and there being 2 normal copies, what's the purity?
    t=median_af
    n=(1-median_af)/2
    p=t/(t+n)
    
    purity=purecn_stat
    purity$mutation_Purity=round(p,2)
    
    exportJson <- toJSON(purity)
    write(exportJson, opts$purity_json)
    
}

## Collect copy number for selected genes ####
cn_calls=rep('NO_CALL',nrow(genes))
names(cn_calls) = genes$label
## Get the AR copy number
g=genes[genes$label=='AR',]
data_line=data.frame(g)
cbins=bins[bins$chromosome==g$chromosome,] # this chromosome
left=cbins$end<g$start-3e6 # left control
left.median=median(cbins$log2[left],na.rm=T)
target=cbins$start>g$start-.1e6 & cbins$end<g$end+.1e6
target.median=median(cbins$log2[target],na.rm=T)
right=cbins$start>g$end+3e6 # right control
right.median=median(cbins$log2[right],na.rm=T)
target.diff=target.median-c(left.median,right.median) # difference compared to both controls
control.toUse=which.min(abs(target.diff)) # which control is most similar? left:1 and right:2
controlpoints=list(cbins$log2[left],cbins$log2[right]) # list of left and right control points
pval=wilcox.test(x = cbins$log2[target],y = controlpoints[[control.toUse]])$p.value # p value with correct control points
data_line$cnvkit.logr=target.median
data_line$cnvkit.to.control=round(target.diff[control.toUse],3)
data_line$cnvkit.pval=round(pval,5)
# For AR amp, require the AR signal to exceed control by 0.5 (there is no PureCN data for AR)
if (data_line$cnvkit.to.control > 0.5) cn_calls['AR']='AMPLIFIED'

## Check for amplifications
for (g in genes$label) {
    t=purecn_genes[purecn_genes$gene.symbol==g,]
    if (all(is.na(t))) next()  # skip if only mocked purecn data
    if (t$focal & t$C>3) cn_calls[g]='AMPLIFIED' #   calls amplified if focal and >3 copies
    if (t$C >= 7) cn_calls[g]='AMPLIFIED'            #  also if ≥7 copies
    if (!is.na(t$type)) if (t$type=='AMPLIFIED') cn_calls[g]='AMPLIFIED'   #           or if >5 copies regardless of focal status
}
## Check for focal deletion and LOH of remaining
for (g in genes$label) {
    t=purecn_genes[purecn_genes$gene.symbol==g,]
    if (all(is.na(t)) | cn_calls[g] == 'AMPLIFIED') next()  # skip if only mocked purecn data or if gene already marked as amplified
    if (nrow(t)>0){
        if (isTRUE(t$loh)) cn_calls[g]='LOSS_OF_HETEROZYGOSITY'
        if (t$focal & t$C==1) cn_calls[g]='FOCAL_DELETION'
        if (t$C==0) cn_calls[g]='HOMOZYGOUS_DELETION'
    }
}
## TMPRSS2 fusion is special case
g=genes[genes$label=='TMPRSS2',]
data_line=data.frame(g)
cbins=bins[bins$chromosome==g$chromosome,] # this chromosome
left=cbins$end<g$start # left control
left.median=median(cbins$log2[left],na.rm=T)
target=cbins$start>g$start-.1e6 & cbins$end<g$end+.1e6
target.median=median(cbins$log2[target],na.rm=T)
right=cbins$start>g$end # right control
right.median=median(cbins$log2[right],na.rm=T)
target.diff=target.median-c(left.median,right.median) # difference compared to both controls
control.toUse=which.min(abs(target.diff)) # which control is most similar? left:1 and right:2
controlpoints=list(cbins$log2[left],cbins$log2[right]) # list of left and right control points
pval=wilcox.test(x = cbins$log2[target],y = controlpoints[[control.toUse]])$p.value # p value with correct control points
data_line$cnvkit.logr=target.median
data_line$cnvkit.to.control=round(target.diff[control.toUse],3)
data_line$cnvkit.pval=round(pval,5)
# For TMPRSS-ERG del, require the signal to be 0.1 below nearest control
if (data_line$cnvkit.to.control < -.1) cn_calls['TMPRSS2']='IMPLIED_FUSION'

# Write copy numbers to JSON file
exportJson <- toJSON(cn_calls)
write(exportJson, opts$cna_json)

save.image('fr_T.Rdata')
### Frankenplot code - TUMOR  ####
{
    try( {
        # Prepare the per-SNP logratio
        alf[,log2:=as.numeric(NA)]
        # delta <- 3e6; for (i in 1:nrow(alf)) {
        #   ix <- bins$chromosome==alf$chromosome[i] & bins$start>alf$start[i]-delta & bins$end<alf$start[i]+delta
        #   alf$log2[i]=median(bins$log2[ix],na.rm=T)
        # }
        
        
        snpranges <- makeGRangesFromDataFrame(alf[,.(chromosome,start=start-1e2,end=end+1e2)])
        binranges <- makeGRangesFromDataFrame(bins)
        overlap <- as.data.table(findOverlaps(snpranges,binranges))
        bins[,smoothedLog2:=runmed(log2,k=1001),by='chromosome']
        alf[overlap$queryHits,log2:=bins[overlap$subjectHits]$smoothedLog2]
        
        # Prepare the smoothed AI
        ai <- smoothedAi <- alf$ai
        # for (i in 1:length(ai)) {
        #   ss=max(1,i-5); e=min(i+5,length(ai))
        #   smoothedAi[i]=median(ai[ss:e])
        # }
        smoothedAi <- runmed(ai,11)
    }, silent=T)
    
    ## File name here:
    destination_folder <- dirname(opts$plot_png)
    png(filename = opts$plot_png,width=11.7,height=8.3,units="in",res=600)
    
    ## set screens
    split.screen(figs=c(2,1))
    split.screen(as.matrix(data.frame(
        left=c(0,0.25),
        right=c(0.2,.9635),
        bottom=c(0,0.097),
        top=c(1,0.903))),1)
    split.screen(figs=c(2,1),3)
    split.screen(as.matrix(data.frame(left=c(rep(0.015,3)),
                                      right=c(rep(1,3)),
                                      bottom=c(0.10,0.5,0.55),
                                      top=c(0.5,0.55,1))),2)
    split.screen(figs=c(4,6),4)
    
    
    screen(1)
    #### Top margin with text #####
    try({
        plot(1,type='n',axes=F,xlab='',ylab='')
        snpcov='NA'; if (!is.null(alf)) if (nrow(alf)>0)
            snpcov=paste0('T',round(median(alf$td,na.rm=T)),'/N',round(median(alf$nd,na.rm=T)))
        mtext(
            text = paste0(purity$Sampleid,' ----- ', format(Sys.time(), "%F %H:%M:%S"),'',
                          '  SNPcov:',snpcov,
                          '  purCn.Ploidy:',round(purity$Ploidy,2),
                          '  purCn.Purity:',purity$Purity,
                          '  smut.Purity:',purity$mutation_Purity),
            side = 3,padj=-5.5,adj=1,cex=0.7)
    },silent=T)
    
    cex.axis <- .6
    cex.mtext <- 1.5
    cex.main <- 2
    cex.text <- .7
    axis1padj <- -2.35
    axis2hadj <- 0.3
    text1padj <- 1.7
    text2padj<- -3
    padj <- -3
    lwd=1
    color <- '#00000010'
    
    
    
    screen(5)
    #### Plot snp cov vs alleleratio ####
    xlim=c(1,5000)
    ylim=(0:1)
    try( {
        par(mar=c(2,3,2,.5),xaxs='i',yaxs='i',las=1)
        plot(1,type='n',xlim=xlim,ylim=ylim,axes=F,main='Variants',cex.main=0.5,log='x')
        text(x=c(2,4,6),y=0.48,labels = c('1 alt','2','3'),cex=0.6,col='grey',srt=-55)
        points(x = 7:1000,y=3/(7:1000),type='l',col='grey')
        axis(1,c(1,10,100,1000),lwd=lwd,lend=1,cex.axis=cex.axis,padj=axis1padj,tck=-0.03)
        axis(2,seq(0,1,.25),lwd=lwd,lend=1,cex.axis=cex.axis,hadj=axis2hadj,tck=-0.03)
        mtext('Coverage',1,cex=cex.text,padj=text1padj)
        mtext('Alt allele ratio',2,cex=cex.text,las=0,padj=text2padj)
        points(alf$td,alf$t,cex=0.3,col='#00000080',xlim=xlim,ylim=ylim,pch=16,lwd=lwd)
        points(alf$nd,alf$n,cex=0.1,col='#60606080',xlim=xlim,ylim=ylim,pch=3,lwd=lwd)
        scol=rep('#C00000CC',nrow(salf))
        scol[salf$FILTER!='PASS']='#500000CC'
        points(salf$DP.T,salf$AF.T,cex=0.4,col=scol,xlim=xlim,ylim=ylim,pch=salf$pch,lwd=lwd)
        segments(x0=median(alf$td,na.rm = T),y0=0,x1=median(alf$td,na.rm = T),y1=1,col='#00000080',lwd=lwd,lty=3)
        segments(x0=median(alf$nd,na.rm = T),y0=0,x1=median(alf$nd,na.rm = T),y1=1,col='#7070FF80',lwd=lwd,lty=3)
    }, silent=T)
    
    
    
    
    screen(6)
    try( {
        #####  Plot AR only  #####
        #xlim=c(55e6,80e6)
        ylim=c(-1.2,2)
        
        ar=genes[which(genes$label=="AR"),]
        par(mar=c(2,3,2,.5),xaxs='i',yaxs='i',las=1)
        ix=bins$chromosome=='X' & bins$gene!='Background' & bins$gene!='Antitarget'
        t=bins[ix,]
        plot(t$log2,ylim=ylim,pch=16,cex=0.3,main='AR',cex.main=0.6,axes=F)
        axis(2,c(-1,0,1,2),lwd=lwd,lend=1,cex.axis=cex.axis,hadj=axis2hadj,tck=-0.03)
        points(t$log2,type='l',col='#00000050')
        # rect(xleft = c(66763873, 66863097, 66905851, 66931243, 66937319, 66941674,
        #                66942668, 66943527),
        #      xright = c(66766604, 66863249, 66905968, 66931531, 66937464, 66941805,
        #                 66942826, 66950461),
        #      ybottom = -0.1,ytop = 0.1)
        mtext('Log ratio',2,cex=cex.text,las=0,padj=text2padj)
    }, silent=T)
    
    
    
    #### Imbalance vs logR panels ######
    
    screen(4)
    try( {
        par(mar=c(0,0,0,0))
        plot(1,type='n',axes=F,xlab='',ylab='')
        mtext('DNA ratio',1,padj=text1padj,cex=cex.text)
        mtext('Allelic imbalance',2,padj=-3,cex=cex.text)
        
        for (c in 1:24)
        {
            cex=0.4
            screen(c+9) # for chrY, screen is 33.
            par(mar = c(0, 0, 0, 0))
            par(oma = c(0,0,0,0))
            par(mgp =c(1,0.5,0))
            #xlim=c(-1.1,1.1)
            xlim=c(0.3,2.1)
            ylim=c(0,1)
            ix <- alf$chromosome %in% c(chrsizes$chr[c],'X','Y')
            ixCurChr <-  alf$chromosome %in% chrsizes$chr[c]
            plot(2^alf$log2[!ix],smoothedAi[!ix],xlim=xlim,ylim=ylim,lwd=lwd,axes=F,ylab='',xlab='',pch=16,
                 col='#B0B0B030',cex=cex)
            points(2^alf$log2[ixCurChr],smoothedAi[ixCurChr],col='#00800070',pch=16,cex=cex)
            
            
            if (c==24 & !is.null(alf)) {
                d=density(2^alf$log2[alf$chromosome %in% as.numeric(1:22)])
                points(d$x,d$y/max(d$y),type='l')
            }
            
            whole=(c(0.5,1,1.5,2))
            
            segments(
                x0=whole,x1=whole,
                y0=-0.05,y1=1.035,
                col='#D3D3D360',
                lwd=1)
            
            segments(
                x0=c(-2),x1=c(5),
                y0=c(1/3,1/2),y1=c(1/3,1/2),
                col='#D3D3D360',
                lwd=1)
            
            if(c == 23)
            {
                mtext("X", side = 3, line = -1, adj = 0.92, cex = 1)
            } else if(c == 24)
            {
                mtext("", side = 3, line = -1, adj = 0.92, cex = 1)
            } else
            {
                mtext(c, side = 3, line = -1, adj = 0.92, cex = 1)
            }
            
            if(c %in% 19:24) axis(side=1,cex.axis=0.5,at=c(0.5,1,1.5,2),tck=-0.05,padj=-1.6,#las=3,
                                  labels=c('-50%','±0','+50%','+100%'),
                                  col='white',col.ticks='black',lend=1)
            if(c %in% c(1,7,13,19)) axis(side=2,cex.axis=0.6,tck=-0.05,
                                         at=c(1/3,1/2),
                                         labels=c('2:1','3:1'),
                                         las=1,col='white',col.ticks='black',lend=1)
            if(c %in% c(6,12,18,24)) axis(side=4,cex.axis=0.6,tck=-0.05,
                                          at=c(1/3,1/2),
                                          labels=c('2:1','3:1'),
                                          las=1,col='white',col.ticks='black',lend=1)
            box(lwd=1)
        }
    }, silent=T)
    
    
    
    
    #### LogR track ####
    screen(9)
    try( {
        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(1,0.5,0))
        par(lend=1)
        
        ymin = -2
        ymax = 2
        #Plot signal over complete genome
        plot(NA,NA,#mpos[ix],mval[ix],
             pch=16,
             cex=0.3,
             main='',
             xlab = "",
             ylab = "",
             col = '#00000003',
             xaxt="n",
             axes=F,
             ylim = c(ymin,ymax),
             xlim = c(0,3095370729)
        )
        seqminmax <- seq(ymin,ymax,by=0.5)
        #Add axis to the left & right of signal
        axis(side=2,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=0,las=1)
        axis(side=4,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=3095370729,las=1)
        mtext("Log ratio",side=2,line=0,cex=0.7,padj=1.15)
        
        whole=(c(.5,1,1.5,2))
        #Add grey segments
        segments(
            y0=log2(whole),y1=log2(whole),
            x0=0,x1=3e9,
            col='#00000020',
            lwd=1)
        
        
        #Add genes as lines
        segments(
            x0=(genes$cumstart+genes$cumend)/2,
            y0=-100,y1=100,
            col='#0000C020',
            lwd=2)
        
        if (!is.null(bins)) {
            ix=bins$gene=='Background' | bins$gene=='Antitarget'
            points((bins$cumstart[!ix]+bins$cumend[!ix])/2,bins$log2[!ix], #ontargets
                   pch=1,
                   cex=0.7,
                   type='l',
                   col = '#00000020'
            )
            points((bins$cumstart[!ix]+bins$cumend[!ix])/2,bins$log2[!ix], #ontargets
                   pch=16,
                   cex=0.5,
                   #type='l',
                   col = '#00000080'
            )
            
            #Add segments
            col='#00C000CC'
            segments(x0=segments$cumstart,x1=segments$cumend,
                     y0=segments$log2,y1=segments$log2,
                     col=col,
                     lwd=3)
            #Add segments from PureCN
            segments(x0=0,x1=3e9,
                     y0=-1.8,y1=-1.8,
                     col='#00000010',
                     lwd=5)
            ix=purecn_loh$M==0 | purecn_loh$C==1
            if (sum(ix, na.rm=TRUE)>0) segments(x0=purecn_loh$cumstart[ix],x1=purecn_loh$cumend[ix],
                                                y0=-1.8,y1=-1.8,
                                                col='cyan',
                                                lwd=5)
            col=rep(NA,nrow(purecn_loh))
            col[purecn_loh$C==1]='blue'
            col[purecn_loh$C==0]='violet'
            col[purecn_loh$C==3]='red'
            col[purecn_loh$C>=4]='orange'
            segments(x0=purecn_loh$cumstart,x1=purecn_loh$cumend,
                     y0=-1.8,y1=-1.8,
                     col=col,
                     lwd=3)
            ix=purecn_loh$C==0
            ixSum=sum(ix, na.rm=TRUE)
            if (ixSum>0) points(x = (purecn_loh$cumstart[ix]+purecn_loh$cumend[ix])/2,y = rep(-1.8, ixSum),pch=24,bg='lightblue')
        }
        
        #Add a bar between chromosomes to distinguish them
        segments(
            x0=chrsizes$cumstart+chrsizes$size,
            y0=-100,y1=100,
            col='#00000099',
            lwd=1)
    }, silent=T)
    
    ## add any structural variants
    try ( {
        text(x = (t_strvs$cumstart+t_strvs$cumend)/2,y = 1.8,labels = t_strvs$type,col='red',cex=0.4,srt=90)
        text(x = (n_strvs$cumstart+n_strvs$cumend)/2,y = 1.8,labels = n_strvs$type,col='blue',cex=0.4,srt=90)
    }, silent=T)
    
    
    
    #Genes
    screen(8)
    try( {
        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(1,0.5,0))
        
        #Create a empty plot with the same ylim and xlim as the plot directly above it (signal) and below (AI)
        plot(0,0,xlab="",ylab="",main="",type="n",axes=F,xaxt="n",ylim=c(0,1),xlim=c(0,3095370729))
        
        chrsizes$labelpos[chrsizes$chr=='MT']=NA
        text(x = chrsizes$cumstart+0.5*chrsizes$size,y = 0.5,labels = chrsizes$chr,cex=.5)
        
        mtext(text='',side=2,las=1,line=-1.2)
        
        # Gene names
        text(x = (genes$cumstart+genes$cumend)[order(genes$chromosome, genes$start)]/2,y = rep(c(0.3,0.7), length.out=nrow(genes)),labels = genes$label[order(genes$chromosome, genes$start)],srt=45,cex=0.3,col='#0000C0CC')
        
    }, silent=T)
    
    
    
    
    #### Allele ratio  track  #####
    screen(7)
    
    try( {
        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(1,0.5,0))
        
        
        #Plot normal BAF over whole genome
        plot(NA,#alf$cumstart,alf$n,
             pch=16,
             cex=0.6,
             cex.axis=1,
             main='',
             xlab = "",
             ylab = "",
             col = "#80808080",
             xaxt="n",
             axes=F,
             ylim = c(-0.1,1.1),
             xlim = c(0,3095370729)
        )
        #Add axis to the left,right and below of AI. The below axis is the chromosome numbers 1-24.
        axis(side=2,tck=-0.04,at=seq(0,1,.2),cex.axis=0.6,pos=0,las=1) #at=c(0,0.25,0.33,0.5,0.67,0.75,1),labels=c('0','1/4','1/3','1/2','2/3','3/4','1')
        #axis(side=1,at=pre,pos=0,labels=c(seq(from=1,to=22),"X",'Y'),cex.axis=0.50,lty=0)#,tck=0,col.ticks='#00000000')
        axis(side=4,tck=-0.04,at=seq(0,1,.2),cex.axis=0.6,pos=3095370729,las=1) #at=c(0,0.25,0.33,0.5,0.67,0.75,1),labels=c('0','1/4','1/3','1/2','2/3','3/4','1')
        mtext("SNP allele ratio",side=2,line=0,cex=0.7,padj = 1.15)
        
        #Add a bar between chromosomes to distinguish them
        segments(
            x0=chrsizes$cumstart+chrsizes$size,
            y0=-100,y1=100,
            col='#00000099',
            lwd=1)
        
        #Add genes as lines
        segments(
            x0=(genes$cumstart+genes$cumend)/2,
            y0=-100,y1=100,
            col='#0000C020',
            lwd=2)
        
        segments(
            y0=c(0,1/4,1/3,2/3,3/4,1),y1=c(0,1/4,1/3,2/3,3/4,1),
            x0=0,x1=3095370729,
            # col='#D3D3D360',
            col='#00000020',
            lwd=1)
        
        points(alf$cumstart,alf$t,
               pch=16,
               cex=0.6,
               col = "#00000040")
        
        
        
    }, silent=T)
    try( {
        
        
        ## Add somatic mutations
        if (nrow(salf)>0) {
            scol=rep('#C00000CC',nrow(salf))
            ix=salf$AF.T<0.02 | salf$AO.T<6 | salf$FILTER!="PASS"
            scol[ix]='#500000CC'
            points(salf$cumpos[ix],salf$AF.T[ix],
                   pch=salf$pch[ix],
                   cex=0.6,
                   bg=scol[ix]
            )
            points(salf$cumpos[!ix],salf$AF.T[!ix],
                   pch=salf$pch[!ix],
                   cex=0.6,
                   bg=scol[!ix]
            )
            
            g <- pos <- aa <- rep('',nrow(salf))
            for (j in 1:nrow(salf)) {
                g[j]=as.character(salf$SYMBOL[j])
                g[is.na(g)]=''
                pos[j]=strsplit(as.character(salf$Protein_position[j]),'/')[[1]][1]
                aa[j]=strsplit(as.character(salf$Amino_acids[j]),'/')[[1]][2]
                pos[is.na(pos)]=''; pos[pos=='-']=''; aa[is.na(aa)]=''
            }
            ix=pos!='' & g %in% genes$label  & aa!=''
            if (sum(ix)>0) text(x = salf$cumpos[ix],y=salf$AF.T[ix]-0.07,labels = paste0(g[ix],' ',pos[ix],aa[ix]),cex=0.4,srt=30,col=scol[ix])
        }
        
    }, silent=T)
    try( {
        
        ## Add germline mutations
        scol=rep('#0000C0CC',nrow(galf))
        
        ix=(galf$hasConsequence) & galf$AO>=12 & galf$AF>=0.2
        
        points(galf$cumpos[ix],galf$AF[ix],
               pch=galf$pch[ix],
               cex=0.6,
               bg=scol[ix]
        )
        if (nrow(galf)>0) {
            g <- pos <- aa <- rep('',nrow(galf))
            for (j in 1:nrow(galf)) {
                g[j]=as.character(galf$SYMBOL[j])
                g[is.na(g)]=''
                pos[j]=strsplit(as.character(galf$Protein_position[j]),'/')[[1]][1]
                aa[j]=strsplit(as.character(galf$Amino_acids[j]),'/')[[1]][2]
                pos[is.na(pos)]=''; pos[pos=='-']=''; aa[is.na(aa)]=''
            }
            ix=ix & pos!='' & g %in% genes$label & aa!=''
            text(x = galf$cumpos[ix],y=galf$AF[ix]+0.07,labels = paste0(g[ix],' ',pos[ix],aa[ix]),cex=0.4,srt=30,col=scol[ix])
        }
        
    }, silent=T)
    
    #Close all the opened split.screens and release the figure
    
    close.screen(all.screens=T)
    dev.off()
}

save.image('fr_N.Rdata')

### Frankenplot code - NORMAL  ####
{
    try( {
        # Prepare the per-SNP logratio (NORMAL)   <<--- Overwrites the alf$log2 with normal, but is not used anymore.
        alf$log2 <- NA
        # delta <- 3e6; for (i in 1:nrow(alf)) {
        #   ix <- bins_n$chromosome==alf$chromosome[i] & bins_n$start>alf$start[i]-delta & bins_n$end<alf$start[i]+delta
        #   alf$log2[i]=median(bins_n$log2[ix],na.rm=T)
        # }
        snpranges <- makeGRangesFromDataFrame(alf[,.(chromosome,start=start-1e2,end=end+1e2)])
        binranges <- makeGRangesFromDataFrame(bins_n)
        overlap <- as.data.table(findOverlaps(snpranges,binranges))
        bins_n[,smoothedLog2:=runmed(log2,k=1001),by='chromosome']
        alf[overlap$queryHits,log2:=bins_n[overlap$subjectHits]$smoothedLog2]
        
        # Prepare the smoothed AI (NORMAL)     <<--- Overwrites the smoothedAi with normal, but is not used anymore.
        ai <- smoothedAi <- alf$ai_n
        # for (i in 1:length(ai)) {
        #   ss=max(1,i-5); e=min(i+5,length(ai))
        #   smoothedAi[i]=median(ai[ss:e])
        # }
        smoothedAi <- runmed(ai,11)
    }, silent=T)
    
    ## File name here:
    destination_folder <- dirname(opts$plot_png)
    png(filename = opts$plot_png_normal,width=11.7,height=8.3,units="in",res=600)
    ## set screens
    split.screen(figs=c(2,1))
    split.screen(as.matrix(data.frame(
        left=c(0,0.25),
        right=c(0.2,.9635),
        bottom=c(0,0.097),
        top=c(1,0.903))),1)
    split.screen(figs=c(2,1),3)
    split.screen(as.matrix(data.frame(left=c(rep(0.015,3)),
                                      right=c(rep(1,3)),
                                      bottom=c(0.10,0.5,0.55),
                                      top=c(0.5,0.55,1))),2)
    split.screen(figs=c(4,6),4)
    
    
    screen(1)
    #### Top margin with text #####
    try({
        plot(1,type='n',axes=F,xlab='',ylab='')
        snpcov='NA'; if (!is.null(alf)) if (nrow(alf)>0)
            snpcov=paste0('T',round(median(alf$td,na.rm=T)),'/N',round(median(alf$nd,na.rm=T)))
        mtext(
            text = paste0(purity$Sampleid,' --NORMAL-- ', format(Sys.time(), "%F %H:%M:%S"),'',
                          '  SNPcov:',snpcov),
            side = 3,padj=-5.5,adj=1,cex=0.7)
    },silent=T)
    
    cex.axis <- .6
    cex.mtext <- 1.5
    cex.main <- 2
    cex.text <- .7
    axis1padj <- -2.35
    axis2hadj <- 0.3
    text1padj <- 1.7
    text2padj<- -3
    padj <- -3
    lwd=1
    color <- '#00000010'
    
    
    
    screen(5)
    #### Plot snp cov vs alleleratio ####
    xlim=c(1,5000)
    ylim=(0:1)
    try( {
        par(mar=c(2,3,2,.5),xaxs='i',yaxs='i',las=1)
        plot(1,type='n',xlim=xlim,ylim=ylim,axes=F,main='Variants',cex.main=0.5,log='x')
        text(x=c(2,4,6),y=0.48,labels = c('1 alt','2','3'),cex=0.6,col='grey',srt=-55)
        points(x = 7:1000,y=3/(7:1000),type='l',col='grey')
        axis(1,c(1,10,100,1000),lwd=lwd,lend=1,cex.axis=cex.axis,padj=axis1padj,tck=-0.03)
        axis(2,seq(0,1,.25),lwd=lwd,lend=1,cex.axis=cex.axis,hadj=axis2hadj,tck=-0.03)
        mtext('Coverage',1,cex=cex.text,padj=text1padj)
        mtext('Alt allele ratio',2,cex=cex.text,las=0,padj=text2padj)
        points(alf$td,alf$t,cex=0.3,col='#00000080',xlim=xlim,ylim=ylim,pch=16,lwd=lwd)
        points(alf$nd,alf$n,cex=0.1,col='#60606080',xlim=xlim,ylim=ylim,pch=3,lwd=lwd)
        #scol=rep('#C00000CC',nrow(salf))
        #scol[salf$FILTER!='PASS']='#500000CC'
        #points(salf$DP.T,salf$AF.T,cex=0.4,col=scol,xlim=xlim,ylim=ylim,pch=salf$pch,lwd=lwd)
        segments(x0=median(alf$td,na.rm = T),y0=0,x1=median(alf$td,na.rm = T),y1=1,col='#00000080',lwd=lwd,lty=3)
        segments(x0=median(alf$nd,na.rm = T),y0=0,x1=median(alf$nd,na.rm = T),y1=1,col='#7070FF80',lwd=lwd,lty=3)
    }, silent=T)
    
    
    
    
    screen(6)
    try( {
        #####  Plot AR only  #####
        #xlim=c(55e6,80e6)
        ylim=c(-1.2,2)
        
        ar=genes[which(genes$label=="AR"),]
        par(mar=c(2,3,2,.5),xaxs='i',yaxs='i',las=1)
        ix=bins_n$chromosome=='X' & bins_n$gene!='Background' & bins_n$gene!='Antitarget'
        t=bins_n[ix,]
        plot(t$log2,ylim=ylim,pch=16,cex=0.3,main='AR',cex.main=0.6,axes=F)
        axis(2,c(-1,0,1,2),lwd=lwd,lend=1,cex.axis=cex.axis,hadj=axis2hadj,tck=-0.03)
        points(t$log2,type='l',col='#00000050')
        # rect(xleft = c(66763873, 66863097, 66905851, 66931243, 66937319, 66941674,
        #                66942668, 66943527),
        #      xright = c(66766604, 66863249, 66905968, 66931531, 66937464, 66941805,
        #                 66942826, 66950461),
        #      ybottom = -0.1,ytop = 0.1)
        mtext('Normal Log ratio',2,cex=cex.text,las=0,padj=text2padj)
    }, silent=T)
    
    
    
    #### Imbalance vs logR panels ######
    
    screen(4)
    try( {
        par(mar=c(0,0,0,0))
        plot(1,type='n',axes=F,xlab='',ylab='')
        mtext('Normal DNA ratio',1,padj=text1padj,cex=cex.text)
        mtext('Normal Allelic imbalance',2,padj=-3,cex=cex.text)
        
        for (c in 1:24)
        {
            cex=0.4
            screen(c+9) # for chrY, screen is 33.
            par(mar = c(0, 0, 0, 0))
            par(oma = c(0,0,0,0))
            par(mgp =c(1,0.5,0))
            #xlim=c(-1.1,1.1)
            xlim=c(0.3,2.1)
            ylim=c(0,1)
            ix <- alf$chromosome %in% c(chrsizes$chr[c],'X','Y')
            ixCurChr <-  alf$chromosome %in% chrsizes$chr[c]
            plot(2^alf$log2[!ix],smoothedAi[!ix],xlim=xlim,ylim=ylim,lwd=lwd,axes=F,ylab='',xlab='',pch=16,
                 col='#B0B0B030',cex=cex)
            points(2^alf$log2[ixCurChr],smoothedAi[ixCurChr],col='#00800070',pch=16,cex=cex)
            
            
            if (c==24 & !is.null(alf)) {
                d=density(2^alf$log2[alf$chromosome %in% as.numeric(1:22)])
                points(d$x,d$y/max(d$y),type='l')
            }
            
            whole=(c(0.5,1,1.5,2))
            
            segments(
                x0=whole,x1=whole,
                y0=-0.05,y1=1.035,
                col='#D3D3D360',
                lwd=1)
            
            segments(
                x0=c(-2),x1=c(5),
                y0=c(1/3,1/2),y1=c(1/3,1/2),
                col='#D3D3D360',
                lwd=1)
            
            if(c == 23)
            {
                mtext("X", side = 3, line = -1, adj = 0.92, cex = 1)
            } else if(c == 24)
            {
                mtext("", side = 3, line = -1, adj = 0.92, cex = 1)
            } else
            {
                mtext(c, side = 3, line = -1, adj = 0.92, cex = 1)
            }
            
            if(c %in% 19:24) axis(side=1,cex.axis=0.5,at=c(0.5,1,1.5,2),tck=-0.05,padj=-1.6,#las=3,
                                  labels=c('-50%','±0','+50%','+100%'),
                                  col='white',col.ticks='black',lend=1)
            if(c %in% c(1,7,13,19)) axis(side=2,cex.axis=0.6,tck=-0.05,
                                         at=c(1/3,1/2),
                                         labels=c('2:1','3:1'),
                                         las=1,col='white',col.ticks='black',lend=1)
            if(c %in% c(6,12,18,24)) axis(side=4,cex.axis=0.6,tck=-0.05,
                                          at=c(1/3,1/2),
                                          labels=c('2:1','3:1'),
                                          las=1,col='white',col.ticks='black',lend=1)
            box(lwd=1)
        }
    }, silent=T)
    
    
    
    
    #### LogR track ####
    screen(9)
    try( {
        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(1,0.5,0))
        par(lend=1)
        
        ymin = -2
        ymax = 2
        #Plot signal over complete genome
        plot(NA,NA,#mpos[ix],mval[ix],
             pch=16,
             cex=0.3,
             main='',
             xlab = "",
             ylab = "",
             col = '#00000003',
             xaxt="n",
             axes=F,
             ylim = c(ymin,ymax),
             xlim = c(0,3095370729)
        )
        seqminmax <- seq(ymin,ymax,by=0.5)
        #Add axis to the left & right of signal
        axis(side=2,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=0,las=1)
        axis(side=4,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=3095370729,las=1)
        mtext("Normal Log ratio",side=2,line=0,cex=0.7,padj=1.15)
        
        whole=(c(.5,1,1.5,2))
        #Add grey segments
        segments(
            y0=log2(whole),y1=log2(whole),
            x0=0,x1=3e9,
            col='#00000020',
            lwd=1)
        
        
        #Add genes as lines
        segments(
            x0=(genes$cumstart+genes$cumend)/2,
            y0=-100,y1=100,
            col='#0000C020',
            lwd=2)
        
        if (!is.null(bins_n)) {
            ix=bins_n$gene=='Background' | bins_n$gene=='Antitarget'
            points((bins_n$cumstart[!ix]+bins_n$cumend[!ix])/2,bins_n$log2[!ix], #ontargets
                   pch=1,
                   cex=0.7,
                   type='l',
                   col = '#00000020'
            )
            points((bins_n$cumstart[!ix]+bins_n$cumend[!ix])/2,bins_n$log2[!ix], #ontargets
                   pch=16,
                   cex=0.5,
                   #type='l',
                   col = '#00000080'
            )
            
            #Add segments
            col='#00C000CC'
            segments(x0=segments_n$cumstart,x1=segments_n$cumend,
                     y0=segments_n$log2,y1=segments_n$log2,
                     col=col,
                     lwd=3)
            #   #Add segments from PureCN
            #   segments(x0=0,x1=3e9,
            #            y0=-1.8,y1=-1.8,
            #            col='#00000010',
            #            lwd=5)
            #   ix=purecn_loh$M==0 | purecn_loh$C==1
            #   if (sum(ix, na.rm=TRUE)>0) segments(x0=purecn_loh$cumstart[ix],x1=purecn_loh$cumend[ix],
            #                                       y0=-1.8,y1=-1.8,
            #                                       col='cyan',
            #                                       lwd=5)
            #   col=rep(NA,nrow(purecn_loh))
            #   col[purecn_loh$C==1]='blue'
            #   col[purecn_loh$C==0]='violet'
            #   col[purecn_loh$C==3]='red'
            #   col[purecn_loh$C>=4]='orange'
            #   segments(x0=purecn_loh$cumstart,x1=purecn_loh$cumend,
            #            y0=-1.8,y1=-1.8,
            #            col=col,
            #            lwd=3)
            #   ix=purecn_loh$C==0
            #   ixSum=sum(ix, na.rm=TRUE)
            #   if (ixSum>0) points(x = (purecn_loh$cumstart[ix]+purecn_loh$cumend[ix])/2,y = rep(-1.8, ixSum),pch=24,bg='lightblue')
        }
        
        #Add a bar between chromosomes to distinguish them
        segments(
            x0=chrsizes$cumstart+chrsizes$size,
            y0=-100,y1=100,
            col='#00000099',
            lwd=1)
    }, silent=T)
    
    ## add any structural variants  (NORMALS only)
    try ( {
        #text(x = (t_strvs$cumstart+t_strvs$cumend)/2,y = 1.8,labels = t_strvs$type,col='red',cex=0.4,srt=90)
        text(x = (n_strvs$cumstart+n_strvs$cumend)/2,y = 1.8,labels = n_strvs$type,col='blue',cex=0.4,srt=90)
    }, silent=T)
    
    
    
    #Genes
    screen(8)
    try( {
        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(1,0.5,0))
        
        #Create a empty plot with the same ylim and xlim as the plot directly above it (signal) and below (AI)
        plot(0,0,xlab="",ylab="",main="",type="n",axes=F,xaxt="n",ylim=c(0,1),xlim=c(0,3095370729))
        
        chrsizes$labelpos[chrsizes$chr=='MT']=NA
        text(x = chrsizes$cumstart+0.5*chrsizes$size,y = 0.5,labels = chrsizes$chr,cex=.5)
        
        mtext(text='',side=2,las=1,line=-1.2)
        
        # Gene names
        text(x = (genes$cumstart+genes$cumend)[order(genes$chromosome, genes$start)]/2,y = rep(c(0.3,0.7), length.out=nrow(genes)),labels = genes$label[order(genes$chromosome, genes$start)],srt=45,cex=0.3,col='#0000C0CC')
        
    }, silent=T)
    
    
    
    
    #### Allele ratio  track  #####
    screen(7)
    
    try( {
        #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
        par(mar = c(0, 0, 0, 0))
        par(oma = c(0,0,0,0))
        par(mgp =c(1,0.5,0))
        
        
        #Plot normal BAF over whole genome
        plot(NA,#alf$cumstart,alf$n,
             pch=16,
             cex=0.6,
             cex.axis=1,
             main='',
             xlab = "",
             ylab = "",
             col = "#80808080",
             xaxt="n",
             axes=F,
             ylim = c(-0.1,1.1),
             xlim = c(0,3095370729)
        )
        #Add axis to the left,right and below of AI. The below axis is the chromosome numbers 1-24.
        axis(side=2,tck=-0.04,at=seq(0,1,.2),cex.axis=0.6,pos=0,las=1) #at=c(0,0.25,0.33,0.5,0.67,0.75,1),labels=c('0','1/4','1/3','1/2','2/3','3/4','1')
        #axis(side=1,at=pre,pos=0,labels=c(seq(from=1,to=22),"X",'Y'),cex.axis=0.50,lty=0)#,tck=0,col.ticks='#00000000')
        axis(side=4,tck=-0.04,at=seq(0,1,.2),cex.axis=0.6,pos=3095370729,las=1) #at=c(0,0.25,0.33,0.5,0.67,0.75,1),labels=c('0','1/4','1/3','1/2','2/3','3/4','1')
        mtext("Normal SNP allele ratio",side=2,line=0,cex=0.7,padj = 1.15)
        
        #Add a bar between chromosomes to distinguish them
        segments(
            x0=chrsizes$cumstart+chrsizes$size,
            y0=-100,y1=100,
            col='#00000099',
            lwd=1)
        
        #Add genes as lines
        segments(
            x0=(genes$cumstart+genes$cumend)/2,
            y0=-100,y1=100,
            col='#0000C020',
            lwd=2)
        
        segments(
            y0=c(0,1/4,1/3,2/3,3/4,1),y1=c(0,1/4,1/3,2/3,3/4,1),
            x0=0,x1=3095370729,
            # col='#D3D3D360',
            col='#00000020',
            lwd=1)
        
        points(alf$cumstart,alf$n,
               pch=16,
               cex=0.6,
               col = "#00000040")
        
        
        
    }, silent=T)
    try( {
        
        
        # ## Add somatic mutations
        # if (nrow(salf)>0) {
        #   scol=rep('#C00000CC',nrow(salf))
        #   ix=salf$AF.T<0.02 | salf$AO.T<6 | salf$FILTER!="PASS"
        #   scol[ix]='#500000CC'
        #   points(salf$cumpos[ix],salf$AF.T[ix],
        #          pch=salf$pch[ix],
        #          cex=0.6,
        #          bg=scol[ix]
        #   )
        #   points(salf$cumpos[!ix],salf$AF.T[!ix],
        #          pch=salf$pch[!ix],
        #          cex=0.6,
        #          bg=scol[!ix]
        #   )
        #
        #   g <- pos <- aa <- rep('',nrow(salf))
        #   for (j in 1:nrow(salf)) {
        #     g[j]=as.character(salf$SYMBOL[j])
        #     g[is.na(g)]=''
        #     pos[j]=strsplit(as.character(salf$Protein_position[j]),'/')[[1]][1]
        #     aa[j]=strsplit(as.character(salf$Amino_acids[j]),'/')[[1]][2]
        #     pos[is.na(pos)]=''; pos[pos=='-']=''; aa[is.na(aa)]=''
        #   }
        #   ix=pos!='' & g %in% genes$label  & aa!=''
        #   if (sum(ix)>0) text(x = salf$cumpos[ix],y=salf$AF.T[ix]-0.07,labels = paste0(g[ix],' ',pos[ix],aa[ix]),cex=0.4,srt=30,col=scol[ix])
        # }
        
    }, silent=T)
    # try( {
    #
    #   ## Add germline mutations
    #   scol=rep('#0000C0CC',nrow(galf))
    #
    #   ix=(galf$hasConsequence) & galf$AO>=12 & galf$AF>=0.2
    #
    #   points(galf$cumpos[ix],galf$AF_n[ix],
    #          pch=galf$pch[ix],
    #          cex=0.6,
    #          bg=scol[ix]
    #   )
    #   if (nrow(galf)>0) {
    #     g <- pos <- aa <- rep('',nrow(galf))
    #     for (j in 1:nrow(galf)) {
    #       g[j]=as.character(galf$SYMBOL[j])
    #       g[is.na(g)]=''
    #       pos[j]=strsplit(as.character(galf$Protein_position[j]),'/')[[1]][1]
    #       aa[j]=strsplit(as.character(galf$Amino_acids[j]),'/')[[1]][2]
    #       pos[is.na(pos)]=''; pos[pos=='-']=''; aa[is.na(aa)]=''
    #     }
    #     ix=ix & pos!='' & g %in% genes$label & aa!=''
    #     text(x = galf$cumpos[ix],y=galf$AF[ix]+0.07,labels = paste0(g[ix],' ',pos[ix],aa[ix]),cex=0.4,srt=30,col=scol[ix])
    #   }
    #
    # }, silent=T)
    
    #Close all the opened split.screens and release the figure
    
    close.screen(all.screens=T)
    dev.off()
}


## Normal Plot, first version ####
if (F) try( {
    ## File name here:
    png(filename = opts$plot_png_normal,width=11.7,height=2.1,units="in",res=600)
    
    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(1, 1, 1, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(1,0.5,0))
    par(lend=1)
    
    ymin = -2
    ymax = 2
    #Plot signal over complete genome
    plot(NA,NA,#mpos[ix],mval[ix],
         pch=16,
         cex=0.3,
         main=sub(".cnr","",basename(opts$normal_cnr)),
         xlab = "",
         ylab = "",
         col = '#00000003',
         xaxt="n",
         axes=F,
         ylim = c(ymin,ymax),
         xlim = c(0,3095370729)
    )
    seqminmax <- seq(ymin,ymax,by=0.5)
    #Add axis to the left & right of signal
    axis(side=2,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=0,las=1)
    axis(side=4,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=3095370729,las=1)
    mtext("Log ratio",side=2,line=0,cex=0.7,padj=1.15)
    
    whole=(c(.5,1,1.5,2))
    #Add grey segments
    segments(
        y0=log2(whole),y1=log2(whole),
        x0=0,x1=3e9,
        col='#00000020',
        lwd=1)
    
    
    #Add genes as lines
    segments(
        x0=(genes$cumstart+genes$cumend)/2,
        y0=-100,y1=100,
        col='#0000C020',
        lwd=2)
    
    if (!is.null(bins_n)) {
        ix=bins_n$gene=='Background' | bins_n$gene=='Antitarget'
        points((bins_n$cumstart[!ix]+bins_n$cumend[!ix])/2,bins_n$log2[!ix], #ontargets
               pch=1,
               cex=0.7,
               type='l',
               col = '#00000020'
        )
        points((bins_n$cumstart[!ix]+bins_n$cumend[!ix])/2,bins_n$log2[!ix], #ontargets
               pch=16,
               cex=0.5,
               #type='l',
               col = '#00000080'
        )
        
        #Add segments
        col='#00C000CC'
        segments(x0=segments_n$cumstart,x1=segments_n$cumend,
                 y0=segments_n$log2,y1=segments_n$log2,
                 col=col,
                 lwd=3)
    }
    
    #Add a bar between chromosomes to distinguish them
    segments(
        x0=chrsizes$cumstart+chrsizes$size,
        y0=-100,y1=100,
        col='#00000099',
        lwd=1)
    
    # chromosomes
    text(x = chrsizes$cumstart+0.5*chrsizes$size,y = ymin,labels = chrsizes$chr,cex=.5)
    
    mtext(text='',side=2,las=1,line=-1.2)
    
    # Gene names
    text(x = (genes$cumstart+genes$cumend)[order(genes$chromosome, genes$start)]/2,y = ymin+rep(c(0.1,0.4), length.out=nrow(genes)),labels = genes$label[order(genes$chromosome, genes$start)],srt=45,cex=0.4,col='#0000C0CC')
    
    ## add any structural variants
    try ( {
        text(x = (n_strvs$cumstart+n_strvs$cumend)/2,y = ymax-0.2,labels = n_strvs$type,col='blue',cex=0.4,srt=90)
    }, silent=T)
    
    dev.off()
    
}, silent=T)

save.image('ws.Rdata')

## INTERACTIVE plot code from here: #############################
#p<-ggplot()+
#  #ggtitle("Plot of length by dose")   +
#  xlim(0,3095370729) + theme_bw() +
#  theme(
#    panel.border = element_blank(),
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    axis.line = element_line(colour = "black"),
#    #axis.text.y = element_text(margin=margin(5,5,10,5,"pt")),     # Change y axis tick labels only
#    #axis.text.y.right = element_text(),# y  axis tick labels on top axis
#    axis.ticks = element_line(),
#    #axis.title.y = element_text("Log ratio"),
#    axis.title=element_text(face="bold.italic",
#                            size="12", color="blue"),
#    axis.title.x=element_blank(),
#    axis.text.x=element_blank(),
#    axis.ticks.x=element_blank(),
#    axis.line.x = element_blank(),
#    legend.position="top",
#    axis.ticks.length = unit(0.25, "cm")
#  ) +
#  labs(y = "Log ratio", x = "", colour = "blue")+
#  scale_y_continuous(
#    breaks = seq(-2,2, by = 0.5),
#    labels = seq(-2,2, by = 0.5),
#    expand = c(0,0),
#    limits = c(-2,2),
#    sec.axis = dup_axis()
#  ) +
#  geom_hline(yintercept=log2(c(.5,1,1.5,2)), linetype="dashed", color='grey' ) +
#  geom_vline(aes(xintercept =(genes$cumstart+genes$cumend)/2, text=sprintf("gene:  %s", genes$label[order(genes$chromosome, genes$start)]), colour='blue'), color='violet' )

#p<-p+
#  geom_line(aes((bins_n$cumstart[!ix]+bins_n$cumend[!ix])/2,bins_n$log2[!ix]))+
#  geom_point(aes((bins_n$cumstart[!ix]+bins_n$cumend[!ix])/2,bins_n$log2[!ix])) +
#  geom_segment(
#    aes(x = segments_n$cumstart, y = segments_n$log2, xend = segments_n$cumend, yend = segments_n$log2,
#        text=sprintf("CNVs:<br>: start: %f <br> end: %f <br> length: %f", segments_n$cumstart, segments_n$cumend, segments_n$cumend-segments_n$cumstart + 1 )
#        ),
#    size=0.80,
#    show.legend=F,
#    colour="green"
#  ) +
#  geom_vline(aes(xintercept = chrsizes$cumstart+chrsizes$size, text=sprintf("chromosome:  %s", chrsizes$chr))) +
#
#  geom_text(aes(x = chrsizes$cumstart+0.5*chrsizes$size,
#                y = ymin+0.2,label = chrsizes$chr, text=sprintf("chromosome:  %s", chrsizes$chr)),
#               show.legend=F, size=2.5) #+
#
#
# p <- p + geom_text(aes(x = (n_strvs$cumstart+n_strvs$cumend)/2,y = ymax-0.2,label = n_strvs$type),
#            col='blue',cex=2,srt=90, size=2)
#
# p <- p +  annotate("text",
#                 x = (genes$cumstart+genes$cumend)[order(genes$chromosome, genes$start)]/2,
#                 y = ymin+rep(c(0.1,0.8),length.out=nrow(genes)),
#                 label = genes$label[order(genes$chromosome, genes$start)] ,
#                 color="orange",
#                size=2.5 , angle=45)
#
#p<-ggplotly(p, tooltip = "text")
#x<-plotly_json(p, FALSE)
#filename<-paste(destination_folder,'/',purity$Sampleid,"_normal.json", sep="")
#write(x, filename)

#############end################################


# downsample bins and alf to 300,000

set.seed(25)
if (exists('bins')) if (!is.null(bins)) if (nrow(bins)>300e3) {
    bins <- bins[order(rnorm(nrow(bins)))[1:300e3]]
}
if (exists('alf')) if (!is.null(alf)) if (nrow(alf)>300e3) {
    alf <- alf[order(rnorm(nrow(alf)))[1:300e3]]
}



######################## alt allele ratio #################
xlim=c(1,3000)
ylim=(0:1)

destination_folder <- dirname(opts$plot_png)

p<-ggplot()+
    ggtitle("Variants")   +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5)
    )+
    xlim(xlim) +
    ylim(ylim) +
    geom_text(aes(x=c(2,4,6),y=0.48,label = c('1 alt', '2','3'), angle=-55),colour='grey')+
    geom_point(aes(x = 7:1000,y=3/(7:1000)), cex=0.1, colour='grey', type="line") +
    geom_line(aes(x = 7:1000,y=3/(7:1000), size=3), cex=0.1, colour='grey')+
    scale_x_continuous(trans = "log2",name="Coverage", limit=c(1,3000), breaks=c(1,10,100,1000),
                       label=c("1","10","100","1000")) +
    scale_y_continuous(name="Alt allele ration" , limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
    geom_point(aes(x=alf$td,y=alf$t),cex=0.3,fill='#00000080', color='black') +
    geom_point(aes(x=alf$nd, y=alf$n),cex=0.1, color='grey')

scol=rep('#C00000CC',nrow(salf))
scol[salf$FILTER!="PASS"]='#500000CC'

p<- p + geom_point(aes(salf$DP.T,salf$AF.T),size=1,colour=scol, fill=scol, pch=salf$pch) +
    geom_vline(aes(xintercept =(median(alf$td,na.rm = T))   ), linetype="dashed", colour="grey") +
    geom_vline(aes(xintercept =(median(alf$nd,na.rm = T))   ), linetype="dashed", colour="grey")

#m<-plotly_json(p)
x<-plotly_json(p, FALSE)
filename<-paste(destination_folder, '/' , purity$Sampleid,"_variant_allelic.json", sep="")
write(x, filename)

###### Plot AR #################
ylim=c(-1.2,2)
ar=genes[which(genes$label=="AR"),]
ix=bins$chromosome=='X' & bins$gene!='Background' & bins$gene!='Antitarget'
t=bins[ix,]
p<-ggplot(t)+
    ggtitle("AR")   +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x.bottom = element_blank(),
        plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(name="Log ratio" , limits=ylim ) +
    geom_line(aes(x=seq(1,length(t$log2)), y=t$log2), colour='grey',show.legend = FALSE)+
    geom_point(aes(x=seq(1,length(t$log2)), y=t$log2, text=sprintf(
        "**AR**<br>
  chr   : %s <br>
  start : %f <br>
  end   : %f <br>
  depth : %f <br>
  gene  : %s <br>
  ", t$chromosome, t$start, t$end, t$depth, t$gene)), cex=0.8,  show.legend = FALSE)

p<-ggplotly(p, tooltip = "text")
x<-plotly_json(p, FALSE)
filename<-paste(destination_folder, '/' , purity$Sampleid,"_plot_AR.json", sep="")
write(x, filename)

###########################################SNP Allelic ratio####################################################

p<-ggplot()+
    xlim(0,3095370729) + theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(),
        axis.title=element_text(face="bold.italic", size="12", color="blue"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),
        legend.position="top",
        axis.ticks.length = unit(0.25, "cm")
    ) +
    labs(y = "Log ratio", x = "", colour = "blue")+
    scale_y_continuous(
        breaks = seq(-2,2, by = 0.5),
        labels = seq(-2,2, by = 0.5),
        expand = c(0,0),
        sec.axis = dup_axis()
    )
# limits = c(-2,2),
p <- p + geom_hline(yintercept=log2(c(.5,1,1.5,2)),  color='grey', size=1 , alpha=0.4) +
    geom_vline(aes(xintercept =(genes$cumstart+genes$cumend)/2,
                   text=sprintf("chr: %s <br>gene:  %s <br>start:%2f <br>end:%2f",
                                genes$chromosome,
                                genes$label,
                                genes$start,
                                genes$end
                   ),
                   colour=genes$chromosome), col='#0000C020', show.legend = FALSE )

if (!is.null(bins)) {
    ix=bins$gene=='Background' | bins$gene=='Antitarget'
    p<-p+
        geom_line(aes((bins$cumstart[!ix]+bins$cumend[!ix])/2,bins$log2[!ix]),  col = '#00000020')+
        geom_point(aes((bins$cumstart[!ix]+bins$cumend[!ix])/2,bins$log2[!ix]), cex=0.5) +
        geom_segment(
            aes(x = segments$cumstart, y = segments$log2, xend = segments$cumend, yend = segments$log2,
                text=sprintf("CNVs:<br>chr: %s <br> log2: %f <br> start: %f <br> end: %f <br> depth: %f <br> probes: %f <br>",
                             segments$chromosome,
                             segments$log2,
                             segments$start,
                             segments$end,
                             segments$depth,
                             segments$probes)
            ),
            size=0.80,
            show.legend=F,
            col='#00C000CC'
        )+
        geom_hline(yintercept=-1.8,  color='grey', size=4 , alpha=0.2)
    
    
    ##Add segments from PureCN
    ix=purecn_loh$M==0 | purecn_loh$C==1
    if (sum(ix, na.rm=TRUE)>0){
        
        xstart_n <- purecn_loh$cumstart[ix]
        xend_n   <- purecn_loh$cumend[ix]
        chr      <- purecn_loh$chr[ix]
        type     <-  purecn_loh$type[ix]
        arm      <- purecn_loh$arm[ix]
        cumstart_n <- purecn_loh$cumstart[ix]
        cumend_nd   <- purecn_loh$cumend[ix]
        c <-  purecn_loh$C[ix]
        df = data.frame(xstart_n, xend_n, chr,  type, arm , cumstart_n, cumend_nd, c  )
        
        col_now=rep(NA,nrow(df))
        col_now[df$c==1]= '#4169E1' #'blue'
        col_now[df$c==0]= '#EE82EE' #'violet'
        col_now[df$c==3]= '#FF0000' #'red'
        col_now[df$c>=4]= '#FFA500' #'orange'
        purecn_loh$C[purecn_loh$C>=4]
        p <- p + geom_segment(data=df,
                              aes(x = df$xstart_n, y = -1.8,
                                  xend = df$xend_n, yend = -1.8,
                                  text=sprintf("chr:%s <br> start: %2f <br> end: %2f <br>type: %s <br> arm:%s  ",
                                               chr,
                                               xstart_n,
                                               xend_n,
                                               type,
                                               arm
                                  )
                              ),
                              size=4,
                              show.legend=F,
                              col='cyan'
        )
        
    }
    
    
    col=rep(NA,nrow(purecn_loh))
    col[purecn_loh$C==1]='blue'
    col[purecn_loh$C==0]='violet'
    col[purecn_loh$C==3]='red'
    col[purecn_loh$C>=4]='orange'
    if(sum(purecn_loh$cumstart, na.rm=TRUE)>0){
        p <- p + geom_segment(data=purecn_loh,
                              aes(x = purecn_loh$cumstart, y = -1.8,
                                  xend = purecn_loh$cumend , yend = -1.8,
                                  text=sprintf("chr:%s <br> start: %2f <br> end: %2f <br>type: %s <br> arm:%s  ",
                                               purecn_loh$chr,
                                               purecn_loh$cumstart,
                                               purecn_loh$cumend,
                                               purecn_loh$type,
                                               purecn_loh$arm )
                                  
                              ),
                              size=1,
                              col =  col,
                              show.legend=F)
    }
    
    ix=purecn_loh$C==0
    if (sum(ix, na.rm=TRUE)>0){
        #p <- p + geom_point(aes(x=(purecn_loh$cumstart[ix]+purecn_loh$cumend[ix])/2, y=-1.8,
        #                        text=sprintf("chr:%s <br> start: %2f <br> end: %2f <br>type: %s <br> arm:%s  ",
        #                                     purecn_loh$chr,
        #                                     purecn_loh$cumstart,
        #                                     purecn_loh$cumend,
        #                                     purecn_loh$type,
        #                                     purecn_loh$arm )
        #
        #                        ), shape =24, size =2, colour="blue", fill = "blue")
        
        xstart_n <- purecn_loh$cumstart[ix]
        xend_n   <- purecn_loh$cumend[ix]
        chr      <- purecn_loh$chr[ix]
        type     <-  purecn_loh$type[ix]
        arm      <- purecn_loh$arm[ix]
        cumstart <- purecn_loh$cumstart[ix]
        cumend   <- purecn_loh$cumend[ix]
        c <-  purecn_loh$C[ix]
        df_c_new = data.frame(xstart_n, xend_n, chr,  type, arm , cumstart, cumend, c )
        
        p <- p + geom_point(data = df_c_new, aes(x=(df_c_new$cumstart+df_c_new$cumend)/2, y=-1.8,
                                                 text=sprintf("chr:%s <br> start: %2f <br> end: %2f <br>type: %s <br> arm:%s  ",
                                                              df_c_new$chr,
                                                              df_c_new$cumstart,
                                                              df_c_new$cumend,
                                                              df_c_new$type,
                                                              df_c_new$arm )
                                                 
        ), shape =24, size =2, colour="blue", fill = "blue")
        
    }
    
}

# chromosome lines
#p <- p + geom_vline(aes(xintercept = chrsizes$cumstart+chrsizes$size,
#                        text=sprintf("chromosome:  %s", chrsizes$chr))) +
#  geom_text(aes(x = (n_strvs$cumstart+n_strvs$cumend)/2,y = ymax-0.2,label = n_strvs$type,
#                text=sprintf("chr:%s, <br> gene_feature: %s <br> strand:%s", n_strvs$chr, n_strvs$V3, n_strvs$V7)
#                ),
#            col='blue',cex=2,srt=90, size=4) +
# add any stuctural variants
#
#  geom_text(aes(x = (t_strvs$cumstart+t_strvs$cumend)/2,y = ymax-0.2,label = t_strvs$type,
#                text=sprintf("chr:%s, <br> gene_feature: %s <br> strand:%s", t_strvs$chr, t_strvs$V3, t_strvs$V7)
#                ),
#          col='red',cex=2,srt=90, size=4)


try(  {
    p <- p + geom_vline(aes(xintercept = chrsizes$cumstart+chrsizes$size,
                            text=sprintf("chromosome:  %s", chrsizes$chr)))
    if(!is.null(dim(n_strvs))) {
        p <- p + geom_text(aes(x = (n_strvs$cumstart+n_strvs$cumend)/2,y = ymax-0.2,label = n_strvs$type,
                               text=sprintf("chr:%s, <br> gene_feature: %s <br> strand:%s", n_strvs$chr, n_strvs$V3, n_strvs$V7)
        ),
        col='blue',cex=2,srt=90, size=4)
    }
    # add any stuctural variants
    if(!is.null(dim(t_strvs))) {
        p <- p + geom_text(aes(x = (t_strvs$cumstart+t_strvs$cumend)/2,y = ymax-0.2,label = t_strvs$type,
                               text=sprintf("chr:%s, <br> gene_feature: %s <br> strand:%s", t_strvs$chr, t_strvs$V3, t_strvs$V7)
        ),
        col='red',cex=2,srt=90, size=4)
    }
    
}, silent=T)


p<-ggplotly(p, tooltip = "text", height = 1024) # this p will be used in subplot
#####################snp plot with germline and somatic#######################################################
p1<-ggplot()+
    xlim(0,3095370729) + theme_bw() +
    ylim(-0.1,1.1 ) +
    theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(),
        axis.title=element_text(face="bold.italic", size="12", color="blue"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),
        legend.position="top",
        #plot.margin = c(0, 0, 0, 0),
        axis.ticks.length = unit(0.25, "cm")
    ) +
    labs(y = "SNP Allelic Ratio", x = "", colour = "blue")+
    scale_y_continuous(
        breaks = seq(-0.1, 1.1, by = 0.2),
        labels = seq(-0.1, 1.1, by = 0.2),
        #expand = c(0,0),
        limits = c(-0.1,1.1),
        sec.axis = dup_axis()
    ) +
    geom_vline(aes(xintercept = chrsizes$cumstart+chrsizes$size,
                   text=sprintf("chromosome:  %s", chrsizes$chr))) +
    geom_hline(yintercept=c(0,0.25,0.33,0.66,0.75,1),  color='grey', size=1 , alpha=0.4) +
    geom_vline(aes(xintercept =(genes$cumstart+genes$cumend)/2,
                   text=sprintf("chr: %s <br>gene:  %s <br>start:%2f <br>end:%2f",
                                genes$chromosome,
                                genes$label,
                                genes$start,
                                genes$end
                   ),
                   colour=genes$chromosome), col='#0000C020', show.legend = FALSE ) +
    
    geom_point(aes(x=alf$cumstart, y=alf$t),pch=16,cex=1,col = "#00000040")

## Add somatic mutations

if (nrow(salf)>0) {
    scol=rep('#C00000CC',nrow(salf))
    ix=salf$AF.T<0.02 | salf$AO.T<6 | salf$FILTER!="PASS"
    scol[ix]='#500000CC'
    
    p1 <- p1 +  geom_point(aes(x=salf$cumpos[ix], y=salf$AF.T[ix],
                               text=sprintf("chr: %s <br>pos:  %f <br>Ref:%s <br>Alt:%s <br> Symbol: %s <br> Biotype: %s",
                                            salf$chromosome[ix],
                                            salf$pos[ix],
                                            salf$REF[ix],
                                            salf$ALT[ix],
                                            salf$SYMBOL[ix],
                                            salf$BIOTYPE[ix]
                               )
    ),pch=salf$pch[ix],cex=1,bg=scol[ix]) +
        geom_point(aes(x=salf$cumpos[!ix], y=salf$AF.T[!ix],
                       text=sprintf("chr: %s <br>pos:  %f <br>Ref:%s <br>Alt:%s <br> Symbol: %s <br> Biotype: %s",
                                    salf$chromosome[!ix],
                                    salf$pos[!ix],
                                    salf$REF[!ix],
                                    salf$ALT[!ix],
                                    salf$SYMBOL[!ix],
                                    salf$BIOTYPE[!ix]
                       )
        ),pch=salf$pch[!ix],cex=1,bg=scol[!ix])
    
    g <- pos <- aa <- rep('',nrow(salf))
    for (j in 1:nrow(salf)) {
        g[j]=as.character(salf$SYMBOL[j])
        g[is.na(g)]=''
        pos[j]=strsplit(as.character(salf$Protein_position[j]),'/')[[1]][1]
        aa[j]=strsplit(as.character(salf$Amino_acids[j]),'/')[[1]][2]
        pos[is.na(pos)]=''; pos[pos=='-']=''; aa[is.na(aa)]=''
    }
    #ix=pos!='' & g %in% genes$label  & aa!=''
    #if (sum(ix)>0){
    #
    #  p1 <- p1 + geom_text(aes(x = salf$cumpos[ix] , y = salf$AF.T[ix]-0.07, paste0(g[ix],' ',pos[ix],aa[ix])),
    #                     cex=0.4,srt=30,col=scol[ix])
    #}
    
}

## Add germline mutations
g_scol=rep('#0000C0CC',nrow(galf))
g_ix=(galf$hasConsequence) & galf$AO>=12 & galf$AF>=0.2
p1 <- p1 + geom_point(aes(x=galf$cumpos[g_ix], y=galf$AF[g_ix],
                          text=sprintf("chr: %s <br>pos:  %f <br>Ref:%s <br>Alt:%s <br> Symbol: %s <br> Biotype: %s",
                                       galf$chromosome[g_ix],
                                       galf$pos[g_ix],
                                       galf$REF[g_ix],
                                       galf$ALT[g_ix],
                                       galf$SYMBOL[g_ix],
                                       galf$BIOTYPE[g_ix]
                          )
),pch=galf$pch[g_ix],cex=1,bg=g_scol[g_ix])

if (nrow(galf)>0) {
    g_g <- pos_g <- aa_g <- rep('',nrow(galf))
    for (j in 1:nrow(galf)) {
        g_g[j]=as.character(galf$SYMBOL[j])
        g_g[is.na(g_g)]=''
        pos_g[j]=strsplit(as.character(galf$Protein_position[j]),'/')[[1]][1]
        aa_g[j]=strsplit(as.character(galf$Amino_acids[j]),'/')[[1]][2]
        pos_g[is.na(pos_g)]=''; pos_g[pos_g=='-']=''; aa_g[is.na(aa_g)]=''
    }
    ix_g=g_ix & pos_g!='' & g_g %in% genes$label & aa_g!=''
    if (sum(ix_g, na.rm = TRUE)>0) {
        p1 <- p1 + geom_text(aes(x = galf$cumpos[ix_g],y=galf$AF[ix_g]+0.07, label = paste0(g_g[ix_g],' ',pos_g[ix_g],
                                                                                            aa_g[ix_g])),cex=1,srt=30,col=g_scol[ix_g])
        #text(x = galf$cumpos[ix],y=galf$AF[ix]+0.07,labels = paste0(g[ix],' ',pos[ix],aa[ix]),cex=0.4,srt=30,col=scol[ix])
    }
}

p1<-ggplotly(p1, tooltip = "text",  height = 1024)

### draw gene track with ggplot#####
p2<-ggplot()+
    xlim(0,3095370729) + theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(),
        axis.title=element_text(face="bold.italic", size="12", color="blue"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),
        legend.position="top",
        axis.ticks.length = unit(0.25, "cm")
    ) +
    labs(y = "Log ratio", x = "", colour = "blue")+
    scale_y_continuous(
        breaks = seq(-2,2, by = 0.5),
        labels = seq(-2,2, by = 0.5),
        expand = c(0,0),
        limits = c(-2,2),
        sec.axis = dup_axis()
    )

p2 <- p2 + geom_rect(data= data2_mRNA,
                     mapping=aes(xmin=data2_mRNA$cumstart,
                                 xmax=data2_mRNA$cumend,
                                 ymin=data2_mRNA$yaxis+0.05,
                                 ymax=data2_mRNA$yaxis-0.05,
                                 text=sprintf("chr: %s <br>start: %f <br>end: %f <br>gene_id:%s <br> transcript_id: %s <br> Symbol: %s  <br>strand: %s<br>type:gene",
                                              data2_mRNA$Chromosome,
                                              data2_mRNA$Start,
                                              data2_mRNA$End,
                                              data2_mRNA$gene_id,
                                              data2_mRNA$transcript_id,
                                              data2_mRNA$gene_name,
                                              data2_mRNA$Strand)))





p2 <- p2 + geom_segment(data=data2_3p, mapping=aes(x = data2_3p$cumstart, y = data2_3p$yaxis, xend = data2_3p$cumend, yend = data2_3p$yaxis,
                                                   
                                                   text=sprintf("chr: %s <br>start: %f <br>end: %f <br>gene_id:%s <br> transcript_id: %s <br> Symbol: %s  <br>strand: %s",
                                                                data2_3p$Chromosome,
                                                                data2_3p$Start,
                                                                data2_3p$End,
                                                                data2_3p$gene_id,
                                                                data2_3p$transcript_id,
                                                                data2_3p$gene_name,
                                                                data2_3p$Strand)
                                                   
), size = 9, color='red' )

p2 <- p2 + geom_segment(data=data2_5p, mapping=aes(x = data2_5p$cumstart, y = data2_5p$yaxis, xend = data2_5p$cumend, yend = data2_5p$yaxis
                                                   
                                                   ,text=sprintf("chr: %s <br>start: %f <br>end: %f <br>gene_id:%s <br> transcript_id: %s <br> Symbol: %s  <br>strand: %s",
                                                                 data2_5p$Chromosome,
                                                                 data2_5p$Start,
                                                                 data2_5p$End,
                                                                 data2_5p$gene_id,
                                                                 data2_5p$transcript_id,
                                                                 data2_5p$gene_name,
                                                                 data2_5p$Strand)
                                                   
), size = 9, color='blue' )



p2 <- p2 + geom_segment(data=data2_exon, mapping=aes(x = data2_exon$cumstart, y = data2_exon$yaxis, xend = data2_exon$cumend, yend = data2_exon$yaxis
                                                     ,            text=sprintf("chr: %s <br>start: %f <br>end: %f <br>gene_id:%s <br> transcript_id: %s <br> Symbol: %s  <br>strand: %s",
                                                                               data2_exon$Chromosome,
                                                                               data2_exon$Start,
                                                                               data2_exon$End,
                                                                               data2_exon$gene_id,
                                                                               data2_exon$transcript_id,
                                                                               data2_exon$gene_name,
                                                                               data2_exon$Strand)
                                                     
), size = 7, color='green' )


x <- subplot(p, p1, p2, shareX = TRUE, nrows=3, which_layout = 1)
x<-plotly_json(x, FALSE)
filename<-paste(destination_folder,'/',purity$Sampleid,"_snpratio.json", sep="")
write(x, filename)

# ###########scatter plot with 24 boxes allelic frequency####################
# #draw title
# library(tibble)
# library(gridExtra)
# xlim=c(0.3,2.1)
# ylim=c(0,1)
# snpcov='NA'; if (!is.null(alf)) if (nrow(alf)>0)
#   snpcov=paste0('T',round(median(alf$td)),'/N',round(median(alf$nd)))
# purity$Ploidy = ifelse(is.na(purity$Ploidy), 0, purity$Ploidy)
# text = paste0(purity$Sampleid,' ----- ', format(Sys.time(), "%F %H:%M:%S"),'',
#               '  SNPcov:',snpcov,
#               '  purCn.Ploidy:',round(purity$Ploidy,2),
#               '\n purCn.Purity:',purity$Purity,
#               '  smut.Purity:',purity$mutation_Purity)
# 
# newalf = tibble(chrom=alf$chromosome, log2 = 2^alf$log2, smoothedAi_new = smoothedAi)
# newalf$chrom <- factor(newalf$chrom,levels=c(1:22,"X"))
# d_bg <- newalf[, c('log2','smoothedAi_new') ]
# whole=(c(0.5,1,1.5,2))
# p4 <- ggplot( newalf, aes(x=log2, y=smoothedAi_new)) + #, aes(x = 2^alf$log2, y = smoothedAi)) +
#   #ggtitle(text) +
#   theme_bw() +
#   theme(
#     panel.border = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black"),
#     plot.subtitle = element_text("|",hjust = 0.5, size = 8)
#     #plot.title = element_text(hjust = 0.5, size = 8)
#   )+
#   scale_x_continuous(name="DNA ratio", limit=c(0.3,2.1), breaks=c(0.5,1.0,1.5,2.0),
#                      label=c("-50%","±0","+50","+100"))+
#   scale_y_continuous(name="Allelic imbalance", limit=c(0,1), breaks=c(1/3, 1/2),
#                      label=c("2:1","3:1"))+
#   geom_vline(xintercept=(c(0.5,1,1.5,2)), colour='grey', size=0.1) +
#   geom_hline(yintercept=(c(1/3,1/2)), colour='grey', size=0.1) +
#   geom_point(data = d_bg, colour = "grey", cex=0.6)+
#   geom_point(colour="green", cex=0.4) +
#   facet_wrap(~ chrom, nrow=4, ncol = 6, strip.position = "bottom", as.table = TRUE)
# 
# p4<-ggplotly(p4, width = 1024, height = 800)
# x<-plotly_json(p4, FALSE)
# filename<-paste(destination_folder,'/', purity$Sampleid,"_scatter_plot.json", sep="")
# write(x, filename)
# 
# ########draw density plot ###################
# ix <- alf$chromosome %in% c(chrsizes$chr[24],'X','Y')
# ixCurChr <-  alf$chromosome %in% chrsizes$chr[24]
# d=density(2^alf$log2[alf$chromosome %in% as.numeric(1:22)])
# p5 <- ggplot()+
#   ggtitle(text) +
#   theme_bw() +
#   theme(
#     panel.border = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black"),
#     plot.subtitle = element_text("|",hjust = 0.5, size = 8),
#     plot.title = element_text(hjust = 0.5, size = 8)
#   )+
#   scale_x_continuous(name="DNA ratio", limit=c(0.3,2.1), breaks=c(0.5,1.0,1.5,2.0),
#                      label=c("-50%","±0","+50","+100"))+
#   scale_y_continuous(name="Allelic imbalance", limit=c(0,1), breaks=c(1/3, 1/2),
#                      label=c("2:1","3:1"))+
#   geom_vline(xintercept=(c(0.5,1,1.5,2)), colour='grey', size=0.1) +
#   geom_hline(yintercept=(c(1/3,1/2)), colour='grey', size=0.1, cex=0.4) +
#   geom_point(aes(2^alf$log2[!ix],smoothedAi[!ix]), colour='grey', cex=0.4)+
#   geom_line(aes(x=d$x,y=d$y/max(d$y)))
# 
# p5 <-ggplotly(p5)
# x<-plotly_json(p5, FALSE)
# filename<-paste(destination_folder, '/', purity$Sampleid,"_scatter_density_plot.json", sep="")
# write(x, filename)
# #####redraw the density plot ######
# ix <- alf$chromosome %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X', 'Y', 'MT')
# ixCurChr <-  alf$chromosome %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22, 'X', 'Y', 'MT')
# d=density(alf$log2[alf$chromosome %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X', 'Y', 'MT')])
# p3 <- ggplot()+
#   ggtitle(text) +
#   theme_bw() +
#   theme(
#     panel.border = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black"),
#     plot.subtitle = element_text("|",hjust = 0.5, size = 8),
#     plot.title = element_text(hjust = 0.5, size = 8),
#     axis.title.x=element_blank(),
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank(),
#     axis.line.x = element_blank()
# 
#   ) +
#   scale_y_continuous(
#     breaks = seq(-2,2, by = 0.5),
#     labels = seq(-2,2, by = 0.5),
#     expand = c(0,0),
#     limits = c(-2,2)
#   ) +
#   geom_path(aes(x=d$y/max(d$y), y=d$x))
# 
#   p3 <-ggplotly(p3, width = 200, height = 1024)
#   x1<-plotly_json(p3, FALSE)
#   filename<-paste(destination_folder, '/', purity$Sampleid,"_snpratio_density.json", sep="")
#   write(x1, filename)
# ####################################
