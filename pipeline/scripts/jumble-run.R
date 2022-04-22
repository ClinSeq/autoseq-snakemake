#!/usr/bin/env Rscript

# Markus Mayrhofer 2022
# Dependencies ------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bamsignals))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(PSCBS))
suppressPackageStartupMessages(library(ggplot2)); theme_set(theme_bw())
suppressPackageStartupMessages(library(patchwork))


# Options ------------------------------------------------------------
option_list <- list(
    make_option(c("-r", "--reference_file"), action = "store", type = "character",default = NULL, 
                help = "reference file e.g. <design>.bed.references.RDS"),
    make_option(c("-b", "--input_bam"), action = "store", type = "character",default = NULL, 
                help = "input bam or bam.counts.RDS file"),
    make_option(c("-n", "--normal_bam"), action = "store", type = "character",default = NULL, 
                help = "optional matched normal bam or bam.counts.RDS file"),
    make_option(c("-v", "--snp_vcf"), action = "store", type = "character",default = NULL, 
                help = "optional vcf file with SNPs, including any matched normal"),
    make_option(c("-o", "--output_dir"), action = "store", type = "character",default = '.', 
                help = "directory for output")
)

opt <- parse_args(OptionParser(option_list = option_list))


# #with RDS
# setwd('~/Data/CLINSEQ/')
# opt <- list(reference_file='comprehensive3_baits_twist.bed.reference.RDS',
#             input_bam='sample1.bam.counts.RDS',
#             normal_bam='sample1.bam.counts.RDS',
#             snp_vcf='sample1n.vcf.gz',
#             output_dir='.'
# )


# #with bam
# setwd('~/Data/CLINSEQ/')
# opt <- list(reference_file='comprehensive3_baits_twist.bed.reference.RDS',
#             input_bam='sample1.bam',
#             normal_bam=NULL,
#             snp_vcf='sample1.vcf',
#             output_dir='.'
# )


# #PSFF with RDS and vcf
# setwd('~/Data/CLINSEQ/')
# opt <- list(reference_file='pancancer2_baits_twist.bed.reference.RDS',
#             input_bam='psff_PN/iPCM-P-00436099-T-MOL2594i-KH20211123-PN20211124_nodups.bam.counts.RDS',
#             snp_vcf='psff_vcf/iPCM-P-00436099-T-MOL2594i-KH-PN-iPCM-P-00436099-N-03998767-KH-PN.vardict-somatic-purecn.vcf.gz',
#             output_dir='psff_pdf/'
# )


# Reference file ------------------------------------------------------------

reference <- readRDS(opt$reference_file)
counts_template <- reference[c("target_bed_file","chromlength","target_ranges","background_ranges")]






# Fragment counts ------------------------------------------------------------

# function for bam process
countsFromBam <- function(counts,bampath) {
    counts$date_count <- date()
    counts$input_bam_file <- bampath
    target_ranges <- counts$target_ranges
    background_ranges <- counts$background_ranges
    counts[['target_count']] <- bamCount(bampath, target_ranges, paired.end="midpoint", 
                                         mapq=20, filteredFlag=1024, verbose=F)
    counts[['target_short']] <- bamCount(bampath, target_ranges, paired.end="midpoint", 
                                         mapq=20, filteredFlag=1024, tlenFilter=c(0,160), verbose=F)
    counts[['background_count']] <- bamCount(bampath, background_ranges, paired.end="midpoint", 
                                             mapq=20, filteredFlag=1024, verbose=F)
    counts[['background_short']] <- bamCount(bampath, background_ranges, paired.end="midpoint", 
                                             mapq=20, filteredFlag=1024, tlenFilter=c(0,160), verbose=F)
    return(counts)
    }

# parse main sample
input <- opt$input_bam
if (str_detect(input,'.RDS$')) {
    if (!file.exists(input)) stop(paste("Cannot find",input))
    counts <- readRDS(input)
} else if (str_detect(input,'.[bB][aA][mM]$')) {
    if (!file.exists(input)) stop(paste("Cannot find",input))
    counts <- countsFromBam(counts_template,input)
} else stop("No input file?")

# parse optional normal
input <- opt$normal_bam
counts_normal <- NULL
normal <- FALSE
if (!is.null(input)) {
    normal <- TRUE
    if (str_detect(input,'.RDS$')) {
        if (!file.exists(input)) stop(paste("Cannot find",input))
        counts_normal <- readRDS(input)
    } else if (str_detect(input,'.[bB][aA][mM]$')) {
        if (!file.exists(input)) stop(paste("Cannot find",input))
        counts_normal <- countsFromBam(counts_template,input)
    } else stop("No normal input file?")
}



# Tables of bins ------------------------------------------------------------


targets <- reference$targets

background <- reference$background
background$gene='Background'

name <- str_remove(opt$input_bam,'.*/')
name <- str_remove(name,'\\.counts.RDS')

targets$sample <- name
background$sample <- name

if (normal) {
    name_n <- str_remove(opt$normal_bam,'.*/')
    name_n <- str_remove(name_n,'\\.counts.RDS')
    targets$sample_n <- name_n
}

# Add counts
targets[,count:=counts$target_count]
sum_by_bg <- targets[,sum(count),by=background]

background[,raw_count:=counts$background_count]
background[,count:=raw_count][sum_by_bg$background,count:=raw_count-sum_by_bg$V1]

# same, short
targets[,count_short:=counts$target_short]
sum_by_bg <- targets[,sum(count_short),by=background]
background[,raw_short:=counts$background_short]
background[,count_short:=raw_short][sum_by_bg$background,count_short:=raw_short-sum_by_bg$V1]


# median correct to backbone and set raw log ratio
targets[,rawLR:=log2(count+1)][,rawLR:=rawLR-median(rawLR[is_backbone]),by=sample]
targets[,rawLR_short:=log2(count_short+1)][,rawLR_short:=rawLR_short-median(rawLR_short[is_backbone]),by=sample]

background[,rawLR:=log2(count+1)][,rawLR:=rawLR-median(rawLR[is_backbone]),by=sample]
background[,rawLR_short:=log2(count_short+1)][,rawLR_short:=rawLR_short-median(rawLR_short[is_backbone]),by=sample]

if (normal) {
    # add counts & median correct (standard only) for the normal
    targets[,count_n:=counts_normal$target_count]
    
    sum_by_bg <- targets[,sum(count_n),by=background]
    
    background[,raw_n:=counts_normal$background_count]
    background[,count_n:=raw_n][sum_by_bg$background,count_n:=raw_n-sum_by_bg$V1]
    
    targets[,rawLR_n:=log2(count_n+1)][,rawLR_n:=rawLR_n-median(rawLR_n[is_backbone]),by=sample]
    background[,rawLR_n:=log2(count_n+1)][,rawLR_n:=rawLR_n-median(rawLR_n[is_backbone]),by=sample]
}

# keep only bins "ok" in this reference
targets <- targets[reference$keep_targets]
background <- background[reference$keep_background]

target_ranges <- makeGRangesFromDataFrame(targets)
background_ranges <- makeGRangesFromDataFrame(background)



# Noise correction ------------------------------------------------------------

# correct using reference
jcorrect <- function(temp,train_ix=NULL) {
    
    if (is.null(train_ix)) train_ix <- rep(TRUE,nrow(temp))
    
    loess_temp=loess(lr ~ PC1+PC2+PC3, data = temp, 
                     subset = train_ix, 
                     family="symmetric", control = loess.control(surface = "direct"))
    temp[,lr:=lr-predict(loess_temp,temp)]
    
    # loess_temp=loess(lr ~ gc, data = temp,
    #                  subset = train_ix,
    #                  family="symmetric", control = loess.control(surface = "direct"))
    # temp[,lr:=lr-predict(loess_temp,temp)]
    
    if (F) for (i in 1:(ncol(temp)-1)) {
        temp$thispc=temp[[paste0('PC',i)]]
        loess_temp=rlm(lr ~ thispc, data=temp,
                       subset = targets$is_backbone)
        temp[,lr:=lr-predict(loess_temp,temp)]
    }
    return(temp$lr)
}

# standard 
temp <- cbind(data.table(lr=targets$rawLR),reference$targets_ref)
targets[,log2:=jcorrect(temp,targets$is_backbone)]

# short
temp <- cbind(data.table(lr=targets$rawLR_short),reference$targets_ref_short)
targets[,log2_short:=jcorrect(temp,targets$is_backbone)]

# standard bg
temp <- cbind(data.table(lr=background$rawLR),reference$background_ref)
background[,log2:=jcorrect(temp,is_backbone)]

# bg short
temp <- cbind(data.table(lr=background$rawLR_short),reference$background_ref_short)
background[,log2_short:=jcorrect(temp,is_backbone)]

if (normal) {
    # correct the normal
    temp <- cbind(data.table(lr=targets$rawLR_n),reference$targets_ref)
    targets[,log2_n:=jcorrect(temp,targets$is_backbone)]
    
    temp <- cbind(data.table(lr=background$rawLR_n),reference$background_ref)
    background[,log2_n:=jcorrect(temp,is_backbone)]
}


# ggplot(background) + geom_point(aes(x=background,y=log2)) +
#     geom_point(data=targets,mapping = aes(x=background,y=log2),col='blue')


# SNPs ------------------------------------------------------------
#save.image('ws.Rdata')

targets$allele_ratio <- as.double(NA)
targets$allele_ratio_n <- as.double(NA)

snp_allele_ratio <- FALSE
input <- opt$snp_vcf
if (!is.null(input)) {
    snp_allele_ratio <- TRUE
    if (!str_detect(input,'.[vV][cC][fF]$') & !str_detect(input,'.[vV][cC][fF].[gG][zZ]$')) stop("SNP vcf file appears incorrect")
    
    vcf=readVcf(input)
    
    # Get alleles
    ref_allele <- as.character(ref(vcf))
    alt_allele <- as.data.table(alt(vcf))[, .(values = list(value)), by = group][,values]
    n_alt_alleles <- sapply(alt_allele, length)
    alt_allele <- sapply(alt_allele, "[[", 1)
    
    # Keep only simple substitutions
    vcf <- vcf[str_length(ref_allele)==1 & str_length(alt_allele)==1 & n_alt_alleles==1]
    
    # Get alleles again
    ref_allele <- as.character(ref(vcf))
    alt_allele <- as.data.table(alt(vcf))[, .(values = list(value)), by = group][,values]
    alt_allele <- sapply(alt_allele, "[[", 1)
    from <- rep('C/G',length(vcf)); from[ref_allele %in% c('A','T')] <- 'A/T'
    to <- rep('C/G',length(vcf)); to[alt_allele %in% c('A','T')] <- 'A/T'
    type <- paste(from,'>',to)
    
    # Prepare table of SNPs for this patient
    snp_table <- data.table(snp=1:length(vcf),target=as.integer(NA),id=names(vcf),type)
    
    # Used for filtering
    min_DP <- rep(Inf,length(vcf))
    max_AD <- rep(-Inf,length(vcf))
    max_RD <- rep(-Inf,length(vcf))
    
    # If there is a matched normal, which sample is the primary?
    t_ix <- 1 # assumes first is "tumor" or main sample
    n_ix <- NA 
    names <- colnames(vcf)
    if (normal) {
        n_ix <- 2 # assumes second one is the matched normal
        if (str_detect(counts$input_bam_file,names[2])) { # unless the second sample name matches the main/tumor bam
            n_ix <- 1
            t_ix <- 2
        }
    } else if (length(names)==2) {
        if (str_detect(counts$input_bam_file,names[2])) { # unless the second sample name matches the main bam
            t_ix <- 2
        }
    }
    
    
    # Iterate samples in this vcf file
    for (j in 1:ncol(vcf)) {
        
        
        g <- geno(vcf[,j])
        
        
        AD <- sapply(g$AD, "[[", 2)
        RD <- sapply(g$AD, "[[", 1)
        DP=AD+RD
        
        
        raw_allele_ratio <- unname(round(AD/DP,4))
        raw_allele_ratio[is.nan(raw_allele_ratio)] <- NA
        
        # for use in filtering
        min_DP <- apply(data.table(DP,min_DP),1,min,na.rm=T)
        max_AD <- apply(data.table(AD,max_AD),1,max,na.rm=T)
        max_RD <- apply(data.table(RD,max_RD),1,max,na.rm=T)
        
        # Selection to use for alt/ref correction model
        ix <- AD>50 & RD > 50 & raw_allele_ratio > .4 & raw_allele_ratio < .6 & DP < quantile(DP,.95)*1.1
        
        # ggplot() + geom_point(aes(x=RD[ix]+AD[ix], y=AD[ix]/RD[ix],fill=type[ix]),shape=21) +
        #     scale_y_log10() + 
        #     geom_smooth(aes(x=RD[ix]+AD[ix], y=AD[ix]/RD[ix],group=type[ix],col=type[ix]))
        # ggplot() + geom_point(aes(x=1:sum(ix), y=AD[ix]/RD[ix],fill=type[ix]),shape=21)
        
        allele_ratio <- raw_allele_ratio
        
        # compute error factors
        
        if (F) {
            temp <- data.table(AD,RD,depth=AD+RD,ratio=AD/RD,type=type,error_factor=1.0,mapd=1)
            # m <- rlm(ratio ~ depth, temp[ix])
            # error_factor <- predict(m,temp)
            # # avoid extrapolating into positive error factor
            # error_factor[error_factor>1] <- 1
            # ggplot() + geom_point(aes(x=RD+AD, y=AD/RD)) + scale_y_log10() +
            #     geom_point(aes(x=RD+AD, y=error_factor),col='red')
            
            all <- NULL
            current <- temp[ix]
            ix_g <- current$type=='A/T > C/G'
            ix_b <- current$type=='C/G > A/T'
            ix <- !ix_g & !ix_b
            i=0
            for (other in 90:100) {
                current[ix,error_factor:=other]
                current[ix,ratio:=AD/(AD+RD*error_factor/100)]
                for (green in 90:100) {
                    current[ix_g,error_factor:=green]
                    current[ix_g,ratio:=AD/(AD+RD*error_factor/100)]
                    for (blue in 90:100) {
                        current[ix_b,error_factor:=blue]
                        current[ix_b,ratio:=AD/(AD+RD*error_factor/100)]
                        
                        current[,mapd:=median(abs(diff(abs(ratio-.5))))]
                        i <- i+1
                        all[[i]] <- copy(current)
                    }
                } 
            } 
            bound <- rbindlist(all)
            best <- bound[mapd==min(mapd)]
            #plot(best$ratio)
            best <- unique(best[,.(type,error_factor)])
            
            
            all <- NULL
            current <- temp[AD>10 & RD>10]
            ix_g <- current$type=='A/T > C/G'
            ix_b <- current$type=='C/G > A/T'
            ix <- !ix_g & !ix_b
            i=0
            best_other <- best[type=='C/G > C/G']$error_factor[1]
            best_green <- best[type=='A/T > C/G']$error_factor[1]
            best_blue <- best[type=='C/G > A/T']$error_factor[1]
            for (other in -5:5) {
                current[ix,error_factor:=best_other+other/10]
                current[ix,ratio:=AD/(AD+RD*error_factor/100)]
                for (green in -5:5) {
                    current[ix_g,error_factor:=best_green+green/10]
                    current[ix_g,ratio:=AD/(AD+RD*error_factor/100)]
                    for (blue in -5:5) {
                        current[ix_b,error_factor:=best_blue+blue/10]
                        current[ix_b,ratio:=AD/(AD+RD*error_factor/100)]
                        
                        current[,mapd:=median(abs(diff(abs(ratio-.5))))]
                        i <- i+1
                        all[[i]] <- copy(current)
                    }
                } 
            } 
            bound <- rbindlist(all)
            best <- bound[mapd==min(mapd)]
            
            # ggplot(best) + geom_point(aes(x=1:nrow(best), y=AD/(AD+RD),fill=type),shape=21)
            # ggplot(best) + geom_point(aes(x=1:nrow(best), y=AD/(AD+RD*error_factor/100),fill=type),shape=21)
            best <- unique(best[,.(type,error_factor)])
            
            
            
            #for (t in unique(temp$type)) temp[type==t]$error_factor <- median(temp[ix & type==t]$ratio,na.rm=T)
            for (t in unique(temp$type)) temp[type==t]$error_factor <- best[type==t]$error_factor[1]
            
            #AD_corrected <- AD/temp$error_factor
            allele_ratio <- unname(round(AD/(AD+RD*temp$error_factor/100),4))
            
        }
        
        #ggplot() + geom_point(aes(x=RD+AD, y=allele_ratio))
        
        # if main sample is this sample
        if (j==t_ix) snp_table$allele_ratio <- allele_ratio
        
        # if main sample is not this sample, and there is a normal sample
        if (j!=t_ix) if (normal) snp_table$allele_ratio_n <- allele_ratio
    }
    
    
    overlaps <- findOverlaps(target_ranges,rowRanges(vcf))
    snp_table[subjectHits(overlaps)]$target <- queryHits(overlaps)
    
    #best <- min_DP > 100 & max_AD > 10 & max_RD > 10 # 300, 50, 50
    #overlaps <- findOverlaps(target_ranges,rowRanges(vcf[best]))

    
    targets[queryHits(overlaps)]$allele_ratio <-  snp_table[subjectHits(overlaps)]$allele_ratio
    
    
    if (normal) targets[queryHits(overlaps)]$allele_ratio_n <- snp_table[subjectHits(overlaps)]$allele_ratio_n
    
}



# Segmentation ------------------------------------------------------------



getsegs <- function(targets, logratio) {
    segments <- segmentByCBS(y=logratio,avg='median',
                             chromosome=as.numeric(str_replace(targets$chromosome,'X','23')),
                             alpha = 0.3,undo=1)
    segments <- as.data.table(segments)[!is.na(chromosome),-1]
    segments[,start_pos:=targets$start[ceiling(start)]]
    segments[,end_pos:=targets$end[floor(end)]]
    segments[,genes:=''] 
    for (i in 1:nrow(segments)) {
        ix <- ceiling(segments[i]$start):floor(segments[i]$end)
        #suppressWarnings(segments[i]$mean <- round(mean(targets[ix]$log2,trim=.1,na.rm=T)),3)
        genes <- unique(targets$gene[ix])
        genes <- genes[!genes %in% c('','Background')]
        if (length(genes)>0) segments[i]$genes <- paste(genes,collapse = ', ')
    }
    segments[,chromosome:=str_replace(as.character(chromosome),'23','X')]
    return(segments)
}


# Merge targets and background for segmentation
# both <- rbind(targets[,.(chromosome,start,mid,end,log2,type='t')],
#               background[,.(chromosome,start,mid,end,log2,type='b')])

both <- rbind(targets,background,fill=T)
both[,chromosome:=str_replace(chromosome,'Y','24')][,chromosome:=str_replace(chromosome,'X','23')][,chromosome:=as.numeric(chromosome)]

both <- both[order(chromosome,mid)][,bin:=1:.N]

#segments <- getsegs(targets, targets$log2)
segments <- getsegs(both, both$log2)

if (normal) {
    
    # calculate values for normal sample and tumor segments
    segments$mean_n <- as.numeric(NA)
    for (i in 1:nrow(segments)) {
        ix <- ceiling(segments[i]$start):floor(segments[i]$end)
        segments[i]$mean_n <- round(median(both[ix]$log2_n,na.rm=T),3)
    }
    
    # segment the normal...
    segments_n <- getsegs(both, both$log2_n)
    
}


# Gene/segment tables ------------------------------------------------------------



segments_ <- segments[,.(segment=paste(1:.N),type='segment',chromosome,start=start_pos,end=end_pos,
                         length=end_pos-start_pos,
                         bins=nbrOfLoci,genes,mean)]
genes <- targets[,.(segment='',type='gene',chromosome,start,end,length=NA,bins=0,
                    genes=gene,log2,
                    mean=0)]
if (normal) {
    segments_$mean_n <- segments$mean_n
    genes$log2_n <- targets$log2_n
    genes$mean_n <- 0
} 
for (i in 1:nrow(segments)) {
    ix <- ceiling(segments[i]$start):floor(segments[i]$end)
    genes[ix]$segment <- i
}
genes <- genes[genes!='']
genes[,segment:=paste(unique(segment),collapse = ','),by=genes]
genes[,start:=min(start),by=genes]
genes[,end:=max(end),by=genes]
genes[,length:=end-start]
genes[,bins:=.N,by=genes]
suppressWarnings(genes[,mean:=round(median(log2,na.rm=T),3),by=genes])
if (normal) suppressWarnings(genes[,mean_n:=round(median(log2_n,na.rm=T),3),by=genes])
genes[,log2:=NULL]
if (normal) genes[,log2_n:=NULL]

suppressWarnings(
    segments_tumor <- rbind(segments_,unique(genes[bins>20]))[order(as.numeric(chromosome),start)]
    )




if (normal) {
    
    segments_ <- segments_n[,.(segment=paste(1:.N),type='segment',chromosome,start=start_pos,end=end_pos,
                             length=end_pos-start_pos,
                             bins=nbrOfLoci,genes,mean)]
    genes <- targets[,.(segment='',type='gene',chromosome,start,end,length=NA,bins=0,
                        genes=gene,log2=log2_n,
                        mean=0)]
    for (i in 1:nrow(segments_)) {
        ix <- ceiling(segments[i]$start):floor(segments[i]$end)
        genes[ix]$segment <- i
    }
    genes <- genes[genes!='']
    genes[,segment:=paste(unique(segment),collapse = ','),by=genes]
    genes[,start:=min(start),by=genes]
    genes[,end:=max(end),by=genes]
    genes[,length:=end-start]
    genes[,bins:=.N,by=genes]
    suppressWarnings(genes[,mean:=round(median(log2,na.rm=T),3),by=genes])
    #if (normal) suppressWarnings(genes[,mean_n:=round(median(log2_n,na.rm=T),3),by=genes])
    genes[,log2:=NULL]
    #if (normal) genes[,log2_n:=NULL]
    
    suppressWarnings(
        segments_normal <- rbind(segments_,unique(genes[bins>20]))[order(as.numeric(chromosome),start)]
    )
    

}


# Table output ------------------------------------------------------------


clinbarcode <- str_remove(name, "_nodups.bam")

fwrite(x = segments_tumor,file = paste0(opt$output_dir,'/',clinbarcode,'.segments.csv'))
fwrite(x = both,file = paste0(opt$output_dir,'/',clinbarcode,'.bins.csv'))
fwrite(x = targets,file = paste0(opt$output_dir,'/',clinbarcode,'.targets.csv'))
# fwrite(x = background,file = paste0(opt$output_dir,'/',name,'.background.csv'))


# for compatibility with CNVkit. Currently with targets only to allow PureCN to work.
# cnr:  chromosome      start   end     gene    depth   log2    weight
cnr <- targets[,.(chromosome=as.character(chromosome),start,end,gene,depth=round(count/width*160,3),log2,weight=1)][gene=='',gene:='-']
cnr[,chromosome:=str_replace(chromosome,'23','X')][,chromosome:=str_replace(chromosome,'24','Y')]
fwrite(x = cnr,file = paste0(opt$output_dir,'/',clinbarcode,'.cnr'),sep = '\t')

# cns:  chromosome      start   end     gene    log2    depth   probes  weight
cns <- segments_tumor[type=='segment',.(chromosome,start,end,gene=genes,log2=mean,depth=mean,probes=bins)]
fwrite(x = cns,file = paste0(opt$output_dir,'/',clinbarcode,'.cns'),sep = '\t')

if (normal) {
    n_clinbarcode <- str_remove(name_n, "_nodups.bam")
    fwrite(x = segments_normal,file = paste0(opt$output_dir,'/',n_clinbarcode,'.segments.csv'))
    #fwrite(x = background,file = paste0(opt$output_dir,'/',name,'-',name_n,'.background.csv'))

    # for compatibility
    # cnr:  chromosome  start   end     gene    depth   log2    weight
    cnr <- both[,.(chromosome=as.character(chromosome),start,end,gene,depth=round(count_n/width*160,3),log2=log2_n,weight=1)][gene=='',gene:='-']
    cnr[,chromosome:=str_replace(chromosome,'23','X')][,chromosome:=str_replace(chromosome,'24','Y')]
    fwrite(x = cnr,file = paste0(opt$output_dir,'/',n_clinbarcode,'.cnr'),sep = '\t')

    # cns:  chromosome  start   end     gene    log2    depth   probes  weight
    cns <- segments_normal[type=='segment',.(chromosome,start,end,gene=genes,log2=mean,depth=mean,probes=bins)]
    fwrite(x = cns,file = paste0(opt$output_dir,'/',n_clinbarcode,'.cns'),sep = '\t')
}

# DNAcopy segment file for main sample only:
# ID    chrom   loc.start       loc.end num.mark        seg.mean        C
seg <- segments[,.(ID=name,chrom=chromosome,loc.start=start_pos,loc.end=end_pos,num.mark=nbrOfLoci,seg.mean=mean,C=NA)]
fwrite(x = seg,file = paste0(opt$output_dir,'/',clinbarcode,'_dnacopy.seg'),sep = '\t')


# Count file output ------------------------------------------------------------

saveRDS(counts,paste0(opt$output_dir,'/',clinbarcode,'.counts.RDS'))
if (normal) saveRDS(counts_normal,paste0(opt$output_dir,'/',n_clinbarcode,'.counts.RDS'))


# Plot 1: overview ------------------------------------------------------------

both$type <- 'Target'
both[gene=='Background']$type <- 'Background'


pdf(file = paste0(opt$output_dir,'/',clinbarcode,'.pdf'),width = 20,height=14,pointsize = 1)

stats <- paste('Coverage:', 
               paste(round(quantile(targets[is_backbone==T]$count,c(.01,.99))),collapse = '-')
               )

# make chroms object for plots
chroms=data.table(chromosome=unique(both$chromosome),
                  start=0,
                  end=0,
                  mid=0)
for (chr in chroms$chromosome) {
    chroms$start[chr==chroms$chromosome]=both[chromosome==chr,min(bin)]
    chroms$end[chr==chroms$chromosome]=both[chromosome==chr,max(bin)]
    chroms$mid[chr==chroms$chromosome]=both[chromosome==chr,mean(bin)]
}


p1 <- ggplot(both[count/width*160 > 100]) + 
    geom_point(mapping = aes(x=bin,y=count/width*160,fill=label),
               shape=21,col='#00000050',size=1) +
    facet_wrap(facets = vars(sample),ncol = 1) + ylab('coverage') +
    scale_fill_hue() + scale_y_log10() + 
    guides(fill=guide_legend(override.aes=list(shape=21,size=2))) +
    scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                       expand = c(.01,.01),labels = chroms$chromosome) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_line(),
          panel.grid.minor.x = element_line(color = 'black'),
          axis.line = element_line(),
          axis.ticks = element_line())

p2 <- ggplot(both) + 
    geom_point(mapping = aes(x=bin,y=log2,fill=label),shape=21,col='#00000050',size=1) +
    facet_wrap(facets = vars(sample),ncol = 1) +
    scale_fill_hue() + ylim(c(-3,3)) +
    geom_segment(data=segments,col='green',size=1,
                 mapping = aes(x=start,xend=end,y=mean,yend=mean)) +
    guides(fill=guide_legend(override.aes=list(shape=21,size=2))) +
    scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                       expand = c(.01,.01),labels = chroms$chromosome) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_line(),
          panel.grid.minor.x = element_line(color = 'black'),
          axis.line = element_line(),
          axis.ticks = element_line())

if (snp_allele_ratio) p3 <- ggplot(both) + 
    geom_point(mapping = aes(x=bin,y=allele_ratio,fill=label),
               shape=21,col='#00000050',size=1) +
    facet_wrap(facets = vars(sample),ncol = 1) +
    scale_fill_hue() + ylim(0:1) +
    guides(fill=guide_legend(override.aes=list(shape=21,size=2))) +
    scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                       expand = c(.01,.01),labels = chroms$chromosome) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_line(),
          panel.grid.minor.x = element_line(color = 'black'),
          axis.line = element_line(),
          axis.ticks = element_line())

q1 <- ggplot(targets[count/width*160 > 100],mapping = aes(x=gc,y=count/width*160)) +
    xlim(c(0.1,0.85)) + ylab('Coverage') +
    geom_point(aes(fill=label),shape=21,col='#00000050',size=1) + 
    scale_fill_hue() + scale_y_log10() + 
    facet_wrap(facets = vars(label),ncol=3) +
    guides(fill=guide_legend(override.aes=list(shape=21,size=2))) +
    theme(panel.spacing = unit(0, "lines"))

q2 <- ggplot(targets,mapping = aes(x=gc,y=log2)) +
    xlim(c(0.1,0.85)) + ylim(c(-3,3)) +
    geom_point(aes(fill=label),shape=21,col='#00000050',size=1) + 
    scale_fill_hue() +
    facet_wrap(facets = vars(label),ncol=3) +
    guides(fill=guide_legend(override.aes=list(shape=21,size=2))) +
    theme(panel.spacing = unit(0, "lines"))


temp <- targets[!is.na(label),median(log2),by=label]
targets[,smooth_log2:=runmed(log2,k=9),by=chromosome]
if (snp_allele_ratio) q3 <- ggplot() +
    geom_point(data=targets,aes(x=smooth_log2,y=allele_ratio),col='#00000050',shape=21,fill='grey',size=1,show.legend = F) + 
    geom_point(data=targets[!is.na(label)],aes(x=smooth_log2,y=allele_ratio,fill=label),shape=21,col='#00000050',size=1,show.legend=F) + 
    geom_point(mapping=aes(x=temp$V1,y=0,fill=temp$label),size=2,shape=24,show.legend=F) +
    scale_fill_hue() +
    scale_x_continuous(limits = c(min(temp$V1-.5,-2),max(temp$V1+.5,2))) +
    scale_y_continuous(limits = c(0,1)) + #,breaks = c(seq(0,1,.1)),minor_breaks = NULL,labels = c(seq(0,1,.1))) +
    guides(fill=guide_legend(override.aes=list(shape=21,size=2))) +
    theme(panel.spacing = unit(0, "lines"))



if (snp_allele_ratio) {
    pa <- plot_annotation(
        title = paste(name,date(), stats),
    )
    layout <- "AAAAADD\nBBBBBEE\nCCCCCFF"
    fig=p1+p2+p3+q1+q2+q3+
        plot_layout(design = layout,guides = 'collect')
} else {
    pa <- plot_annotation(
        title = paste(name,date(), stats),
    )
    layout <- "AAAAACC\nBBBBBDD"
    fig=p1+p2+q1+q2+
        plot_layout(design = layout,guides = 'collect')
}

print(fig+pa)



# Plot 2: T/N overview ------------------------------------------------------------


if (normal) {
    p2_n <- ggplot(both) + 
        geom_point(mapping = aes(x=bin,y=log2_n,fill=label),shape=21,col='#00000050',size=1) +
        facet_wrap(facets = vars(sample_n),ncol = 1) +
        scale_fill_hue() + ylim(c(-3,3)) +
        geom_segment(data=segments_n,col='green',size=1,
                     mapping = aes(x=start,xend=end,y=mean,yend=mean)) +
        guides(fill=guide_legend(override.aes=list(shape=21,size=2))) +
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line())
    
    if (snp_allele_ratio) p3_n <- ggplot(both) + 
        geom_point(mapping = aes(x=bin,y=allele_ratio_n,fill=label),
                   shape=21,col='#00000050',size=1) +
        facet_wrap(facets = vars(sample_n),ncol = 1) +
        scale_fill_hue() + ylim(0:1) +
        guides(fill=guide_legend(override.aes=list(shape=21,size=2))) +
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line())
    
    q2_n <- ggplot(targets,mapping = aes(x=gc,y=log2_n)) +
        xlim(c(0.1,0.85)) + ylim(c(-3,3)) +
        geom_point(aes(fill=label),shape=21,col='#00000050',size=1) + 
        scale_fill_hue() +
        facet_wrap(facets = vars(label),ncol=3) +
        guides(fill=guide_legend(override.aes=list(shape=21,size=2))) +
        theme(panel.spacing = unit(0, "lines"))
    
    temp_n <- targets[!is.na(label),median(log2_n),by=label]
    targets[,smooth_log2_n:=runmed(log2_n,k=9),by=chromosome]
    if (snp_allele_ratio) q3_n <- ggplot() +
        geom_point(data=targets,aes(x=smooth_log2_n,y=allele_ratio_n),col='#00000050',shape=21,fill='grey',size=1,show.legend = F) + 
        geom_point(data=targets[!is.na(label)],aes(x=smooth_log2_n,y=allele_ratio_n,fill=label),shape=21,col='#00000050',size=1,show.legend=F) + 
        geom_point(mapping=aes(x=temp_n$V1,y=0,fill=temp_n$label),size=2,shape=24,show.legend=F) +
        scale_fill_hue() +
        scale_x_continuous(limits = c(min(temp$V1-.5,-2),max(temp$V1+.5,2))) +
        scale_y_continuous(limits = c(0,1)) + #,breaks = c(seq(0,1,.1)),minor_breaks = NULL,labels = c(seq(0,1,.1))) +
        guides(fill=guide_legend(override.aes=list(shape=21,size=2))) +
        theme(panel.spacing = unit(0, "lines"))
    
    if (snp_allele_ratio) {
        layout <- "AAAAAAEE\nBBBBBBFF\nCCCCCCGG\nDDDDDDHH"
        fig=p2+p3+p2_n+p3_n+q2+q3+q2_n+q3_n+
            plot_layout(design = layout,guides = 'collect')
        print(fig+pa)
    } else {
        layout <- "AAAACC\nBBBBDD"
        fig=p2+p2_n+q2+q2_n+
            plot_layout(design = layout,guides = 'collect')
        print(fig+pa)
    }
}




# Plot 3: fragment size effect ------------------------------------------------------------

if (F) {
    short <- both
    short[,smooth_log2:=runmed(log2_short,21),by='chromosome']
    short[,smooth_rawLR:=runmed(rawLR_short,21),by='chromosome']
    short$fragments <- '<160'
    all <- both
    all[,smooth_log2:=runmed(log2,21),by='chromosome']
    all[,smooth_rawLR:=runmed(rawLR,21),by='chromosome']
    all$fragments <- 'all'
    
    p1_ <- ggplot(both) + 
        geom_point(mapping = aes(x=bin,y=rawLR,fill=label),shape=21,col='#00000050',size=2) +
        geom_line(data=rbind(all,short),mapping = aes(x=bin,y=smooth_rawLR,col=fragments),size=2) +
        facet_wrap(facets = vars(sample),ncol = 1) + ylim(c(-3,3)) +
        scale_fill_hue() + scale_color_manual(values = c('all'='black','<160'='lightgrey')) + 
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line())
    
    
    p2_ <- ggplot(both) + ylab('<160 log2') +
        geom_point(mapping = aes(x=bin,y=log2_short,fill=label),shape=21,col='#00000050',size=1) +
        geom_line(data=rbind(all,short),mapping = aes(x=bin,y=smooth_log2,col=fragments),size=1.5) +
        facet_wrap(facets = vars(sample),ncol = 1) + ylim(c(-3,3)) +
        scale_fill_hue() + scale_color_manual(values = c('all'='black','<160'='lightgrey')) + 
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line())
    
    q1_ <- ggplot(targets,mapping = aes(x=rawLR,y=rawLR_short)) +
        xlab('All fragments') + ylab('Fragments < 160 bases') +
        geom_point(aes(fill=label),shape=21,col='#00000050',size=1) + 
        scale_fill_hue() + geom_abline() +
        facet_wrap(facets = vars(label),ncol=3) +
        guides(fill=guide_legend(override.aes=list(shape=21,size=2))) +
        theme(panel.spacing = unit(0, "lines"))
    
    q2_ <- ggplot(targets,mapping = aes(x=log2,y=log2_short)) +
        xlab('All fragments') + ylab('Fragments < 160 bases') +
        geom_point(aes(fill=label),shape=21,col='#00000050',size=1) + 
        scale_fill_hue() + geom_abline() +
        facet_wrap(facets = vars(label),ncol=3) +
        guides(fill=guide_legend(override.aes=list(shape=21,size=2))) +
        theme(panel.spacing = unit(0, "lines"))
    
    pa=plot_annotation(
        title = paste(name,date()),
    )
    # layout <- "AAAAACC\nBBBBBDD"
    # fig=p1+p2+q1+q2+
    #     plot_layout(design = layout,guides = 'collect')
    # 
    # print(fig+pa)
    
    if (snp_allele_ratio) {
        pa <- plot_annotation(
            title = paste(name,date(), stats),
        )
        layout <- "AAAAAEE\nBBBBBFF\nCCCCCGG\nDDDDDHH"
        fig=p1+p2+p2_+p3+q1+q2+q2_+q3+
            plot_layout(design = layout,guides = 'collect')
    } else {
        pa <- plot_annotation(
            title = paste(name,date(), stats),
        )
        layout <- "AAAAADD\nBBBBBEE\nCCCCCFF"
        fig=p1+p2+p2_+q1+q2+q2_+
            plot_layout(design = layout,guides = 'collect')
    }
    
    print(fig+pa)
}

# Close pdf ------------------------------------------------------------
dev.off()
