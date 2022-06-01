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
#             snp_vcf='sample1n.vcf.gz',
#             output_dir='.'
# )


# #with bam
# setwd('~/Data/CLINSEQ/')
# opt <- list(reference_file='comprehensive3_baits_twist.bed.reference.RDS',
#             input_bam='sample1.bam',
#             snp_vcf='sample1.vcf',
#             output_dir='.'
# )


# #PSFF with RDS and vcf
# setwd('~/Data/CLINSEQ/')
# opt <- list(reference_file='references/pancancer2_baits_twist.bed.reference.RDS',
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



# Tables of bins ------------------------------------------------------------


targets <- reference$targets

background <- reference$background
background$gene='Background'

name <- str_remove(opt$input_bam,'.*/')
name <- str_remove(name,'\\.counts.RDS')

targets$sample <- name
background$sample <- name


# Add counts
targets[,count:=counts$target_count]
#sum_by_bg <- targets[,sum(count),by=background]

background[,count:=counts$background_count]
#background[,count:=raw_count][sum_by_bg$background,count:=raw_count-sum_by_bg$V1]
#background[count<0,count:=0]

# same, short
targets[,count_short:=counts$target_short]
#sum_by_bg <- targets[,sum(count_short),by=background]
background[,count_short:=counts$background_short]
#background[,count_short:=raw_short][sum_by_bg$background,count_short:=raw_short-sum_by_bg$V1]
#background[count_short<0,count_short:=0]


# keep only bins "ok" in this reference set
targets <- targets[reference$keep_targets]
background <- background[reference$keep_background]


target_ranges <- makeGRangesFromDataFrame(targets)
background_ranges <- makeGRangesFromDataFrame(background)



# SNP allele ratio ------------------------------------------------------------
save.image('ws.Rdata')

snp_allele_ratio <- FALSE
input <- opt$snp_vcf
if (!is.null(input)) {
    snp_allele_ratio <- TRUE
    if (!str_detect(input,'.[vV][cC][fF]$') & !str_detect(input,'.[vV][cC][fF].[gG][zZ]$')) stop("SNP vcf file appears incorrect")
    
    vcf=readVcf(input)
    
    # remove SNPs that do not have 2 alleles
    alleles <- as.data.table(table(as.data.table(alt(vcf))$group))$N
    vcf <- vcf[alleles==1]
    
    # If there is more than one sample in the VCF, which one to use?
    ix <- 1 # base assumption: first is this sample
    names <- colnames(vcf)
    if (length(names)>1) {
        # but if one name in the VCF fits the sample name better, use it:
        check <- which(str_detect(name,names))
        if (length(check)==1) ix <- check
    }
    
    # Keep only this sample
    vcf <- vcf[,ix]
    
    # Prepare table of SNPs for this sample
    snp_table <- data.table(snp=1:length(vcf),target_row=as.integer(NA),id=names(vcf))
    
    # Get alleles
    snp_table$ref_allele <- as.character(ref(vcf))
    alt_allele <- as.data.table(alt(vcf))[, .(values = list(value)), by = group][,values]
    snp_table$n_alt_alleles <- sapply(alt_allele, length)
    snp_table$alt_allele <- sapply(alt_allele, "[[", 1)
    
    # Get SNP type
    from <- rep('C/G',length(vcf)); from[snp_table$ref_allele %in% c('A','T')] <- 'A/T'
    to <- rep('C/G',length(vcf)); to[snp_table$alt_allele %in% c('A','T')] <- 'A/T'
    snp_table$type <- paste(from,'>',to)
    snp_table[str_length(ref_allele)!=1 | str_length(alt_allele)!=1]$type <- 'other'
    
    # Get read counts
    g <- geno(vcf)
    snp_table$AD <- sapply(g$AD, "[[", 2)
    snp_table$RD <- sapply(g$AD, "[[", 1)
    snp_table$DP=snp_table$AD+snp_table$RD
    
    # Compute allele ratio
    raw_allele_ratio <- unname(round(snp_table$AD/snp_table$DP,4))
    raw_allele_ratio[is.nan(raw_allele_ratio)] <- 0
    snp_table$allele_ratio <- raw_allele_ratio
    snp_table[allele_ratio==0,type:='none']
    
    # Assign to correct targets
    overlaps <- findOverlaps(target_ranges,rowRanges(vcf))
    snp_table[subjectHits(overlaps)]$target_row <- queryHits(overlaps)

    # compute error factors by SNP type (in development, and paused...)
    if (FALSE) {
        
        # Selection to use for alt/ref correction model
        ix <- AD>50 & RD > 50 & raw_allele_ratio > .4 & raw_allele_ratio < .6 & DP < quantile(DP,.95)*1.1
        
        # ggplot() + geom_point(aes(x=RD[ix]+AD[ix], y=AD[ix]/RD[ix],fill=type[ix]),shape=21) +
        #     scale_y_log10() +
        #     geom_smooth(aes(x=RD[ix]+AD[ix], y=AD[ix]/RD[ix],group=type[ix],col=type[ix]))
        # ggplot() + geom_point(aes(x=1:sum(ix), y=AD[ix]/RD[ix],fill=type[ix]),shape=21)
        
        
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
    
    
}



# LogR and genotype correction ------------------------------------------------------------


peakx <- function(data) {
    d <- density(data)
    maxy <- max(d$y)
    maxx <- d$x[d$y==maxy][1]
    return(maxx)
} 
mapd <- function(data) {
    return(median(abs(diff(data))))
} 

targets$snp <- 'none'
targets$allele_ratio <- 0


# set raw log ratio (overwritten below if SNPs available)
targets[,rawLR:=log2(count+1)]
targets[,rawLR_short:=log2(count_short+1)]


# put allele ratio in target table, keeping the highest AF where multiple SNPs map to same target
if (snp_allele_ratio) { 
    temp <- snp_table[type!='other'][order(allele_ratio)]
    
    targets[temp$target_row]$snp <-  temp$type
    targets[temp$target_row]$allele_ratio <-  temp$allele_ratio
    
    
    
    
    
    
    # plot(density(targets[snp!='other' & is_backbone & gene=='' & allele_ratio==0]$count))
    # points(density(targets[snp!='other' & is_backbone & gene=='' & allele_ratio==1]$count),col='red')
    # ggplot(snp_table[type!='other' ]) + geom_point(aes(x=DP,y=allele_ratio,col=type)) + facet_wrap(facets = vars(type)) + ylim(c(.4,.6))
    # ggplot(snp_table[type!='other' ]) + geom_point(aes(x=target_row,y=allele_ratio,col=type)) + facet_wrap(facets = vars(type))
    # ggplot(targets) + geom_point(aes(x=count,y=allele_ratio,col=snp)) + facet_wrap(facets = vars(snp))
    # ggplot(targets[is_backbone==T]) + geom_density(aes(x=count,col=snp)) + facet_wrap(facets = vars(snp))
    # ggplot(targets[T==is_backbone]) + geom_histogram(aes(x=count,col=snp),bins=100) + facet_wrap(facets = vars(snp))
    
    # compute SNP allele bias
    targets[,ref_bias:=0]
    targets[,ref_bias:=peakx(count[snp!='other' & is_backbone & gene=='' & allele_ratio<.01]/
                                 peakx(count[snp!='other' & is_backbone & gene=='' & allele_ratio>.99]))-1]
    
    # bias should not be negative and not too high
    targets[ref_bias<0,ref_bias:=0]
    targets[ref_bias>.10,ref_bias:=.05]
    
    
    
    
    # SNP corrected logR, targets
    targets[,rawLR:=log2(count+ref_bias*allele_ratio*count+1)]
    targets[,rawLR_short:=log2(count_short+ref_bias*allele_ratio*count_short+1)]
    
    
}

# median correct to backbone
targets[,rawLR:=rawLR-median(rawLR[is_backbone])]
targets[,rawLR_short:=rawLR_short-median(rawLR_short[is_backbone])]

# median correct by target
targets[,rawLR:=rawLR-reference$targets_median]
targets[,rawLR_short:=rawLR_short-reference$targets_median_short]


# background
background[,rawLR:=log2(count+1)]
background[,rawLR_short:=log2(count_short+1)]

background[,rawLR:=rawLR-median(rawLR[is_backbone])]
background[,rawLR_short:=rawLR_short-median(rawLR_short[is_backbone])]

background[,rawLR:=rawLR-reference$background_median]
background[,rawLR_short:=rawLR_short-reference$background_median_short]


# Reference data correction ------------------------------------------------------------

# correct using reference
jcorrect <- function(temp,train_ix=NULL,gc=NULL) {
    
    if (is.null(train_ix)) train_ix <- rep(TRUE,nrow(temp))
    
    loess_temp=loess(lr ~ PC1+PC2+PC3, data = temp,
                     subset = train_ix,
                     family="symmetric", control = loess.control(surface = "direct"))
    temp[,lr:=lr-predict(loess_temp,temp)]
    
    # if (!is.null(gc)) { 
    # loess_temp=loess(lr ~ gc, data = temp,
    #                  subset = train_ix,
    #                  family="symmetric", control = loess.control(surface = "direct"))
    # temp[,lr:=lr-predict(loess_temp,temp)]
    # }
    
    runs <- 0
    best_lr <- temp$lr
    best_mapd <- mapd(temp$lr)
    for (i in 1:(ncol(temp)-2)) {
        
        temp$thispc=temp[[paste0('PC',i)]]
        loess_temp=rlm(lr ~ thispc, data=temp,
                       subset = train_ix)
        temp[,lr:=lr-predict(loess_temp,temp)]
        
        mapd <- mapd(temp$lr)
        cat(round(mapd,4),'>>')
        if (mapd < best_mapd) {
            best_lr <- temp$lr
            best_mapd <- mapd
            runs <- i
        }
        
    }
    cat('ran ',runs,'/',ncol(temp)-3, ' times\n')
    return(best_lr)
}

# standard
temp <- cbind(data.table(lr=targets$rawLR),reference$targets_ref,gc=targets$gc)
targets[,log2:=jcorrect(temp,targets$is_backbone)]

# short
temp <- cbind(data.table(lr=targets$rawLR_short),reference$targets_ref_short,gc=targets$gc)
targets[,log2_short:=jcorrect(temp,targets$is_backbone)]

# standard bg
temp <- cbind(data.table(lr=background$rawLR),reference$background_ref)
background[,log2:=jcorrect(temp,is_backbone)]

# bg short
temp <- cbind(data.table(lr=background$rawLR_short),reference$background_ref_short)
background[,log2_short:=jcorrect(temp,is_backbone)]



# ggplot(background) + geom_point(aes(x=background,y=log2)) +
#     geom_point(data=targets,mapping = aes(x=background,y=log2),col='blue')



# Segmentation ------------------------------------------------------------



getsegs <- function(targets, logratio) {
    segments <- segmentByCBS(y=logratio,avg='median',
                             chromosome=as.numeric(str_replace(targets$chromosome,'X','23')),
                             alpha = 0.01,undo=1)
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
# bins <- rbind(targets[,.(chromosome,start,mid,end,log2,type='t')],
#               background[,.(chromosome,start,mid,end,log2,type='b')])

bins <- rbind(targets,background,fill=T)
bins[,chromosome:=str_replace(chromosome,'Y','24')][,chromosome:=str_replace(chromosome,'X','23')][,chromosome:=as.numeric(chromosome)]

bins <- bins[order(chromosome,mid)][,bin:=1:.N]

segments <- getsegs(bins, bins$log2)



# Gene/segment tables ------------------------------------------------------------



segments_temp <- segments[,.(segment=paste(1:.N),type='segment',chromosome,start=start_pos,end=end_pos,
                         length=end_pos-start_pos,
                         bins=nbrOfLoci,genes,mean)]
genes <- targets[,.(segment='',type='gene',chromosome,start,end,length=NA,bins=0,
                    genes=gene,log2,
                    mean=0)]


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

genes[,log2:=NULL]


suppressWarnings(
    segments_genes <- rbind(segments_temp,unique(genes[bins>20]))[order(as.numeric(chromosome),start)]
)




# Table output ------------------------------------------------------------


clinbarcode <- str_remove(name, "_nodups.bam")
bins[,chromosome:=as.character(chromosome)][chromosome=='23',chromosome:='X'][chromosome=='24',chromosome:='Y']

fwrite(x = segments_genes,file = paste0(opt$output_dir,'/',clinbarcode,'.segments.csv'))
fwrite(x = bins,file = paste0(opt$output_dir,'/',clinbarcode,'.bins.csv'))
fwrite(x = targets,file = paste0(opt$output_dir,'/',clinbarcode,'.targets.csv'))
# fwrite(x = background,file = paste0(opt$output_dir,'/',name,'.background.csv'))


# for compatibility with CNVkit. Currently with targets only to allow PureCN to work (does not allow overlaps)
# cnr:  chromosome      start   end     gene    depth   log2    weight
cnr <- targets[,.(chromosome=as.character(chromosome),start,end,gene,depth=round(count/width*160,3),log2,weight=1)][gene=='',gene:='-']
cnr[,chromosome:=str_replace(chromosome,'23','X')][,chromosome:=str_replace(chromosome,'24','Y')]
fwrite(x = cnr,file = paste0(opt$output_dir,'/',clinbarcode,'.cnr'),sep = '\t')

# cns:  chromosome      start   end     gene    log2    depth   probes  weight
cns <- segments_genes[type=='segment',.(chromosome,start,end,gene=genes,log2=mean,depth=mean,probes=bins)]
fwrite(x = cns,file = paste0(opt$output_dir,'/',clinbarcode,'.cns'),sep = '\t')


# DNAcopy segment file:
# ID    chrom   loc.start       loc.end num.mark        seg.mean        C
seg <- segments[,.(ID=name,chrom=chromosome,loc.start=start_pos,loc.end=end_pos,num.mark=nbrOfLoci,seg.mean=mean,C=NA)]
fwrite(x = seg,file = paste0(opt$output_dir,'/',clinbarcode,'_dnacopy.seg'),sep = '\t')


# Count file output ------------------------------------------------------------
# (not overwrite)
if (!file.exists(paste0(opt$output_dir,'/',clinbarcode,'.counts.RDS')))
    saveRDS(counts,paste0(opt$output_dir,'/',clinbarcode,'.counts.RDS'))


# Save workspace ------------------------------------------------------------
#save.image(paste0(opt$output_dir,'/',clinbarcode,'.jumble_workspace.RDS'))


bins$type <- 'Target'
bins[gene=='Background']$type <- 'Background'




# Jumble frankenplot ------------------------------------------------------------

noise <- function(data) {
    m <- mapd(data)
    f <- 2^m-1
    return(round(100*f,2))
}

if (T) {
    
    
    p <- NULL
    t <- bins[chromosome==13 & str_detect(gene,'RB1')]$target; bins[target %in% min(t):max(t),gene:='RB1']
    bins[,label:=NA][gene %in% c('AR','ATM','BRCA2','PTEN','RB1','NTRK3','ERG','CDK12','TMPRSS2'),label:=gene]
    bins[,smooth_log2:=runmed(log2,k=21),by=chromosome]
    ylims <- c(.4,2) #c(min(.5,min(2^bins$smooth_log2)),max(2,max(2^bins$smooth_log2)))
    #if (ylims[1]<.25) ylims[1] <- .25
    #if (ylims[2]>4) ylims[1] <- 4

    if (snp_allele_ratio) { 
        bins[snp=='other',allele_ratio:=NA]
        #bins[allele_ratio>.95 | allele_ratio<.05,allele_ratio:=NA]
        bins[,maf:=abs(allele_ratio-.5)+.5]
        bins[!is.na(maf)][maf<.95,maf:=runmed(maf,9)]
        # snp (grid) smooth-to-alleleratio plot
        p$grid <- ggplot(bins[maf<.95]) + xlim(c(.2,1.8)) + ylim(c(.5,1)) + xlab('Corrected depth (smooth)') + ylab('Major allele ratio (smooth)') +
            geom_point(data=bins[maf<.95,.(log2,maf)],aes(x=2^log2,y=maf),col='lightgrey') +
            geom_point(aes(x=2^log2,y=maf),fill='#60606090',col='#20202090',shape=21) +
            geom_point(data=bins[maf<.95 & label!=''],aes(x=2^log2,y=maf,fill=label),shape=21,col='#00000050',size=1) +
            facet_wrap(facets = vars(factor(chromosome,levels=unique(chromosome),ordered=T)),ncol = 8) +
            theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 8))
        # snp (all) smooth-to-alleleratio plot
        temp <- bins[!is.na(label),median(log2),by=label]
        p$nogrid <- ggplot(bins[maf<.95]) + xlim(c(0,2.5)) + ylim(c(.5,1)) + xlab('Corrected depth (smooth)') + ylab('Major allele ratio (smooth)') +
            geom_point(data=bins[maf<.95,.(log2,maf)],aes(x=2^log2,y=maf),fill='#60606090',col='#20202090',shape=21) +
            geom_point(data=bins[label!=''&maf<.95],aes(x=2^log2,y=abs(allele_ratio-.5)+.5,fill=label),shape=21,col='#00000050',size=1) +
            geom_point(data=temp,mapping=aes(x=2^V1,y=1,fill=label),size=2,shape=25,show.legend=F)
    }
    # chroms object by genomic pos
    chroms <- data.table(chromosome=names(reference$chromlength),length=reference$chromlength)
    chroms[,start:=as.double(0)] 
    chroms[,stop:=as.double(length)] 
    chroms[,mid:=as.double(round(length/2))] 
    for (i in 2:nrow(chroms)) {
        chroms[i,start:=chroms$stop[i-1]]
        chroms[i,stop:=chroms$stop[i-1]+length]
        chroms[i,mid:=chroms$stop[i-1]+round(length/2)]
    }    
    # fix positions by genome
    bins[,gpos:=mid]
    for (chr in unique(bins$chromosome)[-1]) bins[chromosome==chr,gpos:=gpos+sum(chroms[1:(which(chromosome==chr)-1)]$length)]
    segments[,gstart:=as.double(start_pos)][,gstop:=as.double(end_pos)]
    for (chr in unique(bins$chromosome)[-1]) {
        segments[chromosome==chr,gstart:=gstart+sum(chroms[1:(which(chromosome==chr)-1)]$length)]
        segments[chromosome==chr,gstop:=gstop+sum(chroms[1:(which(chromosome==chr)-1)]$length)]
    }
    # logR by pos + segments (2nd left)
    p$pos_log2 <- ggplot(bins) + xlab('Genomic position') + ylab('Corrected depth') +
        geom_point(data=bins[is.na(label)],mapping = aes(x=gpos,y=2^log2),fill='#60606070',col='#20202070',shape=21,size=1) +
        geom_point(data=bins[!is.na(label)],mapping = aes(x=gpos,y=2^log2,fill=label),shape=21,col='#00000050',size=1) +
        scale_fill_hue() + scale_y_log10(limits=ylims) +
        geom_segment(data=segments,col='green',size=1,
                     mapping = aes(x=gstart,xend=gstop,y=2^mean,yend=2^mean)) +
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line()) 
    
    if (snp_allele_ratio) {
        # allele ratio by pos
        p$pos_alleleratio <- ggplot(bins) + xlab('Genomic position') + ylab('Allele ratio') +
            geom_point(data=bins[is.na(label)],mapping = aes(x=gpos,y=allele_ratio),fill='#60606080',col='#20202080',shape=21,size=1) +
            geom_point(data=bins[!is.na(label)],mapping = aes(x=gpos,y=allele_ratio,fill=label),shape=21,col='#00000050',size=1) +
            scale_fill_hue() + ylim(0:1) +
            scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                               expand = c(.01,.01),labels = chroms$chromosome) +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_line(),
                  panel.grid.minor.x = element_line(color = 'black'),
                  axis.line = element_line(),
                  axis.ticks = element_line()) 
    }
    
    # chroms object by order
    chroms=data.table(chromosome=unique(bins$chromosome),
                      start=0,
                      end=0,
                      mid=0)
    for (chr in chroms$chromosome) {
        chroms$start[chr==chroms$chromosome]=bins[chromosome==chr,min(bin)]
        chroms$end[chr==chroms$chromosome]=bins[chromosome==chr,max(bin)]
        chroms$mid[chr==chroms$chromosome]=bins[chromosome==chr,mean(bin)]
    }
    # logR by order
    p$order_log2 <- ggplot(bins) + xlab('Order of genomic position') + ylab('Corrected depth') +
        geom_point(data=bins[is.na(label)],mapping = aes(x=bin,y=2^log2),fill='#60606070',col='#20202070',shape=21,size=1) +
        geom_point(data=bins[!is.na(label)],mapping = aes(x=bin,y=2^log2,fill=label),shape=21,col='#00000050',size=1) +
        scale_fill_hue() + scale_y_log10(limits=ylims) +
        geom_segment(data=segments,col='green',size=1,
                     mapping = aes(x=start,xend=end,y=2^mean,yend=2^mean)) +
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line()) 
    # logR by gc
    p$gc_log2 <- ggplot(bins) + xlab('Target GC content') + ylab('Corrected depth') +
        geom_point(data=bins,mapping = aes(x=gc,y=2^log2),fill='#60606040',col='#20202040',shape=21,size=1) + # fill='#60606040'
        geom_smooth(data=bins[!is.na(label)],mapping = aes(x=gc,y=2^log2,col=label),size=.5,se=F,show.legend = F) +
        scale_fill_hue() + scale_y_log10(limits=ylims) 
    m <- bins[!is.na(target),median(count*160/width)]
    if (snp_allele_ratio) {
        # allele ratio by order 
        p$order_alleleratio <- ggplot(bins) + xlab('Order of genomic position') + ylab('Allele ratio') +
            geom_point(data=bins[is.na(label)],mapping = aes(x=bin,y=allele_ratio),fill='#60606050',col='#20202050',shape=21,size=1) +
            geom_point(data=bins[!is.na(label)],mapping = aes(x=bin,y=allele_ratio,fill=label),shape=21,col='#00000050',size=1) +
            scale_fill_hue() + ylim(0:1) +
            scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                               expand = c(.01,.01),labels = chroms$chromosome) +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_line(),
                  panel.grid.minor.x = element_line(color = 'black'),
                  axis.line = element_line(),
                  axis.ticks = element_line()) 
        # allele ratio by depth
        p$depth_alleleratio <- ggplot(bins) + xlab('Sequence depth') + ylab('Allele ratio') +
            geom_point(data=bins[is.na(label)],mapping = aes(x=count,y=allele_ratio),fill='#60606070',col='#20202070',shape=21,size=1) +
            geom_point(data=bins[!is.na(label)],mapping = aes(x=count,y=allele_ratio,fill=label),shape=21,col='#00000050',size=1) +
            scale_fill_hue() + ylim(0:1) + scale_x_log10(limits=c(m/3,m*3))
    }
    # depth by order 
    p$order_rawdepth <- ggplot(bins) + xlab('Order of genomic position') + ylab('Sequence depth') +
        geom_point(data=bins[is.na(background)],mapping = aes(x=bin,y=count),fill='#60606050',col='#20202050',shape=21,size=1) +
        geom_point(data=bins[!is.na(label)],mapping = aes(x=bin,y=count,fill=label),shape=21,col='#00000050',size=1) +
        scale_fill_hue() + scale_y_log10(limits=c(750,3000)) +
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line()) 
    # depth by GC 
    p$gc_rawdepth <- ggplot(bins) + xlab('Target GC content') + ylab('Sequence depth') +
        geom_point(data=bins,mapping = aes(x=gc,y=count),fill='#60606040',col='#20202040',shape=21,size=1) + # 
        geom_smooth(data=bins[!is.na(label)],mapping = aes(x=gc,y=count,col=label),size=.5,se=F,show.legend = F) +
        scale_fill_hue() + scale_y_log10(limits=c(750,3000))
    
    
    for (i in 1:length(p)) p[[i]] <- p[[i]] + guides(fill=guide_legend(override.aes=list(shape=21,size=3)))
    
    stats <- paste0('Coverage: ',
                    paste(round(quantile(targets[is_backbone==T]$count,c(.01,.99))),collapse = '-'),
                    ', Noise: ',
                    noise(targets$log2),'% / ', noise(background$log2),'%'
    )
    
    pa <- plot_annotation(
        title = paste(clinbarcode,'         ',date(),'         ',stats),
    )

    if (snp_allele_ratio) {
        
        layout <-  "ABBBB
                CDDDD
                EFFFF
                GGGGG
                HHHHH
                IJJJJ
                IJJJJ"
        fig <- 
            p$gc_rawdepth+p$order_rawdepth+
            p$gc_log2+p$order_log2+
            p$depth_alleleratio+p$order_alleleratio+
            p$pos_log2+
            p$pos_alleleratio+
            p$nogrid+p$grid+
            plot_layout(design = layout,guides = 'collect')
    } else {

        layout <-  "ABBBB
                CDDDD
                EEEEE
                "
        fig <- 
            p$gc_rawdepth+p$order_rawdepth+
            p$gc_log2+p$order_log2+
            p$pos_log2+
            plot_layout(design = layout,guides = 'collect')
    }
    
    png(file = paste0(opt$output_dir,'/',clinbarcode,'.png'),width = 1800,height=1400,res=100)
    print(fig+pa)
    
}


# Close image

dev.off()


# Save image object ------------------------------------------------------------

p$bins <- bins
p$segments <- segments
saveRDS(p,paste0(opt$output_dir,'/',clinbarcode,'.plots.RDS'))


