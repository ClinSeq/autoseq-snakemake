#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(StructuralVariantAnnotation))
suppressPackageStartupMessages(library(stringr))
options(warn=-1)

simpleEventType <- function(gr) {
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "TRA", # inter-chromosomosal
    ifelse(strand(gr) == strand(pgr), "INV",
      ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
        ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
          "DUP")))))
}

# Options ------------------------------------------------------------
option_list <- list(
    make_option(c("-v", "--vcf"), action = "store", type = "character",default = NULL,
                help = "GRIDSS sv vcf file"),
    make_option(c("-o", "--output_vcf"), action = "store", type = "character",default = '.',
                help = "Annotated output vcf file")
)

opt <- parse_args(OptionParser(option_list = option_list))

# using the example in the GRIDSS /example directory
vcf <- readVcf(opt$vcf, "hg19")
info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
	row.names=c("SIMPLE_TYPE"),
	Number=c("1"),
	Type=c("String"),
	Description=c("Simple event type annotation based purely on breakend position and orientation."))), "DataFrame"))

gr <- breakpointRanges(vcf)
svtype <- simpleEventType(gr)
info(vcf)$SIMPLE_TYPE <- NA_character_
info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype
# info(vcf[gr$sourceId])$SVLEN <- gr$svLen
writeVcf(vcf, opt$output_vcf)
