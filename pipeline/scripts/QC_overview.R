#!/usr/bin/env Rscript

# This is a script creating plots showing the distribution of some QC values, 
# and where the current normal and tumor sample lies within those distributions.
# QC metrics shown are read count, duplication rate, coverage, fold enrichment, on-bait rate, insert size & contamination


library(ggplot2)
library(optparse)
library(data.table)
library(rjson)

# read in command line options
option_list <- list(
  make_option(c("-s", "--samples"), action = "store", type = "character",
              default = NULL, help = "Samples of interest to match with qc file names. If more than one they must be separated with colon, e.g. tumor:normal"),
  make_option(c("-d", "--analysisdir"), action = "store", type = "character",
              default = NULL, help = "Analysis directory for samples of interest, in order to separate between the same samples that possibly are being run in multiple analyses"),
  make_option(c("-o", "--outfile"), action = "store", type = "character",
              default = "QC_overview.pdf", help = "Path to output pdf file"),
  make_option(c("-m", "--mainpath"), action = "store", type = "character",
              default = "/nfs/PROBIO/autoseq-output", help = "Path to output pdf file")
)

opt <- parse_args(OptionParser(option_list = option_list))
sample_string = opt$samples
analysis_dir = opt$analysisdir
outfile = opt$outfile
main_path = opt$mainpath


# find the qc files for all samples
cat("Find all available qc files...\n")
hsmetrics_files = dir(Sys.glob(paste0(main_path, "/*/*/qc/picard")), 
                pattern = "picard-hsmetrics.txt$", recursive = TRUE, full.names = TRUE)
markduplicates_files = dir(Sys.glob(paste0(main_path, "/*/*/qc/picard")), 
                           pattern = "markdup.metrics.txt$", recursive = TRUE, full.names = TRUE)
insertsize_files = dir(Sys.glob(paste0(main_path, "/*/*/qc/picard")), 
                           pattern = "picard-insertsize.txt$", recursive = TRUE, full.names = TRUE)
contest_files = dir(Sys.glob(paste0(main_path, "/*/*/contamination")), 
                    pattern = "contest.txt$", recursive = TRUE, full.names = TRUE)
msings_files = dir(path = Sys.glob(paste0(main_path, "/*/*/msings*")),
                   pattern = paste0("MSI_Analysis.txt$"), full.names = TRUE, recursive = TRUE)
flagstat_files = dir(path = Sys.glob(paste0(main_path, "/*/*/qc/samtools")),
                   pattern = "flagstats.json$", full.names = TRUE, recursive = TRUE)

# remove files from old pipeline (modified before 2019-04-30)
hsmetrics_files = hsmetrics_files[file.mtime(hsmetrics_files) > as.POSIXct("2019-04-30")]
markduplicates_files = markduplicates_files[file.mtime(markduplicates_files) > as.POSIXct("2019-04-30")]
insertsize_files = insertsize_files[file.mtime(insertsize_files) > as.POSIXct("2019-04-30")]
contest_files = contest_files[file.mtime(contest_files) > as.POSIXct("2019-04-30")]
msings_files = msings_files[file.mtime(msings_files) > as.POSIXct("2019-04-30")]

# remove files from WGS samples (will otherwise break the script due to different columns etc)
hsmetrics_files = hsmetrics_files[grep("WGS", hsmetrics_files, invert = TRUE)]
markduplicates_files = markduplicates_files[grep("WGS", markduplicates_files, invert = TRUE)]
insertsize_files = insertsize_files[grep("WGS", insertsize_files, invert = TRUE)]
contest_files = contest_files[grep("WGS", contest_files, invert = TRUE)]
msings_files = msings_files[grep("WGS", msings_files, invert = TRUE)]


# read in the files
cat("Read in the files...\n")
HsMetrics = data.frame()
for (f in hsmetrics_files) {
  tryCatch({
    SAMP = strsplit(basename(f), split = "\\.")[[1]][1]
    DIR = dirname(dirname(dirname(f)))
    HsMetrics = rbind(HsMetrics, cbind(SAMP, DIR, read.table(f, skip = 6, nrow = 1, sep = "\t", 
                                                        header = TRUE, stringsAsFactors = FALSE), 
                                                        stringsAsFactors = FALSE))
  }, error = function(err) {
      print(paste("Sample: ", SAMP, " QC: HsMetrics"))
      print(paste("ERROR: ", err))
  })
}

MarkDuplicates = data.frame()
for (f in markduplicates_files) {
    tryCatch({
        SAMP = sub("-picard-markdup.metrics.txt", "", basename(f))
        DIR = dirname(dirname(dirname(f)))
        MarkDuplicates = rbind(MarkDuplicates, cbind(SAMP, DIR, read.table(f, skip = 6, nrow = 1, sep = "\t", 
                                                                    header = TRUE, stringsAsFactors = FALSE), 
                                                                    stringsAsFactors = FALSE))
    }, error = function(err) {
      print(paste("Sample: ", SAMP, " QC: MarkDuplicates"))
      print(paste("ERROR: ", err))
    }) 
}

InsertSize = data.frame()
InsertSize_histogram = data.frame(stringsAsFactors = FALSE)
for (f in insertsize_files) {
    tryCatch({
        SAMP = strsplit(basename(f), split = "\\.")[[1]][1]
        DIR = dirname(dirname(dirname(f)))
        InsertSize = rbind(InsertSize, cbind(SAMP, DIR, read.table(f, skip = 6, nrow = 1, sep = "\t", 
                                                             header = TRUE, stringsAsFactors = FALSE), 
                                       stringsAsFactors = FALSE))
        InsertSize_histogram = rbind(InsertSize_histogram, 
                               cbind(SAMP, DIR, read.table(f, skip = 10, sep = "\t",
                                                           header = TRUE, stringsAsFactors = FALSE),
                                     stringsAsFactors = FALSE))
    }, error = function(err) {
      print(paste("Sample: ", SAMP, " QC: InsertSize"))
      print(paste("ERROR: ", err))
    }) 
}

ContEst = data.frame()
for (f in contest_files) {
  tryCatch({
    fname = strsplit(basename(f), split = "\\.")[[1]][1]
    clinseq_barcodes = rev(strsplit(f, split = "/")[[1]])[3]
    samples = unlist(strsplit(clinseq_barcodes, split = "_"))
    if(grepl('-N-', fname, fixed=TRUE)){
      SAMP = samples[1]
    } else {
      SAMP = samples[2]
    }
    DIR = dirname(dirname(f))
    ContEst = rbind(ContEst, cbind(SAMP, DIR, read.table(f, nrow = 1, sep = "\t", header = TRUE, stringsAsFactors = FALSE), 
                                  stringsAsFactors = FALSE))
  }, error = function(err) {
    print(paste("Sample: ", SAMP, " QC: contamination test"))
    print(paste("ERROR: ", err))
  })    
}

msings = data.frame()
for (f in msings_files) {
  tryCatch({
    SAMP = sub("_nodups.MSI_Analysis.txt", "", basename(f))
    DIR = dirname(dirname(dirname(f)))
    msings = rbind(msings, cbind(SAMP, DIR, t(read.table(f, header = FALSE, nrows = 5, sep = "\t", stringsAsFactors = FALSE)))[2,], stringsAsFactors=FALSE)
  }, error = function(err) {
    print(paste("Sample: ", SAMP, " QC: mSINGs"))
    print(paste("ERROR: ", err))
  })
}

colnames(msings) = c("SAMP", "DIR", read.table(f, header = FALSE, nrows = 5, sep = "\t", stringsAsFactors = FALSE)[,1])
msings$`msi status`[which(msings$msing_score<0.2)] = "NEG"  # use cut-off 0.2 for MSI-H (default 0.1 is too low)
msings$msing_score = as.numeric(msings$msing_score)

# samtools flagstats
flagstat_data = data.frame(matrix(ncol = 8, nrow = 0))
colnames(flagstat_data) = c("SAMP", "DIR", "mapped_reads", "paired_reads", "properly_paired_reads", 
                       "with_itself_and_mate_mapped_reads", "singletons_reads", "mate_mapped_diff_chr_reads")

for (f in flagstat_files){
  json_data = fromJSON(file = f)
  read_total = as.numeric(json_data[['QC-passed reads']]['total'] )
  SAMP = sub("-flagstats.json", "", basename(f))
  DIR = dirname(dirname(dirname(f)))
  mapped = as.numeric(json_data[['QC-passed reads']]['mapped']) / read_total 
  paired = as.numeric(json_data[['QC-passed reads']]['paired in sequencing']) / read_total 
  properly_paired = as.numeric(json_data[['QC-passed reads']]['properly paired']) / read_total
  with_itself_mate_mapped = as.numeric(json_data[['QC-passed reads']]['with itself and mate mapped']) / read_total 
  matemapped_diff_chr = as.numeric(json_data[['QC-passed reads']]['with mate mapped to a different chr']) / read_total 
  singletons = as.numeric(json_data[['QC-passed reads']]['singletons']) / read_total
  nrows = nrow(flagstat_data)
  flagstat_data[nrows+1, ] = c(SAMP, DIR, mapped, paired, properly_paired, 
                          with_itself_mate_mapped, singletons, matemapped_diff_chr)
}


# add missing values in the insert size histogram table
for (d in unique(InsertSize_histogram$DIR)) {
  d1 = subset(InsertSize_histogram, DIR == d)
  for (s in unique(d1$SAMP)) {
    d2 = subset(d1, SAMP == s)
    missing = which(!seq(2,max(d2$insert_size)) %in% d2$insert_size) + 1
    if (length(missing)>0) {
      InsertSize_histogram = rbind(InsertSize_histogram,
                                   data.frame(SAMP = s, DIR = d, insert_size = missing, All_Reads.fr_count = 0))
    }
    # print(paste0("Sample ", s, " missing sizes: ", paste0(missing, collapse = ", ")))
  }
}


# add normalized counts for insert size histogram
InsertSize_histogram$count_norm = NA
for (d in unique(InsertSize_histogram$DIR)) {
  for (s in unique(subset(InsertSize_histogram, DIR == d)$SAMP)) {
    tot_reads = MarkDuplicates$READ_PAIRS_EXAMINED[which(MarkDuplicates$DIR == d & MarkDuplicates$SAMP == s)]
    idx = which(InsertSize_histogram$DIR == d & InsertSize_histogram$SAMP == s)
    InsertSize_histogram$count_norm[idx] = InsertSize_histogram$All_Reads.fr_count[idx]/tot_reads
  }
}


process_samp <- function(sample) {
    pre = unlist(strsplit(sample, split = "-"))
    tmp = paste(pre[1:5], collapse="-")
    libkit = substr(pre[6], 1, 2)
    capture = substr(pre[7], 1, 2)
    return(paste(c(tmp, libkit, capture), collapse="-"))
}


# merge the QC tables 
qc_merge = merge(merge(merge(merge(merge(HsMetrics, MarkDuplicates, by = c("SAMP", "DIR")), 
                InsertSize, by = c("SAMP", "DIR")), ContEst, by = c("SAMP", "DIR")),
                msings, by = c("SAMP", "DIR"), all.x = TRUE), flagstat_data, by = c("SAMP", "DIR"), all.x = TRUE)

qc_merge$mapped_reads = as.numeric(qc_merge$mapped_reads)
qc_merge$paired_reads = as.numeric(qc_merge$paired_reads)
qc_merge$properly_paired_reads = as.numeric(qc_merge$properly_paired_reads)
qc_merge$with_itself_and_mate_mapped_reads = as.numeric(qc_merge$with_itself_and_mate_mapped_reads)
qc_merge$singletons_reads = as.numeric(qc_merge$singletons_reads)
qc_merge$mate_mapped_diff_chr_reads = as.numeric(qc_merge$mate_mapped_diff_chr_reads)
# add sample type and capture kit
qc_merge$SAMP = sapply(qc_merge$SAMP, process_samp)
qc_merge$sample_type = sapply(strsplit(qc_merge$SAMP, split = "-"), "[", 4)
qc_merge$capture = sapply(strsplit(qc_merge$SAMP, split = "-"), "[", 7)
qc_merge$capture = substr(qc_merge$capture, 1, 2)
InsertSize_histogram$sample_type = sapply(strsplit(InsertSize_histogram$SAMP, split = "-"), "[", 4)
InsertSize_histogram$capture = sapply(strsplit(InsertSize_histogram$SAMP, split = "-"), "[", 7)
InsertSize_histogram$capture = substr(InsertSize_histogram$capture, 1, 2)
InsertSize_histogram$SAMP = sapply(InsertSize_histogram$SAMP, process_samp)

# label the samples of interest (soi) and analysis directory of interest (doi)
samples = strsplit(sample_string, split = ":")[[1]]
qc_merge$soi = qc_merge$SAMP %in% samples
qc_merge$doi = qc_merge$DIR == Sys.glob(analysis_dir)
InsertSize_histogram$soi = InsertSize_histogram$SAMP %in% samples
InsertSize_histogram$doi = InsertSize_histogram$DIR == Sys.glob(analysis_dir)


# create an ouput table for the samples of interest
soi_table = data.table(qc_merge)[i = soi&doi, 
                                 j =list(SAMP, MEAN_TARGET_COVERAGE, FOLD_ENRICHMENT, dedupped_on_bait_rate=ON_BAIT_BASES/PF_BASES_ALIGNED, FOLD_80_BASE_PENALTY,
                                         READ_PAIRS_EXAMINED, PERCENT_DUPLICATION, "contamination_%"=contamination, MEDIAN_INSERT_SIZE, msing_score)]
table_outfile = sub("pdf$", "txt", outfile)
write.table(x = soi_table, file = table_outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




##################################################################################
#   ____   __       ____   _______
#  | o  \ | |     /     \ |__  __|
#  | __/  | |    |  / \ |   | |
#  ||     | |___ | \ /  |   | |
#  ||     |____| \____ /    |_|
#################################################################################

# Center plot titles (must be specified as of ggplot v 2.2.0)
theme_update(plot.title = element_text(hjust = 0.5))
theme_update(plot.subtitle = element_text(hjust = 0.5))

# plotting function to create desired histograms 
my_barplot = function(x, ybreaks, x_string, title_string) {
  p = ggplot(qc_merge, aes_string(x = x)) +
    geom_bar(aes(group = capture, fill = capture), width = 0.5, alpha = 0.7, color = "black") +
    geom_bar(data = subset(qc_merge, soi), width = 0.5, fill = "blue", show.legend = FALSE) + # the sample of interest
    geom_bar(data = subset(qc_merge, soi&doi), width = 0.5, fill = "red", show.legend = FALSE) + # the sample of interest
    scale_fill_manual(values = c("antiquewhite", "aliceblue", "lightpink", "palegreen")) +
    scale_y_continuous(breaks = ybreaks) +
    scale_x_discrete(name = x_string) +
    facet_wrap(~sample_type, ncol = 1) +
    ggtitle(title_string, subtitle = "Red piles show present samples, blue piles show these samples if run earlier")
  print(p)
}

# plotting function to create desired scatter plots
my_scatter = function(x, y, xbreaks, ybreaks, x_string, y_string, title_string) {
  ggplot(qc_merge, aes(shape = capture)) +
    geom_point(aes_string(x = x, y = y), size = 3) +
    geom_point(data = subset(qc_merge, soi), aes_string(x = x, y = y), fill = "blue", size = 3, show.legend = FALSE) +
    geom_point(data = subset(qc_merge, soi&doi), aes_string(x = x, y = y), fill = "red", size = 3, show.legend = FALSE) +
    scale_alpha_manual(values = c(0.7, 1)) +
    scale_shape_manual(values = c(24, 25, 21, 22), guide = guide_legend(override.aes = list(fill = NA))) +
    scale_x_continuous(name = x_string, breaks = xbreaks) +
    scale_y_continuous(name = y_string, breaks = ybreaks) +
    facet_wrap(~sample_type, ncol = 1) +
    ggtitle(title_string, subtitle = "Red symbols show present samples, blue symbols show these samples if run earlier")
}

# Plot histograms and scatter plots
cat(paste0("Create plots (saved in ", outfile, ") ...\n"))
pdf(file = outfile, width=14)

# duplication vs read count scatter plot
my_scatter(x = "READ_PAIRS_EXAMINED", y = "PERCENT_DUPLICATION", xbreaks = seq(0, 1e12, 1e7), ybreaks = seq(0, 1, 0.1),
           x_string = "number of read pairs", y_string = "duplication rate", title_string = "Duplication rate vs Read count")

# coverage vs read count scatter plot
my_scatter(x = "READ_PAIRS_EXAMINED", y = "MEAN_TARGET_COVERAGE", xbreaks = seq(0, 1e12, 1e7), ybreaks = seq(0, 5000, 500),
           x_string = "number of read pairs", y_string = "mean target coverage", title_string = "Coverage vs Read count")

# duplication vs fold enrichment scatter plot
my_scatter(x = "FOLD_ENRICHMENT", y = "PERCENT_DUPLICATION", xbreaks = seq(0, 5000, 50), ybreaks = waiver(),
           x_string = "fold enrichment", y_string = "duplication rate", title_string = "Duplication rate vs Fold enrichment")

# duplication vs on-bait rate scatter plot
my_scatter(x = "ON_BAIT_BASES/PF_BASES_ALIGNED", y = "PERCENT_DUPLICATION", xbreaks = waiver(), ybreaks = waiver(),
           x_string = "dedupped on-bait rate", y_string = "duplication rate", title_string = "Duplication rate vs Dedupped on-bait rate")

# fold80 base penalty vs coverage scatter plot
my_scatter(x = "MEAN_TARGET_COVERAGE", y = "FOLD_80_BASE_PENALTY", xbreaks = seq(0, 5000, 500), ybreaks = waiver(),
           x_string = "mean target coverage", y_string = "fold 80 base penalty", title_string = "Fold 80 base penalty vs Coverage")

# insert size histogram
p = ggplot(InsertSize_histogram, aes(x = insert_size, y = count_norm*1e6, group = interaction (SAMP, DIR),
                                 linetype = capture)) +
  geom_line(color = "black", alpha = 0.5) +
  geom_line(data = subset(InsertSize_histogram, soi&!doi), color = "blue", show.legend = FALSE) +
  geom_line(data = subset(InsertSize_histogram, soi&doi), color = "red", show.legend = FALSE) +
  # scale_linetype(guide = guide_legend(override.aes = list(color = "black"))) +
  scale_x_continuous(name = "insert size", breaks = seq(0,2000,100)) +
  scale_y_continuous(name = "count per million reads in total for each sample") +
  facet_wrap(~sample_type, ncol = 1, scales = "free_y") +
  ggtitle("Insert size", subtitle = "Red lines show present samples, blue lines show these samples if run earlier.")
# add a line marking the length of one chromatosome if any cfDNA sample is present, since this is the typical insert size for cfDNA
if (any(grepl("CFDNA", InsertSize_histogram$SAMP))) {
  p = p + geom_vline(xintercept = 167, alpha = 0.8, linetype = 3) +
    geom_label(data=data.frame(sample_type="CFDNA"), aes(x = 167, y = Inf), label = "167 bases, length of one chromatosome", 
               inherit.aes = FALSE, vjust = "top", hjust = 0, nudge_x = 2)
}
print(p)

# contamination bar plot
my_barplot(x = "factor(contamination, levels = sort(unique(contamination)))", ybreaks = waiver(),
             x_string = "contamination, %", title_string = "Contamination")

# msings score vs read count scatter plot
my_scatter(x = "READ_PAIRS_EXAMINED", y = "msing_score", xbreaks = seq(0, 1e12, 1e7), ybreaks = waiver(),
             x_string = "number of read pairs", y_string = "mSINGS score", title_string = "mSINGS score vs Read count")

my_scatter(x = "mapped_reads", y= "paired_reads", xbreaks = waiver(), ybreaks = waiver(),
             x_string = "fraction of mapped reads", y_string = "fraction of paired reads", title_string = "Paired Reads vs Mapped reads")

my_scatter(x = "mapped_reads", y= "properly_paired_reads", xbreaks = waiver(), ybreaks = waiver(),
             x_string = "fraction of mapped reads", y_string = "fraction of properly paired reads", title_string = "Properly Paired Reads vs Mapped reads")

my_scatter(x = "mapped_reads", y= "with_itself_and_mate_mapped_reads", xbreaks = waiver(), ybreaks = waiver(),
             x_string = "fraction of mapped reads", y_string = "fraction of with itself and mate mapped reads", title_string = "WithItself and MateMappedReads vs Mapped reads")

my_scatter(x = "mapped_reads", y= "singletons_reads", xbreaks = waiver(), ybreaks = waiver(),
             x_string = "fraction of mapped reads", y_string = "fraction of singletons reads", title_string = "Singletons Reads vs Mapped reads")

my_scatter(x = "mapped_reads", y= "mate_mapped_diff_chr_reads", xbreaks = waiver(), ybreaks = waiver(),
             x_string = "fraction of mapped reads", y_string = "fraction of mate mapped diff chr reads", title_string = "Mate Mapped Diff chr Reads vs Mapped reads")

dev.off()

cat("Completed.\n")


