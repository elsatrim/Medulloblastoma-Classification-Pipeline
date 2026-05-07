#!/usr/bin/env Rscript
# ichorCNA Wrapper - Fixed WIG format (no scientific notation)

suppressPackageStartupMessages({
  library(optparse)
  library(ichorCNA)
  library(HMMcopy)
  library(Rsamtools)
  library(GenomicRanges)
})

option_list <- list(
  make_option(c("--id"), type="character", default=NULL, help="Sample ID"),
  make_option(c("--bam"), type="character", default=NULL, help="Path to BAM file"),
  make_option(c("--outdir"), type="character", default=".", help="Output directory"),
  make_option(c("--threads"), type="integer", default=1, help="Number of threads")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$id) || is.null(opt$bam)) {
  print_help(opt_parser)
  stop("Both --id and --bam required", call.=FALSE)
}

if (!file.exists(opt$bam)) {
  stop(paste("BAM not found:", opt$bam), call.=FALSE)
}

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

cat("=== Running ichorCNA for MB:", opt$id, "===\n")

ichorCNA_path <- system.file(package="ichorCNA")
gc_wig <- file.path(ichorCNA_path, "extdata", "gc_hg38_1000kb.wig")
map_wig <- file.path(ichorCNA_path, "extdata", "map_hg38_1000kb.wig")
centromere_file <- file.path(ichorCNA_path, "extdata", "GRCh38.GCA_000001405.2_centromere_acen.txt")

if (!file.exists(gc_wig)) {
  stop("ichorCNA reference files not found", call.=FALSE)
}

cat("[1/4] Detecting chromosomes...\n")

bam_header <- scanBamHeader(opt$bam)
chr_names <- names(bam_header[[1]]$targets)
has_chr_prefix <- any(grepl("^chr", chr_names))

if (has_chr_prefix) {
  chrs_to_use <- paste0("chr", c(1:22, "X"))
} else {
  chrs_to_use <- c(1:22, "X")
}

chrs_available <- intersect(chrs_to_use, chr_names)
cat("  Found", length(chrs_available), "chromosomes\n")

cat("[2/4] Generating WIG file (fixedStep format)...\n")

wig_file <- file.path(opt$outdir, paste0(opt$id, ".wig"))

bin_size <- 1000000
chr_lengths <- bam_header[[1]]$targets[chrs_available]

bins <- GRanges()
for (chr in names(chr_lengths)) {
  starts <- seq(1, chr_lengths[chr], by = bin_size)
  ends <- pmin(starts + bin_size - 1, chr_lengths[chr])
  bins <- c(bins, GRanges(seqnames = chr, ranges = IRanges(start = starts, end = ends)))
}

cat("  Counting reads in", length(bins), "bins (10-20 min)...\n")

param <- ScanBamParam(
  which = bins,
  flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE, isNotPassingQualityControls = FALSE),
  mapqFilter = 20
)

counts <- countBam(opt$bam, param = param)

wig_data <- data.frame(
  chr = as.character(seqnames(bins)),
  start = start(bins),
  count = counts$records,
  stringsAsFactors = FALSE
)

if (!has_chr_prefix) {
  wig_data$chr <- paste0("chr", wig_data$chr)
}

cat("  Writing fixedStep WIG format (no scientific notation)...\n")

wig_conn <- file(wig_file, "w")

for (chr in unique(wig_data$chr)) {
  chr_data <- wig_data[wig_data$chr == chr, ]
  chr_data <- chr_data[order(chr_data$start), ]
  
  # CRITICAL: Use format() to avoid scientific notation
  # Write as plain integers
  header <- sprintf("fixedStep chrom=%s start=1 step=%d span=%d", 
                    chr, 
                    as.integer(bin_size), 
                    as.integer(bin_size))
  
  writeLines(header, wig_conn)
  writeLines(as.character(chr_data$count), wig_conn)
}

close(wig_conn)

total_reads <- sum(wig_data$count)
cat("  WIG created:", wig_file, "\n")
cat("  Total reads:", total_reads, "\n")
cat("  Estimated coverage:", round(total_reads * 150 / 3e9, 2), "x\n")

cat("[3/4] Running ichorCNA (30-60 min)...\n")

ichorCNA_script <- file.path(ichorCNA_path, "scripts", "runIchorCNA.R")

cmd <- paste(
  "Rscript", ichorCNA_script,
  "--WIG", wig_file,
  "--gcWig", gc_wig,
  "--mapWig", map_wig,
  "--centromere", centromere_file,
  "--id", opt$id,
  "--outDir", opt$outdir,
  "--includeHOMD FALSE",
  "--chrs 'c(1:22, \"X\")'",
  "--chrTrain 'c(1:22)'",
  "--estimateNormal TRUE",
  "--estimatePloidy TRUE",
  "--maxCN 5",
  "--scStates 'NULL'",
  "--txnE 0.9999",
  "--txnStrength 10000",
  "--genomeStyle UCSC",
  "--genomeBuild hg38"
)

result <- system(cmd, intern = FALSE)

if (result != 0) {
  stop("ichorCNA failed with exit code ", result, call.=FALSE)
}

cat("[4/4] Validating...\n")

params_file <- file.path(opt$outdir, paste0(opt$id, ".params.txt"))
seg_file <- file.path(opt$outdir, paste0(opt$id, ".seg.txt"))

if (!file.exists(params_file)) {
  params_file <- file.path(opt$outdir, opt$id, paste0(opt$id, ".params.txt"))
  seg_file <- file.path(opt$outdir, opt$id, paste0(opt$id, ".seg.txt"))
}

if (!file.exists(params_file)) {
  stop("ichorCNA output not found in:", opt$outdir, call.=FALSE)
}

params <- read.table(params_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

cat("\n=== ichorCNA Results ===\n")
cat("Tumor Fraction:", params$Tumor_Fraction[1], "\n")
cat("Ploidy:", params$Ploidy[1], "\n")

tf <- as.numeric(params$Tumor_Fraction[1])
if (!is.na(tf)) {
  if (tf < 0.08) {
    cat("\nWARNING: Low tumor fraction (", round(tf*100, 1), "%)\n", sep="")
  } else {
    cat("\nTumor fraction adequate\n")
  }
}

cat("\nOutput files:\n")
cat("  Parameters:", params_file, "\n")
cat("  Segments:", seg_file, "\n")

cat("\n=== Complete! ===\n")