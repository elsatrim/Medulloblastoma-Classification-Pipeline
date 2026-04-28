#!/usr/bin/env Rscript
# Run ichorCNA for MB CNV detection
# Based on Liu 2021 parameters

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option(c("-i", "--id"), type="character", help="Sample ID"),
  make_option(c("-b", "--bam"), type="character", help="Input BAM"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("--ichorcna"), type="character", help="ichorCNA path"),
  make_option(c("-t", "--threads"), type="integer", default=4)
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("=== Running ichorCNA for MB:", opt$id, "===\n")

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(opt$outdir, "wig"), showWarnings=FALSE)

# Step 1: Create WIG file
cat("[1/3] Generating WIG file...\n")
wig_file <- file.path(opt$outdir, "wig", paste0(opt$id, ".wig"))

system(paste(
  "readCounter",
  "--window 1000000",
  "--quality 20",
  "--chromosome '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X'",
  opt$bam,
  ">", wig_file
))

# Step 2: Run ichorCNA
# Parameters optimized for MB (based on Liu 2021)
cat("[2/3] Running ichorCNA...\n")

ichorcna_script <- file.path(opt$ichorcna, "scripts", "runIchorCNA.R")
gc_wig <- file.path(opt$ichorcna, "inst", "extdata", "gc_hg38_1000kb.wig")
map_wig <- file.path(opt$ichorcna, "inst", "extdata", "map_hg38_1000kb.wig")
centromere <- file.path(opt$ichorcna, "inst", "extdata", "GRCh38.GCA_000001405.2_centromere_acen.txt")

# MB-specific parameters (Liu 2021)
# - Pediatric tumors often near-diploid
# - CSF samples can have low tumor fraction
# - Need to detect focal amplifications (MYC, MYCN)

ichorcna_cmd <- paste(
  "Rscript", ichorcna_script,
  "--id", opt$id,
  "--WIG", wig_file,
  "--ploidy 'c(2,3)'",  # MB typically diploid or triploid
  "--normal 'c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)'",  # Wide range for CSF
  "--maxCN 7",  # Higher to detect amplifications
  "--gcWig", gc_wig,
  "--mapWig", map_wig,
  "--centromere", centromere,
  "--includeHOMD FALSE",
  "--chrs 'c(1:22, \"X\")'",
  "--chrTrain 'c(1:22)'",
  "--estimateNormal TRUE",
  "--estimatePloidy TRUE",
  "--estimateScPrevalence TRUE",
  "--scStates 'c(1,3)'",
  "--txnE 0.9999",
  "--txnStrength 10000",
  "--fracReadsInChrYForMale 0.001",
  "--minMapScore 0.75",
  "--outDir", opt$outdir
)

system(ichorcna_cmd)

# Step 3: Parse and validate results
cat("[3/3] Validating results...\n")

params_file <- file.path(opt$outdir, paste0(opt$id, ".params.txt"))
seg_file <- file.path(opt$outdir, paste0(opt$id, ".seg.txt"))

if (file.exists(params_file) && file.exists(seg_file)) {
  params <- fread(params_file)
  segs <- fread(seg_file)
  
  tumor_fraction <- params$Tumor_Fraction[1]
  ploidy <- params$Ploidy[1]
  
  cat("\n=== ichorCNA Results ===\n")
  cat("Tumor Fraction:", round(tumor_fraction * 100, 2), "%\n")
  cat("Ploidy:", ploidy, "\n")
  cat("Segments detected:", nrow(segs), "\n")
  
  # QC based on Markowitz 2025
  pass_qc <- tumor_fraction >= 0.08
  
  cat("QC Status:", ifelse(pass_qc, "PASS", "LOW TUMOR FRACTION"), "\n")
  
  if (!pass_qc) {
    cat("\nWARNING: Tumor fraction below 8%\n")
    cat("Classification confidence may be reduced.\n")
  }
  
  # Save summary
  summary_data <- data.frame(
    Sample_ID = opt$id,
    Tumor_Fraction = tumor_fraction,
    Ploidy = ploidy,
    N_Segments = nrow(segs),
    Pass_QC = pass_qc,
    Timestamp = Sys.time()
  )
  
  fwrite(summary_data, 
         file.path(opt$outdir, paste0(opt$id, "_ichorCNA_summary.csv")))
  
} else {
  cat("\nERROR: ichorCNA failed!\n")
  quit(status=1)
}

cat("\n=== ichorCNA complete ===\n")
