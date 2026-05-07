#!/usr/bin/env Rscript
# MB Subgroup Classification with Fragmentomics Support
# Based on Liu et al. 2021 (CNV) + Markowitz et al. 2025 (fragmentomics)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(yaml)
  library(jsonlite)
})

# Parse command line arguments
option_list <- list(
  make_option("--id", type="character", help="Sample ID"),
  make_option("--seg", type="character", help="Path to ichorCNA seg file"),
  make_option("--params", type="character", help="Path to ichorCNA params file"),
  make_option("--config", type="character", help="Path to CNV signatures config"),
  make_option("--fragmentomics", type="character", default="", 
              help="Path to fragmentomics metrics JSON (optional)"),
  make_option("--output", type="character", help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("=====================================\n")
cat("MB SUBGROUP CLASSIFICATION\n")
cat("=====================================\n")
cat(sprintf("Sample: %s\n", opt$id))

# Load CNV signatures configuration
cat("\nLoading CNV signatures...\n")
cnv_signatures <- read_yaml(opt$config)

# Load ichorCNA results
cat("Loading ichorCNA results...\n")
seg <- fread(opt$seg)

# Normalize column names (ichorCNA uses different naming)
if ("chrom" %in% names(seg)) {
  setnames(seg, old = "chrom", new = "Chromosome")
}
if ("start" %in% names(seg)) {
  setnames(seg, old = "start", new = "Start")
}
if ("end" %in% names(seg)) {
  setnames(seg, old = "end", new = "End")
}
if ("seg.median.logR" %in% names(seg)) {
  setnames(seg, old = "seg.median.logR", new = "Median_logR")
}

cat(sprintf("  Loaded %d segments\n", nrow(seg)))

# Read params file - ichorCNA format is key-value pairs
params_raw <- readLines(opt$params)
params_list <- list()
for (line in params_raw) {
  if (grepl(":", line) && !grepl("^#", line)) {
    parts <- strsplit(line, ":")[[1]]
    if (length(parts) == 2) {
      key <- trimws(parts[1])
      value <- trimws(parts[2])
      params_list[[key]] <- value
    }
  }
}

# Extract tumor fraction and ploidy
tumor_fraction <- as.numeric(params_list[["Tumor Fraction"]])
ploidy <- as.numeric(params_list[["Ploidy"]])

if (is.na(tumor_fraction) || is.na(ploidy)) {
  cat("ERROR: Could not extract tumor fraction or ploidy from params file\n")
  cat("Available parameters:\n")
  print(names(params_list))
  quit(status = 1)
}

cat(sprintf("  Tumor fraction: %.2f%%\n", tumor_fraction * 100))
cat(sprintf("  Ploidy: %.2f\n", ploidy))

# Define function to check CNV alterations
check_cnv <- function(seg_data, chrom, type, threshold = 0.3) {
  chrom_seg <- seg_data[Chromosome == chrom | Chromosome == paste0("chr", chrom)]
  if (nrow(chrom_seg) == 0) return(FALSE)
  
  if (type == "gain") {
    return(any(chrom_seg$Median_logR > threshold))
  } else if (type == "loss") {
    return(any(chrom_seg$Median_logR < -threshold))
  } else if (type == "amplification") {
    return(any(chrom_seg$Median_logR > 0.8))
  }
  return(FALSE)
}

# Check for isochromosome 17q
check_iso17q <- function(seg_data) {
  chr17_seg <- seg_data[Chromosome == "17" | Chromosome == "chr17"]
  if (nrow(chr17_seg) == 0) {
    cat("  [DEBUG] No chr17 segments found\n")
    return(FALSE)
  }
  
  cat(sprintf("  [DEBUG] Found %d chr17 segments\n", nrow(chr17_seg)))
  
  # Sort by position
  chr17_seg <- chr17_seg[order(Start)]
  
  # Chromosome 17 centromere is around 25 Mb (GRCh38)
  # p arm: 0-25 Mb, q arm: 25-83 Mb
  centromere <- 25000000
  
  p_arm <- chr17_seg[End < centromere]
  q_arm <- chr17_seg[Start > centromere]
  
  cat(sprintf("  [DEBUG] Chr17p segments: %d, Chr17q segments: %d\n", 
              nrow(p_arm), nrow(q_arm)))
  
  if (nrow(p_arm) > 0) {
    cat(sprintf("  [DEBUG] Chr17p logR range: %.2f to %.2f\n", 
                min(p_arm$Median_logR), max(p_arm$Median_logR)))
  }
  if (nrow(q_arm) > 0) {
    cat(sprintf("  [DEBUG] Chr17q logR range: %.2f to %.2f\n", 
                min(q_arm$Median_logR), max(q_arm$Median_logR)))
  }
  
  # Check for loss in p arm and gain in q arm
  p_loss <- nrow(p_arm) > 0 && any(p_arm$Median_logR < -0.2)
  q_gain <- nrow(q_arm) > 0 && any(q_arm$Median_logR > 0.2)
  
  cat(sprintf("  [DEBUG] p_loss=%s, q_gain=%s\n", p_loss, q_gain))
  
  return(p_loss && q_gain)
}

# Detect CNV alterations
cat("\nDetecting CNV alterations...\n")
monosomy_6 <- check_cnv(seg, "6", "loss")
chr9q_loss <- check_cnv(seg, "9", "loss")
myc_amp <- check_cnv(seg, "8", "amplification")
iso17q <- check_iso17q(seg)
chr8_loss <- check_cnv(seg, "8", "loss")

cat(sprintf("  Monosomy 6: %s\n", ifelse(monosomy_6, "DETECTED", "NOT DETECTED")))
cat(sprintf("  Isochromosome 17q: %s\n", ifelse(iso17q, "DETECTED", "NOT DETECTED")))
cat(sprintf("  Chromosome 8 loss: %s\n", ifelse(chr8_loss, "DETECTED", "NOT DETECTED")))

# Classification logic
cat("\nApplying classification logic...\n")
classification <- "INDETERMINATE"
confidence <- "LOW"
rationale <- character()

if (monosomy_6) {
  classification <- "WNT"
  confidence <- "HIGH"
  rationale <- c("Monosomy 6 detected (100% specific for WNT)")
} else if (myc_amp) {
  classification <- "GROUP_3"
  confidence <- "HIGH"
  rationale <- c("MYC amplification detected")
} else if (chr9q_loss) {
  classification <- "SHH"
  confidence <- "MEDIUM"
  rationale <- c("Chromosome 9q loss detected")
} else if (iso17q && chr8_loss) {
  classification <- "GROUP_4"
  confidence <- "HIGH"
  rationale <- c("Isochromosome 17q detected", "Chromosome 8 loss detected")
} else if (iso17q) {
  classification <- "GROUP_3_OR_4"
  confidence <- "MEDIUM"
  rationale <- c("Isochromosome 17q detected")
}

cat(sprintf("\nClassification: %s\n", classification))
cat(sprintf("Confidence: %s\n", confidence))

# Fragmentomics integration
fragmentomics_enabled <- FALSE
sl_ratio <- NA
fragmentomics_concordance <- "N/A"

if (opt$fragmentomics != "" && file.exists(opt$fragmentomics)) {
  cat("\n--- Fragmentomics Integration ---\n")
  tryCatch({
    frag_metrics <- fromJSON(opt$fragmentomics)
    fragmentomics_enabled <- TRUE
    sl_ratio <- frag_metrics$short_long_ratio
    
    cat(sprintf("  S/L ratio: %.3f\n", sl_ratio))
    
    if (classification == "GROUP_4" && sl_ratio >= 0.8 && sl_ratio <= 1.2) {
      fragmentomics_concordance <- "CONCORDANT"
      confidence <- "VERY_HIGH"
      rationale <- c(rationale, "Fragmentomics concordant with Group 4")
      cat("  ✓ Fragmentomics CONCORDANT\n")
    }
  }, error = function(e) {
    cat(sprintf("  ✗ Error: %s\n", e$message))
  })
}

# Save results
output <- data.frame(
  sample_id = opt$id,
  tumor_fraction = tumor_fraction,
  ploidy = ploidy,
  classification = classification,
  confidence = confidence,
  monosomy_6 = monosomy_6,
  isochromosome_17q = iso17q,
  chromosome_8_loss = chr8_loss,
  fragmentomics_sl_ratio = sl_ratio,
  fragmentomics_concordance = fragmentomics_concordance,
  rationale = paste(rationale, collapse = " | ")
)

output_file <- file.path(opt$output, paste0(opt$id, "_classification.csv"))
fwrite(output, output_file)
cat(sprintf("\n✓ Saved: %s\n", output_file))