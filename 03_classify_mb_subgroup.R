#!/usr/bin/env Rscript
# Classify MB using SIMPLE PRESENCE/ABSENCE logic
# Based on Liu et al. 2021 approach

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option(c("-i", "--id"), type="character", help="Sample ID"),
  make_option(c("-s", "--seg"), type="character", help="ichorCNA seg file"),
  make_option(c("-p", "--params"), type="character", help="ichorCNA params file"),
  make_option(c("-o", "--output"), type="character", help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("=== MB Subgroup Classification (Liu 2021 Method) ===\n")

# Load data
seg_data <- fread(opt$seg)
params <- fread(opt$params)

tumor_fraction <- params$Tumor_Fraction[1]
ploidy <- params$Ploidy[1]

cat("Sample:", opt$id, "\n")
cat("Tumor Fraction:", round(tumor_fraction * 100, 2), "%\n\n")

# ==========================================
# FEATURE DETECTION FUNCTIONS
# Based on Escudero 2020 + Schwalbe 2017
# ==========================================

get_chr_cn <- function(chr, arm = NULL) {
  chr_segs <- seg_data[chr == as.character(chr)]
  if (nrow(chr_segs) == 0) return(NA)
  
  if (!is.null(arm)) {
    chr_mid <- max(chr_segs$end) / 2
    if (arm == "p") {
      chr_segs <- chr_segs[end < chr_mid]
    } else if (arm == "q") {
      chr_segs <- chr_segs[start > chr_mid]
    }
  }
  
  if (nrow(chr_segs) == 0) return(NA)
  median(chr_segs$copy.number, na.rm = TRUE)
}

# Monosomy 6 (WNT diagnostic feature)
# Escudero 2020: 100% of WNT cases
check_monosomy_6 <- function() {
  chr6_cn <- get_chr_cn(6)
  if (is.na(chr6_cn)) return(FALSE)
  chr6_cn < 1.5  # Monosomy = 1 copy
}

# 9q loss (SHH characteristic)
# Escudero 2020: 100% of SHH cases
check_chr9q_loss <- function() {
  chr9q_cn <- get_chr_cn(9, arm = "q")
  if (is.na(chr9q_cn)) return(FALSE)
  chr9q_cn < 1.5
}

# Isochromosome 17q (Group 3/4)
# Schwalbe 2017: 42% Group 3, 67% Group 4
check_iso17q <- function() {
  chr17p_cn <- get_chr_cn(17, arm = "p")
  chr17q_cn <- get_chr_cn(17, arm = "q")
  
  if (is.na(chr17p_cn) || is.na(chr17q_cn)) return(FALSE)
  
  # i17q = loss of 17p + gain of 17q
  (chr17p_cn < 1.5) && (chr17q_cn > 2.5)
}

# MYC amplification (Group 3 specific)
# Schwalbe 2017: 17% of Group 3
check_myc_amplification <- function() {
  # MYC at chr8:q24 (approximately 128-129 Mb)
  chr8q_segs <- seg_data[chr == "8" & start > 100000000]
  if (nrow(chr8q_segs) == 0) return(FALSE)
  
  max_cn_8q <- max(chr8q_segs$copy.number, na.rm = TRUE)
  max_cn_8q > 4  # Amplification = 4+ copies
}

# MYCN amplification (SHH/Group 4)
check_mycn_amplification <- function() {
  # MYCN at chr2:p24.3
  chr2p_segs <- seg_data[chr == "2" & start < 50000000]
  if (nrow(chr2p_segs) == 0) return(FALSE)
  
  max_cn_2p <- max(chr2p_segs$copy.number, na.rm = TRUE)
  max_cn_2p > 4
}

# Chr 8 loss (Group 4 characteristic)
check_chr8_loss <- function() {
  chr8_cn <- get_chr_cn(8)
  if (is.na(chr8_cn)) return(FALSE)
  chr8_cn < 1.5
}

# 10q loss (SHH)
check_chr10q_loss <- function() {
  chr10q_cn <- get_chr_cn(10, arm = "q")
  if (is.na(chr10q_cn)) return(FALSE)
  chr10q_cn < 1.5
}

# ==========================================
# CLASSIFICATION LOGIC
# Based on Liu et al. 2021 decision tree
# ==========================================

# Detect all features
features <- list(
  monosomy_6 = check_monosomy_6(),
  chr9q_loss = check_chr9q_loss(),
  chr10q_loss = check_chr10q_loss(),
  iso17q = check_iso17q(),
  myc_amp = check_myc_amplification(),
  mycn_amp = check_mycn_amplification(),
  chr8_loss = check_chr8_loss()
)

cat("Detected Features:\n")
for (feature in names(features)) {
  if (features[[feature]]) {
    cat(sprintf("  ✓ %s\n", gsub("_", " ", feature)))
  }
}
cat("\n")

# Classification decision tree (Liu 2021 logic)

if (features$monosomy_6) {
  # Monosomy 6 is DIAGNOSTIC for WNT
  # Escudero 2020: 100% specific
  classification <- list(
    subgroup = "WNT",
    confidence = "High",
    rationale = "Monosomy 6 detected (diagnostic for WNT-activated MB)",
    key_features = "Monosomy 6"
  )
  
} else if (features$myc_amp) {
  # MYC amplification is SPECIFIC to Group 3
  # Schwalbe 2017: 17% of Group 3, 0% of others
  classification <- list(
    subgroup = "Group 3",
    confidence = "High",
    rationale = "MYC amplification detected (specific to Group 3)",
    key_features = "MYC amplification"
  )
  
} else if (features$chr9q_loss || (features$chr9q_loss && features$chr10q_loss)) {
  # 9q loss characteristic of SHH
  # Escudero 2020: 100% of SHH in their cohort
  classification <- list(
    subgroup = "SHH",
    confidence = "High",
    rationale = "9q loss detected (characteristic of SHH-activated MB)",
    key_features = "9q deletion"
  )
  
} else if (features$iso17q && features$chr8_loss) {
  # i17q + chr8 loss suggests Group 4
  # Schwalbe 2017: Chr 8 loss in 20% of Group 4
  classification <- list(
    subgroup = "Group 4",
    confidence = "High",
    rationale = "Isochromosome 17q + chromosome 8 loss (characteristic of Group 4)",
    key_features = "i17q, chr8 loss"
  )
  
} else if (features$iso17q && !features$chr8_loss) {
  # i17q alone - could be Group 3 or Group 4
  # Schwalbe 2017: 42% Group 3, 67% Group 4
  classification <- list(
    subgroup = "Group 3 or Group 4",
    confidence = "Medium",
    rationale = "Isochromosome 17q present but cannot distinguish between Group 3 and Group 4",
    key_features = "i17q"
  )
  
} else if (features$mycn_amp) {
  # MYCN amplification - SHH or Group 4
  # Schwalbe 2017: 7% SHH, 4% Group 4
  classification <- list(
    subgroup = "SHH or Group 4",
    confidence = "Medium",
    rationale = "MYCN amplification detected",
    key_features = "MYCN amplification"
  )
  
} else {
  # No diagnostic features
  classification <- list(
    subgroup = "Unclassified",
    confidence = "Low",
    rationale = "No diagnostic CNV pattern detected. May represent low tumor fraction or non-MB tumor.",
    key_features = "None"
  )
}

# Print results
cat("=== CLASSIFICATION RESULT ===\n")
cat("MB Subgroup:", classification$subgroup, "\n")
cat("Confidence:", classification$confidence, "\n")
cat("Rationale:", classification$rationale, "\n")
cat("Key Features:", classification$key_features, "\n")

# Save results
output_data <- data.frame(
  Sample_ID = opt$id,
  Tumor_Fraction = tumor_fraction,
  Ploidy = ploidy,
  MB_Subgroup = classification$subgroup,
  Confidence = classification$confidence,
  Rationale = classification$rationale,
  Monosomy_6 = features$monosomy_6,
  Chr9q_Loss = features$chr9q_loss,
  Iso17q = features$iso17q,
  MYC_amp = features$myc_amp,
  MYCN_amp = features$mycn_amp,
  Chr8_loss = features$chr8_loss,
  Timestamp = Sys.time()
)

fwrite(output_data, file.path(opt$output, paste0(opt$id, "_classification.csv")))

cat("\n✓ Classification complete\n")# Based on CNV signatures from published papers
# ============================================

classify_medulloblastoma <- function(seg_data, config, ploidy) {
  
  # Initialize scoring system
  scores <- list(
    WNT = 0,
    SHH = 0,
    Group3 = 0,
    Group4 = 0
  )
  
  features_detected <- list()
  
  # Helper function to calculate median CN for chromosome/region
  get_chr_cn <- function(chr, region=NULL, arm=NULL) {
    if (!is.null(arm)) {
      # Get centromere position for arm-level analysis
      # Simplified: p-arm = start of chr, q-arm = middle to end
      chr_segs <- seg_data[chr == as.character(chr)]
      if (nrow(chr_segs) == 0) return(NA)
      
      chr_mid <- max(chr_segs$end) / 2
      if (arm == "p") {
        chr_segs <- chr_segs[end < chr_mid]
      } else if (arm == "q") {
        chr_segs <- chr_segs[start > chr_mid]
      }
    } else {
      chr_segs <- seg_data[chr == as.character(chr)]
    }
    
    if (nrow(chr_segs) == 0) return(NA)
    median(chr_segs$copy.number, na.rm=TRUE)
  }
  
  # =====================================
  # WNT SIGNATURE DETECTION
  # =====================================
  
  # Monosomy 6 (DIAGNOSTIC for WNT)
  chr6_cn <- get_chr_cn(6)
  if (!is.na(chr6_cn) && chr6_cn < 1.5) {
    scores$WNT <- scores$WNT + 10  # Strong evidence
    features_detected$monosomy_6 <- TRUE
  }
  
  # Chr 7 gain (common in WNT)
  chr7_cn <- get_chr_cn(7)
  if (!is.na(chr7_cn) && chr7_cn > 2.5) {
    scores$WNT <- scores$WNT + 2
    features_detected$chr7_gain <- TRUE
  }
  
  # =====================================
  # SHH SIGNATURE DETECTION
  # =====================================
  
  # 9q loss (PTCH1 region)
  chr9q_cn <- get_chr_cn(9, arm="q")
  if (!is.na(chr9q_cn) && chr9q_cn < 1.5) {
    scores$SHH <- scores$SHH + 5
    features_detected$chr9q_loss <- TRUE
  }
  
  # 10q loss
  chr10q_cn <- get_chr_cn(10, arm="q")
  if (!is.na(chr10q_cn) && chr10q_cn < 1.5) {
    scores$SHH <- scores$SHH + 3
    features_detected$chr10q_loss <- TRUE
  }
  
  # 17p loss (TP53)
  chr17p_cn <- get_chr_cn(17, arm="p")
  if (!is.na(chr17p_cn) && chr17p_cn < 1.5) {
    scores$SHH <- scores$SHH + 3
    scores$Group3 <- scores$Group3 + 1  # Also seen in G3
    features_detected$chr17p_loss <- TRUE
  }
  
  # MYCN amplification (chr 2p24)
  chr2p_segs <- seg_data[chr == "2" & start < 50000000]
  if (nrow(chr2p_segs) > 0) {
    max_cn_2p <- max(chr2p_segs$copy.number, na.rm=TRUE)
    if (max_cn_2p > 4) {
      scores$SHH <- scores$SHH + 4
      scores$Group4 <- scores$Group4 + 2
      features_detected$mycn_amplification <- TRUE
    }
  }
  
  # =====================================
  # ISOCHROMOSOME 17q (Group 3/4)
  # =====================================
  
  chr17q_cn <- get_chr_cn(17, arm="q")
  
  # i17q = 17p loss + 17q gain
  if (!is.na(chr17p_cn) && !is.na(chr17q_cn)) {
    if (chr17p_cn < 1.5 && chr17q_cn > 2.5) {
      scores$Group3 <- scores$Group3 + 6
      scores$Group4 <- scores$Group4 + 6
      features_detected$iso17q <- TRUE
    }
  }
  
  # =====================================
  # GROUP 3 SPECIFIC
  # =====================================
  
  # MYC amplification (chr 8q24) - HIGH RISK
  chr8q_segs <- seg_data[chr == "8" & start > 100000000]
  if (nrow(chr8q_segs) > 0) {
    max_cn_8q <- max(chr8q_segs$copy.number, na.rm=TRUE)
    if (max_cn_8q > 4) {
      scores$Group3 <- scores$Group3 + 8  # Strong evidence
      features_detected$myc_amplification <- TRUE
    }
  }
  
  # OTX2 amplification (chr 14q22)
  chr14q_segs <- seg_data[chr == "14" & start > 50000000]
  if (nrow(chr14q_segs) > 0) {
    max_cn_14q <- max(chr14q_segs$copy.number, na.rm=TRUE)
    if (max_cn_14q > 3) {
      scores$Group3 <- scores$Group3 + 3
      features_detected$otx2_amplification <- TRUE
    }
  }
  
  # =====================================
  # GROUP 4 SPECIFIC
  # =====================================
  
  # Chr 8 loss (distinguishes G4 from G3)
  chr8_cn <- get_chr_cn(8)
  if (!is.na(chr8_cn) && chr8_cn < 1.5) {
    scores$Group4 <- scores$Group4 + 4
    scores$Group3 <- scores$Group3 - 2  # Less likely G3
    features_detected$chr8_loss <- TRUE
  }
  
  # Chr X loss (in females)
  chrX_cn <- get_chr_cn("X")
  if (!is.na(chrX_cn) && chrX_cn < 1.5) {
    scores$Group4 <- scores$Group4 + 3
    features_detected$chrX_loss <- TRUE
  }
  
  # CDK6 amplification (chr 7q21)
  chr7q_segs <- seg_data[chr == "7" & start > 80000000]
  if (nrow(chr7q_segs) > 0) {
    max_cn_7q <- max(chr7q_segs$copy.number, na.rm=TRUE)
    if (max_cn_7q > 3) {
      scores$Group4 <- scores$Group4 + 2
      features_detected$cdk6_amplification <- TRUE
    }
  }
  
  # =====================================
  # DETERMINE FINAL CLASSIFICATION
  # =====================================
  
  # Find subgroup with highest score
  max_score <- max(unlist(scores))
  
  if (max_score == 0) {
    classification <- "Unclassified"
    confidence <- "Low"
    second_best <- "None"
  } else {
    # Primary classification
    classification <- names(which.max(unlist(scores)))
    
    # Calculate confidence based on score separation
    sorted_scores <- sort(unlist(scores), decreasing=TRUE)
    score_diff <- sorted_scores[1] - sorted_scores[2]
    
    if (sorted_scores[1] >= 10 && score_diff >= 3) {
      confidence <- "High"
    } else if (sorted_scores[1] >= 5 && score_diff >= 2) {
      confidence <- "Medium"
    } else {
      confidence <- "Low"
    }
    
    # Second best
    second_best <- names(sorted_scores)[2]
  }
  
  # Return comprehensive result
  list(
    classification = classification,
    confidence = confidence,
    scores = scores,
    second_best = second_best,
    features = features_detected,
    cnv_data = list(
      chr6_cn = chr6_cn,
      chr7_cn = chr7_cn,
      chr8_cn = chr8_cn,
      chr17p_cn = chr17p_cn,
      chr17q_cn = chr17q_cn
    )
  )
}

# ============================================
# RUN CLASSIFICATION
# ============================================

result <- classify_medulloblastoma(seg_data, config, ploidy)

# ============================================
# PRINT RESULTS
# ============================================

cat("\n=== CLASSIFICATION RESULTS ===\n\n")
cat("PRIMARY CLASSIFICATION:", result$classification, "\n")
cat("Confidence:", result$confidence, "\n")

if (result$second_best != "None") {
  cat("Second possibility:", result$second_best, "\n")
}

cat("\nSubgroup Scores:\n")
for (subgroup in names(result$scores)) {
  cat(sprintf("  %s: %d\n", subgroup, result$scores[[subgroup]]))
}

cat("\nKey Features Detected:\n")
if (length(result$features) > 0) {
  for (feature in names(result$features)) {
    cat(sprintf("  ✓ %s\n", gsub("_", " ", feature)))
  }
} else {
  cat("  (No diagnostic features detected)\n")
}

cat("\nChromosome Copy Numbers:\n")
cat(sprintf("  Chr 6: %.2f\n", result$cnv_data$chr6_cn))
cat(sprintf("  Chr 7: %.2f\n", result$cnv_data$chr7_cn))
cat(sprintf("  Chr 8: %.2f\n", result$cnv_data$chr8_cn))
cat(sprintf("  Chr 17p: %.2f\n", result$cnv_data$chr17p_cn))
cat(sprintf("  Chr 17q: %.2f\n", result$cnv_data$chr17q_cn))

# ============================================
# SAVE RESULTS
# ============================================

output_data <- data.frame(
  Sample_ID = opt$id,
  Tumor_Fraction = tumor_fraction,
  Ploidy = ploidy,
  MB_Subgroup = result$classification,
  Confidence = result$confidence,
  Second_Best = result$second_best,
  WNT_Score = result$scores$WNT,
  SHH_Score = result$scores$SHH,
  Group3_Score = result$scores$Group3,
  Group4_Score = result$scores$Group4,
  Monosomy_6 = "monosomy_6" %in% names(result$features),
  Iso17q = "iso17q" %in% names(result$features),
  MYC_amp = "myc_amplification" %in% names(result$features),
  MYCN_amp = "mycn_amplification" %in% names(result$features),
  Chr8_loss = "chr8_loss" %in% names(result$features),
  Timestamp = Sys.time()
)

fwrite(output_data, 
       file.path(opt$output, paste0(opt$id, "_classification.csv")))

cat("\n=== Classification Complete ===\n")
cat("Results saved to:", file.path(opt$output, paste0(opt$id, "_classification.csv")), "\n")
