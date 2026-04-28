#!/usr/bin/env Rscript
# MB Risk Stratification
# Based on Schwalbe et al. 2017 (Lancet Oncology)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option(c("-i", "--id"), type="character"),
  make_option(c("-c", "--classification"), type="character"),
  make_option(c("-s", "--seg"), type="character"),
  make_option(c("-o", "--output"), type="character"),
  make_option(c("--metastasis"), type="character", default="unknown"),
  make_option(c("--residual"), type="character", default="unknown"),
  make_option(c("--histology"), type="character", default="unknown")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Load classification
class_data <- fread(opt$classification)
seg_data <- fread(opt$seg)

subgroup <- class_data$MB_Subgroup[1]
tumor_fraction <- class_data$Tumor_Fraction[1]

# Schwalbe 2017 risk stratification
# 4 groups: Favourable, Standard, High Risk, Very High Risk

risk_stratify <- function(subgroup, features, clinical) {
  
  # WNT = Favourable (unless very rare exceptions)
  if (subgroup == "WNT") {
    return(list(
      risk_group = "Favourable",
      pfs_5year = ">90%",
      rationale = "WNT-activated MB has excellent prognosis"
    ))
  }
  
  # SHH risk depends on TP53 and other factors
  if (subgroup == "SHH") {
    if (features$tp53_mut || features$mycn_amp || features$lca_histology) {
      return(list(
        risk_group = "High Risk",
        pfs_5year = "~50-70%",
        rationale = "SHH with TP53 mutation, MYCN amplification, or LCA histology"
      ))
    } else {
      return(list(
        risk_group = "Standard",
        pfs_5year = "~75-85%",
        rationale = "SHH without high-risk features"
      ))
    }
  }
  
  # Group 3 - generally high risk
  if (subgroup == "Group3") {
    if (features$myc_amp || clinical$metastasis) {
      return(list(
        risk_group = "Very High Risk",
        pfs_5year = "<40%",
        rationale = "Group 3 with MYC amplification and/or metastasis"
      ))
    } else {
      return(list(
        risk_group = "High Risk",
        pfs_5year = "~50-60%",
        rationale = "Group 3 MB"
      ))
    }
  }
  
  # Group 4 - intermediate risk
  if (subgroup == "Group4") {
    if (features$mycn_amp || clinical$metastasis || clinical$residual) {
      return(list(
        risk_group = "High Risk",
        pfs_5year = "~60-70%",
        rationale = "Group 4 with high-risk features"
      ))
    } else if (features$chr13_loss) {
      return(list(
        risk_group = "Favourable",
        pfs_5year = "~85-90%",
        rationale = "Group 4 with chr 13 loss (protective)"
      ))
    } else {
      return(list(
        risk_group = "Standard",
        pfs_5year = "~75-80%",
        rationale = "Group 4 without risk modifiers"
      ))
    }
  }
  
  # Unclassified
  return(list(
    risk_group = "Unknown",
    pfs_5year = "N/A",
    rationale = "Unable to classify MB subgroup"
  ))
}

# Extract features from classification
features <- list(
  myc_amp = class_data$MYC_amp[1],
  mycn_amp = class_data$MYCN_amp[1],
  iso17q = class_data$Iso17q[1],
  tp53_mut = FALSE,  # Would need sequencing
  lca_histology = opt$histology == "LCA",
  chr13_loss = FALSE  # Check from seg data
)

# Check for chr 13 loss
chr13_segs <- seg_data[chr == "13"]
if (nrow(chr13_segs) > 0) {
  chr13_cn <- median(chr13_segs$copy.number)
  features$chr13_loss <- chr13_cn < 1.5
}

# Clinical features
clinical <- list(
  metastasis = opt$metastasis == "yes",
  residual = opt$residual == "yes"
)

# Perform risk stratification
risk_result <- risk_stratify(subgroup, features, clinical)

# Print results
cat("\n=== RISK STRATIFICATION ===\n\n")
cat("MB Subgroup:", subgroup, "\n")
cat("Risk Group:", risk_result$risk_group, "\n")
cat("Expected 5-year PFS:", risk_result$pfs_5year, "\n")
cat("Rationale:", risk_result$rationale, "\n\n")

# Save
output_data <- data.frame(
  Sample_ID = opt$id,
  MB_Subgroup = subgroup,
  Risk_Group = risk_result$risk_group,
  PFS_5year = risk_result$pfs_5year,
  Rationale = risk_result$rationale,
  MYC_amplification = features$myc_amp,
  MYCN_amplification = features$mycn_amp,
  Chr13_loss = features$chr13_loss,
  Metastasis = clinical$metastasis,
  Residual_disease = clinical$residual,
  Timestamp = Sys.time()
)

fwrite(output_data, opt$output)
cat("Results saved to:", opt$output, "\n")
