#!/usr/bin/env Rscript
# MB Risk Stratification
# Based on Schwalbe et al. 2017 + Taylor et al. 2012
# Updated: 2024

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# Parse command line options
option_list <- list(
  make_option("--id", type="character", help="Sample ID"),
  make_option("--classification", type="character", 
              help="Classification CSV file"),
  make_option("--seg", type="character", help="ichorCNA seg file"),
  make_option("--output", type="character", help="Output CSV file"),
  make_option("--metastasis", type="character", default="unknown",
              help="Metastasis status: yes/no/unknown"),
  make_option("--residual", type="character", default="unknown",
              help="Residual disease: yes/no/unknown"),
  make_option("--histology", type="character", default="unknown",
              help="Histology: LCA/classic/DN/unknown")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate inputs
if (is.null(opt$id) || is.null(opt$classification) || 
    is.null(opt$seg) || is.null(opt$output)) {
  stop("Missing required arguments")
}

cat("=" , rep("=", 68), "=\n", sep="")
cat("MB RISK STRATIFICATION\n")
cat("=" , rep("=", 68), "=\n", sep="")
cat("Sample:", opt$id, "\n")
cat("Classification:", opt$classification, "\n")
cat("Metastasis:", opt$metastasis, "\n")
cat("Residual disease:", opt$residual, "\n")
cat("Histology:", opt$histology, "\n\n")

# Load classification results
classification <- fread(opt$classification)

subgroup <- classification$classification[1]
confidence <- classification$confidence[1]
tumor_fraction <- classification$tumor_fraction[1]

cat("Molecular subgroup:", subgroup, "\n")
cat("Classification confidence:", confidence, "\n")
cat("Tumor fraction:", round(tumor_fraction * 100, 2), "%\n\n")

# Load seg data for additional markers
seg <- fread(opt$seg)

# ============================================
# RISK STRATIFICATION LOGIC
# Based on Schwalbe et al. 2017
# ============================================

risk_group <- "UNKNOWN"
risk_factors <- c()
notes <- c()

# WNT - typically low risk
if (subgroup == "WNT") {
  risk_group <- "LOW"
  risk_factors <- c("WNT subgroup (excellent prognosis, >95% survival)")
  notes <- c("Consider treatment de-escalation")
}

# SHH - variable risk
if (subgroup == "SHH") {
  # Check for TP53 status (would need mutation data)
  # Check for age (would need clinical data)
  # Default to intermediate
  risk_group <- "INTERMEDIATE"
  risk_factors <- c("SHH subgroup (variable prognosis)")
  
  # High-risk factors
  if (opt$metastasis == "yes") {
    risk_group <- "HIGH"
    risk_factors <- c(risk_factors, "Metastatic disease at diagnosis")
  }
  
  if (opt$histology == "LCA") {
    risk_group <- "HIGH"
    risk_factors <- c(risk_factors, "Large cell/anaplastic histology")
  }
  
  notes <- c("TP53 status needed for final risk assignment")
}

# Group 3 - typically high risk
if (subgroup == "GROUP_3") {
  # Check for MYC amplification
  if (classification$myc_amp[1] == TRUE) {
    risk_group <- "HIGH"
    risk_factors <- c("Group 3 with MYC amplification (poor prognosis)")
  } else {
    risk_group <- "INTERMEDIATE_TO_HIGH"
    risk_factors <- c("Group 3 subgroup")
  }
  
  if (opt$metastasis == "yes") {
    risk_group <- "HIGH"
    risk_factors <- c(risk_factors, "Metastatic disease at diagnosis")
  }
}

# Group 4 - typically intermediate risk
if (subgroup == "GROUP_4") {
  risk_group <- "INTERMEDIATE"
  risk_factors <- c("Group 4 subgroup (intermediate prognosis)")
  
  # High-risk factors
  if (opt$metastasis == "yes") {
    risk_group <- "HIGH"
    risk_factors <- c(risk_factors, "Metastatic disease at diagnosis")
  }
  
  if (opt$residual == "yes") {
    risk_group <- "HIGH"
    risk_factors <- c(risk_factors, "Residual disease post-surgery")
  }
  
  if (opt$histology == "LCA") {
    risk_factors <- c(risk_factors, "Large cell/anaplastic histology (increases risk)")
  }
}

# Group 3/4 uncertain
if (subgroup == "GROUP_3_OR_4") {
  risk_group <- "INTERMEDIATE"
  risk_factors <- c("Group 3/4 ambiguous - chromosome 8 status needed")
  notes <- c("Complete classification needed for accurate risk stratification")
}

# Indeterminate
if (subgroup == "INDETERMINATE") {
  risk_group <- "UNKNOWN"
  risk_factors <- c("Molecular classification incomplete")
  notes <- c("Cannot assign risk without molecular subgroup")
}

# Additional factors
if (tumor_fraction < 0.10) {
  notes <- c(notes, "Low tumor fraction may affect classification accuracy")
}

# Check for fragmentomics concordance
if ("fragmentomics_concordance" %in% names(classification)) {
  if (!is.na(classification$fragmentomics_concordance[1]) && 
      classification$fragmentomics_concordance[1] == "CONCORDANT") {
    notes <- c(notes, "Fragmentomics concordant with classification (Markowitz 2025)")
  }
}

# ============================================
# OUTPUT
# ============================================
cat("\n")
cat("=" , rep("=", 68), "=\n", sep="")
cat("RISK STRATIFICATION RESULT\n")
cat("=" , rep("=", 68), "=\n", sep="")
cat("Risk group:", risk_group, "\n")
cat("Risk factors:\n")
for (i in seq_along(risk_factors)) {
  cat("  ", i, ". ", risk_factors[i], "\n", sep="")
}

if (length(notes) > 0) {
  cat("\nNotes:\n")
  for (i in seq_along(notes)) {
    cat("  • ", notes[i], "\n", sep="")
  }
}
cat("\n")

# Create output dataframe
results <- data.frame(
  sample_id = opt$id,
  molecular_subgroup = subgroup,
  risk_group = risk_group,
  metastasis = opt$metastasis,
  residual_disease = opt$residual,
  histology = opt$histology,
  tumor_fraction = tumor_fraction,
  risk_factors = paste(risk_factors, collapse = "; "),
  notes = paste(notes, collapse = "; "),
  stringsAsFactors = FALSE
)

# Save
write.csv(results, opt$output, row.keys = FALSE, quote = TRUE)

cat("Risk stratification saved to:", opt$output, "\n")
cat("\nIMPORTANT LIMITATIONS:\n")
cat("• Risk stratification requires integration of:\n")
cat("  - Molecular subgroup (✓ provided by cfDNA)\n")
cat("  - Histology (requires tissue)\n")
cat("  - Metastatic status (requires MRI)\n")
cat("  - Age and resection extent (requires clinical data)\n")
cat("  - Specific mutations like TP53 (requires sequencing)\n")
cat("• This analysis provides molecular component only\n")
cat("• Final risk category requires multidisciplinary integration\n\n")