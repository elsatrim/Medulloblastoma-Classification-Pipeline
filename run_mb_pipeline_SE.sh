#!/bin/bash
# Medulloblastoma Classification Pipeline - SINGLE-END VERSION
# For use with single-end (SE) sequencing data (one FASTQ file)
# Modified from paired-end pipeline for compatibility with older datasets

set -euo pipefail

# ==========================================
# CONFIGURATION
# ==========================================

# Reference genome (update this path if needed)
REFERENCE="${HOME}/mb_reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# ichorCNA installation
ICHORCNA_DIR="${HOME}/ichorCNA"

# Default parameters
THREADS=4
METASTASIS="unknown"
RESIDUAL="unknown"
HISTOLOGY="unknown"

# ==========================================
# USAGE
# ==========================================

usage() {
    cat << EOF
Medulloblastoma Classification Pipeline - SINGLE-END VERSION

Usage: $0 SAMPLE_ID FASTQ OUTPUT_DIR [OPTIONS]

Required Arguments:
  SAMPLE_ID       Unique sample identifier
  FASTQ           Path to single-end FASTQ file (.fastq.gz or .fastq)
  OUTPUT_DIR      Directory for output files

Optional Arguments:
  --threads N     Number of threads (default: 4)
  --metastasis    Metastasis status: yes/no/unknown (default: unknown)
  --residual      Residual disease: yes/no/unknown (default: unknown)
  --histology     Histology: LCA/classic/DN/unknown (default: unknown)
  --help          Show this help message

Example:
  $0 MB001 data/MB001.fastq.gz results/MB001 --threads 8 --metastasis no

Note: This is the SINGLE-END version for samples with only one FASTQ file.
      For paired-end data (R1 + R2), use run_mb_pipeline.sh instead.

EOF
    exit 1
}

# ==========================================
# PARSE ARGUMENTS
# ==========================================

if [[ $# -lt 3 ]]; then
    echo "Error: Missing required arguments"
    usage
fi

SAMPLE_ID=$1
FASTQ=$2
OUTPUT_DIR=$3
shift 3

while [[ $# -gt 0 ]]; do
    case $1 in
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --metastasis)
            METASTASIS="$2"
            shift 2
            ;;
        --residual)
            RESIDUAL="$2"
            shift 2
            ;;
        --histology)
            HISTOLOGY="$2"
            shift 2
            ;;
        --help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# ==========================================
# VALIDATION
# ==========================================

echo "=========================================="
echo "MB Classification Pipeline (Single-End)"
echo "=========================================="
echo ""
echo "Sample ID: ${SAMPLE_ID}"
echo "FASTQ: ${FASTQ}"
echo "Output: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"
echo "Metastasis: ${METASTASIS}"
echo "Residual: ${RESIDUAL}"
echo "Histology: ${HISTOLOGY}"
echo ""

# Check if FASTQ exists
if [[ ! -f "$FASTQ" ]]; then
    echo "Error: FASTQ file not found: $FASTQ"
    exit 1
fi

# Check if reference exists
if [[ ! -f "$REFERENCE" ]]; then
    echo "Error: Reference genome not found: $REFERENCE"
    echo "Please update REFERENCE variable in this script"
    exit 1
fi

# Check if reference is indexed
if [[ ! -f "${REFERENCE}.bwt" ]]; then
    echo "Error: Reference genome not indexed for BWA"
    echo "Run: bwa index $REFERENCE"
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# ==========================================
# STEP 1: PREPROCESSING (SINGLE-END)
# ==========================================

echo "=========================================="
echo "STEP 1: Preprocessing (Single-End Mode)"
echo "=========================================="
echo ""

PREPROCESS_DIR="${OUTPUT_DIR}/01_preprocessing"
mkdir -p "${PREPROCESS_DIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running FastQC..."

fastqc \
    -o "${PREPROCESS_DIR}" \
    -t ${THREADS} \
    "${FASTQ}" \
    > "${PREPROCESS_DIR}/fastqc.log" 2>&1

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Aligning reads (BWA-MEM, single-end mode)..."

# Single-end alignment (no R2 file)
bwa mem \
    -t ${THREADS} \
    -M \
    -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA" \
    "${REFERENCE}" \
    "${FASTQ}" \
    2> "${PREPROCESS_DIR}/bwa.log" \
    | samtools view -Sb - \
    | samtools sort -@ ${THREADS} -o "${PREPROCESS_DIR}/${SAMPLE_ID}.sorted.bam"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Indexing BAM..."

samtools index "${PREPROCESS_DIR}/${SAMPLE_ID}.sorted.bam"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Calculating alignment statistics..."

samtools flagstat "${PREPROCESS_DIR}/${SAMPLE_ID}.sorted.bam" \
    > "${PREPROCESS_DIR}/${SAMPLE_ID}.flagstat.txt"

samtools stats "${PREPROCESS_DIR}/${SAMPLE_ID}.sorted.bam" \
    > "${PREPROCESS_DIR}/${SAMPLE_ID}.stats.txt"

# Calculate coverage
samtools depth "${PREPROCESS_DIR}/${SAMPLE_ID}.sorted.bam" \
    | awk '{sum+=$3} END {print "Average coverage:", sum/NR}' \
    > "${PREPROCESS_DIR}/${SAMPLE_ID}.coverage.txt"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Preprocessing complete!"
echo ""

# ==========================================
# STEP 2: CNV DETECTION (ichorCNA)
# ==========================================

echo "=========================================="
echo "STEP 2: Copy Number Variation Detection"
echo "=========================================="
echo ""

ICHORCNA_DIR_OUT="${OUTPUT_DIR}/02_ichorCNA"
mkdir -p "${ICHORCNA_DIR_OUT}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running ichorCNA..."

# Check if ichorCNA R script exists
if [[ ! -f "scripts/02_run_ichorCNA.R" ]]; then
    echo "Error: ichorCNA script not found: scripts/02_run_ichorCNA.R"
    exit 1
fi

Rscript scripts/02_run_ichorCNA.R \
    --id "${SAMPLE_ID}" \
    --bam "${PREPROCESS_DIR}/${SAMPLE_ID}.sorted.bam" \
    --outdir "${ICHORCNA_DIR_OUT}" \
    --threads ${THREADS} \
    > "${ICHORCNA_DIR_OUT}/ichorCNA.log" 2>&1

# Check if ichorCNA completed successfully
if [[ ! -f "${ICHORCNA_DIR_OUT}/${SAMPLE_ID}.params.txt" ]]; then
    echo "Error: ichorCNA failed. Check log: ${ICHORCNA_DIR_OUT}/ichorCNA.log"
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ichorCNA complete!"
echo ""

# Extract tumor fraction
TF=$(grep "Tumor Fraction" "${ICHORCNA_DIR_OUT}/${SAMPLE_ID}.params.txt" | awk '{print $3}')
echo "Tumor Fraction: ${TF}"

# Check if tumor fraction is sufficient
if (( $(echo "$TF < 0.05" | bc -l) )); then
    echo "WARNING: Low tumor fraction (${TF}). Results may be unreliable."
    echo "Recommended minimum: 8% (0.08)"
fi

echo ""

# ==========================================
# STEP 3: CLASSIFICATION
# ==========================================

echo "=========================================="
echo "STEP 3: MB Subgroup Classification"
echo "=========================================="
echo ""

CLASS_DIR="${OUTPUT_DIR}/03_classification"
mkdir -p "${CLASS_DIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Classifying MB subgroup..."

# Check if classification script exists
if [[ ! -f "scripts/03_classify_mb_subgroup_v2.R" ]]; then
    echo "Error: Classification script not found: scripts/03_classify_mb_subgroup_v2.R"
    exit 1
fi

Rscript scripts/03_classify_mb_subgroup_v2.R \
    --sample "${SAMPLE_ID}" \
    --ichor_dir "${ICHORCNA_DIR_OUT}" \
    --outdir "${CLASS_DIR}" \
    > "${CLASS_DIR}/classification.log" 2>&1

# Check if classification completed
if [[ ! -f "${CLASS_DIR}/${SAMPLE_ID}_classification.csv" ]]; then
    echo "Error: Classification failed. Check log: ${CLASS_DIR}/classification.log"
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Classification complete!"

# Display results
echo ""
echo "Classification Results:"
cat "${CLASS_DIR}/${SAMPLE_ID}_classification.csv"
echo ""

# ==========================================
# STEP 4: RISK STRATIFICATION
# ==========================================

echo "=========================================="
echo "STEP 4: Risk Stratification"
echo "=========================================="
echo ""

RISK_DIR="${OUTPUT_DIR}/04_risk_stratification"
mkdir -p "${RISK_DIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Calculating risk group..."

# Check if risk script exists
if [[ ! -f "scripts/05_risk_stratification.R" ]]; then
    echo "Error: Risk stratification script not found"
    exit 1
fi

Rscript scripts/05_risk_stratification.R \
    --sample "${SAMPLE_ID}" \
    --classification "${CLASS_DIR}/${SAMPLE_ID}_classification.csv" \
    --metastasis "${METASTASIS}" \
    --residual "${RESIDUAL}" \
    --outdir "${RISK_DIR}" \
    > "${RISK_DIR}/risk.log" 2>&1

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Risk stratification complete!"

if [[ -f "${RISK_DIR}/${SAMPLE_ID}_risk.csv" ]]; then
    echo ""
    echo "Risk Stratification:"
    cat "${RISK_DIR}/${SAMPLE_ID}_risk.csv"
    echo ""
fi

# ==========================================
# STEP 5: REPORT GENERATION
# ==========================================

echo "=========================================="
echo "STEP 5: Clinical Report Generation"
echo "=========================================="
echo ""

REPORT_DIR="${OUTPUT_DIR}/05_report"
mkdir -p "${REPORT_DIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generating clinical report..."

# Check if report script exists
if [[ ! -f "scripts/06_generate_report.Rmd" ]]; then
    echo "Warning: Report template not found. Skipping report generation."
else
    Rscript -e "
    rmarkdown::render(
        'scripts/06_generate_report.Rmd',
        params = list(
            sample_id = '${SAMPLE_ID}',
            results_dir = '${OUTPUT_DIR}',
            metastasis = '${METASTASIS}',
            residual = '${RESIDUAL}',
            histology = '${HISTOLOGY}'
        ),
        output_file = '${SAMPLE_ID}_MB_Report.pdf',
        output_dir = '${REPORT_DIR}'
    )
    " > "${REPORT_DIR}/report.log" 2>&1

    if [[ -f "${REPORT_DIR}/${SAMPLE_ID}_MB_Report.pdf" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Report generated successfully!"
        echo "Report: ${REPORT_DIR}/${SAMPLE_ID}_MB_Report.pdf"
    else
        echo "Warning: Report generation failed. Check log: ${REPORT_DIR}/report.log"
    fi
fi

echo ""

# ==========================================
# COMPLETION
# ==========================================

echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
echo ""
echo "Sample: ${SAMPLE_ID}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Key Results:"
echo "  - Alignment: ${PREPROCESS_DIR}/${SAMPLE_ID}.sorted.bam"
echo "  - CNV calls: ${ICHORCNA_DIR_OUT}/${SAMPLE_ID}.seg.txt"
echo "  - Classification: ${CLASS_DIR}/${SAMPLE_ID}_classification.csv"
echo "  - Risk group: ${RISK_DIR}/${SAMPLE_ID}_risk.csv"
if [[ -f "${REPORT_DIR}/${SAMPLE_ID}_MB_Report.pdf" ]]; then
    echo "  - Clinical report: ${REPORT_DIR}/${SAMPLE_ID}_MB_Report.pdf"
fi
echo ""
echo "Total runtime: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

exit 0
