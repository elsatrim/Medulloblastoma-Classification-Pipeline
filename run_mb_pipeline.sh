#!/bin/bash
# Complete MB Classification Pipeline
# Based on published methods (Liu 2021, Markowitz 2025, Escudero 2020)

set -e

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <sample_id> <fastq_r1> <fastq_r2> <output_dir> [options]"
    echo ""
    echo "Options:"
    echo "  --threads N          Number of threads (default: 8)"
    echo "  --metastasis yes/no  Metastasis at diagnosis"
    echo "  --residual yes/no    Residual disease post-surgery"
    echo "  --histology LCA/classic/DN  Histology type"
    echo ""
    echo "Example:"
    echo "  $0 MB001 sample_R1.fq.gz sample_R2.fq.gz results/MB001 \\"
    echo "     --threads 8 --metastasis no --residual no --histology classic"
    exit 1
fi

SAMPLE_ID=$1
FASTQ_R1=$2
FASTQ_R2=$3
OUTPUT_DIR=$4
shift 4

# Parse options
THREADS=8
METASTASIS="unknown"
RESIDUAL="unknown"
HISTOLOGY="unknown"

while [[ $# -gt 0 ]]; do
    case $1 in
        --threads) THREADS="$2"; shift 2 ;;
        --metastasis) METASTASIS="$2"; shift 2 ;;
        --residual) RESIDUAL="$2"; shift 2 ;;
        --histology) HISTOLOGY="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# MODIFY THESE PATHS FOR YOUR SYSTEM
REFERENCE="${HOME}/mb_reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
ICHORCNA_PATH="${HOME}/ichorCNA"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config"

echo "=========================================="
echo "MB CLASSIFICATION PIPELINE"
echo "=========================================="
echo "Sample: ${SAMPLE_ID}"
echo "Start: $(date)"
echo ""

# Create directory structure
mkdir -p ${OUTPUT_DIR}
PREPROCESS_DIR="${OUTPUT_DIR}/01_preprocessing"
ICHORCNA_DIR="${OUTPUT_DIR}/02_ichorCNA"
CLASSIFICATION_DIR="${OUTPUT_DIR}/03_classification"
RISK_DIR="${OUTPUT_DIR}/04_risk_stratification"
REPORT_DIR="${OUTPUT_DIR}/05_report"

# STEP 1: Preprocessing
echo "[STEP 1/6] Preprocessing..."
bash ${SCRIPT_DIR}/01_preprocessing.sh \
  ${SAMPLE_ID} \
  ${FASTQ_R1} \
  ${FASTQ_R2} \
  ${REFERENCE} \
  ${PREPROCESS_DIR} \
  ${THREADS}

# STEP 2: ichorCNA
echo ""
echo "[STEP 2/6] Running ichorCNA..."
Rscript ${SCRIPT_DIR}/02_run_ichorCNA.R \
  --id ${SAMPLE_ID} \
  --bam ${PREPROCESS_DIR}/bam/${SAMPLE_ID}.sorted.bam \
  --outdir ${ICHORCNA_DIR} \
  --ichorcna ${ICHORCNA_PATH} \
  --threads ${THREADS}

# STEP 3: MB Subgroup Classification
echo ""
echo "[STEP 3/6] Classifying MB subgroup..."
mkdir -p ${CLASSIFICATION_DIR}

Rscript ${SCRIPT_DIR}/03_classify_mb_subgroup.R \
  --id ${SAMPLE_ID} \
  --seg ${ICHORCNA_DIR}/${SAMPLE_ID}.seg.txt \
  --params ${ICHORCNA_DIR}/${SAMPLE_ID}.params.txt \
  --config ${CONFIG_DIR}/mb_cnv_signatures.yaml \
  --output ${CLASSIFICATION_DIR}

# STEP 4: Risk Stratification
echo ""
echo "[STEP 4/6] Risk stratification..."
mkdir -p ${RISK_DIR}

Rscript ${SCRIPT_DIR}/05_risk_stratification.R \
  --id ${SAMPLE_ID} \
  --classification ${CLASSIFICATION_DIR}/${SAMPLE_ID}_classification.csv \
  --seg ${ICHORCNA_DIR}/${SAMPLE_ID}.seg.txt \
  --output ${RISK_DIR}/${SAMPLE_ID}_risk.csv \
  --metastasis ${METASTASIS} \
  --residual ${RESIDUAL} \
  --histology ${HISTOLOGY}

# STEP 5: Generate Clinical Report
echo ""
echo "[STEP 5/6] Generating clinical report..."
mkdir -p ${REPORT_DIR}

Rscript -e "rmarkdown::render(
  '${SCRIPT_DIR}/06_generate_report.Rmd',
  output_file='${REPORT_DIR}/${SAMPLE_ID}_MB_Report.pdf',
  params=list(
    sample_id='${SAMPLE_ID}',
    classification_file='${CLASSIFICATION_DIR}/${SAMPLE_ID}_classification.csv',
    risk_file='${RISK_DIR}/${SAMPLE_ID}_risk.csv',
    ichorcna_dir='${ICHORCNA_DIR}',
    seg_file='${ICHORCNA_DIR}/${SAMPLE_ID}.seg.txt'
  )
)"

# STEP 6: Generate Summary
echo ""
echo "[STEP 6/6] Pipeline summary..."

SUBGROUP=$(awk -F',' 'NR==2 {print $4}' ${CLASSIFICATION_DIR}/${SAMPLE_ID}_classification.csv)
CONFIDENCE=$(awk -F',' 'NR==2 {print $5}' ${CLASSIFICATION_DIR}/${SAMPLE_ID}_classification.csv)
RISK=$(awk -F',' 'NR==2 {print $3}' ${RISK_DIR}/${SAMPLE_ID}_risk.csv)
TF=$(awk -F',' 'NR==2 {print $2}' ${CLASSIFICATION_DIR}/${SAMPLE_ID}_classification.csv)

echo ""
echo "=========================================="
echo "PIPELINE COMPLETE"
echo "=========================================="
echo "Sample: ${SAMPLE_ID}"
echo "MB Subgroup: ${SUBGROUP}"
echo "Confidence: ${CONFIDENCE}"
echo "Risk Group: ${RISK}"
echo "Tumor Fraction: $(awk "BEGIN {printf \"%.2f%%\", ${TF}*100}")"
echo ""
echo "Results:"
echo "  Classification: ${CLASSIFICATION_DIR}/${SAMPLE_ID}_classification.csv"
echo "  Risk: ${RISK_DIR}/${SAMPLE_ID}_risk.csv"
echo "  Report: ${REPORT_DIR}/${SAMPLE_ID}_MB_Report.pdf"
echo ""
echo "End: $(date)"
