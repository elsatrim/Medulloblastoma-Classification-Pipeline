#!/bin/bash
# Complete MB Classification Pipeline with Fragmentomics
# Based on Liu 2021 (CNV), Markowitz 2025 (fragmentomics), Escudero 2020 (signatures)

set -e

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <sample_id> <fastq_r1> <fastq_r2> <output_dir> [options]"
    echo ""
    echo "Options:"
    echo "  --threads N               Number of threads (default: 8)"
    echo "  --metastasis yes/no       Metastasis at diagnosis"
    echo "  --residual yes/no         Residual disease post-surgery"
    echo "  --histology LCA/classic/DN  Histology type"
    echo "  --no-fragmentomics        Skip fragmentomics analysis (CNV-only)"
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
RUN_FRAGMENTOMICS=true

while [[ $# -gt 0 ]]; do
    case $1 in
        --threads) THREADS="$2"; shift 2 ;;
        --metastasis) METASTASIS="$2"; shift 2 ;;
        --residual) RESIDUAL="$2"; shift 2 ;;
        --histology) HISTOLOGY="$2"; shift 2 ;;
        --no-fragmentomics) RUN_FRAGMENTOMICS=false; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# MODIFY THESE PATHS FOR YOUR SYSTEM
REFERENCE="${HOME}/mb_reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
ICHORCNA_PATH="${HOME}/ichorCNA"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config"

echo "=========================================="
echo "MB CLASSIFICATION PIPELINE v2.0"
echo "With Fragmentomics Support"
echo "=========================================="
echo "Sample: ${SAMPLE_ID}"
echo "Fragmentomics: ${RUN_FRAGMENTOMICS}"
echo "Start: $(date)"
echo ""

# Create directory structure
mkdir -p ${OUTPUT_DIR}
PREPROCESS_DIR="${OUTPUT_DIR}/01_preprocessing"
ICHORCNA_DIR="${OUTPUT_DIR}/02_ichorCNA"
FRAGMENTOMICS_DIR="${OUTPUT_DIR}/03_fragmentomics"
CLASSIFICATION_DIR="${OUTPUT_DIR}/04_classification"
RISK_DIR="${OUTPUT_DIR}/05_risk_stratification"
REPORT_DIR="${OUTPUT_DIR}/06_report"

# STEP 1: Preprocessing
echo "[STEP 1/7] Preprocessing..."
bash ${SCRIPT_DIR}/01_preprocessing.sh \
  ${SAMPLE_ID} \
  ${FASTQ_R1} \
  ${FASTQ_R2} \
  ${REFERENCE} \
  ${PREPROCESS_DIR} \
  ${THREADS}

BAM_FILE="${PREPROCESS_DIR}/bam/${SAMPLE_ID}.sorted.bam"

# STEP 2: ichorCNA
echo ""
echo "[STEP 2/7] Running ichorCNA..."
Rscript ${SCRIPT_DIR}/02_run_ichorCNA.R \
  --id ${SAMPLE_ID} \
  --bam ${BAM_FILE} \
  --outdir ${ICHORCNA_DIR} \
  --ichorcna ${ICHORCNA_PATH} \
  --threads ${THREADS}

# STEP 3: Fragmentomics Analysis (NEW!)
echo ""
if [ "$RUN_FRAGMENTOMICS" = true ]; then
  echo "[STEP 3/7] Running fragmentomics analysis..."
  mkdir -p ${FRAGMENTOMICS_DIR}
  
  # Check if fragmentomics script exists
  if [ -f "${SCRIPT_DIR}/fragmentomics_analysis.py" ]; then
    python ${SCRIPT_DIR}/fragmentomics_analysis.py \
      ${BAM_FILE} \
      --reference ${REFERENCE} \
      --output-dir ${FRAGMENTOMICS_DIR}
    
    echo "  ✓ Fragmentomics complete"
    FRAGMENTOMICS_METRICS="${FRAGMENTOMICS_DIR}/${SAMPLE_ID}_fragmentomics_metrics.json"
  else
    echo "  ⚠ Fragmentomics script not found: ${SCRIPT_DIR}/fragmentomics_analysis.py"
    echo "  ⚠ Continuing with CNV-only classification"
    RUN_FRAGMENTOMICS=false
    FRAGMENTOMICS_METRICS=""
  fi
else
  echo "[STEP 3/7] Fragmentomics analysis skipped (--no-fragmentomics)"
  FRAGMENTOMICS_METRICS=""
fi

# STEP 4: MB Subgroup Classification
echo ""
echo "[STEP 4/7] Classifying MB subgroup..."
mkdir -p ${CLASSIFICATION_DIR}

Rscript ${SCRIPT_DIR}/03_classify_mb_subgroup.R \
  --id ${SAMPLE_ID} \
  --seg ${ICHORCNA_DIR}/${SAMPLE_ID}.seg.txt \
  --params ${ICHORCNA_DIR}/${SAMPLE_ID}.params.txt \
  --config ${CONFIG_DIR}/mb_cnv_signatures.yaml \
  --fragmentomics ${FRAGMENTOMICS_METRICS} \
  --output ${CLASSIFICATION_DIR}

# STEP 5: Risk Stratification
echo ""
echo "[STEP 5/7] Risk stratification..."
mkdir -p ${RISK_DIR}

Rscript ${SCRIPT_DIR}/05_risk_stratification.R \
  --id ${SAMPLE_ID} \
  --classification ${CLASSIFICATION_DIR}/${SAMPLE_ID}_classification.csv \
  --seg ${ICHORCNA_DIR}/${SAMPLE_ID}.seg.txt \
  --output ${RISK_DIR}/${SAMPLE_ID}_risk.csv \
  --metastasis ${METASTASIS} \
  --residual ${RESIDUAL} \
  --histology ${HISTOLOGY}

# STEP 6: Generate Clinical Report
echo ""
echo "[STEP 6/7] Generating clinical report..."
mkdir -p ${REPORT_DIR}

Rscript -e "rmarkdown::render(
  '${SCRIPT_DIR}/06_generate_report.Rmd',
  output_file='${REPORT_DIR}/${SAMPLE_ID}_MB_Report.pdf',
  params=list(
    sample_id='${SAMPLE_ID}',
    classification_file='${CLASSIFICATION_DIR}/${SAMPLE_ID}_classification.csv',
    risk_file='${RISK_DIR}/${SAMPLE_ID}_risk.csv',
    ichorcna_dir='${ICHORCNA_DIR}',
    seg_file='${ICHORCNA_DIR}/${SAMPLE_ID}.seg.txt',
    fragmentomics_dir='${FRAGMENTOMICS_DIR}',
    fragmentomics_enabled='${RUN_FRAGMENTOMICS}'
  )
)"

# STEP 7: Generate Summary
echo ""
echo "[STEP 7/7] Pipeline summary..."

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

if [ "$RUN_FRAGMENTOMICS" = true ] && [ -f "${FRAGMENTOMICS_METRICS}" ]; then
  SL_RATIO=$(python -c "import json; d=json.load(open('${FRAGMENTOMICS_METRICS}')); print(f\"{d['short_long_ratio']:.3f}\")")
  echo "S/L Ratio: ${SL_RATIO}"
  echo "Fragmentomics: ENABLED"
else
  echo "Fragmentomics: DISABLED"
fi

echo ""
echo "Results:"
echo "  Classification: ${CLASSIFICATION_DIR}/${SAMPLE_ID}_classification.csv"
echo "  Risk: ${RISK_DIR}/${SAMPLE_ID}_risk.csv"
echo "  Report: ${REPORT_DIR}/${SAMPLE_ID}_MB_Report.pdf"

if [ "$RUN_FRAGMENTOMICS" = true ]; then
  echo "  Fragmentomics Plot: ${FRAGMENTOMICS_DIR}/${SAMPLE_ID}_fragmentomics.png"
  echo "  Fragmentomics Metrics: ${FRAGMENTOMICS_METRICS}"
fi

echo ""
echo "End: $(date)"