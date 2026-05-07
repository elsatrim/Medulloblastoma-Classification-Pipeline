#!/bin/bash
# MB Classification Pipeline
# Supports both paired-end and single-end sequencing

set -e

# Show usage
show_usage() {
    cat << EOF
MB CLASSIFICATION PIPELINE
==========================

Usage:
  $0 <sample_id> <fastq_files...> <output_dir> [options]

Paired-end example:
  $0 SAMPLE1 reads_R1.fq.gz reads_R2.fq.gz output/ --threads 8

Single-end example:
  $0 SAMPLE2 reads.fq.gz output/ --threads 8 --single-end

Options:
  --threads N              CPU threads (default: 8)
  --single-end            Single-end mode (default: paired-end)
  --no-fragmentomics      Skip fragmentomics analysis
  --metastasis yes/no     Metastasis status (default: unknown)
  --residual yes/no       Residual disease (default: unknown)
  --histology TYPE        Histology: LCA/classic/DN (default: unknown)

Full example:
  $0 MB001 R1.fq.gz R2.fq.gz results/MB001 \\
     --threads 8 \\
     --metastasis no \\
     --residual no \\
     --histology classic

For help:
  $0 --help

EOF
}

# Parse arguments
if [ "$#" -lt 3 ] || [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    show_usage
    exit 0
fi

SAMPLE_ID=$1
shift

# Detect SE vs PE based on number of FASTQ files
FASTQ_FILES=()
while [[ $# -gt 0 ]] && [[ ! $1 == --* ]] && [[ ! -d $1 ]]; do
    FASTQ_FILES+=("$1")
    shift
done

OUTPUT_DIR=$1
shift

# Defaults
THREADS=8
SINGLE_END=false
RUN_FRAGMENTOMICS=true
METASTASIS="unknown"
RESIDUAL="unknown"
HISTOLOGY="unknown"

# Auto-detect SE/PE
if [ ${#FASTQ_FILES[@]} -eq 1 ]; then
    SINGLE_END=true
    FASTQ_R1="${FASTQ_FILES[0]}"
    FASTQ_R2=""
elif [ ${#FASTQ_FILES[@]} -eq 2 ]; then
    FASTQ_R1="${FASTQ_FILES[0]}"
    FASTQ_R2="${FASTQ_FILES[1]}"
else
    echo "ERROR: Expected 1 (SE) or 2 (PE) FASTQ files, got ${#FASTQ_FILES[@]}"
    exit 1
fi

# Parse options
while [[ $# -gt 0 ]]; do
    case $1 in
        --threads) THREADS="$2"; shift 2 ;;
        --single-end) SINGLE_END=true; shift ;;
        --no-fragmentomics) RUN_FRAGMENTOMICS=false; shift ;;
        --metastasis) METASTASIS="$2"; shift 2 ;;
        --residual) RESIDUAL="$2"; shift 2 ;;
        --histology) HISTOLOGY="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Paths (modify for your system)
REFERENCE="${HOME}/mb_reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
ICHORCNA_PATH="${HOME}/ichorCNA"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config"

# Validation
if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference genome not found: $REFERENCE"
    echo "Run: bash setup_reference.sh"
    exit 1
fi

echo "=========================================="
echo "MB CLASSIFICATION PIPELINE"
echo "=========================================="
echo "Sample: ${SAMPLE_ID}"
echo "Mode: $([ "$SINGLE_END" = true ] && echo "Single-end" || echo "Paired-end")"
echo "FASTQ: ${FASTQ_R1} $([ -n "$FASTQ_R2" ] && echo "$FASTQ_R2")"
echo "Threads: ${THREADS}"
echo "Fragmentomics: ${RUN_FRAGMENTOMICS}"
echo "Start: $(date)"
echo ""

# Create directories
mkdir -p ${OUTPUT_DIR}
PREPROCESS_DIR="${OUTPUT_DIR}/01_preprocessing"
ICHORCNA_DIR="${OUTPUT_DIR}/02_ichorCNA"
FRAGMENTOMICS_DIR="${OUTPUT_DIR}/03_fragmentomics"
CLASSIFICATION_DIR="${OUTPUT_DIR}/04_classification"
RISK_DIR="${OUTPUT_DIR}/05_risk_stratification"
REPORT_DIR="${OUTPUT_DIR}/06_report"

# Run pipeline steps...
# (rest of your pipeline code here - same as before)

echo ""
echo "=========================================="
echo "PIPELINE COMPLETE"
echo "=========================================="
echo "Results in: ${OUTPUT_DIR}"
echo "End: $(date)"
