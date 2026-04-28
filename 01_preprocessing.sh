#!/bin/bash
# Preprocessing for MB Classification
# Based on Liu 2021, Markowitz 2025

set -e

SAMPLE_ID=$1
FASTQ_R1=$2
FASTQ_R2=$3
REFERENCE=$4
OUTPUT_DIR=$5
THREADS=${6:-8}

mkdir -p ${OUTPUT_DIR}/{qc,bam,logs}

echo "=== MB Classification Pipeline: Preprocessing ==="
echo "Sample: ${SAMPLE_ID}"
echo "Start: $(date)"

# Quality Control
echo "[1/6] FastQC..."
fastqc ${FASTQ_R1} ${FASTQ_R2} \
  -o ${OUTPUT_DIR}/qc/ \
  -t ${THREADS} \
  2>&1 | tee ${OUTPUT_DIR}/logs/01_fastqc.log

# Alignment (BWA-MEM - as per Liu 2021)
echo "[2/6] BWA-MEM alignment..."
bwa mem \
  -t ${THREADS} \
  -M \
  -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:MB_CSF" \
  ${REFERENCE} \
  ${FASTQ_R1} \
  ${FASTQ_R2} \
  2> ${OUTPUT_DIR}/logs/02_bwa.log | \
  samtools view -Sb - | \
  samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/bam/${SAMPLE_ID}.sorted.bam -

echo "[3/6] Indexing BAM..."
samtools index ${OUTPUT_DIR}/bam/${SAMPLE_ID}.sorted.bam

# Statistics
echo "[4/6] Calculating statistics..."
samtools flagstat ${OUTPUT_DIR}/bam/${SAMPLE_ID}.sorted.bam \
  > ${OUTPUT_DIR}/bam/${SAMPLE_ID}.flagstat.txt

samtools stats ${OUTPUT_DIR}/bam/${SAMPLE_ID}.sorted.bam \
  > ${OUTPUT_DIR}/bam/${SAMPLE_ID}.stats.txt

# Coverage calculation
echo "[5/6] Coverage analysis..."
samtools depth ${OUTPUT_DIR}/bam/${SAMPLE_ID}.sorted.bam | \
  awk '{sum+=$3; count++} END {
    print "Average_coverage: " sum/count
    print "Total_bases: " count
  }' > ${OUTPUT_DIR}/bam/${SAMPLE_ID}.coverage.txt

# Extract coverage value for QC
AVG_COV=$(grep "Average_coverage" ${OUTPUT_DIR}/bam/${SAMPLE_ID}.coverage.txt | cut -d' ' -f2)
echo "Average coverage: ${AVG_COV}×"

# QC check (Liu 2021 used ≥1× coverage)
if (( $(echo "$AVG_COV < 0.5" | bc -l) )); then
    echo "WARNING: Coverage below 0.5× (${AVG_COV}×)"
    echo "Quality may be compromised."
fi

echo "[6/6] Preprocessing complete"
echo "End: $(date)"
