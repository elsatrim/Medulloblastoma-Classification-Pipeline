#!/bin/bash
# Download and index reference genome (GRCh38)
# Run this AFTER install.sh completes

set -e

echo "======================================"
echo "REFERENCE GENOME SETUP"
echo "======================================"
echo ""
echo "This will download and index GRCh38 (~8.2 GB)"
echo "Estimated time: 20-30 minutes"
echo "Estimated space: ~15 GB total"
echo ""

# Prompt user
read -p "Continue? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Setup cancelled."
    exit 0
fi

# Create reference directory
mkdir -p ~/mb_reference
cd ~/mb_reference

# Check if already downloaded
if [ -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
    echo "✓ Reference genome already exists"
else
    echo "Downloading reference genome..."
    wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    
    echo "Extracting..."
    gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    
    echo "✓ Download complete"
fi

# Index with BWA
if [ -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa.bwt" ]; then
    echo "✓ BWA index exists"
else
    echo "Creating BWA index (this takes ~30 minutes)..."
    bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
    echo "✓ BWA index complete"
fi

# Index with samtools
if [ -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai" ]; then
    echo "✓ SAMtools index exists"
else
    echo "Creating SAMtools index..."
    samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
    echo "✓ SAMtools index complete"
fi

echo ""
echo "======================================"
echo "REFERENCE GENOME READY!"
echo "======================================"
echo "Location: ~/mb_reference/"
echo "File: Homo_sapiens.GRCh38.dna.primary_assembly.fa"
echo ""
echo "You can now run the pipeline:"
echo "  bash run_mb_pipeline.sh --help"
echo "======================================"
