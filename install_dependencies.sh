#!/bin/bash
# Installation for MB Classification Pipeline
# Based on Liu 2021, Markowitz 2025, Escudero 2020

echo "=== Installing MB Classification Pipeline ==="

# Create conda environment
conda create -n mb_classifier -y python=3.9 r-base=4.2

conda activate mb_classifier

# Install bioinformatics tools
echo "Installing bioinformatics tools..."
conda install -c bioconda -y \
  bwa=0.7.17 \
  samtools=1.17 \
  fastqc=0.12.1 \
  multiqc=1.14 \
  hmmcopy_utils \
  bedtools=2.30.0

# Install Python packages (for Markowitz fragmentomics)
pip install pandas==1.5.3 \
  numpy==1.24.3 \
  scikit-learn==1.3.0 \
  matplotlib==3.7.1 \
  seaborn==0.12.2

# Install R packages
echo "Installing R packages..."
conda install -c conda-forge -c bioconda -y \
  r-optparse \
  r-data.table \
  r-ggplot2 \
  r-dplyr \
  r-reshape2 \
  r-rcolorbrewer \
  r-scales \
  r-caret \
  r-randomforest \
  r-e1071

# Additional R packages
Rscript -e "install.packages(c('knitr', 'rmarkdown', 'gridExtra', 'pROC', 'glmnet'), repos='http://cran.us.r-project.org')"

# Install ichorCNA
echo "Installing ichorCNA..."
cd ~
git clone https://github.com/broadinstitute/ichorCNA.git

# Clone Markowitz et al. fragmentomics code
echo "Cloning Markowitz MB classification code..."
cd ~
git clone https://github.com/amarkowitzchla/project-mb-lbseq-classification.git

# Download reference genome
echo "Downloading reference genome..."
mkdir -p ~/mb_reference
cd ~/mb_reference

# GRCh38
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Index
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

echo "=== Installation complete! ==="
echo "Activate with: conda activate mb_classifier"
