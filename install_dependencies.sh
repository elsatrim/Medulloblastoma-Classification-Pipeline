#!/bin/bash
# FIXED Installation for MB Classification Pipeline (macOS compatible)
# Based on Liu 2021, Markowitz 2025, Escudero 2020

echo "=== Installing MB Classification Pipeline (macOS) ==="

# Initialize conda if needed
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda not found. Please install Miniconda or Anaconda first."
    echo "Download from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Initialize conda for bash (if not already done)
conda init bash
source ~/.bash_profile 2>/dev/null || source ~/.bashrc 2>/dev/null

# Remove existing environment if it exists
conda env remove -n mb_classifier -y 2>/dev/null || true

# Create fresh environment with FLEXIBLE versioning (macOS compatible)
echo "Creating conda environment..."
conda create -n mb_classifier -y \
  python=3.9 \
  -c conda-forge

# Activate environment
eval "$(conda shell.bash hook)"
conda activate mb_classifier

# Install R separately (macOS compatible)
echo "Installing R..."
conda install -y -c conda-forge r-base=4.3

# Install bioinformatics tools (with flexible versions for macOS)
echo "Installing bioinformatics tools..."
conda install -y -c bioconda -c conda-forge \
  bwa \
  samtools \
  fastqc \
  multiqc \
  bedtools

# HMMcopy (may need special handling on macOS)
echo "Installing HMMcopy utilities..."
conda install -y -c bioconda hmmcopy_utils || {
  echo "WARNING: hmmcopy_utils failed to install via conda"
  echo "Will install from source later if needed"
}

# Install Python packages
echo "Installing Python packages..."
pip install --no-cache-dir \
  pandas==1.5.3 \
  numpy==1.24.3 \
  scikit-learn==1.3.0 \
  matplotlib==3.7.1 \
  seaborn==0.12.2 \
  pyyaml

# Install R packages
echo "Installing R packages (this may take 10-15 minutes)..."

Rscript -e "
options(repos = c(CRAN = 'https://cloud.r-project.org'))

# Core packages
install.packages(c(
  'optparse',
  'data.table',
  'ggplot2',
  'dplyr',
  'reshape2',
  'RColorBrewer',
  'scales',
  'knitr',
  'rmarkdown',
  'gridExtra',
  'pROC',
  'glmnet',
  'yaml'
), dependencies = TRUE, quiet = TRUE)

# Bioconductor packages
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager', quiet = TRUE)

BiocManager::install(c(
  'GenomeInfoDb',
  'GenomicRanges'
), update = FALSE, ask = FALSE)

cat('\nR packages installed successfully!\n')
"

# Install ichorCNA
echo "Installing ichorCNA..."
cd ~
if [ -d "ichorCNA" ]; then
  echo "ichorCNA directory exists, updating..."
  cd ichorCNA
  git pull
else
  git clone https://github.com/broadinstitute/ichorCNA.git
  cd ichorCNA
fi

# Install ichorCNA R package
Rscript -e "
if (!requireNamespace('devtools', quietly = TRUE)) {
  install.packages('devtools', repos = 'https://cloud.r-project.org')
}
devtools::install('.', dependencies = TRUE, upgrade = 'never')
"

# Clone Markowitz et al. fragmentomics code (reference only)
echo "Cloning Markowitz MB classification code (for reference)..."
cd ~
if [ ! -d "project-mb-lbseq-classification" ]; then
  git clone https://github.com/amarkowitzchla/project-mb-lbseq-classification.git
else
  echo "Markowitz code already cloned"
fi

# Create reference directory structure
echo "Setting up reference directories..."
mkdir -p ~/mb_reference

echo ""
echo "=========================================="
echo "Installation Complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Activate environment:"
echo "   conda activate mb_classifier"
echo ""
echo "2. Download reference genome manually:"
echo "   cd ~/mb_reference"
echo "   wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
echo "   gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
echo "   bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa"
echo "   samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa"
echo ""
echo "3. Test installation:"
echo "   conda activate mb_classifier"
echo "   bwa"
echo "   samtools"
echo "   Rscript -e 'library(ichorCNA)'"
echo ""  r-rcolorbrewer \
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
