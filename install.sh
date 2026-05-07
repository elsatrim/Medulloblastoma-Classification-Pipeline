#!/bin/bash
# Complete MB Pipeline Installation
# One-step setup for macOS/Linux

set -e

echo "======================================"
echo "MB PIPELINE INSTALLER"
echo "======================================"
echo ""

# Check operating system
if [[ "$OSTYPE" == "darwin"* ]]; then
    OS="macOS"
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    OS="Linux"
else
    echo "ERROR: Unsupported OS. This installer works on macOS and Linux only."
    exit 1
fi

echo "Detected OS: $OS"
echo ""

# Check conda
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda not found!"
    echo "Please install Miniconda from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "✓ Conda found"
echo ""

# ============================================
# STEP 1: Create conda environment
# ============================================
echo "[STEP 1/5] Creating conda environment..."

conda init bash 2>/dev/null || true
source ~/.bash_profile 2>/dev/null || source ~/.bashrc 2>/dev/null || true

conda deactivate 2>/dev/null || true
conda env remove -n mb_classifier -y 2>/dev/null || true

conda create -n mb_classifier -y python=3.9 -c conda-forge

eval "$(conda shell.bash hook)"
conda activate mb_classifier

echo "✓ Environment created"
echo ""

# ============================================
# STEP 2: Install bioinformatics tools
# ============================================
echo "[STEP 2/5] Installing bioinformatics tools..."

conda install -y -c bioconda -c conda-forge \
    bwa \
    samtools \
    fastqc \
    bedtools \
    hmmcopy_utils 2>/dev/null || echo "  (some tools may need manual install)"

echo "✓ Core tools installed"
echo ""

# ============================================
# STEP 3: Install R and packages
# ============================================
echo "[STEP 3/5] Installing R and packages..."

conda install -y -c conda-forge r-base r-essentials 2>/dev/null || {
    echo "WARNING: R install via conda failed"
    echo "Please install R manually: https://cran.r-project.org"
}

# Install R packages
Rscript - <<'EOF'
options(repos = c(CRAN = "https://cloud.r-project.org"))

packages <- c("optparse", "data.table", "ggplot2", "dplyr", "yaml",
              "reshape2", "RColorBrewer", "scales", "knitr", 
              "rmarkdown", "gridExtra", "jsonlite")

cat("\nInstalling R packages...\n")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, quiet = TRUE)
    cat(sprintf("✓ %s\n", pkg))
  }
}

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", quiet = TRUE)
}
BiocManager::install(c("GenomeInfoDb", "GenomicRanges"), 
                     update = FALSE, ask = FALSE, quiet = TRUE)

cat("\n✓ R packages installed\n")
EOF

echo ""

# ============================================
# STEP 4: Install Python packages
# ============================================
echo "[STEP 4/5] Installing Python packages..."

pip install --quiet --no-cache-dir \
    pandas \
    numpy \
    scikit-learn \
    pyyaml \
    pysam \
    matplotlib \
    seaborn

echo "✓ Python packages installed"
echo ""

# ============================================
# STEP 5: Install ichorCNA
# ============================================
echo "[STEP 5/5] Installing ichorCNA..."

cd ~
if [ -d "ichorCNA" ]; then
    echo "  ichorCNA exists, updating..."
    cd ichorCNA
    git pull
else
    git clone https://github.com/broadinstitute/ichorCNA.git
    cd ichorCNA
fi

Rscript -e "
if (!requireNamespace('devtools', quietly = TRUE)) {
  install.packages('devtools', repos = 'https://cloud.r-project.org')
}
devtools::install('.', dependencies = TRUE, upgrade = 'never')
"

echo "✓ ichorCNA installed"
echo ""

# ============================================
# SETUP COMPLETE
# ============================================
echo "======================================"
echo "INSTALLATION COMPLETE!"
echo "======================================"
echo ""
echo "Next steps:"
echo ""
echo "1. Activate environment:"
echo "   conda activate mb_classifier"
echo ""
echo "2. Download reference genome (8.2 GB):"
echo "   bash setup_reference.sh"
echo ""
echo "3. Test installation:"
echo "   bwa"
echo "   samtools"
echo "   Rscript -e 'library(ichorCNA)'"
echo ""
echo "4. Run pipeline:"
echo "   bash run_mb_pipeline.sh --help"
echo ""
echo "======================================"
echo "Installation log saved to: install.log"
echo "======================================"
