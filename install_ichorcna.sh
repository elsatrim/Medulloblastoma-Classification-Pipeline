# Activate environment
conda activate mb_classifier

# Install devtools now (libgit2 is installed, so this should work)
Rscript -e 'install.packages("devtools", repos="https://cloud.r-project.org")'

# Install HMMcopy from Bioconductor
Rscript -e '
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")
BiocManager::install("HMMcopy")
'

# Clone ichorCNA
cd ~
git clone https://github.com/broadinstitute/ichorCNA.git

# Install ichorCNA
cd ~/ichorCNA
Rscript -e 'devtools::install(".", dependencies=TRUE, upgrade="never")'

# Test installation
Rscript -e 'library(ichorCNA); cat("✓ ichorCNA loaded successfully!\n")'
