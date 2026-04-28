# Medulloblastoma CSF Liquid Biopsy Classification Pipeline


---

## 📖 Table of Contents

1. [Overview](#overview)
2. [Scientific Background](#scientific-background)
3. [Pipeline Architecture](#pipeline-architecture)
4. [Installation](#installation)
5. [Quick Start](#quick-start)
6. [Detailed Usage](#detailed-usage)
7. [Output Files](#output-files)
8. [Interpretation Guide](#interpretation-guide)
9. [Clinical Applications](#clinical-applications)
10. [Troubleshooting](#troubleshooting)
11. [References](#references)
12. [Citation](#citation)

---

## 🎯 Overview

### What is this pipeline?

This bioinformatics pipeline classifies **medulloblastoma (MB)** tumors into molecular subgroups using **cell-free DNA (cfDNA)** from **cerebrospinal fluid (CSF)** samples.

**Input:** CSF cfDNA sequencing data (low-coverage whole genome sequencing, ~0.5-1× coverage)

**Output:** 
- Molecular subgroup (WNT, SHH, Group 3, Group 4)
- Risk stratification (Favourable, Standard, High, Very High)
- Clinical report (PDF)
- Copy number profile
- Quality metrics

**Timeline:** 2-3 days from CSF sample to clinical report

---

### Why use cfDNA instead of tissue?

**Advantages:**
- ✅ **Non-invasive:** CSF routinely collected during EVD placement (no additional procedure)
- ✅ **Pre-operative:** Results available BEFORE surgery
- ✅ **Fast:** 2-3 days vs. 7-10 days for tissue methylation
- ✅ **Serial monitoring:** Track tumor burden over time
- ✅ **MRD detection:** Detect relapse 3+ months before imaging

**Clinical impact:**
- Inform surgical approach (WNT → conservative resection)
- Guide treatment intensity (WNT → reduce therapy; Group 3 → intensify)
- Early family counseling (prognosis from 40% to 95% depending on subgroup)

---

### Who should use this pipeline?

**Intended users:**
- 🏥 Pediatric oncology centers
- 🧬 Molecular diagnostics laboratories
- 🔬 Cancer research groups
- 📊 Bioinformatics core facilities

**Requirements:**
- Experience with NGS data analysis
- Access to low-coverage WGS (Illumina or similar)
- Computational resources (8+ CPU cores, 16GB+ RAM)
- Clinical collaboration (neurosurgeons, oncologists, pathologists)

---

## 🧬 Scientific Background

### Medulloblastoma Molecular Subgroups

Medulloblastoma has **4 distinct molecular subgroups**, each with different:
- Biology
- Prognosis
- Treatment response
- Survival

| Subgroup | Frequency | Key Alterations | 5-Year Survival | Treatment |
|----------|-----------|-----------------|-----------------|-----------|
| **WNT** | 10% | Monosomy 6, CTNNB1 mutation | >90% | **Reduce** intensity |
| **SHH** | 30% | 9q loss, PTCH1/SMO mutations | 50-85% (variable) | Standard ± SMO inhibitors |
| **Group 3** | 25% | MYC amplification, iso17q | <50% | **Intensify** therapy |
| **Group 4** | 35% | Isochromosome 17q, chr8 loss | ~75% | Standard |

**Source:** Schwalbe et al. 2017, Lancet Oncology

---

### Copy Number Alterations (CNV) Signatures

Each subgroup has characteristic CNV patterns detectable by low-coverage WGS:

#### WNT Subgroup
```
Diagnostic features:
├─ Monosomy 6 (100% of cases) 
├─ Chromosome 7 gain (80%)
└─ CTNNB1 exon 3 mutation (requires sequencing)
```

#### SHH Subgroup
```
Characteristic features:
├─ Chromosome 9q loss (PTCH1 region) (32-100%)
├─ Chromosome 10q loss (32%)
├─ 17p loss (TP53) (20%)
├─ MYCN amplification (7%)
└─ GLI2 amplification
```

#### Group 3
```
Characteristic features:
├─ Isochromosome 17q (42%)
├─ MYC amplification (17%) ⚠️ HIGH RISK
├─ OTX2 amplification
└─ GFI1/GFI1B activation
```

#### Group 4
```
Characteristic features:
├─ Isochromosome 17q (67%) MOST COMMON
├─ Chromosome 8 loss (20%)
├─ Chromosome X loss (females)
├─ MYCN amplification (4%)
└─ CDK6 amplification
```

**Sources:** 
- Escudero et al. 2020, Nature Communications
- Schwalbe et al. 2017, Lancet Oncology

---

### Published Evidence for CSF cfDNA

This pipeline is based on validated methods from:

**1. Liu et al. 2021 (Cancer Cell)**
- Serial CSF liquid biopsy for MB
- lcWGS (0.5× coverage) + ichorCNA for CNV detection
- MRD detection with **3.8 month lead time** before imaging
- https://doi.org/10.1016/j.ccell.2021.08.012

**2. Markowitz et al. 2025 (NPJ Precision Oncology)**
- Fragmentomics + CNV for MB classification
- Meta-classifier AUC 0.94
- Tumor fraction ≥8% threshold
- https://doi.org/10.1038/s41698-025-01067-5

**3. Escudero et al. 2020 (Nature Communications)**
- CSF ctDNA detects 98.9% of tumor mutations
- 100% concordance with tissue for subgroup classification
- CSF superior to plasma for brain tumors
- https://doi.org/10.1038/s41467-020-19175-0

**4. Christodoulou et al. 2023 (NPJ Precision Oncology)**
- LBSeq4Kids clinical assay (CHLA)
- Combined lcWGS + targeted sequencing
- Clinical validation in pediatric solid tumors
- https://doi.org/10.1038/s41698-023-00361-0

---

## 🏗️ Pipeline Architecture

### Workflow Overview

```
CSF Sample (from EVD)
    ↓
DNA Extraction
    ↓
Library Preparation
    ↓
Low-Coverage WGS (0.5-1×)
    ↓
FASTQ Files
    ↓
┌─────────────────────────────────────────┐
│  BIOINFORMATICS PIPELINE                │
│                                         │
│  [1] Preprocessing                      │
│      ├─ Quality Control (FastQC)       │
│      ├─ Alignment (BWA-MEM)            │
│      ├─ BAM Processing (SAMtools)      │
│      └─ Coverage Calculation           │
│                                         │
│  [2] CNV Detection (ichorCNA)          │
│      ├─ Tumor Fraction Estimation     │
│      ├─ Copy Number Calling            │
│      └─ Ploidy Estimation              │
│                                         │
│  [3] Subgroup Classification           │
│      ├─ CNV Signature Matching        │
│      ├─ Decision Tree Logic           │
│      └─ Confidence Assessment          │
│                                         │
│  [4] Risk Stratification               │
│      ├─ Schwalbe 2017 Criteria        │
│      ├─ Prognostic Features           │
│      └─ Survival Estimates            │
│                                         │
│  [5] Clinical Report Generation        │
│      ├─ CNV Visualization             │
│      ├─ Quality Metrics               │
│      └─ Treatment Recommendations      │
│                                         │
└─────────────────────────────────────────┘
    ↓
PDF Clinical Report
```

---

### Technical Components

**Core Tools:**
- **BWA-MEM v0.7.17:** Read alignment
- **SAMtools v1.17:** BAM processing
- **ichorCNA:** CNV detection + tumor fraction estimation
- **FastQC:** Quality control
- **R/Bioconductor:** Data analysis and visualization
- **RMarkdown:** Report generation

**Analysis Method:**
- **Simple presence/absence logic** (matches Liu et al. 2021)
- Based on diagnostic CNV features from published literature
- Transparent, clinically interpretable
- No machine learning required

---

### Directory Structure

```
mb-classifier/
├── README.md                          # This file
├── install_dependencies_fixed.sh      # Installation script (macOS/Linux)
├── install_ichorcna.sh                # ichorCNA installation
├── run_mb_pipeline.sh                 # Main pipeline wrapper
│
├── config/
│   ├── mb_cnv_signatures.yaml         # Subgroup CNV definitions
│   ├── reference_paths.yaml           # Reference genome paths
│   └── clinical_thresholds.yaml       # QC thresholds
│
├── scripts/
│   ├── 01_preprocessing.sh            # FastQC, BWA, SAMtools
│   ├── 02_run_ichorCNA.R              # CNV detection
│   ├── 03_classify_mb_subgroup_v2.R   # Classification (simple logic)
│   ├── 04_risk_stratification.R       # Risk stratification
│   └── 05_generate_report.Rmd         # Clinical report
│
├── reference_data/
│   ├── Homo_sapiens.GRCh38.fa         # Reference genome
│   ├── Homo_sapiens.GRCh38.fa.bwt     # BWA index
│   └── ...                            # Other reference files
│
├── data/
│   └── test_samples/                  # Example data
│
└── results/
    └── [SAMPLE_ID]/
        ├── 01_preprocessing/
        ├── 02_ichorCNA/
        ├── 03_classification/
        ├── 04_risk_stratification/
        └── 05_report/
```

---

## 💻 Installation

### System Requirements

**Hardware:**
- CPU: 8+ cores recommended (minimum 4)
- RAM: 16GB minimum (32GB recommended)
- Storage: 50GB per sample (temporary), 10GB per sample (final)
- OS: macOS, Linux (Ubuntu 20.04+), or Windows WSL2

**Software Prerequisites:**
- Conda/Miniconda (https://docs.conda.io/en/latest/miniconda.html)
- Git
- Internet connection (for installation only)

---

### Installation Steps

#### Step 1: Clone Repository

```bash
# Clone the repository
git clone https://github.com/your-org/mb-classifier.git
cd mb-classifier

# Make scripts executable
chmod +x install_dependencies_fixed.sh
chmod +x install_ichorcna.sh
chmod +x run_mb_pipeline.sh
chmod +x scripts/*.sh
```

---

#### Step 2: Install Dependencies

**For macOS/Linux:**

```bash
# Run main installation script
bash install_dependencies_fixed.sh

# This will:
# - Create conda environment 'mb_classifier'
# - Install Python 3.9
# - Install R 4.3
# - Install bioinformatics tools (BWA, SAMtools, FastQC)
# - Install R packages
# - Clone ichorCNA repository

# Installation takes 10-20 minutes
```

**Troubleshooting installation issues:**

If conda installation fails, see [Troubleshooting](#troubleshooting) section.

---

#### Step 3: Install ichorCNA

```bash
# Activate environment
conda activate mb_classifier

# Install ichorCNA separately
bash install_ichorcna.sh

# This will:
# - Install HMMcopy dependencies
# - Install ichorCNA R package
```

---

#### Step 4: Download Reference Genome

```bash
# Create reference directory
mkdir -p ~/mb_reference
cd ~/mb_reference

# Download GRCh38 reference genome
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Decompress
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Index for BWA (takes ~60 minutes)
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Index for SAMtools (takes ~5 minutes)
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

**Alternative:** Use pre-indexed reference if available from your institution.

---

#### Step 5: Test Installation

```bash
# Activate environment
conda activate mb_classifier

# Test core tools
bwa
samtools --version
fastqc --version

# Test R packages
Rscript -e "library(ichorCNA); library(data.table); library(ggplot2)"

# Should all complete without errors
```

**Expected output:**
```
Program: bwa (alignment via Burrows-Wheeler transformation)
...
samtools 1.17
...
FastQC v0.12.1
...
```

---

## 🚀 Quick Start

### Minimal Example

Process a single sample with default parameters:

```bash
# Activate environment
conda activate mb_classifier

# Run pipeline
bash run_mb_pipeline.sh \
  MB001 \
  data/MB001_R1.fastq.gz \
  data/MB001_R2.fastq.gz \
  results/MB001 \
  --threads 8

# Timeline:
# - Preprocessing: 1-2 hours
# - ichorCNA: 1-2 hours
# - Classification: 5 minutes
# - Risk stratification: 2 minutes
# - Report generation: 3 minutes
# Total: ~4-6 hours
```

**Output:**
```
results/MB001/
├── 01_preprocessing/
│   ├── qc/                    # FastQC reports
│   ├── bam/                   # Aligned BAM file
│   └── logs/
├── 02_ichorCNA/
│   ├── MB001.seg.txt          # CNV segments
│   ├── MB001.params.txt       # Tumor fraction, ploidy
│   └── MB001_genomeWide.pdf   # CNV plot
├── 03_classification/
│   └── MB001_classification.csv   # Subgroup result
├── 04_risk_stratification/
│   └── MB001_risk.csv         # Risk group
└── 05_report/
    └── MB001_MB_Report.pdf    # Clinical report ⭐
```

---

## 📖 Detailed Usage

### Full Command Syntax

```bash
bash run_mb_pipeline.sh \
  <sample_id> \
  <fastq_r1> \
  <fastq_r2> \
  <output_dir> \
  [OPTIONS]
```

**Required Arguments:**
- `sample_id`: Unique sample identifier (e.g., MB001, CSF_20260428_01)
- `fastq_r1`: Path to R1 FASTQ file (gzipped or uncompressed)
- `fastq_r2`: Path to R2 FASTQ file
- `output_dir`: Directory for results

**Optional Arguments:**
- `--threads N`: Number of CPU threads (default: 8)
- `--metastasis yes|no|unknown`: Metastatic disease at diagnosis (default: unknown)
- `--residual yes|no|unknown`: Residual disease post-surgery (default: unknown)
- `--histology LCA|classic|DN|unknown`: Histology type (default: unknown)
  - LCA = Large cell/anaplastic (high risk)
  - DN = Desmoplastic/nodular
  - classic = Classic medulloblastoma

---

### Example with Clinical Parameters

```bash
# Patient with complete clinical information
bash run_mb_pipeline.sh \
  MB_Patient_2026_042 \
  /path/to/CSF_R1.fastq.gz \
  /path/to/CSF_R2.fastq.gz \
  results/MB_Patient_2026_042 \
  --threads 16 \
  --metastasis no \
  --residual no \
  --histology classic
```

This provides more accurate risk stratification by incorporating clinical features.

---

### Batch Processing Multiple Samples

Create a sample sheet (`samples.txt`):

```
sample_id,fastq_r1,fastq_r2,metastasis,residual,histology
MB001,data/MB001_R1.fq.gz,data/MB001_R2.fq.gz,no,no,classic
MB002,data/MB002_R1.fq.gz,data/MB002_R2.fq.gz,yes,no,LCA
MB003,data/MB003_R1.fq.gz,data/MB003_R2.fq.gz,no,yes,classic
```

Run batch script:

```bash
# Create batch processing script
cat > run_batch.sh << 'EOF'
#!/bin/bash

while IFS=',' read -r SAMPLE R1 R2 MET RES HIST; do
  # Skip header
  if [[ "$SAMPLE" == "sample_id" ]]; then continue; fi
  
  echo "Processing ${SAMPLE}..."
  
  bash run_mb_pipeline.sh \
    ${SAMPLE} \
    ${R1} \
    ${R2} \
    results/${SAMPLE} \
    --threads 8 \
    --metastasis ${MET} \
    --residual ${RES} \
    --histology ${HIST}
    
done < samples.txt
EOF

chmod +x run_batch.sh
bash run_batch.sh
```

---

### Serial Monitoring (MRD Detection)

For patients undergoing treatment, serial CSF sampling can detect minimal residual disease:

```bash
# Baseline (pre-treatment)
bash run_mb_pipeline.sh \
  MB001_Baseline \
  data/MB001_baseline_R1.fq.gz \
  data/MB001_baseline_R2.fq.gz \
  results/MB001_Baseline

# Month 1 (during treatment)
bash run_mb_pipeline.sh \
  MB001_Month1 \
  data/MB001_month1_R1.fq.gz \
  data/MB001_month1_R2.fq.gz \
  results/MB001_Month1

# Month 3
bash run_mb_pipeline.sh \
  MB001_Month3 \
  data/MB001_month3_R1.fq.gz \
  data/MB001_month3_R2.fq.gz \
  results/MB001_Month3

# Compare tumor fractions over time
grep "Tumor_Fraction" results/MB001_*/03_classification/*.csv
```

**Expected pattern:**

```
Baseline: 45% (high tumor burden)
Month 1:  12% (responding to treatment)
Month 3:  2%  (minimal residual disease)
```

**Clinical significance:**
- **Decreasing:** Good response to treatment
- **Stable/Increasing:** Possible treatment resistance or relapse
- **Undetectable → Reappearing:** Early relapse detection (3+ months before imaging)

---

## 📂 Output Files

### Complete Output Structure

```
results/SAMPLE_ID/
│
├── 01_preprocessing/
│   ├── qc/
│   │   ├── SAMPLE_ID_R1_fastqc.html
│   │   ├── SAMPLE_ID_R2_fastqc.html
│   │   └── multiqc_report.html
│   ├── bam/
│   │   ├── SAMPLE_ID.sorted.bam
│   │   ├── SAMPLE_ID.sorted.bam.bai
│   │   ├── SAMPLE_ID.flagstat.txt
│   │   ├── SAMPLE_ID.stats.txt
│   │   └── SAMPLE_ID.coverage.txt
│   └── logs/
│       ├── 01_fastqc.log
│       └── 02_bwa.log
│
├── 02_ichorCNA/
│   ├── SAMPLE_ID.seg.txt              # CNV segments 
│   ├── SAMPLE_ID.params.txt           # Tumor fraction, ploidy 
│   ├── SAMPLE_ID.cna.seg
│   ├── SAMPLE_ID_genomeWide.pdf       # CNV visualization
│   ├── SAMPLE_ID_genomeWide_all_sols.pdf
│   ├── SAMPLE_ID_ichorCNA_summary.csv
│   └── wig/
│       └── SAMPLE_ID.wig
│
├── 03_classification/
│   └── SAMPLE_ID_classification.csv   # Subgroup classification 
│
├── 04_risk_stratification/
│   └── SAMPLE_ID_risk.csv             # Risk stratification 
│
└── 05_report/
    └── SAMPLE_ID_MB_Report.pdf        # Clinical report 
```

---

### Key Output Files Explained

#### 1. Classification Results (`SAMPLE_ID_classification.csv`)

```csv
Sample_ID,Tumor_Fraction,Ploidy,MB_Subgroup,Confidence,Rationale,Monosomy_6,Chr9q_Loss,Iso17q,MYC_amp,MYCN_amp,Chr8_loss,Timestamp
MB001,0.32,2,WNT,High,"Monosomy 6 detected (diagnostic for WNT-activated MB)",TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,2026-04-28 14:23:15
```

**Columns:**
- `Tumor_Fraction`: Proportion of tumor DNA in CSF (0-1 scale)
- `MB_Subgroup`: WNT, SHH, Group 3, Group 4, or Unclassified
- `Confidence`: High, Medium, Low
- `Rationale`: Explanation for classification
- Feature columns: TRUE/FALSE for each diagnostic CNV

---

#### 2. Risk Stratification (`SAMPLE_ID_risk.csv`)

```csv
Sample_ID,MB_Subgroup,Risk_Group,PFS_5year,Rationale,MYC_amplification,MYCN_amplification,Chr13_loss,Metastasis,Residual_disease
MB001,WNT,Favourable,>90%,"WNT-activated MB has excellent prognosis",FALSE,FALSE,FALSE,FALSE,FALSE
```

**Columns:**
- `Risk_Group`: Favourable, Standard, High, Very High
- `PFS_5year`: Expected 5-year progression-free survival
- `Rationale`: Clinical explanation
- Risk modifiers: MYC/MYCN amplification, metastasis, etc.

---

#### 3. ichorCNA Parameters (`SAMPLE_ID.params.txt`)

```
Tumor_Fraction	Ploidy	Sample
0.32	2	MB001
```

**Key metrics:**
- `Tumor_Fraction`: Critical for QC (≥8% recommended)
- `Ploidy`: Tumor ploidy (diploid = 2, triploid = 3)

---

#### 4. Clinical Report (`SAMPLE_ID_MB_Report.pdf`)

**Contents:**
1. **Executive Summary:** Subgroup, risk, tumor fraction
2. **Classification:** Scores, detected features
3. **Risk Stratification:** Risk group, survival estimate, modifiers
4. **CNV Profile:** Genome-wide copy number plot
5. **Quality Metrics:** Coverage, tumor fraction, confidence
6. **Clinical Interpretation:** Subgroup-specific treatment implications
7. **Methods:** Analysis pipeline description
8. **Limitations:** Caveats and recommended confirmatory testing

**Example interpretation:**

> **WNT-activated medulloblastoma** has excellent prognosis (>90% 5-year survival).
>
> **Treatment Implications:**
> - Consider reduced-intensity therapy to minimize late effects
> - Gross total resection may not be necessary (survival excellent regardless)
> - Conservative surgical approach to reduce posterior fossa syndrome risk
> - Reduced craniospinal radiation dose may be appropriate

---

## 📊 Interpretation Guide

### Classification Confidence Levels

#### High Confidence ✅

**Criteria:**
- Diagnostic CNV detected (e.g., monosomy 6 for WNT)
- Tumor fraction ≥8%
- Clear CNV pattern matching one subgroup

**Action:** Use for clinical decision-making

---

#### Medium Confidence ⚠️

**Criteria:**
- Characteristic but non-diagnostic CNVs
- Tumor fraction 5-8%
- Pattern suggests multiple possible subgroups

**Action:** Consider confirmatory testing (tissue methylation)

**Example:**
> "Isochromosome 17q present but cannot distinguish between Group 3 and Group 4"

---

#### Low Confidence ❌

**Criteria:**
- No diagnostic CNVs detected
- Tumor fraction <5%
- Ambiguous pattern

**Action:** Repeat CSF sampling or perform tissue testing

**Possible causes:**
- Low tumor burden
- Early disease stage
- Sampling error
- Non-MB tumor

---

### Quality Control Checkpoints

#### PASS ✅

```
Tumor Fraction: 32%     ✅ >8%
Average Coverage: 0.8×  ✅ >0.5×
Alignment Rate: 94%     ✅ >80%
Classification: WNT     ✅ High confidence
```

**Action:** Proceed with clinical interpretation

---

#### WARNING ⚠️

```
Tumor Fraction: 6%      ⚠️ Below 8% threshold
Average Coverage: 0.7×  ✅ Adequate
Alignment Rate: 92%     ✅ Good
Classification: Group 4 ⚠️ Medium confidence
```

**Action:**
- Report results with caveat
- Consider repeat sampling
- Recommend tissue confirmation

---

#### FAIL ❌

```
Tumor Fraction: 2%      ❌ Too low
Average Coverage: 0.3×  ❌ Insufficient
Alignment Rate: 65%     ❌ Poor quality
Classification: Unclassified
```

**Action:**
- DO NOT use for clinical decisions
- Repeat sequencing with more input DNA
- Check sample quality

---

### Subgroup-Specific Interpretations

#### WNT Subgroup

**Key Features:**
- Monosomy 6 (diagnostic)
- Chromosome 7 gain
- Excellent prognosis (>90% survival)

**Clinical Implications:**
- **Surgery:** Conservative resection acceptable
- **Radiation:** Consider dose reduction (23.4 Gy vs 36 Gy)
- **Chemotherapy:** May de-escalate in clinical trials
- **Prognosis:** Excellent
- **Counseling:** Reassure family, focus on minimizing late effects

---

#### SHH Subgroup

**Key Features:**
- 9q loss (PTCH1)
- 10q loss
- MYCN amplification (poor prognosis)
- TP53 mutation (poor prognosis - requires sequencing)

**Clinical Implications:**
- **Surgery:** Maximal safe resection
- **Therapy:** Consider SMO inhibitors (vismodegib) for recurrent disease
- **Risk:** Variable (depends on TP53 status)
- **Counseling:** If MYCN+ or TP53+, discuss higher risk

---

#### Group 3

**Key Features:**
- Isochromosome 17q
- MYC amplification (high risk)
- OTX2 amplification

**Clinical Implications:**
- **Surgery:** Maximal safe resection critical
- **Therapy:** Intensified chemotherapy
- **MYC amplification:** Very high risk, consider experimental therapies
- **Prognosis:** Poor (<50% survival, worse if MYC+)
- **Counseling:** Prepare family for intensive treatment

---

#### Group 4

**Key Features:**
- Isochromosome 17q (most common)
- Chromosome 8 loss
- Chromosome X loss (females)
- MYCN or CDK6 amplification

**Clinical Implications:**
- **Surgery:** Standard approach
- **Therapy:** Standard risk-adapted
- **Chr 13 loss:** Protective (better prognosis)
- **Prognosis:** Intermediate (~75% survival)
- **Counseling:** Guardedly optimistic

---

## ✅ Validation & Quality Control

### Pre-Analytical QC

**CSF Sample Quality:**
- Volume: ≥1 mL recommended
- Appearance: Clear (bloody samples may reduce quality)
- Processing: <2 hours from collection to extraction
- Storage: -80°C if not processed immediately

**DNA Quality:**
- Concentration: ≥10 ng/μL recommended
- Fragment size: Peak 150-180 bp (cfDNA characteristic)
- 260/280 ratio: 1.8-2.0
- Quantification: Qubit or similar (accurate for low concentrations)

---

### Analytical QC

**Sequencing Metrics:**
- Total reads: ≥50 million paired-end reads
- Read length: ≥75 bp (100-150 bp optimal)
- Q30: ≥80%
- Coverage: ≥0.5× average
- Alignment rate: ≥80%

**ichorCNA QC:**
- Tumor fraction: ≥8% (Markowitz 2025 threshold)
- Ploidy estimate: Reasonable (2 or 3 typical)
- Segments: Should show clear CNV calls

---


## 🔧 Troubleshooting

### Installation Issues

#### Problem: Conda package conflicts

```
LibMambaUnsatisfiableError: Encountered problems while solving
```

**Solution:**
```bash
# Remove old environment
conda deactivate
conda env remove -n mb_classifier -y

# Create minimal environment
conda create -n mb_classifier python=3.9 -y
conda activate mb_classifier

# Install tools one at a time
conda install -c bioconda bwa -y
conda install -c bioconda samtools -y
# etc.
```

---

#### Problem: R package installation fails

```
Error: package 'X' is not available
```

**Solution:**
```bash
# Update R packages
conda activate mb_classifier
Rscript -e "update.packages(ask=FALSE)"

# Install manually
Rscript -e "install.packages('PACKAGE_NAME', repos='https://cloud.r-project.org')"
```

---

#### Problem: ichorCNA installation fails

```
Error: HMMcopy not found
```

**Solution:**
```bash
# Install HMMcopy from Bioconductor
Rscript -e "
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('HMMcopy')
"

# Then reinstall ichorCNA
cd ~/ichorCNA
Rscript -e "devtools::install('.', dependencies=TRUE)"
```

---

### Runtime Issues

#### Problem: Low coverage

```
Average coverage: 0.2×
```

**Causes:**
- Insufficient sequencing depth
- Low input DNA
- High adapter contamination

**Solutions:**
1. Increase sequencing depth (target 100M reads)
2. Check DNA concentration before library prep
3. Run FastQC to check adapter content

---

#### Problem: Low tumor fraction

```
Tumor fraction: 3%
Classification: Unclassified
```

**Causes:**
- Early disease stage
- Low tumor burden
- Sampling location (CSF may not be tumor-adjacent)
- Non-tumor sample

**Solutions:**
1. Repeat CSF sampling (different location if possible)
2. Increase sequencing depth
3. Perform tissue biopsy for confirmation
4. Consider alternative diagnosis

---

#### Problem: Pipeline crashes

```
Error: Cannot allocate memory
```

**Causes:**
- Insufficient RAM
- Too many parallel processes

**Solutions:**
```bash
# Reduce thread count
bash run_mb_pipeline.sh SAMPLE R1 R2 OUT --threads 4

# Close other applications
# Consider running on HPC cluster
```

---

#### Problem: ichorCNA fails to converge

```
Warning: ichorCNA did not converge
```

**Causes:**
- Very low coverage
- Highly aneuploid tumor
- Poor quality data

**Solutions:**
1. Check alignment rate (should be >80%)
2. Check coverage (should be >0.5×)
3. Try different ploidy parameters in ichorCNA
4. Visual inspection of WIG file

---

### Classification Issues

#### Problem: Ambiguous classification

```
MB_Subgroup: Group 3 or Group 4
Confidence: Medium
```

**Interpretation:**
- Both subgroups share iso17q
- Need additional features to distinguish
- Chr 8 loss → Group 4
- MYC amplification → Group 3

**Solutions:**
1. Check for chr 8 loss carefully
2. Consider tissue methylation for definitive classification
3. Proceed with treatment appropriate for either subgroup

---

#### Problem: Conflicting CNVs

```
Features detected:
✓ monosomy 6 (WNT)
✓ MYC amplification (Group 3)
```

**Interpretation:**
- Very rare
- May indicate:
  - Two tumor populations (intratumoral heterogeneity)
  - Artifact
  - Mixed sample

**Solutions:**
1. Visual inspection of CNV plots
2. Check tumor fraction and quality metrics
3. Repeat CSF sampling
4. Tissue analysis mandatory

---

## 📚 References

### Primary Publications

1. **Liu et al. 2021** - Serial assessment of measurable residual disease in medulloblastoma liquid biopsies  
   *Cancer Cell* 39(11):1519-1530  
   https://doi.org/10.1016/j.ccell.2021.08.012

2. **Markowitz et al. 2025** - Genome-wide cfDNA fragmentation patterns in cerebrospinal fluid reflect medulloblastoma groups  
   *NPJ Precision Oncology* 9:12  
   https://doi.org/10.1038/s41698-025-01067-5

3. **Escudero et al. 2020** - Circulating tumour DNA from the cerebrospinal fluid allows the characterisation and monitoring of medulloblastoma  
   *Nature Communications* 11:5376  
   https://doi.org/10.1038/s41467-020-19175-0

4. **Schwalbe et al. 2017** - Novel molecular subgroups for clinical classification and outcome prediction in childhood medulloblastoma  
   *Lancet Oncology* 18(7):958-971  
   https://doi.org/10.1016/S1470-2045(17)30434-6

5. **Christodoulou et al. 2023** - Combined low-pass whole genome and targeted sequencing in liquid biopsies for pediatric solid tumors  
   *NPJ Precision Oncology* 7:21  
   https://doi.org/10.1038/s41698-023-00361-0

---

### Technical Methods

6. **Adalsteinsson et al. 2017** - Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors (ichorCNA)  
   *Nature Communications* 8:1324  
   https://doi.org/10.1038/s41467-017-00965-y

7. **Li & Durbin 2009** - Fast and accurate short read alignment with Burrows-Wheeler transform (BWA)  
   *Bioinformatics* 25(14):1754-1760  
   https://doi.org/10.1093/bioinformatics/btp324

---

### Clinical Context

8. **Northcott et al. 2017** - The whole-genome landscape of medulloblastoma subtypes  
   *Nature* 547:311-317  
   https://doi.org/10.1038/nature22973

9. **Cavalli et al. 2017** - Intertumoral heterogeneity within medulloblastoma subgroups  
   *Cancer Cell* 31(6):737-754  
   https://doi.org/10.1016/j.ccell.2017.05.005

10. **Ramaswamy et al. 2016** - Risk stratification of childhood medulloblastoma in the molecular era  
    *Journal of Clinical Oncology* 34(5):405-414  
    https://doi.org/10.1200/JCO.2015.65.1031

---

## 📝 Citation

If you use this pipeline in your research, please cite:

**This pipeline:**
```
[Your Name] et al. (2026). Medulloblastoma CSF Liquid Biopsy Classification Pipeline.
GitHub: https://github.com/your-org/mb-classifier
```

**Key methods:**
```
Liu AP et al. (2021). Serial assessment of measurable residual disease 
in medulloblastoma liquid biopsies. Cancer Cell 39(11):1519-1530.

Adalsteinsson VA et al. (2017). Scalable whole-exome sequencing of 
cell-free DNA reveals high concordance with metastatic tumors. 
Nat Commun 8:1324.

Schwalbe EC et al. (2017). Novel molecular subgroups for clinical 
classification and outcome prediction in childhood medulloblastoma. 
Lancet Oncol 18(7):958-971.
```


---

### For Clinical Interpretation

**Consult with:**
- Pediatric neuro-oncologist
- Molecular pathologist
- Genetic counselor

**This is a research tool. Clinical decisions should incorporate:**
- Tissue diagnosis (gold standard)
- Imaging findings
- Clinical presentation
- Multidisciplinary team input

---

## ⚠️ Disclaimer

**IMPORTANT:**

This pipeline is provided for **research use only**. It has not been cleared or approved by regulatory agencies (FDA, CE-IVD, etc.) for clinical diagnostic use.


**Not a Replacement for:**
- Tissue histopathology
- Tissue molecular profiling (methylation arrays)
- Standard-of-care diagnostics
- Clinical judgment

**Use in clinical care requires:**
- Institutional validation
- Regulatory compliance
- Quality management system
- Appropriate clinical oversight

---

## 🎓 Acknowledgments

**Scientific Foundation:**
- St. Jude Children's Research Hospital (Liu et al.)
- Children's Hospital Los Angeles (Markowitz et al., Christodoulou et al.)
- Hospital Sant Joan de Déu, Barcelona (Escudero et al.)
- Newcastle University (Schwalbe et al.)

**Technical Tools:**
- Broad Institute (ichorCNA)
- Heng Li (BWA)
- SAMtools developers
- R/Bioconductor community


## 🔄 Version History

**Version 1.0 (April 2026):**
- Initial release
- CNV-based classification
- Risk stratification
- Clinical report generation

**Planned Features (v2.0):**
- Fragmentomics integration
- Machine learning classifier
- Multi-site validation results
- Docker containerization

---

## 🌐 Related Resources

**DKFZ Methylation Classifier (Gold Standard):**  
https://www.molecularneuropathology.org/mnp

**St. Jude Cloud (MB Data):**  
https://stjude.cloud/

**PeCan Data Portal:**  
https://pecan.stjude.cloud/

**Markowitz Fragmentomics Code:**  
https://github.com/amarkowitzchla/project-mb-lbseq-classification

---

**END OF README**
