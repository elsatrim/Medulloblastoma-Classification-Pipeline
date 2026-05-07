# Medulloblastoma Classification Pipeline

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.9+-green.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)

Complete bioinformatics pipeline for medulloblastoma subgroup classification from cerebrospinal fluid (CSF) cell-free DNA using low-coverage whole genome sequencing.

Based on:
- **Liu et al. 2021** - Cancer Cell (CNV-based classification)
- **Escudero et al. 2020** - Acta Neuropathologica (CNV signatures)
- **Markowitz et al. 2025** - npj Precision Oncology (fragmentomics, optional)

## 📋 Overview

This pipeline classifies medulloblastoma into four molecular subgroups (WNT, SHH, Group 3, Group 4) using copy number variation (CNV) patterns from low-coverage whole genome sequencing (lcWGS) of CSF cfDNA.

### Features

- ✅ **CNV-based classification** (85% accuracy, Liu 2021)
- ✅ **Low-coverage sequencing** (0.5-1× coverage, cost-effective)
- ✅ **Tumor fraction estimation** (ichorCNA integration)
- ✅ **Risk stratification** (integrates clinical factors)
- ✅ **Fragmentomics-ready** (optional enhancement, 94% accuracy)
- ✅ **Automated reporting** (PDF generation with RMarkdown)

## 🚀 Quick Start

### Prerequisites

- macOS or Linux
- Conda/Miniconda
- 16 GB RAM minimum
- 4-8 CPU cores

### Installation

```bash
# Clone repository
git clone https://github.com/yourusername/Medulloblastoma-Classification-Pipeline.git
cd Medulloblastoma-Classification-Pipeline

# Install dependencies
bash install_dependencies.sh

# Activate environment
conda activate mb_classifier

# Download reference genome (GRCh38)
cd ~/mb_reference
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

### Usage

```bash
# Run complete pipeline
bash run_mb_pipeline.sh \
  SAMPLE_ID \
  reads_R1.fastq.gz \
  reads_R2.fastq.gz \
  output_directory \
  --threads 8 \
  --metastasis no \
  --residual no \
  --histology classic
```

## 📊 Pipeline Workflow

```
INPUT: FASTQ files (paired-end)
  ↓
1. PREPROCESSING
   • Quality control (FastQC)
   • Alignment (BWA-MEM)
   • Duplicate marking
  ↓
2. CNV DETECTION
   • ichorCNA analysis
   • Tumor fraction estimation
   • Copy number calling
  ↓
3. FRAGMENTOMICS (Optional)
   • Fragment length distribution
   • Short-to-long ratios
   • End motif analysis
  ↓
4. CLASSIFICATION
   • CNV pattern matching
   • Liu 2021 decision tree
   • Confidence scoring
  ↓
5. RISK STRATIFICATION
   • Molecular + clinical integration
   • Risk group assignment
  ↓
6. REPORT GENERATION
   • PDF report (RMarkdown)
   • Plots and visualizations
   • Literature validation
  ↓
OUTPUT: Classification + Risk + Report
```

## 📂 Directory Structure

```
Medulloblastoma-Classification-Pipeline/
├── LICENSE
├── README.md
├── install_dependencies.sh       # Environment setup
├── run_mb_pipeline.sh           # Main pipeline script
├── config/
│   └── mb_cnv_signatures.yaml   # CNV signatures definition
├── scripts/
│   ├── 01_preprocessing.sh      # Alignment & QC
│   ├── 02_run_ichorCNA.R        # CNV detection
│   ├── 03_classify_mb_subgroup.R    # Classification logic
│   ├── 04_fragmentomics_analysis.py # Fragmentomics (optional)
│   ├── 05_risk_stratification.R     # Risk assignment
│   └── 06_generate_report.Rmd       # PDF report
└── test_data/                   # Example data (not included)
```

## 🧬 Classification Logic

Based on hierarchical CNV signatures (Escudero et al. 2020):

```
IF monosomy_6:
    → WNT (100% specific)

ELSE IF myc_amplification:
    → GROUP_3 (17% frequency)

ELSE IF chr9q_loss:
    → SHH (32-100% frequency)

ELSE IF isochromosome_17q AND chromosome_8_loss:
    → GROUP_4 (67% + 20% frequencies)

ELSE:
    → INDETERMINATE
```

## 📊 Validation Results

Test sample: ERR550408_MB_POOL (pooled CSF from 5 patients)

| Metric | Value |
|--------|-------|
| **Classification** | GROUP_4 |
| **Confidence** | HIGH |
| **Tumor Fraction** | 38.33% |
| **Ploidy** | 2.02 |
| **Isochromosome 17q** | ✓ Detected |
| **Chromosome 8 loss** | ✓ Detected |
| **Literature concordance** | ✓ Validated |

## 🔬 Technical Specifications

### System Requirements

- **CPU:** 4-8 cores
- **RAM:** 16 GB minimum, 32 GB recommended
- **Storage:** 50 GB for reference + results
- **Runtime:** ~90-150 minutes per sample

### Software Dependencies

- **Alignment:** BWA-MEM 0.7.17
- **CNV calling:** ichorCNA 0.3.2
- **Python:** 3.9+ (pysam, matplotlib, seaborn)
- **R:** 4.0+ (data.table, ggplot2, rmarkdown)

### Input Requirements

- **Sequencing:** Paired-end lcWGS
- **Coverage:** 0.5-1× (for CNV)
- **Coverage:** >30× (for mutations, not used)
- **Sample type:** CSF cell-free DNA
- **Format:** FASTQ (gzipped)

## 📚 References

1. **Liu et al. 2021** - Cancer Cell  
   "Liquid biopsy using cell-free DNA enables medulloblastoma classification"  
   DOI: 10.1016/j.ccell.2021.04.006

2. **Escudero et al. 2020** - Acta Neuropathologica  
   "Molecular subgroup-specific copy number alterations"  
   DOI: 10.1007/s00401-020-02146-w

3. **Markowitz et al. 2025** - npj Precision Oncology  
   "Genome-wide cfDNA fragmentation patterns"  
   DOI: 10.1038/s41698-025-01067-5

4. **Schwalbe et al. 2017** - Lancet Oncology  
   "Molecular subgroups for clinical classification"  
   DOI: 10.1016/S1470-2045(17)30243-7

## ⚠️ Limitations

1. **Pooled data validation:** Test sample is pooled (5 patients), accuracy cannot be calculated
2. **Risk stratification:** Molecular classification only; complete risk requires histology + imaging
3. **Fragmentomics:** Requires paired-end sequencing with >90bp fragments
4. **Clinical use:** Research tool only; not validated for clinical diagnostics

## 🚀 Future Work

### Phase 1: Clinical Validation
- Individual patient samples (n=20-30)
- Sensitivity/specificity calculation
- Comparison with tissue diagnosis

### Phase 2: Fragmentomics Integration
- Paired-end sequencing dataset
- Markowitz 2025 methodology validation
- Target: 94% accuracy

### Phase 3: Clinical Implementation
- Hospital EHR integration
- Automated reporting pipeline
- Cost-effectiveness analysis

## 🤝 Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## 👤 Author

**Elsa Sanchez Fernandez**  
Master's Thesis Project  
University of Geneva, 2024

## 🙏 Acknowledgments

- ichorCNA developers (Broad Institute)
- Liu lab (Harvard/Stanford) for methodology
- Markowitz lab (UCLA) for fragmentomics approach

## 📧 Contact

For questions or issues, please open a GitHub issue or contact [your email].

---

**Note:** This is a research tool. Clinical decisions should not be made based solely on this pipeline's output without proper validation and integration with clinical data.
