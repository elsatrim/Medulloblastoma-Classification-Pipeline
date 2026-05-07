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

## 🧬 Fragmentomics Analysis (Optional)

### What is Fragmentomics?

Fragmentomics analyzes cell-free DNA fragment patterns - the **size distribution** and **end motifs** of DNA fragments released from tumor cells. Different cancer types and molecular subtypes produce characteristic fragmentation patterns.

Based on **Markowitz et al. 2025**, fragmentomics can improve medulloblastoma classification accuracy from 85% (CNV-only) to **94%** (CNV + fragmentomics).

### Biological Basis

Cell-free DNA fragments are produced by:
1. **Apoptosis** (programmed cell death) - produces ~167 bp fragments (nucleosome-sized)
2. **Necrosis** (uncontrolled death) - produces irregular fragments
3. **Active secretion** - produces shorter fragments

Tumor cells have:
- Different chromatin accessibility patterns
- Distinct nucleosome positioning
- Altered DNA packaging
- **Result:** Unique fragmentation "fingerprints" for each MB subgroup

### Key Fragmentomics Features

#### 1. **Fragment Length Distribution**

**Short fragments (90-150 bp):**
- Represent highly accessible chromatin regions
- Enriched in active genes
- Tumor-specific patterns

**Long fragments (151-220 bp):**
- Represent nucleosome-protected regions  
- More stable chromatin structure
- Background/normal tissue contribution

**Short-to-Long (S/L) Ratio:**
- **Group 4:** S/L = 0.8-1.2 (balanced)
- **Group 3:** S/L = 0.9-1.3 (slightly elevated)
- **WNT:** S/L = 0.7-1.0 (more long fragments)
- **SHH:** S/L = 1.0-1.4 (more short fragments)

#### 2. **Fragment End Motifs**

The 4-nucleotide sequences at fragment ends reveal:
- Nuclease cleavage preferences
- Chromatin accessibility
- Transcription factor binding

**Characteristic Patterns:**
- **Group 4:** CCCA, CCAG enrichment
- **Group 3:** Different motif preference (TGCA, AGCC)
- **WNT/SHH:** Distinct nucleosome positioning signatures

### MB Subgroup-Specific Patterns

Based on Markowitz et al. 2025 analysis of 73 CSF cfDNA samples:

| Feature | WNT | SHH | Group 3 | Group 4 |
|---------|-----|-----|---------|---------|
| **S/L Ratio** | 0.75 ± 0.12 | 1.15 ± 0.18 | 1.05 ± 0.15 | 0.95 ± 0.10 |
| **Mean Length** | 172 bp | 158 bp | 165 bp | 168 bp |
| **Top Motif** | CCAG | AGCC | TGCA | CCCA |
| **Fragment Diversity** | Low | High | Medium | Medium |

### How to Use Fragmentomics

#### Requirements

**Data requirements:**
- ✅ **Paired-end sequencing** (required)
- ✅ Fragment length >90 bp (required)
- ✅ Coverage: 0.5-1× sufficient
- ❌ Single-end sequencing (incompatible)

#### Running with Fragmentomics

```bash
# Standard run (fragmentomics enabled by default)
bash run_pipeline.sh \
  SAMPLE_ID \
  reads_R1.fastq.gz \
  reads_R2.fastq.gz \
  output/

# Skip fragmentomics (CNV-only)
bash run_pipeline.sh \
  SAMPLE_ID \
  reads_R1.fastq.gz \
  reads_R2.fastq.gz \
  output/ \
  --no-fragmentomics
```

#### Output Files

```
output/03_fragmentomics/
├── SAMPLE_fragmentomics.png          # 4-panel visualization
├── SAMPLE_fragmentomics_metrics.json # Quantitative metrics
└── SAMPLE_end_motifs.csv            # Top end motifs (if reference provided)
```

#### Fragmentomics Plot (4 panels)

1. **Fragment Length Distribution**
   - Histogram of fragment sizes (50-600 bp)
   - Shows short vs long fragment balance

2. **Cumulative Distribution**
   - Cumulative percentage by length
   - Helps identify size cutoffs

3. **Short/Long Ratio by Chromosome**
   - S/L ratio for each chromosome
   - Identifies chromosomal patterns

4. **Top Fragment End Motifs**
   - Most frequent 4bp end sequences
   - Shows tumor-specific signatures

### Integration with CNV Classification

The pipeline uses fragmentomics to **enhance confidence** when patterns agree:

```
CNV Classification: GROUP_4 (based on iso17q + chr8 loss)
Fragmentomics S/L: 0.95 (within Group 4 range: 0.8-1.2)
→ Confidence upgraded: HIGH → VERY_HIGH
```

**Concordance states:**
- **CONCORDANT:** Fragmentomics pattern matches CNV classification
- **DISCORDANT:** Patterns disagree (flag for review)
- **N/A:** Fragmentomics not run or inconclusive

### Clinical Interpretation

**When fragmentomics is concordant:**
- Higher confidence in classification
- Can guide treatment decisions with more certainty
- Reduces risk of misclassification

**When fragmentomics is discordant:**
- May indicate mixed tumors
- Sample quality issues
- Consider additional testing (tissue biopsy)

### Limitations

1. **Data requirements:** Needs paired-end sequencing with adequate fragment lengths
2. **Reference data:** Limited to Markowitz 2025 dataset (n=73)
3. **Computational cost:** Adds ~30-60 minutes to pipeline runtime
4. **Not validated:** For clinical diagnostics (research use only)

### Technical Notes

**Fragment extraction:**
- Filters for properly paired reads
- Fragment length range: 90-600 bp
- Excludes PCR duplicates and low-quality reads

**Metrics calculated:**
- Mean fragment length
- Short-to-long ratio (90-150bp / 151-220bp)
- Fragment length distribution
- End motif frequencies (if reference provided)

**End motif extraction** (optional):
- Requires reference genome
- Extracts 4bp sequences at fragment ends
- Compares to published MB subgroup signatures

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
| **Fragmentomics** | Not applicable* |

*Test dataset consists of single-end 36bp reads, which are incompatible with fragmentomics analysis (requires paired-end reads with >90bp fragments). CNV-based classification alone successfully identified the correct molecular subgroup.

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


## 🤝 Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## 📄 License

This project is licensed under the MIT License.

## 👤 Author

**Elsa Sanchez Fernandez**  

## 🙏 Acknowledgments

- ichorCNA developers (Broad Institute)
- Liu lab (Harvard/Stanford) for methodology
- Markowitz lab (UCLA) for fragmentomics approach

## 📧 Contact

For questions or issues, please open a GitHub issue or contact the author.

---

**Note:** This is a research tool. Clinical decisions should not be made based solely on this pipeline's output without proper validation and integration with clinical data.
