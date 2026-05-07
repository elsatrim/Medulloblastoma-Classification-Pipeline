# 🚀 Quick Start Guide

## Installation (5 minutes + download time)

```bash
# 1. Clone repository
git clone https://github.com/YOUR_USERNAME/Medulloblastoma-Classification-Pipeline.git
cd Medulloblastoma-Classification-Pipeline

# 2. Run installer (installs conda environment + all tools)
bash install.sh

# 3. Activate environment
conda activate mb_classifier

# 4. Download reference genome (~30 minutes)
bash setup_reference.sh
```

## Running the Pipeline

### Paired-end data (typical):
```bash
bash run_pipeline.sh \
  SAMPLE_ID \
  reads_R1.fastq.gz \
  reads_R2.fastq.gz \
  output_directory \
  --threads 8
```

### Single-end data:
```bash
bash run_pipeline.sh \
  SAMPLE_ID \
  reads.fastq.gz \
  output_directory \
  --threads 8 \
  --single-end
```

### With clinical info:
```bash
bash run_pipeline.sh \
  MB001 \
  R1.fq.gz R2.fq.gz \
  results/MB001 \
  --threads 8 \
  --metastasis no \
  --residual no \
  --histology classic
```

## Output

```
output_directory/
├── 01_preprocessing/
│   └── SAMPLE.sorted.bam
├── 02_ichorCNA/
│   ├── SAMPLE.seg.txt
│   ├── SAMPLE_genomeWide.pdf
│   └── SAMPLE.params.txt
├── 04_classification/
│   ├── SAMPLE_classification.csv
│   └── SAMPLE_classification.json
├── 05_risk_stratification/
│   └── SAMPLE_risk.csv
└── 06_report/
    └── SAMPLE_MB_Report.pdf
```

## Troubleshooting

**"conda: command not found"**
→ Install Miniconda: https://docs.conda.io/en/latest/miniconda.html

**"Reference genome not found"**
→ Run: `bash setup_reference.sh`

**"ichorCNA not found"**
→ Re-run: `bash install.sh`

**Need help?**
→ Open an issue on GitHub
