# Medulloblastoma-Classification-Pipeline

📚 SCIENTIFIC FOUNDATION:
This pipeline directly implements methods from:

Liu et al. 2021 (Cancer Cell) - St. Jude's MB liquid biopsy study

https://www.cell.com/cancer-cell/fulltext/S1535-6108(21)00466-3
lcWGS for MRD detection, 3+ month lead time


Markowitz et al. 2025 (NPJ Precision Oncology) - CHLA fragmentomics

https://www.nature.com/articles/s41698-025-01067-5
Code available: https://github.com/amarkowitzchla/project-mb-lbseq-classification
Machine learning classification (AUC 0.94)


Escudero et al. 2020 (Nature Communications) - Subgroup classification

https://www.nature.com/articles/s41467-020-19175-0
WES of CSF ctDNA for MB subtyping


Christodoulou et al. 2023 (NPJ Precision Oncology) - LBSeq4Kids

https://www.nature.com/articles/s41698-023-00361-0
Clinical validation at CHLA




🎯 PIPELINE OVERVIEW:
INPUT: CSF cfDNA from MB patient (collected during EVD placement)
  ↓
STEP 1: Low-coverage WGS (0.5-1× coverage)
  ↓
STEP 2: CNV detection (ichorCNA)
  ↓
STEP 3: MB subgroup classification (WNT, SHH, Group 3, Group 4)
  ↓
STEP 4: Risk stratification
  ↓
OUTPUT: Clinical report with subgroup + prognosis
  
TIMELINE: 3-5 days (CSF → Report)

📁 COMPLETE PROJECT STRUCTURE:
bashmb-classifier/
├── install_dependencies.sh
├── config/
│   ├── mb_cnv_signatures.yaml
│   ├── reference_paths.yaml
│   └── clinical_thresholds.yaml
├── scripts/
│   ├── 01_preprocessing.sh
│   ├── 02_run_ichorCNA.R
│   ├── 03_classify_mb_subgroup.R
│   ├── 04_fragmentomics_features.R    # Markowitz et al.
│   ├── 05_risk_stratification.R
│   └── 06_generate_report.Rmd
├── reference_data/
│   ├── download_references.sh
│   └── mb_training_data/
├── data/
│   └── test_samples/
├── results/
└── run_mb_pipeline.sh
