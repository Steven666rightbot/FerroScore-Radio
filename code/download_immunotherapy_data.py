#!/usr/bin/env python3
"""
Download real immunotherapy datasets for FerroScore-Immuno validation
"""

import os
import pandas as pd
import requests
from tqdm import tqdm

EXTERNAL_DIR = '../data/external'
os.makedirs(EXTERNAL_DIR, exist_ok=True)

print("=" * 60)
print("Downloading Immunotherapy Datasets")
print("=" * 60)

# Dataset 1: IMvigor210 (Urothelial carcinoma, anti-PD-L1)
# Available from R package: IMvigor210CoreBiologies
print("\n1. IMvigor210 Dataset")
print("   Source: Urothelial carcinoma anti-PD-L1 trial")
print("   Note: Download from https://research-pub.gene.com/IMvigor210CoreBiologies/")
print("   Or use R: BiocManager::install('IMvigor210CoreBiologies')")

# Dataset 2: GSE78220 (Melanoma, anti-PD-1) - already have
print("\n2. GSE78220 Dataset")
print("   Status: Already downloaded")
print("   Location: data/external/GSE78220_series_matrix.txt.gz")

# Dataset 3: GSE91061 (Melanoma, anti-PD-1) - already have phenotype
print("\n3. GSE91061 Dataset")
print("   Status: Phenotype downloaded")
print("   Location: data/external/GSE91061_phenotype.csv")
print("   Note: Expression data needs SRA download")

# Dataset 4: TIDE database (Tumor Immune Dysfunction and Exclusion)
print("\n4. TIDE Database")
print("   URL: http://tide.dfci.harvard.edu/")
print("   Contains multiple immunotherapy cohorts")

# Dataset 5: CheckMate trials (Nivolumab)
print("\n5. CheckMate Trials")
print("   CheckMate 025, 275, etc.")
print("   Available via dbGaP or publications")

# Create a summary file
summary = """
# Immunotherapy Datasets for FerroScore-Immuno

## Available Datasets

### 1. GSE78220 - Melanoma Anti-PD-1
- **Cancer**: Melanoma
- **Treatment**: Anti-PD-1 (Pembrolizumab)
- **Samples**: ~28
- **Status**: Downloaded
- **Use**: Validation cohort

### 2. GSE91061 - Melanoma Anti-PD-1
- **Cancer**: Melanoma  
- **Treatment**: Anti-PD-1
- **Samples**: ~109
- **Status**: Phenotype only
- **Action**: Need to download expression from SRA

### 3. IMvigor210 - Urothelial Carcinoma Anti-PD-L1
- **Cancer**: Urothelial carcinoma (BLCA)
- **Treatment**: Atezolizumab (anti-PD-L1)
- **Samples**: ~348
- **Status**: Need to download
- **Source**: Gene.com research portal

### 4. TIDE Database
- **Cancers**: Multiple
- **Treatments**: Anti-PD-1, Anti-CTLA-4
- **Status**: Web resource
- **URL**: http://tide.dfci.harvard.edu/

## Next Steps

1. Process GSE78220 (already have data)
2. Download IMvigor210 dataset
3. Use TCGA as training set (with proxy labels)
4. Validate on immunotherapy cohorts

## TCGA Proxy Label Strategy

Since TCGA doesn't have immunotherapy data, use:
- High FerroImmuno Score + High immune infiltration = Likely responder
- TMB (Tumor Mutation Burden) as proxy for immunogenicity
- CIBERSORT immune cell fractions
"""

with open(f'{EXTERNAL_DIR}/immunotherapy_datasets_summary.md', 'w') as f:
    f.write(summary)

print("\n" + "=" * 60)
print("Summary saved to: data/external/immunotherapy_datasets_summary.md")
print("=" * 60)

# Process GSE78220 if available
print("\nProcessing GSE78220...")
gse78220_file = f'{EXTERNAL_DIR}/GSE78220_series_matrix.txt.gz'

if os.path.exists(gse78220_file):
    print(f"  Found: {gse78220_file}")
    print("  Run: python analyze_gse7820.py to process")
else:
    print(f"  Not found: {gse78220_file}")

print("\nRecommended next action:")
print("  1. Process GSE78220 data")
print("  2. Download IMvigor210 dataset")
print("  3. Re-train model with TCGA + validate on immunotherapy cohorts")
