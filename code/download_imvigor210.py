#!/usr/bin/env python3
"""
Download and process IMvigor210 dataset
Urothelial carcinoma anti-PD-L1 (Atezolizumab)
"""

import os
import pandas as pd
import requests
from tqdm import tqdm

EXTERNAL_DIR = '../data/external'
os.makedirs(EXTERNAL_DIR, exist_ok=True)

print("=" * 60)
print("IMvigor210 Dataset Download")
print("=" * 60)

print("""
IMvigor210 is a Phase 2 study of Atezolizumab in urothelial carcinoma.
Dataset includes:
- 348 patients
- RNA-seq expression data
- Response labels (CR, PR, SD, PD)
- Survival data

Download options:
1. R package: IMvigor210CoreBiologies
   install.packages('BiocManager')
   BiocManager::install('IMvigor210CoreBiologies')

2. Web portal: https://research-pub.gene.com/IMvigor210CoreBiologies/

3. Direct download (if available):
""")

# Try to download from known URLs
urls_to_try = [
    "https://research-pub.gene.com/IMvigor210CoreBiologies/data/IMvigor210.csv",
]

print("Attempting download...")
for url in urls_to_try:
    print(f"  Trying: {url}")
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            output_file = f'{EXTERNAL_DIR}/IMvigor210.csv'
            with open(output_file, 'wb') as f:
                f.write(response.content)
            print(f"  [OK] Downloaded: {output_file}")
            break
        else:
            print(f"  Status: {response.status_code}")
    except Exception as e:
        print(f"  Error: {e}")

# Create instructions file
instructions = """
# IMvigor210 Dataset Download Instructions

## Option 1: R Package (Recommended)

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install IMvigor210CoreBiologies
BiocManager::install("IMvigor210CoreBiologies")

# Load and extract data
library(IMvigor210CoreBiologies)
data(cds)

# Extract expression matrix
expr_matrix = counts(cds)
write.csv(expr_matrix, "IMvigor210_expression.csv")

# Extract phenotype data
pheno_data = pData(cds)
write.csv(pheno_data, "IMvigor210_phenotype.csv")
```

## Option 2: Manual Download

1. Visit: https://research-pub.gene.com/IMvigor210CoreBiologies/
2. Download the dataset
3. Place files in: data/external/

## Dataset Information

- **Study**: IMvigor210 Phase 2
- **Cancer**: Urothelial carcinoma (BLCA)
- **Treatment**: Atezolizumab (anti-PD-L1)
- **Samples**: 348 patients
- **Data types**: RNA-seq, clinical, response, survival

## Response Definitions

- CR: Complete Response
- PR: Partial Response  
- SD: Stable Disease
- PD: Progressive Disease

## Use in FerroScore-Immuno

This dataset can be used as:
1. External validation cohort
2. Training data for immunotherapy-specific model
3. Benchmark for performance comparison
"""

with open(f'{EXTERNAL_DIR}/IMvigor210_instructions.md', 'w') as f:
    f.write(instructions)

print(f"\n[OK] Saved: IMvigor210_instructions.md")
print("\nPlease follow the instructions to download the dataset.")
