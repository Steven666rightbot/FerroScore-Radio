#!/usr/bin/env python3
"""
Download IMvigor210 dataset
Urothelial carcinoma anti-PD-L1 (Atezolizumab)
"""

import os
import requests

EXTERNAL_DIR = '../data/external'
os.makedirs(EXTERNAL_DIR, exist_ok=True)

print("=" * 60)
print("Downloading IMvigor210 Dataset")
print("=" * 60)

print("""
IMvigor210 is a Phase 2 study of Atezolizumab in urothelial carcinoma.
- 348 patients
- RNA-seq expression data
- Response labels (CR, PR, SD, PD)
- Survival data

Download options:
1. R package: IMvigor210CoreBiologies (recommended)
2. cBioPortal: http://www.cbioportal.org/study?id=blca_iatlas_imvigor210_2017
3. EGA: https://ega-archive.org/dacs/EGAC00001001611
""")

# Try to download from cBioPortal
print("\n1. Attempting download from cBioPortal...")

# cBioPortal data hub URL
# Note: cBioPortal requires clicking through UI, but has API
# For now, provide manual download instructions

print("""
   cBioPortal requires manual download:
   
   Steps:
   1. Visit: http://www.cbioportal.org/study?id=blca_iatlas_imvigor210_2017
   2. Click "Download" tab
   3. Download:
      - data_mrna_seq_v2_rsem.txt (expression)
      - data_clinical_sample.txt (sample info)
      - data_clinical_patient.txt (patient info)
   4. Save to: data/external/imvigor210/
""")

# Try R package installation script
print("\n2. Creating R script for IMvigor210CoreBiologies...")

r_script = """# IMvigor210 Download Script
# Run this in R

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install IMvigor210CoreBiologies
BiocManager::install("IMvigor210CoreBiologies")

# Load package
library(IMvigor210CoreBiologies)

# Load data
data(cds)

# Extract expression matrix (TPM)
expr_matrix = counts(cds, normalized=TRUE)
write.csv(expr_matrix, "IMvigor210_expression.csv")

# Extract phenotype data
pheno_data = pData(cds)
write.csv(pheno_data, "IMvigor210_phenotype.csv")

# Extract feature data (gene info)
feature_data = fData(cds)
write.csv(feature_data, "IMvigor210_genes.csv")

print("Download complete!")
print(paste("Samples:", ncol(expr_matrix)))
print(paste("Genes:", nrow(expr_matrix)))
"""

r_script_path = f'{EXTERNAL_DIR}/download_imvigor210.R'
with open(r_script_path, 'w') as f:
    f.write(r_script)

print(f"   [OK] Created: {r_script_path}")
print("\n   Run this script in R to download IMvigor210")

# Alternative: Try to find direct download links
print("\n3. Checking for direct download links...")

# Try to download from known academic mirrors
urls_to_try = [
    "https://raw.githubusercontent.com/cran/IMvigor210CoreBiologies/master/data/cds.rda",
]

for url in urls_to_try:
    print(f"   Trying: {url}")
    try:
        response = requests.head(url, timeout=10)
        print(f"   Status: {response.status_code}")
    except Exception as e:
        print(f"   Error: {e}")

print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print("""
To get IMvigor210 data:

Option 1 (Recommended): Use R script
   Run: Rscript download_imvigor210.R
   
Option 2: Manual download from cBioPortal
   Visit: http://www.cbioportal.org/study?id=blca_iatlas_imvigor210_2017
   Download expression and clinical data

Option 3: EGA (requires application)
   https://ega-archive.org/dacs/EGAC00001001611

Expected files:
- IMvigor210_expression.csv (RNA-seq TPM)
- IMvigor210_phenotype.csv (sample info + response)
- ~348 samples
""")
