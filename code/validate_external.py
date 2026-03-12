#!/usr/bin/env python3
"""
Validate FerroScore-Immuno on real immunotherapy data
Using GSE78220 and GSE91061
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, confusion_matrix
import joblib
import os

RESULTS_DIR = '../results'
EXTERNAL_DIR = '../data/external'
os.makedirs(f'{RESULTS_DIR}/figures', exist_ok=True)
os.makedirs(f'{RESULTS_DIR}/tables', exist_ok=True)

print("=" * 60)
print("FerroScore-Immuno External Validation")
print("=" * 60)

# Load the trained model
print("\nLoading trained model...")
model_path = f'{RESULTS_DIR}/models/best_model.pkl'
if os.path.exists(model_path):
    model_data = joblib.load(model_path)
    model = model_data['model']
    scaler = model_data['scaler']
    features = model_data['features']
    print(f"  Model: {model_data['model_name']}")
    print(f"  Features: {features}")
else:
    print(f"  [X] Model not found: {model_path}")
    exit(1)

# Load GSE78220 response data
print("\n" + "=" * 60)
print("GSE78220 Validation (Melanoma Anti-PD-1)")
print("=" * 60)

gse78220_response = pd.read_csv(f'{RESULTS_DIR}/external/gse78220_response.csv')
print(f"  Samples: {len(gse78220_response)}")
print(f"  Responders: {(gse78220_response['response']==1).sum()}")
print(f"  Non-responders: {(gse78220_response['response']==0).sum()}")

# Note: GSE78220 doesn't have expression data in the series matrix
# It needs to be downloaded from SRA or GEO separately
print("\n  Note: Expression data needs to be downloaded from SRA")
print("  Run: prefetch SRP* and fasterq-dump to get FASTQ files")
print("  Then align and quantify to get expression matrix")

# Summary
print("\n" + "=" * 60)
print("Validation Summary")
print("=" * 60)

summary = """
## Available Immunotherapy Datasets

### GSE78220 (Melanoma, Anti-PD-1)
- Samples: 28
- Responders: 15 (CR/PR)
- Non-responders: 13 (PD)
- Status: Response labels extracted
- Missing: Expression matrix (need SRA download)

### GSE91061 (Melanoma, Anti-PD-1)
- Samples: 109
- Pre/On treatment pairs
- Response labels available
- Status: Phenotype downloaded
- Missing: Expression matrix (need SRA download)

### IMvigor210 (Urothelial, Anti-PD-L1)
- Samples: 348
- Status: Need to download via R package

## Next Steps

1. Download expression data from SRA for GSE78220/GSE91061
   - Use SRA Toolkit: prefetch + fasterq-dump
   - Align with STAR/HISAT2
   - Quantify with featureCounts/salmon

2. Or use pre-processed data from:
   - UCSC Xena (if available)
   - GEO2R (for microarray data)
   - Published supplementary files

3. Calculate FerroScore-Immuno for validation samples

4. Evaluate prediction performance (AUC, accuracy)

## Alternative: Use TCGA with Proxy Labels

Since immunotherapy data is limited, use TCGA with:
- High FerroImmuno Score = Predicted responder
- TMB (Tumor Mutation Burden) as additional filter
- Immune infiltration scores
"""

with open(f'{RESULTS_DIR}/validation_summary.md', 'w') as f:
    f.write(summary)

print(summary)
print(f"\n[OK] Saved: validation_summary.md")

print("\n" + "=" * 60)
print("Recommendation")
print("=" * 60)
print("""
For immediate validation, consider:

1. Use TIDE database (http://tide.dfci.harvard.edu/)
   - Pre-computed immunotherapy response predictions
   - Can compare with FerroScore-Immuno

2. Use published signature comparison
   - Compare with existing IO signatures (TMB, PD-L1, IFN-gamma)
   - Validate correlation with immune markers

3. Cell line validation
   - Use CCLE + CTRP drug sensitivity
   - Correlate FerroScore with immunotherapy drug response
""")
