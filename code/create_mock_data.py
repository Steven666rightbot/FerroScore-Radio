#!/usr/bin/env python3
"""
Create Mock Data for FerroScore-Radio Testing
WARNING: This is synthetic data for code testing only!
DO NOT use for actual research or publication.
"""

import os
import pandas as pd
import numpy as np

# Set random seed for reproducibility
np.random.seed(42)

# Paths
PROC_DIR = '../data/processed'
os.makedirs(PROC_DIR, exist_ok=True)

print("=" * 70)
print("WARNING: Generating MOCK DATA for testing only!")
print("This data is synthetic and should NOT be used for research.")
print("=" * 70)

# Load gene sets
gene_sets = {
    'Ferroptosis Driver': ['ACSL4', 'LPCAT3', 'ALOX15', 'ALOX5', 'NOX1', 'NOX4', 'P53', 'SAT1', 'CARS1'],
    'Ferroptosis Suppressor': ['GPX4', 'SLC7A11', 'SLC3A2', 'NFE2L2', 'HMOX1', 'FTH1', 'FTL', 'SOD1', 'SOD2', 'GCLC', 'GCLM'],
    'DNA Repair': ['BRCA1', 'BRCA2', 'ATM', 'ATR', 'CHEK1', 'CHEK2', 'TP53', 'RAD51', 'PARP1'],
    'ROS-related': ['NOX1', 'NOX2', 'NOX4', 'SOD1', 'SOD2', 'CAT', 'GPX1', 'PRDX1']
}

all_genes = []
for genes in gene_sets.values():
    all_genes.extend(genes)
all_genes = list(set(all_genes))

print(f"\nGene set: {len(all_genes)} unique genes")

# Create mock expression data
print("\nGenerating mock expression data...")

n_samples = 500  # 500 mock samples
n_genes = len(all_genes)

# Generate realistic expression values (log2 TPM)
# Most genes have low expression, some have high
expression_matrix = np.random.lognormal(mean=2, sigma=2, size=(n_genes, n_samples))

# Create sample IDs
sample_ids = [f'TCGA-TEST-{i:04d}' for i in range(n_samples)]

# Create DataFrame
expr_df = pd.DataFrame(
    expression_matrix,
    index=all_genes,
    columns=sample_ids
)

print(f"  Expression matrix: {expr_df.shape}")

# Save expression data
expr_df.to_csv(f'{PROC_DIR}/tcga_ferro_radio_expression.csv')
print(f"  ✓ Saved: tcga_ferro_radio_expression.csv")

# Create mock clinical data
print("\nGenerating mock clinical data...")

cancer_types = ['BRCA', 'LUAD', 'LUSC', 'HNSC', 'GBM', 'LGG', 'CESC', 'ESCA', 'BLCA']
stages = ['Stage I', 'Stage II', 'Stage III', 'Stage IV']

clinical_data = []
for i, sample_id in enumerate(sample_ids):
    received_rt = np.random.choice([True, False], p=[0.4, 0.6])
    
    clinical_data.append({
        'sample_id': sample_id,
        'cancer_type': np.random.choice(cancer_types),
        'age_at_initial_pathologic_diagnosis': np.random.randint(30, 85),
        'gender': np.random.choice(['male', 'female']),
        'stage': np.random.choice(stages),
        'radiotherapy': 'yes' if received_rt else 'no',
        'received_radiotherapy': received_rt,
        'targeted_therapy': np.random.choice(['yes', 'no'], p=[0.2, 0.8])
    })

clinical_df = pd.DataFrame(clinical_data)
print(f"  Clinical data: {clinical_df.shape}")
print(f"    - With RT: {clinical_df['received_radiotherapy'].sum()}")
print(f"    - Without RT: {(~clinical_df['received_radiotherapy']).sum()}")

# Save clinical data
clinical_df.to_csv(f'{PROC_DIR}/tcga_clinical.csv', index=False)
print(f"  ✓ Saved: tcga_clinical.csv")

# Create mock survival data
print("\nGenerating mock survival data...")

survival_data = []
for sample_id in sample_ids:
    # Generate survival time (days)
    survival_time = np.random.exponential(scale=1000) + 100
    
    # Generate event status (1 = death, 0 = censored)
    event = np.random.choice([0, 1], p=[0.6, 0.4])
    
    survival_data.append({
        'sample': sample_id,
        'OS.time': survival_time,
        'OS': event,
        'PFI.time': survival_time * 0.8,
        'PFI': event
    })

survival_df = pd.DataFrame(survival_data)
print(f"  Survival data: {survival_df.shape}")
print(f"    - Events (deaths): {survival_df['OS'].sum()}")
print(f"    - Censored: {(~survival_df['OS'].astype(bool)).sum()}")

# Save survival data
survival_df.to_csv(f'{PROC_DIR}/tcga_survival.csv', index=False)
print(f"  ✓ Saved: tcga_survival.csv")

# Save gene list
gene_list_df = pd.DataFrame({'gene': all_genes})
gene_list_df.to_csv(f'{PROC_DIR}/ferro_radio_gene_list.csv', index=False)
print(f"  ✓ Saved: ferro_radio_gene_list.csv ({len(all_genes)} genes)")

print("\n" + "=" * 70)
print("MOCK DATA GENERATION COMPLETE!")
print("=" * 70)
print("\n⚠️  IMPORTANT REMINDERS:")
print("   1. This data is SYNTHETIC (fake)")
print("   2. For CODE TESTING only")
print("   3. DO NOT use for research or publication")
print("   4. Replace with real TCGA data when available")
print("\nNext step: Run 03_ferroscore_algorithm.py")
print("=" * 70)
