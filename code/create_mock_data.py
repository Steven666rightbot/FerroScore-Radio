#!/usr/bin/env python3
"""
Create mock expression data for testing FerroScore-Immuno
"""

import pandas as pd
import numpy as np
import os

PROC_DIR = '../data/processed'
os.makedirs(PROC_DIR, exist_ok=True)

# Load gene sets
gene_file = '../gene_sets/ferro_immuno_genes.txt'
gene_sets = {}
current_set = None

with open(gene_file, 'r', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if line.startswith('##'):
            current_set = line.split('(', 1)[0].strip()
            current_set = current_set.replace('## ', '').replace('##', '').strip()
            if '. ' in current_set:
                current_set = current_set.split('. ', 1)[1]
            gene_sets[current_set] = []
        elif line and not line.startswith('#') and current_set:
            gene_sets[current_set].append(line)

# Get all unique genes
all_genes = []
for genes in gene_sets.values():
    all_genes.extend(genes)
all_genes = list(set(all_genes))

print(f"Total unique genes: {len(all_genes)}")
print(f"Gene sets: {list(gene_sets.keys())}")

# Create mock expression matrix (genes x samples)
n_samples = 200  # Mock 200 samples
np.random.seed(42)

# Generate realistic expression values (log2 TPM-like)
expr_matrix = pd.DataFrame(
    np.random.lognormal(mean=2, sigma=1.5, size=(len(all_genes), n_samples)),
    index=all_genes,
    columns=[f'SAMPLE_{i:03d}' for i in range(n_samples)]
)

print(f"\nMock expression matrix: {expr_matrix.shape}")
print(f"Expression range: [{expr_matrix.min().min():.2f}, {expr_matrix.max().max():.2f}]")

# Save
output_file = f'{PROC_DIR}/tcga_ferro_immuno_expression.csv'
expr_matrix.to_csv(output_file)
print(f"\nSaved: {output_file}")

# Also create mock clinical data
clinical = pd.DataFrame({
    'sample_id': expr_matrix.columns,
    'cancer_type': np.random.choice(['BRCA', 'LUAD', 'LUSC', 'SKCM', 'BLCA'], n_samples),
    'age': np.random.randint(30, 80, n_samples),
    'gender': np.random.choice(['male', 'female'], n_samples),
    'received_immunotherapy': np.random.choice([True, False], n_samples, p=[0.3, 0.7])
})
clinical.to_csv(f'{PROC_DIR}/tcga_clinical.csv', index=False)
print(f"Saved: tcga_clinical.csv")

# Create mock survival data
survival = pd.DataFrame({
    'sample': expr_matrix.columns,
    'OS': np.random.choice([0, 1], n_samples, p=[0.4, 0.6]),
    'OS.time': np.random.exponential(1000, n_samples).astype(int)
})
survival.to_csv(f'{PROC_DIR}/tcga_survival.csv', index=False)
print(f"Saved: tcga_survival.csv")

print("\nMock data creation complete!")
