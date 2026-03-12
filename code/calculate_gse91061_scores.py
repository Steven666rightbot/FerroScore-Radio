#!/usr/bin/env python3
"""
Calculate FerroScore-Immuno for GSE91061
"""

import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import MinMaxScaler

RESULTS_DIR = '../results/external'
os.makedirs(RESULTS_DIR, exist_ok=True)

print("=" * 60)
print("Calculating FerroScore-Immuno for GSE91061")
print("=" * 60)

# Load mapped FPKM data
print("\n1. Loading GSE91061 mapped data...")
fpkm_file = f'{RESULTS_DIR}/gse91061_fpkm_mapped.csv'
expr_matrix = pd.read_csv(fpkm_file, index_col=0)
print(f"   Shape: {expr_matrix.shape}")

# Load gene sets
print("\n2. Loading gene sets...")
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

for name, genes in gene_sets.items():
    print(f"   {name}: {len(genes)} genes")

# Get all unique genes
all_genes = []
for genes in gene_sets.values():
    all_genes.extend(genes)
all_genes = list(set(all_genes))
print(f"\n   Total unique genes: {len(all_genes)}")

# Check available genes
available_genes = [g for g in all_genes if g in expr_matrix.index]
missing_genes = [g for g in all_genes if g not in expr_matrix.index]
print(f"   Available: {len(available_genes)}/{len(all_genes)}")
print(f"   Missing: {len(missing_genes)}")
if missing_genes:
    print(f"   Missing genes: {missing_genes}")

# Calculate scores
def calculate_ferroscore(expr_matrix, driver_genes, suppressor_genes):
    print("\n3. Calculating FerroScore...")
    
    driver_in_data = [g for g in driver_genes if g in expr_matrix.index]
    suppressor_in_data = [g for g in suppressor_genes if g in expr_matrix.index]
    
    print(f"   Driver genes: {len(driver_in_data)}/{len(driver_genes)}")
    print(f"   Suppressor genes: {len(suppressor_in_data)}/{len(suppressor_genes)}")
    
    # Log transform FPKM
    expr_log = np.log2(expr_matrix + 1)
    
    # Calculate scores
    driver_scores = expr_log.loc[driver_in_data].mean() if driver_in_data else pd.Series(0, index=expr_matrix.columns)
    suppressor_scores = expr_log.loc[suppressor_in_data].mean() if suppressor_in_data else pd.Series(0, index=expr_matrix.columns)
    
    # Combined score
    ferroscore = driver_scores - suppressor_scores
    
    # Normalize
    scaler = MinMaxScaler()
    ferroscore_scaled = scaler.fit_transform(ferroscore.values.reshape(-1, 1)).flatten()
    
    return pd.Series(ferroscore_scaled, index=expr_matrix.columns, name='FerroScore')

def calculate_immune_score(expr_matrix, immune_genes):
    print("\n4. Calculating Immune Score...")
    
    immune_in_data = [g for g in immune_genes if g in expr_matrix.index]
    print(f"   Immune genes: {len(immune_in_data)}/{len(immune_genes)}")
    
    expr_log = np.log2(expr_matrix + 1)
    immune_score = expr_log.loc[immune_in_data].mean() if immune_in_data else pd.Series(0, index=expr_matrix.columns)
    
    scaler = MinMaxScaler()
    immune_scaled = scaler.fit_transform(immune_score.values.reshape(-1, 1)).flatten()
    
    return pd.Series(immune_scaled, index=expr_matrix.columns, name='Immune_Score')

# Get gene lists
driver_genes = gene_sets.get('Ferroptosis Driver Genes', [])
suppressor_genes = gene_sets.get('Ferroptosis Suppressor Genes', [])
immune_genes = (gene_sets.get('T Cell Infiltration Markers', []) + 
                gene_sets.get('Antigen Presentation Genes', []) + 
                gene_sets.get('Immune Checkpoint Genes', []))

# Calculate scores
ferroscore = calculate_ferroscore(expr_matrix, driver_genes, suppressor_genes)
immune_score = calculate_immune_score(expr_matrix, immune_genes)

# Combined score
ferro_immuno_score = 0.6 * ferroscore + 0.4 * immune_score
ferro_immuno_score = pd.Series(ferro_immuno_score, index=expr_matrix.columns, name='FerroImmuno_Score')

# Combine
scores_df = pd.DataFrame({
    'FerroScore': ferroscore,
    'Immune_Score': immune_score,
    'FerroImmuno_Score': ferro_immuno_score
})

print(f"\n5. Score statistics:")
print(f"   FerroScore: [{ferroscore.min():.3f}, {ferroscore.max():.3f}]")
print(f"   Immune Score: [{immune_score.min():.3f}, {immune_score.max():.3f}]")
print(f"   FerroImmuno Score: [{ferro_immuno_score.min():.3f}, {ferro_immuno_score.max():.3f}]")

# Save
scores_df.to_csv(f'{RESULTS_DIR}/gse91061_ferro_immuno_scores.csv')
print(f"\n   [OK] Saved: gse91061_ferro_immuno_scores.csv")

print("\n" + "=" * 60)
print("Complete! Next: Extract response labels and combine with GSE78220")
print("=" * 60)
