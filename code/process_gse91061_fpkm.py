#!/usr/bin/env python3
"""
Process GSE91061 FPKM data and calculate FerroScore-Immuno
"""

import pandas as pd
import numpy as np
import gzip
import os

EXTERNAL_DIR = '../data/external'
RESULTS_DIR = '../results/external'
os.makedirs(RESULTS_DIR, exist_ok=True)

print("=" * 60)
print("Processing GSE91061 FPKM Data")
print("=" * 60)

# Load FPKM data
fpkm_file = f'{EXTERNAL_DIR}/GSE91061_fpkm.csv.gz'
print(f"\nLoading: {fpkm_file}")

# Read CSV (tab or comma separated)
with gzip.open(fpkm_file, 'rt', encoding='utf-8', errors='ignore') as f:
    # Check first line to determine separator
    first_line = f.readline()
    f.seek(0)
    
    if '\t' in first_line:
        fpkm_data = pd.read_csv(f, sep='\t')
    else:
        fpkm_data = pd.read_csv(f)

print(f"  FPKM data shape: {fpkm_data.shape}")
print(f"  Columns: {list(fpkm_data.columns[:5])}...")
print(f"  First few rows:")
print(fpkm_data.head())

# Load gene sets
print("\n" + "=" * 60)
print("Loading Gene Sets")
print("=" * 60)

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
    print(f"  {name}: {len(genes)} genes")

# Get all unique genes
all_genes = []
for genes in gene_sets.values():
    all_genes.extend(genes)
all_genes = list(set(all_genes))
print(f"\n  Total unique genes: {len(all_genes)}")

# Check which genes are in FPKM data
print("\n" + "=" * 60)
print("Matching Genes")
print("=" * 60)

# FPKM data should have gene symbols in first column
gene_col = fpkm_data.columns[0]
print(f"  Gene column: {gene_col}")

available_genes = [g for g in all_genes if g in fpkm_data[gene_col].values]
missing_genes = [g for g in all_genes if g not in fpkm_data[gene_col].values]

print(f"  Available genes: {len(available_genes)}/{len(all_genes)}")
print(f"  Missing genes: {len(missing_genes)}")
if missing_genes:
    print(f"    Missing: {missing_genes[:10]}...")

# Extract expression for available genes
print("\n" + "=" * 60)
print("Extracting Expression Matrix")
print("=" * 60)

expr_matrix = fpkm_data[fpkm_data[gene_col].isin(available_genes)]
expr_matrix = expr_matrix.set_index(gene_col)

# Remove gene symbol column, keep only sample columns
sample_cols = [c for c in expr_matrix.columns if c != gene_col]
expr_matrix = expr_matrix[sample_cols]

print(f"  Expression matrix: {expr_matrix.shape}")
print(f"  Samples: {len(sample_cols)}")
print(f"  Genes: {len(available_genes)}")

# Save extracted expression
expr_matrix.to_csv(f'{RESULTS_DIR}/gse91061_fpkm_extracted.csv')
print(f"  [OK] Saved: gse91061_fpkm_extracted.csv")

# Calculate FerroScore-Immuno
print("\n" + "=" * 60)
print("Calculating FerroScore-Immuno")
print("=" * 60)

from sklearn.preprocessing import MinMaxScaler

def calculate_ferroscore(expr_matrix, driver_genes, suppressor_genes):
    """Calculate FerroScore"""
    print("\n  Calculating FerroScore...")
    
    driver_in_data = [g for g in driver_genes if g in expr_matrix.index]
    suppressor_in_data = [g for g in suppressor_genes if g in expr_matrix.index]
    
    print(f"    Driver genes: {len(driver_in_data)}")
    print(f"    Suppressor genes: {len(suppressor_in_data)}")
    
    # Log transform FPKM (add 1 to avoid log(0))
    expr_log = np.log2(expr_matrix + 1)
    
    # Calculate scores
    driver_scores = expr_log.loc[driver_in_data].mean() if driver_in_data else pd.Series(0, index=expr_matrix.columns)
    suppressor_scores = expr_log.loc[suppressor_in_data].mean() if suppressor_in_data else pd.Series(0, index=expr_matrix.columns)
    
    # Combined score (high driver + low suppressor = high ferroptosis)
    ferroscore = driver_scores - suppressor_scores
    
    # Normalize to [0, 1]
    scaler = MinMaxScaler()
    ferroscore_scaled = scaler.fit_transform(ferroscore.values.reshape(-1, 1)).flatten()
    
    return pd.Series(ferroscore_scaled, index=expr_matrix.columns, name='FerroScore')

def calculate_immune_score(expr_matrix, immune_genes):
    """Calculate Immune Score"""
    print("\n  Calculating Immune Score...")
    
    immune_in_data = [g for g in immune_genes if g in expr_matrix.index]
    print(f"    Immune genes: {len(immune_in_data)}")
    
    # Log transform
    expr_log = np.log2(expr_matrix + 1)
    
    # Mean expression of immune genes
    immune_score = expr_log.loc[immune_in_data].mean() if immune_in_data else pd.Series(0, index=expr_matrix.columns)
    
    # Normalize
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

# Calculate combined score
ferro_immuno_score = 0.6 * ferroscore + 0.4 * immune_score
ferro_immuno_score = pd.Series(ferro_immuno_score, index=expr_matrix.columns, name='FerroImmuno_Score')

# Combine results
scores_df = pd.DataFrame({
    'FerroScore': ferroscore,
    'Immune_Score': immune_score,
    'FerroImmuno_Score': ferro_immuno_score
})

print(f"\n  FerroScore range: [{ferroscore.min():.3f}, {ferroscore.max():.3f}]")
print(f"  Immune Score range: [{immune_score.min():.3f}, {immune_score.max():.3f}]")
print(f"  FerroImmuno Score range: [{ferro_immuno_score.min():.3f}, {ferro_immuno_score.max():.3f}]")

# Save scores
scores_df.to_csv(f'{RESULTS_DIR}/gse91061_ferro_immuno_scores.csv')
print(f"\n  [OK] Saved: gse91061_ferro_immuno_scores.csv")

print("\n" + "=" * 60)
print("GSE91061 Processing Complete!")
print("=" * 60)
print("\nNext: Extract response labels and combine with GSE78220")
