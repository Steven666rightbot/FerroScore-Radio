#!/usr/bin/env python3
"""
Step 2: Data Preprocessing for FerroScore-Radio
Process TCGA/GTEx expression data and clinical information
"""

import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import re

# Paths
RAW_DIR = '../data/raw'
PROC_DIR = '../data/processed'
os.makedirs(PROC_DIR, exist_ok=True)

def load_gene_sets():
    """Load Ferro-Radio gene sets"""
    gene_file = '../gene_sets/ferro_radio_genes.txt'
    
    gene_sets = {}
    current_set = None
    
    with open(gene_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('## '):
                current_set = line.replace('## ', '').split('(')[0].strip()
                gene_sets[current_set] = []
            elif line and not line.startswith('#') and current_set:
                gene_sets[current_set].append(line)
    
    # Get unique genes
    all_genes = []
    for genes in gene_sets.values():
        all_genes.extend(genes)
    unique_genes = list(set(all_genes))
    
    print(f"Loaded {len(unique_genes)} unique genes from {len(gene_sets)} categories")
    print(f"Categories: {list(gene_sets.keys())}")
    
    return gene_sets, unique_genes

def load_expression_data():
    """Load and merge TCGA + GTEx expression data"""
    print("\nLoading expression data...")
    
    # TCGA
    tcga_file = f'{RAW_DIR}/tcga_RSEM_gene_tpm'
    if os.path.exists(tcga_file):
        tcga_expr = pd.read_csv(tcga_file, sep='\t', index_col=0)
        print(f"  TCGA: {tcga_expr.shape}")
    else:
        print(f"  ✗ TCGA file not found: {tcga_file}")
        tcga_expr = None
    
    # GTEx
    gtex_file = f'{RAW_DIR}/gtex_RSEM_gene_tpm'
    if os.path.exists(gtex_file):
        gtex_expr = pd.read_csv(gtex_file, sep='\t', index_col=0)
        print(f"  GTEx: {gtex_expr.shape}")
    else:
        print(f"  ✗ GTEx file not found: {gtex_file}")
        gtex_expr = None
    
    return tcga_expr, gtex_expr

def extract_gene_symbol(index):
    """Extract gene symbol from Ensembl ID or full name"""
    # Handle different formats
    if isinstance(index, str):
        # Remove Ensembl ID part if present
        if '|' in index:
            return index.split('|')[0]
        # Remove version number
        if index.startswith('ENSG'):
            return index.split('.')[0]
    return index

def preprocess_expression(expr_df, gene_list):
    """Preprocess expression matrix"""
    print(f"\nPreprocessing expression data...")
    print(f"  Original shape: {expr_df.shape}")
    
    # Extract gene symbols
    expr_df.index = [extract_gene_symbol(idx) for idx in expr_df.index]
    
    # Remove duplicates (keep first)
    expr_df = expr_df[~expr_df.index.duplicated(keep='first')]
    print(f"  After removing duplicates: {expr_df.shape}")
    
    # Filter to Ferro-Radio genes
    available_genes = [g for g in gene_list if g in expr_df.index]
    missing_genes = [g for g in gene_list if g not in expr_df.index]
    
    print(f"  Available Ferro-Radio genes: {len(available_genes)}/{len(gene_list)}")
    if missing_genes:
        print(f"  Missing genes: {missing_genes[:10]}...")
    
    # Subset expression
    expr_subset = expr_df.loc[available_genes]
    
    # Remove low expression genes (median < 0.1 TPM)
    median_expr = expr_subset.median(axis=1)
    expr_filtered = expr_subset[median_expr >= 0.1]
    print(f"  After filtering low expression: {expr_filtered.shape}")
    
    # Log2 transform (if not already)
    # Xena data is usually already log2(TPM+0.001)
    
    return expr_filtered, available_genes

def load_clinical_data():
    """Load and process clinical data"""
    print("\nLoading clinical data...")
    
    # TCGA phenotype
    tcga_clin_file = f'{RAW_DIR}/TCGA_phenotype.tsv'
    if os.path.exists(tcga_clin_file):
        tcga_clin = pd.read_csv(tcga_clin_file, sep='\t')
        print(f"  TCGA clinical: {tcga_clin.shape}")
        
        # Key columns to extract
        key_cols = ['sample', '_primary_site', '_primary_disease', 
                   'age_at_initial_pathologic_diagnosis', 'gender',
                   'race', 'ethnicity', 'ajcc_pathologic_tumor_stage',
                   'radiation_therapy', 'targeted_molecular_therapy']
        
        available_cols = [c for c in key_cols if c in tcga_clin.columns]
        tcga_clin_subset = tcga_clin[available_cols].copy()
        
        # Rename columns
        col_mapping = {
            'sample': 'sample_id',
            '_primary_site': 'tissue',
            '_primary_disease': 'cancer_type',
            'ajcc_pathologic_tumor_stage': 'stage',
            'radiation_therapy': 'radiotherapy',
            'targeted_molecular_therapy': 'targeted_therapy'
        }
        tcga_clin_subset = tcga_clin_subset.rename(columns=col_mapping)
        
    else:
        print(f"  ✗ TCGA clinical file not found")
        tcga_clin_subset = None
    
    # Survival data
    surv_file = f'{RAW_DIR}/TCGA_survival.txt'
    if os.path.exists(surv_file):
        surv_data = pd.read_csv(surv_file, sep='\t')
        print(f"  Survival data: {surv_data.shape}")
    else:
        print(f"  ✗ Survival file not found")
        surv_data = None
    
    return tcga_clin_subset, surv_data

def identify_radiotherapy_patients(clinical_df):
    """Identify patients who received radiotherapy"""
    if clinical_df is None or 'radiotherapy' not in clinical_df.columns:
        return None
    
    rt_col = 'radiotherapy'
    
    # Check unique values
    unique_vals = clinical_df[rt_col].dropna().unique()
    print(f"\nRadiotherapy values: {unique_vals}")
    
    # Define positive indicators
    rt_positive = ['yes', 'YES', 'Yes', 'y', 'Y', '1', 'true', 'TRUE', True]
    
    clinical_df['received_radiotherapy'] = clinical_df[rt_col].isin(rt_positive)
    rt_count = clinical_df['received_radiotherapy'].sum()
    print(f"  Patients with radiotherapy: {rt_count}/{len(clinical_df)}")
    
    return clinical_df

def main():
    """Main preprocessing function"""
    print("=" * 60)
    print("FerroScore-Radio: Data Preprocessing")
    print("=" * 60)
    
    # 1. Load gene sets
    gene_sets, unique_genes = load_gene_sets()
    
    # 2. Load expression data
    tcga_expr, gtex_expr = load_expression_data()
    
    # 3. Preprocess TCGA
    if tcga_expr is not None:
        tcga_processed, available_genes = preprocess_expression(tcga_expr, unique_genes)
        tcga_processed.to_csv(f'{PROC_DIR}/tcga_ferro_radio_expression.csv')
        print(f"  ✓ Saved: tcga_ferro_radio_expression.csv")
    
    # 4. Preprocess GTEx
    if gtex_expr is not None:
        gtex_processed, _ = preprocess_expression(gtex_expr, unique_genes)
        gtex_processed.to_csv(f'{PROC_DIR}/gtex_ferro_radio_expression.csv')
        print(f"  ✓ Saved: gtex_ferro_radio_expression.csv")
    
    # 5. Load and process clinical data
    tcga_clin, surv_data = load_clinical_data()
    
    if tcga_clin is not None:
        tcga_clin = identify_radiotherapy_patients(tcga_clin)
        tcga_clin.to_csv(f'{PROC_DIR}/tcga_clinical.csv', index=False)
        print(f"  ✓ Saved: tcga_clinical.csv")
    
    if surv_data is not None:
        surv_data.to_csv(f'{PROC_DIR}/tcga_survival.csv', index=False)
        print(f"  ✓ Saved: tcga_survival.csv")
    
    # 6. Save gene list
    pd.DataFrame({'gene': available_genes}).to_csv(
        f'{PROC_DIR}/ferro_radio_gene_list.csv', index=False
    )
    print(f"  ✓ Saved: ferro_radio_gene_list.csv ({len(available_genes)} genes)")
    
    print("\n" + "=" * 60)
    print("Preprocessing Complete!")
    print("=" * 60)
    print("\nNext step: Run 03_ferroscore_algorithm.py")

if __name__ == "__main__":
    main()
