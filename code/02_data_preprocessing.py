#!/usr/bin/env python3
"""
Step 2: Data Preprocessing for FerroScore-Radio (Real Data Version)
Process real TCGA data from Xena
"""

import os
import pandas as pd
import numpy as np
import gzip

# Paths
RAW_DIR = '../data/raw'
PROC_DIR = '../data/processed'
os.makedirs(PROC_DIR, exist_ok=True)

def load_gene_sets():
    """Load Ferro-Radio gene sets"""
    gene_file = '../gene_sets/ferro_radio_genes.txt'
    
    gene_sets = {}
    current_set = None
    
    with open(gene_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('##'):
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
    
    return gene_sets, unique_genes

def load_expression_data():
    """Load TCGA expression data"""
    print("\nLoading expression data...")
    
    # Try different possible filenames
    possible_files = [
        'tcga_RSEM_gene_tpm',
        'tcga_RSEM_gene_tpm.gz',
        'TCGA_HiSeqV2_PANCAN',
        'TCGA_HiSeqV2_PANCAN.gz'
    ]
    
    expr_file = None
    for f in possible_files:
        path = f'{RAW_DIR}/{f}'
        if os.path.exists(path):
            expr_file = path
            break
    
    if expr_file is None:
        print(f"  Expression file not found in {RAW_DIR}")
        print(f"  Available files: {os.listdir(RAW_DIR)}")
        return None
    
    print(f"  Loading: {expr_file}")
    
    # Load expression matrix
    if expr_file.endswith('.gz'):
        expr_df = pd.read_csv(expr_file, sep='\t', index_col=0, compression='gzip')
    else:
        expr_df = pd.read_csv(expr_file, sep='\t', index_col=0)
    
    print(f"  Expression matrix: {expr_df.shape}")
    
    return expr_df

def extract_gene_symbol(index):
    """Extract gene symbol from Ensembl ID or full name"""
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
        print(f"  Missing genes ({len(missing_genes)}): {missing_genes[:10]}...")
    
    if len(available_genes) == 0:
        print("  ERROR: No genes found!")
        return None, []
    
    # Subset expression
    expr_subset = expr_df.loc[available_genes]
    
    # Remove low expression genes (median < 0.1 TPM)
    median_expr = expr_subset.median(axis=1)
    expr_filtered = expr_subset[median_expr >= 0.1]
    print(f"  After filtering low expression: {expr_filtered.shape}")
    
    return expr_filtered, available_genes

def load_clinical_data():
    """Load clinical and survival data"""
    print("\nLoading clinical data...")
    
    # Try to find clinical file
    possible_clinical_files = [
        'TCGA_phenotype_denseDataOnlyDownload.tsv',
        'TCGA_phenotype_denseDataOnlyDownload.tsv.gz',
        'Survival_SupplementalTable_S1_20171025_xena_sp',
        'survivalTCGA_pancan.txt'
    ]
    
    clin_file = None
    for f in possible_clinical_files:
        path = f'{RAW_DIR}/{f}'
        if os.path.exists(path):
            clin_file = path
            break
    
    if clin_file is None:
        print(f"  Clinical file not found")
        return None
    
    print(f"  Loading: {clin_file}")
    
    # Load clinical data
    if clin_file.endswith('.gz'):
        clinical_df = pd.read_csv(clin_file, sep='\t', compression='gzip')
    else:
        clinical_df = pd.read_csv(clin_file, sep='\t')
    
    print(f"  Clinical data: {clinical_df.shape}")
    print(f"  Columns: {list(clinical_df.columns[:10])}...")
    
    return clinical_df

def main():
    """Main preprocessing function"""
    print("=" * 60)
    print("FerroScore-Radio: Data Preprocessing (Real Data)")
    print("=" * 60)
    
    # 1. Load gene sets
    gene_sets, unique_genes = load_gene_sets()
    
    # 2. Load expression data
    tcga_expr = load_expression_data()
    if tcga_expr is None:
        print("\nERROR: No expression data found!")
        print("Please download TCGA expression data first.")
        return
    
    # 3. Preprocess expression
    tcga_processed, available_genes = preprocess_expression(tcga_expr, unique_genes)
    if tcga_processed is None:
        return
    
    # Save processed expression
    tcga_processed.to_csv(f'{PROC_DIR}/tcga_ferro_radio_expression.csv')
    print(f"  Saved: tcga_ferro_radio_expression.csv")
    
    # 4. Load clinical data
    clinical = load_clinical_data()
    
    if clinical is not None:
        # Rename columns to standard names
        col_mapping = {}
        
        # Find sample ID column
        sample_col = None
        for col in clinical.columns:
            if 'sample' in col.lower() or 'bcr_patient_barcode' in col.lower():
                sample_col = col
                col_mapping[col] = 'sample_id'
                break
        
        # Find other important columns
        for col in clinical.columns:
            col_lower = col.lower()
            if 'age' in col_lower:
                col_mapping[col] = 'age'
            elif 'gender' in col_lower or 'sex' in col_lower:
                col_mapping[col] = 'gender'
            elif 'stage' in col_lower:
                col_mapping[col] = 'stage'
            elif 'cancer' in col_lower or 'type' in col_lower:
                col_mapping[col] = 'cancer_type'
            elif 'os.time' in col_lower or 'os_time' in col_lower:
                col_mapping[col] = 'OS.time'
            elif col_lower == 'os':
                col_mapping[col] = 'OS'
            elif 'pfi.time' in col_lower:
                col_mapping[col] = 'PFI.time'
            elif col_lower == 'pfi':
                col_mapping[col] = 'PFI'
        
        # Rename columns
        clinical_renamed = clinical.rename(columns=col_mapping)
        
        # Save
        clinical_renamed.to_csv(f'{PROC_DIR}/tcga_clinical.csv', index=False)
        print(f"  Saved: tcga_clinical.csv")
        
        # Also save as survival data if survival columns exist
        survival_cols = ['sample_id', 'OS', 'OS.time', 'PFI', 'PFI.time']
        available_surv_cols = [c for c in survival_cols if c in clinical_renamed.columns]
        if len(available_surv_cols) >= 3:
            survival_df = clinical_renamed[available_surv_cols].copy()
            survival_df.to_csv(f'{PROC_DIR}/tcga_survival.csv', index=False)
            print(f"  Saved: tcga_survival.csv")
    
    # 5. Save gene list
    pd.DataFrame({'gene': available_genes}).to_csv(
        f'{PROC_DIR}/ferro_radio_gene_list.csv', index=False
    )
    print(f"  Saved: ferro_radio_gene_list.csv ({len(available_genes)} genes)")
    
    print("\n" + "=" * 60)
    print("Preprocessing Complete!")
    print("=" * 60)
    print("\nNext step: Run 03_ferroscore_algorithm.py")

if __name__ == "__main__":
    main()
