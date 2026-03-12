#!/usr/bin/env python3
"""
Parse GEO series matrix files to extract expression data
"""

import pandas as pd
import numpy as np
import gzip
import os

EXTERNAL_DIR = '../data/external'
RESULTS_DIR = '../results/external'
os.makedirs(RESULTS_DIR, exist_ok=True)

def parse_geo_series_matrix(filepath):
    """Parse GEO series matrix file to extract expression data"""
    
    print(f"Parsing: {filepath}")
    
    metadata = {}
    expression_data = []
    in_expression_section = False
    sample_ids = []
    
    with gzip.open(filepath, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            
            if not line:
                continue
            
            # Parse metadata
            if line.startswith('!Series_'):
                parts = line.split('\t')
                key = parts[0].replace('!Series_', '')
                values = [p.strip('"') for p in parts[1:]]
                metadata[key] = values
            
            # Parse sample metadata
            elif line.startswith('!Sample_'):
                parts = line.split('\t')
                key = parts[0].replace('!Sample_', '')
                values = [p.strip('"') for p in parts[1:]]
                if key not in metadata:
                    metadata[key] = {}
                metadata[key] = values
            
            # Start of expression data
            elif line == '!series_matrix_table_begin':
                in_expression_section = True
                # Next line should be header
                header_line = next(f).strip()
                headers = header_line.split('\t')
                sample_ids = [h.strip('"') for h in headers[1:]]
                print(f"  Found {len(sample_ids)} samples")
                continue
            
            # End of expression data
            elif line == '!series_matrix_table_end':
                in_expression_section = False
                break
            
            # Parse expression values
            elif in_expression_section:
                parts = line.split('\t')
                gene_id = parts[0].strip('"')
                values = []
                for p in parts[1:]:
                    try:
                        values.append(float(p))
                    except:
                        values.append(np.nan)
                expression_data.append([gene_id] + values)
    
    # Create DataFrame
    if expression_data and sample_ids:
        df = pd.DataFrame(expression_data, columns=['Gene'] + sample_ids)
        df = df.set_index('Gene')
        print(f"  Expression matrix: {df.shape}")
        return df, metadata
    else:
        print("  No expression data found")
        return None, metadata

print("=" * 60)
print("Parsing GEO Expression Data")
print("=" * 60)

# Parse GSE78220
print("\n1. GSE78220")
expr_78220, meta_78220 = parse_geo_series_matrix(f'{EXTERNAL_DIR}/GSE78220_series_matrix_full.txt.gz')

if expr_78220 is not None:
    # Save expression
    expr_78220.to_csv(f'{RESULTS_DIR}/gse78220_expression_matrix.csv')
    print(f"  [OK] Saved: gse78220_expression_matrix.csv")
    
    # Print sample info
    print(f"  Genes: {expr_78220.shape[0]}")
    print(f"  Samples: {expr_78220.shape[1]}")
    print(f"  Expression range: [{expr_78220.min().min():.2f}, {expr_78220.max().max():.2f}]")

# Parse GSE91061
print("\n2. GSE91061")
expr_91061, meta_91061 = parse_geo_series_matrix(f'{EXTERNAL_DIR}/GSE91061_series_matrix_full.txt.gz')

if expr_91061 is not None:
    # Save expression
    expr_91061.to_csv(f'{RESULTS_DIR}/gse91061_expression_matrix.csv')
    print(f"  [OK] Saved: gse91061_expression_matrix.csv")
    
    # Print sample info
    print(f"  Genes: {expr_91061.shape[0]}")
    print(f"  Samples: {expr_91061.shape[1]}")
    print(f"  Expression range: [{expr_91061.min().min():.2f}, {expr_91061.max().max():.2f}]")

print("\n" + "=" * 60)
print("Parsing Complete!")
print("=" * 60)
print("\nNext: Map gene symbols and calculate FerroScore-Immuno")
