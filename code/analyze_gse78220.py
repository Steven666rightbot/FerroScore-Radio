#!/usr/bin/env python3
"""
Analyze GSE78220 - Melanoma Anti-PD-1 (Pembrolizumab)
This is the real immunotherapy dataset
"""

import pandas as pd
import numpy as np
import gzip
import os

EXTERNAL_DIR = '../data/external'
RESULTS_DIR = '../results/external'
os.makedirs(RESULTS_DIR, exist_ok=True)

print("=" * 60)
print("GSE78220 Analysis - Melanoma Anti-PD-1")
print("=" * 60)

# Read GSE78220 data
file_path = f'{EXTERNAL_DIR}/GSE78220_series_matrix.txt.gz'

if not os.path.exists(file_path):
    print(f"File not found: {file_path}")
    exit(1)

print(f"\nReading: {file_path}")

# Parse series matrix file
sample_metadata = {}
expression_rows = []
in_data_section = False
sample_ids = []

with gzip.open(file_path, 'rt', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        
        if not line:
            continue
            
        # Parse metadata
        if line.startswith('!Sample_'):
            parts = line.split('\t')
            key = parts[0].replace('!Sample_', '')
            values = [p.strip('"') for p in parts[1:]]
            sample_metadata[key] = values
            
        # Start of expression data
        elif line.startswith('"ID_REF"'):
            in_data_section = True
            header = line.split('\t')
            sample_ids = [h.strip('"') for h in header[1:]]
            
        elif in_data_section and not line.startswith('!'):
            parts = line.split('\t')
            gene_id = parts[0].strip('"')
            values = [float(p) if p != 'null' else np.nan for p in parts[1:]]
            expression_rows.append([gene_id] + values)

print(f"Samples: {len(sample_ids)}")
print(f"Metadata fields: {list(sample_metadata.keys())[:10]}")

# Create expression DataFrame
expr_df = pd.DataFrame(expression_rows, columns=['Gene'] + sample_ids)
expr_df = expr_df.set_index('Gene')
print(f"Expression matrix: {expr_df.shape}")

# Extract sample characteristics and response
print("\n" + "=" * 60)
print("Extracting Response Labels")
print("=" * 60)

response_labels = {}
characteristics = sample_metadata.get('characteristics_ch1', [])

print(f"Sample characteristics (first 5):")
for i, char in enumerate(characteristics[:5]):
    print(f"  {sample_ids[i]}: {char}")

# Try to identify response from source_name (contains response info)
if 'source_name_ch1' in sample_metadata:
    sources = sample_metadata['source_name_ch1']
    for i, source in enumerate(sources):
        sample_id = sample_ids[i]
        source_lower = source.lower()
        
        # Look for response indicators in source name
        if 'complete response' in source_lower or 'partial response' in source_lower:
            response_labels[sample_id] = 1  # Responder (CR/PR)
        elif 'progressive disease' in source_lower or 'stable disease' in source_lower:
            response_labels[sample_id] = 0  # Non-responder (PD/SD)
        else:
            response_labels[sample_id] = -1  # Unknown

responder_count = sum(1 for v in response_labels.values() if v == 1)
non_responder_count = sum(1 for v in response_labels.values() if v == 0)
unknown_count = sum(1 for v in response_labels.values() if v == -1)

print(f"\nResponder: {responder_count}")
print(f"Non-responder: {non_responder_count}")
print(f"Unknown: {unknown_count}")

# Check sample titles for more info
if 'title' in sample_metadata:
    titles = sample_metadata['title']
    print(f"\nSample titles (first 10):")
    for i, title in enumerate(titles[:10]):
        print(f"  {sample_ids[i]}: {title}")

# Check source name
if 'source_name_ch1' in sample_metadata:
    sources = sample_metadata['source_name_ch1']
    print(f"\nSource names (first 10):")
    for i, source in enumerate(sources[:10]):
        print(f"  {sample_ids[i]}: {source}")

# Save processed data
if response_labels:
    response_df = pd.DataFrame([
        {'sample_id': k, 'response': v} 
        for k, v in response_labels.items()
    ])
    response_df.to_csv(f'{RESULTS_DIR}/gse78220_response.csv', index=False)
    print(f"\n[OK] Saved: gse78220_response.csv")

# Save expression
expr_df.to_csv(f'{RESULTS_DIR}/gse78220_expression.csv')
print(f"[OK] Saved: gse78220_expression.csv")

print("\n" + "=" * 60)
print("GSE78220 Processing Complete!")
print("=" * 60)
