#!/usr/bin/env python3
"""
Convert Affymetrix Probe ID to Gene Symbol for GSE35452
"""

import pandas as pd
import numpy as np
import gzip
import os

EXTERNAL_DIR = '../data/external'

print("=" * 60)
print("Converting Probe ID to Gene Symbol")
print("=" * 60)

# Read GPL570 annotation
print("\nReading GPL570 annotation...")
gpl_file = f'{EXTERNAL_DIR}/GPL570-55999.txt'

# Read first few lines to check format
with open(gpl_file, 'r') as f:
    for i in range(10):
        line = f.readline().strip()
        print(line)

print("\nLoading full annotation...")
# Skip comment lines
annotation = pd.read_csv(gpl_file, sep='\t', skiprows=16)
print(f"Annotation shape: {annotation.shape}")
print(f"Columns: {list(annotation.columns)}")

# Check columns
if 'ID' in annotation.columns and 'Gene Symbol' in annotation.columns:
    probe_to_gene = dict(zip(annotation['ID'], annotation['Gene Symbol']))
    print(f"Created mapping for {len(probe_to_gene)} probes")
    
    # Check some examples
    print("\nExample mappings:")
    for probe, gene in list(probe_to_gene.items())[:10]:
        if pd.notna(gene) and gene != '':
            print(f"  {probe} -> {gene}")
else:
    print("ERROR: Required columns not found!")
    print("Available columns:", annotation.columns.tolist())
