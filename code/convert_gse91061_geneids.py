#!/usr/bin/env python3
"""
Convert Entrez Gene ID to Gene Symbol for GSE91061
Using NCBI gene_info or alternative mapping
"""

import pandas as pd
import numpy as np
import gzip
import os
import requests
from io import StringIO

EXTERNAL_DIR = '../data/external'
RESULTS_DIR = '../results/external'
os.makedirs(RESULTS_DIR, exist_ok=True)

print("=" * 60)
print("Converting GSE91061 Entrez ID to Gene Symbol")
print("=" * 60)

# Load GSE91061 FPKM data
print("\n1. Loading GSE91061 FPKM data...")
fpkm_file = f'{EXTERNAL_DIR}/GSE91061_fpkm.csv.gz'

with gzip.open(fpkm_file, 'rt', encoding='utf-8', errors='ignore') as f:
    fpkm_data = pd.read_csv(f)

print(f"   Shape: {fpkm_data.shape}")
print(f"   First few gene IDs: {fpkm_data.iloc[:5, 0].tolist()}")

# Get unique Entrez IDs
entrez_ids = fpkm_data.iloc[:, 0].astype(str).tolist()
print(f"\n   Total genes: {len(entrez_ids)}")

# Method 1: Try to download gene_info from NCBI
print("\n2. Attempting to download gene ID mapping...")

# Alternative: Use mygene.info API (free, no registration)
print("\n   Using mygene.info API...")

def batch_query_mygene(entrez_ids, batch_size=1000):
    """Query mygene.info in batches"""
    import requests
    import json
    
    base_url = "https://mygene.info/v3/gene"
    mapping = {}
    
    for i in range(0, len(entrez_ids), batch_size):
        batch = entrez_ids[i:i+batch_size]
        ids_str = ','.join(batch)
        
        url = f"{base_url}?ids={ids_str}&fields=symbol,name"
        try:
            response = requests.get(url, timeout=60)
            if response.status_code == 200:
                data = response.json()
                for item in data:
                    query = item.get('query', '')
                    symbol = item.get('symbol', '')
                    if symbol:
                        mapping[query] = symbol
            print(f"   Processed {min(i+batch_size, len(entrez_ids))}/{len(entrez_ids)} genes...")
        except Exception as e:
            print(f"   Error in batch {i}: {e}")
    
    return mapping

# Query mygene.info (use subset for testing first)
print("\n   Querying mygene.info for gene symbols...")
test_size = min(100, len(entrez_ids))  # Start with 100 genes
test_ids = entrez_ids[:test_size]

mapping_test = batch_query_mygene(test_ids, batch_size=100)
print(f"\n   Mapped {len(mapping_test)}/{test_size} test genes")
print(f"   Examples: {list(mapping_test.items())[:5]}")

# If successful, do all genes
if len(mapping_test) > test_size * 0.5:  # If >50% success rate
    print("\n   Success rate good, processing all genes...")
    full_mapping = batch_query_mygene(entrez_ids, batch_size=1000)
    print(f"\n   Total mapped: {len(full_mapping)}/{len(entrez_ids)}")
    
    # Save mapping
    mapping_df = pd.DataFrame(list(full_mapping.items()), columns=['entrez_id', 'gene_symbol'])
    mapping_df.to_csv(f'{RESULTS_DIR}/gse91061_gene_id_mapping.csv', index=False)
    print(f"   [OK] Saved: gse91061_gene_id_mapping.csv")
    
    # Apply mapping to FPKM data
    print("\n3. Applying gene symbol mapping...")
    fpkm_data['gene_symbol'] = fpkm_data.iloc[:, 0].astype(str).map(full_mapping)
    
    # Check how many genes we can map
    mapped_count = fpkm_data['gene_symbol'].notna().sum()
    print(f"   Genes with symbol: {mapped_count}/{len(fpkm_data)}")
    
    # Filter to genes with symbols
    fpkm_mapped = fpkm_data[fpkm_data['gene_symbol'].notna()].copy()
    fpkm_mapped = fpkm_mapped.set_index('gene_symbol')
    fpkm_mapped = fpkm_mapped.drop(fpkm_mapped.columns[0], axis=1)  # Drop original ID column
    
    print(f"   Final shape: {fpkm_mapped.shape}")
    
    # Save mapped data
    fpkm_mapped.to_csv(f'{RESULTS_DIR}/gse91061_fpkm_mapped.csv')
    print(f"   [OK] Saved: gse91061_fpkm_mapped.csv")
    
    print("\n" + "=" * 60)
    print("Conversion Complete!")
    print("=" * 60)
    print("\nNext: Calculate FerroScore-Immuno for GSE91061")
    
else:
    print("\n   [X] Mapping success rate too low, trying alternative method...")
    print("\n   Alternative: Download NCBI gene_info file")
    print("   Visit: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz")
    print("   Or use: https://www.ncbi.nlm.nih.gov/gene")
