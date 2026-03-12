#!/usr/bin/env python3
"""
Download NCBI gene_info and convert GSE91061 gene IDs
"""

import pandas as pd
import gzip
import os
import requests
from io import BytesIO

EXTERNAL_DIR = '../data/external'
RESULTS_DIR = '../results/external'
os.makedirs(EXTERNAL_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

print("=" * 60)
print("Downloading NCBI Gene Info for ID Conversion")
print("=" * 60)

# Download gene_info.gz from NCBI
print("\n1. Downloading gene_info from NCBI...")
print("   This file is ~500MB compressed, may take a few minutes...")

# Try HTTP mirror instead of FTP
url = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz'
output_file = f'{EXTERNAL_DIR}/gene_info.gz'

if os.path.exists(output_file):
    print(f"   File already exists: {output_file}")
else:
    try:
        print(f"   Downloading from: {url}")
        response = requests.get(url, timeout=600, stream=True)
        
        if response.status_code == 200:
            with open(output_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            file_size = os.path.getsize(output_file) / 1024 / 1024
            print(f"   [OK] Downloaded: {output_file} ({file_size:.1f} MB)")
        else:
            print(f"   [X] Status: {response.status_code}")
    except Exception as e:
        print(f"   [X] Error: {e}")
        print("\n   Please download manually from:")
        print("   ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz")
        exit(1)

# Parse gene_info for human genes only
print("\n2. Parsing gene_info for human genes...")

# Read only necessary columns: tax_id, GeneID, Symbol
gene_info_file = f'{EXTERNAL_DIR}/gene_info.gz'

# Read first few lines to check format
with gzip.open(gene_info_file, 'rt', encoding='utf-8', errors='ignore') as f:
    header = f.readline().strip().split('\t')
    print(f"   Columns: {header[:5]}...")

# Read full file (this may take a while)
print("   Reading gene_info (this may take a few minutes)...")
gene_info = pd.read_csv(gene_info_file, sep='\t', usecols=['#tax_id', 'GeneID', 'Symbol'], 
                        dtype={'#tax_id': int, 'GeneID': int, 'Symbol': str})

# Filter for human genes only (tax_id = 9606)
human_genes = gene_info[gene_info['#tax_id'] == 9606].copy()
human_genes = human_genes.rename(columns={'#tax_id': 'tax_id', 'GeneID': 'entrez_id', 'Symbol': 'gene_symbol'})

print(f"   Total genes: {len(gene_info)}")
print(f"   Human genes: {len(human_genes)}")
print(f"   First few: {human_genes.head()}")

# Create mapping dictionary
print("\n3. Creating ID mapping dictionary...")
id_to_symbol = dict(zip(human_genes['entrez_id'].astype(str), human_genes['gene_symbol']))
print(f"   Mapping size: {len(id_to_symbol)}")

# Load GSE91061 and apply mapping
print("\n4. Loading GSE91061 and applying mapping...")

fpkm_file = f'{EXTERNAL_DIR}/GSE91061_fpkm.csv.gz'
with gzip.open(fpkm_file, 'rt', encoding='utf-8', errors='ignore') as f:
    fpkm_data = pd.read_csv(f)

print(f"   GSE91061 shape: {fpkm_data.shape}")

# Get gene IDs (first column)
gene_ids = fpkm_data.iloc[:, 0].astype(str)
print(f"   First few IDs: {gene_ids.head().tolist()}")

# Map to symbols
fpkm_data['gene_symbol'] = gene_ids.map(id_to_symbol)

# Check mapping success
mapped_count = fpkm_data['gene_symbol'].notna().sum()
print(f"\n   Mapped genes: {mapped_count}/{len(fpkm_data)} ({mapped_count/len(fpkm_data)*100:.1f}%)")

# Filter to mapped genes only
fpkm_mapped = fpkm_data[fpkm_data['gene_symbol'].notna()].copy()
fpkm_mapped = fpkm_mapped.set_index('gene_symbol')
fpkm_mapped = fpkm_mapped.drop(fpkm_mapped.columns[0], axis=1)  # Drop original ID column

print(f"   Final shape: {fpkm_mapped.shape}")
print(f"   First few genes: {fpkm_mapped.index[:10].tolist()}")

# Save mapped data
fpkm_mapped.to_csv(f'{RESULTS_DIR}/gse91061_fpkm_mapped.csv')
print(f"\n   [OK] Saved: gse91061_fpkm_mapped.csv")

# Save mapping for reference
mapping_df = pd.DataFrame({
    'entrez_id': gene_ids,
    'gene_symbol': fpkm_data['gene_symbol']
})
mapping_df.to_csv(f'{RESULTS_DIR}/gse91061_gene_id_mapping.csv', index=False)
print(f"   [OK] Saved: gse91061_gene_id_mapping.csv")

print("\n" + "=" * 60)
print("Conversion Complete!")
print("=" * 60)
print("\nNext: Calculate FerroScore-Immuno for GSE91061")
