#!/usr/bin/env python3
"""
Process downloaded gene_info and convert GSE91061 gene IDs
"""

import pandas as pd
import gzip
import os

EXTERNAL_DIR = '../data/external'
RESULTS_DIR = '../results/external'

print("=" * 60)
print("Processing Gene Info and Converting GSE91061")
print("=" * 60)

# Parse gene_info for human genes only
print("\n1. Parsing gene_info for human genes...")

gene_info_file = f'{EXTERNAL_DIR}/gene_info.gz'
file_size = os.path.getsize(gene_info_file) / 1024 / 1024 / 1024
print(f"   File size: {file_size:.2f} GB")

# Read first few lines to check format
with gzip.open(gene_info_file, 'rt', encoding='utf-8', errors='ignore') as f:
    header = f.readline().strip().split('\t')
    print(f"   Columns: {header[:5]}...")

# Read full file (this may take a few minutes)
print("\n   Reading gene_info (this may take a few minutes)...")
print("   Processing in chunks to save memory...")

chunk_size = 100000
human_genes_list = []

for chunk in pd.read_csv(gene_info_file, sep='\t', usecols=['#tax_id', 'GeneID', 'Symbol'], 
                         dtype={'#tax_id': int, 'GeneID': int, 'Symbol': str},
                         chunksize=chunk_size):
    # Filter for human genes only (tax_id = 9606)
    human_chunk = chunk[chunk['#tax_id'] == 9606].copy()
    if len(human_chunk) > 0:
        human_genes_list.append(human_chunk)
    print(f"   Processed chunk...", end='\r')

print(f"\n   Concatenating {len(human_genes_list)} chunks...")
human_genes = pd.concat(human_genes_list, ignore_index=True)
human_genes = human_genes.rename(columns={'#tax_id': 'tax_id', 'GeneID': 'entrez_id', 'Symbol': 'gene_symbol'})

print(f"   Total human genes: {len(human_genes)}")
print(f"   First few: {human_genes.head()}")

# Create mapping dictionary
print("\n2. Creating ID mapping dictionary...")
id_to_symbol = dict(zip(human_genes['entrez_id'].astype(str), human_genes['gene_symbol']))
print(f"   Mapping size: {len(id_to_symbol)}")

# Load GSE91061 and apply mapping
print("\n3. Loading GSE91061 and applying mapping...")

fpkm_file = f'{EXTERNAL_DIR}/GSE91061_fpkm.csv.gz'
with gzip.open(fpkm_file, 'rt', encoding='utf-8', errors='ignore') as f:
    fpkm_data = pd.read_csv(f)

print(f"   GSE91061 shape: {fpkm_data.shape}")

# Get gene IDs (first column)
gene_ids = fpkm_data.iloc[:, 0].astype(str)
print(f"   First few IDs: {gene_ids.head().tolist()}")

# Map to symbols
print("\n   Mapping gene IDs to symbols...")
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
print("\n4. Saving mapped data...")
fpkm_mapped.to_csv(f'{RESULTS_DIR}/gse91061_fpkm_mapped.csv')
print(f"   [OK] Saved: gse91061_fpkm_mapped.csv")

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
