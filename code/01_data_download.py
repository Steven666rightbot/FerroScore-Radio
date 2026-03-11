#!/usr/bin/env python3
"""
Step 1: Data Download for FerroScore-Radio
Download TCGA and GTEx data from UCSC Xena
"""

import os
import pandas as pd
import numpy as np
from urllib.request import urlretrieve
from urllib.error import URLError
import gzip
import shutil

# Create directories
os.makedirs('../data/raw', exist_ok=True)
os.makedirs('../data/processed', exist_ok=True)

# UCSC Xena data hub URLs
XENA_HUB = "https://toil-xenagoatskins.s3.amazonaws.com"

# TCGA Pan-Cancer (PANCAN) datasets
DATASETS = {
    'tcga_expression': {
        'url': f'{XENA_HUB}/tcga_RSEM_gene_tpm.gz',
        'filename': 'tcga_RSEM_gene_tpm.gz',
        'description': 'TCGA RNA-seq TPM (log2)'
    },
    'tcga_phenotype': {
        'url': f'{XENA_HUB}/TCGA_phenotype_denseDataOnlyDownload.tsv.gz',
        'filename': 'TCGA_phenotype.tsv.gz',
        'description': 'TCGA clinical data'
    },
    'gtex_expression': {
        'url': f'{XENA_HUB}/gtex_RSEM_gene_tpm.gz',
        'filename': 'gtex_RSEM_gene_tpm.gz',
        'description': 'GTEx RNA-seq TPM (log2)'
    },
    'gtex_phenotype': {
        'url': f'{XENA_HUB}/GTEX_phenotype.tsv.gz',
        'filename': 'GTEX_phenotype.tsv.gz',
        'description': 'GTEx clinical data'
    },
    'survival_data': {
        'url': 'https://xenagoatskins.s3.amazonaws.com/survivalTCGA_pancan.txt',
        'filename': 'TCGA_survival.txt',
        'description': 'TCGA Pan-Cancer survival'
    }
}

def download_file(url, filename, description):
    """Download file with progress"""
    filepath = f'../data/raw/{filename}'
    
    if os.path.exists(filepath):
        print(f"✓ {description} already exists")
        return filepath
    
    print(f"Downloading {description}...")
    print(f"  URL: {url}")
    
    try:
        urlretrieve(url, filepath)
        print(f"  ✓ Saved to {filepath}")
        return filepath
    except URLError as e:
        print(f"  ✗ Error: {e}")
        return None

def decompress_gz(filepath):
    """Decompress .gz file"""
    if filepath.endswith('.gz'):
        output_path = filepath[:-3]
        if os.path.exists(output_path):
            print(f"  ✓ Decompressed file already exists")
            return output_path
        
        print(f"  Decompressing...")
        with gzip.open(filepath, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"  ✓ Decompressed to {output_path}")
        return output_path
    return filepath

def main():
    """Main download function"""
    print("=" * 60)
    print("FerroScore-Radio: Data Download")
    print("=" * 60)
    
    downloaded_files = {}
    
    for key, info in DATASETS.items():
        print(f"\n[{key}]")
        filepath = download_file(info['url'], info['filename'], info['description'])
        if filepath:
            downloaded_files[key] = filepath
            # Decompress if needed
            if filepath.endswith('.gz'):
                decompressed = decompress_gz(filepath)
                downloaded_files[key + '_decompressed'] = decompressed
    
    print("\n" + "=" * 60)
    print("Download Summary")
    print("=" * 60)
    for key, filepath in downloaded_files.items():
        if os.path.exists(filepath):
            size = os.path.getsize(filepath) / (1024**3)  # GB
            print(f"✓ {key}: {filepath} ({size:.2f} GB)")
    
    print("\nNext step: Run 02_data_preprocessing.py")
    
    return downloaded_files

if __name__ == "__main__":
    main()
