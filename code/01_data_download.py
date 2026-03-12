#!/usr/bin/env python3
"""
Step 1: Data Download for FerroScore-Immuno
Download TCGA and GTEx data from UCSC Xena
"""

import os
import pandas as pd
import numpy as np
import requests
from tqdm import tqdm
import gzip
import shutil

# Create directories
os.makedirs('../data/raw', exist_ok=True)
os.makedirs('../data/processed', exist_ok=True)

# UCSC Xena data hub - updated URLs
XENA_HUB = "https://gdc-xenagoatskins.s3.amazonaws.com"

# TCGA Pan-Cancer datasets - updated URLs
DATASETS = {
    'tcga_expression': {
        'url': f'{XENA_HUB}/download/TCGA.PANCAN.sampleMap%2FHiSeqV2_PANCAN.gz',
        'filename': 'TCGA_HiSeqV2_PANCAN.gz',
        'description': 'TCGA RNA-seq (HiSeqV2, log2)'
    },
    'tcga_phenotype': {
        'url': 'https://xena.ucsc.edu/public-hubs/TCGA/TCGA_phenotype_denseDataOnlyDownload.tsv.gz',
        'filename': 'TCGA_phenotype.tsv.gz',
        'description': 'TCGA clinical data'
    },
    'gtex_expression': {
        'url': f'{XENA_HUB}/download/GTEX_phenotype%2FGTEX_expression.gz',
        'filename': 'GTEX_expression.gz',
        'description': 'GTEx RNA-seq'
    },
    'gtex_phenotype': {
        'url': 'https://xena.ucsc.edu/public-hubs/GTEX/GTEX_phenotype.gz',
        'filename': 'GTEX_phenotype.tsv.gz',
        'description': 'GTEx clinical data'
    },
    'survival_data': {
        'url': 'https://xena.ucsc.edu/public-hubs/TCGA/survivalTCGA_pancan.txt',
        'filename': 'TCGA_survival.txt',
        'description': 'TCGA Pan-Cancer survival'
    }
}

def download_file(url, filename, description):
    """Download file with progress bar"""
    filepath = f'../data/raw/{filename}'
    
    if os.path.exists(filepath):
        print(f"[OK] {description} already exists")
        return filepath
    
    print(f"\nDownloading {description}...")
    print(f"  URL: {url}")
    
    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        
        with open(filepath, 'wb') as f:
            if total_size > 0:
                with tqdm(total=total_size, unit='B', unit_scale=True) as pbar:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
            else:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
        
        size_mb = os.path.getsize(filepath) / (1024**2)
        print(f"  [OK] Saved to {filepath} ({size_mb:.1f} MB)")
        return filepath
        
    except Exception as e:
        print(f"  [ERROR] {e}")
        return None

def decompress_gz(filepath):
    """Decompress .gz file"""
    if filepath.endswith('.gz'):
        output_path = filepath[:-3]
        if os.path.exists(output_path):
            print(f"  [OK] Decompressed file already exists")
            return output_path
        
        print(f"  Decompressing...")
        try:
            with gzip.open(filepath, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f"  [OK] Decompressed to {output_path}")
            return output_path
        except Exception as e:
            print(f"  [ERROR] Decompression failed: {e}")
            return filepath
    return filepath

def alternative_download_methods():
    """Provide alternative download methods if direct download fails"""
    print("\n" + "=" * 60)
    print("Alternative Download Methods")
    print("=" * 60)
    print("""
If automatic download fails, please manually download from:

1. UCSC Xena Browser: https://xena.ucsc.edu/
   - Go to "Data Hubs" -> "TCGA" 
   - Download: "TCGA Pan-Cancer (PANCAN)"
   
2. Specific files needed:
   - tcga_RSEM_gene_tpm (or HiSeqV2_PANCAN)
   - TCGA_phenotype_denseDataOnlyDownload.tsv
   - survivalTCGA_pancan.txt
   - gtex_RSEM_gene_tpm (optional)
   
3. Place files in: D:\Research\FerroScore-Radio\data\raw\

4. Then run: python 02_data_preprocessing.py
""")

def main():
    """Main download function"""
    print("=" * 60)
    print("FerroScore-Radio: Data Download")
    print("=" * 60)
    print("\nNote: TCGA data is ~2-3 GB and may take 10-30 minutes to download")
    print("      depending on your internet speed.\n")
    
    downloaded_files = {}
    failed_downloads = []
    
    for key, info in DATASETS.items():
        filepath = download_file(info['url'], info['filename'], info['description'])
        if filepath:
            downloaded_files[key] = filepath
            # Decompress if needed
            if filepath.endswith('.gz'):
                decompressed = decompress_gz(filepath)
                if decompressed != filepath:
                    downloaded_files[key + '_decompressed'] = decompressed
        else:
            failed_downloads.append(key)
    
    # Summary
    print("\n" + "=" * 60)
    print("Download Summary")
    print("=" * 60)
    
    if downloaded_files:
        print(f"\nSuccessfully downloaded {len(downloaded_files)} files:")
        for key, filepath in downloaded_files.items():
            if os.path.exists(filepath):
                size_mb = os.path.getsize(filepath) / (1024**2)
                print(f"  [OK] {key}: {size_mb:.1f} MB")
    
    if failed_downloads:
        print(f"\nFailed downloads ({len(failed_downloads)}):")
        for key in failed_downloads:
            print(f"  [FAIL] {key}")
        
        print("\n" + "=" * 60)
        alternative_download_methods()
        return
    
    print("\n" + "=" * 60)
    print("All downloads complete!")
    print("=" * 60)
    print("\nNext step: Run 02_data_preprocessing.py")
    
    return downloaded_files

if __name__ == "__main__":
    main()
