#!/usr/bin/env python3
"""
Download processed expression data from GEO
"""

import os
import pandas as pd
import requests
import gzip
import shutil

EXTERNAL_DIR = '../data/external'
os.makedirs(EXTERNAL_DIR, exist_ok=True)

print("=" * 60)
print("Downloading Processed GEO Data")
print("=" * 60)

# GSE78220 - Try to download series matrix with expression
print("\n1. GSE78220 - Checking for expression data...")

# GEO series matrix URLs
geo_urls = {
    'GSE78220': [
        'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/matrix/GSE78220_series_matrix.txt.gz',
        'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/',
    ],
    'GSE91061': [
        'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/matrix/GSE91061_series_matrix.txt.gz',
        'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/',
    ]
}

# Try to download GSE78220 series matrix
url = geo_urls['GSE78220'][0]
output_file = f'{EXTERNAL_DIR}/GSE78220_series_matrix_full.txt.gz'

print(f"   Downloading: {url}")
try:
    response = requests.get(url, timeout=60, stream=True)
    if response.status_code == 200:
        with open(output_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"   [OK] Saved: {output_file}")
        
        # Check if it has expression data
        with gzip.open(output_file, 'rt', encoding='utf-8') as f:
            lines = f.readlines()
            has_expression = any('!series_matrix_table_begin' in line for line in lines)
            if has_expression:
                print("   [OK] Contains expression data!")
            else:
                print("   [INFO] No expression table found")
    else:
        print(f"   [X] Status: {response.status_code}")
except Exception as e:
    print(f"   [X] Error: {e}")

# Try to download GSE91061 series matrix
print("\n2. GSE91061 - Checking for expression data...")
url = geo_urls['GSE91061'][0]
output_file = f'{EXTERNAL_DIR}/GSE91061_series_matrix_full.txt.gz'

print(f"   Downloading: {url}")
try:
    response = requests.get(url, timeout=60, stream=True)
    if response.status_code == 200:
        with open(output_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"   [OK] Saved: {output_file}")
        
        # Check if it has expression data
        with gzip.open(output_file, 'rt', encoding='utf-8') as f:
            lines = f.readlines()
            has_expression = any('!series_matrix_table_begin' in line for line in lines)
            if has_expression:
                print("   [OK] Contains expression data!")
            else:
                print("   [INFO] No expression table found")
    else:
        print(f"   [X] Status: {response.status_code}")
except Exception as e:
    print(f"   [X] Error: {e}")

print("\n" + "=" * 60)
print("Alternative Sources")
print("=" * 60)
print("""
If GEO series matrix doesn't have expression data, try:

1. GEO Supplementary Files:
   - Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220
   - Check "Supplementary file" section
   - Look for: *_counts.txt, *_TPM.txt, *_expression.txt

2. cBioPortal (IMvigor210):
   - Visit: http://www.cbioportal.org/study?id=blca_iatlas_imvigor210_2017
   - Download RNA-seq data

3. Published Papers:
   - GSE78220: https://www.cell.com/cell/fulltext/S0092-8674(16)30215-X
   - Check supplementary materials

4. Gemma Database:
   - https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=10486
   - Pre-processed expression data
""")

print("\n[OK] Download attempt complete!")
