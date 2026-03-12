#!/usr/bin/env python3
"""
Download GSE91061 supplementary expression data
"""

import os
import requests

EXTERNAL_DIR = '../data/external'
os.makedirs(EXTERNAL_DIR, exist_ok=True)

print("=" * 60)
print("Downloading GSE91061 Expression Data")
print("=" * 60)

# Try to download supplementary files
urls_to_try = [
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_Riaz_TPM.txt.gz',
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_TPM.txt.gz',
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_counts.txt.gz',
]

for url in urls_to_try:
    filename = url.split('/')[-1]
    output_file = f'{EXTERNAL_DIR}/{filename}'
    
    print(f"\nTrying: {url}")
    try:
        response = requests.get(url, timeout=120, stream=True)
        if response.status_code == 200:
            with open(output_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"  [OK] Downloaded: {filename}")
            print(f"  Size: {os.path.getsize(output_file) / 1024 / 1024:.2f} MB")
            break
        else:
            print(f"  [X] Status: {response.status_code}")
    except Exception as e:
        print(f"  [X] Error: {e}")

print("\n" + "=" * 60)
print("Alternative: Manual Download")
print("=" * 60)
print("""
If automatic download fails:

1. Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061
2. Scroll to "Supplementary file" section
3. Download files like:
   - GSE91061_Riaz_TPM.txt.gz
   - GSE91061_gene_counts.txt.gz
   - Any .txt or .csv with expression data

4. Save to: data/external/

Note: GSE91061 has 109 samples with pre/on treatment pairs
""")
