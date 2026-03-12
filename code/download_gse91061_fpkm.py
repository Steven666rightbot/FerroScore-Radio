#!/usr/bin/env python3
"""
Download GSE91061 FPKM data
"""

import requests
import os

EXTERNAL_DIR = '../data/external'
os.makedirs(EXTERNAL_DIR, exist_ok=True)

print("=" * 60)
print("Downloading GSE91061 FPKM Data")
print("=" * 60)

url = 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE91061&format=file&file=GSE91061%5FBMS038109Sample%2Ehg19KnownGene%2Efpkm%2Ecsv%2Egz'
output_file = f'{EXTERNAL_DIR}/GSE91061_fpkm.csv.gz'

print(f"\nDownloading from: {url}")
print(f"Saving to: {output_file}")

try:
    response = requests.get(url, timeout=300, stream=True)
    print(f"Status: {response.status_code}")
    
    if response.status_code == 200:
        with open(output_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        
        file_size = os.path.getsize(output_file) / 1024 / 1024
        print(f"\n[OK] Downloaded: {output_file}")
        print(f"Size: {file_size:.2f} MB")
        
        if file_size < 1:
            print("\n[WARNING] File size is small, may be an error page")
            print("Please download manually from:")
            print("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061")
    else:
        print(f"\n[X] Download failed: Status {response.status_code}")
        print("Please download manually")
        
except Exception as e:
    print(f"\n[X] Error: {e}")
    print("Please download manually from GEO website")
