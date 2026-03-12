#!/usr/bin/env python3
"""
Download supplementary expression files from GEO
"""

import os
import requests

EXTERNAL_DIR = '../data/external'
os.makedirs(EXTERNAL_DIR, exist_ok=True)

print("=" * 60)
print("Downloading GEO Supplementary Expression Files")
print("=" * 60)

# GSE78220 supplementary file
print("\n1. GSE78220 - Patient FPKM data")
url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.xlsx'
output_file = f'{EXTERNAL_DIR}/GSE78220_PatientFPKM.xlsx'

print(f"   URL: {url}")
print(f"   Output: {output_file}")
print("\n   Note: FTP download may not work via requests")
print("   Alternative: Use browser to download from:")
print(f"   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220")
print("   Then click 'Supplementary file' -> GSE78220_PatientFPKM.xlsx")

# Try HTTP instead of FTP
http_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.xlsx'
print(f"\n   Trying HTTP: {http_url}")

try:
    response = requests.get(http_url, timeout=120, stream=True)
    if response.status_code == 200:
        with open(output_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"   [OK] Downloaded: {output_file}")
        print(f"   Size: {os.path.getsize(output_file) / 1024 / 1024:.2f} MB")
    else:
        print(f"   [X] Status: {response.status_code}")
except Exception as e:
    print(f"   [X] Error: {e}")

# GSE91061 supplementary files
print("\n2. GSE91061 - Check for supplementary files")
print("   Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061")
print("   Look for: counts, TPM, or expression matrix files")

print("\n" + "=" * 60)
print("Manual Download Instructions")
print("=" * 60)
print("""
If automatic download fails, please manually download:

1. GSE78220:
   - Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220
   - Scroll to "Supplementary file"
   - Download: GSE78220_PatientFPKM.xlsx
   - Save to: data/external/

2. GSE91061:
   - Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061
   - Check supplementary files for expression data
   - Download and save to: data/external/

3. Alternative - cBioPortal (IMvigor210):
   - Go to: http://www.cbioportal.org/study?id=blca_iatlas_imvigor210_2017
   - Download RNA-seq data
""")

print("\n[OK] Download script complete!")
